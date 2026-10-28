"""One-stream annotation-only index builder. Root owns source/files/once/release."""
import hashlib,json,os,re,sqlite3,time
from pathlib import Path
from accepted_manifest_receiver import manifest_plan,SUPPORTED,SOURCE_ONLY,HOLD
from index_protocol import SCHEMA,BIN,DDL,need,encoded,manifest_hash,rowframe,limits_check,file_hash,Guard,configure,no_sort

class IndexHold(ValueError):
 def __init__(self,msg,state):super().__init__(msg);self.state=state

def parse_id(attrs):
 found=False;value=b''
 for item in attrs.split(b';'):
  k,sep,v=item.partition(b'=')
  if sep and k==b'ID':
   need(not found,'duplicate retained ID attribute');found=True;value=v
 return value

def fields(line,prefix,keys):
 parts=line.split(b'\t');need(parts[0]==prefix,'footer tag');d={}
 for item in parts[1:]:
  k,sep,v=item.partition(b'=');need(sep and k not in d,'footer duplicate/malformed key');d[k]=v
 need(set(d)==set(keys),'exact source footer fields')
 need(all(v.isdigit() for v in d.values()),'source integer fields');return {k:int(v) for k,v in d.items()}

def feature(line):
 f=line.split(b'\t');need(len(f)==9 and f[0] and f[2] and f[3].isdigit() and f[4].isdigit(),'GFF3 feature syntax')
 a,b=int(f[3]),int(f[4]);need(1<=a<=b<2**53 and f[6] in [b'+',b'-',b'.',b'?'],'exact integer coordinates/source strand')
 return f,a-1,b

def build(stream,manifest,db_path,temp_dir,limits):
 p,_,windows,nas,mh,_=manifest_plan(manifest);limits_check(limits,manifest);need(p['distinct_feature_types']<=4096,'bounded type counter dictionary')
 db=Path(db_path);temp=Path(temp_dir);g=Guard(db,temp,limits);c=None
 count={'read_bytes':0,'physical_lines':0,'feature_records':0,'retained_rows':0,'sequence_region_directives':0,'comment_lines':0,'blank_lines':0,'bin_links':0,'query_index_hits':0,'candidate_rows':0,'logical_feature_window_pairs':0,'logical_event_keys':0,'logical_gene_window_keys':0,'merged_intervals':0,'complete_EOF':False,'physical_EOF_seen':False,'source_terminal_validated':False}
 digest=hashlib.sha256();rd=hashlib.sha256();types={};footers={};sourcecounts=None;upstreamEOF=None;headers=0;eof=False
 try:
  fd=os.open(db,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600);os.close(fd)
  c=sqlite3.connect(db);configure(c,temp,limits,g);c.executescript(DDL)
  c.execute('PRAGMA user_version=1');c.execute('INSERT INTO source VALUES(1,?,?,NULL)',('PENDING',mh));c.commit()
  for i,t in enumerate(manifest['targets'],1):
   c.execute('INSERT INTO targets VALUES(?,?,?)',(i,t['target_key'],encoded(t)))
   if t['coordinate_status']!=SUPPORTED:c.execute('INSERT INTO typed_na VALUES(?,?)',(i,encoded(t)))
  c.commit()
  while True:
   g.check();raw=stream.readline(min(limits['max_line_bytes']+1,limits['expected_input_bytes']-count['read_bytes']+1))
   need(type(raw) is bytes,'binary source stream')
   if not raw:count['physical_EOF_seen']=True;break
   count['read_bytes']+=len(raw);count['physical_lines']+=1;digest.update(raw)
   need(count['read_bytes']<=limits['expected_input_bytes'] and len(raw)<=limits['max_line_bytes'] and raw.endswith(b'\n') and b'\r' not in raw,'exact bounded complete-LF source')
   need(not eof,'source content after EOF');line=raw[:-1]
   if re.match(br'^##FASTA(?:\s|$)',line):raise ValueError('FASTA not structural projection')
   if re.match(br'^##gff-version(?:\s|$)',line):
    headers+=1;need(line==b'##gff-version 3' and headers==1,'exact unique GFF3 header');count['comment_lines']+=1;continue
   if re.match(br'^##sequence-region(?:\s|$)',line):
    v=line.split();need(len(v)==4 and v[2].isdigit() and v[3].isdigit() and 1<=int(v[2])<=int(v[3]),'sequence region syntax');count['sequence_region_directives']+=1;count['comment_lines']+=1;continue
   if line.startswith(b'#TYPE_COUNT\t'):
    v=line.split(b'\t');need(len(v)==4 and v[2].isdigit() and v[3].isdigit() and v[1] not in footers,'type footer');footers[v[1]]=(int(v[2]),int(v[3]));count['comment_lines']+=1;continue
   if line.startswith(b'#SOURCE_COUNTS\t'):
    need(sourcecounts is None,'duplicate source counts');sourcecounts=fields(line,b'#SOURCE_COUNTS',[b'source_lines',b'feature_records',b'emitted_feature_records',b'comment_lines',b'blank_lines',b'sequence_region_directives',b'fasta_directives',b'fasta_tail_lines']);count['comment_lines']+=1;continue
   if line.startswith(b'#EOF_OK\t'):
    upstreamEOF=fields(line,b'#EOF_OK',[b'version',b'source_lines',b'feature_records',b'emitted_feature_records']);eof=True;count['comment_lines']+=1;continue
   if line.startswith(b'#'):count['comment_lines']+=1;continue
   if not line.strip():count['blank_lines']+=1;continue
   f,a,b=feature(line);count['feature_records']+=1
   if f[2] not in types:need(sum(len(k) for k in types)+len(f[2])<=2**20,'bounded source type-label bytes')
   types[f[2]]=types.get(f[2],0)+1
   need(len(types)<=p['distinct_feature_types'] and count['feature_records']<=p['feature_records'],'source feature/type bounds')
   if f[2] not in [b'gene',b'pseudogene',b'exon']:continue
   rid=parse_id(f[8]);need(f[2]==b'exon' or rid,'gene/pseudogene raw ID required')
   count['retained_rows']+=1;need(count['retained_rows']<=limits['max_retained_rows'],'retained-row cap')
   fid=count['feature_records'];phy=count['physical_lines'];c.execute('INSERT INTO features VALUES(?,?,?,?,?,?,?,?,?)',(fid,phy,f[0],f[2],rid,a,b,f[6],raw));rd.update(rowframe(fid,phy,raw))
   lo=a//BIN;hi=(b-1)//BIN;need(hi-lo+1<=limits['max_bins_per_feature'],'per-feature bin fanout cap')
   for bi in range(lo,hi+1):
    count['bin_links']+=1;need(count['bin_links']<=limits['max_bin_links'],'total index fanout cap');c.execute('INSERT INTO bins VALUES(?,?,?)',(f[0],bi,fid))
   if count['retained_rows']%limits['transaction_rows']==0:c.commit();g.check()
  need(eof and upstreamEOF is not None and sourcecounts is not None and headers==1,'missing complete unique source terminal/header')
  need(count['read_bytes']==p['bytes'] and digest.hexdigest()==limits['expected_input_sha256'],'exact source bytes/SHA')
  for k in ['physical_lines','feature_records','sequence_region_directives']:need(count[k]==p[k],'complete projection counter '+k)
  need(count['comment_lines']==p['sequence_region_directives']+p['distinct_feature_types']+6 and count['blank_lines']==0,'F+Q+D+6 exact physical closure')
  need(set(types)==set(footers) and len(types)==p['distinct_feature_types'],'all original types present')
  for kind,n in types.items():need(footers[kind]==(n,n),'all type counts equal emitted');c.execute('INSERT INTO feature_types VALUES(?,?,?,?)',(kind,n,n,n))
  F=p['feature_records'];src=p['source_lines']
  need(sourcecounts[b'source_lines']==src and sourcecounts[b'feature_records']==sourcecounts[b'emitted_feature_records']==F and sourcecounts[b'source_lines']==F+sourcecounts[b'comment_lines']+sourcecounts[b'blank_lines']+sourcecounts[b'fasta_tail_lines'] and sourcecounts[b'sequence_region_directives']==p['sequence_region_directives'],'upstream original counts closure')
  need(sourcecounts[b'fasta_directives'] in [0,1] and (not sourcecounts[b'fasta_tail_lines'] or sourcecounts[b'fasta_directives']==1) and sourcecounts[b'comment_lines']>=sourcecounts[b'sequence_region_directives']+sourcecounts[b'fasta_directives'],'upstream directives')
  need(upstreamEOF=={b'version':2,b'source_lines':src,b'feature_records':F,b'emitted_feature_records':F},'exact upstream EOF')
  count['source_terminal_validated']=True
  c.commit();query='SELECT f.* FROM work_candidates q CROSS JOIN features f ON f.id=q.id ORDER BY q.id';no_sort(c,query)
  for w in windows:
   g.check();t=next(x for x in manifest['targets'] if x['target_key']==w['target'])
   c.execute('DELETE FROM work_candidates');c.execute('DELETE FROM work_ids');c.execute('DELETE FROM work_events')
   seq=t['seqid'].encode('ascii');lo=w['requested_start0'];hi=w['requested_end0'];target_id=t['source_gene_ID'].encode('ascii')
   records=0;bad=False
   for sq,a,b in c.execute("SELECT seqid,start0,end0 FROM features INDEXED BY feature_rawid WHERE raw_id=? AND kind IN (x'67656e65',x'70736575646f67656e65') ORDER BY id",(target_id,)):
    records+=1;bad=bad or sq!=seq or a!=t['full_original_gene_start1']-1 or b!=t['full_original_gene_end1']
   candidates=0
   for (fid,) in c.execute('SELECT feature_id FROM bins WHERE seqid=? AND bin BETWEEN ? AND ?',(seq,lo//BIN,(hi-1)//BIN)):
    count['query_index_hits']+=1;need(count['query_index_hits']<=limits['max_query_index_hits'],'query index-hit cap')
    if c.execute('INSERT OR IGNORE INTO work_candidates VALUES(?)',(fid,)).rowcount:candidates+=1;need(candidates<=limits['max_query_candidates'],'bounded query-candidate cap')
   count['candidate_rows']+=candidates;genes=inner=pairs=events=0
   for fid,phy,sq,kind,rid,a,b,strand,raw in c.execute(query):
    g.check()
    if a>=hi or b<=lo:continue
    pairs+=1;count['logical_feature_window_pairs']+=1;need(count['logical_feature_window_pairs']<=limits['max_feature_window_pairs'],'logical row/window work cap')
    if kind==b'exon':
     for pos,delta in [(max(a,lo),1),(min(b,hi),-1)]:
      if c.execute('INSERT OR IGNORE INTO work_events VALUES(?,0)',(pos,)).rowcount:events+=1
      c.execute('UPDATE work_events SET delta=delta+? WHERE position=?',(delta,pos))
    else:
     inside=int(a<hi-(hi-lo)//4 and b>lo+(hi-lo)//4)
     if c.execute('INSERT OR IGNORE INTO work_ids VALUES(?,?)',(rid,inside)).rowcount:genes+=1;inner+=inside
   active=0;start=None;bp=segments=0
   for pos,delta in c.execute('SELECT position,delta FROM work_events ORDER BY position'):
    before=active;active+=delta;need(active>=0,'negative ALLexon event state')
    if before==0 and active>0:start=pos
    if before>0 and active==0:
     segments+=1;bp+=pos-start;count['merged_intervals']+=1;need(count['merged_intervals']<=limits['max_merged_intervals'],'merged interval cap');c.execute('INSERT INTO exons VALUES(?,?,?,?)',(w['window_index'],segments,start,pos))
   need(active==0,'unclosed ALLexon union')
   status=SOURCE_ONLY if records==1 and not bad else HOLD
   c.execute('INSERT INTO windows VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,0)',(w['window_index'],w['target'],w['radius'],seq,lo,hi,t['RNA_start1'],t['RNA_end1'],t['full_original_gene_start1'],t['full_original_gene_end1'],target_id,t['strand'].encode(),status,genes,inner,w['full_original_gene_length'],bp,segments,records,pairs,events))
   count['logical_event_keys']+=events;count['logical_gene_window_keys']+=genes;c.commit();g.check()
  for i in range(2,len(windows)+1,2):
   a=c.execute('SELECT inner_half_gene_count FROM windows WHERE window_index=?',(i-1,)).fetchone()[0];b=c.execute('SELECT gene_count_including_target FROM windows WHERE window_index=?',(i,)).fetchone()[0];need(a==b,'INNER50 must equal GENES25; no duplicate-ID rescue')
  for name in ['work_candidates','work_ids','work_events']:c.execute('DELETE FROM '+name)
  count['complete_EOF']=True
  receipt=dict(count,schema=SCHEMA,status='BUILT_UNACCEPTED',manifest_sha256=mh,input_sha256=digest.hexdigest(),logical_retained_row_sha256=rd.hexdigest(),upstream_source_counts={k.decode():v for k,v in sourcecounts.items()},upstream_EOF={k.decode():v for k,v in upstreamEOF.items()},windows=len(windows),typed_NA=len(nas),primary_callable_or_model_eligible=False)
  c.execute('UPDATE source SET state=?,receipt_json=? WHERE id=1',('CLOSED',encoded(receipt)));c.commit();g.check();c.close();c=None
  receipt['db_file']=file_hash(db,limits['main_bytes']);g.check();receipt['resource_peak_sampled']=g.peak;receipt['seconds']=time.monotonic()-g.started
  return receipt
 except BaseException as e:
  count.update(input_sha256=digest.hexdigest(),logical_retained_row_sha256=rd.hexdigest(),automatic_retry=False,db_path=str(db),complete_EOF=False)
  if c:
   try:c.rollback();c.execute('UPDATE source SET state=?,receipt_json=? WHERE id=1',('HOLD',encoded(dict(count,reason=str(e)))));c.commit()
   except BaseException:pass
   c.close()
  raise IndexHold(type(e).__name__+': '+str(g.fault or e),count) from e
