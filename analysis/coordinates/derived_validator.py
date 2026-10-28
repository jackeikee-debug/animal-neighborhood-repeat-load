"""Independent bounded derived-SQLite validation. Never reads original stream.

Count/union implementation is independent: raw rows -> indexed candidates ->
SQL first-rawID grouping and sorted interval merge (not producer event sweep).
"""
import hashlib,json,os,sqlite3,time
from pathlib import Path
from accepted_manifest_receiver import manifest_plan,SUPPORTED,SOURCE_ONLY,HOLD
from index_protocol import SCHEMA,BIN,need,encoded,rowframe,limits_check,file_hash,Guard,configure,no_sort,readonly


def validate(db_path,manifest,expected,temp_dir,limits):
 p,_,ws,nas,mh,_=manifest_plan(manifest);limits_check(limits,manifest)
 need(expected.get('schema')==SCHEMA and expected.get('status')=='BUILT_UNACCEPTED' and expected.get('complete_EOF') is True,'closed unaccepted build receipt required')
 g=Guard(db_path,temp_dir,limits)
 need(file_hash(db_path,limits['main_bytes'])==expected['db_file'],'independent closed DB physical hash')
 g.check()
 c=readonly(db_path,cache_kib=4096);scratch=None
 try:
  need(c.execute('PRAGMA user_version').fetchone()[0]==1,'exact index user version')
  row=c.execute('SELECT state,manifest_hash,receipt_json FROM source WHERE id=1').fetchone();need(row and row[0]=='CLOSED' and row[1]==mh,'source schema/complete state')
  src=json.loads(row[2]);core={k:v for k,v in expected.items() if k not in ['db_file','resource_peak_sampled','seconds']};need(src==core and src['input_sha256']==limits['expected_input_sha256'] and src['physical_EOF_seen'] is True and src['source_terminal_validated'] is True,'source receipt/pins')
  need(src['physical_lines']==p['physical_lines']==p['feature_records']+p['sequence_region_directives']+p['distinct_feature_types']+6 and src['feature_records']==p['feature_records'] and src['read_bytes']==p['bytes'] and src['sequence_region_directives']==p['sequence_region_directives'] and src['blank_lines']==0,'full source F+Q+D+6 facts')
  types={};typeSum=0
  for kind,n,s,e in c.execute('SELECT * FROM feature_types'):
   need(n==s==e and n>0,'complete source types');types[kind]=n;typeSum+=n
  need(typeSum==p['feature_records'] and len(types)==p['distinct_feature_types'],'all type count closure')
  rd=hashlib.sha256();retained=links=0;previous=0;observed={}
  for fid,physical,seq,kind,rid,a,b,strand,raw in c.execute('SELECT * FROM features ORDER BY id'):
   g.check();need(fid>previous and physical>=fid and physical<=p['physical_lines'] and type(raw) is bytes and raw.endswith(b'\n'),'raw source ordinal/LF');previous=fid
   cols=raw[:-1].split(b'\t');need(len(cols)==9 and cols[0]==seq and cols[2]==kind and cols[3].isdigit() and cols[4].isdigit() and int(cols[3])-1==a and int(cols[4])==b and cols[6]==strand and 0<=a<b<2**53,'raw row exact BLOB/coordinate fields')
   found=[]
   for item in cols[8].split(b';'):
    k,sep,v=item.partition(b'=')
    if sep and k==b'ID':found.append(v)
   need(len(found)<=1 and rid==(found[0] if found else b'') and (kind==b'exon' or rid) and kind in [b'gene',b'pseudogene',b'exon'],'raw ID/type identity, no duplicate attributes')
   retained+=1;need(retained<=limits['max_retained_rows'],'derived row cap');observed[kind]=observed.get(kind,0)+1;rd.update(rowframe(fid,physical,raw))
   lo=a//BIN;hi=(b-1)//BIN;need(hi-lo+1<=limits['max_bins_per_feature'],'index fanout bound')
   cur=c.execute('SELECT seqid,bin FROM bins INDEXED BY bins_by_feature WHERE feature_id=? ORDER BY bin',(fid,))
   for bi in range(lo,hi+1):
    need(next(cur,None)==(seq,bi),'complete per-source-row bin identity');links+=1;need(links<=limits['max_bin_links'],'total bin fanout bound')
   need(next(cur,None) is None,'extra bin link')
  need(retained==src['retained_rows']==sum(types.get(k,0) for k in [b'gene',b'pseudogene',b'exon']) and rd.hexdigest()==src['logical_retained_row_sha256'] and links==src['bin_links']==c.execute('SELECT count(*) FROM bins').fetchone()[0],'logical raw-row/bin digest closure')
  cur=c.execute('SELECT ordinal,target,input_json FROM targets ORDER BY ordinal')
  for i,t in enumerate(manifest['targets'],1):need(next(cur,None)==(i,t['target_key'],encoded(t)),'exact all-target input ledger')
  need(next(cur,None) is None,'extra target')
  nc=c.execute('SELECT input_json FROM typed_na ORDER BY ordinal')
  for t in nas:need(next(nc,None)==(encoded(t),),'unchanged complete typed NA')
  need(next(nc,None) is None,'extra NA')
  # Own bounded validation scratch, never mutate the retained source/index.
  sp=Path(temp_dir)/'validation-scratch.sqlite';fd=os.open(sp,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600);os.close(fd)
  scratch=sqlite3.connect(sp);configure(scratch,temp_dir,limits,g,main=False,cache_kib=4096)
  scratch.executescript('CREATE TABLE candidates(id INTEGER PRIMARY KEY);CREATE TABLE genes(raw_id BLOB,fid INTEGER,a INTEGER,b INTEGER,PRIMARY KEY(raw_id,fid)) WITHOUT ROWID;CREATE TABLE intervals(fid INTEGER PRIMARY KEY,a INTEGER,b INTEGER);CREATE INDEX interval_order ON intervals(a,b);CREATE TABLE points(pos INTEGER PRIMARY KEY);')
  no_sort(scratch,'SELECT raw_id,min(fid) FROM genes GROUP BY raw_id');no_sort(scratch,'SELECT a,b FROM intervals INDEXED BY interval_order ORDER BY a,b')
  totals={'pairs':0,'events':0,'genes':0,'segments':0};targetMap={t['target_key']:t for t in manifest['targets']};hits=0
  columns=[x[1] for x in c.execute('PRAGMA table_info(windows)')]
  for w in ws:
   g.check();rawAudit=c.execute('SELECT * FROM windows WHERE window_index=?',(w['window_index'],)).fetchone();need(rawAudit is not None,'missing window');a=dict(zip(columns,rawAudit));t=targetMap[w['target']];seq=t['seqid'].encode('ascii');lo=w['requested_start0'];hi=w['requested_end0']
   for key in ['target','window_index','radius','requested_start0','requested_end0','full_original_gene_length']:need(a[key]==w[key],'frozen requested window binding')
   need(a['seqid_raw']==seq and a['RNA_start1']==t['RNA_start1'] and a['RNA_end1']==t['RNA_end1'] and a['full_original_gene_start1']==t['full_original_gene_start1'] and a['full_original_gene_end1']==t['full_original_gene_end1'] and a['source_gene_id_raw']==t['source_gene_ID'].encode() and a['strand_raw']==t['strand'].encode() and a['primary_callable_assigned']==0,'RNA/fullgene/rawID carry')
   for table in ['candidates','genes','intervals','points']:scratch.execute('DELETE FROM '+table)
   candidates=0
   for (fid,) in c.execute('SELECT feature_id FROM bins WHERE seqid=? AND bin BETWEEN ? AND ?',(seq,lo//BIN,(hi-1)//BIN)):
    hits+=1;need(hits<=limits['max_query_index_hits'],'validation query-hit cap')
    if scratch.execute('INSERT OR IGNORE INTO candidates VALUES(?)',(fid,)).rowcount:candidates+=1;need(candidates<=limits['max_query_candidates'],'validation candidate cap')
   pairCount=0
   for (fid,) in scratch.execute('SELECT id FROM candidates ORDER BY id'):
    sq,kind,rid,x,y=c.execute('SELECT seqid,kind,raw_id,start0,end0 FROM features WHERE id=?',(fid,)).fetchone()
    if x>=hi or y<=lo:continue
    pairCount+=1;totals['pairs']+=1;need(totals['pairs']<=limits['max_feature_window_pairs'],'validator logical pair cap')
    if kind==b'exon':
     l=max(x,lo);r=min(y,hi);scratch.execute('INSERT INTO intervals VALUES(?,?,?)',(fid,l,r));scratch.execute('INSERT OR IGNORE INTO points VALUES(?)',(l,));scratch.execute('INSERT OR IGNORE INTO points VALUES(?)',(r,))
    else:scratch.execute('INSERT INTO genes VALUES(?,?,?,?)',(rid,fid,x,y))
   genes=inner=0
   for rid,fid in scratch.execute('SELECT raw_id,min(fid) FROM genes GROUP BY raw_id'):
    x,y=scratch.execute('SELECT a,b FROM genes WHERE raw_id=? AND fid=?',(rid,fid)).fetchone();genes+=1;inner+=int(x<hi-(hi-lo)//4 and y>lo+(hi-lo)//4)
   expectedIntervals=c.execute('SELECT start0,end0 FROM exons WHERE window_index=? ORDER BY segment_index',(w['window_index'],));start=end=None;bp=segments=0
   for x,y in scratch.execute('SELECT a,b FROM intervals INDEXED BY interval_order ORDER BY a,b'):
    if end is None:start,end=x,y
    elif x<=end:end=max(end,y)
    else:
     need(next(expectedIntervals,None)==(start,end),'independent ALLexon union');bp+=end-start;segments+=1;start,end=x,y
   if end is not None:need(next(expectedIntervals,None)==(start,end),'independent last ALLexon interval');bp+=end-start;segments+=1
   need(next(expectedIntervals,None) is None,'extra retained union segment')
   events=scratch.execute('SELECT count(*) FROM points').fetchone()[0];records=0;bad=False
   for sq,x,y in c.execute("SELECT seqid,start0,end0 FROM features INDEXED BY feature_rawid WHERE raw_id=? AND kind IN (x'67656e65',x'70736575646f67656e65') ORDER BY id",(t['source_gene_ID'].encode(),)):
    records+=1;bad=bad or sq!=seq or x!=t['full_original_gene_start1']-1 or y!=t['full_original_gene_end1']
   status=SOURCE_ONLY if records==1 and not bad else HOLD
   for key,value in [('status',status),('target_gene_records',records),('gene_count_including_target',genes),('inner_half_gene_count',inner),('exon_union_bp',bp),('union_segments',segments),('logical_feature_window_pairs',pairCount),('logical_event_keys',events)]:need(a[key]==value,'independent scientific audit: '+key)
   totals['events']+=events;totals['genes']+=genes;totals['segments']+=segments;need(totals['segments']<=limits['max_merged_intervals'],'validation segment cap');scratch.commit();g.check()
  need(c.execute('SELECT count(*) FROM windows').fetchone()[0]==len(ws) and c.execute('SELECT count(*) FROM typed_na').fetchone()[0]==len(nas) and c.execute('SELECT count(*) FROM exons').fetchone()[0]==totals['segments'],'complete windows/NA/segments')
  for i in range(2,len(ws)+1,2):
   big=c.execute('SELECT requested_start0,requested_end0,inner_half_gene_count FROM windows WHERE window_index=?',(i-1,)).fetchone();small=c.execute('SELECT requested_start0,requested_end0,gene_count_including_target FROM windows WHERE window_index=?',(i,)).fetchone();need(big[2]==small[2],'INNER50 = GENES25')
   actual=c.execute('SELECT start0,end0 FROM exons WHERE window_index=? ORDER BY segment_index',(i,))
   for x,y in c.execute('SELECT start0,end0 FROM exons WHERE window_index=? ORDER BY segment_index',(i-1,)):
    l=max(x,small[0]);r=min(y,small[1])
    if l<r:need(next(actual,None)==(l,r),'A25 equals clipped A50')
   need(next(actual,None) is None,'extra A25 interval')
  need(totals['pairs']==src['logical_feature_window_pairs'] and totals['events']==src['logical_event_keys'] and totals['genes']==src['logical_gene_window_keys'] and totals['segments']==src['merged_intervals'],'global logical counters')
  for table in ['work_candidates','work_ids','work_events']:need(c.execute('SELECT count(*) FROM '+table).fetchone()[0]==0,'closed empty producer scratch')
  scratch.close();scratch=None;g.check();need(file_hash(db_path,limits['main_bytes'])==expected['db_file'],'DB unchanged during independent QA');g.check()
  return {'schema':SCHEMA,'status':'DERIVED_INDEX_VALIDATED_SOURCE_ONLY','manifest_sha256':mh,'db_file':expected['db_file'],'source_receipt_sha256':hashlib.sha256(encoded(src)).hexdigest(),'windows':len(ws),'typed_NA':len(nas),'merged_intervals':totals['segments'],'primary_callable_or_model_eligible':False,'original_source_reads':0,'resource_peak_sampled':g.peak,'seconds':time.monotonic()-g.started}
 finally:
  if scratch:scratch.close()
  c.close()
