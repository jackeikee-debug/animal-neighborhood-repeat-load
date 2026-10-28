"""One new, bounded occupancy/candidate-closure audit; no phylogeny execution."""
import collections
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

def mcl_clusters(lines):
    active=False; index=None; members=[]
    for line in lines:
        s=line.strip()
        if not active:
            if s=='begin': active=True
            continue
        if not s or s.startswith('#'): continue
        if s==')':
            if index is not None: raise ValueError('unterminated MCL cluster')
            continue
        tokens=s.split()
        if index is None:
            if not tokens[0].isdigit(): raise ValueError('unexpected MCL index')
            index=int(tokens.pop(0))
        if '$' in tokens:
            p=tokens.index('$')
            if p!=len(tokens)-1: raise ValueError('MCL trailing fields')
            members.extend(tokens[:p])
            if len(members)!=len(set(members)): raise ValueError('duplicate member in cluster')
            yield f'OG{index:07d}',members
            index=None;members=[]
        else: members.extend(tokens)
    if index is not None: raise ValueError('truncated MCL source')

def bind_anchor_ids(lines,anchors):
    by_anchor={a:set() for a in anchors};by_uid={}
    for line in lines:
        if not line.strip(): continue
        uid,value=line.rstrip('\n').split(': ',1)
        sid=int(uid.split('_',1)[0]);key=(sid,value.split()[0])
        if key in anchors:
            if uid in by_uid: raise ValueError('duplicate native UID')
            by_anchor[key].add(uid);by_uid[uid]=set(anchors[key])
    return {'by_anchor':by_anchor,'by_uid':by_uid}

def anchor_closure(clusters,by_uid):
    wanted=set(by_uid);selected=[];seen={u:[] for u in wanted};total=0
    for og,members in clusters:
        total+=1;hit=wanted.intersection(members)
        if not hit: continue
        groups=set()
        for u in hit: seen[u].append(og);groups.update(by_uid[u])
        selected.append({'native_OG':og,'members':members,'old_groups':sorted(groups),'anchor_uids':sorted(hit),'species_count':len({u.split('_',1)[0] for u in members})})
    return selected,seen

def alignment_stats(lines):
    canonical=set('ACDEFGHIKLMNPQRSTVWY');record=None;width=None;seen=set()
    def finish(r):
        nonlocal width
        if width is None: width=r['columns']
        if width!=r['columns']: raise ValueError('nonrectangular concatenation')
        r['noncanonical_characters']=dict(sorted(r['noncanonical_characters'].items()))
        return r
    for line in lines:
        s=line.strip()
        if not s: continue
        if s.startswith('>'):
            if record is not None: yield finish(record)
            sid=int(s[1:].split()[0])
            if sid in seen: raise ValueError('duplicate alignment species')
            seen.add(sid);record={'species_id':sid,'columns':0,'gap':0,'canonical20':0,'noncanonical':0,'noncanonical_characters':collections.Counter()}
        else:
            if record is None: raise ValueError('alignment data before header')
            for ch,n in collections.Counter(s.upper()).items():
                record['columns']+=n
                if ch in '-.': record['gap']+=n
                elif ch in canonical: record['canonical20']+=n
                else: record['noncanonical']+=n;record['noncanonical_characters'][ch]+=n
    if record is not None: yield finish(record)

def snap(p):
    s=Path(p).stat();return dict(dev=s.st_dev,ino=s.st_ino,mode=s.st_mode,size=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)

class HashedRaw(io.RawIOBase):
    def __init__(self,fp,check,progress):
        super().__init__();self.fp=fp;self.check=check;self.progress=progress;self.n=0;self.h=hashlib.sha256();self.mark=0
    def readable(self):return True
    def readinto(self,b):
        data=self.fp.read(len(b));n=len(data)
        if n:self.h.update(data);self.n+=n;b[:n]=data
        self.check()
        if self.n-self.mark>=4*1024**2:self.mark=self.n;self.progress(self.n)
        return n

class Sources:
    def __init__(self,cfg,check,events):
        self.sources={s['name']:s for s in cfg['sources']};self.check=check;self.events=events;self.done={};self.used=set()
    def lines(self,name):
        if name in self.used:raise ValueError('source purpose attempted twice')
        self.used.add(name);s=self.sources[name];p=Path(s['path']);assert snap(p)==s['closed_stat']
        self.events(name,'STARTED_BEFORE_SOURCE_BODY',{'path':str(p),'durable_compressed_bytes':0})
        count=0;raw=None
        try:
            with p.open('rb') as fp:
                raw=HashedRaw(fp,self.check,lambda n:self.events(name,'PROGRESS_DURABLE_LOWER_BOUND',{'compressed_bytes_lower_bound':n,'exact_final_bytes_known':False}))
                buffered=io.BufferedReader(raw,buffer_size=65536)
                content=gzip.GzipFile(fileobj=buffered) if s.get('gzip') else buffered
                text=io.TextIOWrapper(content,encoding='utf-8',newline='')
                for line in text:count+=1;yield line
                # Drain the raw envelope as well; gzip verifies CRC and EOF above.
                if raw.n!=s['bytes']:raise ValueError('compressed byte accounting')
                digest=raw.h.hexdigest()
                if s.get('sha256') and digest!=s['sha256']:raise ValueError('source hash mismatch')
                assert snap(p)==s['closed_stat']
                r={'compressed_bytes':raw.n,'inline_prefix_sha256':digest,'detail':{'full_EOF':True,'text_lines':count},'source_stat_unchanged':True}
                self.done[name]=r;self.events(name,'COMPLETE_EOF_SHA_STAT_ACCEPTED',r)
        except BaseException as e:
            self.events(name,'UNKNOWN_FAILED_PARTIAL_NO_RETRY',{'compressed_bytes_lower_bound':raw.n if raw else 0,'exact_final_bytes_known':False,'reason':type(e).__name__+':'+str(e)})
            raise

def write_json(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')

def build(cfg,out,check,events):
    out=Path(out);out.mkdir();src=Sources(cfg,check,events)
    species={}
    for line in src.lines('species'):
        sid,name=line.rstrip().split(': ',1);species[int(sid)]=name.removesuffix('.fa')
    assert len(species)==cfg['expected_species']
    roles={k:int(v) for k,v in cfg['old_role_species_ids'].items()}
    anchors=collections.defaultdict(set);old_rows=[]
    for r in csv.DictReader(src.lines('old_bridge'),delimiter='\t'):
        key=(roles[r['species']],r['protein_id']);anchors[key].add(r['original_group']);old_rows.append(r)
    assert len(old_rows)==cfg['expected_old_rows'] and len({r['original_group'] for r in old_rows})==cfg['expected_old_groups']
    selected={l.strip() for l in src.lines('selected_OGs') if l.strip()};assert len(selected)==cfg['expected_selected_OGs']
    stats=list(alignment_stats(src.lines('concat')));assert {r['species_id'] for r in stats}==set(species)
    for r in stats:r['species_name']=species[r['species_id']]
    write_json(out/'species-alignment-occupancy.json',stats)
    reader=csv.DictReader(src.lines('gene_count'),delimiter='\t');assert reader.fieldnames[1:-1]==list(species.values()) and reader.fieldnames[-1]=='Total'
    acc={sid:{'absent':0,'single':0,'multi':0,'copies':0} for sid in species};found=set()
    with (out/'selected-OG-copy-occupancy.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['native_OG','species_id','copies','copy_state'])
        for r in reader:
            og=r['Orthogroup']
            if og not in selected:continue
            assert og not in found;found.add(og)
            for sid,name in species.items():
                n=int(r[name]);state='absent' if n==0 else 'single' if n==1 else 'multi';acc[sid][state]+=1;acc[sid]['copies']+=n;w.writerow([og,sid,n,state])
    assert found==selected
    with (out/'species-occupancy.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['species_id','species','absent','single','multi','copies','columns','gap','canonical20','noncanonical'])
        statby={r['species_id']:r for r in stats}
        for sid in species:w.writerow([sid,species[sid],*[acc[sid][x] for x in ['absent','single','multi','copies']],*[statby[sid][x] for x in ['columns','gap','canonical20','noncanonical']]])
    bound=bind_anchor_ids(src.lines('sequence_ids'),anchors)
    clusters_seen=0
    def indexed():
        nonlocal clusters_seen
        for og,members in mcl_clusters(src.lines('clusters')):
            assert int(og[2:])==clusters_seen,'native cluster indices contiguous';clusters_seen+=1;yield og,members
    closure,seen=anchor_closure(indexed(),bound['by_uid']);assert clusters_seen==cfg['expected_native_OGs']
    by_og={x['native_OG']:x for x in closure};trees=Path(cfg['trees_directory']);nbytes=0;missing=[]
    with (out/'anchor-family-closure.tsv').open('x') as f,(out/'full-OG-closure-members.jsonl').open('x') as full:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['native_OG','member_count','species_count','old_group_count','excluded725','tree_status','tree_bytes'])
        for x in closure:
            og=x['native_OG'];p=trees/(og+'_tree_id.txt');exists=p.exists();size=p.stat().st_size if exists else 0;nbytes+=size
            status='EXPLICIT_OG725_EXCLUSION' if og=='OG0000725' else 'EXISTING_UNROOTED_TREE' if exists else 'NATIVE_SMALL_OR_SINGLE_SPECIES_PATH' if len(x['members'])<=3 or x['species_count']==1 else 'MISSING_REQUIRED_TREE'
            if status=='MISSING_REQUIRED_TREE':missing.append(og)
            x.update(tree_path=str(p) if exists else None,tree_bytes=size,tree_status=status)
            w.writerow([og,len(x['members']),x['species_count'],len(x['old_groups']),og=='OG0000725',status,size]);full.write(json.dumps(x,separators=(',',':'))+'\n')
    states=collections.Counter()
    with (out/'all-original-anchor-membership.tsv').open('x') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['original_group','role','species_id','protein_id','native_uids','native_OGs','state'])
        for r in old_rows:
            sid=roles[r['species']];uids=bound['by_anchor'][(sid,r['protein_id'])];ogs=sorted({og for u in uids for og in seen[u]})
            state='NO_NATIVE_SEQUENCE_ID' if not uids else 'NO_NATIVE_CLUSTER' if not ogs else 'AMBIGUOUS_SEQUENCE_ID' if len(uids)>1 else 'MULTIPLE_NATIVE_CLUSTERS' if len(ogs)>1 else 'EXPLICIT_OG725_EXCLUSION' if ogs==['OG0000725'] else 'BOUND'
            states[state]+=1;w.writerow([r['original_group'],r['species'],sid,r['protein_id'],','.join(sorted(uids)),','.join(ogs),state])
    # OG partitioning is species-tree independent. All candidates for an old
    # protein lie in its complete native cluster(s); forward/reverse collision
    # tests later use all members, never only the current successful anchors.
    report={'status':'INPUT_AUDIT_COMPLETE','species':len(species),'selected_OGs':len(selected),'alignment_columns':stats[0]['columns'],'original_groups':len({r['original_group'] for r in old_rows}),'original_rows':len(old_rows),'native_clusters_scanned':clusters_seen,'complete_OG_closure':len(closure),'closure_members':sum(len(x['members']) for x in closure),'closure_existing_tree_bytes':nbytes,'required_tree_missing':missing,'anchor_states':dict(states),'closure_complete_for_all_bound_anchors':not missing and all(len(v)==1 for v in seen.values()),'all1896_including_unresolved_candidates_examined':True,'OG725_excluded':True,'per_locus_trimmed_site_occupancy':'UNESTABLISHED_NO_FINAL_PARTITION_MAP','source_receipts':src.done,'new_scientific_tools':0}
    write_json(out/'closure-completeness.json',report)
    write_json(out/'terminal.json',report)
    return {'numeric_call_counts':{'input_occupancy_and_closure_audit':1,'phylogeny_HOG_models':0},'source_receipts':src.done}
