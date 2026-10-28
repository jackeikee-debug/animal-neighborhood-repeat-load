#!/usr/bin/env python3
"""Prepare real-protein evidence; no inference, repeats, network, or automatic retries.

The CLI can run a real join ONLY with a Root exact release and fresh resource gates.
Input large tables are read solely by run_join(), never at import or test discovery.
"""
from pathlib import Path
from dataclasses import dataclass, field
import argparse, collections, csv, datetime, gzip, hashlib, itertools, json, math, os, re, sqlite3, sys, time
sys.dont_write_bytecode=True
CODES=('G','P','S')
COORD_STATES={'ACCEPTED_COORDINATE','ANNOTATION_ONLY','COORDINATE_PENDING','COORDINATE_CONFLICT','COORDINATE_EVIDENCE_NOT_SUPPLIED'}
DECISION_KEYS=('root_primary_conditional','mixed_N0_reference_only','N26_old_members_required','all_partner_unique_1to1','resolved_relevant_copy_filter','preserve_convergent_groups','no_optional_phylogenetic_fit_release')

def require_decisions(d):
 if any(d.get(k) is not True for k in DECISION_KEYS):raise ValueError('Missing explicit Root preoutcome node/filter decision')
 if d.get('animal_reference_node')!='N0' or d.get('old_reference_node')!='N26':raise ValueError('Root node contract differs; review a successor')

def tokens(s):
 if not s:return []
 x=s.split(', ')
 if any(not t or ',' in t or t!=t.strip() for t in x) or len(x)!=len(set(x)):raise ValueError('Malformed/duplicate native gene list')
 return x

def pair_status(seen,row_bad,forward_count,reverse_bad):
 if not seen:return 'UNOBSERVED'
 return 'COORTHOLOG' if row_bad or forward_count!=1 or reverse_bad else 'ONE_TO_ONE'

def classify_candidate(f):
 reasons=[]
 if f['source_identity']!='EXACT_UNIQUE':reasons.append(f['source_identity'])
 if not f['hog_same']:reasons.append('HOG_CONFLICT')
 if f['hog_multicopy']:reasons.append('HOG_MULTICOPY')
 for p in sorted(set(f['pairwise'])):
  if p!='ONE_TO_ONE':reasons.append('PAIRWISE_'+p)
 if f['tree']!='COMPATIBLE':reasons.append('GENE_TREE_'+f['tree'])
 if reasons:return 'UNQUALIFIED',sorted(set(reasons))
 c=f['coordinate']
 if c=='ACCEPTED_COORDINATE':return 'QUALIFIED_CANDIDATE',[]
 if c=='ANNOTATION_ONLY':return 'SUPPORTED_ANNOTATION_ONLY',[c]
 if c=='COORDINATE_CONFLICT':return 'SUPPORTED_COORDINATE_CONFLICT',[c]
 if c=='COORDINATE_EVIDENCE_NOT_SUPPLIED':return 'SUPPORTED_COORDINATE_EVIDENCE_NOT_SUPPLIED',[c]
 return 'SUPPORTED_COORDINATE_PENDING',[c]

def cell_status(statuses,g0_blocked):
 if g0_blocked:return 'G0_UNRESOLVED',[]
 passing=[x for x in statuses if x!='UNQUALIFIED']
 if len(passing)>1:return 'MULTIPLE_SUPPORTED_CANDIDATES',[]
 if len(passing)==1:return passing[0],[]
 return ('UNOBSERVED' if not statuses else 'NO_UNIQUE_COMPATIBLE_CANDIDATE'),[]

@dataclass(eq=False)
class Node:
 name:str=''
 length:float=0.0
 children:list=field(default_factory=list)
 parent:object=None
 def leaves(self):
  return [self] if not self.children else [x for c in self.children for x in c.leaves()]
 def walk(self):
  yield self
  for c in self.children:yield from c.walk()

def parse_newick(text):
 """Strict unquoted grammar of accepted native outputs, no tree fitting."""
 s=text.strip();i=0
 def label():
  nonlocal i
  start=i
  while i<len(s) and s[i] not in ',():;':i+=1
  x=s[start:i]
  if any(c in x for c in "[]'\"\n\r\t "):raise ValueError('Unsupported Newick label')
  return x
 def subtree():
  nonlocal i
  n=Node()
  if i<len(s) and s[i]=='(':
   i+=1
   while True:
    ch=subtree();ch.parent=n;n.children.append(ch)
    if i>=len(s):raise ValueError('Unclosed Newick')
    c=s[i];i+=1
    if c==')':break
    if c!=',':raise ValueError('Malformed Newick')
  n.name=label()
  if not n.children and not n.name:raise ValueError('Empty tip')
  if i<len(s) and s[i]==':':
   i+=1;n.length=float(label())
   if not math.isfinite(n.length):raise ValueError('Nonfinite branch')
  return n
 root=subtree()
 if s[i:]!=';':raise ValueError('Trailing Newick or missing terminator')
 tips=[x.name for x in root.leaves()]
 if len(set(tips))!=len(tips):raise ValueError('Duplicate tips')
 return root

def mrca(root,names):
 want=set(names);tips={x.name:x for x in root.leaves()}
 if not want or not want<=tips.keys():raise ValueError('Missing MRCA tips')
 n=tips[next(iter(want))]
 while not want<={x.name for x in n.leaves()}:n=n.parent
 return n

def tree_consistency(gene,species,names,tip_species):
 """Only the frozen visible relevant-copy constraint; no topology fitting/test."""
 leaves={x.name for x in gene.leaves()}
 if not set(names)<=leaves:return 'UNOBSERVED'
 if any(x not in tip_species for x in leaves):return 'IDENTITY_CONFLICT'
 species_names={x.name for x in species.leaves()}
 if not set(tip_species.values())<=species_names:return 'IDENTITY_CONFLICT'
 clade=mrca(gene,names);counts=collections.Counter(tip_species[x.name] for x in clade.leaves())
 relevant={tip_species[x] for x in names}
 return 'MULTICOPY_CONFLICT' if any(counts[s]>1 for s in relevant) else 'COMPATIBLE'

def text_open(p,mode='rt'):
 return gzip.open(p,mode,encoding='utf-8',newline='') if str(p).endswith('.gz') else open(p,mode,encoding='utf-8',newline='')

def load_old(path,expected=1896):
 result={};seen=set()
 with text_open(path) as f:
  reader=csv.DictReader(f,delimiter='\t')
  for row in reader:
   g=row['original_group'];c=row['species']
   if c not in CODES or c in result.setdefault(g,{}):raise ValueError('Old bridge duplicate/malformed anchor')
   key=(c,row['protein_id'])
   if key in seen:raise ValueError('Original species/protein key appears in multiple groups')
   seen.add(key)
   result[g][c]=row
 if len(result)!=expected or any(set(x)!=set(CODES) for x in result.values()):raise ValueError('Old bridge exact group union mismatch')
 return result

def load_species(path):
 result={}
 with open(path) as f:
  for line in f:
   key,name=line.rstrip('\n').split(': ',1);result[int(key)]=name.rsplit('.',1)[0]
 if len(result)!=54 or set(result)!=set(range(54)):raise ValueError('Species54 identity gate')
 return result

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def json_bytes(x):return (json.dumps(x,ensure_ascii=False,separators=(',',':'))+'\n').encode()
def save_new(p,x):
 with open(p,'xb') as f:f.write(json_bytes(x));f.flush();os.fsync(f.fileno())
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()

class Facts:
 def __init__(self,path):
  self.db=sqlite3.connect(path)
  self.db.executescript('''PRAGMA journal_mode=DELETE; PRAGMA cache_size=-8192;
CREATE TABLE hog(groupid TEXT,node TEXT,hog TEXT,og TEXT,parent TEXT,members TEXT,PRIMARY KEY(groupid,node,hog));
CREATE TABLE candidate(groupid TEXT,sp TEXT,gene TEXT,og TEXT,origin TEXT,PRIMARY KEY(groupid,sp,gene,og));
CREATE INDEX candidate_key ON candidate(sp,gene);
CREATE TABLE pair(groupid TEXT,ref TEXT,sp TEXT,gene TEXT,seen INT DEFAULT 0,bad INT DEFAULT 0,reverse_bad INT DEFAULT 0,PRIMARY KEY(groupid,ref,sp,gene));
CREATE TABLE evidence(groupid TEXT,ref TEXT,sp TEXT,gene TEXT,rownum INT,kind TEXT,leftn INT,rightn INT,og TEXT,PRIMARY KEY(groupid,ref,sp,gene,rownum,kind));
CREATE TABLE needed(sp TEXT,gene TEXT,PRIMARY KEY(sp,gene));
CREATE TABLE identity(sp TEXT,gene TEXT,internal TEXT,excluded INT,PRIMARY KEY(sp,gene,internal));
CREATE TABLE coordinate(sp TEXT,gene TEXT,status TEXT,record TEXT,PRIMARY KEY(sp,gene));''')
  self.db.executescript('CREATE TABLE supported(groupid TEXT,sp TEXT,gene TEXT,PRIMARY KEY(groupid,sp,gene)); CREATE INDEX supported_target ON supported(sp,gene);')
 def close(self):self.db.commit();self.db.close()
 def candidates(self,g):return list(self.db.execute('SELECT sp,gene,og,origin FROM candidate WHERE groupid=? ORDER BY sp,gene,og',(g,)))
 def hogs(self,g,node):
  return [(h,og,parent,json.loads(mem)) for h,og,parent,mem in self.db.execute('SELECT hog,og,parent,members FROM hog WHERE groupid=? AND node=?',(g,node))]
 def pair_mark(self,g,ref,sp,gene,row,kind,left,right,og):
  self.db.execute('INSERT OR IGNORE INTO pair(groupid,ref,sp,gene) VALUES(?,?,?,?)',(g,ref,sp,gene))
  if kind=='forward':self.db.execute('UPDATE pair SET seen=1,bad=max(bad,?) WHERE groupid=? AND ref=? AND sp=? AND gene=?',(int(len(left)!=1 or len(right)!=1),g,ref,sp,gene))
  else:self.db.execute('UPDATE pair SET reverse_bad=1 WHERE groupid=? AND ref=? AND sp=? AND gene=?',(g,ref,sp,gene))
  self.db.execute('INSERT OR IGNORE INTO evidence VALUES(?,?,?,?,?,?,?,?,?)',(g,ref,sp,gene,row,kind,len(left),len(right),og))
 def pair_state(self,g,ref,sp,gene):
  row=self.db.execute('SELECT seen,bad,reverse_bad FROM pair WHERE groupid=? AND ref=? AND sp=? AND gene=?',(g,ref,sp,gene)).fetchone()
  if row is None:return 'UNOBSERVED'
  count=self.db.execute('SELECT count(*) FROM pair WHERE groupid=? AND ref=? AND sp=? AND seen=1',(g,ref,sp)).fetchone()[0]
  return pair_status(bool(row[0]),bool(row[1]),count,bool(row[2]))
 def id_state(self,sp,gene):
  rows=list(self.db.execute('SELECT internal,excluded FROM identity WHERE sp=? AND gene=?',(sp,gene)))
  if any(x[1] for x in rows):return 'EXCLUDED_OG725'
  return 'EXACT_UNIQUE' if len(rows)==1 else ('SOURCE_ID_UNOBSERVED' if not rows else 'SOURCE_ID_AMBIGUOUS')

def collect_hogs(facts,paths,anchors,ref_species,tick=lambda:None,expected_species=None):
 seeds={(ref_species[c],row['protein_id']):g for g,refs in anchors.items() for c,row in refs.items()}
 counts={}
 for node,path in paths.items():
  count=0
  with text_open(path) as f:
   rd=csv.DictReader(f,delimiter='\t')
   if rd.fieldnames[:3]!=['HOG','OG','Gene Tree Parent Clade']:raise ValueError('Native HOG schema mismatch')
   if expected_species is not None and set(rd.fieldnames[3:])!=set(expected_species):raise ValueError('HOG species header differs')
   for rownum,row in enumerate(rd,1):
    if rownum%1024==0:tick()
    count+=1;members={s:tokens(row[s]) for s in rd.fieldnames[3:]}
    groups={seeds[s,gene] for s,genes in members.items() for gene in genes if (s,gene) in seeds}
    for group in groups:
     tick()
     facts.db.execute('INSERT INTO hog VALUES(?,?,?,?,?,?)',(group,node,row['HOG'],row['OG'],row['Gene Tree Parent Clade'],json.dumps(members,separators=(',',':'))))
     for sp,genes in members.items():
      for gene in genes:facts.db.execute('INSERT OR IGNORE INTO candidate VALUES(?,?,?,?,?)',(group,sp,gene,row['OG'],'HOG'))
  facts.db.commit();counts[node]=count
 return counts

def collect_pairs(facts,paths,anchors,ref_species,tick=lambda:None,expected_species=None):
 seeds={c:{refs[c]['protein_id']:g for g,refs in anchors.items()} for c in CODES};counts={}
 # All HOG-derived candidates already exist: reverse badness from ANY source
 # partner is caught, even when that source gene is not an old anchor.
 for ref,path in paths.items():
  count=0
  with text_open(path) as f:
   rd=csv.DictReader(f,delimiter='\t');source=ref_species[ref]
   if rd.fieldnames!=['Orthogroup','Species',source,'Orthologs']:raise ValueError('Aggregate header mismatch')
   for rownum,row in enumerate(rd,1):
    if rownum%1024==0:tick()
    count+=1;sp=row['Species'];left=tokens(row[source]);right=tokens(row['Orthologs']);og=row['Orthogroup']
    if expected_species is not None and sp not in expected_species:raise ValueError('Foreign aggregate species')
    for gene in left:
     if gene not in seeds[ref]:continue
     group=seeds[ref][gene]
     for target in right:
      facts.db.execute('INSERT OR IGNORE INTO candidate VALUES(?,?,?,?,?)',(group,sp,target,og,'PAIRWISE'))
      facts.pair_mark(group,ref,sp,target,rownum,'forward',left,right,og)
    for target in right:
     for (group,) in facts.db.execute("SELECT DISTINCT groupid FROM candidate WHERE sp=? AND gene=? AND origin='HOG'",(sp,target)).fetchall():
      expected=anchors[group][ref]['protein_id']
      if left!=[expected] or len(right)!=1:facts.pair_mark(group,ref,sp,target,rownum,'reverse_other_partner',left,right,og)
  facts.db.commit();counts[ref]=count
 return counts

def collect_identities(facts,path,species,anchors,ref_species,excluded,tick=lambda:None):
 facts.db.execute('INSERT OR IGNORE INTO needed SELECT DISTINCT sp,gene FROM candidate')
 for refs in anchors.values():
  for c,row in refs.items():facts.db.execute('INSERT OR IGNORE INTO needed VALUES(?,?)',(ref_species[c],row['protein_id']))
 count=0
 with open(path) as f:
  for line in f:
   if count%1024==0:tick()
   count+=1;internal,description=line.rstrip('\n').split(': ',1);spid=int(internal.split('_',1)[0]);protein=description.split()[0];sp=species[spid]
   if facts.db.execute('SELECT 1 FROM needed WHERE sp=? AND gene=?',(sp,protein)).fetchone():facts.db.execute('INSERT INTO identity VALUES(?,?,?,?)',(sp,protein,internal,int(internal in excluded)))
 facts.db.commit();return count

def collect_coordinates(facts,path,receipt,tick=lambda:None):
 if path is None:return 0
 if not receipt:raise ValueError('Coordinate ledger needs accepted source receipt')
 count=0
 with text_open(path) as f:
  rd=csv.DictReader(f,delimiter='\t')
  for row in rd:
   if count%1024==0:tick()
   count+=1;status=row['coordinate_status']
   if status not in COORD_STATES:raise ValueError('Unknown coordinate state')
   if status=='ACCEPTED_COORDINATE':
    for k in ['seqid','gene_id','transcript_id','strand','gene_start0','gene_end0','window_center0']:
     if not row.get(k):raise ValueError('Incomplete accepted coordinate evidence')
    start,end=int(row['gene_start0']),int(row['gene_end0'])
    if start<0 or end<=start or row['strand'] not in ('+','-'):raise ValueError('Coordinate schema conflict')
   facts.db.execute('INSERT INTO coordinate VALUES(?,?,?,?)',(row['species_name'],row['protein_id'],status,json.dumps({**row,'accepted_source_receipt':receipt},separators=(',',':'))))
 facts.db.commit();return count

def shared_hog(facts,g,node,required):
 rows=[r for r in facts.hogs(g,node) if all(protein in r[3].get(sp,[]) for sp,protein in required)]
 return rows

def tree_tip_map(tree,species_names):
 prefixes=sorted([(s.replace('.','_').replace(' ','_')+'_',s) for s in species_names],key=lambda x:-len(x[0]));result={}
 for leaf in tree.leaves():
  matches=[(p,s) for p,s in prefixes if leaf.name.startswith(p)]
  if len(matches)!=1:raise ValueError('Public tree species prefix mismatch')
  result[leaf.name]=matches[0][1]
 return result

def gene_label(sp,protein):return sp.replace('.','_').replace(' ','_')+'_'+protein

def read_gene_trees(og,new,maxbytes):
 result=[]
 for directory in ['Gene_Trees','Resolved_Gene_Trees']:
  p=Path(new)/directory/(og+'_tree.txt')
  if not p.is_file():return None
  if p.stat().st_size>maxbytes:raise ValueError('Released single-tree byte bound exceeded')
  result.append(parse_newick(p.read_text()))
 return result

def evaluate_group(facts,g,refs,ref_species,species_tree,new,max_tree_bytes):
 original=[(ref_species[c],refs[c]['protein_id']) for c in CODES];g0=[]
 for sp,gene in original:
  state=facts.id_state(sp,gene)
  if state!='EXACT_UNIQUE':g0.append(state)
 base=shared_hog(facts,g,'N0',original);old=shared_hog(facts,g,'N26',original)
 if len(base)!=1 or len(old)!=1:g0.append('G0_HOG_CONFLICT')
 if len(base)==1 and len(old)==1 and base[0][1]!=old[0][1]:g0.append('G0_HOG_OG_CONFLICT')
 if any(og=='OG0000725' for _,og,_,_ in base+old):g0.append('EXCLUDED_OG725')
 for c in CODES:
  for other in CODES:
   if c!=other and facts.pair_state(g,c,ref_species[other],refs[other]['protein_id'])!='ONE_TO_ONE':g0.append('G0_PAIRWISE_CONFLICT')
 # Distinct original labels sharing a root HOG are retained as competing groups.
 converging=[]
 for hog,_,_,_ in base:
  converging += [x[0] for x in facts.db.execute('SELECT DISTINCT groupid FROM hog WHERE node=? AND hog=? AND groupid!=?',('N0',hog,g))]
 # Shared root HOG is a diagnostic, not an added whole-group exclusion.
 # True partner ambiguity is already enforced by complete forward/reverse facts.
 by_og={};out=[]
 for sp,gene,og,origin in facts.candidates(g):
  req=original+[(sp,gene)];localnode=mrca(species_tree,[x[0] for x in req]).name
  n0=shared_hog(facts,g,'N0',req);local=shared_hog(facts,g,localnode,req)
  hog_same=(len(base)==1 and len(n0)==1 and n0[0][0]==base[0][0] and len(local)==1 and n0[0][1]==local[0][1]==og and og!='OG0000725')
  multi=any(len(row[3].get(sp,[]))>1 for row in local)
  pairs=[facts.pair_state(g,c,sp,gene) for c in CODES if ref_species[c]!=sp]
  if og not in by_og:by_og[og]=read_gene_trees(og,new,max_tree_bytes)
  trees=by_og[og];tree_state='UNOBSERVED'
  if trees:
   labels=list(dict.fromkeys(gene_label(s,p) for s,p in req))
   if all(set(labels)<={x.name for x in t.leaves()} for t in trees):tree_state=tree_consistency(trees[1],species_tree,labels,tree_tip_map(trees[1],[x.name for x in species_tree.leaves()]))
  co=facts.db.execute('SELECT status,record FROM coordinate WHERE sp=? AND gene=?',(sp,gene)).fetchone()
  coordinate=co[0] if co else 'COORDINATE_EVIDENCE_NOT_SUPPLIED'
  if any(sp==s and gene==p for s,p in original):coordinate='ACCEPTED_COORDINATE'
  f={'hog_same':hog_same,'hog_multicopy':multi,'pairwise':pairs,'tree':tree_state,'source_identity':facts.id_state(sp,gene),'coordinate':coordinate}
  status,reasons=classify_candidate(f)
  if not pairs and sp not in ref_species.values():status='UNQUALIFIED';reasons.append('PAIRWISE_UNOBSERVED')
  evidence=[dict(zip(['reference','row','kind','source_cardinality','target_cardinality','OG'],x)) for x in facts.db.execute('SELECT ref,rownum,kind,leftn,rightn,og FROM evidence WHERE groupid=? AND sp=? AND gene=? ORDER BY ref,rownum,kind',(g,sp,gene))]
  overlaps=[x[0] for x in facts.db.execute('SELECT DISTINCT groupid FROM candidate WHERE sp=? AND gene=? AND groupid!=? ORDER BY groupid',(sp,gene,g))]
  out.append({'original_group':g,'species_name':sp,'protein_id':gene,'new_OG':og,'origin':origin,'local_species_MRCA':localnode,'N0_HOGs':[x[0] for x in n0],'local_HOGs':[x[0] for x in local],'status':status,'reasons':sorted(set(reasons)),'G0_reasons':sorted(set(g0)),'member_supported':not g0 and status!='UNQUALIFIED','collision_overlay_required':True,'measurement_or_model_eligible':False,'shared_root_HOG_old_labels':sorted(set(converging)),'other_original_groups_same_candidate':overlaps,'competition_is_diagnostic_only':True,'facts':f,'pairwise_evidence':evidence,'coordinate_record':json.loads(co[1]) if co else None,'inference_scope':'accepted-native-primary-conditional','repeat_effect_used':False})
 # Cache is per original group only. Shared OGs across groups are flagged above;
 # this reads immutable trees, never reinfers or runs native reconciliation.
 return sorted(set(g0)),sorted(set(converging)),out

def export_supported_qualification(facts,out,tick=lambda:None):
 """Mandatory collision overlay; only actual supported targets are competing."""
 counts={'supported_keys':0,'collision_targets':0,'ambiguous_supported_keys':0}
 with gzip.open(Path(out)/'supported-collisions.jsonl.gz','wt',encoding='utf-8') as col:
  for sp,gene in facts.db.execute('SELECT sp,gene FROM supported GROUP BY sp,gene HAVING count(*)>1 ORDER BY sp,gene'):
   tick();groups=[x[0] for x in facts.db.execute('SELECT groupid FROM supported WHERE sp=? AND gene=? ORDER BY groupid',(sp,gene))]
   col.write(json_bytes({'species_name':sp,'protein_id':gene,'original_groups':groups,'status':'AMBIGUOUS_SUPPORTED_TARGET','all_original_labels_retained':True}).decode());counts['collision_targets']+=1
 with gzip.open(Path(out)/'supported-qualification.jsonl.gz','wt',encoding='utf-8') as dst:
  for group,sp,gene in facts.db.execute('SELECT groupid,sp,gene FROM supported ORDER BY groupid,sp,gene'):
   tick();groups=[x[0] for x in facts.db.execute('SELECT groupid FROM supported WHERE sp=? AND gene=? ORDER BY groupid',(sp,gene))];ambiguous=len(groups)>1
   counts['supported_keys']+=1;counts['ambiguous_supported_keys']+=int(ambiguous)
   dst.write(json_bytes({'original_group':group,'species_name':sp,'protein_id':gene,'status':'AMBIGUOUS_SUPPORTED_TARGET' if ambiguous else 'SUPPORTED_UNIQUE_ORIGINAL_LABEL','supported_original_labels':groups,'post_collision_member_eligible':not ambiguous,'measurement_or_model_eligible':False,'coordinate_and_G3_qualification_still_required':True}).decode())
 return counts

def resource_guard(release,out,start):
 """Root also supervises the whole process; this checks streamed metadata steps."""
 import psutil
 proc=psutil.Process();lim=release['limits'];last=[0.0]
 guardian=release['external_guard_identity'];parent=psutil.Process(proc.ppid())
 if parent.pid!=guardian['pid'] or parent.create_time()!=guardian['create_time']:raise ValueError('Exact Root resource guardian identity differs')
 def tick():
  now_m=time.monotonic()
  if now_m-last[0]<0.1:return
  last[0]=now_m
  if now_m-start>lim['wall_seconds']:raise RuntimeError('Wall bound exceeded')
  if proc.memory_info().rss>lim['process_rss_bytes']:raise RuntimeError('RSS bound exceeded')
  if psutil.virtual_memory().available<lim['host_reserve_bytes']:raise RuntimeError('Host reserve bound exceeded')
  if psutil.disk_usage(out).free<lim['disk_floor_bytes']:raise RuntimeError('Disk floor exceeded')
  if sum(p.stat().st_size for p in Path(out).iterdir() if p.is_file())>lim['new_artifact_bytes']:raise RuntimeError('New artifact bound exceeded')
 return tick

def run_join(config,release,outdir):
 """Production entry. Root prepares source/resource/once release; not called in prep."""
 require_decisions(config['root_decisions'])
 if sha(config['Root_disposition_path'])!=config['Root_disposition_sha256']:raise ValueError('Root preoutcome disposition binding differs')
 if release.get('status')!='ROOT_EXACT_ONCE_G2_MEMBER_BRIDGE' or release.get('approved_by')!='current-conversation-root':raise ValueError('Root exact release required')
 if datetime.datetime.now(datetime.timezone.utc)>=datetime.datetime.fromisoformat(release['first_action_valid_until']):raise ValueError('Expired fresh gate')
 if release.get('config_sha256')!=hashlib.sha256(json_bytes(config)).hexdigest():raise ValueError('Config binding mismatch')
 # Code/small contracts pinned here. Huge accepted bodies must be bound by Root
 # receipts, not silently rehashed by this adapter.
 for p,h in release['small_file_pins'].items():
  if sha(p)!=h:raise ValueError('Changed code/contract pin')
 if not release.get('fresh_resources_accepted') or not release.get('accepted_large_source_identities_bound'):raise ValueError('Fresh Root resource/source gates absent')
 out=Path(outdir);out.mkdir(exist_ok=False);save_new(out/'once-intent.json',{'at_utc':now(),'source':release,'automatic_retry':False})
 csv.field_size_limit(8*1024**2);start=time.monotonic();facts=None
 try:
  tick=resource_guard(release,out,start);tick()
  species=load_species(config['species_ids']);tree=parse_newick(Path(config['rooted_node_tree']).read_text())
  if {x.name for x in tree.leaves()}!=set(species.values()):raise ValueError('Tree54 names differ')
  refs={c:species[int(i)] for c,i in config['original_species_ids'].items()}
  if set(refs)!=set(CODES) or mrca(tree,list(refs.values())).name!='N26':raise ValueError('Old member node differs')
  nodes={mrca(tree,list(refs.values())+[sp]).name for sp in species.values()}|{'N0','N26'}
  if nodes!=set(config['hog_paths']):raise ValueError('Required ancestor HOG nodes not exact')
  anchors=load_old(config['old_bridge']);facts=Facts(out/'scratch.sqlite')
  counts={'HOG_rows_read':collect_hogs(facts,config['hog_paths'],anchors,refs,tick,set(species.values())),'aggregate_rows_read':collect_pairs(facts,config['aggregate_paths'],anchors,refs,tick,set(species.values()))}
  excluded=set(json.loads(Path(config['excluded_members']).read_text())['genes'])
  counts['SequenceIDs_rows_read']=collect_identities(facts,config['sequence_ids'],species,anchors,refs,excluded,tick)
  counts['coordinate_rows_read']=collect_coordinates(facts,config.get('coordinate_ledger'),config.get('coordinate_acceptance_receipt'),tick)
  with gzip.open(out/'status-matrix.tsv.gz','wt',newline='',encoding='utf-8') as matrix,gzip.open(out/'cell-states.jsonl.gz','wt',encoding='utf-8') as cells,gzip.open(out/'candidates.jsonl.gz','wt',encoding='utf-8') as candidates:
   writer=csv.writer(matrix,delimiter='\t');writer.writerow(['original_group']+[species[i] for i in range(54)]);cellcount=0;candidatecount=0
   for g in sorted(anchors):
    tick()
    g0,converging,rows=evaluate_group(facts,g,anchors[g],refs,tree,config['new_results'],release['limits']['max_single_tree_bytes']);by_sp=collections.defaultdict(list)
    for row in rows:
     candidates.write(json_bytes(row).decode());by_sp[row['species_name']].append(row);candidatecount+=1
     if row['member_supported']:facts.db.execute('INSERT OR IGNORE INTO supported VALUES(?,?,?)',(g,row['species_name'],row['protein_id']))
    wide=[g]
    for i in range(54):
     sp=species[i];items=by_sp[sp];status,_=cell_status([x['status'] for x in items],bool(g0));original_code=next((c for c in CODES if refs[c]==sp),None)
     if original_code:status='ORIGINAL_MEMBER_RETAINED'
     cell={'original_group':g,'species_id':i,'species_name':sp,'status':status,'G0_reasons':g0,'shared_root_HOG_old_labels':converging,'candidate_overlap_original_groups':sorted({z for x in items for z in x['other_original_groups_same_candidate']}),'candidate_count':len(items),'candidate_proteins':[x['protein_id'] for x in items],'candidate_reasons':sorted({z for x in items for z in x['reasons']}),'collision_overlay_required':True,'primary_conditional':True,'measurement_or_model_eligible':False,'original_coordinate':anchors[g][original_code] if original_code else None}
     cells.write(json_bytes(cell).decode());wide.append(status);cellcount+=1
    writer.writerow(wide)
    if time.monotonic()-start>release['limits']['wall_seconds']:raise RuntimeError('Released wall bound exceeded')
  if cellcount!=1896*54:raise ValueError('Matrix exact accounting mismatch')
  counts['supported_qualification']=export_supported_qualification(facts,out,tick)
  facts.close();facts=None
  save_new(out/'terminal.json',{'at_utc':now(),'status':'G2_MEMBER_BRIDGE_JOIN_RETURNED_NOT_SCIENTIFIC_ACCEPTANCE','groups':1896,'cells':cellcount,'candidate_rows':candidatecount,'counts':counts,'automatic_retry':False,'source_science_repeated':0,'new_repeat_outcomes':0,'Root_output_acceptance_required':True})
 except BaseException as exc:
  if facts:facts.close()
  if not (out/'HOLD.json').exists():save_new(out/'HOLD.json',{'at_utc':now(),'error':type(exc).__name__+': '+str(exc),'automatic_retry':False})
  raise

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--config',required=True);p.add_argument('--release',required=True);p.add_argument('--out',required=True);a=p.parse_args()
 run_join(json.loads(Path(a.config).read_text()),json.loads(Path(a.release).read_text()),a.out)
if __name__=='__main__':main()
