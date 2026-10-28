"""Root new downstream row projection; accepted source reads only, no QA or model."""
from pathlib import Path
import contextlib,collections,gzip,hashlib,importlib.util,json,os,sqlite3,sys,time

def need(x,msg):
 if not x:raise ValueError('HOLD:'+msg)
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def encoded(x):return json.dumps(x,ensure_ascii=False,sort_keys=True,separators=(',',':')).encode()
def snapshot(p):
 s=Path(p).stat();return dict(dev=s.st_dev,ino=s.st_ino,mode=s.st_mode,size=s.st_size,mtime_ns=s.st_mtime_ns,ctime_ns=s.st_ctime_ns)
def build(cfg,output,check,event):
 start=time.monotonic();out=Path(output);out.mkdir();sources={x['name']:x for x in cfg['sources']};receipts={};code=cfg['codes']
 srcmod=module('root_source_recording',code['source']);Source=srcmod.Source
 def source_json(name):
  with Source(sources[name],check,receipts,event) as src:
   raw=src.bytes_all();return srcmod.decode(raw)
 def hash_database(path,cap):
  need(str(path)==sources['Genome_DB']['path'] and sources['Genome_DB']['bytes']<=cap,'only exact accepted GenomeDB and originalcap')
  with Source(sources['Genome_DB'],check,receipts,event) as src:
   while True:
    check();x=src.f.read(65536);src.n+=len(x);need(src.n<=cap,'original derivedDB cap');src.h.update(x);src.progress()
    if not x:break
   src.eof=True
   return {'bytes':src.n,'sha256':src.h.hexdigest()}
 # All metadata below was independently pinned by Root; no scientific prehash.
 cache_proof=cfg['cache_acceptance'];cp=cfg['cache_database'];need(snapshot(cp['path'])==cfg['cache_closed_stat'],'accepted globalcache unchanged')
 connection=sqlite3.connect(Path(cp['path']).resolve().as_uri()+'?mode=ro&immutable=1',uri=True);connection.row_factory=sqlite3.Row;connection.execute('PRAGMA query_only=ON');connection.execute('PRAGMA cache_size=-4096');connection.execute('PRAGMA mmap_size=0')
 native=cfg['catalogue_record']['native_species_name'];assembly=cfg['catalogue_record']['assembly'];rec=cfg['catalogue_record'];n=rec['targets'];supported=rec['source_supported'];na=rec['typed_NA'];refs=cfg['refs'];need(n==supported+na,'full source target accounting')
 def rows(q,args=()):
  for row in connection.execute(q,args):check();yield dict(row)
 cache={'target_rows':list(rows('SELECT native,protein,assembly,ord FROM targets WHERE native=? ORDER BY ord',(native,))), 'signature_rows':list(rows('SELECT g,xa,xk,source FROM signatures ORDER BY g')), 'state_rows':list(rows('SELECT g,native,status FROM cells WHERE native=? ORDER BY g',(native,))), 'member_rows':list(rows('SELECT native,protein,g,body,candidate_ok,qualification_ok FROM member WHERE native=? ORDER BY protein',(native,))), 'collision_rows':[]}
 need(cache_proof['actual_supported_collision_targets']==0 and cache_proof['complete_source_streams']==7,'not guessed empty collision view')
 ids=[x['protein'] for x in cache['target_rows']];need(len(ids)==n,'full cache targets forspecies')
 context={'schema':'ROOT_ACCEPTED_MEMBER_SIGNATURE_BINDING_INPUT_VIEW_V1','approved_by':'current-conversation-root','complete_input_source_pins_EOF_verified':True,'full_signature_groups':1896,'signature_ref':refs['signature_ref'],'refs':{'member':refs['member_acceptance_ref'],'evidence':refs['member_evidence_ref'],'collision':refs['collision_ref'],'canonical':refs['canonical_ref'],'report':refs['report_ref'],'NA_authority':refs['source_flag_authority_ref']},'limits':{'max_targets':2000,'max_evidence_rows':100000,'max_signature_rows':2000,'max_collision_rows':100000,'max_stored_evidence_bytes':16777216,'max_matches_per_target':8},'expected_rows':{'canonical':n,'evidence':n,'signatures':1896,'collisions':0},'frozen_species':{'assembly':assembly,'native_species_name':native,'target_ids':ids,'target_count':n,'target_set_sha256':hashlib.sha256(('\n'.join(sorted(ids))+'\n').encode()).hexdigest()}}
 canonical=[]
 with Source(sources['canonical'],check,receipts,event) as source:
  for line in source.lines():check();canonical.append(srcmod.decode(line));need(len(canonical)<=n,'canonical cap')
 need(len(canonical)==n and collections.Counter(x['coordinate_status'] for x in canonical)['SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY']==supported,'canonicalfulltarget/sourceNAcounts')
 report=source_json('report');need(report['schema']=='ncbi-sequence-report-state-v1' and report['assembly']==assembly and report['complete'] is True and report['next_page_token'] is None and len(report['records'])==report['total_count']==rec['report_metadata_declared_owners'],'full acceptedreport exactassembly/EOF')
 owners=report['records'];token_owners={}
 for i,r in enumerate(owners):
  need(r['assembly_accession']==assembly,'exact report ownerassembly')
  for k in ('refseq_accession','genbank_accession','sequence_name'):
   if r.get(k):token_owners.setdefault(r[k].encode('ascii'),set()).add(i)
 old_assemblies=set(cfg['old3_assemblies']);need(assembly not in old_assemblies and cfg['role_record']['fixed_study_role']!='old_discovery','new51 distinctnamespace only; old3neverevaluate')
 policy_ref=refs['source_flag_authority_ref'];canonical_policy={'schema':'ROOT_LITERAL_CANONICAL_VIEW_POLICY_V1','approved_by':'current-conversation-root','marker':'COPY_ACCEPTED_CANONICAL_STATUS_NO_RECLASSIFICATION','ref':policy_ref}
 role={'study_role':cfg['role_record']['fixed_study_role'],'major_layer':cfg['role_record']['fixed_major_layer'],'native_species_name':native,'assembly':assembly}
 builder=module('root_binding',code['builder']);cv=module('root_literal_canonical',code['canonical']);fa=module('root_fact_adapter',code['facts']);driver=module('root_existing_row_driver',code['driver']);bridge=module('root_existing_row_bridge',code['bridge']);integration=module('root_row_integration',code['integration'])
 sys.path.insert(0,str(Path(code['protocol']).parent));proto=module('genome_protocol',code['protocol']);sys.modules['genome_protocol']=proto
 factory={'schema':'ROOT_SINGLE_SPECIES_PREOUTCOME_ROW_FACTORY_V1','approved_by':'current-conversation-root',**refs,'species_role':role,'facts_mode':'ACCEPTED_NA_ONLY_NO_DEFINED_WINDOWS' if supported==0 else 'ACCEPTED_GENOME_FACTS'}
 fact_p={'schema':'ROOT_FACT_VIEW_PROVENANCE_V1','kind':rec['GENOME']['kind'],'actual_projection_released_by_this_adapter':False,'catalogue_ref':cfg['catalogue_ref'],'catalogue_acceptance_ref':cfg['catalogue_acceptance_ref'],'canonical_acceptance_ref':refs['canonical_acceptance_ref'],'source_disposition_ref':refs['all_QC_dispositions_closed_ref'],'code_freeze_ref':cfg['code_freeze_ref']}
 report_bindings=[];protection=[];source_eof=False
 def protection_row(t,radius):
  return {'target_key':t,'radius':radius,'status':'LEGACY_COORDINATE_EXCLUSION_NOT_APPLICABLE_DISTINCT_ASSEMBLY','legacy_overlap_status':'NOT_COMPARABLE_OUTSIDE_FROZEN_COORDINATE_DOMAIN','new_species_Hox_exclusion_status':'UNEVALUATED_NO_COMPLETE_NEW_HOX_MAP_CLAIM','policy_ref':policy_ref,'assembly':assembly}
 ordered=[]
 for x in canonical:
  target=builder.key(assembly,x['key']['protein_id']);ordered.append({'target_key':target,'coordinate_status':x['coordinate_status'],'NA_reason':x['NA_reason']})
  for radius in ((50000,25000) if x['coordinate_status']=='SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY' else (None,)):protection.append(protection_row(target,radius))
 root={'binding_context':context,'canonical_records':canonical,'canonical_source_authority_refs':[refs['canonical_acceptance_ref'],policy_ref],'canonical_policy':canonical_policy,'canonical_limits':cv.CAPS,'cache_acceptance_ref':cfg['cache_acceptance_ref'],'root_policy_ref':policy_ref,'fact_kind':rec['GENOME']['kind'],'fact_provenance':fact_p,'expected_output_counts':{'targets':n,'windows':supported*2,'typed_NA':na},'report_bindings':report_bindings,'protection_records':protection}
 admitted={(x['path'],x['sha256']) for x in refs.values() if isinstance(x,dict) and 'sha256' in x}
 for x in [cfg['cache_acceptance_ref'],policy_ref,cfg['catalogue_ref'],cfg['catalogue_acceptance_ref']]:admitted.add((x['path'],x['sha256']))
 class Hooks:
  def verify_accepted_control(self,r,role):need((r['path'],r['sha256']) in admitted,'Root admitted exactreference '+role)
  def verify_exact_release_and_source_roles(self,f):need(cfg['Root_actual_projection_released'] is True,'actualsourcegate');check()
  def check_same_resource_deadline(self):check()
  def verify_final_counts_source_EOF_hashes(self,f,**counts):
   need(counts==root['expected_output_counts'] and set(receipts)==set(sources),'allnewinputEOF source receipts/counts');check()
 def attach_report(f):
  key=f['target'];radius=f['radius']
  if f['geometry_status']=='SEQID_UNOBSERVED':return
  cand=token_owners.get(f['actual_token_raw'],set()) & token_owners.get(f['canonical_seqid_raw'],set());need(len(cand)==1,'sameexactuniqueactualtoken and canonicaltoken reportowner')
  idx=next(iter(cand));r=owners[idx];need(r['length']==f['sequence_length'],'actual accepted header reportlength')
  required=('role','assigned_molecule_location_type','assembly_unit','chr_name')
  if not all(k in r for k in required):report_bindings.append({'target_key':key,'radius':radius,'owner':None,'location':None,'typed_pending_reason':'ACCEPTED_REPORT_OWNER_FIELDS_MISSING','raw_report_owner':r});return
  real=r['role']=='assembled-molecule' and r['assigned_molecule_location_type']=='Chromosome' and r['chr_name'] not in ('Un','',None)
  loc={'schema':'ROOT_ACCEPTED_REPORT_LOCATION_DISPOSITION_V1','approved_by':'current-conversation-root','owner_index':idx,'report_ref':refs['report_ref'],'source_fields':{k:r[k] for k in required},'accepted_assembly_unit':r['assembly_unit'],'kind':'REAL_CHROMOSOME' if real else 'LOCAL_DESCRIPTIVE','location':r['chr_name'] if real else 'REPORT_OWNER:'+str(idx)}
  need(type(r['assembly_unit']) is str,'explicitliteralassemblyunit')
  report_bindings.append({'target_key':key,'radius':radius,'owner':{'owner_index':idx,'record':r},'location':loc})
 @contextlib.contextmanager
 def fact_context():
  if rec['GENOME']['kind']=='TYPED_NA_ONLY_NO_FACTS_OR_FAKE_QA':
   need(supported==0 and rec['GENOME']['QA_ref'] is None,'honestNAonly');root['fact_inputs']={'original_targets':ordered};factory['NA_only_disposition_ref']=rec['GENOME']['Root_acceptance_ref'];yield;return
  need(rec['GENOME']['kind']=='ACCEPTED_READONLY_SQLITE','THIS_FIRST_IO_RELEASE_SUPPORTS_SQLITE_AND_HONEST_NA_ONLY; no defaultformatconversion')
  qa=source_json('Genome_QA');ri=module('root_original_Genome_readonly',code['read_interface']);ri.file_hash=hash_database
  fact_p.update(Genome_QA_ref=rec['GENOME']['QA_ref'],Genome_acceptance_ref=rec['GENOME']['Root_acceptance_ref'],scientific_artifact_ref=rec['GENOME']['scientific_artifact']['accepted_descriptor']);context['refs']['genome_QA']=rec['GENOME']['QA_ref']
  factory['genome_QA_ref']=rec['GENOME']['QA_ref'];factory['genome_acceptance_ref']=rec['GENOME']['Root_acceptance_ref']
  with ri.open_view(sources['Genome_DB']['path'],qa) as original:
   facts=[];nas=[]
   for f in original.facts():check();facts.append(f);need(len(facts)<=2*supported,'factsbound')
   for f in original.typed_na():check();nas.append(f);need(len(nas)<=na,'NAbound')
   need(len(facts)==2*supported and len(nas)==na,'allfactandNAphysicalcursorsEOF')
   for f in facts:attach_report(f)
   class Buffered:
    def __init__(self):self.used=set()
    def facts(self):need('facts' not in self.used,'noDBfactreplay');self.used.add('facts');yield from facts
    def typed_na(self):need('na' not in self.used,'noDBNAreplay');self.used.add('na');yield from nas
   root['fact_inputs']={'accepted_view':Buffered(),'ordered_targets':ordered};yield
 parts={'driver':driver.prepared_rows,'builder':builder.build_bindings,'map_fact':bridge.map_fact,'canonical_view':cv.prepare_canonical_view,'fact_adapter':fa,'target_key':builder.key}
 counts=collections.Counter();mapped=collections.Counter();reasons=collections.Counter();width_counts=collections.Counter();outpath=out/'row-events.jsonl.gz';hh=hashlib.sha256();bytes_raw=0;count=0
 try:
  with fact_context():
   with outpath.open('xb') as fd:
    with gzip.GzipFile(fileobj=fd,mode='wb',mtime=0,compresslevel=6) as gz:
     for e in integration.iter_species_events(factory,cache,root,parts,Hooks()):
      check();raw=integration.encode_event(e,proto.pack)+b'\n';need(len(raw)<=4*2**20,'eventcap');hh.update(raw);bytes_raw+=len(raw);count+=1;counts[e['kind']]+=1;gz.write(raw)
      if e['kind']=='TARGET_ROW_DRAFT':
       m=e['draft']['row_preparation'];mapped[m['status']]+=1;reasons.update(m['reasons']);width_counts[str(e['draft']['radius'])]+=1
    fd.flush();os.fsync(fd.fileno())
  need(snapshot(cp['path'])==cfg['cache_closed_stat'],'globalcache statclosure');need(counts['FULL_GROUP_SOURCE_STATE']==1896 and counts['TARGET_ROW_DRAFT']==2*supported+na and counts['SPECIES_DRAFT_COMPLETE']==1,'fullspeciescountclosure')
  hout=hashlib.sha256()
  with outpath.open('rb') as f:
   for chunk in iter(lambda:f.read(1048576),b''):check();hout.update(chunk)
  ledger=out/'report-and-protection-bindings.json.gz'
  with ledger.open('xb') as fd:
   with gzip.GzipFile(fileobj=fd,mode='wb',mtime=0) as gz:gz.write(proto.pack({'assembly':assembly,'report_ref':refs['report_ref'],'root_policy_ref':policy_ref,'report_bindings':report_bindings,'protection_records':protection})+b'\n')
   fd.flush();os.fsync(fd.fileno())
  terminal={'schema':'ROOT_REAL_SINGLE_SPECIES_PREOUTCOME_ROWS_V1','status':'FULL_SPECIES_ROW_PROJECTION_COMPLETE_PENDING_ROOT_ACCEPTANCE','packet':rec['packet'],'assembly':assembly,'native_species_name':native,'source_receipts':receipts,'source_counts':root['expected_output_counts'],'row_event_counts':dict(counts),'mapped_status_counts':dict(mapped),'typed_NA_reason_counts':dict(reasons),'window_counts':dict(width_counts),'output':{'path':str(outpath),'bytes':outpath.stat().st_size,'sha256':hout.hexdigest(),'uncompressed_bytes':bytes_raw,'uncompressed_sha256':hh.hexdigest(),'records':count},'role':role,'J_weights_S_or_models_assigned':False,'old_results_recomputed':False,'seconds':time.monotonic()-start}
  with (out/'terminal.json').open('xb') as f:f.write(encoded(terminal)+b'\n');f.flush();os.fsync(f.fileno())
  return terminal
 finally:connection.close()
