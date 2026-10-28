"""First permitted species finite fits from accepted observations, no new RM IO."""
from pathlib import Path
import hashlib,importlib.util,json,os,time

def module(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);import sys;sys.modules[n]=m;s.loader.exec_module(m);return m
def need(x,m):
 if not x:raise ValueError('HOLD:'+m)
def syncdir(p):
 fd=os.open(p,os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def new(p,raw):
 p=Path(p);q=p.with_name(p.name+'.publish-'+str(os.getpid()))
 with q.open('xb') as f:f.write(raw+b'\n');f.flush();os.fsync(f.fileno())
 os.link(q,p,follow_symlinks=False);syncdir(p.parent);q.unlink();syncdir(p.parent)
def ref(p):p=Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def project_saved_observation(x,physical,rawline_sha256,source_ref,primary_ids,local_ids,seen,obs,lineage,local_refs):
 radius=x['radius'];g=x['old_group'];key=g,radius
 need(radius in (50000,25000) and key not in seen,'unique exact saved target/radius')
 seen.add(key)
 reference=dict(old_group=g,radius=radius,physical_record_index=physical,rawline_sha256=rawline_sha256,source_ref=source_ref,status=x['status'])
 if g not in primary_ids:
  need(g in local_ids and x['frozen_J_member'] is False,'only frozen localdescriptive rows outside mainJ, never effectchosen')
  local_refs.append(reference);return
 need(x['frozen_J_member'] is True,'literal mainJ source flag')
 if x['status']=='MEASURED_DRAFT':
  c=x['coverage'];need(type(c['A']) is int and type(c['U']) is int and 0<=c['U']<=c['A']<=c['D'],'accepted integer counts')
  obs[radius][g]=dict(status='MEASURED',A=c['A'],U=c['U'])
 else:
  need(x['status']=='TYPED_NA' and x['NA_reasons'],'honest accepted sourceNA')
  obs[radius][g]=dict(status='TYPED_NA',reason=';'.join(x['NA_reasons']),reasons=x['NA_reasons'])
 lineage[radius][g]=reference

def build(cfg,output,check,event):
 start=time.monotonic();out=Path(output);out.mkdir();runtime=out.parent;(runtime/'finite-call-claims').mkdir();syncdir(runtime)
 source=module('firstspecies_Source',cfg['codes']['source']);io=module('firstspecies_JSON',cfg['codes']['source_json']);caller=module('accepted_finite_caller',cfg['codes']['caller']);backend=module('fixed_original_backend',cfg['codes']['backend']);sources={x['name']:x for x in cfg['sources']};receipts={}
 with source.Source(sources['panel'],check,receipts,event) as src:panel=io.json_from_source(src,sources['panel'],source.decode)
 obs={50000:{},25000:{}};lineage={50000:{},25000:{}};physical=0;local_refs=[];primary_ids=set(panel['panel']['ids']);local_ids=set(panel['panel']['local_descriptive_ids']);seen_measurements=set()
 with source.Source(sources['measurements'],check,receipts,event) as src:
  for raw in src.lines():
   check();physical+=1;x=source.decode(raw);project_saved_observation(x,physical,hashlib.sha256(raw).hexdigest(),cfg['source_measurements'],primary_ids,local_ids,seen_measurements,obs,lineage,local_refs)
 need(physical==cfg['source_rows'] and len(local_refs)==2*len(local_ids) and len(seen_measurements)==2*(len(primary_ids)+len(local_ids)) and all(set(v)==primary_ids for v in obs.values()) and set(receipts)=={'panel','measurements'} and all(r['full_EOF'] for r in receipts.values()),'complete source files/counts before firstmodel')
 panelref={k:sources['panel'][k] for k in ('path','bytes','sha256')};base=dict(schema='ROOT_ACCEPTED_FINITE_MODEL_INPUT_BINDING_V1',approved_by='current-conversation-root',complete_source_hash_EOF_verified=True,backend_sha256=cfg['backend_sha256'],numerical_supplement_sha256=cfg['numerical_supplement_sha256'],preoutcome_acceptance_ref=cfg['preoutcome_gate_ref'],panel_ref=panelref,panel_acceptance_ref=cfg['J_acceptance_ref'],native_species_name=panel['native_species_name'],assembly=panel['assembly'],expected_main_J_count=cfg['expected_main_J_count'])
 species=caller.restore_species_panel(panel,base);observations={}
 for radius in (50000,25000):
  binding=dict(schema='ROOT_ACCEPTED_FINITE_OBSERVATION_BINDING_V1',approved_by='current-conversation-root',complete_source_hash_EOF_verified=True,source_ref=cfg['source_measurements'],acceptance_ref=cfg['measurement_acceptance_ref'],radius=radius,native_species_name=species.native_species_name,assembly=species.assembly,panel_sha256=panelref['sha256']);observations[radius]=caller.restore_observations(obs[radius],species,radius,binding)
 new(out/'observation-record-lineage.json',caller.encode_result(dict(source_file_ref=cfg['source_measurements'],physical_rows=physical,radius_maps={str(k):v for k,v in lineage.items()},NA_records_not_zero=True,frozen_local_descriptive_records_not_primary_coefficients=local_refs)))
 global_claims=Path(cfg['global_numeric_claims']);global_claims.mkdir(exist_ok=True);syncdir(global_claims.parent)
 def claim(operation,inputs):
  check();need(operation.startswith('SPATIAL:'+panelref['sha256']+':') and operation.rsplit(':',1)[1] in cfg['allowed_versions'],'exact Root allowed3speciesversions only');key=hashlib.sha256(operation.encode()).hexdigest();obj=dict(schema='ROOT_ONCE_FINITE_CALL_CLAIM_V1',approved_by='current-conversation-root',operation_id=operation,status='CLAIMED_ONCE',inputs_immutable_bound=True,inputs=inputs,first_actual_generation=str(runtime),automatic_retry=False)
  new(global_claims/(key+'.json'),caller.encode_result(obj));new(runtime/'finite-call-claims'/(key+'.START.json'),caller.encode_result(obj));return obj
 outputs=[];calls=0;returns=[]
 for result in caller.iter_species_versions(species,observations[50000],observations[25000],backend,check,claim):
  check();variant=result['variant'];p=out/(variant+'-complete-result.json');new(p,caller.encode_result(result));calls+=result['spatial_jackknife_calls'];returns.append(dict(variant=variant,status=result['status'],backend_status=result.get('backend_result',{}).get('status') if isinstance(result.get('backend_result'),dict) else None));outputs.append(ref(p))
  if result.get('claim_receipt'):
   key=hashlib.sha256(result['claim_receipt']['operation_id'].encode()).hexdigest();new(runtime/'finite-call-claims'/(key+'.RETURNED.json'),caller.encode_result(dict(operation_id=result['claim_receipt']['operation_id'],source_result_ref=outputs[-1],saved_complete_return=True,no_retry=True)))
 need([x['variant'] for x in returns]==cfg['allowed_versions'],'exactfinite3versions retained')
 summary=dict(status='FIRST_SPECIES_FINITE_A50_K50_A25_RETURNED_PENDING_ROOT',packet=cfg['packet'],native_species_name=species.native_species_name,assembly=species.assembly,source_receipts=receipts,source_measurement_acceptance_ref=cfg['measurement_acceptance_ref'],frozen_preoutcome_ref=cfg['preoutcome_gate_ref'],fixed_J_count=len(species.ids),results=outputs,result_statuses=returns,numeric_call_counts=dict(spatial_jackknife=calls,separate_full_fit=0,new_working_REML=0,Ctop_species_refit=0,prepare_primary_or_J_weights_C=0,new_RM=0,API=0),seconds=time.monotonic()-start,all_NA_constant_failure_details_preserved=True,working_S_not_reselected=True)
 new(out/'terminal.json',caller.encode_result(summary));check();return summary
