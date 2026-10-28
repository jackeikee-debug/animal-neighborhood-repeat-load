"""First new unique numerical inputs -> original spatial returns, durable once per fingerprint."""
from pathlib import Path
import datetime,gzip,hashlib,importlib.util,json,os,sys,time
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def syncdir(p):
 fd=os.open(p,os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)
def atomic_new(p,raw):
 p=Path(p);tmp=p.with_name(p.name+'.publish-'+str(os.getpid()))
 with tmp.open('xb') as f:f.write(raw);f.flush();os.fsync(f.fileno())
 os.link(tmp,p,follow_symlinks=False);syncdir(p.parent);tmp.unlink();syncdir(p.parent)
def build(cfg,out,check,events):
 start=time.monotonic();source=load('effect_fit_Sources',cfg['codes']['sources']);src=source.Sources(cfg,check,events);helper=load('effect_fit_jsonsafe',cfg['codes']['json_safe']);restore=load('effect_fit_restore',cfg['codes']['restore_bytes']);adapter=load('effect_fit_adapter',cfg['codes']['effect_inputs']);backend=load('effect_fit_original_backend',cfg['codes']['backend']);caller=load('effect_fit_original_once_caller',cfg['codes']['original_caller']);out=Path(out);out.mkdir();(out/'results').mkdir();runtime=out.parent;(runtime/'model-call-claims').mkdir()
 def encode(x):return json.dumps(helper.json_safe(x),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()+b'\n'
 def read(n):
  parts=[];length=0
  for line in src.lines(n):length+=len(line.encode());assert length<=32*2**20;parts.append(line)
  return restore.restore_exact_bytes(json.loads(''.join(parts)))
 manifest=read('prepared_manifest');assert manifest['status']=='NEW_FINITE_NUMERIC_INPUTS_PREPARED_NO_FITS_PENDING_ROOT' and manifest['new_or_old_model_calls']==0
 accepted=read('input_acceptance');contracts=load('effect_compute_contracts',Path(__file__).with_name('contract_checks.py'));contracts.check_descriptor(src.sources['prepared_manifest'],cfg['input_manifest_ref']);contracts.check_descriptor(src.sources['input_acceptance'],cfg['input_acceptance_ref']);contracts.verify_compute_acceptance(accepted,src.sources['prepared_manifest'],manifest);assert set(cfg['expected_numeric_fingerprints'])==set(manifest['input_refs']) and len(cfg['expected_numeric_fingerprints'])==len(manifest['input_refs'])
 claims=Path(cfg['global_claim_root']);assert claims.is_dir();calls=0;outputs=[]
 for fp,desc in sorted(manifest['input_refs'].items()):
  check();name='input_'+fp;d=src.sources[name];assert d['sha256']==desc['sha256'] and d['bytes']==desc['bytes'] and d['path']==desc['path']
  inp=read(name);assert name in src.done and inp['status']=='PREPARED_NUMERIC_INPUT' and inp['schema']=='FROZEN_SENSITIVITY_NUMERIC_INPUT_V1'
  assert inp['numeric_fingerprint']==fp==adapter.input_fingerprint(inp['effective'],helper.json_safe)
  assert inp['native_species_name']==manifest['native_species_name'] and inp['assembly']==manifest['assembly'] and not inp['active_NA_ids']
  variant=inp['variant'];assert variant in ('A50','K50','A25');op='SPATIAL_EFFECT_NUMERIC:'+fp
  claimref=None
  def claim(operation,inputs):
   nonlocal claimref
   check();assert operation==op and inputs['numeric_fingerprint']==fp
   obj=dict(schema='ROOT_ONCE_FINITE_CALL_CLAIM_V1',approved_by='current-conversation-root',operation_id=operation,status='CLAIMED_ONCE',inputs_immutable_bound=True,inputs=inputs,first_actual_generation=str(runtime),automatic_retry=False,at_utc=now())
   raw=encode(obj);atomic_new(claims/(fp+'.json'),raw);atomic_new(runtime/'model-call-claims'/(fp+'.START.json'),raw);claimref=obj;return obj
  inputbinding={'numeric_fingerprint':fp,'prepared_input_ref':desc,'input_acceptance_ref':cfg['input_acceptance_ref'],'backend_sha256':adapter.BACKEND_SHA,'variant':variant,'native_species_name':inp['native_species_name'],'assembly':inp['assembly']}
  kwargs={'narrow_rows':inp['narrow_rows']} if variant=='A25' else {}
  try:result,token=caller._invoke(op,inputbinding,lambda:backend.spatial_jackknife(inp['main_panel'],inp['observations'],variant=variant,**kwargs),check,claim)
  except BaseException as error:
   if claimref is not None:atomic_new(runtime/'model-call-claims'/(fp+'.UNKNOWN.json'),encode({'operation_id':op,'at_utc':now(),'reason':type(error).__name__+':'+str(error),'automatic_retry':False,'backend_return_or_complete_output_not_accepted':True}))
   raise
  calls+=1;assert result['excluded']==inp['excluded'];full=result.get('full',{})
  if full.get('status')=='OK':assert full['ids']==inp['active_ids'] and full['weights']==inp['effective']['weights']
  record={'schema':'ORIGINAL_BACKEND_SENSITIVITY_COMPLETE_RETURN_V1','numeric_fingerprint':fp,'variant':variant,'native_species_name':inp['native_species_name'],'assembly':inp['assembly'],'prepared_input_ref':desc,'input_acceptance_ref':cfg['input_acceptance_ref'],'claim_receipt':token,'backend_result':result,'old_results_not_recomputed':True}
  raw=encode(record);assert len(raw)<=4*2**20;path=out/'results'/(fp+'.json.gz')
  with path.open('xb') as f:
   with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,compresslevel=6) as gz:gz.write(raw)
   f.flush();os.fsync(f.fileno())
  packed=path.read_bytes();ref={'path':str(path),'bytes':len(packed),'sha256':hashlib.sha256(packed).hexdigest(),'gzip':True,'uncompressed_bytes':len(raw),'uncompressed_sha256':hashlib.sha256(raw).hexdigest()}
  summary={'numeric_fingerprint':fp,'variant':variant,'source_result_ref':ref,'status':result['status'],'reason':result.get('reason'),'full_status':full.get('status'),'n':full.get('n'),'beta':full.get('beta'),'vJ':result.get('vJ'),'spatial_eligible':result.get('information',{}).get('spatial_eligible',False),'deletions':len(result.get('deletions',[])),'deletion_min':result.get('deletion_min'),'deletion_max':result.get('deletion_max'),'sign_changes':result.get('sign_changes'),'p_value':None,'calibrated_confidence_interval':None,'NA_negative_uncertain_complete_return_retained':True}
  returned={'operation_id':op,'at_utc':now(),'source_result_ref':ref,'summary':summary,'saved_complete_return':True,'automatic_retry':False};atomic_new(runtime/'model-call-claims'/(fp+'.RETURNED.json'),encode(returned));atomic_new(claims/(fp+'.RETURNED.json'),encode(returned));outputs.append(summary);del inp,result,record,raw
 assert calls==len(manifest['input_refs']) and set(src.done)==set(src.sources)
 terminal={'status':'NEW_UNIQUE_EFFECT_INPUTS_ALL_ORIGINAL_BACKEND_RETURNS_SAVED_PENDING_ROOT','native_species_id':manifest['native_species_id'],'native_species_name':manifest['native_species_name'],'assembly':manifest['assembly'],'input_manifest_ref':cfg['input_manifest_ref'],'results':outputs,'source_receipts':src.done,'numeric_call_counts':{'new_spatial_jackknife':calls,'separate_full_fit':0,'old_or_baseline_refits':0,'primary_reweighting':0,'new_RM_QC_API':0},'seconds':time.monotonic()-start,'all_NA_originalS_OG725_retained':True}
 source.write_json(out/'terminal.json',terminal);return {'numeric_call_counts':terminal['numeric_call_counts'],'source_receipts':src.done}
