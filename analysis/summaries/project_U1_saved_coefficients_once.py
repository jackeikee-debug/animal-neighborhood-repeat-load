"""U1 only: saved B0 full coefficients -> new fixed score projection, no refit."""
from pathlib import Path
import collections,gc,hashlib,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();pure=load('U1pure',cfg['adapter']);source=load('U1sources',cfg['sources_code']);restore=load('U1restore',cfg['restore_code']);consumer=load('U1consumer',cfg['consumer_code']);src=source.Sources(cfg,check,events);out=Path(out);out.mkdir()
 def read(n,restore_all=False):
  parts=[];size=0
  for line in src.lines(n):size+=len(line.encode());assert size<=32*2**20;parts.append(line)
  obj=json.loads(''.join(parts));return restore.restore_exact_bytes(obj)if restore_all else obj
 policy=read('policy');graphproof=read('graph_acceptance');graph=read('graph');index=read('cache_index');known=read('accepted_consumer_index');assert graphproof['status']=='ROOT_ACCEPTED_U0_RESULTBLIND_42_GLOBAL_BLOCKS_P31_ALLSTOP'and graphproof['all_owned_exactSTOP'];assert all(src.sources['graph'][k]==graphproof['graph_ref'][k]for k in('path','bytes','sha256'));assert graph['P_species']==policy['P_species']and graph['P_native_ids']==policy['P_native_ids']and graph['component_count']==42;index={x['native_species_id']:x for x in index};b0={k:v for k,v in known.items()if v['kind']=='ACCEPTED_B0_RETURN'and v['variant']=='A50'};assert len(b0)==34;del known;gc.collect()
 scores={};records={};witness=[]
 for sid,sp in zip(policy['P_native_ids'],policy['P_species']):
  check();desc=index[sid]['cache'];assert all(src.sources['cache'+sid][k]==desc[k]for k in('path','bytes','sha256'));raw=read('cache'+sid);assert raw['native_species_id']==sid and raw['native_species_name']==sp and raw['unchanged_J_weights_S'];panel=restore.restore_exact_bytes(raw['base']['panel']);observations=restore.restore_exact_bytes(raw['observations']['50000']);assembly=raw['assembly'];del raw;gc.collect()
  active=dict(baseline='B0',species=sp,assembly=assembly,variant='A50',rows=panel['rows'],weights=panel['weights'],x=[r['xA']for r in panel['rows']],observations=observations);assert [r['id']for r in active['rows']]==panel['ids'];sig=consumer.consumer_signature(active);old=b0[sig['sha256']];assert all(old['source_result_ref'][k]==src.sources['result'+sid][k]for k in('path','bytes','sha256'));saved=read('result'+sid,restore_all=True);score=pure.accepted_full_to_scores(active,saved,graph['block_map'],policy)
  assert score['species']==sp and score['old_fit_calls']==0 and score['new_nuisance_projection_calls']<=1
  # U3 later uses this complete-file-bound full beta only, never reevaluates a fit.
  full=saved['backend_result']['full'];records[sp]=dict(native_species_name=sp,assembly=assembly,variant='A50',backend_result=dict(full={k:full.get(k)for k in('status','beta','reason')}));scores[sp]=score
  witness.append(dict(native_species_id=sid,species=sp,assembly=assembly,source_cache_ref=desc,source_result_ref=old['source_result_ref'],accepted_original_root_ref=old['accepted_root_ref'],consumer_signature=sig,projection_status=score['status'],projection_reason=score['reason'],complete_old_return_body_consumed_once=True))
  del saved,full,panel,observations,active,score;gc.collect()
 assert set(scores)==set(policy['P_species'])and len(scores)==31 and set(src.done)==set(src.sources)
 source.write_json(out/'all31-fixed-B0-A50-block-scores-and-numeric-NA.json',scores);source.write_json(out/'all31-accepted-B0-full-beta-only-for-U3.json',records)
 source.write_json(out/'terminal.json',dict(status='U1_FIXED31_SAVED_COEFFICIENT_SCORE_PROJECTIONS_COMPLETE_PENDING_ROOT',P_species=policy['P_species'],original_S35=policy['original_S35'],projection_status_counts=dict(collections.Counter(x['status']for x in scores.values())),projection_NA_reasons={sp:r['reason']for sp,r in scores.items()if r['status']!='SCORES_PREPARED'},new_nuisance_projections=sum(x['new_nuisance_projection_calls']for x in scores.values()),old_or_new_model_fits=0,spatial_REML_optimizer_calls=0,source_receipts=src.done,source_witness=witness,policy_id=policy['policy_id'],all_original_NA_S35_OG725_retained=True,seconds=time.monotonic()-start))
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20
 return dict(numeric_call_counts=dict(new_nuisance_score_projections=sum(x['new_nuisance_projection_calls']for x in scores.values()),old_model_fits=0,spatial_REML_optimizer_calls=0),source_receipts=src.done)
