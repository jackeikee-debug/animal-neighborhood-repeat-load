"""First result-blind U0 geometry/native-OG closure; no beta/model/random access."""
from pathlib import Path
import collections,csv,gc,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();source=load('u0sources',cfg['sources_code']);src=source.Sources(cfg,check,events);restore=load('u0restore',cfg['restore_code']);pure=load('u0pure',cfg['adapter']);out=Path(out);out.mkdir()
 read=lambda n:restore.restore_exact_bytes(json.loads(''.join(src.lines(n))))
 policy=read('policy');S=read('fixedS');index=read('cache_index');assert S['S_order']==policy['original_S35'] and len(S['S_order'])==35 and policy['selection_uses_beta_or_sign']is False;index={r['native_species_id']:r for r in index};assert len(index)==51
 mapping=collections.defaultdict(set);rowcount=0;states=collections.Counter();roles=collections.defaultdict(set)
 for r in csv.DictReader(src.lines('full_native_OG_membership'),delimiter='\t'):
  check();rowcount+=1;assert r['state']in('BOUND','EXPLICIT_OG725_EXCLUSION') and r['native_OGs'] and r['native_uids'];mapping[r['original_group']].update(r['native_OGs'].split(','));roles[r['original_group']].add(r['role']);states[r['state']]+=1
 assert rowcount==5688 and len(mapping)==1896 and all(len(x)==3 for x in roles.values()) and states=={'BOUND':5682,'EXPLICIT_OG725_EXCLUSION':6}
 records=[];counts={}
 for sid in policy['P_native_ids']:
  check();desc=index[sid]['cache'];assert all(src.sources['cache'+sid][k]==desc[k]for k in('path','bytes','sha256'));c=read('cache'+sid);assert c['native_species_id']==sid and c['unchanged_J_weights_S'] and c['native_species_name']in policy['P_species'];p=c['base']['panel'];assert len(p['ids'])==len(p['rows'])
  for row,g in zip(p['rows'],p['ids']):
   geometry=c['accepted_nested_geometry'][g];assert row['id']==g and row['real_chromosome']is True and (row['start0'],row['end0'])==(geometry['start50'],geometry['end50']) and geometry['same_mask_clip_exact']is True
   records.append(dict(species_id=sid,old_group=g,chromosome=row['location'],start0=row['start0'],end0=row['end0'],real_chromosome=True))
  counts[sid]=len(p['ids']);del c;gc.collect()
 frozen=dict(policy_id=policy['policy_id'],species_ids=policy['P_native_ids'],original_group_ids=sorted(mapping),window_bp=100000,strictly_positive_overlap=True,distance_band_bp=0,same_complete_native_OG=True,shared_old_group_across_species=True,result_blind=True)
 ogmap=dict(mapping_kind='COMPLETE_NATIVE_OG_MEMBERSHIP',complete=True,source_known=True,source_ref={k:src.sources['full_native_OG_membership'][k]for k in('path','bytes','sha256')},groups={g:dict(native_OGs=sorted(v),complete=True,source_known=True)for g,v in mapping.items()})
 result=pure.dependency_blocks(records,ogmap,frozen);assert len(counts)==31 and result['active_species_count']==31 and set(src.done)==set(src.sources)
 result.update(status='U0_GLOBAL_DEPENDENCY_CLOSURE_COMPLETE_PENDING_ROOT',original_S35=policy['original_S35'],P_species=policy['P_species'],P_native_ids=policy['P_native_ids'],sourceNA_four_original_records=policy['sourceNA_four_original_records'],per_species_window_counts=counts,original_membership_states=dict(states),new_models=0,new_beta_or_result_body_reads=0,source_receipts=src.done,seconds=time.monotonic()-start,all_original_NA_OG725_S35_preserved=True)
 result['U1_U2_graph_gate']='STOP_TOTAL_COMPONENTS_BELOW20_THUS_CANNOT_HAVE20_NONZERO_SCORE_BLOCKS' if result['component_count']<20 else 'SCORE_INFORMATION_GATE_STILL_REQUIRED_NO_GRAPH_SIZE_SUBSTITUTION'
 source.write_json(out/'global-dependency-blocks-and-fixed-P31.json',result);source.write_json(out/'terminal.json',{k:v for k,v in result.items()if k not in('components','block_map')})
 return dict(numeric_call_counts=dict(resultblind_dependency_projection=1,new_or_old_models=0,random_draws=0),source_receipts=src.done)
