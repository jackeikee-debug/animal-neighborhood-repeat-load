"""One999-vector shared-block working perturbation plus five arithmetic summaries."""
from pathlib import Path
import collections,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();pure=load('U2pure',cfg['adapter']);source=load('U2sources',cfg['sources_code']);src=source.Sources(cfg,check,events);out=Path(out);out.mkdir();read=lambda n:json.loads(''.join(src.lines(n)))
 policy=read('policy');proof=read('U1_acceptance');assert proof['status']=='ROOT_ACCEPTED_U1_ALL31_SAVED_COEFFICIENT_SCORE_PROJECTIONS_ALLSTOP'and proof['all_owned_exactSTOP']and proof['old_or_new_model_fits']==0
 for name,key in [('scores','scores_ref'),('accepted_betas','accepted_betas_ref')]:assert all(src.sources[name][k]==proof[key][k]for k in('path','bytes','sha256'))
 scores=read('scores');records=read('accepted_betas');check();joint=pure.joint_score_envelope(scores,policy);check();layers=pure.leave_layer_out_means(records,policy);assert len(joint['P_order'])==31 and len(joint['original_S35_rows'])==35 and joint['signs']['vectors']==999 and joint['model_fits']==layers['model_fits']==0 and set(src.done)==set(src.sources)
 source.write_json(out/'all-original35-and-P31-joint-block-working-ranges.json',joint);source.write_json(out/'five-leave-layer-out-descriptive-means.json',layers)
 source.write_json(out/'terminal.json',dict(status='U2_SHARED999_SIGNS_AND_U3_FIVE_ARITHMETIC_SUMMARIES_COMPLETE_PENDING_ROOT',P_species=policy['P_species'],original_S35=policy['original_S35'],species_status_counts=dict(collections.Counter(r['status']for r in joint['species_targets'].values())),theta_P=joint['P_conditional_equalweight_mean'],working_range_failed_reasons=dict(collections.Counter(reason for r in joint['species_targets'].values()for reason in r.get('gate',{}).get('reasons',[]))),original_S35_mean=joint['original_S35_mean_and_range'],signs=joint['signs'],source_receipts=src.done,old_or_new_models=0,new_joint_sign_vectors=999,new_layer_arithmetic_summaries=5,interpretation=joint['interpretation'],all_NA_OG725_preserved=True,seconds=time.monotonic()-start))
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20
 return dict(numeric_call_counts=dict(new_shared_globalblock_sign_vectors=999,new_leave_layer_out_arithmetic=5,old_or_new_model_fits=0),source_receipts=src.done)
