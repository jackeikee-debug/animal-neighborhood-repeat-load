"""Saved accepted summaries -> true2x2/common-scale tables. No model/body reexecution."""
from pathlib import Path
import collections,csv,importlib.util,json,sys,time
PURPOSES={'OFFDIAGONAL_xA50_TO_K50':('2x2','xA50','K50','B0_ORIGINAL_J'),'OFFDIAGONAL_xK50_TO_A50':('2x2','xK50','A50','B0_ORIGINAL_J'),'FIXED_COMMON_SCALE_50000':('scale_common','xA50','A50','B0_A25_ACTIVE_IDS_STRUCTURE50'),'FIXED_COMMON_SCALE_25000':('scale_common','xA50','A25','B0_A25_ACTIVE_IDS_STRUCTURE50')}
def tables(index,baseline,logical,summaries,number,baseline_ref):
 base={(r['species'],r['assembly'],r['variant']):r for r in baseline};assert len(base)==len(baseline);subjects={r['native_species_id']:r for r in index};assert len(subjects)==len(index);rows=[]
 def row(sp,analysis,purpose,discovery,target,cohort,z,kind,record=None,resultref=None):
  old=base[(sp['species'],sp['assembly'],target)];assert old['role']==sp['role'];record=record or {};beta=number(z.get('beta'));full=z['full_status'];assert (beta is None) or full=='OK'
  return dict(native_species_id=sp['native_species_id'],species=sp['species'],assembly=sp['assembly'],role=old['role'],baseline_analysis_layer=old['layer'],fixed_taxon_layer=sp['layer'],in_frozen_working_S=old['in_frozen_working_S']=='True',analysis=analysis,purpose=purpose,discovery_predictor=discovery,target_response=target,cohort=cohort,original_B0_J=int(old['frozen_J']),planned_active_n=record.get('planned_active_n'),n_fit=None if z.get('n') in (None,'NA') else int(z['n']),beta=beta,full_status=full,spatial_status=z['status'],full_NA_reason=z.get('full_reason') if full=='NA' else None,spatial_NA_reason=z.get('reason') if z['status']=='NA' else None,typed_source_NA_reason=record.get('reason') if record.get('status')=='TYPED_NA_NO_NUMERIC_CALL' else None,original_B0_same_response_beta=number(old['beta']),original_B0_NA_reason_literal=old['reason'],return_kind=kind,numeric_fingerprint=z.get('numeric_fingerprint'),consumer_signature_sha=(record.get('consumer_signature') or {}).get('sha256'),selection_excluded=record.get('selection_excluded',[]),typed_NA_ids=record.get('typed_NA_ids',[]),source_result_ref=resultref,logical_labels=record.get('logical_labels',{}),original_baseline_table_ref=baseline_ref,original_NA_record_preserved=True)
 for sid,sp in subjects.items():
  for v,x in [('A50','xA50'),('K50','xK50')]:
   old=base[(sp['species'],sp['assembly'],v)];z=dict(beta=old['beta'],n=old['n_fit'],full_status=old['full_status'],status=old['spatial_status'],reason=old['reason'],full_reason=old['reason']);rows.append(row(sp,'2x2','BASELINE_'+x+'_TO_'+v,x,v,'B0_ORIGINAL_J',z,'REUSE_ACCEPTED_B0_DIAGONAL',resultref=baseline_ref))
 seen=set()
 for r in logical:
  sid=r['native_species_id'];sp=subjects[sid];assert (r['species'],r['assembly'])==(sp['species'],sp['assembly']);assert (sid,r['purpose']) not in seen;seen.add((sid,r['purpose']));analysis,x,target,cohort=PURPOSES[r['purpose']];state=r['status'];resultref=None
  if state=='NEW_UNIQUE_NUMERIC_INPUT':z=summaries[r['numeric_fingerprint']];resultref=z['source_result_ref'];assert z['native_species_id']==sid and z['consumer_signature']==r['consumer_signature'];kind='NEW_COMPLETE_RETURN_BODY_QA_ACCEPTED'
  elif state in ('REUSE_EXACT_CONSUMER_INPUT_RETURN','ALIAS_TO_NEW_UNIQUE_INPUT_PENDING_RETURN'):
   reused=r['reused'];assert reused['consumer_signature']==r['consumer_signature']
   if reused['kind']=='NEW_PHASE2_INPUT_PENDING_NUMERIC_RETURN':z=summaries[reused['numeric_fingerprint']];kind='ALIAS_TO_NEW_NOW_COMPLETE_VERIFIED_RETURN';resultref=z['source_result_ref']
   else:z=reused['summary'];kind=reused['kind'];resultref=reused['source_result_ref'];assert kind in ('ACCEPTED_B0_RETURN','ACCEPTED_PC_RETURN')
  else:
   assert state=='TYPED_NA_NO_NUMERIC_CALL';z=dict(beta=None,n=None,full_status='NA',status='NA',reason=r['reason'],full_reason=r['reason']);kind='TYPED_NA_NO_NUMERIC_CALL'
  rows.append(row(sp,analysis,r['purpose'],x,target,cohort,z,kind,r,resultref))
 assert len(seen)==4*len(index) and len(rows)==6*len(index)
 keys={(r['native_species_id'],r['analysis'],r['discovery_predictor'],r['target_response']) for r in rows};assert len(keys)==len(rows);return rows
def build(cfg,out,check,events):
 start=time.monotonic();s=importlib.util.spec_from_file_location('p2report_sources',cfg['sources_code']);source=importlib.util.module_from_spec(s);s.loader.exec_module(source);src=source.Sources(cfg,check,events);s=importlib.util.spec_from_file_location('p2report_metrics',cfg['metrics_code']);metrics=importlib.util.module_from_spec(s);s.loader.exec_module(metrics);out=Path(out);out.mkdir()
 def read(n):return json.loads(''.join(src.lines(n)))
 qa=read('QA_acceptance');assert qa['status']=='ROOT_ACCEPTED_ALL114_COMPLETE_SAVED_RETURN_BODIES_AND_EXACTSTOP' and qa['complete_new_pairs']==114
 for name,ref in cfg['required_refs'].items():assert all(src.sources[name][k]==ref[k] for k in ('path','bytes','sha256'))
 assert all(src.sources['verified_summaries'][k]==qa['verified_summary_ref'][k] for k in ('path','bytes','sha256'))
 summaries=read('verified_summaries');assert len(summaries)==114
 index=read('cache_index');assert len(index)==51;baseline=list(csv.DictReader(src.lines('baseline'),delimiter='\t'));assert len(baseline)==153;logical=[json.loads(l) for l in src.lines('logical')];assert len(logical)==204
 rows=tables(index,baseline,logical,summaries,metrics.number,{k:src.sources['baseline'][k] for k in ('path','bytes','sha256')});normal=[r for r in rows if r['role']=='new_normal_candidate'];lookup={(r['native_species_id'],r['analysis'],r['discovery_predictor'],r['target_response']):r for r in rows};normal_ids=sorted({r['native_species_id'] for r in normal},key=int);assert len(normal_ids)==45
 contrasts=[('TARGET_A50_PREDICTOR_K_VS_A',('2x2','xA50','A50'),('2x2','xK50','A50')),('TARGET_K50_PREDICTOR_A_VS_K',('2x2','xK50','K50'),('2x2','xA50','K50')),('PREDICTOR_A_TARGET_K_VS_A',('2x2','xA50','A50'),('2x2','xA50','K50')),('PREDICTOR_K_TARGET_K_VS_A',('2x2','xK50','A50'),('2x2','xK50','K50')),('FIXED_COMMON_SCALE_Y25_VS_Y50',('scale_common','xA50','A50'),('scale_common','xA50','A25')),('FIXED_COMMON_Y50_VS_B0_A50',('2x2','xA50','A50'),('scale_common','xA50','A50'))]
 comparison=[];pairedrows=[]
 for label,oldkey,newkey in contrasts:
  pairs=[]
  for sid in normal_ids:
   a=lookup[(sid,*oldkey)];b=lookup[(sid,*newkey)];pairs.append(dict(old=a['beta'],new=b['beta']));pairedrows.append(dict(contrast=label,native_species_id=sid,species=a['species'],in_frozen_working_S=a['in_frozen_working_S'],old_beta=a['beta'],new_beta=b['beta'],old_NA_reason=a['full_NA_reason'],new_NA_reason=b['full_NA_reason'],signed_delta=None if a['beta'] is None or b['beta'] is None else b['beta']-a['beta']))
  comparison.append(dict(contrast=label,population='fixed_original45_normal_candidates',**metrics.paired(pairs)))
 pairs=[]
 for sid in normal_ids:
  b=lookup[(sid,'scale_common','xA50','A25')];old=next(r for r in baseline if (r['species'],r['assembly'],r['variant'])==(b['species'],b['assembly'],'A25'));pairs.append(dict(old=metrics.number(old['beta']),new=b['beta']))
 comparison.append(dict(contrast='Y25_STRUCTURE50_VS_ORIGINAL_A25_STRUCTURE25',population='fixed_original45_normal_candidates',**metrics.paired(pairs)))
 cells=[]
 for key in sorted({(r['analysis'],r['discovery_predictor'],r['target_response']) for r in normal}):
  z=[r for r in normal if (r['analysis'],r['discovery_predictor'],r['target_response'])==key];v=[r['beta'] for r in z if r['beta'] is not None];assert len(z)==45
  import statistics
  cells.append(dict(analysis=key[0],discovery_predictor=key[1],target_response=key[2],fixed_normal_candidates=45,finite=len(v),NA=45-len(v),positive=sum(x>0 for x in v),negative=sum(x<0 for x in v),zero=sum(x==0 for x in v),median_beta=statistics.median(v) if v else None,minimum_beta=min(v) if v else None,maximum_beta=max(v) if v else None))
 def tsv(name,records):
  with (out/name).open('x') as f:
   w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader()
   for r in records:w.writerow({k:'NA' if v is None else json.dumps(v,ensure_ascii=False,sort_keys=True,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in r.items()})
 tsv('all51-true-discovery-target-2x2.tsv',[r for r in rows if r['analysis']=='2x2']);tsv('all51-fixed-common-member-scale.tsv',[r for r in rows if r['analysis']=='scale_common']);tsv('fixed45-descriptive-paired-contrasts.tsv',comparison);tsv('fixed45-paired-species-ledger.tsv',pairedrows);tsv('fixed45-six-cell-effect-distributions.tsv',cells)
 assert len(rows)==306 and sum(r['return_kind']=='REUSE_ACCEPTED_B0_DIAGONAL' for r in rows)==102 and set(src.done)==set(src.sources)
 terminal=dict(status='ALL51_TRUE2X2_AND_COMMON_SCALE_TABLES_COMPLETE_PENDING_ROOT',all51_cells=306,original_diagonals_reused=102,prepared_logical_states=dict(collections.Counter(r['status'] for r in logical)),new_body_accepted_returns=114,coefficient_distributions=cells,contrasts=comparison,all_NA_OG725_S35_preserved=True,original_fixed_S_not_reselected=True,original_6644UNKNOWN1='untouched',source_receipts=src.done,new_models=0,old_refits=0,p_value=None,calibrated_confidence_interval=None,seconds=time.monotonic()-start)
 source.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(saved_summary_table_join=1,new_or_old_fits=0),source_receipts=src.done)
