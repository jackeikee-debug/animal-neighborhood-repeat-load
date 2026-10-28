"""Saved R15 accepted summaries -> all48 finite/typedNA cells, no refits."""
from pathlib import Path
import collections,csv,importlib.util,json,statistics,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();source=load('r15reportsource',cfg['sources_code']);src=source.Sources(cfg,check,events);metrics=load('r15metrics',cfg['metrics_code']);out=Path(out);out.mkdir()
 read=lambda n:json.loads(''.join(src.lines(n)))
 proof=read('QA_acceptance');assert proof['complete_pairs']==23 and proof['all_owned_exactSTOP'];assert all(src.sources['summaries'][k]==proof['verified_summary_ref'][k]for k in ('path','bytes','sha256'));summaries=read('summaries');assert len(summaries)==23
 logical=[json.loads(l)for l in src.lines('logical')];baseline=list(csv.DictReader(src.lines('baseline'),delimiter='\t'));base={(r['species'],r['assembly']):r for r in baseline if r['variant']=='A50'};assert len(base)==51 and len(logical)==48
 rows=[]
 for r in logical:
  check();b=base[(r['species'],r['assembly'])];state=r['status'];z=None
  if state=='NEW_UNIQUE_NUMERIC_INPUT':z=summaries[r['numeric_fingerprint']]
  elif state=='ALIAS_TO_NEW_UNIQUE_INPUT_PENDING_RETURN':z=summaries[r['reused']['numeric_fingerprint']]
  else:assert state=='TYPED_NA_NO_NUMERIC_CALL'
  if z is not None:assert z['consumer_signature']==r['consumer_signature'] and (z['species'],z['assembly'])==(r['species'],r['assembly'])
  beta=metrics.number(z['beta'])if z else None;old=metrics.number(b['beta']);s=r['summary']
  row=dict(native_species_id=r['native_species_id'],species=r['species'],assembly=r['assembly'],layer=r['layer'],rule=r['rule'],role=b['role'],in_original_S=b['in_frozen_working_S'],original_J=s['original_J_count'],common_J=s['common_count'],retention_fraction=s['retention_fraction'],original_xA_min=s['original_J_xA50_min'],original_xA_max=s['original_J_xA50_max'],common_xA_min=s['common_xA50_min'],common_xA_max=s['common_xA50_max'],beta=beta,baseline_A50_beta=old,delta_from_baseline=None if beta is None or old is None else beta-old,full_status=z['full_status']if z else 'NA',spatial_status=z['status']if z else 'NA',reason=(z.get('reason')if z else r['reason']),selected_sourceNA=r['selected_sourceNA_ids'],original_sourceNA_ids=list(r['original_sourceNA_records']),original_sourceNA_outside_common=r['sourceNA_outside_predefined_common_ids'],conditional_subset_not_original_NA_repair=True,return_kind=state,source_result_ref=z['source_result_ref']if z else None)
  rows.append(row)
 group=[]
 for key in sorted({(r['layer'],r['rule'])for r in rows}):
  rs=[r for r in rows if (r['layer'],r['rule'])==key];vals=[r['beta']for r in rs if r['beta']is not None];pairs=[dict(old=r['baseline_A50_beta'],new=r['beta'])for r in rs]
  group.append(dict(layer=key[0],rule=key[1],species=len(rs),finite_full_beta=len(vals),typedNA=len(rs)-len(vals),spatial_variance_OK=sum(r['spatial_status']=='OK'for r in rs),common_J=rs[0]['common_J'],median_retention=statistics.median(r['retention_fraction']for r in rs),minimum_retention=min(r['retention_fraction']for r in rs),maximum_retention=max(r['retention_fraction']for r in rs),median_beta=statistics.median(vals)if vals else None,positive=sum(x>0 for x in vals),negative=sum(x<0 for x in vals),baseline_pair_comparison=metrics.paired(pairs),meaning='Conditional within-frozen-clade shared oldgroups; eachspecies own accepted protein, originalweights/structure50; no tip independence or wholeclade sampling claim'))
 def tsv(n,data):
  with(out/n).open('x',newline='')as f:
   w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows({k:'NA'if v is None else json.dumps(v,ensure_ascii=False,separators=(',',':'))if isinstance(v,(list,dict))else v for k,v in r.items()}for r in data)
 tsv('all48-clade-common-member-effects-and-NA.tsv',rows);tsv('all10-frozen-clade-rule-descriptions.tsv',group)
 assert set(src.done)==set(src.sources)
 terminal=dict(status='ALL48_CLADE_COMMON_EFFECTS_AND_NA_REPORT_COMPLETE_PENDING_ROOT',logical_rows=48,finite_beta=sum(r['beta']is not None for r in rows),typed_NA=sum(r['beta']is None for r in rows),unique_complete_return_pairs=23,identical_input_aliases=2,group_summaries=group,source_receipts=src.done,new_models=0,old_refits=0,original_sourceNA_S35_OG725_unchanged=True,seconds=time.monotonic()-start)
 source.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(saved_report_join=1,new_or_old_models=0),source_receipts=src.done)
