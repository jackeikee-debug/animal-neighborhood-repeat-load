"""First R19/R22 saved nested-count descriptive projection, not a model."""
from pathlib import Path
import csv,gc,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();source=load('r19sources',cfg['sources_code']);src=source.Sources(cfg,check,events);pure=load('r19pure',cfg['adapter']);out=Path(out);out.mkdir()
 def read(n):return json.loads(''.join(src.lines(n)))
 proof=read('cache_acceptance');index=read('cache_index');assert proof['accepted_profiles']==38 and proof['all_owned_exactSTOP'];ctx=dict(accepted_nested_mask_source_proof=True,species={},expected_total_rows=29328,expected_source_NA_rows=9,source_ref={k:src.sources['counts'][k]for k in ('path','bytes','sha256')})
 expected={r['native_species_id']:r for r in proof['cache_refs']}
 for sid in sorted(expected,key=int):
  check();s=src.sources['cache'+sid];assert all(s[k]==expected[sid][k]for k in ('path','bytes','sha256'));c=read('cache'+sid);panel=c['base']['panel'];assert c['native_species_id']==sid and c['unchanged_J_weights_S'];ctx['species'][sid]=dict(J_count=len(panel['ids']),ordered_J_ids=panel['ids'],weights=panel['weights'],species=c['native_species_name'],assembly=c['assembly']);del c;gc.collect()
 result=pure.summarize_stream(csv.DictReader(src.lines('counts'),delimiter='\t'),ctx,check)
 assert len(result['species'])==38 and set(src.done)==set(src.sources)
 result.update(status='ALL38_FROZEN_WINDOW_LOAD_AND_NESTED_DESCRIPTIONS_COMPLETE_PENDING_ROOT',all51_input_status=index,unavailable_profiles_retained=[r for r in index if not r['cache']],source_receipts=src.done,new_RM_DNA_QC_or_fits=0,original_NA_records_preserved=True,seconds=time.monotonic()-start)
 source.write_json(out/'all51-window-based-load-and-nested-descriptions.json',result)
 flat=[]
 for sp in result['species']:
  for scale,metrics in sp['scales'].items():
   for metric,z in metrics.items():
    r=dict(native_species_id=sp['native_species_id'],species=sp['species'],assembly=sp['assembly'],scale=scale,metric=metric,J_count=sp['J_count'],sourceNA_n=sp['source_NA_n'])
    r.update({'frozen_'+k:v for k,v in z['frozen_J'].items()});r.update({'observed_subset_'+k:v for k,v in z['conditional_observed_subset'].items()});flat.append(r)
 with(out/'per-species-inner-outer-load-summary.tsv').open('x',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(flat[0]),delimiter='\t');w.writeheader();w.writerows({k:'NA'if v is None else v for k,v in r.items()}for r in flat)
 source.write_json(out/'terminal.json',dict(status='SAVED_NESTED_COUNTS_WINDOW_DESCRIPTIONS_COMPLETE_PENDING_ROOT',species_profiles=38,all51_statuses=len(index),target_rows=result['total_rows'],original_NA_rows=result['source_NA_rows'],summary_rows=len(flat),source_receipts=src.done,new_models=0,original_NA_OG725_preserved=True,seconds=time.monotonic()-start))
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20
 return dict(numeric_call_counts=dict(saved_window_count_descriptions=1,new_or_old_models=0),source_receipts=src.done)
