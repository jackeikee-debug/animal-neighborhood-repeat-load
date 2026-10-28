"""New purpose: exact saved-family link-base evidence accounting, no source recompute."""
from pathlib import Path
import csv,hashlib,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();source=load('r20sources',cfg['sources_code']);src=source.Sources(cfg,check,events);adapter=load('r20pure',cfg['adapter']);out=Path(out);out.mkdir()
 def rows(n):return csv.DictReader(src.lines(n),delimiter='\t')
 ctx=cfg['context'];ctx['source_refs']={n:dict(id=n,sha256=src.sources[n]['sha256'])for n in ('family_links','evidence','selected36','final36')}
 result=adapter.join_accumulate(rows('family_links'),rows('evidence'),rows('selected36'),rows('final36'),ctx,check)
 assert set(src.done)==set(src.sources) and result['new_search_or_model_calls']==0
 with (out/'family-link-base-evidence-strata.tsv').open('x',newline='') as f:
  fields=['species','axis','stratum','family_count','link_count','family_neighborhood_link_base_weight_bp'];w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(result['strata'])
 with (out/'all5540-family-weight-and-source-evidence.jsonl').open('x') as f:
  for r in result['family_audit']:f.write(json.dumps(r,ensure_ascii=False,sort_keys=True,separators=(',',':'))+'\n')
 summary={}
 for r in result['strata']:
  if r['axis']=='evidence_status':summary.setdefault(r['species'],[]).append(r)
 for sp,rs in summary.items():
  total=sum(r['family_neighborhood_link_base_weight_bp']for r in rs)
  for r in rs:r['within_species_link_base_fraction']=r['family_neighborhood_link_base_weight_bp']/total if total else None
 result.pop('family_audit');result.update(status='SAVED_UNKNOWN_COVERAGE_WEIGHTED_EVIDENCE_AVAILABILITY_COMPLETE_PENDING_ROOT',species_evidence_summaries=summary,full_source_descriptors={n:{k:s[k]for k in ('path','bytes','sha256')}for n,s in src.sources.items()},source_receipts=src.done,new_fits=0,source_RM_DNA_classifier_repeats=0,seconds=time.monotonic()-start,original_NA_and_OG725_untouched=True)
 source.write_json(out/'terminal.json',result)
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20,'8MiB finite output bound'
 return dict(numeric_call_counts=dict(saved_family_evidence_projection=1,new_or_old_models=0),source_receipts=src.done)
