"""Resolve the original7495 universe by saved availability ledger, not reselection."""
from pathlib import Path
import collections,csv,importlib.util,json,time
def build(cfg,out,check,events):
 start=time.monotonic();s=importlib.util.spec_from_file_location('r17_sources',cfg['sources_code']);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);src=m.Sources(cfg,check,events);s=importlib.util.spec_from_file_location('r17_stats',cfg['statistics_code']);stats=importlib.util.module_from_spec(s);s.loader.exec_module(stats);out=Path(out);out.mkdir()
 authority=json.loads(''.join(src.lines('freeze')));assert authority['status']=='ROOT_FROZEN_R17_EXACT7495_SAVED_LEDGER_JOIN_NO_RESELECTION';available={};eligible=set()
 for r in csv.DictReader(src.lines('availability'),delimiter='\t'):
  g=r['orthogroup'];assert g not in available and r['geometric_candidate'] in ('1','0');available[g]=r
  if r['geometric_candidate']=='1':eligible.add(g)
 assert len(available)==9214 and len(eligible)==7495
 groups=collections.defaultdict(list);counts=collections.Counter();seen=set();selected=set();ont=collections.Counter();fields=None
 with (out/'exact7495-universe-coordinate-property-rows.tsv').open('x') as f:
  for r in csv.DictReader(src.lines('properties'),delimiter='\t'):
   check();g=r['orthogroup'];key=(g,r['species']);assert g in available and key not in seen;seen.add(key);a=available[g]
   chosen=r['cohort']=='selected1896'
   if chosen:selected.add(g);assert g in eligible
   cohort='selected1896' if chosen else 'geometric_eligible_not_selected5599' if g in eligible else 'geometrically_excluded1719';r['cohort']=cohort;r['original_geometric_eligible']=a['geometric_candidate'];r['original_geometric_exclusion_reasons']=a['reasons'];counts[(cohort,r['species'])]+=1
   if fields is None:fields=list(r);writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
   writer.writerow(r)
   for metric in ('gene_span_bp','transcript_span_bp','nearby_strict_anchor_representative_count'):groups[(cohort,r['species'],metric)].append(None if r[metric]=='NA' else int(r[metric]))
   if r['species']=='Glandiceps_talaboti':ont[(cohort,r['source_ontology_group'],r['source_ontology_review_clean'])]+=1
 assert len(seen)==3*9214 and len(selected)==1896 and selected<=eligible
 summaries=[dict(cohort=c,species=sp,metric=metric,**stats.describe(v)) for (c,sp,metric),v in sorted(groups.items())]
 with (out/'exact7495-selected-versus-unselected-distributions.tsv').open('x') as f:
  w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader()
  for r in summaries:w.writerow({k:'NA' if v is None else v for k,v in r.items()})
 with (out/'exact7495-ontology-availability.tsv').open('x') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['cohort','source_ontology_group','source_review_clean','orthogroups'])
  for (c,g,r),n in sorted(ont.items()):w.writerow([c,g,r,n])
 assert set(src.done)==set(src.sources)
 terminal=dict(status='EXACT_ORIGINAL7495_SELECTION_PROPERTY_JOIN_COMPLETE_PENDING_ROOT',original_frame=9214,original_geometric_eligible=7495,selected=1896,eligible_not_selected=5599,geometrically_excluded=1719,selection_fraction=1896/7495,source_receipts=src.done,prior_eligible9044_is_distinct_gate=True,old_protected_outcomes_read=0,new_models_QC_RM_API=0,seconds=time.monotonic()-start,all_old_results_preserved=True)
 m.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(new_exact_universe_property_join=1,models=0),source_receipts=src.done)
