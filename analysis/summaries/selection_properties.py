"""New bounded descriptive selection audit from a frozen coordinate frame; no refits."""
from pathlib import Path
import collections,csv,importlib.util,json,statistics,time
SP=('Glandiceps_talaboti','Ptychodera_flava','Saccoglossus_kowalevskii')
def span(start,end):
 if start in ('','NA') or end in ('','NA'):return None
 s=int(start);e=int(end)
 return e-s+1 if s>=1 and e>=s else None
def describe(values):
 x=sorted(v for v in values if v is not None)
 def quantile(q):
  z=(len(x)-1)*q;lo=int(z);hi=min(lo+1,len(x)-1);return x[lo]+(z-lo)*(x[hi]-x[lo])
 return dict(total=len(values),finite=len(x),NA=len(values)-len(x),median=statistics.median(x) if x else None,q25=quantile(.25) if x else None,q75=quantile(.75) if x else None,minimum=min(x) if x else None,maximum=max(x) if x else None)
def build(cfg,out,check,events):
 start=time.monotonic();s=importlib.util.spec_from_file_location('selection_sources',cfg['sources_code']);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);src=m.Sources(cfg,check,events);out=Path(out);out.mkdir();freeze=json.loads(''.join(src.lines('freeze')));assert freeze['status']=='FROZEN_R17_GENE_COORDINATE_PROPERTIES_BEFORE_NEW_FRAME_PROJECTION'
 selected={r['orthogroup'] for r in csv.DictReader(src.lines('selected'),delimiter='\t')};assert len(selected)==1896
 receipt=json.loads(''.join(src.lines('historical_receipt')));assert receipt['code']['sha256']==cfg['original_producer_sha256'];counts=collections.Counter();groups=collections.defaultdict(list);ontology=collections.Counter();seen=set();eligible=set()
 fields=['orthogroup','species','cohort','original_prior_eligible','gene_span_bp','transcript_span_bp','nearby_strict_anchor_representative_count','coordinate_system','coordinate_feature','source_ontology_group','source_ontology_review_clean','source_prior_exclusion_reasons','window_error']
 with (out/'all9214-selected-status-and-coordinate-properties.tsv').open('x') as f:
  writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
  for row in csv.DictReader(src.lines('frame'),delimiter='\t'):
   check();g=row['orthogroup'];assert g not in seen;seen.add(g);assert row['prior_eligible'] in ('true','false');prior=row['prior_eligible']=='true'
   if prior:eligible.add(g)
   cohort='selected1896' if g in selected else 'initially_eligible_not_selected' if prior else 'initially_excluded';counts[cohort]+=1;ontology[(cohort,row['ontology_group_from_valid_terms'],row['ontology_review_clean'])]+=1
   neighbor=json.loads(row['nearby_anchor_representative_gene_counts_by_species'])
   for sp in SP:
    p=sp+'_';gl=span(row[p+'gene_start'],row[p+'gene_end']);tl=span(row[p+'transcript_start'],row[p+'transcript_end']);n=neighbor.get(sp);assert n is None or type(n) is int and n>=0
    if gl is not None:assert row[p+'coordinate_system']=='GFF_1_based_closed'
    values=dict(orthogroup=g,species=sp,cohort=cohort,original_prior_eligible=prior,gene_span_bp=gl,transcript_span_bp=tl,nearby_strict_anchor_representative_count=n,coordinate_system=row[p+'coordinate_system'],coordinate_feature=row[p+'coordinate_feature'],source_ontology_group=row['ontology_group_from_valid_terms'],source_ontology_review_clean=row['ontology_review_clean'],source_prior_exclusion_reasons=row['prior_exclusion_reasons'],window_error=row['window_error']);writer.writerow({k:'NA' if v is None else v for k,v in values.items()})
    for metric,val in [('gene_span_bp',gl),('transcript_span_bp',tl),('nearby_strict_anchor_representative_count',n)]:groups[(cohort,sp,metric)].append(val)
 assert len(seen)==9214 and selected<=seen and counts['selected1896']==1896
 exact7495=len(eligible)==7495 and selected<=eligible
 summaries=[dict(cohort=c,species=sp,metric=metric,**describe(v)) for (c,sp,metric),v in sorted(groups.items())]
 with (out/'selected-versus-unselected-property-distributions.tsv').open('x') as f:
  w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t',lineterminator='\n');w.writeheader()
  for r in summaries:w.writerow({k:'NA' if v is None else v for k,v in r.items()})
 with (out/'selection-ontology-availability.tsv').open('x') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['cohort','source_ontology_group','source_review_clean','orthogroups'])
  for (c,g,r),n in sorted(ontology.items()):w.writerow([c,g,r,n])
 assert set(src.done)==set(src.sources)
 terminal=dict(status='DESCRIPTIVE_SELECTION_PROPERTY_PROJECTION_COMPLETE_PENDING_ROOT' if exact7495 else 'FRAME_PROJECTION_COMPLETE_EXACT7495_UNIVERSE_NOT_MATCHED_NO_SUBSTITUTION',frame=9214,initial_eligible=len(eligible),selected=1896,selected_in_initial_eligible=len(selected&eligible),exact7495_universe_bound=exact7495,cohort_counts=dict(counts),source_coverage='original three-species coordinate frame; gene/transcript spans and nearby strict-anchor representative counts; not new full gene or GC/exon/repeat measurements',unselected_window_repeat_GC_ALL_exon_fields='UNMEASURED_NOT_IMPUTED',selection_fraction=1896/7495 if exact7495 else None,new_models=0,new_QC_or_source_API=0,source_receipts=src.done,seconds=time.monotonic()-start,all_old_results_OG725_NA_retained=True)
 m.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(new_descriptive_selection_audit=1,new_models=0),source_receipts=src.done)
