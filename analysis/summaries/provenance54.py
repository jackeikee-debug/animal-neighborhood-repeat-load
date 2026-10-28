"""Lossless exact-assembly join of saved provenance tables, no external calls."""
from pathlib import Path
import csv,importlib.util,json,time
def join(roster,new51,historical54,old_db):
 def index(rows):
  out={r['assembly']:r for r in rows};assert len(out)==len(rows);return out
 modern=index(new51);historical=index(historical54);assert len(modern)==51 and len(historical)==54 and len(roster)==54
 db={r['species']:r for r in old_db};assert set(db)=={'G','P','S'}
 oldnames={'Glandiceps_talaboti__GCF_964340395_1':'G','Ptychodera_flava__GCF_041260155_1':'P','Saccoglossus_kowalevskii__GCF_978046115_1':'S'};out=[]
 for r in roster:
  assembly=r['assembly'];old=r['fixed_study_role']=='old_discovery';table=historical if old else modern;source=table[assembly];assert source['species']==r['species_display_name']
  row=dict(native_species_name=r['native_species_name'],species_display_name=r['species_display_name'],assembly=assembly,fixed_role=r['fixed_study_role'],fixed_taxon_layer=r['fixed_major_layer'],provenance_table='historical54_original3' if old else 'accepted_v11_expansion51',historical_is_not_current_callability=True,source_row=dict(source))
  if old:
   code=oldnames[r['native_species_name']];row['old3_database_provenance']=db[code];row['classifier_reference_database_status']='UNESTABLISHED_RELEASE_AND_CONTENT_HASH_NOT_RECORDED';row['current_expansion_J']='NOT_APPLICABLE_DISCOVERY1896_SCOPE';row['current_response_scope']='original_discovery_saved_results_retained'
  else:
   assert source['current_role']==r['fixed_study_role'];row['old3_database_provenance']=None;row['classifier_reference_database_status']=source['classifier_database_version_hash'];row['current_expansion_J']=source['current_frozen_J'];row['current_response_scope']=source['current_A50_status']
  out.append(row)
 assert len({r['assembly'] for r in out})==54 and sum(r['fixed_role']=='old_discovery' for r in out)==3;return out
def build(cfg,out,check,events):
 start=time.monotonic();spec=importlib.util.spec_from_file_location('provenance_sources',cfg['sources_code']);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);src=m.Sources(cfg,check,events);out=Path(out);out.mkdir()
 roster=json.loads(''.join(src.lines('roster')))['roster54']['records'];new=list(csv.DictReader(src.lines('new51'),delimiter='\t'));old=list(csv.DictReader(src.lines('historical54'),delimiter='\t'));db=list(csv.DictReader(src.lines('old_db'),delimiter='\t'));rows=join(roster,new,old,db)
 with (out/'annotation-provenance-by-all54-species.jsonl').open('x') as f:
  for r in rows:f.write(json.dumps(r,ensure_ascii=False,sort_keys=True)+'\n')
 # Flat table retains every original column, with source versions unchanged.
 fields=list(dict.fromkeys(k for r in rows for k in r['source_row']))
 with (out/'annotation-provenance-by-all54-species.tsv').open('x') as f:
  columns=['native_species_name','species_display_name','fixed_role','fixed_taxon_layer','provenance_table','classifier_reference_database_status','current_expansion_J','current_response_scope','old3_database_provenance']+fields;w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow(dict({k:(json.dumps(r[k],ensure_ascii=False,sort_keys=True) if isinstance(r[k],dict) else ('NA_NOT_APPLICABLE' if r[k] is None else r[k])) for k in columns if k in r},**{k:r['source_row'].get(k,'NA_NOT_RECORDED_IN_THIS_SOURCE') for k in fields}))
 assert set(src.done)==set(src.sources)
 terminal=dict(status='ALL54_EXACT_ASSEMBLY_PROVENANCE_JOIN_COMPLETE_PENDING_ROOT',species=54,new51=51,old3=3,reference_database_unknown_preserved=True,rows_not_current_eligibility=True,new_API_or_scientific_calls=0,source_receipts=src.done,seconds=time.monotonic()-start)
 m.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(provenance_join=1,models=0),source_receipts=src.done)
