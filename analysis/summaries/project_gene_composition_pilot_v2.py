"""R23 first exact readonly source join and16-target pilot. No annotation/QC/model."""
from pathlib import Path
import collections,csv,gc,gzip,hashlib,importlib.util,json,sqlite3,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def build(cfg,out,check,events):
 start=time.monotonic();source=load('r23sources',cfg['sources_code']);src=source.Sources(cfg,check,events);restore=load('r23restore',cfg['restore_code']);safe=load('r23safe',cfg['safe_code']);pure=load('r23pure',cfg['adapter']);out=Path(out);out.mkdir()
 def read(n):
  obj=json.loads(''.join(src.lines(n)))
  # Only lossless source-bound fields consumed by this NEW purpose remain
  # resident; all original source bytes are still SHA/EOF/stat checked once.
  if n=='cache':
   obj={k:obj[k]for k in('native_species_id','native_species_name','assembly','unchanged_J_weights_S','source_lineage','accepted_nested_geometry')} | {'base':{'panel':{'ids':obj['base']['panel']['ids']}}}
  elif n=='mask':
   fields=('id','target_key','wid','radius','seqid','start0','end0','D','complete','masks','mask_basis')
   obj={'ordered_primary_ids':obj['ordered_primary_ids'],'windows':[{k:w[k]for k in fields}for w in obj['windows']]}
  return restore.restore_exact_bytes(obj)
 cert=read('cache_acceptance');c=read('cache');mask=read('mask');ann=read('ANN_acceptance');gen=read('GENOME_config');assert cert['accepted_profiles']==38 and cert['all_owned_exactSTOP'];assert c['native_species_id']==cfg['native_species_id'] and c['unchanged_J_weights_S'] and mask['ordered_primary_ids']==c['base']['panel']['ids'];assert all(c['source_lineage']['accepted_mask_source'][k]==src.sources['mask'][k]for k in('path','bytes','sha256'))
 db=cfg['readonly_annotation_DB'];assert source.snap(db['path'])==db['closed_stat']==ann['closed_DB_stat'];assert all(db[k]==ann['closed_DB'][k]for k in('path','bytes','sha256'));assert gen['annotation_db_sha256']==db['sha256'] and gen['annotation_manifest_sha256']==ann['manifest_canonical_sha256'] and gen['native_species_name']==ann['native_species_name']==c['native_species_name'] and gen['assembly']==ann['assembly']==c['assembly'];ids=c['base']['panel']['ids'];assert len(ids)<=2000
 pins=dict(annotation_db_sha256=db['sha256'],genome_annotation_db_sha256=gen['annotation_db_sha256'],annotation_manifest_sha256=ann['manifest_canonical_sha256'],genome_annotation_manifest_sha256=gen['annotation_manifest_sha256'],p2_mask_source_sha256=c['source_lineage']['accepted_mask_source']['sha256'],loaded_mask_source_sha256=src.sources['mask']['sha256'],p2_assembly=c['assembly'],annotation_assembly=ann['assembly'],p2_species=c['native_species_name'],annotation_species=ann['native_species_name'])
 windows={(w['id'],w['radius']):w for w in mask['windows']};assert len(windows)==len(mask['windows']);budget=pure.QueryBudget(max_targets=2000);plans=[];rows=[];source.write_json(out/'SQL-query-once-intent.json',dict(operation='R23_EXACT_SOURCE_JOIN_ALL_FROZEN_WINDOWS_AND_FIRST16_TARGET_FEATURES',readonly_annotation_DB=db,source_pins=pins,selection='first16 originalJ IDs, no outcome selection',targets=ids[:16],all_targets_bound_once=len(ids),radii=[50000,25000],max_index_feature_gene_hits_each=200000,queries_schema=cfg['adapter_sha256'],no_full_DB_hash_or_scan=True))
 con=sqlite3.connect(Path(db['path']).as_uri()+'?mode=ro&immutable=1',uri=True);con.execute('PRAGMA query_only=ON');con.execute('PRAGMA cache_size=-4096');con.execute('PRAGMA mmap_size=0');con.execute('PRAGMA temp_store=MEMORY')
 try:
  sqlplans=pure.query_plan_gate(con);s=con.execute('SELECT state,manifest_hash,receipt_json FROM source WHERE id=1').fetchone();assert s is not None and s[1]==ann['manifest_canonical_sha256'];source.write_json(out/'actual-SQL-index-source-metadata-gate.json',dict(plans=sqlplans,source_state=s[0],manifest_hash=s[1],receipt_json_sha256=hashlib.sha256(s[2].encode()if isinstance(s[2],str)else s[2]).hexdigest(),readonly=True,cache_limit=4*2**20,mmap=0))
  focal_proofs={}
  for i,g in enumerate(ids):
   check();geometry=c['accepted_nested_geometry'][g];target=None;pair=[]
   for radius in(50000,25000):
    mw=windows[(g,radius)];aw=pure.read_annotation_window(con,mw['wid'])
    if target is None:target=pure.read_annotation_target(con,geometry['target_key'])
    binding=pure.bind_focal_window(geometry,radius,mw,aw,target,pins);pair.append(binding)
    if i<16:
     if g not in focal_proofs:focal_proofs[g]=pure.exact_focal_feature(con,binding['focal_id'],binding['seqid'],*binding['focal'],budget=budget,check=check)
     features=pure.iter_gene_features(con,binding['seqid'],*binding['window'],budget=budget,window_key=(binding['target_key'],radius),check=check)
     result=pure.compose(window=binding['window'],combined_mask=binding['combined_mask'],D=binding['D'],focal=binding['focal'],focal_id=binding['focal_id'],seqid=binding['seqid'],gene_features=features)
     assert result['status']=='CALLER_BOUND_GEOMETRIC_COMPOSITION_NOT_FUNCTIONAL_INTRON';counts=result['counts'];rows.append(dict(native_species_id=c['native_species_id'],species=c['native_species_name'],assembly=c['assembly'],old_group=g,target_key=binding['target_key'],wid=binding['window_index'],radius=radius,**counts,other_gene_features_seen=result['full_gene_features_seen'],source_binding_accepted=True,not_functional_intron=True))
   plans.extend(pair)
  assert source.snap(db['path'])==db['closed_stat'] and set(src.done)==set(src.sources)
 finally:con.close()
 raw=json.dumps(safe.json_safe(dict(schema='R23_ALL_FROZEN_QUERY_BINDINGS_SOURCE_ACCEPTED_NO_FEATURE_REQUERY_FOR_PILOT',native_species_id=c['native_species_id'],species=c['native_species_name'],assembly=c['assembly'],plans=plans,pilot_completed_old_groups=ids[:16],full_J_ids=ids,readonly_annotation_DB=db,source_pins=pins)),ensure_ascii=False,separators=(',',':'),allow_nan=False).encode()+b'\n'
 with(out/'all-frozen-window-query-bindings.json.gz').open('xb')as f:
  with gzip.GzipFile(fileobj=f,mode='wb',mtime=0)as z:z.write(raw)
 with(out/'first16-targets-two-radii-composition.tsv').open('x',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 source.write_json(out/'terminal.json',dict(status='R23_FIRST16_TARGET_FEATURE_PILOT_AND_ALL_FROZEN_SOURCE_BINDINGS_COMPLETE_PENDING_ROOT',native_species_id=c['native_species_id'],all_frozen_targets_bound=len(ids),all_window_bindings=len(plans),pilot_targets=16,pilot_window_rows=len(rows),query_hit_counts=budget.counts,source_closed_stat_unchanged=True,all_target_raw_fullgene_maskwid_join_exact=True,new_ANN_Genome_RM_QC_models=0,source_receipts=src.done,source_pins=pins,seconds=time.monotonic()-start,all_NA_OG725_S35_unchanged=True))
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20
 return dict(numeric_call_counts=dict(exact_source_join=len(plans),new_bounded_target_gene_composition=32,new_or_old_models=0),source_receipts=src.done)
