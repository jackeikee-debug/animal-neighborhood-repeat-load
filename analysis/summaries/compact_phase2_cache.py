"""One accepted compact-source pass; exact reusable inputs and count projections, no fits."""
from pathlib import Path
import collections,csv,gzip,hashlib,importlib.util,json,sys,time
def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m
def valid_counts(c):
 assert all(type(c[k]) is int and c[k]>=0 for k in ('A','C','U','I','K','D'))
 assert c['A']==c['C']+c['U']-c['I'] and c['K']==c['A']-c['U'] and c['A']<=c['D'] and c['I']<=min(c['C'],c['U']) and c['K']<=c['C']
def nested_geometry(w50,w25):
 assert all(w50[k]==w25[k] for k in ('target_key','id','seqid')) and w50['complete'] and w25['complete']
 assert w50['end0']-w50['start0']==100000 and w25['start0']==w50['start0']+25000 and w25['end0']==w50['end0']-25000
 for w in (w50,w25):
  right=w['start0'];bp=0
  for a,b in w['masks']:assert right<=a<b<=w['end0'];right=b;bp+=b-a
  assert w['D']==w['end0']-w['start0']-bp
 clipped=[[max(a,w25['start0']),min(b,w25['end0'])] for a,b in w50['masks'] if max(a,w25['start0'])<min(b,w25['end0'])]
 assert clipped==w25['masks'],'same accepted combined mask restricted to nested window'
 return dict(target_key=w50['target_key'],id=w50['id'],seqid=w50['seqid'],start50=w50['start0'],end50=w50['end0'],start25=w25['start0'],end25=w25['end0'],masks50=w50['masks'],masks25=w25['masks'],D50=w50['D'],D25=w25['D'],same_mask_clip_exact=True)
def count_pair(c50,c25):
 if c50 is None or c25 is None:return dict(status='TYPED_NA',reason='AT_LEAST_ONE_RADIUS_SOURCE_NA',outer=None)
 for c in (c50,c25):valid_counts(c)
 outer={k:c50[k]-c25[k] for k in ('A','C','U','I','K','D')};valid_counts(outer)
 return dict(status='TRUSTED_NESTED_COUNTS',reason=None,outer=outer)
def build(cfg,out,check,events):
 start=time.monotonic();source=load('phase2_sources',cfg['codes']['sources']);src=source.Sources(cfg,check,events);safe=load('phase2_safe',cfg['codes']['json_safe']);restore=load('phase2_restore',cfg['codes']['restore_bytes']);proj=load('phase2_original_projection',cfg['codes']['projection']);out=Path(out);out.mkdir();(out/'species').mkdir()
 def read(n):
  parts=[];total=0
  for line in src.lines(n):total+=len(line.encode());assert total<32*2**20;parts.append(line)
  return restore.restore_exact_bytes(json.loads(''.join(parts)))
 def gznew(path,obj):
  raw=json.dumps(safe.json_safe(obj),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()+b'\n';assert len(raw)<32*2**20
  with path.open('xb') as f:
   with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,compresslevel=6) as z:z.write(raw)
  packed=path.read_bytes();return dict(path=str(path),bytes=len(packed),sha256=hashlib.sha256(packed).hexdigest(),gzip=True,uncompressed_bytes=len(raw),uncompressed_sha256=hashlib.sha256(raw).hexdigest())
 contract=read('contract');assert contract['status']=='ROOT_FROZEN_BOUNDED_PHASE2_EXTENSION_UNDER_DIRECT_CONTINUATION_AND_PROSPECTIVE_CONSENT';rolecfg=read('roster');roles={(r['native_species_name'],r['assembly']):r for r in rolecfg['roster54']['records']};binding=read('binding');assert len(binding['packets'])==51
 refs=[];summaries=[];Jsets={};windows=0;oldNA=collections.Counter();profilecount=0
 loadfile=(out/'nested-window-counts-and-NA.tsv').open('x');fields=['native_species_id','species','assembly','old_group','target_key','weight','status','reason','D50','D25','A50','A25','U50','U25','K50','K25','outer_A','outer_C','outer_U','outer_I','outer_K','outer_D'];writer=csv.DictWriter(loadfile,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
 for descriptor in binding['packets']:
  sid=descriptor['native_species_id'];plan=read('plan'+sid);assert plan['native_species_id']==sid;role=roles[(plan['native_species_name'],plan['assembly'])]
  if plan['profile'] is None:
   summaries.append(dict(native_species_id=sid,species=plan['native_species_name'],assembly=plan['assembly'],role=role['fixed_study_role'],layer=role['fixed_major_layer'],profile='NO_ACCEPTED_MEASURED_PROFILE',J=plan['recipes'][0]['main_J_count'],cache=None,NA_preserved=True));continue
  check();profilecount+=1;profile=plan['profile'];base=read('panel'+sid);assert (base['native_species_name'],base['assembly'])==(plan['native_species_name'],plan['assembly'])
  mask=read('mask'+sid);assert mask['native_species_name']==base['native_species_name'] and mask['assembly']==base['assembly'] and mask['ordered_primary_ids']==base['panel']['ids']
  bywindow={(w['id'],w['radius']):w for w in mask['windows']};mask_events={(e['old_group_from_cache'],e['radius']):e for e in mask['source_event_ledger'] if e['radius'] in (25000,50000)}
  ids=base['panel']['ids'];primary=set(ids);local=set(base['panel']['local_descriptive_ids']);obs={50000:{},25000:{}};lineage={50000:{},25000:{}};seen=set();localrefs=[];coverage={50000:{},25000:{}};targetkeys={};physical=0;geometry={}
  for line in src.lines('measurement'+sid):
   check();physical+=1;r=restore.restore_exact_bytes(json.loads(line));proj.project_saved_observation(r,physical,hashlib.sha256(line.encode()).hexdigest(),profile['source_measurement_ref'],primary,local,seen,obs,lineage,localrefs)
   if r['old_group'] in primary:
    g=r['old_group'];radius=r['radius'];coverage[radius][g]=r['coverage'];targetkeys.setdefault(g,r['target_key']);assert targetkeys[g]==r['target_key'] and r['target_key']==bywindow[g,radius]['target_key'] and r['window_raw_event_ref']==mask_events[g,radius]
  assert physical==profile['source_rows'] and len(seen)==2*(len(primary)+len(local)) and len(localrefs)==2*len(local) and all(set(z)==primary for z in obs.values()) and len(ids)==profile['root_metadata']['fixed_J_count']
  counts=collections.Counter();byid={r['id']:r for r in base['panel']['rows']};narrow=base['A25'].get('narrow_rows',{})
  for i,g in enumerate(ids):
   geometry[g]=nested_geometry(bywindow[g,50000],bywindow[g,25000]);c50=coverage[50000][g];c25=coverage[25000][g];pair=count_pair(c50,c25);counts[pair['status']]+=1;windows+=1
   if c50 is not None:assert c50['D']==byid[g]['nonexonic_acgt']==bywindow[g,50000]['D']
   if c25 is not None:assert c25['D']==bywindow[g,25000]['D']
   if c25 is not None and g in narrow:assert c25['D']==narrow[g]['nonexonic_acgt']
   if pair['status']=='TYPED_NA':oldNA[sid]+=1
   row=dict(native_species_id=sid,species=base['native_species_name'],assembly=base['assembly'],old_group=g,target_key=targetkeys[g],weight=base['panel']['weights'][i],status=pair['status'],reason=pair['reason'])
   for k in ('D','A','U','K'):
    for radius in (50000,25000):c=coverage[radius][g];row[k+str(radius//1000)]=c[k] if c is not None else None
   for k in ('A','C','U','I','K','D'):row['outer_'+k]=pair['outer'][k] if pair['outer'] else None
   writer.writerow({k:'NA' if v is None else v for k,v in row.items()})
  cache=dict(schema='ROOT_COMPACT_PHASE2_ACCEPTED_BASELINE_AND_PAIRED_OBSERVATIONS_V1',native_species_id=sid,native_species_name=base['native_species_name'],assembly=base['assembly'],base=base,observations=obs,coverage=coverage,target_keys=targetkeys,accepted_nested_geometry=geometry,source_lineage=dict(plan_ref=cfg['plan_refs'][sid],baseline_panel_ref=plan['source_panel'],measurement_source_ref=profile['source_measurement_ref'],measurement_acceptance_ref=profile['measurement_acceptance_ref'],accepted_mask_source=cfg['mask_refs'][sid],old_root_ref=profile['accepted_root_ref']),unchanged_J_weights_S=True,new_models=0)
  ref=gznew(out/'species'/('native-'+sid+'.json.gz'),cache);refs.append(dict(native_species_id=sid,**ref));summaries.append(dict(native_species_id=sid,species=base['native_species_name'],assembly=base['assembly'],role=role['fixed_study_role'],layer=role['fixed_major_layer'],profile='ACCEPTED_SAVED_PROFILE_REUSED',J=len(ids),cache=ref,nested_count_states=dict(counts),NA_preserved=True));Jsets[sid]=set(ids);del cache,base,obs,coverage,lineage,mask,bywindow,mask_events,geometry
 loadfile.close();assert profilecount==38 and len(refs)==38 and sum(oldNA.values())==9 and set(oldNA)=={'2','7','16','32'}
 # Rules were frozen before new count projection: all nonempty normal J and
 # highest-QC-J representatives; no effect values used for member selection.
 clades=collections.defaultdict(list)
 for s in summaries:
  if s['role']=='new_normal_candidate' and s['J']>0:clades[s['layer']].append(s)
 cladeplans=[]
 for layer,records in sorted(clades.items()):
  for label,selected in [('ALL_NONEMPTY_NORMAL_J',records),('TWO_LARGEST_FROZEN_J',sorted(records,key=lambda x:(-x['J'],int(x['native_species_id'])))[:2])]:
   assert all(x['native_species_id'] in Jsets for x in selected);common=set.intersection(*(Jsets[x['native_species_id']] for x in selected));ids=sorted(common)
   cladeplans.append(dict(layer=layer,rule=label,species=[x['species'] for x in selected],native_species_ids=[x['native_species_id'] for x in selected],common_old_group_ids=ids,common_count=len(ids),original_J_counts={x['native_species_id']:x['J'] for x in selected},sourceNA_retained=sorted(set(oldNA)&{x['native_species_id'] for x in selected},key=int),no_outcome_or_sign_selection=True,meaning='same homologous old group with each species own accepted protein, not identical cross-species protein string'))
 source.write_json(out/'frozen-clade-common-member-plans.json',cladeplans);source.write_json(out/'species-cache-index.json',summaries);assert set(src.done)==set(src.sources)
 terminal=dict(status='PHASE2_COMPACT_INPUT_CACHE_AND_NESTED_COUNTS_COMPLETE_PENDING_ROOT',species=51,accepted_profiles=38,cache_refs=refs,main_J_windows=windows,sourceNA_units=dict(oldNA),clade_common_rules=len(cladeplans),new_models=0,new_primary_or_QC=0,raw_RM_DNA_API_reads=0,source_receipts=src.done,seconds=time.monotonic()-start,all_NA_OG725_preserved=True)
 source.write_json(out/'terminal.json',terminal);return dict(numeric_call_counts=dict(compact_source_projection=38,new_models=0),source_receipts=src.done)
