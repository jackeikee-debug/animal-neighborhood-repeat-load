"""First compact assembly of NEW disjoint R23 rows; preserves all old typed NA."""
from pathlib import Path
import collections,csv,gzip,importlib.util,json,math,statistics,time

def build(cfg,out,check,events):
 start=time.monotonic();s=importlib.util.spec_from_file_location('r23collectsource',cfg['sources_code']);source=importlib.util.module_from_spec(s);s.loader.exec_module(source);src=source.Sources(cfg,check,events);out=Path(out);out.mkdir();read=lambda n:json.loads(''.join(src.lines(n)));index=read('cache_index');admission=read('admission');supp=read('large_metadata_supplement');assert len(index)==51 and admission['profile_count']==38;subject={x['native_species_id']:x for x in index};reason={x['native_species_id']:x for x in admission['records']};data={};expected={};proofs={}
 for sid in cfg['explicit_supported_native_ids']:
  p=read('proof'+sid);assert p['all_owned_exactSTOP']and p['native_species_id']==sid;J=p.get('targets',p.get('all_frozen_targets'));assert J==subject[sid]['J'];expected[sid]=J;proofs[sid]=p
  for n in cfg['count_source_names'][sid]:
   for r in csv.DictReader(src.lines(n),delimiter='\t'):
    check();assert r['native_species_id']==sid and(r['species'],r['assembly'])==(subject[sid]['species'],subject[sid]['assembly']);g=r['old_group'];rad=int(r['radius']);assert rad in(50000,25000);key=(sid,g,rad);assert key not in data;counts={k:int(r[k])for k in('F','O','I','D')};assert all(v>=0 for v in counts.values())and sum(counts[k]for k in('F','O','I'))==counts['D'];data[key]=dict(target_key=r['target_key'],wid=int(r['wid']),**counts)
  assert sum(k[0]==sid for k in data)==2*J
 assert len(expected)==32
 aggregate=collections.defaultdict(list);allrows=0;sourceNA=0;keys=set();newNA=collections.Counter()
 fields=['native_species_id','species','assembly','old_group','target_key','radius','weight','composition_status','composition_NA_reason','F','O','I','D','F_fraction','O_fraction','I_fraction','original_repeat_pair_status','original_repeat_pair_reason','repeat_D_literal_at_radius']
 with gzip.open(out/'all51-frozen-J-gene-composition-and-distinct-NA.tsv.gz','xt',encoding='utf-8',newline='')as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader()
  for old in csv.DictReader(src.lines('nested_counts'),delimiter='\t'):
   check();sid,g=old['native_species_id'],old['old_group'];assert (sid,g)not in keys;keys.add((sid,g));sourceNA+=old['status']=='TYPED_NA';weight=float(old['weight']);assert weight>0 and math.isfinite(weight);pair={}
   for radius in(50000,25000):
    r=dict(native_species_id=sid,species=old['species'],assembly=old['assembly'],old_group=g,target_key=old['target_key'],radius=radius,weight=weight,original_repeat_pair_status=old['status'],original_repeat_pair_reason=old['reason'],repeat_D_literal_at_radius=old['D'+str(radius//1000)])
    if sid in expected:
     c=data.pop((sid,g,radius));assert c['target_key']==old['target_key']
     if r['repeat_D_literal_at_radius']!='NA':assert c['D']==int(r['repeat_D_literal_at_radius'])
     pair[radius]=c;r.update(composition_status='SOURCE_BOUND_COMPLETE_GEOMETRIC_PARTITION',composition_NA_reason=None,**{k:c[k]for k in('F','O','I','D')});r.update({k+'_fraction':c[k]/c['D']if c['D'] else None for k in('F','O','I')});aggregate[sid,radius].append((weight,c))
    else:
     rec=reason[sid];r.update(composition_status='TYPED_NA_DIAGNOSTIC_SOURCE_UNAVAILABLE',composition_NA_reason=rec['status']+':'+','.join(rec.get('reason_codes',[])),F=None,O=None,I=None,D=None,F_fraction=None,O_fraction=None,I_fraction=None);newNA[sid]+=1
    w.writerow({k:'NA'if r[k]is None else r[k]for k in fields});allrows+=1
   if pair:
    outer={k:pair[50000][k]-pair[25000][k]for k in('F','O','I','D')};assert all(v>=0 for v in outer.values())and sum(outer[k]for k in('F','O','I'))==outer['D'];aggregate[sid,0].append((weight,outer))
 assert not data and len(keys)==29328 and sourceNA==9 and allrows==58656 and set(src.done)==set(src.sources)
 rows=[];status=[]
 for p in index:
  sid=p['native_species_id'];ok=sid in expected
  status.append(dict(native_species_id=sid,species=p['species'],assembly=p['assembly'],role=p['role'],layer=p['layer'],frozen_J=p['J'],prior_profile_status=p['profile'],composition_status='ALL_FROZEN_J_TWO_RADII_ACCEPTED'if ok else reason[sid]['status']if sid in reason else 'NO_ACCEPTED_MEASURED_PROFILE',new_other_gene_coverage_assumed_zero=False,original_RM_NA_unchanged=True))
  for radius in(50000,25000,0):
   values=aggregate.get((sid,radius),[]);swD=math.fsum(weight*c['D']for weight,c in values);positive=[(weight,c)for weight,c in values if c['D']>0];sw=math.fsum(weight for weight,c in positive)
   r=dict(native_species_id=sid,species=p['species'],assembly=p['assembly'],scale='annulus25to50'if radius==0 else str(radius//1000),J=p['J'],measured_n=len(values),D_zero_n=len(values)-len(positive),status='COMPLETE_WINDOW_BASED_COMPOSITION'if ok else 'TYPED_NA',reason=None if ok else status[-1]['composition_status'],sum_weighted_D=swD if values else None)
   for k in('F','O','I'):
    r[k+'_ratio_weighted_counts_over_weighted_D']=math.fsum(weight*c[k]for weight,c in values)/swD if swD else None;r[k+'_weighted_mean_fraction_Dpositive']=math.fsum(weight*c[k]/c['D']for weight,c in positive)/sw if sw else None;r[k+'_median_fraction_Dpositive']=statistics.median(c[k]/c['D']for weight,c in positive)if positive else None
   rows.append(r)
 def tsv(n,data):
  with(out/n).open('x',newline='')as f:
   writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');writer.writeheader();writer.writerows({k:'NA'if v is None else v for k,v in row.items()}for row in data)
 tsv('all51-gene-composition-source-and-NA-status.tsv',status);tsv('all51-inner-annulus-gene-composition-summary.tsv',rows)
 source.write_json(out/'terminal.json',dict(status='R23_ALL51_DISPOSITION_AND32_COMPLETE_FROZEN_J_GENE_COMPOSITION_PENDING_ROOT',all51_species_statuses=51,fully_composed_species=32,composed_target_pairs=sum(expected.values()),composed_two_radius_windows=2*sum(expected.values()),all38_target_rows=29328,all38_two_radius_rows=58656,original_repeat_pair_NA_units_retained=9,composition_NA_window_counts_by_native_id=dict(newNA),unprofiled_original_states_retained=13,fullANNledger_NA=admission['ANN_ledger_profile_typedNA_ids'],legacy_specific_binding_gap=admission['prelaunch_metadata_HOLD_ids'],inner_annulus_partition_closure=True,window_based_not_disjoint_genome_coverage=True,not_functional_intron_or_repeat_compartment_allocation=True,new_SQL_old_QC_models=0,source_receipts=src.done,seconds=time.monotonic()-start,all_NA_OG725_preserved=True))
 assert sum(p.stat().st_size for p in out.iterdir())<=8*2**20
 return dict(numeric_call_counts=dict(new_saved_composition_assembly=1,new_SQL_QC_model=0),source_receipts=src.done)
