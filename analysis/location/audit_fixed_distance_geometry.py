"""First geometry-only audit; join saved statistics without any numerical replay."""
from pathlib import Path
import csv,gzip,hashlib,json
import numpy as np
from location_geometry import allowed_point,admit_direct,canonical_mode
Q=Path(__import__('os').environ['ANRL_R4_ROOT']);R=Q.parent
gfile=R/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz';raw=gfile.read_bytes()
bind=json.loads((R/'geometry-design-preparation-syntaxfix-once/Root-paired-geometry-X-row-order-binding.json').read_text());assert hashlib.sha256(raw).hexdigest()==bind['geometry_sha256']
geo=json.loads(gzip.decompress(raw));assert len(geo['ids'])==1896 and geo['ids']==bind['oldgroup_ids']
phase=R/'all1998-offsets-geometry-freeze-first';meta=json.loads((phase/'return.json').read_text());plans=json.loads((phase/'allowed-phase-plans.json').read_text())
streams=[]
for sp in ['G','S']:
 s=geo['species'][sp]
 for ch in s['chromosomes']:
  W=[(t['start0'],t['end0'])for t in s['targets']if t['seqid']==ch];L=s['lengths'][ch]
  assert W and L>=2*max(b-a for a,b in W)
  streams.append(dict(species=sp,chromosome=ch,L=L,W=W,B=s['forbidden'].get(ch,[]),n=len(W),widths=np.array([b-a for a,b in W],dtype=np.int64)))
assert len(streams)==45 and sum(x['n']for x in streams)==3792
cases=[];audit=[];grid=[250,500,1000,2500,5000,10000,25000,50000,100000]
for mode in ['circ','linear']:
 for distance in grid:
  for Gsign,Ssign in [(1,1),(1,-1),(-1,1),(-1,-1)]:
   case=f'{mode}_d{distance}_G{Gsign:+d}_S{Ssign:+d}';blocked=[];no_wrap=True
   for x in streams:
    sp=x['species'];ch=x['chromosome'];delta=distance*(Gsign if sp=='G'else Ssign)
    membership=allowed_point(plans[mode][sp][ch]['allowed'],delta,x['L'],mode=='circ')
    direct,reason=admit_direct(x['W'],x['B'],delta,x['L'],mode=='circ')
    assert membership==direct,('independent admission mismatch',case,sp,ch,membership,direct)
    lo,hi=plans['linear'][sp][ch]['domain'];no_wrap=no_wrap and lo<=delta<hi
    row=dict(case=case,mode=mode,distance_bp=distance,species=sp,chromosome=ch,targets=x['n'],signed_delta=delta,admissible=int(membership),reason=reason or'PASS')
    audit.append(row)
    if not membership:blocked.append(dict(species=sp,chromosome=ch,reason=reason))
   admitted=not blocked
   key=f'{canonical_mode(mode,bool(no_wrap))}:G{distance*Gsign}:S{distance*Ssign}'if admitted else None
   cases.append(dict(case=case,mode=mode,distance_bp=distance,G_delta=distance*Gsign,S_delta=distance*Ssign,admissible=admitted,canonical_key=key,all_whole_windows_no_wrap=bool(no_wrap),blocked_chromosomes=blocked,per_species_members=1896))
assert len(cases)==72 and len(audit)==3240
unique={}
for x in cases:
 if x['admissible']:unique.setdefault(x['canonical_key'],x)
payload=dict(status='GEOMETRY_ONLY_FIXED72_SCENARIOS_FROZEN_BEFORE_NEW_RESPONSES',grid_bp=grid,patterns=[[1,1],[1,-1],[-1,1],[-1,-1]],cases=cases,unique_admitted_cases=list(unique.values()),zero_reference='REUSE_ORIGINAL_OBSERVED_NO_QUERY_NO_REFIT',new_repeat_track_reads=0,new_response_projections=0,original_fits=0)
Path('fixed-distance-admission-freeze.json').write_text(json.dumps(payload,indent=2)+'\n')
with Path('per-chromosome-fixed-distance-admission.tsv').open('x')as f:
 w=csv.DictWriter(f,fieldnames=list(audit[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(audit)
summary=[];joined=[];base_total=sum(int(x['widths'].sum())for x in streams);weights=np.array([x['n']for x in streams],dtype=np.int64)
for mode in ['circ','linear']:
 item=next(x for x in meta['outputs']if x['mode']==mode);p=phase/item['path'];assert hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
 offsets=np.load(p,allow_pickle=False);assert offsets.shape==(999,45)and offsets.dtype==np.dtype('int64')and item['columns']==[[x['species'],x['chromosome']]for x in streams]
 with(R/(mode+'-all999-first')/'all999-conditional-background-statistics.tsv').open()as f:old=list(csv.DictReader(f,delimiter='\t'))
 assert len(old)==999 and[int(x['draw'])for x in old]==list(range(1,1000))
 local=[]
 for i,row in enumerate(old):
  delta=offsets[i];dist=np.array([min(int(v)%x['L'],(-int(v))%x['L'])if mode=='circ'else abs(int(v))for v,x in zip(delta,streams)],dtype=np.int64)
  retained=sum(int(np.maximum(x['widths']-int(dd),0).sum())for dd,x in zip(dist,streams))
  order=np.argsort(dist);cw=np.cumsum(weights[order]);v=dist[order];med=(int(v[np.searchsorted(cw,1896,side='left')])+int(v[np.searchsorted(cw,1897,side='left')]))/2
  rec=dict(mode=mode,draw=i+1,T_reused=float(row['T']),target_weighted_mean_abs_distance_bp=float(np.dot(dist,weights)/3792),target_weighted_median_abs_distance_bp=med,whole_window_self_overlap_bp=retained,whole_window_total_bp=base_total,whole_window_self_overlap_fraction=retained/base_total,target_fraction_moved_at_least100kb=float(weights[dist>=100000].sum()/3792))
  assert np.isfinite(rec['T_reused'])and 0<=retained<=base_total
  local.append(rec);joined.append(rec)
 a=np.array([x['whole_window_self_overlap_fraction']for x in local]);d=np.array([x['target_weighted_median_abs_distance_bp']for x in local]);summary.append(dict(mode=mode,draws=999,self_overlap_fraction_min=float(a.min()),self_overlap_fraction_median=float(np.median(a)),self_overlap_fraction_max=float(a.max()),target_median_displacement_min=float(d.min()),target_median_displacement_median=float(np.median(d)),target_median_displacement_max=float(d.max())))
with Path('saved1998-distance-and-footprint-T.tsv').open('x')as f:
 w=csv.DictWriter(f,fieldnames=list(joined[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(joined)
counts={mode:{str(d):sum(x['admissible']for x in cases if x['mode']==mode and x['distance_bp']==d)for d in grid}for mode in ['circ','linear']}
Path('return.json').write_text(json.dumps(dict(status='R4_GEOMETRY_AND_SAVED1998_DESCRIPTIVE_JOIN_COMPLETE',logical_scenarios=72,admitted=sum(x['admissible']for x in cases),unique_admitted=len(unique),admission_counts=counts,saved_draw_geometry=summary,new_response_projections=0,new_source_reads=0,old_models_refitted=0,old_nulls_reexecuted=0,full_1896_retained=True),indent=2)+'\n')
