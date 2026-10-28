"""Exactly one released frozen 999-draw mode; no scheduling or retry."""
from pathlib import Path
import csv,gzip,hashlib,json,sys,time
import numpy as np
from physical_shift_root_numeric_gates import AnnotationUnion,prepare_callable_fragments,query_counts,fixed_qr,response_residuals,pearson
R3=Path(__import__('os').environ['ANRL_R3_ROOT'])
P=Path(__import__('os').environ['ANRL_PROJECT_ROOT'])
O=P/'results/migration/2026-10-04-mac-handoff/new-mac-root-v1/three-species-neighborhood-input-statistics-20261013-v1'
cfg=json.loads(Path(sys.argv[1]).read_text());mode=cfg['mode'];assert mode in ['circ','linear']
geo_path=R3/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz'
geo_raw=geo_path.read_bytes();geo=json.loads(gzip.decompress(geo_raw))
binding=json.loads((R3/'geometry-design-preparation-syntaxfix-once/Root-paired-geometry-X-row-order-binding.json').read_text())
assert hashlib.sha256(geo_raw).hexdigest()==binding['geometry_sha256']
phase_root=R3/'all1998-offsets-geometry-freeze-first';phase_meta=json.loads((phase_root/'return.json').read_text());phase=next(x for x in phase_meta['outputs']if x['mode']==mode)
offsets=np.load(phase_root/phase['path'],allow_pickle=False);assert offsets.shape==(999,45)and offsets.dtype==np.dtype('int64')
plans=json.loads((phase_root/'allowed-phase-plans.json').read_text())[mode]
ids=geo['ids'];assert len(ids)==1896 and ids==sorted(set(ids));N=len(ids)
assert binding['oldgroup_ids']==ids and hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()==binding['ids_sha256']
assert all([t['orthogroup']for t in geo['species'][sp]['targets']]==ids for sp in ['G','S'])
old=json.loads((O/'Stage-B-main-comparison-v1/main_50kb-summary.json').read_text());T0=old['T'];assert np.isfinite(T0)
assert old['n_final']==1896 and old['version']=='main_50kb'and old['half_window_bp']==50000
npz_path=R3/'geometry-design-preparation-syntaxfix-once/original-fixed-design-and-saved-residuals.npz'
assert hashlib.sha256(npz_path.read_bytes()).hexdigest()==binding['npz_sha256']
z=np.load(npz_path,allow_pickle=False)
P_fixed=z['P_fixed_residual'];assert P_fixed.shape==(N,)and np.isfinite(P_fixed).all()
Q={sp:fixed_qr(z[sp+'_X'])for sp in ['G','S']}
denominators={sp:np.array([t['D']for t in geo['species'][sp]['targets']],dtype=np.int64)for sp in ['G','S']}
tracks={};column=0
for sp in ['G','S']:
    s=geo['species'][sp];stage=R3/(sp+'-original-annotation-projection-first')
    meta=json.loads((stage/'track-metadata.json').read_text());assert meta['original_zero_phase_numerators_exact']==1896 and meta['species']==sp
    assert [x['chromosome']for x in meta['chromosomes']]==s['chromosomes']
    payload=np.load(stage/'accepted-A-unions.npz',allow_pickle=False);parts=[];seen=set();sum_D=np.zeros(N,dtype=np.int64)
    for cm in meta['chromosomes']:
        ch=cm['chromosome'];L=cm['L'];key=cm['key'];assert phase['columns'][column]==[sp,ch]
        A=AnnotationUnion(payload[key+'_starts'],payload[key+'_ends'],payload[key+'_cumulative'],L)
        assert np.all(A.starts<A.ends)and np.all(A.ends[:-1]<A.starts[1:])and A.cumulative[0]==0 and np.array_equal(np.diff(A.cumulative),A.ends-A.starts)
        targets={t['orthogroup']:t for t in s['targets']if t['seqid']==ch}
        assert not seen.intersection(targets);seen.update(targets)
        expected_D=[targets[g]['D']if g in targets else 0 for g in ids]
        fr=prepare_callable_fragments({g:targets[g]['callable']if g in targets else []for g in ids},L,target_ids=ids,expected_D=expected_D);sum_D+=fr.denominators
        allowed=plans[sp][ch]['allowed'];assert plans[sp][ch]['integer_width']>1
        assert all(any(a<=int(delta)<b for a,b in allowed)for delta in offsets[:,column]),('unreleased_offset',mode,sp,ch)
        parts.append((A,fr,L,column));column+=1
    assert seen==set(ids)and np.array_equal(sum_D,denominators[sp]);tracks[sp]=parts
assert column==45
count_store=np.empty((999,2,N),dtype=np.uint32);records=[];start=time.monotonic()
progress=Path('draw-progress.jsonl').open('x')
for draw in range(999):
    residuals={};diag={}
    progress.write(json.dumps(dict(draw=draw+1,state='STARTED'))+'\n');progress.flush()
    for spi,sp in enumerate(['G','S']):
        counts=np.zeros(N,dtype=np.int64)
        for A,fr,L,col in tracks[sp]:counts+=query_counts(A,fr,int(offsets[draw,col]),L,circular=mode=='circ')
        assert np.all(counts>=0)and np.all(counts<=denominators[sp])and np.max(counts)<=100000
        count_store[draw,spi]=counts
        e=response_residuals(counts,denominators[sp],Q[sp]);assert np.isfinite(e).all()
        residuals[sp]=e;diag[sp+'_residual_norm']=float(np.linalg.norm(e))
    GP=pearson(residuals['G'],P_fixed);GS=pearson(residuals['G'],residuals['S']);PS=pearson(P_fixed,residuals['S']);T=float((GP+GS+PS)/3)
    record=dict(draw=draw+1,r_GP=GP,r_GS=GS,r_PS=PS,T=T,exceeds_observed=int(T>=T0),**diag);records.append(record)
    progress.write(json.dumps(dict(draw=draw+1,state='RETURNED',T=T))+'\n');progress.flush()
progress.close();assert len(records)==999
with Path('all999-conditional-background-statistics.tsv').open('x',newline='')as f:
    w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
np.savez_compressed('all999-two-species-frozen-window-counts.npz',counts=count_store)
null=np.array([r['T']for r in records]);exceed=int(sum(r['exceeds_observed']for r in records));summary=dict(status='ALL999_NEW_TRUE_BP_PHASE_DRAWS_COMPLETE',mode=mode,observed_T_reused=T0,draws=999,new_rank_QR_projections=1998,old_models_refitted=0,per_draw_panel=1896,conditional_nominal_upper_tail=(1+exceed)/1000,exceedances=exceed,T_min=float(null.min()),T_q025=float(np.quantile(null,.025)),T_median=float(np.median(null)),T_q975=float(np.quantile(null,.975)),T_max=float(null.max()),T_mean=float(null.mean()),undefined_draws=0,redraws=0,phase_plan_precedes_tracks=True,source_body_reads=0,Genome_GFF_RM_execution=0,protected_outcomes_opened=0,calibrated_population_CI_claimed=False,exact_group_randomization_claimed=False,interpretation='Conditional annotation-track phase background with geometric protection restrictions; linear is same-panel boundary sensitivity, not independent replication',counts_shape=list(count_store.shape),count_array_sha256=hashlib.sha256(count_store.tobytes()).hexdigest(),seconds=time.monotonic()-start)
Path('return.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
