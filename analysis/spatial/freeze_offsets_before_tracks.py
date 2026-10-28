"""Geometry-only exact integer phase plans and draw freeze; no repeat tracks."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
from physical_shift import offset_plan,sample_offset
R3=Path(__import__('os').environ['ANRL_R3_ROOT'])
geo=json.loads(gzip.decompress((R3/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz').read_bytes()))
modes=['circ','linear'];streams=np.random.SeedSequence(2026102803).spawn(2)
plans={};summary=[];outputs=[]
for mode,seed in zip(modes,streams):
    rng=np.random.Generator(np.random.PCG64(seed));plans[mode]={};offsets=[];keys=[]
    for sp in ['G','S']:
        s=geo['species'][sp];plans[mode][sp]={}
        for ch in s['chromosomes']:
            L=s['lengths'][ch];W=[(t['start0'],t['end0'])for t in s['targets']if t['seqid']==ch];B=s['forbidden'].get(ch,[])
            plan=offset_plan(W,B,L,circular=mode=='circ');assert plan['status']=='ELIGIBLE_CONDITIONAL_PHASE_SET',(mode,sp,ch,plan['status'])
            plans[mode][sp][ch]=plan;keys.append([sp,ch]);offsets.append([sample_offset(plan['allowed'],rng)for _ in range(999)])
            delta=np.array(offsets[-1],dtype=np.int64);distance=np.minimum(delta,L-delta)if mode=='circ'else np.abs(delta)
            summary.append(dict(mode=mode,species=sp,chromosome=ch,L=L,targets=len(W),forbidden_intervals=len(B),allowed_integer_phases=plan['integer_width'],domain=plan['domain'],draws=999,min_distance=int(distance.min()),median_distance=float(np.median(distance)),max_distance=int(distance.max()),fraction_distance_at_least100kb=float(np.mean(distance>=100000)),zero_draws=int(np.sum(delta==0)),duplicates_retained=True))
    array=np.asarray(offsets,dtype=np.int64).T;assert array.shape==(999,45)
    f=Path(mode+'-all999-frozen-offsets.npy');np.save(f,array);outputs.append(dict(mode=mode,path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),columns=keys,spawn_key=list(seed.spawn_key)))
Path('allowed-phase-plans.json').write_text(json.dumps(plans,separators=(',',':'))+'\n')
Path('offset-distance-geometry-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
Path('return.json').write_text(json.dumps(dict(status='GEOMETRY_ONLY_ALL1998_DRAW_PLANS_FROZEN_BEFORE_TRACKS',seed=2026102803,mode_order=modes,outputs=outputs,real_repeat_source_reads=0,real_null_statistics=0,original_observed_refits=0,phase_model_conditional_not_exact_group_test=True),indent=2)+'\n')
