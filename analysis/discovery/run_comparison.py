"""One fixed main comparison and S1, from accepted compact Stage-A inputs only.

No API, source rescanning, new candidates, optional regressions, or retries.
Output creation is exclusive; a claim prevents an accidental second analysis.
"""
import csv
import gzip
import hashlib
import json
import os
import resource
import signal
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
STAGE = OUT.parent
SPECIES = ('G', 'P', 'S')
PAIRS = ((0, 1), (0, 2), (1, 2))
PAIR_NAMES = ('GP', 'GS', 'PS')
ASSEMBLIES = dict(G='GCF_964340395.1', P='GCF_041260155.1', S='GCF_978046115.1')
COVS = ('nonexonic_GC_fraction', 'overlapping_gene_count', 'target_gene_length_bp', 'exon_fraction', 'nonACGT_fraction')
COUNTS = dict(fits=0, shifts=0)


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(name, value):
    with (OUT/name).open('x') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def table(name, rows):
    if not rows:
        raise ValueError('empty_output_table:' + name)
    path = OUT/name
    if name.endswith('.gz'):
        f = gzip.open(path, 'xt', newline='')
    else:
        f = path.open('x', newline='')
    with f:
        w = csv.DictWriter(f, list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def average_rank(x):
    return pd.Series(np.asarray(x, dtype=float)).rank(method='average').to_numpy()


def pair_correlations(x):
    c = np.corrcoef(x, rowvar=False)
    result = np.array([c[i,j] for i,j in PAIRS])
    if not np.isfinite(result).all():
        raise ValueError('undefined_correlation')
    return result


def trim_singletons(chrom, ids):
    keep = np.ones(len(ids), dtype=bool)
    exclusions = []
    iteration = 0
    while True:
        iteration += 1
        counts = {s: Counter(chrom[keep,j]) for j,s in enumerate(SPECIES)}
        remove = []
        for i in np.flatnonzero(keep):
            reasons = [s+':'+str(chrom[i,j]) for j,s in enumerate(SPECIES) if counts[s][chrom[i,j]]==1]
            if reasons:
                remove.append(i)
                exclusions.append(dict(orthogroup=ids[i], iteration=iteration, reason=';'.join(reasons)))
        if not remove:
            return keep, exclusions
        keep[remove] = False


def fit_one(y, covariates, chrom, version, sp):
    cols, names, dropped, transforms = [np.ones(len(y))], ['intercept'], [], {}
    for name in COVS:
        z = average_rank(covariates[name])
        sd = float(np.std(z))
        if sd == 0:
            dropped.append(name)
            continue
        mean = float(np.mean(z))
        cols.append((z-mean)/sd)
        names.append('rank_z:'+name)
        transforms[name] = dict(rank_mean=mean, rank_sd_ddof0=sd)
    levels = sorted(set(chrom))
    for level in levels[1:]:
        cols.append((chrom==level).astype(float))
        names.append('chromosome:'+str(level))
    X = np.column_stack(cols)
    # Each invocation is counted before numerical work. One SVD solve, no refit.
    COUNTS['fits'] += 1
    beta, _, rank, singular = np.linalg.lstsq(X, y, rcond=None)
    residual = y-X@beta
    finite = bool(np.isfinite(beta).all() and np.isfinite(residual).all())
    orthogonality = float(np.linalg.norm(X.T@residual)/(np.linalg.norm(X)*np.linalg.norm(y)))
    checks = dict(finite=finite, design_full_rank=bool(rank==X.shape[1]),
                  relative_normal_equation_error=orthogonality,
                  residual_positive_variance=bool(np.std(residual)>0))
    fit = dict(version=version, species=sp, n=len(y), columns=X.shape[1], rank=int(rank),
        response='average_rank_of_annotated_interspersed_bp_over_nonexonic_ACGT_bp',
        continuous_rank_normalization='centered and population-SD scaled; same column space',
        chromosome_reference=levels[0], chromosome_levels=levels, dropped_constant_covariates=dropped,
        covariate_transforms=transforms, coefficients=dict(zip(names,map(float,beta))),
        singular_values=singular.tolist(), condition_number=float(singular[0]/singular[-1]),
        residual_mean=float(np.mean(residual)), residual_SD=float(np.std(residual)),
        rank_response_R2=float(1-np.sum(residual**2)/np.sum((y-y.mean())**2)), checks=checks)
    save(version+'-'+sp+'-fixed-adjustment.json', fit)
    if not (finite and rank==X.shape[1] and orthogonality<1e-10 and np.std(residual)>0):
        raise ValueError('unidentifiable_or_invalid_fixed_adjustment:'+version+':'+sp)
    return residual, fit


def describe(v):
    q = np.quantile(v, [0, .25, .5, .75, 1], method='linear')
    return dict(n=len(v), mean=float(np.mean(v)), minimum=float(q[0]), q25=float(q[1]), median=float(q[2]), q75=float(q[3]), maximum=float(q[4]))


def run_version(frame, half):
    version = 'main_50kb' if half==50000 else 'S1_25kb'
    all_ids = sorted(frame.orthogroup.unique())
    sub = frame[frame.half_window_bp==half]
    good_ids = sorted(og for og,g in sub.groupby('orthogroup') if len(g)==3 and (g.triad_QC_pass==1).all())
    qc_dropped = sorted(set(all_ids)-set(good_ids))
    aligned = [sub[sub.species==s].set_index('orthogroup').loc[good_ids] for s in SPECIES]
    chrom0 = np.column_stack([g.seqid.to_numpy() for g in aligned])
    keep, singletons = trim_singletons(chrom0,good_ids)
    ids = np.array(good_ids)[keep].tolist()
    aligned = [g.loc[ids].copy() for g in aligned]
    chrom = np.column_stack([g.seqid.to_numpy() for g in aligned])
    nchr = {s:len(set(chrom[:,j])) for j,s in enumerate(SPECIES)}
    if len(ids)<500 or min(nchr.values())<10:
        raise ValueError('predefined_panel_minimum_failed:'+version)
    if any(min(Counter(chrom[:,j]).values())<2 for j in range(3)):
        raise ValueError('singleton_remains')
    save(version+'-panel.json', dict(n_input_groups=len(all_ids),n_QC_groups=len(good_ids),
        QC_excluded=qc_dropped, singleton_exclusions=singletons, final_groups=ids, n_final=len(ids),chromosomes=nchr))
    load, yrank, residuals, fits, descriptions = [], [], [], [], {}
    for sp,g in zip(SPECIES,aligned):
        denom = g.nonexonic_ACGT_bp.to_numpy(dtype=float)
        values = g.annotated_interspersed_bp.to_numpy(dtype=float)/denom
        cov = dict(nonexonic_GC_fraction=g.nonexonic_GC_bp.to_numpy(dtype=float)/denom,
            overlapping_gene_count=g.overlapping_gene_count.to_numpy(dtype=float),
            target_gene_length_bp=g.target_gene_length_bp.to_numpy(dtype=float),
            exon_fraction=g.all_isoform_exon_union_bp.to_numpy(dtype=float)/g.window_bp.to_numpy(dtype=float),
            nonACGT_fraction=g.nonACGT_bp.to_numpy(dtype=float)/g.window_bp.to_numpy(dtype=float))
        if not np.isfinite(np.column_stack([values]+list(cov.values()))).all():
            raise ValueError('nonfinite_measurement:'+version+':'+sp)
        y = average_rank(values)
        residual, fit = fit_one(y,cov,g.seqid.to_numpy(),version,sp)
        load.append(values);yrank.append(y);residuals.append(residual);fits.append(fit)
        descriptions[sp] = dict(repeat_load=describe(values),
            classified_TE_fraction=describe(g.classified_TE_bp.to_numpy(dtype=float)/denom),
            Unknown_fraction=describe(g.Unknown_bp.to_numpy(dtype=float)/denom),
            broad_class_conflict_fraction=describe(g.broad_class_conflict_bp.to_numpy(dtype=float)/denom),
            covariates={k:describe(v) for k,v in cov.items()})
    raw = np.column_stack(load)
    ranks = np.column_stack(yrank)
    res = np.column_stack(residuals)
    raw_r = pair_correlations(ranks)
    r = pair_correlations(res)
    T = float(np.mean(r))
    saved = []
    for i,og in enumerate(ids):
        row = dict(orthogroup=og)
        for j,s in enumerate(SPECIES):
            row.update({s+'_seqid':str(chrom[i,j]),s+'_tss0':int(aligned[j].iloc[i].tss0),
                s+'_repeat_load':float(raw[i,j]),s+'_response_average_rank':float(ranks[i,j]),s+'_residual':float(res[i,j])})
        saved.append(row)
    table(version+'-saved-values.tsv.gz',saved)
    # Freeze sort order and chromosome membership once. No refits in background.
    blocks = []
    for j in (0,2):
        for ch in sorted(set(chrom[:,j])):
            idx = np.flatnonzero(chrom[:,j]==ch)
            idx = np.array(sorted(idx,key=lambda i:(int(aligned[j].iloc[i].tss0),ids[i])))
            blocks.append((j,str(ch),idx))
    rng = np.random.Generator(np.random.PCG64(20261013))
    null, shift_log = [], []
    for b in range(999):
        shifted = res.copy()
        for j,ch,idx in blocks:
            k = int(rng.integers(0,len(idx)))
            shifted[idx,j] = np.roll(res[idx,j],k)
            shift_log.append(dict(draw=b+1,species=SPECIES[j],chromosome=ch,n=len(idx),offset=k))
        rb = pair_correlations(shifted)
        COUNTS['shifts'] += 1
        null.append(dict(draw=b+1,**dict(zip(PAIR_NAMES,map(float,rb))),T=float(np.mean(rb))))
    table(version+'-null.tsv.gz',null)
    table(version+'-cyclic-shifts.tsv.gz',shift_log)
    null_t = np.array([row['T'] for row in null])
    exceed = int(np.sum(null_t>=T))
    influence=[]
    for j,s in enumerate(SPECIES):
        for ch in sorted(set(chrom[:,j])):
            use=chrom[:,j]!=ch
            rr=pair_correlations(res[use])
            influence.append(dict(deleted_species=s,deleted_chromosome=str(ch),n_removed=int(np.sum(~use)),
                n_retained=int(np.sum(use)),**dict(zip(PAIR_NAMES,map(float,rr))),T=float(np.mean(rr)),
                delta_T=float(np.mean(rr)-T),all_pairs_positive=bool(np.all(rr>0))))
    table(version+'-fixed-residual-chromosome-influence.tsv',influence)
    summary=dict(version=version,half_window_bp=half,n_QC=len(good_ids),n_final=len(ids),
        n_singleton_excluded=len(singletons),chromosomes=nchr,
        raw_Spearman=dict(zip(PAIR_NAMES,map(float,raw_r))),residual_Pearson=dict(zip(PAIR_NAMES,map(float,r))),
        T=T,null_T_q025=float(np.quantile(null_t,.025)),null_T_q975=float(np.quantile(null_t,.975)),
        null_T_median=float(np.median(null_t)),null_exceedances=exceed,empirical_upper_tail_p=(1+exceed)/1000,
        all_pairs_positive=bool(np.all(r>0)),nominal_primary_rule_pass=bool(np.all(r>0) and (1+exceed)/1000<=.05),
        fixed_residual_chromosome_influence=dict(n_deletions=len(influence),T_min=min(x['T'] for x in influence),
            T_max=max(x['T'] for x in influence),all_deletions_all_pairs_positive=all(x['all_pairs_positive'] for x in influence),
            pair_min={p:min(x[p] for x in influence) for p in PAIR_NAMES},pair_max={p:max(x[p] for x in influence) for p in PAIR_NAMES},
            largest_absolute_delta=max(influence,key=lambda x:abs(x['delta_T']))),
        numerical_fits=[{k:f[k] for k in ['species','rank','columns','condition_number','rank_response_R2','checks','dropped_constant_covariates']} for f in fits],
        annotation_distributions=descriptions,seed=20261013,draws=999,
        residual_correlation_type='Pearson on fixed rank-regression residuals, not reranked',
        background_caveat='approximate cyclic exchangeability; unequal physical distances not preserved',
        influence_caveat='fixed-residual descriptive diagnostic, no refit and not confidence interval')
    save(version+'-summary.json',summary)
    return summary


def main():
    scope=json.loads((OUT/'execution-scope.json').read_text())
    for entry in scope['inputs']:
        if digest(entry['path'])!=entry['sha256']:
            raise ValueError('input_identity_changed:'+entry['path'])
    with (OUT/'analysis-once-claim.json').open('x') as f:
        json.dump(dict(time_utc=now(),scope_sha256=digest(OUT/'execution-scope.json'),code_sha256=digest(__file__),
                       automatic_retry=False,max_fits=6,max_shifts=1998),f)
    start=time.monotonic()
    try:
        signal.alarm(300)
        frame=pd.read_csv(STAGE/'all-1896-triplets-two-scales-input-statistics.tsv.gz',sep='\t')
        if len(frame)!=11376 or frame.duplicated(['orthogroup','species','half_window_bp']).any():
            raise ValueError('input_key_or_row_count')
        for sp,assembly in ASSEMBLIES.items():
            if set(frame.loc[frame.species==sp,'assembly'])!={assembly}:
                raise ValueError('assembly_mismatch:'+sp)
        summaries=[run_version(frame,half) for half in (50000,25000)]
        save('comparison-summary.json',dict(time_utc=now(),status='TWO_VERSIONS_COMPLETE_STOP',versions=summaries,
            S2=dict(status='IDENTICAL_TO_MAIN;CORE_CLUSTER_INTERSECTIONS_ZERO',additional_fits=0,additional_shifts=0,
                independent_replication=False,all_Hox_related_loci_excluded=False,post_Stage_A_metadata_reconciliation=True),
            actual_counts=COUNTS,main_plus_S1_directions_consistent=all(s['all_pairs_positive'] for s in summaries),
            new_Galaxy_jobs=0,Stage_C_started=False,manuscript_v6_modified=False))
        size=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
        if size>scope['output_cap_B']:
            raise ValueError('compact_output_budget_exceeded')
        save('analysis-terminal.json',dict(time_utc=now(),status='COMPLETE_STOP',actual_counts=COUNTS,
            wall_s=time.monotonic()-start,maxrss_platform_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            platform=sys.platform,output_B_at_terminal=size,automatic_retry=False))
        print(json.dumps({'status':'COMPLETE_STOP','counts':COUNTS,'versions':[
            {k:s[k] for k in ['version','n_final','raw_Spearman','residual_Pearson','T','empirical_upper_tail_p','fixed_residual_chromosome_influence']}
            for s in summaries]},ensure_ascii=False))
    except BaseException as e:
        save('analysis-terminal.json',dict(time_utc=now(),status='HOLD_NO_RETRY',actual_counts=COUNTS,
            wall_s=time.monotonic()-start,error=type(e).__name__+':'+str(e),automatic_retry=False))
        raise
    finally:
        signal.alarm(0)


if __name__=='__main__':
    main()
