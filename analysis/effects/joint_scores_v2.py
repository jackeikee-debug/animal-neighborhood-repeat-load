"""Pure NEW conditional score projection/envelope. No old fit or IO calls.

Runtime/source/roles/accepted results/graph/once/resource authority remain Root's.
No backend module is imported. NumPy solve/QR is NEW nuisance projection only.
"""
from copy import deepcopy
import hashlib
import numpy as np

STRUCTURE_NAMES=('nonexonic_GC','gene_count','target_gene_length','exon_fraction','nonACGT_fraction')
LAYERS=('non-Bilateria','Spiralia','Ecdysozoa','Deuterostomia','Xenacoelomorpha_separate')


def need(ok, reason):
    if not ok: raise ValueError('HOLD: '+reason)


def _policy(p):
    P=tuple(p['P_species']);S=tuple(p['original_S35'])
    need(0<len(P)<=35 and len(set(P))==len(P) and len(set(S))==len(S)
         and len(S)<=35 and set(P)<=set(S),'exact unique fixed P/original S')
    na=p['sourceNA_four_original_records']
    need(len(na)==len(S)-len(P) and {r['species'] for r in na.values()}==set(S)-set(P),
         'all original sourceNA identities retained without duplicates')
    need(p['numerical_atol']==1e-10 and p['numerical_rtol']==1e-10
         and p['score_zero_abs_tolerance']==1e-14,'frozen numeric/score-zero tolerances')
    need(p['minimum_nonzero_score_blocks']==20 and p['max_squared_score_share']==.20
         and p['min_effective_information_blocks']==10,'unchanged score gates')
    need(p['U2_joint_sign_vectors']==999 and p['seed']==2026102801
         and p['quantile_method']=='linear' and p['quantile_probability']==.95,
         'frozen one joint perturbation stream/quantile')
    return P,S


def midranks(values):
    a=np.asarray(values,dtype=np.float64)
    need(a.ndim==1 and np.isfinite(a).all(),'finite rank vector')
    order=np.argsort(a,kind='stable');out=np.empty(len(a),dtype=np.float64);lo=0
    while lo<len(a):
        hi=lo+1
        while hi<len(a) and a[order[hi]]==a[order[lo]]:hi+=1
        out[order[lo:hi]]=(lo+1+hi)/2.;lo=hi
    return out


def rank_z(values):
    r=midranks(values)
    if len(r)<2 or np.all(r==r[0]):return None
    return (r-r.mean())/r.std(ddof=1)


def _number(x):
    need(type(x) in (int,float) and np.isfinite(x),'literal finite numeric value')
    return float(x)


def _accepted_summary(record):
    need(isinstance(record,dict) and isinstance(record.get('backend_result'),dict),
         'accepted saved result outer record.backend_result required')
    full=record['backend_result'].get('full')
    need(isinstance(full,dict),'accepted spatial backend_result.full required')
    return full,{k:deepcopy(full[k]) for k in ('status','reason','beta','coefficients',
        'design_columns','dropped_constant_covariates','ids','weights','n') if k in full}


def accepted_full_to_scores(active_input, accepted_full_return, block_map, policy):
    P,S=_policy(policy);a=active_input;record=accepted_full_return
    species=a['species'];assembly=a['assembly']
    need(species in P and a['baseline']=='B0' and a['variant']=='A50','fixed P original B0/A50 only')
    need(record['native_species_name']==species and record['assembly']==assembly and record['variant']=='A50',
         'saved outer result exact species/assembly/variant')
    rows=a['rows'];ids=[r['id'] for r in rows];n=len(ids)
    need(0<n<=1896 and len(set(ids))==n and set(a['observations'])==set(ids), 'exact active old-group input closure')
    need(set(ids)<=set(block_map) and 0<len(block_map)<=1896,'complete fixed global block map')
    need(all(type(b) is str and b for b in block_map.values()),'literal block IDs')
    blocks=sorted(set(block_map.values()))
    w=np.asarray(a['weights'],dtype=np.float64);x=np.asarray(a['x'],dtype=np.float64)
    need(w.shape==x.shape==(n,) and np.isfinite(w).all() and np.isfinite(x).all() and np.all(w>0), 'active weights/x')
    need(np.array_equal(x,np.array([_number(r['xA']) for r in rows])),'exact unchanged frozen xA50')
    full,accepted=_accepted_summary(record)
    out={'schema':'FIXED_B0_A50_ACCEPTED_FULL_BLOCK_SCORES_V1','species':species,'assembly':assembly,
         'status':'TYPED_NA','reason':None,'accepted_full':accepted,'block_order':blocks,
         'block_scores':None,'diagnostics':{},'old_fit_calls':0,'new_nuisance_projection_calls':0,
         'original_S35_retained':True,'policy_id':policy['policy_id']}
    def na(reason):out.update(reason=reason,status='TYPED_NA');return out
    if 'fixed_J_ids' in record:need(record['fixed_J_ids']==ids,'saved outer original fixed J identity')
    if full.get('status')!='OK':return na('ACCEPTED_FULL_ESTIMATE_UNESTIMABLE')
    need(full['ids']==ids and full['n']==n and np.array_equal(np.asarray(full['weights'],dtype=float),w),
         'accepted full exact IDs/order/weights/n; no reconstruction by beta')
    if n<50:return na('ORIGINAL_MINIMUM_50_ANCHORS_NOT_MET')
    beta=_number(full['beta']);coef=np.asarray(full['coefficients'],dtype=np.float64)
    if not np.isfinite(coef).all():return na('NONFINITE_ACCEPTED_COEFFICIENTS')
    locations=[r['location'] for r in rows]
    need(all(type(loc) is str and loc for loc in locations),'original backend-consumable location identities; no decoding')
    loc_order=sorted(set(locations));loc_index={loc:i for i,loc in enumerate(loc_order)}
    li=np.array([loc_index[loc] for loc in locations],dtype=int);L=len(loc_order)
    Qraw=np.asarray([r['structure'] for r in rows],dtype=np.float64)
    need(Qraw.shape==(n,5) and np.isfinite(Qraw).all(),'exact accepted five structure50 values')
    need(np.all((Qraw[:,[0,3,4]]>=0)&(Qraw[:,[0,3,4]]<=1))
         and np.all(Qraw[:,1]>=0) and np.all(Qraw[:,1]==np.floor(Qraw[:,1]))
         and np.all(Qraw[:,2]>0),'original backend structural-value contract')
    columns=[];names=['intercept'];dropped=[]
    for j,name in enumerate(STRUCTURE_NAMES):
        z=rank_z(Qraw[:,j])
        if z is None:dropped.append(name)
        else:columns.append(z);names.append(name)
    Z=np.column_stack(columns) if columns else np.empty((n,0),dtype=float)
    names+=['location:'+loc for loc in loc_order[1:]];names+=['frozen_signature']
    need(full['design_columns']==names and full['dropped_constant_covariates']==dropped
         and coef.shape==(len(names),),'accepted design/drop/coefficients exactly match original algorithm')
    if not np.isclose(coef[-1],beta,atol=1e-10,rtol=1e-10):return na('ACCEPTED_BETA_COEFFICIENT_MISMATCH')
    y_fraction=[];na_ids=[]
    for r in rows:
        o=a['observations'][r['id']]
        if o.get('status')=='TYPED_NA':na_ids.append(r['id']);continue
        A,U,D=o['A'],o['U'],r['nonexonic_acgt']
        need(o.get('status') in (None,'MEASURED') and all(type(v) is int for v in (A,U,D))
             and 0<=U<=A<=D and D>0,'accepted A50 observation bounds')
        y_fraction.append(A/D)
    if na_ids:out['diagnostics']['active_sourceNA_ids']=na_ids;return na('ACTIVE_SOURCE_NA_NOT_DROPPED')
    y=rank_z(y_fraction)
    if y is None:return na('CONSTANT_RESPONSE_NO_SCORE_INTERVAL')
    k=Z.shape[1];pred=np.full(n,coef[0],dtype=float)
    if k:pred+=Z@coef[1:1+k]
    loc_coef=np.zeros(L)
    if L>1:loc_coef[1:]=coef[1+k:-1]
    pred+=loc_coef[li]+x*coef[-1];e=y-pred
    loc_w=np.bincount(li,weights=w,minlength=L)
    def demean(v):
        return v-(np.bincount(li,weights=w*v,minlength=L)/loc_w)[li]
    xd=demean(x);Zd=np.column_stack([demean(Z[:,j]) for j in range(k)]) if k else Z
    sw=np.sqrt(w);A=sw[:,None]*Zd;b=sw*xd
    out['new_nuisance_projection_calls']=1
    if k:
        try:
            q,r=np.linalg.qr(A,mode='reduced');sing=np.linalg.svd(r,compute_uv=False)
        except np.linalg.LinAlgError:
            return na('NEW_NUISANCE_QR_OR_SVD_NUMERIC_FAILURE')
        tol=max(A.shape)*np.finfo(np.float64).eps*sing[0] if len(sing) else 0.
        rank=int(np.sum(sing>tol))
        out['diagnostics'].update(qr_rank=rank,qr_columns=k,
            qr_condition=float(sing[0]/sing[-1]) if sing[-1]>0 else None)
        if rank!=k:return na('NEW_NUISANCE_PROJECTION_RANK_DEFICIENT')
        try:gamma=np.linalg.solve(r,q.T@b)
        except np.linalg.LinAlgError:return na('NEW_NUISANCE_TRIANGULAR_SOLVE_NUMERIC_FAILURE')
        xr=xd-Zd@gamma
        nuisance_q_leverage=np.sum(q*q,axis=1)
    else:
        xr=xd;nuisance_q_leverage=np.zeros(n)
        out['diagnostics'].update(qr_rank=0,qr_columns=0,qr_condition=None)
    denom=float(np.sum(w*xr*xr))
    if not np.isfinite(denom) or denom<=0:return na('NO_REMAINING_SIGNATURE_INFORMATION')
    l=w*xr/denom
    row_x_share=w*xr*xr/denom
    nuisance_leverage=w/loc_w[li]+nuisance_q_leverage
    full_leverage=nuisance_leverage+row_x_share
    out['diagnostics'].update(maximum_nuisance_leverage=float(nuisance_leverage.max()),
        maximum_full_model_leverage=float(full_leverage.max()),
        maximum_single_row_signature_information_share=float(row_x_share.max()),
        signature_information_denominator=denom)
    # No explicit n x chromosome-dummy matrix: evaluate the same columns as
    # independent vector/group operations for identities and normal equations.
    def moments(v):
        pieces=[float(np.sum(v))]
        pieces.extend((Z.T@v).tolist())
        pieces.extend(np.bincount(li,weights=v,minlength=L)[1:].tolist())
        return np.asarray(pieces,dtype=float)
    orth=moments(l);location_orth=np.bincount(li,weights=l,minlength=L)
    rhs=np.r_[moments(w*y),np.dot(x,w*y)]
    lhs=np.r_[moments(w*pred),np.dot(x,w*pred)]
    reproduced=float(np.dot(l,y));normal_error=lhs-rhs
    checks={'sum_l_x':bool(np.isclose(np.dot(l,x),1,atol=1e-10,rtol=1e-10)),
            'nuisance_orthogonality':bool(np.allclose(orth,0,atol=1e-10,rtol=1e-10)),
            'all_location_orthogonality':bool(np.allclose(location_orth,0,atol=1e-10,rtol=1e-10)),
            'accepted_weighted_normal_equations':bool(np.allclose(lhs,rhs,atol=1e-10,rtol=1e-10)),
            'accepted_beta_reproduced':bool(np.isclose(reproduced,beta,atol=1e-10,rtol=1e-10))}
    out['diagnostics'].update(checks=checks,sum_l_x=float(np.dot(l,x)),beta_from_projection_identity=reproduced,
        maximum_nuisance_orthogonality_error=float(np.max(np.abs(orth))),
        maximum_location_orthogonality_error=float(np.max(np.abs(location_orth))),
        maximum_weighted_normal_equation_error=float(np.max(np.abs(normal_error))),
        maximum_normal_equation_scaled_error=float(np.max(np.abs(normal_error)/(1e-10+1e-10*np.abs(rhs)))))
    if not all(checks.values()):return na('ACCEPTED_COEFFICIENT_OR_PROJECTION_IDENTITY_FAILED')
    # A singleton location can have nuisance leverage1 and signature share0.
    # It does not identify beta. A unit full-design leverage row with positive
    # signature share DOES determine signature identification on deletion,
    # including a two-row contrast whose individual shares are about .5.
    unit=full_leverage>=1-1e-10
    signature_critical=unit & (row_x_share>policy['score_zero_abs_tolerance'])
    out['diagnostics'].update(nuisance_only_unit_leverage_rows=int(np.sum(unit & ~signature_critical)),
        signature_critical_unit_leverage_rows=int(np.sum(signature_critical)),
        signature_specific_unit_gate=True)
    if np.any(signature_critical) or np.any(row_x_share>=1-1e-10):
        return na('UNIT_LEVERAGE_OR_SINGLE_ROW_DETERMINATION')
    h={block:0. for block in blocks}
    for g,influence in zip(ids,l*e):h[block_map[g]]+=float(influence)
    if not all(np.isfinite(v) for v in h.values()):return na('NONFINITE_BLOCK_SCORE')
    out.update(status='SCORES_PREPARED',reason=None,block_scores=h,
        diagnostics={**out['diagnostics'],'row_influence_sum':float(np.dot(l,e)),
                     'residual_min':float(e.min()),'residual_max':float(e.max())})
    return out


def _gate(h,policy):
    h=np.asarray(h,dtype=float);count=int(np.sum(np.abs(h)>policy['score_zero_abs_tolerance']))
    with np.errstate(over='ignore',invalid='ignore'):
        sq=h*h;total=float(sq.sum())
    if not np.isfinite(total):
        return {'passed':False,'reasons':['NONFINITE_SCORE_INFORMATION'],
                'nonzero_score_blocks':count,'score_zero_abs_tolerance':1e-14,
                'maximum_squared_score_share':None,'effective_information_blocks':None,
                'sum_squared_block_scores':None}
    maximum=None;effective=None
    if total>0:
        shares=sq/total;maximum=float(shares.max());effective=float(1/np.sum(shares*shares))
    reasons=[]
    if count<20:reasons.append('FEWER_THAN_20_NONZERO_SCORE_BLOCKS')
    if maximum is None or maximum>.20:reasons.append('MAXIMUM_SQUARED_SCORE_SHARE_ABOVE_POINT20_OR_ZERO_INFORMATION')
    if effective is None or effective<10:reasons.append('EFFECTIVE_INFORMATION_BLOCKS_BELOW_10_OR_ZERO_INFORMATION')
    return {'passed':not reasons,'reasons':reasons,'nonzero_score_blocks':count,
            'score_zero_abs_tolerance':1e-14,'maximum_squared_score_share':maximum,
            'effective_information_blocks':effective,'sum_squared_block_scores':total,
            'tolerance_classifies_nonzero_only_raw_scores_not_thresholded':True}


def joint_score_envelope(scores,policy):
    P,S=_policy(policy)
    need(set(scores)==set(P),'every fixed P projection result present; no species removal')
    orders=[scores[s]['block_order'] for s in P]
    need(all(o==orders[0] for o in orders) and 0<len(orders[0])<=1896,'same full global block order for all P')
    blocks=orders[0];need(len(set(blocks))==len(blocks),'unique global blocks')
    H=np.zeros((len(P),len(blocks)));valid=[];betas={};projection_na={}
    for i,s in enumerate(P):
        r=scores[s];need(r['species']==s and r['policy_id']==policy['policy_id'],'exact projection species/policy')
        beta=r['accepted_full'].get('beta');betas[s]=float(beta) if type(beta) in (int,float) and np.isfinite(beta) else None
        if r['status']!='SCORES_PREPARED':projection_na[s]=deepcopy(r);continue
        need(set(r['block_scores'])==set(blocks),'full block vector, not independent per-species blocks')
        H[i]=[r['block_scores'][b] for b in blocks]
        need(np.isfinite(H[i]).all() and betas[s] is not None,'finite accepted beta/scores')
        valid.append(i)
    # Invalid rows are never used as zero-valued estimates/covariances: the
    # zero scratch rows below are excluded explicitly from every public target.
    rng=np.random.Generator(np.random.PCG64(2026102801))
    signs=rng.integers(0,2,size=(999,len(blocks)),dtype=np.int8)*2-1
    delta=signs@H.T
    targets={}
    def target(estimate,h,perturb):
        gate=_gate(h,policy)
        if not gate['passed']:return {'status':'TYPED_NA','reason':'SCORE_INFORMATION_GATE_FAILED',
            'estimate':estimate,'working_range':None,'gate':gate}
        if not np.isfinite(perturb).all():return {'status':'TYPED_NA',
            'reason':'NONFINITE_JOINT_PERTURBATION','estimate':estimate,'working_range':None,'gate':gate}
        q=float(np.quantile(np.abs(perturb),.95,method='linear'))
        return {'status':'WORKING_RANGE_PREPARED','estimate':estimate,'absolute_delta_q95':q,
                'working_range':[estimate-q,estimate+q],'gate':gate}
    for i,s in enumerate(P):
        if s in projection_na:targets[s]={'status':'TYPED_NA','reason':'PROJECTION_NUMERIC_NA',
            'estimate':betas[s],'working_range':None,'projection':deepcopy(projection_na[s])}
        else:targets[s]=target(betas[s],H[i],delta[:,i])
    if projection_na:
        theta={'status':'TYPED_NA','reason':'AT_LEAST_ONE_FIXED_P_PROJECTION_NA_NO_SPECIES_REMOVED',
            'estimate':float(np.mean(list(betas.values()))) if all(v is not None for v in betas.values()) else None,
            'working_range':None,'projection_NA_species':list(projection_na)}
    else:theta=target(float(np.mean([betas[s] for s in P])),H.mean(axis=0),delta.mean(axis=1))
    covariance=H@H.T
    cov=[[float(covariance[i,j]) if i in valid and j in valid and np.isfinite(covariance[i,j]) else None
          for j in range(len(P))] for i in range(len(P))]
    na_by_species={r['species']:r for r in policy['sourceNA_four_original_records'].values()}
    full_rows=[]
    for s in S:
        full_rows.append({'species':s,'in_fixed_observable_P':s in targets,
            'result':deepcopy(targets[s]) if s in targets else {'status':'ORIGINAL_SOURCE_NA_RETAINED',
                'source_record':deepcopy(na_by_species[s]),'estimate':None,'working_range':None}})
    return {'schema':'SYNCHRONIZED_GLOBAL_BLOCK_SCORE_WORKING_ENVELOPE_V1',
        'policy_id':policy['policy_id'],'P_order':list(P),'original_S35_order':list(S),
        'original_S35_rows':full_rows,'original_S35_mean_and_range':{'status':'TYPED_NA',
            'reason':'ORIGINAL_SOURCE_NA_NOT_IMPUTED','estimate':None,'working_range':None},
        'P_conditional_equalweight_mean':theta,'species_targets':targets,
        'sourceNA_original_records':deepcopy(policy['sourceNA_four_original_records']),
        'joint_working_score_covariance':{'species_order':list(P),'matrix':cov},
        'signs':{'generator':'NumPy Generator/PCG64','seed':2026102801,'vectors':999,
                 'block_order':list(blocks),'sign_array_sha256':hashlib.sha256(signs.tobytes()).hexdigest(),
                 'same_block_same_sign_for_all_species':True},
        'interpretation':'95% synchronized dependency-block score perturbation WORKING range; not calibrated CI; no p',
        'node_component_fraction_not_a_score_gate':True,'model_fits':0,'variance_refits':0}


def leave_layer_out_means(accepted_records,policy):
    """U3 only fixed-P arithmetic of existing betas, no intervals or refits."""
    P,S=_policy(policy);need(set(accepted_records)==set(P),'all fixed P accepted records')
    layers=policy['layer_by_species'];need(set(layers)==set(S) and all(v in LAYERS for v in layers.values()),'fixed all-S layer roster')
    beta={}
    for s in P:
        record=accepted_records[s]
        need(record['native_species_name']==s and record['variant']=='A50','fixed accepted record identity')
        full,_=_accepted_summary(record);v=full.get('beta')
        beta[s]=float(v) if full.get('status')=='OK' and type(v) in (int,float) and np.isfinite(v) else None
    base=float(np.mean(list(beta.values()))) if all(v is not None for v in beta.values()) else None
    results=[]
    for layer in LAYERS:
        keep=[s for s in P if layers[s]!=layer]
        bad=[s for s in keep if beta[s] is None]
        mean=float(np.mean([beta[s] for s in keep])) if keep and not bad else None
        results.append({'deleted_layer':layer,'retained_P_species':keep,'n':len(keep),
            'status':'DESCRIPTIVE_ARITHMETIC_ONLY' if mean is not None else 'TYPED_NA',
            'reason':None if mean is not None else 'NO_REMAINING_SPECIES' if not keep else 'RETAINED_ACCEPTED_BETA_NA',
            'mean':mean,'change_from_complete_P_mean':mean-base if mean is not None and base is not None else None,
            'beta_NA_species':bad,'working_range':None})
    source_na={r['species']:r for r in policy['sourceNA_four_original_records'].values()}
    original_rows=[{'species':s,'layer':layers[s],'in_fixed_observable_P':s in beta,
        'accepted_beta':beta.get(s),'sourceNA_record':deepcopy(source_na[s]) if s in source_na else None}
        for s in S]
    return {'schema':'FIXED_P_LEAVE_LAYER_OUT_DESCRIPTIVE_MEANS_V1','original_S35_order':list(S),
        'original_S35_rows':original_rows,
        'P_order':list(P),'sourceNA_original_records':deepcopy(policy['sourceNA_four_original_records']),
        'P_mean':base,'layers':results,'interpretation':'Descriptive arithmetic, not layer-balanced population inference',
        'model_fits':0,'variance_refits':0}
