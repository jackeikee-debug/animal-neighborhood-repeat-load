"""Frozen-contract numerical preparation. Pure functions, no files/network/CLI.

Real-data execution requires separate Root acceptance and exact release. None of
these functions grants permission to read outcomes or repeat an execution.
"""
from copy import deepcopy
import math
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

EPS=np.finfo(np.float64).eps
STRUCTURE_NAMES=('nonexonic_GC','gene_count','target_gene_length','exon_fraction','nonACGT_fraction')
LAYERS=('non-Bilateria','Spiralia','Ecdysozoa','Deuterostomia','Xenacoelomorpha_separate')
OPT_OPTIONS=dict(maxiter=1000,maxfun=5000,maxls=20,ftol=1e-12,gtol=1e-8)


def _na(reason,**extra):
    return dict(status='NA',reason=reason,beta=None,SE=None,p_value=None,**extra)


def _vector(value,name):
    a=np.asarray(value,dtype=np.float64)
    if a.ndim!=1 or not np.isfinite(a).all(): raise ValueError(name+' must be finite and one-dimensional')
    return a


def _number(value,name):
    if isinstance(value,(bool,np.bool_)) or not np.isscalar(value): raise ValueError(name+' must be numeric')
    x=float(value)
    if not math.isfinite(x): raise ValueError(name+' must be finite')
    return x


def _integer(value,name):
    x=_number(value,name)
    if x!=math.floor(x): raise ValueError(name+' must be an integer')
    return int(x)


def midranks(values):
    """1-based exact-equality midranks; deterministic stable sorting."""
    a=_vector(values,'values')
    order=np.argsort(a,kind='stable'); out=np.empty(len(a),dtype=float)
    lo=0
    while lo<len(a):
        hi=lo+1
        while hi<len(a) and a[order[hi]]==a[order[lo]]: hi+=1
        out[order[lo:hi]]=(lo+1+hi)/2.0
        lo=hi
    return out


def _rank_z(values):
    r=midranks(values)
    if len(r)<2 or np.all(r==r[0]): return None
    return (r-r.mean())/r.std(ddof=1)


def _rank(X):
    singular=np.linalg.svd(X,compute_uv=False)
    if not len(singular) or singular[0]==0: return 0
    return int(np.sum(singular>max(X.shape)*EPS*singular[0]))


def _lstsq(X,y,w):
    sw=np.sqrt(w); A=X*sw[:,None]
    if _rank(A)!=X.shape[1]: return None
    coef,_,_,_=np.linalg.lstsq(A,y*sw,rcond=max(A.shape)*EPS)
    return coef


def overlap_weights(ids,locations,starts,ends):
    """Each location's components equal-weight; touching is not overlap."""
    n=len(ids)
    if any(len(x)!=n for x in (locations,starts,ends)) or len(set(ids))!=n: raise ValueError('coordinate keys/lengths invalid')
    if any(not isinstance(x,str) or not x for x in ids) or any(not isinstance(x,str) or not x for x in locations): raise ValueError('IDs/locations must be nonempty strings')
    starts=[_integer(x,'start0') for x in starts]; ends=[_integer(x,'end0') for x in ends]
    if any(s<0 or e<=s for s,e in zip(starts,ends)): raise ValueError('invalid half-open interval')
    weights=np.empty(n); components=[None]*n; counts={}
    for loc in sorted(set(locations)):
        ind=sorted((i for i in range(n) if locations[i]==loc),key=lambda i:(starts[i],ends[i],ids[i]))
        groups=[]; right=None
        for i in ind:
            if right is None or starts[i]>=right:
                groups.append([]); right=ends[i]
            else: right=max(right,ends[i])
            groups[-1].append(i)
        counts[loc]=len(groups)
        for j,group in enumerate(groups):
            for i in group:
                weights[i]=1.0/(len(groups)*len(group)); components[i]=(loc,j)
    return dict(weights=weights,components=components,components_per_location=counts)


def _geometry(row,minimum):
    """Only geometry/callability, never repeat outcomes."""
    if not isinstance(row['complete'],(bool,np.bool_)): raise ValueError('complete must be boolean')
    a=_number(row['acgt_fraction'],'acgt_fraction'); d=_integer(row['nonexonic_acgt'],'nonexonic_acgt')
    if not 0<=a<=1 or d<0: raise ValueError('invalid callability values')
    reasons=[]
    if not row['complete']: reasons.append('INCOMPLETE_WINDOW')
    if a<.95: reasons.append('INSUFFICIENT_ACGT')
    if d<minimum: reasons.append('INSUFFICIENT_NONEXONIC_ACGT')
    return reasons


def _structure(row):
    q=_vector(row['structure'],'structure')
    if len(q)!=5: raise ValueError('five structural values required')
    if any(not 0<=q[j]<=1 for j in (0,3,4)) or q[1]<0 or q[1]!=math.floor(q[1]) or q[2]<=0: raise ValueError('invalid structural values')
    return q


def _design(rows):
    n=len(rows); Q=np.asarray([_structure(r) for r in rows],dtype=float).reshape(n,5)
    columns=[np.ones(n)]; names=['intercept']; dropped=[]
    for j,name in enumerate(STRUCTURE_NAMES):
        z=_rank_z(Q[:,j])
        if z is None: dropped.append(name)
        else: columns.append(z); names.append(name)
    locations=sorted(set(r['location'] for r in rows))
    for loc in locations[1:]:
        columns.append(np.array([float(r['location']==loc) for r in rows])); names.append('location:'+loc)
    return np.column_stack(columns),names,dropped


def _information(rows,w,x):
    n=len(rows)
    if not n: return dict(status='NA',reason='EMPTY_PANEL',spatial_eligible=False)
    Z,names,dropped=_design(rows)
    X=np.column_stack((Z,x))
    if _rank(Z*np.sqrt(w)[:,None])!=Z.shape[1] or _rank(X*np.sqrt(w)[:,None])!=X.shape[1]:
        return dict(status='NA',reason='RANK_DEFICIENT_OR_NO_X_VARIATION',spatial_eligible=False,
                    dropped_constant_covariates=dropped,design_columns=names)
    coef=_lstsq(Z,x,w); xr=x-Z@coef
    amounts={loc:float(np.sum(w[np.array([r['location']==loc for r in rows])]*xr[np.array([r['location']==loc for r in rows])]**2)) for loc in sorted(set(r['location'] for r in rows))}
    total=sum(amounts.values())
    if not total>0: return dict(status='NA',reason='NO_X_INFORMATION',spatial_eligible=False)
    fractions={k:v/total for k,v in amounts.items()}; eff=1/sum(v*v for v in fractions.values()); maximum=max(fractions.values())
    reasons=[]
    if n<200: reasons.append('FEWER_THAN_200_ANCHORS')
    if len(amounts)<5: reasons.append('FEWER_THAN_5_REAL_CHROMOSOMES')
    if eff<4: reasons.append('EFFECTIVE_INFORMATION_CHROMOSOMES_BELOW_4')
    if maximum>.5: reasons.append('ONE_CHROMOSOME_INFORMATION_ABOVE_HALF')
    return dict(status='OK',spatial_eligible=not reasons,spatial_exclusions=reasons,n=n,real_chromosomes=len(amounts),
                information_fractions=fractions,effective_chromosomes=eff,maximum_fraction=maximum,
                dropped_constant_covariates=dropped,design_columns=names)


def prepare_primary(rows):
    """Input mapping/coordinate bridge must already be Root-accepted.

    Each row has id,location,real_chromosome,start0,end0,complete,
    acgt_fraction,nonexonic_acgt,structure[5],xA,xK. Never reads outcomes.
    Unplaced qualifying measurements are retained as local descriptive IDs,
    not pooled into a real-chromosome coefficient or independence count.
    """
    rows=deepcopy(list(rows)); ids=[r['id'] for r in rows]
    if len(set(ids))!=len(ids) or any(not isinstance(x,str) or not x for x in ids): raise ValueError('duplicate/invalid anchor IDs')
    eligible=[]; local=[]; excluded=[]; location_roles={}
    for r in rows:
        if not isinstance(r['location'],str) or not r['location']: raise ValueError('invalid location')
        if not isinstance(r['real_chromosome'],(bool,np.bool_)): raise ValueError('real_chromosome must be boolean')
        if r['location'] in location_roles and location_roles[r['location']]!=r['real_chromosome']: raise ValueError('location role conflict')
        location_roles[r['location']]=r['real_chromosome']
        _structure(r); _number(r['xA'],'xA'); _number(r['xK'],'xK')
        reasons=_geometry(r,40000)
        if reasons: excluded.append(dict(id=r['id'],reasons=reasons)); continue
        s=_integer(r['start0'],'start0'); e=_integer(r['end0'],'end0')
        if e-s!=100000 or s<0 or r['nonexonic_acgt']>100000: raise ValueError('complete A50 window geometry invalid')
        if not r['real_chromosome']: local.append(r['id']); continue
        eligible.append(r)
    ow=overlap_weights([r['id'] for r in eligible],[r['location'] for r in eligible],
                       [r['start0'] for r in eligible],[r['end0'] for r in eligible])
    info={v:_information(eligible,ow['weights'],np.array([r['xK' if v=='K50' else 'xA'] for r in eligible],dtype=float)) for v in ('A50','K50')}
    return dict(status='OK',ids=[r['id'] for r in eligible],rows=eligible,weights=ow['weights'],components=ow['components'],
                components_per_location=ow['components_per_location'],information=info,excluded=excluded,local_descriptive_ids=local)


def _variant(panel,variant,narrow_rows=None):
    if variant not in ('A50','K50','A25'): raise ValueError('unknown finite version')
    rows=deepcopy(panel['rows']); w=np.array(panel['weights'],dtype=float,copy=True); exclusions=[]
    if variant=='A25':
        if narrow_rows is None or set(narrow_rows)!=set(panel['ids']): raise ValueError('A25 metadata must contain exactly primary J')
        retained=[]; wi=[]
        for i,r in enumerate(rows):
            nr=deepcopy(narrow_rows[r['id']]); reasons=_geometry(nr,20000); _structure(nr)
            if nr['nonexonic_acgt']>50000: raise ValueError('A25 denominator exceeds window')
            if reasons: exclusions.append(dict(id=r['id'],reasons=reasons)); continue
            r.update(complete=nr['complete'],acgt_fraction=nr['acgt_fraction'],nonexonic_acgt=nr['nonexonic_acgt'],structure=nr['structure'])
            retained.append(r); wi.append(w[i])
        rows=retained; w=np.array(wi,dtype=float)
    x=np.array([r['xK' if variant=='K50' else 'xA'] for r in rows],dtype=float)
    return rows,w,x,exclusions


def _fit_rows(rows,w,x,observations,variant):
    n=len(rows)
    if n<50: return _na('INSUFFICIENT_ANCHORS',n=n)
    coverage=[]
    for r in rows:
        o=observations[r['id']]; A=_integer(o['A'],'A'); U=_integer(o['U'],'U'); D=_integer(r['nonexonic_acgt'],'D')
        if not 0<=U<=A<=D: raise ValueError('outcome bounds require 0<=U<=A<=D')
        coverage.append((A-U if variant=='K50' else A)/D)
    y=_rank_z(coverage)
    if y is None: return _na('CONSTANT_RESPONSE',n=n,constant_level=float(coverage[0]))
    Z,names,dropped=_design(rows); X=np.column_stack((Z,x))
    if _rank(X*np.sqrt(w)[:,None])!=X.shape[1]:
        return _na('RANK_DEFICIENT_OR_NO_X_VARIATION',n=n,dropped_constant_covariates=dropped)
    coef=_lstsq(X,y,w)
    if coef is None: return _na('RANK_DEFICIENT_OR_NO_X_VARIATION',n=n)
    return dict(status='OK',beta=float(coef[-1]),n=n,ids=[r['id'] for r in rows],weights=w.tolist(),
                dropped_constant_covariates=dropped,design_columns=names+['frozen_signature'],
                coefficients=coef.tolist(),response_unit='new_species_midrank_sample_SD',SE=None,p_value=None)


def fit_species(panel,observations,variant='A50',narrow_rows=None):
    if set(observations)!=set(panel['ids']): raise ValueError('outcomes must contain exactly fixed primary J; no silent selection')
    rows,w,x,excluded=_variant(panel,variant,narrow_rows)
    result=_fit_rows(rows,w,x,observations,variant); result.update(variant=variant,excluded=excluded)
    return result


def spatial_jackknife(panel,observations,variant='A50',narrow_rows=None):
    if set(observations)!=set(panel['ids']): raise ValueError('outcomes must contain exactly fixed primary J')
    rows,w,x,excluded=_variant(panel,variant,narrow_rows)
    info=_information(rows,w,x)
    full=_fit_rows(rows,w,x,observations,variant)
    if not info.get('spatial_eligible',False):
        return dict(status='NA',reason='SPATIAL_LAYER_INELIGIBLE',vJ=None,full=full,information=info,deletions=[],excluded=excluded)
    if full['status']!='OK':
        return dict(status='NA',reason='FULL_ESTIMATE_UNESTIMABLE',vJ=None,full=full,information=info,deletions=[],excluded=excluded)
    locations=sorted(set(r['location'] for r in rows)); deletions=[]
    for loc in locations:
        keep=np.array([r['location']!=loc for r in rows])
        result=_fit_rows([r for i,r in enumerate(rows) if keep[i]],w[keep],x[keep],observations,variant)
        deletions.append(dict(deleted_chromosome=loc,**result))
    if any(r['status']!='OK' for r in deletions):
        return dict(status='NA',reason='AT_LEAST_ONE_DELETION_UNESTIMABLE',vJ=None,full=full,
                    information=info,deletions=deletions,excluded=excluded)
    betas=np.array([r['beta'] for r in deletions]); L=len(betas)
    variance=float((L-1)/L*np.sum((betas-betas.mean())**2))
    return dict(status='OK',vJ=variance,beta=full['beta'],full=full,information=info,deletions=deletions,
                deletion_min=float(betas.min()),deletion_max=float(betas.max()),
                sign_changes=int(np.sum(np.sign(betas)!=np.sign(full['beta']))),excluded=excluded,
                interpretation='Working chromosome deletion error, not validated biological sampling variance')


def _validate_model(beta,C,v):
    beta=_vector(beta,'beta'); v=_vector(v,'vJ'); C=np.asarray(C,dtype=float)
    n=len(beta)
    if n<3 or len(v)!=n or C.shape!=(n,n) or not np.isfinite(C).all() or np.any(v<0): raise ValueError('model input shape/variance invalid')
    if not np.allclose(C,C.T,rtol=0,atol=1e-12) or not np.allclose(np.diag(C),1,rtol=0,atol=1e-12): raise ValueError('C must be symmetric with diagonal1')
    tol=n*EPS*max(1,float(np.linalg.norm(C,2)))
    if np.linalg.eigvalsh(C).min() < -tol: raise ValueError('C not positive semidefinite')
    return beta,C,v


def reml_eval(theta,beta,C,v):
    """Analytic restricted objective, gradient and observed Hessian.

    This evaluator intentionally raises at singular V, without ridge repair.
    It is also the arithmetic-fixture derivative inspection API.
    """
    beta,C,v=_validate_model(beta,C,v); theta=_vector(theta,'theta')
    if len(theta)!=2 or np.any(theta<0): raise ValueError('two nonnegative components required')
    n=len(beta); I=np.eye(n); V=theta[0]*C+theta[1]*I+np.diag(v)
    try: cf=cho_factor(V,lower=True,check_finite=False)
    except np.linalg.LinAlgError as e: raise ValueError('V is not strictly positive definite') from e
    inv=cho_solve(cf,I,check_finite=False); a=inv@np.ones(n); info=float(a.sum())
    if not info>0: raise ValueError('nonpositive GLS information')
    mu=float(a@beta/info); P=inv-np.outer(a,a)/info; P=(P+P.T)/2
    r=beta-mu; py=inv@r
    f=.5*((n-1)*math.log(2*math.pi)+2*np.log(np.diag(cf[0])).sum()+math.log(info)+float(r@py))
    D=(C,I); grad=np.array([.5*(np.trace(P@d)-float(py@d@py)) for d in D])
    H=np.empty((2,2))
    for i in range(2):
        for j in range(2): H[i,j]=.5*(-np.trace(P@D[i]@P@D[j])+2*float(py@D[i]@P@D[j]@py))
    H=(H+H.T)/2
    return dict(objective=float(f),gradient=grad,hessian=H,mu=mu,conditional_se=math.sqrt(1/info),V=V)


def _terminal(theta,evaluation):
    theta=np.asarray(theta); g=evaluation['gradient']; H=evaluation['hessian']; active=(theta==0)
    projected=np.where(active,np.minimum(g,0),g)
    critical=np.flatnonzero((~active)|(g<=0))
    diagnostics=dict(active_bounds=active.tolist(),gradient=g.tolist(),projected_KKT_gradient=projected.tolist(),
                     hessian=H.tolist(),critical_coordinates=critical.tolist())
    if not len(critical):
        diagnostics.update(hessian_status='NOT_APPLICABLE_STRICT_BOUNDARY',scaled_gradient_max=0.0,accepted=True)
        return diagnostics
    Hc=H[np.ix_(critical,critical)]
    try:
        L=np.linalg.cholesky(Hc)
        scaled=np.linalg.solve(L,projected[critical]); maximum=float(np.max(np.abs(scaled)))
    except np.linalg.LinAlgError:
        diagnostics.update(hessian_status='NOT_POSITIVE_DEFINITE_ON_CRITICAL_SUBSPACE',scaled_gradient_max=None,accepted=False)
        return diagnostics
    diagnostics.update(hessian_status='POSITIVE_DEFINITE_ON_CRITICAL_SUBSPACE',scaled_gradient_max=maximum,
                       accepted=bool(math.isfinite(maximum) and maximum<=.002))
    return diagnostics


def fit_reml(beta,C,v):
    """Numerical core only; caller must enforce frozen species-role/S gates.

    Exactly one optimizer call, one fixed start. No p or CI is produced.
    Public fit_working_summary below enforces the scientific panel gates.
    """
    beta,C,v=_validate_model(beta,C,v); n=len(beta)
    P0=np.eye(n)-np.ones((n,n))/n
    basis=np.column_stack(((P0@C@P0).ravel(),P0.ravel()))
    if _rank(basis)!=2: return dict(status='NA',reason='VARIANCE_COMPONENTS_NOT_IDENTIFIABLE',optimizer_calls=0)
    scale=float(beta.var(ddof=1)+v.mean())
    if not scale>0: return dict(status='NA',reason='ALL_ZERO_VARIANCE_SCALE',optimizer_calls=0)
    offset=float(beta.mean()); y=(beta-offset)/math.sqrt(scale); vv=v/scale
    invalid_domain_evaluations=0
    def objective(phi):
        nonlocal invalid_domain_evaluations
        try:
            e=reml_eval(phi,y,C,vv); return e['objective'],e['gradient']
        except ValueError:
            # An infeasible line-search point is not repaired or retried.
            # NaN derivative means optimizer cannot claim valid convergence there.
            invalid_domain_evaluations+=1
            return np.inf,np.full(2,np.nan)
    result=minimize(objective,np.array([.5,.5]),jac=True,method='L-BFGS-B',bounds=((0,None),(0,None)),options=OPT_OPTIONS)
    base=dict(optimizer_calls=1,optimizer_success=bool(result.success),optimizer_status=int(result.status),
              optimizer_message=str(result.message),iterations=int(result.nit),function_evaluations=int(result.nfev),
              scale=scale,normalized_start=[.5,.5],terminal_phi=np.asarray(result.x).tolist(),
              invalid_domain_evaluations=invalid_domain_evaluations)
    if not result.success or not np.isfinite(result.x).all():
        return dict(status='NA',reason='SINGLE_OPTIMIZER_DID_NOT_CONVERGE',**base)
    try: final=reml_eval(result.x,y,C,vv)
    except ValueError: return dict(status='NA',reason='SINGULAR_TERMINAL_V',**base)
    diagnostics=_terminal(result.x,final)
    if not diagnostics['accepted']: return dict(status='NA',reason='TERMINAL_KKT_OR_HESSIAN_GATE_FAILED',diagnostics=diagnostics,**base)
    values=np.array([final['objective'],final['mu'],final['conditional_se'],*final['gradient'],*final['hessian'].ravel()])
    if not np.isfinite(values).all(): return dict(status='NA',reason='NONFINITE_TERMINAL_RESULT',**base)
    return dict(status='OK',mu_work=offset+math.sqrt(scale)*final['mu'],conditional_se=math.sqrt(scale)*final['conditional_se'],
                tau2=(scale*np.asarray(result.x)).tolist(),V=(scale*final['V']).tolist(),
                objective=final['objective']+.5*(n-1)*math.log(scale),diagnostics=diagnostics,
                report_limit='Conditional working SE only; shared errors uncalibrated; no p or calibrated CI',**base)


def fit_working_summary(records,C):
    """records must be exact outcome-blind final S, never silently filtered."""
    records=list(records)
    if len(set(r['species'] for r in records))!=len(records): raise ValueError('duplicate species')
    if any(r['role']!='new_normal_metazoan' or not r['spatial_eligible'] for r in records): raise ValueError('records not exact accepted new-animal spatial S')
    if any(r['layer'] not in LAYERS for r in records): raise ValueError('unknown fixed layer')
    counts={layer:sum(r['layer']==layer for r in records) for layer in LAYERS}
    if len(records)<20 or sum(counts[x]>=3 for x in LAYERS[:4])<3:
        return dict(status='NA',reason='INSUFFICIENT_SPECIES_OR_MAJOR_LAYER_COVERAGE',layer_counts=counts,optimizer_calls=0)
    result=fit_reml([r['beta'] for r in records],C,[r['vJ'] for r in records])
    result.update(species_order=[r['species'] for r in records],layer_counts=counts)
    return result


def fixed_v_influence(beta,V,layers):
    """Use only a previously accepted saved V, no variance refit."""
    beta=_vector(beta,'beta'); V=np.asarray(V,dtype=float); layers=list(layers); n=len(beta)
    if V.shape!=(n,n) or len(layers)!=n or not np.isfinite(V).all() or not np.allclose(V,V.T,rtol=0,atol=1e-12): raise ValueError('fixed V/labels invalid')
    if any(x not in LAYERS for x in layers): raise ValueError('unknown fixed layer')
    output=[]
    for layer in LAYERS:
        if layer not in layers: continue
        keep=np.array([x!=layer for x in layers])
        if not keep.any(): output.append(dict(deleted_layer=layer,status='NA',reason='NO_REMAINING_SPECIES',mu_fixedV=None)); continue
        A=V[np.ix_(keep,keep)]
        try: cf=cho_factor(A,lower=True,check_finite=False)
        except np.linalg.LinAlgError:
            output.append(dict(deleted_layer=layer,status='NA',reason='SINGULAR_FIXED_SUBMATRIX',mu_fixedV=None)); continue
        ones=np.ones(int(keep.sum())); a=cho_solve(cf,ones,check_finite=False); info=float(a.sum())
        if info<=0: output.append(dict(deleted_layer=layer,status='NA',reason='NONPOSITIVE_INFORMATION',mu_fixedV=None)); continue
        output.append(dict(deleted_layer=layer,status='OK',mu_fixedV=float(a@beta[keep]/info),remaining_species=int(keep.sum())))
    return output
