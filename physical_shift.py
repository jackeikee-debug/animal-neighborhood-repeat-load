"""Pure sparse physical-bp shifts for Root's frozen R3 contract.

NO file/source/API IO, observed regressions, execution loop, RNG creation, seed
selection, mode selection, retries or scientific execution authority. NumPy and
SciPy imports are LAZY, inside numeric functions only; importing this module
loads only the standard library. Root executes NEW fixtures and reviews before
any exact release. All genomic intervals/delta sets are integer half-open.

Positive delta moves annotation A forward; callable fragments/masks/D are fixed.
Counts use the source preimage W-delta. Protection qualification uses WHOLE
original windows, not just their callable fragments. Root must use offset_plan
before sampling, and HOLD when status ZERO_ONLY_NOT_ESTIMABLE rather than run a
misleading null. The sampler itself is generic and does not reapply that gate.
"""
from dataclasses import dataclass
from numbers import Integral
import math

INT64_MAX=(1<<63)-1


class OffsetHold(ValueError):
    pass


class UndefinedStatistic(ValueError):
    """Preserve this draw as undefined/HOLD; never drop, redraw or impute."""
    pass


def need(value,message):
    if not value:raise ValueError(message)


def integer(value):
    return isinstance(value,Integral)and not isinstance(value,bool)


def length_value(L):
    need(integer(L)and 0<int(L)<=INT64_MAX,'positive int64 chromosome length')
    return int(L)


def merge_halfopen(intervals,L):
    """Canonical sparse union within[0,L); empty segments contribute no bases."""
    L=length_value(L);items=[]
    for pair in intervals:
        need(len(pair)==2 and all(integer(v)for v in pair),'integer half-open endpoints')
        a,b=map(int,pair);need(0<=a<=b<=L,'interval outside chromosome or reversed')
        if a<b:items.append((a,b))
    result=[]
    for a,b in sorted(items):
        if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(result[-1][1],b))
        else:result.append((a,b))
    return result


def _merge_clipped(ranges,domain):
    lo,hi=domain;items=sorted((max(a,lo),min(b,hi))for a,b in ranges if max(a,lo)<min(b,hi));result=[]
    for a,b in items:
        if result and a<=result[-1][1]:result[-1]=(result[-1][0],max(result[-1][1],b))
        else:result.append((a,b))
    return result


def delta_domain(windows,L,circular=True):
    """Half-open integer domain, inclusive linear endpoints encoded as hi+1.

    circular:[0,L); linear:[max(Wend)-L,min(Wstart)+1). Windows remain fixed.
    """
    L=length_value(L);need(type(circular)is bool,'explicit circular/linear mode')
    W=merge_halfopen(windows,L);need(W,'at least one nonempty original window')
    return (0,L)if circular else (max(b for a,b in W)-L,min(a for a,b in W)+1)


def forbidden_delta_ranges(windows,blocked,L,circular=True):
    """Exact integer offsets whose source preimage overlaps any blocked base.

    W[a,b), B[u,v) forbids[a-v+1,b-u), NOT[a-v,b-u]. Thus merely touching
    boundaries is allowed. Circular ranges are projected moduloL and split,
    while linear ranges are intersected with the bounded no-wrap domain.
    """
    L=length_value(L);W=merge_halfopen(windows,L);B=merge_halfopen(blocked,L)
    domain=delta_domain(W,L,circular);ranges=[]
    for a,b in W:
        for u,v in B:
            low,stop=a-v+1,b-u
            if circular:
                width=stop-low
                if width>=L:return [(0,L)]
                start=low%L;end=start+width
                if end<=L:ranges.append((start,end))
                else:ranges.extend(((start,L),(0,end-L)))
            else:ranges.append((low,stop))
    return _merge_clipped(ranges,domain)


def _complement(forbidden,domain):
    low,stop=domain;cursor=low;allowed=[]
    for a,b in forbidden:
        if cursor<a:allowed.append((cursor,a))
        cursor=max(cursor,b)
    if cursor<stop:allowed.append((cursor,stop))
    return allowed


def allowed_delta_ranges(windows,blocked,L,circular=True):
    """Exact complement; empty/zero-ineligible whole-window phase is a HOLD."""
    W=merge_halfopen(windows,L);B=merge_halfopen(blocked,L)
    domain=delta_domain(W,L,circular)
    forbidden=forbidden_delta_ranges(W,B,L,circular)
    allowed=_complement(forbidden,domain)
    if not allowed:raise OffsetHold('EMPTY_ALLOWED_DELTAS')
    if not any(a<=0<b for a,b in allowed):raise OffsetHold('ZERO_NOT_ALLOWED_IN_FROZEN_PHASE_SET')
    return allowed


def offset_plan(windows,blocked,L,circular=True):
    """Compact phase-set metadata; no offsets drawn and no data inspected.

    ZERO_ONLY_NOT_ESTIMABLE is an explicit boundary-model disposition. Root
    must not turn it into an estimable999-draw test or silently omit a species.
    """
    # Freeze iterable inputs once so generator callers retain exact semantics.
    W=merge_halfopen(windows,L);B=merge_halfopen(blocked,L)
    domain=delta_domain(W,L,circular);forbidden=forbidden_delta_ranges(W,B,L,circular)
    allowed=allowed_delta_ranges(W,B,L,circular);width=sum(b-a for a,b in allowed)
    return dict(mode='circ'if circular else'linear',domain=domain,forbidden=forbidden,allowed=allowed,
                integer_width=width,zero_eligible=True,status='ZERO_ONLY_NOT_ESTIMABLE'if width==1 else'ELIGIBLE_CONDITIONAL_PHASE_SET',
                masks_or_denominators_moved=False,permutation_group_not_assumed=True)


def sample_offset(ranges,rng):
    """One exact uniform integer ticket across widths, not uniform over ranges.

    Root supplies its frozenPCG64 child stream. Zero/repeated/small offsets are
    retained. No RNG construction, rejection sampling or automatic retries.
    """
    ranges=list(ranges)
    need(ranges,'EMPTY_ALLOWED_DELTAS');previous=None;clean=[];total=0
    for pair in ranges:
        need(len(pair)==2 and all(integer(x)for x in pair),'integer offset ranges')
        a,b=map(int,pair);need(a<b and(previous is None or a>=previous),'sorted disjoint positive-width offset ranges')
        clean.append((a,b));total+=b-a;previous=b
    need(total<=INT64_MAX,'uniform ticket fits int64 generator domain')
    ticket=rng.integers(total);need(integer(ticket)and 0<=int(ticket)<total,'RNG integer ticket outside exact width');ticket=int(ticket)
    for a,b in clean:
        if ticket<b-a:return a+ticket
        ticket-=b-a
    raise AssertionError('unreachable ticket accounting')


@dataclass(frozen=True)
class AnnotationUnion:
    starts:object
    ends:object
    cumulative:object
    L:int


@dataclass(frozen=True)
class CallableFragments:
    starts:object
    ends:object
    owners:object
    denominators:object
    target_ids:tuple
    L:int


def _readonly(array):
    array.setflags(write=False);return array


def canonical_union(intervals,L):
    """Read-only int64 interval arrays with cumulative covered bp; no denseL array."""
    import numpy as np
    L=length_value(L);merged=merge_halfopen(intervals,L)
    starts=np.asarray([a for a,b in merged],dtype=np.int64);ends=np.asarray([b for a,b in merged],dtype=np.int64)
    cumulative=np.empty(len(merged)+1,dtype=np.int64);cumulative[0]=0
    if merged:np.cumsum(ends-starts,dtype=np.int64,out=cumulative[1:])
    return AnnotationUnion(_readonly(starts),_readonly(ends),_readonly(cumulative),L)


def coverage_prefix(A,x):
    """Vectorized exact int64 F(x)=|A∩[0,x)|, preserving x's shape.

    x must be integer within0..L inclusive. Scalar returns a0-dimensional int64
    array. No chromosome-length storage and no float prefix interpolation.
    """
    import numpy as np
    need(isinstance(A,AnnotationUnion),'canonical AnnotationUnion required')
    values=np.asarray(x)
    if values.size==0:return np.zeros(values.shape,dtype=np.int64)
    need(values.dtype.kind in('i','u')and values.dtype.kind!='b','integer prefix coordinates')
    need(bool(np.all(values>=0))and bool(np.all(values<=A.L)),'prefix endpoints within0..L')
    points=values.astype(np.int64,copy=False);flat=points.reshape(-1)
    if len(A.starts)==0:return np.zeros(points.shape,dtype=np.int64)
    idx=np.searchsorted(A.ends,flat,side='right')
    result=A.cumulative[idx].copy();active=idx<len(A.starts)
    if np.any(active):
        selected=idx[active];partial=np.maximum(flat[active]-A.starts[selected],0)
        result[active]+=partial
    return result.reshape(points.shape)


def prepare_callable_fragments(fragments_by_target,L,target_ids=None,expected_D=None):
    """One frozen sparse callable-fragment preparation; preserve ALL target IDs.

    Root passes original1896 order explicitly in production. Empty targets stay
    present withD/count0; no per-draw qualification/mask changes. `expected_D`
    optionally checks the independently accepted frozen denominator vector.
    """
    import numpy as np
    L=length_value(L);ids=tuple(fragments_by_target)if target_ids is None else tuple(target_ids)
    need(ids and len(ids)==len(set(ids))and set(ids)==set(fragments_by_target),'exact unique complete target order')
    need(all(isinstance(g,str)and g for g in ids),'old-group target identities')
    starts=[];ends=[];owners=[];D=[]
    for i,g in enumerate(ids):
        merged=merge_halfopen(fragments_by_target[g],L);D.append(sum(b-a for a,b in merged))
        for a,b in merged:starts.append(a);ends.append(b);owners.append(i)
    if expected_D is not None:
        expected=list(expected_D);need(len(expected)==len(ids)and all(integer(v)and int(v)==d for v,d in zip(expected,D)),'fixed callable denominator mismatch')
    arrays=[np.asarray(v,dtype=np.int64)for v in(starts,ends,owners,D)]
    return CallableFragments(*[_readonly(v)for v in arrays],ids,L)


def query_counts(A,callablefragments,delta,L,circular=True):
    """One fixed-fragment sparse query; positive annotation shift via preimage.

    Circular source intervals split at0/L exactly, including a full-L fragment.
    Linear caller must have admitted whole-window delta; this additionally
    checks every callable fragment preimage in bounds. No masks/D are shifted.
    """
    import numpy as np
    L=length_value(L);need(isinstance(A,AnnotationUnion)and A.L==L,'same canonical annotation chromosome')
    fr=callablefragments
    if not isinstance(fr,CallableFragments):fr=prepare_callable_fragments(fr,L)
    need(fr.L==L and integer(delta)and type(circular)is bool,'same frozen fragment chromosome/integer delta/mode')
    delta=int(delta)
    if circular:need(0<=delta<L,'circular delta in0..L-1')
    elif len(fr.starts):need(int(fr.starts.min())-delta>=0 and int(fr.ends.max())-delta<=L,'linear preimages outside chromosome')
    result=np.zeros(len(fr.target_ids),dtype=np.int64)
    if not len(fr.starts):return result
    if circular:
        low=(fr.starts-delta)%L;width=fr.ends-fr.starts;to_end=L-low
        first_end=low+np.minimum(width,to_end)
        pieces=coverage_prefix(A,first_end)-coverage_prefix(A,low)
        second_end=np.maximum(width-to_end,0)
        pieces+=coverage_prefix(A,second_end)
    else:
        low=fr.starts-delta;end=fr.ends-delta
        pieces=coverage_prefix(A,end)-coverage_prefix(A,low)
    np.add.at(result,fr.owners,pieces)
    need(bool(np.all(result>=0))and bool(np.all(result<=fr.denominators)),'exact count/unchanged denominator closure')
    return result


def _vector(values,name):
    import numpy as np
    result=np.asarray(values,dtype=np.float64)
    need(result.ndim==1 and len(result)>0 and bool(np.all(np.isfinite(result))),'finite one-dimensional '+name)
    return result


def average_ranks(values):
    """Exact ties receive average ranks; no binning, jitter or random tie breaks."""
    from scipy.stats import rankdata
    return rankdata(_vector(values,'response'),method='average')


def fixed_qr(X):
    """One reduced QR of frozen originalX, not any observed-response refit.

    Rank-deficient designs HOLD rather than choose/drop columns or substitute a
    pseudoinverse. Numerical rank uses NumPy's ordinary matrix_rank admission;
    Root must verify this equals its frozen originalX column contract once.
    The returnedQ is read-only and reused unchanged across all null responses.
    """
    import numpy as np
    matrix=np.asarray(X,dtype=np.float64)
    need(matrix.ndim==2 and matrix.shape[0]>=matrix.shape[1]>0 and bool(np.all(np.isfinite(matrix))),'finite tall frozen originalX')
    need(np.linalg.matrix_rank(matrix)==matrix.shape[1],'rank-deficient frozen originalX; no automatic redesign')
    Q,_=np.linalg.qr(matrix,mode='reduced')
    return _readonly(Q)


def residualize(y,Q):
    """Only y-Q@(Q.T@y); never recomputeX, QR, coefficients or observation y."""
    import numpy as np
    values=_vector(y,'new null response');basis=np.asarray(Q,dtype=np.float64)
    need(basis.ndim==2 and basis.shape[0]==len(values)and 0<basis.shape[1]<=len(values)and bool(np.all(np.isfinite(basis))),'frozen original QR basis shape')
    return values-basis@(basis.T@values)


def response_residuals(counts,D,Q):
    """New null A/D -> tie-average ranks -> original fixed QR projection only.

    Root assembles all chromosomes in the fixed original1896 order BEFORE this
    global ranking; do not rank separately within chromosomes.
    Constant ranked response is explicitly undefined, not turned into tiny
    floating QR residuals or silently dropped. Root retains undefined draw.
    """
    import numpy as np
    c=np.asarray(counts);d=np.asarray(D)
    need(c.ndim==d.ndim==1 and c.shape==d.shape and c.dtype.kind in('i','u')and d.dtype.kind in('i','u'),'integer count and frozen denominator vectors')
    need(bool(np.all(d>0))and bool(np.all(c>=0))and bool(np.all(c<=d)),'positive unchangedD/count bounds; no draw-specific deletion')
    ranks=average_ranks(c.astype(np.float64)/d.astype(np.float64))
    if np.all(ranks==ranks[0]):raise UndefinedStatistic('CONSTANT_NEW_RANKED_RESPONSE')
    return residualize(ranks,Q)


def pearson(x,y):
    """Residual Pearson, NO re-ranking. Exact zero variance raises, never omit."""
    import numpy as np
    a=_vector(x,'residual x');b=_vector(y,'residual y')
    need(a.shape==b.shape and len(a)>=2,'paired residual vector length')
    a=a-a.mean();b=b-b.mean();na=np.linalg.norm(a);nb=np.linalg.norm(b)
    if na==0 or nb==0:raise UndefinedStatistic('ZERO_RESIDUAL_VARIANCE')
    value=float(np.dot(a/na,b/nb))
    need(math.isfinite(value),'nonfinite residual Pearson')
    # Roundoff may produce1+epsilon; no scientific threshold or tail filtering.
    return max(-1.,min(1.,value))


def mean_three_pair_pearson(G,P,S):
    """Frozen three-pair residual mean; P must be the unchanged saved vector."""
    return (pearson(G,P)+pearson(G,S)+pearson(P,S))/3.
