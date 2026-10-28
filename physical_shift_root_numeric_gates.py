"""Narrow Root numerical admission around unchanged reviewed interval code."""
import numpy as np
import physical_shift as base
from physical_shift import AnnotationUnion,prepare_callable_fragments,query_counts,pearson,average_ranks,UndefinedStatistic
TOL=1e-10
_FACTORY_BASES={}
def fixed_qr(X):
    Q=base.fixed_qr(X)
    assert np.max(np.abs(Q.T@Q-np.eye(Q.shape[1])))<=TOL,'fixed_Q_orthogonality'
    _FACTORY_BASES[id(Q)]=Q
    return Q
def residualize(y,Q):
    if _FACTORY_BASES.get(id(Q))is not Q:raise ValueError('only_accepted_readonly_fixed_Q_factory')
    y=np.asarray(y,dtype=np.float64);e=base.residualize(y,Q)
    response_norm=float(np.linalg.norm(y));centered_norm=float(np.linalg.norm(y-y.mean()))
    if centered_norm==0:raise UndefinedStatistic('CONSTANT_NEW_RANKED_RESPONSE')
    normal=float(np.linalg.norm(Q.T@e))/max(response_norm,np.finfo(float).tiny)
    assert normal<=TOL,'new_response_projection_normal_equation'
    if float(np.linalg.norm(e))/centered_norm<=TOL:raise UndefinedStatistic('NEW_RESPONSE_NO_NUMERIC_RESIDUAL_INFORMATION')
    return e
def response_residuals(counts,D,Q):
    c=np.asarray(counts);d=np.asarray(D)
    assert c.ndim==d.ndim==1 and c.shape==d.shape and c.dtype.kind in('i','u')and d.dtype.kind in('i','u')
    assert np.all(d>0)and np.all(c>=0)and np.all(c<=d)
    return residualize(average_ranks(c.astype(float)/d.astype(float)),Q)
