"""New pure fixed-distance geometry/alias helpers; no source or outcome IO."""
from numbers import Integral
def integer(x):
 if not isinstance(x,Integral)or isinstance(x,bool):raise ValueError('integer coordinate/offset required')
 return int(x)
def allowed_point(ranges,delta,L,circular):
 delta=integer(delta);L=integer(L);assert L>0 and type(circular)is bool
 v=delta%L if circular else delta
 return any(integer(a)<=v<integer(b)for a,b in ranges)
def inverse_pieces(a,b,delta,L,circular):
 a,b,delta,L=map(integer,(a,b,delta,L));assert 0<=a<b<=L
 if not circular:
  x,y=a-delta,b-delta
  return [(x,y)]if 0<=x<y<=L else None
 x=(a-delta)%L;y=x+(b-a)
 return [(x,y)]if y<=L else[(x,L),(0,y-L)]
def admit_direct(windows,blocked,delta,L,circular):
 for a,b in windows:
  pieces=inverse_pieces(a,b,delta,L,circular)
  if pieces is None:return False,'LINEAR_BOUNDARY'
  for x,y in pieces:
   if any(max(x,u)<min(y,v)for u,v in blocked):return False,'PROTECTED_OVERLAP'
 return True,None
def overlap_self(windows,L,delta,circular):
 total=0
 for a,b in windows:
  pieces=inverse_pieces(a,b,delta,L,circular)
  if pieces is None:raise ValueError('inadmissible linear footprint')
  total+=sum(max(0,min(b,y)-max(a,x))for x,y in pieces)
 return total
def canonical_mode(mode,all_no_wrap):
 assert mode in ['circ','linear']and type(all_no_wrap)is bool
 if mode=='linear':assert all_no_wrap
 return 'linear'if all_no_wrap else'circ'
