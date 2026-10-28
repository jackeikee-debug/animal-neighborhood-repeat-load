"""Pure structure/alias binding checks added for independent-review I1."""
def need(v,why):
 if not v:raise ValueError(why)
def validate_freeze(f,expected_keys):
 rows=f['cases'];unique=f['unique_admitted_cases'];expected=set(expected_keys)
 need(len(rows)==72 and len({x['case']for x in rows})==72,'exact72 unique logicalcases')
 admitted=[x for x in rows if x['admissible']];keys=[x['canonical_key']for x in unique]
 need(len(expected_keys)==len(expected)==8 and len(keys)==len(set(keys))==8,'exact8 unique admittedkeys')
 need(set(keys)==expected and len(admitted)==16,'known16logical/8unique bound admission')
 need({x['canonical_key']for x in admitted}==expected,'logical andunique keys correspond')
 need(all(sum(x['canonical_key']==k for x in admitted)==2 for k in keys),'two modealiases perkey')
def validate_draw_rows(rows,n=999):need(len(rows)==n and[int(x['draw'])for x in rows]==list(range(1,n+1)),'ordered exact saved drawaxis')
def phase_matches(old,new,lengths):
 import numpy as np
 need(old.ndim==2 and new.ndim==lengths.ndim==1 and old.shape[1:]==new.shape==lengths.shape,'same chromosomeaxis')
 need(old.dtype==new.dtype==lengths.dtype==np.dtype('int64'),'int64 phaseaxes')
 return np.where(np.all(np.mod(old,lengths[None,:])==np.mod(new,lengths)[None,:],axis=1))[0]
def validate_pair(pair,den):
 import numpy as np
 need(pair.shape==den.shape and pair.ndim==2 and pair.shape[0]==2,'G/S byfixedmembers countshape')
 need(pair.dtype==np.dtype('uint32')and den.dtype.kind in'iu','acceptedintegercounts')
 need(bool(np.all(pair<=den)),'count bounds')
