"""Complete saved-return checks; arithmetic checks saved deletion summary only, no fits."""
import math
def validate_backend_return(inp,result):
 rows=inp['effective']['rows'];ids=[r['id'] for r in rows];w=inp['effective']['weights'];assert ids==inp['active_ids'] and result['excluded']==inp['excluded']
 full=result['full'];assert full['status'] in ('OK','NA')
 if 'n' in full:assert full['n']==len(ids)
 if full['status']=='OK':assert full['ids']==ids and full['weights']==w and type(full['beta']) in (int,float) and math.isfinite(full['beta'])
 else:assert full.get('beta') is None and full.get('reason')
 info=result['information'];assert info['status'] in ('OK','NA') and type(info['spatial_eligible']) is bool
 if 'n' in info:assert info['n']==len(ids)
 deleted=[]
 for d in result['deletions']:
  loc=d['deleted_chromosome'];assert loc not in deleted;deleted.append(loc);keep=[i for i,r in enumerate(rows) if r['location']!=loc]
  if 'n' in d:assert d['n']==len(keep)
  if d['status']=='OK':assert d['ids']==[ids[i] for i in keep] and d['weights']==[w[i] for i in keep] and math.isfinite(d['beta'])
  else:assert d['status']=='NA' and d.get('beta') is None and d.get('reason')
 if deleted:assert deleted==sorted({r['location'] for r in rows})
 if result['status']=='OK':
  assert full['status']=='OK' and len(deleted)>=2 and all(d['status']=='OK' for d in result['deletions']);b=[d['beta'] for d in result['deletions']];L=len(b);mean=sum(b)/L;v=(L-1)/L*sum((x-mean)**2 for x in b)
  assert math.isclose(v,result['vJ'],rel_tol=1e-10,abs_tol=1e-12) and result['beta']==full['beta'] and result['deletion_min']==min(b) and result['deletion_max']==max(b)
  sign=lambda x:1 if x>0 else -1 if x<0 else 0
  assert result['sign_changes']==sum(sign(x)!=sign(full['beta']) for x in b)
 else:assert result['status']=='NA' and result.get('reason') and result.get('vJ') is None
