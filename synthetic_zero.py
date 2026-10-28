"""New synthetic arrays only; immutable accepted numerical functions, no real inputs."""
from pathlib import Path
import sys,json,hashlib,time,os
import numpy as np
V=Path(os.environ['PFLAVA_SYNTHETIC_ROOT']).resolve()
sys.path.insert(0,str(Path(__file__).resolve().parent))
from physical_shift_root_numeric_gates import fixed_qr,response_residuals,pearson
N=128
def counts(z):return np.rint(32/(1+np.exp(-np.asarray(z)))).astype(np.int64)
def rep_ids(phase):
 if phase=='pilot':return list(range(2))
 if phase=='remainder':return list(range(2,100))
 raise ValueError(phase)
def allowed(scenario):return np.arange(16,112,dtype=np.int64)if scenario=='S3'else np.arange(128,dtype=np.int64)
def design():
 t=2*np.pi*np.arange(N)/N
 return fixed_qr(np.column_stack([np.ones(N),np.sin(t),np.cos(t),np.sin(2*t),np.cos(2*t)]))
def residual(c,Q):return response_residuals(np.asarray(c,dtype=np.int64),np.full(N,32,dtype=np.int64),Q)
def statistic(G,P,S):return float((pearson(G,P)+pearson(G,S)+pearson(P,S))/3)
def upper_tail(obs,inner):
 x=np.asarray(inner);assert x.shape==(99,)and np.isfinite(x).all()and np.isfinite(obs)
 return (1+int(np.count_nonzero(x>=obs-1e-12)))/100
def generate(scenario,rep):
 assert scenario in ['S1','S2','S3']
 rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([20261028,int(scenario[1]),rep])))
 eps=rng.standard_normal((3,N));z=sum(w*np.roll(eps,j-2,axis=1)for j,w in enumerate([1,2,3,2,1]))/np.sqrt(19)
 if scenario=='S2':
  sigma=np.full(N,.2);sigma[48:80]=3.;z=z*sigma
 c=counts(z)
 initial=[]
 for k in range(3):
  phase=int(rng.integers(0,N))if scenario!='S2'else 0;c[k]=np.roll(c[k],phase);initial.append(phase)
 return c,rng,initial
def write_json(path,obj):
 with Path(path).open('x')as f:json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def run(phase):
 Q=design();claims=V/'replicate-claims';assert claims.is_dir();out=Path.cwd()/'replicates';out.mkdir();start=time.monotonic();done=[]
 for scenario in ['S1','S2','S3']:
  for rep in rep_ids(phase):
   label=f'{scenario}-r{rep:03d}';claim=claims/(label+'.CLAIMED.json')
   write_json(claim,dict(scenario=scenario,replicate=rep,seed=[20261028,int(scenario[1]),rep],phase=phase,automatic_retry=False))
   c,rng,initial=generate(scenario,rep);a=allowed(scenario);shifts=rng.choice(a,size=(99,2),replace=True);P=residual(c[1],Q);G=residual(c[0],Q);S=residual(c[2],Q);obs=statistic(G,P,S);values=[]
   for dg,ds in shifts:values.append(statistic(residual(np.roll(c[0],int(dg)),Q),P,residual(np.roll(c[2],int(ds)),Q)))
   p=upper_tail(obs,values);record=dict(label=label,scenario=scenario,replicate=rep,seed=[20261028,int(scenario[1]),rep],initial_phases=initial,counts=c.tolist(),shifts=shifts.tolist(),observed_T=obs,inner_T=values,tail=p,reject05=p<=.05,reject10=p<=.10,rank_QR_projections=201,real_models_refitted=0)
   dest=out/(label+'.json');write_json(dest,record);h=hashlib.sha256(dest.read_bytes()).hexdigest();write_json(claims/(label+'.RETURNED.json'),dict(path=str(dest),sha256=h));done.append(label)
 write_json(Path.cwd()/'return.json',dict(phase=phase,replicates=done,n=len(done),rank_QR_projections=201*len(done),inner_transforms=99*len(done),seconds=time.monotonic()-start,real_source_reads=0,old_result_recomputation=0,undefined_redraws=0))
if __name__=='__main__':run(sys.argv[1])
