"""Pure new sensitivity input adapter. No IO, source selection, fits or primary reweighting."""
from copy import deepcopy
import hashlib,json,math
BACKEND_SHA='36a901dce504df25c8e2d8cc0f9df2810b41ceb004d4d64b2d53e149699f0c26'
READY_A25='PREPARED_FROM_EXACT_ACCEPTED_HELPERS'
def need(x,s):
 if not x:raise ValueError('HOLD: '+s)
def select_panel(base,ids,weights):
 original=base['panel'];need(len(ids)==len(set(ids))==len(weights),'unique selected IDs and same length weights')
 positions={g:i for i,g in enumerate(original['ids'])};need(set(ids)<=set(positions),'selected IDs in accepted J')
 need([positions[g] for g in ids]==sorted(positions[g] for g in ids),'original per-species J order')
 need(all(type(w) in (int,float) and math.isfinite(w) and w>0 for w in weights),'finite positive accepted weights')
 need([r['id'] for r in original['rows']]==original['ids'],'accepted full row order')
 byid={r['id']:r for r in original['rows']}
 # These are exactly the fields consumed by original _variant/spatial_jackknife.
 # Omitted component labels are not needed numerically and are not fabricated.
 return {'ids':list(ids),'rows':[deepcopy(byid[g]) for g in ids],'weights':list(weights)}
def active_NA(observations,ids):
 need(set(ids)<=set(observations),'all active observations supplied including NA')
 bad=[]
 for g in ids:
  o=observations[g]
  if o.get('status')=='TYPED_NA':need(bool(o.get('reason')),'typed NA reason');bad.append(g)
  else:need(o.get('status') in (None,'MEASURED') and type(o['A']) is int and type(o['U']) is int and 0<=o['U']<=o['A'],'measured integer A/U')
 return bad
def input_fingerprint(effective,json_safe):
 identity={'schema':'EXACT_SPECIES_VARIANT_NUMERIC_ARRAYS_V1','backend_sha256':BACKEND_SHA,'runtime':'Py3.12.14/NumPy2.3.5/SciPy1.16.3/one_BLAS_thread','effective':effective}
 return hashlib.sha256(json.dumps(json_safe(identity),sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode()).hexdigest()
def prepare(base,ids,weights,variant,observations,backend,json_safe,*,accepted_baseline=False):
 need(variant in ('A50','K50','A25'),'finite original variant')
 panel=select_panel(base,ids,weights);need(set(observations)==set(ids),'exact main J observations, including inactive A25 slots')
 result={'status':'TYPED_NA_NO_NUMERIC_CALL','native_species_name':base['native_species_name'],'assembly':base['assembly'],'variant':variant,'main_panel':panel,'observations':deepcopy(observations),'narrow_rows':None,'excluded':[],'numeric_fingerprint':None,'active_NA_ids':[],'baseline_input_reconstruction':accepted_baseline,'primary_reweighting_calls':0}
 if not ids:result['reason']='EMPTY_MAIN_J';return result
 if variant=='A25':
  a=base['A25']
  if a['status']!=READY_A25:result.update(reason='FROZEN_A25_TECHNICAL_METADATA_UNAVAILABLE',A25_status=a['status']);return result
  need(set(a['narrow_rows'])==set(base['panel']['ids']),'accepted A25 narrow mainJ closure')
  narrow={g:deepcopy(a['narrow_rows'][g]) for g in ids};result['narrow_rows']=narrow
  if accepted_baseline:
   need(ids==base['panel']['ids'] and weights==base['panel']['weights'],'baseline reconstruction is exact accepted J and weights')
   rows=deepcopy(a['rows']);w=list(a['weights']);x=list(a['xA']);excluded=deepcopy(a['excluded'])
   active=[r['id'] for r in rows];need(active==[g for g in ids if g in set(active)],'accepted A25 active order')
   index={g:i for i,g in enumerate(ids)};need(w==[weights[index[g]] for g in active],'accepted narrow corresponding weights')
  else:rows,w,x,excluded=backend._variant(panel,variant,narrow)
 elif accepted_baseline:
  need(ids==base['panel']['ids'] and weights==base['panel']['weights'],'baseline accepted inputs')
  rows=deepcopy(panel['rows']);w=list(weights);x=[r['xK' if variant=='K50' else 'xA'] for r in rows];excluded=[]
 else:rows,w,x,excluded=backend._variant(panel,variant)
 active=[r['id'] for r in rows];result['excluded']=excluded;result['active_ids']=active
 bad=active_NA(observations,active);result['active_NA_ids']=bad
 if bad:result['reason']='FROZEN_OBSERVATIONS_HAVE_TYPED_NA';return result
 if not active:result['reason']='NO_ACTIVE_VARIANT_ROWS';return result
 for r in rows:need(0<=observations[r['id']]['U']<=observations[r['id']]['A']<=r['nonexonic_acgt'],'outcome bounds in original current denominator')
 effective={'species':base['native_species_name'],'assembly':base['assembly'],'variant':variant,'rows':rows,'weights':list(w),'x':list(x),'observations':{g:deepcopy(observations[g]) for g in active}}
 result.update(status='PREPARED_NUMERIC_INPUT',effective=effective,numeric_fingerprint=input_fingerprint(effective,json_safe));return result
