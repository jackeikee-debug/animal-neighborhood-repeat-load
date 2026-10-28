"""Exact immutable descriptor contracts; purpose labels may differ, content may not."""
def need(x,s):
 if not x:raise ValueError('HOLD: '+s)
def check_descriptor(actual,expected):
 for k in ('path','bytes','sha256'):need(actual.get(k)==expected.get(k),'exact source descriptor '+k)
 for k in ('gzip','closed_stat'):
  if k in expected:need(actual.get(k)==expected[k],'exact source descriptor '+k)
def verify_compute_acceptance(accepted,descriptor,manifest):
 need(accepted.get('status')=='ROOT_ACCEPTED_ONE_SPECIES_NEW_FINITE_INPUTS_NO_FITS','accepted new input stage')
 for k in ('native_species_id','native_species_name','assembly'):need(accepted.get(k)==manifest[k],'accepted input identity '+k)
 check_descriptor(descriptor,accepted['prepared_manifest_ref'])
 need(accepted.get('new_model_calls')==0 and accepted.get('all_owned_exactSTOP') is True,'input stage zero models and exactSTOP')
 need(set(accepted['input_refs'])==set(manifest['input_refs']),'complete accepted numerical fingerprint set')
 for fp,ref in manifest['input_refs'].items():check_descriptor(ref,accepted['input_refs'][fp])
def baseline_seed_disposition(old,prepared):
 if old.get('new_spatial_call')!=1:return 'NO_OLD_CONSUMED_NUMERIC_CALL'
 need(old.get('call_status')=='BACKEND_RETURNED' and bool(old.get('claim_ref')),'old consumed call has accepted return and claim')
 if prepared['status']=='PREPARED_NUMERIC_INPUT':return 'SEED_ACCEPTED_NUMERIC'
 need(prepared.get('reason') in {'EMPTY_MAIN_J','NO_ACTIVE_VARIANT_ROWS'},'old consumed input reconstruction must match original active-NA rules')
 return 'KEEP_OLD_CONSUMED_EMPTY_RETURN_NO_FP'
