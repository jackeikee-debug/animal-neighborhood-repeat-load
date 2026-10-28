"""R15 first finite clade-common input projection; saved cache only, no models."""
from pathlib import Path
import collections,gc,gzip,hashlib,importlib.util,json,sys,time

def load(n,p):
 s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);sys.modules[n]=m;s.loader.exec_module(m);return m

def promote(known,verified,proof):
 out={k:dict(v) for k,v in known.items()};n=0
 for k,v in out.items():
  if v['kind']=='NEW_PHASE2_INPUT_PENDING_NUMERIC_RETURN':
   sm=verified[v['numeric_fingerprint']];assert sm['consumer_signature']['sha256']==k and sm['full_status']=='OK'
   v.update(kind='ACCEPTED_P2_CROSS_SCALE_RETURN',source_result_ref=sm['source_result_ref'],summary=sm,complete_body_acceptance_ref=proof);n+=1
  else:assert v['kind'] in ('ACCEPTED_PC_RETURN','ACCEPTED_B0_RETURN')
 assert n==114 and len(verified)==114
 return out

def build(cfg,out,check,events):
 start=time.monotonic();source=load('r15source',cfg['codes']['sources']);src=source.Sources(cfg,check,events);safe=load('r15safe',cfg['codes']['json_safe']);restore=load('r15restore',cfg['codes']['restore_bytes']);fpmod=load('r15fp',cfg['codes']['effect_inputs']);cross=load('cross_scale_inputs',cfg['cross_adapter']);adapter=load('r15adapter',cfg['adapter']);out=Path(out);out.mkdir();(out/'inputs').mkdir();(out/'manifests').mkdir()
 def read(n):
  parts=[];h=hashlib.sha256();nb=0
  for line in src.lines(n):
   raw=line.encode();nb+=len(raw);assert nb<=32*2**20;h.update(raw);parts.append(line)
  s=src.sources[n]
  if 'uncompressed_bytes' in s:assert nb==s['uncompressed_bytes'] and h.hexdigest()==s['uncompressed_sha256']
  return restore.restore_exact_bytes(json.loads(''.join(parts)))
 def packed(p,obj):
  raw=json.dumps(safe.json_safe(obj),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()+b'\n'
  with p.open('xb') as f:
   with gzip.GzipFile(fileobj=f,mode='wb',mtime=0) as z:z.write(raw)
  b=p.read_bytes();return dict(path=str(p),bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),gzip=True,uncompressed_bytes=len(raw),uncompressed_sha256=hashlib.sha256(raw).hexdigest())
 cert=read('cache_acceptance');index=read('cache_index');plans=read('plans');qa=read('cross_QA');known=read('consumer_cache');verified=read('verified114');assert cert['accepted_profiles']==38 and cert['all_owned_exactSTOP'];assert qa['complete_new_pairs']==114 and qa['all_owned_exactSTOP'] and qa['new_numeric_claim_UNKNOWN']==0
 known=promote(known,verified,cfg['cross_QA_ref']);source.write_json(out/'promoted-B0-PC501-P2cross114-consumer-cache.json',known)
 accepted={r['native_species_id']:r for r in cert['cache_refs']};byid={r['native_species_id']:r for r in index};assert len(plans)==10 and sum(len(p['native_species_ids']) for p in plans)==48
 planby=collections.defaultdict(list)
 for plan in plans:
  assert plan['common_count']==len(plan['common_old_group_ids']) and plan['no_outcome_or_sign_selection']
  for sid in plan['native_species_ids']:planby[sid].append(plan)
 logical=[];manifests=[];counts=collections.Counter();newcount=0
 for sid in sorted(planby,key=int):
  check();ref=byid[sid]['cache'];assert all(ref[k]==accepted[sid][k]==src.sources['cache'+sid][k] for k in ('path','bytes','sha256'))
  c=read('cache'+sid);assert c['native_species_id']==sid and c['schema']==cross.CACHE_SCHEMA;inputs={}
  for plan in planby[sid]:
   p=adapter.prepare_clade_common(c,plan['common_old_group_ids'],plan['layer'],plan['rule']);purpose=plan['layer']+'__'+plan['rule']
   row=dict(native_species_id=sid,species=c['native_species_name'],assembly=c['assembly'],purpose=purpose,layer=plan['layer'],rule=plan['rule'],status=p['status'],reason=p['reason'],logical_labels=p['logical_labels'],summary=p['common_member_summary'],source_cache_ref=ref,original_sourceNA_records=p['original_sourceNA_records'],original_paired_radius_sourceNA_records=p['original_paired_radius_sourceNA_records'],selected_sourceNA_ids=p['selected_sourceNA_ids'],sourceNA_outside_predefined_common_ids=p['sourceNA_outside_predefined_common_ids'],minimum_anchor_gate=p['minimum_anchor_gate'],no_fit_reasons=p['no_fit_reasons'],consumer_signature=p['consumer_signature'])
   if p['status']!='PREPARED_NUMERIC_INPUT':row['status']='TYPED_NA_NO_NUMERIC_CALL'
   else:
    k=p['consumer_signature']['sha256'];fp=fpmod.input_fingerprint(p['effective'],safe.json_safe)
    if k in known:
     v=known[k];row.update(status='ALIAS_TO_NEW_UNIQUE_INPUT_PENDING_RETURN' if v['kind']=='NEW_CLADE_INPUT_PENDING_NUMERIC_RETURN' else 'REUSE_EXACT_CONSUMER_INPUT_RETURN',reused=v)
    else:
     q=dict(schema='FROZEN_SENSITIVITY_NUMERIC_INPUT_V1',status='PREPARED_NUMERIC_INPUT',native_species_name=c['native_species_name'],assembly=c['assembly'],variant='A50',main_panel=p['minimal_panel'],observations=p['observations'],narrow_rows=None,excluded=[],active_ids=p['minimal_panel']['ids'],active_NA_ids=[],effective=p['effective'],numeric_fingerprint=fp,source_lineage=p['source_lineage'],phase2_logical_labels=p['logical_labels'],phase2_selection_excluded=p['selection_excluded'],consumer_signature=p['consumer_signature'])
     r=packed(out/'inputs'/(fp+'.json.gz'),q);inputs[fp]=r;known[k]=dict(kind='NEW_CLADE_INPUT_PENDING_NUMERIC_RETURN',numeric_fingerprint=fp,input_ref=r,native_species_id=sid,consumer_signature=p['consumer_signature']);row.update(status='NEW_UNIQUE_NUMERIC_INPUT',numeric_fingerprint=fp,input_ref=r);newcount+=1
   counts[row['status']]+=1;logical.append(row)
  m=dict(status='NEW_FINITE_NUMERIC_INPUTS_PREPARED_NO_FITS_PENDING_ROOT',native_species_id=sid,native_species_name=c['native_species_name'],assembly=c['assembly'],input_refs=inputs,new_or_old_model_calls=0,phase2_purpose='R15 fixedclade common member A50',source_cache_ref=ref,source_receipts={'cache'+sid:src.done['cache'+sid]})
  p=out/'manifests'/('native-'+sid+'.json');source.write_json(p,m);b=p.read_bytes();manifests.append(dict(native_species_id=sid,path=str(p),bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),new_unique_inputs=len(inputs)));del c;gc.collect()
 assert len(logical)==48 and newcount<=48 and set(src.done)==set(src.sources)
 with (out/'logical-clade-common-dispositions.jsonl').open('x') as f:
  for r in logical:f.write(json.dumps(safe.json_safe(r),ensure_ascii=False,sort_keys=True,allow_nan=False)+'\n')
 source.write_json(out/'clade-inclusive-consumer-index-pending.json',known)
 source.write_json(out/'terminal.json',dict(status='FINITE48_CLADE_LOGICAL_ROWS_INPUTS_COMPLETE_PENDING_ROOT',logical_rows=48,new_unique_inputs=newcount,counts=dict(counts),manifests=manifests,new_or_old_models=0,source_receipts=src.done,seconds=time.monotonic()-start,all_original_NA_S35_OG725_preserved=True))
 return dict(numeric_call_counts={'new_input_preparation':newcount,'new_or_old_fits':0},source_receipts=src.done)
