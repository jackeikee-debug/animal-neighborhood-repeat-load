"""Pure copy of already accepted compact Genome mask intervals, NOT Genome QA.
Caller binds original receiver/raw SHA, positive original manifest acceptance,
preoutcome panel and immutable member cache. No DNA/annotation/count rebuild.
"""
import hashlib,re
KINDS=('ACCEPTED_COMPACT_GFF3_GENOME_QA_JSON','LEGACY32_ACCEPTED_COMPACT_RECEIVER_AND_NAMESPACE_JSON')

def need(x,m):
 if not x:raise ValueError('HOLD:'+m)
def fields(line):
 parts=line.split('\t');out={}
 for part in parts[1:]:
  k,s,v=part.partition('=');need(s and k and k not in out,'literal compact fields');out[k]=v
 return parts[0],out

def copy_compact_masks(raw_lines,qa,panel,member_rows,provenance,*,target_key,report_aliases,check):
 """Consume accepted local GFF3 derived-mask artifact once, through full EOF.
 QA is already accepted. Only source/target identity, literal mask/D/clip closure
 are checked here. Original GFF/DNA counts, GC, annotations and eligibility are
 never recomputed. Missing literal narrow numeric rows stay typed NA.
 """
 check();need(provenance['kind'] in KINDS and provenance['Root_accepted_original_manifest_binding_reused'] is True,'distinct positive compact source kind and originalmanifest binding')
 for k in ('raw_ref','receiver_acceptance_ref','original_manifest_ref','original_manifest_binding_ref','Root_raw_lineage_ref','preoutcome_gate_ref'):
  p=provenance[k];need(type(p) is dict and type(p['path']) is str and p['path'].startswith('/') and re.fullmatch('[0-9a-f]{64}',p['sha256']) is not None,'explicit immutable compact source control '+k)
 need(qa['schema']=='GENOME_WINDOW_RECEIVER_V1' and qa['status']=='SOURCE_GEOMETRY_AND_COUNT_ACCOUNTING_ACCEPTED' and qa['repeat_outcomes_computed'] is False and qa['weights_or_model_computed'] is False,'original accepted receiver schema, not fabricated modernDB')
 need(qa['output_bytes']==provenance['raw_ref']['bytes'] and qa['output_sha256']==provenance['raw_ref']['sha256'],'accepted raw mask output independent receiver byte/SHA binding')
 ids=list(panel['panel']['ids'])+list(panel['panel']['local_descriptive_ids']);need(len(ids)==len(set(ids))<=2000,'frozen primary/local unique ids');rows={r['id']:r for r in panel['literal_rows50']};need(set(ids)<=set(rows),'literal preoutcome row identities')
 proteins={}
 for r in member_rows:
  check()
  if r['g'] in ids:need(r['g'] not in proteins,'unique member group');proteins[r['g']]=r['protein']
 need(set(proteins)==set(ids),'all frozen members from accepted cache');target_groups={target_key(panel['assembly'],proteins[g]):g for g in ids};need(len(target_groups)==len(ids),'unique frozen source targets')
 audits=qa['window_audits'];need(len(audits)<=4000 and len(audits)==qa['terminal']['windows'],'alreadyaccepted source window count')
 selected={};qa_bywid={}
 for wid,a in enumerate(audits,1):
  check();need(a['window_index']==wid and a['radius'] in (50000,25000),'original exact window index/radius');qa_bywid[wid]=a
  if a['target'] in target_groups:
   g=target_groups[a['target']];key=(g,a['radius']);need(key not in selected,'one frozen target/radius');selected[key]=dict(a,masks=[])
 need(len(selected)==2*len(ids),'all exact frozen requested pairs in accepted receiver')
 source_n=source_lines=mask_rows=window_audits_seen=0;terminal=False;digest=hashlib.sha256();ordinals={};lasts={}
 for raw in raw_lines:
  check();need(type(raw) is bytes and raw.endswith(b'\n') and len(raw)<=2**20 and b'\r' not in raw,'literal LF compact row interface');source_n+=len(raw);source_lines+=1;need(source_n<=qa['output_bytes'] and not terminal,'exact accepted raw output bounds/oneEOF');digest.update(raw);line=raw[:-1].decode('ascii')
  if source_lines==1:need(line=='##gff-version 3','originalcompactheader');continue
  if source_lines==2:
   tag,f=fields(line);need(tag=='#GENOME_WINDOW_SCHEMA' and f==dict(version='1',semantics='ALL_original_exons_OR_nonACGT',manifest_sha256=qa['manifest_sha256']),'original receiver-bound schema/manifest semantic identity');continue
  if not line.startswith('#'):
   f=line.split('\t');need(len(f)==9 and f[1:3]==['derived_exon_or_nonACGT','region'] and f[5:8]==['.','.','.'],'literal accepted mask region kind')
   match=re.fullmatch(r'ID=genome_mask_w([1-9][0-9]*)_([1-9][0-9]*);codex_window=([1-9][0-9]*);codex_mask_basis=ALL_original_exons_OR_nonACGT',f[8]);need(match is not None and match[1]==match[3],'literal region window and maskbasis')
   wid=int(match[1]);need(wid in qa_bywid,'mask window index from originalQA');a=qa_bywid[wid];left=int(f[3])-1;right=int(f[4]);need(int(match[2])==ordinals.get(wid,0)+1 and f[0]==a['canonical_seqid'] and 0<=left<right and a['requested_start0']<=left<right<=a['requested_end0'] and (wid not in lasts or left>lasts[wid]),'literal mask interval/order/window')
   ordinals[wid]=ordinals.get(wid,0)+1;lasts[wid]=right;mask_rows+=1
   if a['target'] in target_groups:selected[target_groups[a['target']],a['radius']]['masks'].append([left,right])
  elif line.startswith('#GENOME_WINDOW_AUDIT\t'):
   _,f=fields(line);window_audits_seen+=1;need(window_audits_seen in qa_bywid,'literal source window audit count');a=qa_bywid[window_audits_seen]
   need(all(f[k]==str(a[k]) for k in ('window_index','target','radius','canonical_seqid','requested_start0','requested_end0')),'alreadyaccepted raw/QA target-window identities')
  elif line.startswith('#GENOME_EOF_OK\t'):
   _,f=fields(line);need({k:int(v) for k,v in f.items()}==qa['terminal'],'literal accepted terminal unchanged');terminal=True
  else:need(line.startswith('#TARGET_NA\t') or line.startswith('#FASTA_SEQUENCE_AUDIT\t'),'known original receiver record family')
 need(terminal and source_n==qa['output_bytes'] and source_lines==qa['output_lines'] and digest.hexdigest()==qa['output_sha256'] and window_audits_seen==len(audits) and mask_rows==qa['terminal']['mask_segments'],'complete accepted rawmask EOF/bytes/hash/count identity')
 sequences={};canonical={}
 for z in qa['sequences']:
  check();token=z['token'];canon=z['canonical_seqid'] or None;length=z['length'];need(token not in sequences and token.encode('ascii') in report_aliases and report_aliases[token.encode('ascii')][1]==length,'accepted actual Genome token/length owner')
  if canon is not None:need(canon not in canonical and canon.encode('ascii') in report_aliases and report_aliases[canon.encode('ascii')]==report_aliases[token.encode('ascii')],'accepted actual canonical same reportowner');canonical[canon]=token
  sequences[token]=dict(canonical=canon,owner=report_aliases[token.encode('ascii')][0],length=length)
 windows=[];total_selected_segments=0
 for g in ids:
  for r in (50000,25000):
   check();a=selected[g,r];big=rows[g];nr=panel['literal_narrow_records'][g] if r==25000 else None;expected=big if r==50000 else nr['row'];lo,hi=(big['start0'],big['end0']) if r==50000 else (big['start0']+25000,big['end0']-25000);need((a['requested_start0'],a['requested_end0'])==(lo,hi),'source mask aligned to exact frozen A50/A25')
   if expected is None:need(nr['status']=='TYPED_NA' and nr['reasons'],'literal narrow technicalNA, no fake numericrow')
   else:need(expected['nonexonic_acgt']==a['nonexonic_acgt_bp'] and expected['complete']==bool(a['geometric_complete']),'frozen D/completeness unchanged')
   token=a['actual_token'];canon=a['canonical_seqid'];need(token in sequences and canonical.get(canon)==token and a['sequence_length']==sequences[token]['length'],'positive original actual source sequence binding')
   pairs=a['masks'];need(len(pairs)==a['combined_mask_segments'] and sum(z-l for l,z in pairs)==a['combined_mask_bp']==a['observed_bp']-a['nonexonic_acgt_bp'],'literal selected mask/D/count closure');total_selected_segments+=len(pairs)
   windows.append(dict(id=g,target_key=a['target'],wid=a['window_index'],radius=r,seqid=canon,actual_token=token,owner=sequences[token]['owner'],sequence_length=sequences[token]['length'],start0=lo,end0=hi,D=a['nonexonic_acgt_bp'],complete=bool(a['geometric_complete']),masks=pairs,mask_basis='ALL_original_exons_OR_nonACGT',frozen_narrow_record=nr,frozen_main_J=g in panel['panel']['ids'],frozen_local_descriptive=g in panel['panel']['local_descriptive_ids']))
  w,n=windows[-2:];need(n['masks']==[[max(a,n['start0']),min(z,n['end0'])] for a,z in w['masks'] if max(a,n['start0'])<min(z,n['end0'])],'literal accepted A25 mask clip equality')
 return dict(schema='ROOT_FROZEN_J_AND_LOCAL_ACCEPTED_MASK_INPUT_V1',native_species_name=panel['native_species_name'],assembly=panel['assembly'],ordered_primary_ids=list(panel['panel']['ids']),ordered_local_descriptive_ids=list(panel['panel']['local_descriptive_ids']),windows=windows,sequence_owners=[dict(token=t,canonical=z['canonical'],owner=z['owner'],length=z['length']) for t,z in sequences.items()],mask_segments=total_selected_segments,source_fact_records_observed=len(audits),source_kind=provenance['kind'],original_compact_mask_source_ref=provenance['raw_ref'],original_manifest_binding_reused_ref=provenance['original_manifest_binding_ref'],original_full_source_mask_rows_copied_or_skipped=mask_rows,eligibility_weights_DNA_QC_recomputed=False,new_RM_or_models=False)
