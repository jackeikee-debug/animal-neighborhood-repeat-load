"""Pure bounded accepted-record-view binding builder; no IO/science/driver."""
import hashlib,json,math
from copy import deepcopy
SIGNATURE_SHA='0ecca728c39dec3ea6ffdbe4ee86f2114f5163fa93af1ba6c01eacda9ecffd20'
SUPPORTED='SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY'
def need(ok,why):
 if not ok:raise ValueError(why)
def ref(r):
 need(type(r) is dict and type(r.get('path')) is str and r['path'] and type(r.get('sha256')) is str and len(r['sha256'])==64 and all(x in '0123456789abcdef' for x in r['sha256']),'explicit accepted metadata reference');return deepcopy(r)
def finite(v):
 need(type(v) in (int,float,str),'signature numeric source field');x=float(v);need(math.isfinite(x),'finite fixed signature');return x
def key(assembly,protein):return 't'+hashlib.sha256(json.dumps([assembly,protein],ensure_ascii=False,separators=(',',':')).encode()).hexdigest()
def text(v,name):
 need(type(v) is str and v,name+' exact nonempty string');return v

def build_bindings(canonical_views,evidence_records,signature_records,collision_records,context):
 """Yield draft bindings; accept NONE unless iterator terminates normally.
 Root owns trusted complete source pins/EOF/view-count and actual output receipt.
 No HOG/pair/tree recomputation; member evidence accepted flags are copied only.
 """
 c=context;need(c['schema']=='ROOT_ACCEPTED_MEMBER_SIGNATURE_BINDING_INPUT_VIEW_V1' and c['approved_by']=='current-conversation-root' and c['complete_input_source_pins_EOF_verified'] is True,'Root accepted complete input-view authority')
 need(c['full_signature_groups']==1896 and c['signature_ref']['sha256']==SIGNATURE_SHA,'immutable full1896 signature lineage')
 sr=ref(c['signature_ref']);refs=c['refs']
 for role in ('member','evidence','collision','canonical','report','NA_authority'):ref(refs[role])
 if 'genome_QA' in refs:ref(refs['genome_QA']) # no fabricated default for NA-only
 lim=c['limits'];caps={'max_targets':2000,'max_evidence_rows':100000,'max_signature_rows':2000,'max_collision_rows':100000,'max_stored_evidence_bytes':16777216,'max_matches_per_target':8}
 need(set(lim)==set(caps) and all(type(lim[k]) is int and 0<lim[k]<=v for k,v in caps.items()),'explicit bounded per-species limits')
 exp=c['expected_rows'];need(set(exp)=={'canonical','evidence','signatures','collisions'} and all(type(v) is int and v>=0 for v in exp.values()),'exact Root accepted record-view EOF counts')
 g=c['frozen_species'];assembly=text(g['assembly'],'assembly');native=text(g['native_species_name'],'native');ids=g['target_ids'];need(type(ids) is list and 0<len(ids)<=lim['max_targets'] and g['target_count']==len(ids) and all(type(x) is str and x for x in ids) and len(set(ids))==len(ids),'complete unique frozen target set')
 need(hashlib.sha256(('\n'.join(sorted(ids))+'\n').encode()).hexdigest()==g['target_set_sha256'] and exp['canonical']==len(ids),'accepted target-set membership hash/count')
 allowed=set(ids);sigs={};n=0
 for s in signature_records:
  n+=1;need(n<=lim['max_signature_rows'],'signature view cap');label=text(s['orthogroup'],'original group');need(label not in sigs,'duplicate original-group signature HOLD')
  sigs[label]={'xA':finite(s['x_A50']),'xK':finite(s['x_K50']),'sha256':SIGNATURE_SHA}
 need(n==exp['signatures'],'signature view complete EOF count')
 members={};seen_evidence=set();n=stored=0
 for e in evidence_records:
  n+=1;need(n<=lim['max_evidence_rows'],'evidence stream row cap');en=text(e['species_name'],'evidence native');protein=text(e['protein_id'],'evidence protein')
  if en!=native:continue
  need(protein in allowed,'selected-species evidence outside accepted target set HOLD');group=text(e['original_group'],'original group');ek=(en,protein,group);need(ek not in seen_evidence,'duplicate member evidence HOLD');seen_evidence.add(ek)
  need(type(e['member_supported']) is bool and type(e['post_collision_member_eligible']) is bool and e['repeat_effect_used'] is False and type(e['other_original_groups_same_candidate']) is list and type(e['facts']) is dict,'accepted membership evidence types, never effect selection')
  raw=json.dumps(e,ensure_ascii=False,separators=(',',':')).encode();stored+=len(raw);need(stored<=lim['max_stored_evidence_bytes'],'stored evidence byte cap')
  rows=members.setdefault(protein,[]);rows.append(deepcopy(e));need(len(rows)<=lim['max_matches_per_target'],'per-target evidence fanout HOLD')
 need(n==exp['evidence'],'full member view EOF count')
 collisions={};n=0
 for co in collision_records:
  n+=1;need(n<=lim['max_collision_rows'],'collision stream row cap');cn=text(co['species_name'],'collision native');protein=text(co['protein_id'],'collision protein')
  if cn!=native:continue
  need(protein in allowed and protein not in collisions,'duplicate/outside-target collision HOLD');groups=co['original_groups'];need(type(groups) is list and len(groups)>=2 and all(type(x) is str and x for x in groups) and len(set(groups))==len(groups),'complete actual collision group claims');collisions[protein]=deepcopy(co)
 need(n==exp['collisions'],'complete collision overlay EOF count')
 seen=set();previous=None;n=0
 for view in canonical_views:
  n+=1;need(n<=lim['max_targets'],'canonical stream cap');r=view['record'];protein=text(r['key']['protein_id'],'canonical protein')
  need(r['schema']=='canonical-source-coordinate-v1' and r['key']['assembly']==assembly and r['native_species_name']==native and r['primary_callable_or_model_eligible'] is False,'exact accepted canonical identity/source-only')
  need(protein in allowed and protein not in seen and (previous is None or previous<protein),'complete original canonical sorted unique target order');seen.add(protein);previous=protein
  need(type(view['source_flags']) is list and type(view['blocking_states']) is list and type(view['source_authority_refs']) is list and view['source_authority_refs'],'explicit inherited source flag/NA authority view')
  for authority in view['source_authority_refs']:ref(authority)
  blocks=deepcopy(view['blocking_states']);ev=members.get(protein,[]);group=None;signature=None;status='MEMBER_BINDING_UNRESOLVED'
  if not ev:blocks.append('MEMBER_EVIDENCE_UNOBSERVED')
  else:
   actual_claims={x['original_group'] for x in ev if x['member_supported']}
   if protein in collisions or len(actual_claims)>1:blocks.append('ACTUAL_SUPPORTED_CANDIDATE_COLLISION')
   elif len(ev)!=1:blocks.append('NONUNIQUE_MEMBER_EVIDENCE')
   elif ev[0]['new_OG']=='OG0000725':blocks.append('NATIVE_OG0000725_UNRESOLVED_EXCEPTION')
   elif not (ev[0]['member_supported'] and ev[0]['post_collision_member_eligible']):blocks.append('ACCEPTED_MEMBER_NOT_ELIGIBLE')
   else:
    group=ev[0]['original_group'];signature=deepcopy(sigs.get(group));status='MEMBER_SUPPORTED_UNIQUE_NO_COLLISION'
    if signature is None:blocks.append('OLD_GROUP_SIGNATURE_UNOBSERVED')
  coords={}
  if r['coordinate_status']==SUPPORTED:
   s=r['source_coordinates'];a,z,rs,re_,center=(s[k] for k in ('full_original_gene_start1','full_original_gene_end1','RNA_start1','RNA_end1','window_center0_boundary'))
   need(all(type(x) is int for x in (a,z,rs,re_,center)) and 1<=a<=rs<=re_<=z and s['strand'] in ('+','-') and center==(rs-1 if s['strand']=='+' else re_) and type(s['target_gene_length_bp']) is int and s['target_gene_length_bp']==z-a+1,'accepted fullRNA/fullgene separation, no recanonicalization')
   coords={'RNA_start1':rs,'RNA_end1':re_,'full_original_gene_start1':a,'full_original_gene_end1':z,'source_gene_id_raw':text(s['source_gene_ID'],'source gene').encode('utf-8'),'strand_raw':s['strand'].encode('ascii'),'canonical_seqid_raw':text(s['seqid'],'source seqid').encode('utf-8'),'window_center0':center,'source_RNA_ID':deepcopy(s['source_RNA_ID'])}
  outrefs={k:deepcopy(refs[k]) for k in ('member','canonical','report')}
  if 'genome_QA' in refs:outrefs['genome_QA']=deepcopy(refs['genome_QA'])
  binding={'schema':'ROOT_ACCEPTED_TARGET_GROUP_SIGNATURE_ROW_BINDING_V1','approved_by':'current-conversation-root','assembly':assembly,'native_species_name':native,'target_key':key(assembly,protein),'protein_id':protein,'old_group':group,'member_status':status,'coordinate_status':r['coordinate_status'],'source_flags':deepcopy(view['source_flags']),'blocking_states':blocks,'coordinates':coords,'signature':signature,'refs':outrefs}
  yield {'binding':binding,'lineage':{'source_authority_refs':deepcopy(view['source_authority_refs']),'NA_authority_ref':deepcopy(refs['NA_authority']),'evidence_ref':deepcopy(refs['evidence']),'collision_ref':deepcopy(refs['collision']),'signature_source_ref':sr,'canonical_record':deepcopy(r),'member_evidence':deepcopy(ev),'collision_record':deepcopy(collisions.get(protein)),'binding_is_Root_input_view_not_new_acceptance':True,'no_J_weights_S_tree_model_eligibility':True}}
 need(n==exp['canonical'] and seen==allowed,'complete canonical physical iterator EOF/membership; drafts unaccepted until here')
