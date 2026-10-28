"""Distinct pure legacy32 adapter; Root injects already pinned base helper module.
No IO, DB, source/QC/mask rerun, report-owner selection or scientific assignment.
Input objects must be independently admitted immutable ORIGINAL receiver, manifest
and namespace metadata, never manufactured from the row binding under test.
"""
from copy import deepcopy
import hashlib,json
TOP={'schema','assembly','windows','counts','geometry_statuses','actual_sequence_records','actual_assigned_chromosomes','not_yet_global_J_weights_S_tree_freeze','RepeatMasker_or_models'}
WINDOW={'target_key','radius','source_geometry_ACGT_D_pass','actual_realchr','actual_report_record_index','audit','final_J_primary_model_eligibility_not_frozen'}
ASSEMBLY='GCA_051530055.1'

def exact(a,b):
    """Exact typed equality: bool/int and int/float are not interchangeable."""
    if type(a) is not type(b):return False
    if type(a) is dict:return set(a)==set(b) and all(exact(a[k],b[k]) for k in a)
    if type(a) is list:return len(a)==len(b) and all(exact(x,y) for x,y in zip(a,b))
    return a==b

def encoded(x):return json.dumps(x,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False).encode()

def legacy32_view(receiver_qa, original_manifest, namespace, p, *, base, check=lambda:None):
    """Return base.DraftView with existing driver facts()/typed_na()/EOF contract.
    Root also keeps returned .namespace_summary + per-fact namespace observation.
    Admitted physical FILE EOF/SHA/stats and original lineage are external duties;
    semantic SHA below binds original manifest with frozen producer convention.
    Source rows are drafts until complete external Root output/EOF acceptance.
    """
    need=base.need
    need(p['schema']=='ROOT_LEGACY32_FACT_VIEW_PROVENANCE_V1' and p['kind']==base.LEGACY
         and p['packet']==32 and p['actual_projection_released_by_this_adapter'] is False,
         'distinct explicit legacy32 projection provenance')
    for key in ('catalogue_ref','catalogue_acceptance_ref','canonical_acceptance_ref',
                'source_disposition_ref','code_freeze_ref','base_adapter_code_ref',
                'receiver_QA_ref','first32_Root_acceptance_ref','namespace_QC_ref',
                'namespace_schema_observation_ref','mask_acceptance_ref',
                'original_Genome_manifest_ref','original_manifest_lineage_ref'):
        base.ref(p[key])
    need(p.get('scientific_artifact_ref') is None, 'do not invent modern raw descriptor for legacy receiver route')
    m=original_manifest;q=receiver_qa;n=namespace
    # Finite metadata envelope only; no physical file or source proof inferred.
    check();need(sum(len(encoded(x)) for x in (m,q,n))<=8*2**20,'prior8MiB input metadata envelope')
    need(m['schema']=='GENOME_TARGET_WINDOW_INPUT_V1' and m['assembly']==ASSEMBLY,
         'original legacy32 manifest exact family/assembly')
    mh=hashlib.sha256(encoded(m)).hexdigest()
    need(mh==q['manifest_sha256']==p['original_manifest_semantic_sha256'],
         'independently pinned original manifest semantic binding')
    need(q['schema']=='GENOME_WINDOW_RECEIVER_V1' and q['status']=='SOURCE_GEOMETRY_AND_COUNT_ACCOUNTING_ACCEPTED',
         'original accepted receiver family')
    for key in ('primary_callable_or_model_eligible','real_chromosome_assigned','repeat_outcomes_computed','weights_or_model_computed'):
        need(q[key] is False,'receiver assignment protections unchanged')
    supported,nas=base.targets_index(m['targets'])
    need(len(m['targets'])==len(supported)==322 and nas==[], 'original322supported/0typedNA unchanged')
    need(len(q['window_audits'])==644 and q['typed_NA']==[]
         and type(q['terminal']['targets']) is int and q['terminal']['targets']==322
         and type(q['terminal']['windows']) is int and q['terminal']['windows']==644,
         'full legacy receiver target/window/typedNA closure')
    need(q['terminal']['input_bytes']==m['expected_source']['bytes'], 'original full-source byte accounting')
    need(type(n) is dict and set(n)==TOP and type(n['schema']) is str and n['schema']
         and n['assembly']==ASSEMBLY and type(n['windows']) is list and len(n['windows'])==644,
         'exact admitted namespace schema/assembly/window shape')
    need(type(n['not_yet_global_J_weights_S_tree_freeze']) is bool and n['not_yet_global_J_weights_S_tree_freeze'] is True
         and type(n['RepeatMasker_or_models']) is int and n['RepeatMasker_or_models']==0,
         'historical namespace unfrozen and noRM/model observation retained')
    need(type(n['actual_sequence_records']) is int and n['actual_sequence_records']==len(q['sequences'])==q['terminal']['sequences'],
         'existing complete sequence receipt count equality only')
    need(type(n['actual_assigned_chromosomes']) is int and n['actual_assigned_chromosomes']>=0
         and type(n['counts']) is dict and type(n['geometry_statuses']) is dict,
         'literal namespace diagnostic summary types; no recount/reinterpretation')
    def facts():
        i=0
        for t in supported:
            need(t['assembly']==ASSEMBLY and [w['radius'] for w in t['windows']]==[50000,25000],
                 'original full-target frozen window order')
            for window in t['windows']:
                check();a=q['window_audits'][i];old=n['windows'][i];i+=1
                need(type(old) is dict and set(old)==WINDOW,'exact original namespace window wrapper')
                need(type(old['radius']) is int and old['target_key']==t['target_key']
                     and old['radius']==window['radius'],'namespace exact target/radius/order')
                need(exact(old['audit'],a),'all existing namespace/receiver audit values and types identical')
                need(type(old['source_geometry_ACGT_D_pass']) is bool and type(old['actual_realchr']) is bool
                     and type(old['actual_report_record_index']) is int and old['actual_report_record_index']>=0
                     and old['final_J_primary_model_eligibility_not_frozen'] is True,
                     'existing namespace diagnostics/owner-index/unfrozen types')
                f=base.from_compact_audit(a,t,window,i)
                # Preserve historical owner/realchr/pass as unpromoted observations.
                # Root must independently resolve/bind report owner and current role
                # and protected disposition; these fields cannot assign eligibility.
                f['legacy32_namespace_observation']=deepcopy(old)
                yield f
    result=base.DraftView(facts(),iter(()),644,0,deepcopy(p),check)
    result.namespace_summary=deepcopy({k:v for k,v in n.items() if k!='windows'})
    result.provenance['observed_namespace_schema_literal']=n['schema']
    return result
