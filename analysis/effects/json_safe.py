"""Pure once-per-species preoutcome wrapper; no IO/imported backend/fits.
Root admits event/code/runtime refs, physical EOF/hash/once/resources externally.
All returned objects are drafts until independent Root acceptance.
"""
from copy import deepcopy
import math
import json
import re

LAYERS = ('non-Bilateria','Spiralia','Ecdysozoa','Deuterostomia','Xenacoelomorpha_separate')
ROLES = {'new_normal_candidate':45,'old_discovery':3,'outgroup_background':3,
         'Hemicordylus_method_exception':1,'gene_tree_only':2}
READY_PROTECTION = 'LEGACY_COORDINATE_EXCLUSION_NOT_APPLICABLE_DISTINCT_ASSEMBLY'


def need(ok, why):
    if not ok:
        raise ValueError('HOLD: '+why)


def json_safe(value):
    """Exact finite JSON/NumPy conversion; preserve bytes with frozen hex tag.
    No NaN/inf->null, str(bytes), ID normalization or arbitrary object stringify.
    Caller still owns actual packed-byte/output/RSS/deadline guards.
    """
    active = set();nodes = [0]
    def cv(v, depth):
        nodes[0] += 1
        need(nodes[0] <= 2000000 and depth <= 64, 'JSON node/depth bound')
        t = type(v)
        if v is None or t in (str,bool):return v
        if t is bytes:return {'__bytes_hex__':v.hex()}
        if t is int:
            need(v.bit_length() <= 4096, 'JSON integer bound');return v
        if t is float:
            need(math.isfinite(v), 'nonfinite value must remain an error, never substituted');return v
        if t in (dict,list,tuple):
            need(id(v) not in active, 'cyclic JSON structure')
            active.add(id(v))
            try:
                if t is dict:
                    need(all(type(k) is str for k in v), 'literal string JSON keys')
                    need(set(v) != {'__bytes_hex__'}, 'literal dict collides with frozen byte tag')
                    return {k:cv(z,depth+1) for k,z in v.items()}
                return [cv(z,depth+1) for z in v]
            finally:active.remove(id(v))
        # No NumPy import or actual numerical operation. Accept only NumPy's
        # literal container/scalar conversion, recurse to check finite values.
        if t.__module__.split('.')[0] == 'numpy':
            need(id(v) not in active, 'cyclic NumPy conversion')
            active.add(id(v))
            try:
                if callable(getattr(v,'tolist',None)):return cv(v.tolist(),depth+1)
                if callable(getattr(v,'item',None)):return cv(v.item(),depth+1)
            finally:active.remove(id(v))
        raise ValueError('HOLD: unsupported literal JSON/NumPy/bytes value')
    return cv(value,0)


def _representable_narrow(row):
    if type(row) is not dict or not {'complete','acgt_fraction','nonexonic_acgt','structure'} <= set(row):return False
    if type(row['complete']) is not bool:return False
    if type(row['acgt_fraction']) not in (int,float) or not math.isfinite(row['acgt_fraction']):return False
    if type(row['nonexonic_acgt']) is not int:return False
    q = row['structure']
    return type(q) in (list,tuple) and len(q)==5 and all(type(v) in (int,float) and math.isfinite(v) for v in q)


def _source_ref(r):
    need(type(r) is dict and type(r.get('path')) is str and r['path'].startswith('/')
         and type(r.get('sha256')) is str and re.fullmatch('[0-9a-f]{64}',r['sha256']) is not None,
         'explicit immutable accepted source reference')


def _reference_policy(retain_event, policy):
    if retain_event is None:
        need(policy is None, 'reference policy cannot be silently ignored without callback')
        return None
    need(callable(retain_event) and type(policy) is dict and
         policy.get('schema')=='ROOT_LOSSLESS_EVENT_REFERENCE_POLICY_V1' and
         policy.get('approved_by')=='current-conversation-root' and
         policy.get('mode')=='IMMUTABLE_SOURCE_PHYSICAL_LINE_SHA256' and
         policy.get('complete_source_retained_immutable') is True and
         policy.get('all_target_events_without_downsampling') is True,
         'explicit Root lossless immutable-source retention policy')
    _source_ref(policy['source_ref']);_source_ref(policy['source_acceptance_ref'])
    return deepcopy(policy)


def _retained_reference(event, record_index, retain_event, policy):
    # Snapshot compact identity BEFORE callback; no full-payload deepcopy needed.
    d=event['draft']
    expected={'event_kind':event['kind'],'target_key':d['target_key'],'radius':d['radius'],
        'row_preparation_status':d['row_preparation']['status'],
        'old_group_from_cache':event['old_group_from_cache'],'source_status':event['source_status'],
        'protection_status':d['protected_disposition']['status']}
    r=retain_event(event,record_index)
    need(type(r) is dict and r.get('schema')=='ROOT_LOSSLESS_TARGET_EVENT_REFERENCE_V1'
         and r.get('complete_event_retained') is True, 'callback must retain every complete target event')
    need(type(r.get('physical_record_index')) is int and r['physical_record_index']==record_index,
         'exact original1-based physical event record index')
    _source_ref(r['source_ref'])
    need(r['source_ref']==policy['source_ref'], 'callback immutable source reference must equal Root policy source')
    h=r.get('rawline_sha256')
    need(type(h) is str and re.fullmatch('[0-9a-f]{64}',h) is not None,
         'explicit original raw-line SHA256, not reserialized JSON digest')
    need(all(k in r and type(r[k]) is type(v) and r[k]==v for k,v in expected.items()),
         'retained reference exact target/radius/status/group/source/protection identity')
    # Recheck the compact input identity so callback cannot mutate selection data.
    need(event['kind']==expected['event_kind'] and d['target_key']==expected['target_key']
         and d['radius']==expected['radius'] and d['row_preparation']['status']==expected['row_preparation_status']
         and event['old_group_from_cache']==expected['old_group_from_cache']
         and event['source_status']==expected['source_status']
         and d['protected_disposition']['status']==expected['protection_status'],
         'retention callback must not mutate admitted event identity')
    # This is a bounded metadata receipt, never a disguised full body. Root still
    # owns actual original gzip/line hash/EOF verification and immutable retention.
    raw=json.dumps(r,allow_nan=False,ensure_ascii=False,separators=(',',':')).encode()
    need(len(raw)<=16384, 'compact lossless event-reference receipt cap16KiB')
    return deepcopy(r)


def prepare_species(iter_of_Root_admitted_events, expected_source_counts, backend, check, *,
                    retain_event=None, event_reference_policy=None):
    """Consume one complete admitted event view and call prepare_primary ONCE.
    expected_source_counts preserves all caller metadata, including actual pinned
    row-source refs. Required targets/windows/typed_NA are source counts: window
    TYPED_NA is additional to source typed_NA and is NEVER silently discarded.
    """
    check()
    ledger_policy=_reference_policy(retain_event,event_reference_policy)
    all_na_indices=[]
    need(type(expected_source_counts) is dict, 'explicit Root source count/metadata object')
    exp = {k:expected_source_counts[k] for k in ('targets','windows','typed_NA')}
    need(all(type(v) is int and v>=0 for v in exp.values()) and 0<exp['targets']<=2000
         and exp['windows']==2*(exp['targets']-exp['typed_NA']), 'exact finite source target/window/NA accounting')
    states = {};state_events = [];source_ledger = [];rows50 = [];narrow_by_group = {}
    targets = {};group_targets = {};seen_windows = set();terminal = None;identity = None
    source_na = windows = all_na = record_index = 0
    max_events = 1896+exp['windows']+exp['typed_NA']+1
    for event in iter_of_Root_admitted_events:
        check();record_index += 1
        need(record_index<=max_events and terminal is None and type(event) is dict,
             'bounded complete stream; no trailing event after terminal')
        kind = event['kind']
        if kind == 'FULL_GROUP_SOURCE_STATE':
            need(not source_ledger and len(states)<1896, 'full states precede target drafts')
            native,assembly = event['native_species_name'],event['assembly']
            role,layer = event['study_role'],event['major_layer']
            need(role in ROLES and role!='old_discovery', 'new51 only; original3 reuse, never rerun')
            current=(native,assembly,role,layer)
            need(identity is None or identity==current, 'same exact native/assembly/fixed role/layer')
            identity=current;g=event['old_group']
            need(type(g) is str and g and g not in states and type(event['source_status']) is str,
                 'one original source state per group')
            states[g]=event['source_status'];state_events.append(deepcopy(event))
            continue
        need(len(states)==1896 and identity is not None, 'complete1896 states before target/closure')
        if kind == 'SPECIES_DRAFT_COMPLETE':
            need(event['native_species_name']==identity[0] and event['assembly']==identity[1]
                 and event['full_group_states']==1896 and event['target_window_NA_counts']==exp
                 and event['target_drafts']==exp['windows']+exp['typed_NA'], 'literal source completion counts/identity')
            terminal=deepcopy(event);continue
        need(kind=='TARGET_ROW_DRAFT', 'known frozen event family only')
        d=event['draft'];m=d['row_preparation'];target=d['target_key'];radius=d['radius'];g=event['old_group_from_cache']
        need(type(target) is str and re.fullmatch('t[0-9a-f]{64}',target) is not None,
             'existing exact hashed event target ID, never regenerated')
        need(g in states and event['source_status']==states[g], 'target links to exact full-group source state')
        need(g not in group_targets or group_targets[g]==target, 'no silent target chooser within a group')
        group_targets[g]=target
        need(target not in targets or targets[target]['group']==g, 'same target/group binding across windows')
        if target not in targets:targets[target]={'group':g,'radii':set()}
        need(radius in (None,50000,25000) and (target,radius) not in seen_windows,
             'one literal event per exact target/radius')
        seen_windows.add((target,radius));targets[target]['radii'].add(radius)
        need(type(d['protected_disposition']) is dict and type(d['protected_disposition'].get('status')) is str,
             'explicit source protection record, never a false/default')
        status=m['status'];need(status in ('TYPED_NA','BACKEND_ROW_PREPARED_NOT_ELIGIBLE'), 'frozen row draft status')
        retained = (deepcopy(event) if retain_event is None else
                    _retained_reference(event,record_index,retain_event,ledger_policy))
        source_ledger.append(retained)
        source_index=len(source_ledger)-1
        if status=='TYPED_NA':
            all_na += 1;all_na_indices.append(source_index)
        if radius is None:
            need(status=='TYPED_NA', 'source NA cannot be a fabricated numeric window')
            source_na += 1
        else:
            windows += 1
            if status=='BACKEND_ROW_PREPARED_NOT_ELIGIBLE':
                need(d['protected_disposition']['status']==READY_PROTECTION,
                     'Root protection status must be ready before numeric admission; no default NotProtected')
            if radius==50000 and status=='BACKEND_ROW_PREPARED_NOT_ELIGIBLE':
                row=m['backend_row'];need(type(row) is dict and row['id']==g,
                                          'unchanged original-group50 row identity')
                rows50.append(deepcopy(row))
            if radius==25000:
                need(g not in narrow_by_group, 'duplicate narrow group metadata')
                narrow_by_group[g]={'row':deepcopy(m['narrow_row']), 'status':status,
                    'reasons':deepcopy(m['reasons']), 'target_key':target,
                    'source_event_ledger_index':source_index,'source_event_record_index':record_index}
    check()
    need(terminal is not None and len(states)==1896, 'complete source iterator EOF/terminal required')
    observed={'targets':len(targets),'windows':windows,'typed_NA':source_na}
    need(observed==exp and len(source_ledger)==windows+source_na, 'all source target/window/NA events conserved')
    need(all(t['radii'] in ({None},{50000,25000}) for t in targets.values()),
         'both original source windows or one honest source NA; no missing pair')
    # All input views/counts close before any numerical call. The accepted backend
    # handles exact geometry/callability/local exclusions, weights and information.
    panel=backend.prepare_primary(rows50)
    check()
    ids=panel['ids'];need(len(set(ids))==len(ids) and set(ids)<=set(r['id'] for r in rows50),
                         'accepted backend preserves original row identity/subset')
    missing=[];narrow={}
    for g in ids:
        record=narrow_by_group.get(g)
        if record is None or record['status']!='BACKEND_ROW_PREPARED_NOT_ELIGIBLE' or not _representable_narrow(record['row']):
            missing.append({'id':g,'target_key':group_targets[g],
                'reasons':deepcopy(record['reasons']) if record is not None else ['NARROW_METADATA_UNOBSERVED'],
                'source_event_ledger_index':record['source_event_ledger_index'] if record is not None else None,
                'source_event_record_index':record['source_event_record_index'] if record is not None else None})
        else:narrow[g]=deepcopy(record['row'])
    if missing:
        a25={'status':'NARROW_TECHNICAL_NA_ADAPTATION_NEEDED','unrepresentable':missing,
             'representable_narrow_rows':narrow,'main_J_unchanged':True,'main_weights_unchanged':True,
             'variant_or_information_calls':0,'no_fake_structure_or_zero_GC':True}
    else:
        nr,nw,nx,excluded=backend._variant(panel,'A25',narrow)
        check();info=backend._information(nr,nw,nx);check()
        a25={'status':'PREPARED_FROM_EXACT_ACCEPTED_HELPERS','rows':nr,'weights':nw,'xA':nx,
             'narrow_rows':narrow,'excluded':excluded,'information':info,
             'main_J_unchanged':True,'main_weights_unchanged':True}
    return {'status':'PREOUTCOME_SPECIES_J_DRAFT_NOT_ROOT_ACCEPTANCE',
        'native_species_name':identity[0],'assembly':identity[1],
        'fixed_study_role':identity[2],'fixed_major_layer':identity[3],
        'Root_input_metadata':deepcopy(expected_source_counts),'full_group_states':state_events,
        'source_event_ledger':source_ledger,
        'source_event_ledger_schema':{'schema':'TARGET_EVENT_BODY_OR_LOSSLESS_REFERENCE_LEDGER_V1',
            'mode':'FULL_TARGET_EVENT_BODIES' if retain_event is None else 'LOSSLESS_IMMUTABLE_SOURCE_REFERENCES',
            'retained_target_events':len(source_ledger),'expected_target_events':exp['windows']+exp['typed_NA'],
            'reference_policy':deepcopy(ledger_policy),
            'rawline_sha_or_source_hash_verified_by_helper':False},
        'source_completion':terminal,
        'observed_source_counts':observed,'all_typed_NA_count':all_na,
        'all_typed_NA_ledger_indices':all_na_indices,
        'literal_rows50':rows50,'literal_narrow_records':narrow_by_group,
        'panel':panel,'A25':a25,'actual_source_EOF_or_hash_verified_by_helper':False,
        'prepare_primary_calls':1,'weights_recomputed_separately':False,'fit_or_outcome_calls':0}


def select_fixed_S(full54_order, fixed_roles, prepared_by_species):
    """Fixed-role/mainA50 structural S only. No response/beta/vJ or model call."""
    need(type(full54_order) is list and len(full54_order)==54 and len(set(full54_order))==54
         and set(fixed_roles)==set(full54_order), 'exact full54 order and frozen role roster')
    counts={k:0 for k in ROLES};new=set()
    for name in full54_order:
        r=fixed_roles[name]
        need(r['native_species_name']==name and r['fixed_study_role'] in ROLES, 'literal fixed role identity')
        counts[r['fixed_study_role']] += 1
        if r['fixed_study_role']!='old_discovery':new.add(name)
    need(counts==ROLES and set(prepared_by_species)==new, 'all51 new panels, fixed45/3/3/1/2 roles; old3 reused')
    S=[];ledger=[];layer_counts={k:0 for k in LAYERS}
    for name in full54_order:
        r=fixed_roles[name];role=r['fixed_study_role'];layer=r['fixed_major_layer']
        if role=='old_discovery':
            ledger.append({'species':name,'status':'ORIGINAL_DISCOVERY_REUSE_ONLY','role':role});continue
        p=prepared_by_species[name]
        need(p['native_species_name']==name and p['assembly']==r['assembly'], 'same exact prepared/frozen species assembly')
        info=p['panel']['information']['A50']
        if role!='new_normal_candidate':
            ledger.append({'species':name,'status':'FIXED_SEPARATE_ROLE_NOT_WORKING_S','role':role,
                           'A50_information':deepcopy(info)});continue
        need(layer in LAYERS, 'frozen normal-animal major layer')
        if info.get('spatial_eligible') is True:
            need(info.get('status')=='OK', 'accepted primaryA50 spatial information')
            S.append(name);layer_counts[layer] += 1
            ledger.append({'species':name,'status':'PRIMARY_A50_PREOUTCOME_SPATIAL_S','role':role,'layer':layer})
        else:ledger.append({'species':name,'status':'PRIMARY_A50_SPATIAL_INELIGIBLE','role':role,
                            'layer':layer,'A50_information':deepcopy(info)})
    reasons=[]
    if len(S)<20:reasons.append('FEWER_THAN_20_NEW_NORMAL_ANIMALS')
    if sum(layer_counts[k]>=3 for k in LAYERS[:4])<3:reasons.append('FEWER_THAN_3_FIRST_FOUR_LAYERS_WITH_AT_LEAST_3')
    return {'status':'FIXED_S_DRAFT_NOT_ROOT_ACCEPTANCE','full54_order':deepcopy(full54_order),
        'S_order':S,'layer_counts':layer_counts,'working_cohort_eligible':not reasons,
        'working_exclusions':reasons,'fixed_role_counts':counts,'full54_role_and_spatial_ledger':ledger,
        'working_record_identity_template':[{'species':name,'role':'new_normal_metazoan',
            'layer':fixed_roles[name]['fixed_major_layer'],'spatial_eligible':True} for name in S],
        'same_S_for_K50_A25_C_top':True,'beta_vJ_or_working_fit_assigned':False,
        'sensitivity_ineligibility_must_not_reselect_S':True}
