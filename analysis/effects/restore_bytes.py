"""Pure per-species T5 adapter proposal; no I/O, backend import or fitting.
Root must admit complete immutable inputs and accept outputs independently.
Normalized records are an interface, never a claim that real schema was read.
"""
from copy import deepcopy
import math

LINEAGE = ('coordinate', 'all_exon', 'genome_mask')
VARIANTS = ('A50', 'K50', 'A25')

def need(ok, why):
    if not ok:
        raise ValueError('HOLD: ' + why)

def text(x):
    need(type(x) is str and bool(x), 'exact nonempty identity')
    return x

def finite(x):
    need(type(x) in (int, float) and math.isfinite(x), 'finite fixed numeric field')
    return x

def identity(r):
    """No coordinate/mask normalization; Root supplies hash-bound lineage tokens."""
    return (text(r['assembly']), text(r['protein']),
            *(text(r['lineage'][k]) for k in LINEAGE))

def member_key(r):
    return (text(r['old_group']), *identity(r))

def link_species(cells, targets, signatures, baseline_members, baseline_S, context):
    """Return outcome-blind rows/manifests; no coefficients or backend calls.

    cells: all groups for one scenario/view/species, qualified or typed unresolved.
    targets: normalized accepted measured-target catalog, with backend_row50,
      narrow_record, typed_NA and opaque measurement_ref. No response value used.
    baseline_members: exact old-J membership records and fixed baseline weights.
    signatures: complete global group-> {xA,xK}; never recomputed here.
    """
    need(context['accepted_complete_inputs'] is True and
         context['outcome_used_for_selection'] is False, 'Root complete blind input gate')
    need(context['baseline_mode'] == 'REUSE_NO_REJOIN_NO_REFIT', 'B0 is immutable reuse')
    need(context['scenario'] != 'B0', 'B0 never rebuilt by this adapter')
    groups = context['group_order']
    need(len(groups) == len(set(groups)) == context['expected_groups'] and
         set(signatures) == set(groups), 'complete original signature universe')
    need(context['expected_groups'] <= 1896 and context['synthetic'] is True or
         context['expected_groups'] == 1896 and context['synthetic'] is False,
         'explicit synthetic or production group count')
    need(len(cells) == len(groups) and [c['old_group'] for c in cells] == groups,
         'complete ordered cells, never omit unresolved groups')
    need(len(baseline_S) == len(set(baseline_S)) and
         (context['synthetic'] or len(baseline_S) == 35), 'fixed S preserved')
    catalog = {}
    for t in targets:
        k = (text(t['assembly']), text(t['protein']))
        need(k not in catalog, 'duplicate protein/assembly target')
        identity(t)
        need(type(t['typed_NA']) is list and type(t['measurement_ref']) is dict,
             'explicit typed NA and measurement lineage')
        catalog[k] = t
    need(len(catalog) == context['expected_species_targets'], 'full target catalog count')
    old = {}
    for r in baseline_members:
        k = member_key(r)
        need(k not in old, 'duplicate baseline member key')
        need(finite(r['weight']) > 0, 'positive frozen common baseline weight')
        old[k] = r
    linked = []; own_rows = []; stable = []; narrow = {}
    used = set()
    for c in cells:
        g = c['old_group']; sig = signatures[g]
        finite(sig['xA']); finite(sig['xK'])
        rec = {'old_group': g, 'mapping_status': text(c['status']),
               'typed_NA': deepcopy(c['typed_NA']), 'measurement_ref': None,
               'signature': deepcopy(sig), 'status': 'MAPPING_UNRESOLVED'}
        need(type(c['typed_NA']) is list, 'mapping typed states explicit')
        if c['status'] != 'QUALIFIED':
            linked.append(rec); continue
        k = (text(c['assembly']), text(c['protein']))
        rec.update(protein=c['protein'], assembly=c['assembly'],
                   lineage=deepcopy(c['lineage']))
        if k not in catalog:
            rec['status'] = 'TARGET_NOT_MEASURED'
            rec['typed_NA'].append('TARGET_NOT_MEASURED')
            linked.append(rec); continue
        t = catalog[k]
        if identity(c) != identity(t):
            rec['status'] = 'LINEAGE_MISMATCH'
            rec['typed_NA'].append('COORDINATE_OR_MASK_LINEAGE_MISMATCH')
            linked.append(rec); continue
        need(k not in used, 'qualified target collision must remain unresolved upstream')
        used.add(k)
        rec.update(status='MEASUREMENT_REFERENCE_LINKED', protein=c['protein'],
                   assembly=c['assembly'], lineage=deepcopy(c['lineage']),
                   measurement_ref=deepcopy(t['measurement_ref']))
        rec['typed_NA'].extend(deepcopy(t['typed_NA']))
        # Numeric admission depends only on original preoutcome preparation,
        # never on RM NA/value: repeat-source NA stays in the selected J.
        row = t['backend_row50']
        if row is not None:
            nr = deepcopy(row); nr.update(id=g, xA=sig['xA'], xK=sig['xK'])
            own_rows.append(nr)
        narrow[g] = deepcopy(t['narrow_record'])
        mk = member_key(rec)
        if mk in old:
            stable.append({'member_key': list(mk), 'id': g,
                           'baseline_weight': old[mk]['weight'],
                           'backend_row50': deepcopy(own_rows[-1]) if row is not None else None,
                           'typed_NA': deepcopy(rec['typed_NA']),
                           'measurement_ref': deepcopy(t['measurement_ref'])})
        linked.append(rec)
    return {'status': 'DRAFT_LINK_MANIFEST_NOT_COEFFICIENTS', 'cells': linked,
            'own_prepare_primary_rows': own_rows, 'narrow_records': narrow,
            'stable_baseline_candidates': stable, 'fixed_S_order': list(baseline_S),
            'own_J_next': 'one accepted prepare_primary on original preoutcome rows only',
            'stable_J_next': 'intersect eligible member keys across views with baseline J; fixed baseline weights, no prepare_primary reweight',
            'finite_versions': list(VARIANTS), 'B0_recomputed': False,
            'RM_or_NA_driven_member_selection': False, 'model_calls': 0}


def finite_manifest(species, fixed_S, per_species_states):
    """No silent S omission; exact fixed order even when sensitivity is NA."""
    need(len(fixed_S) == len(set(fixed_S)) and set(per_species_states) >= set(fixed_S),
         'all fixed-S sensitivity states required')
    return {'species': species, 'S_order': list(fixed_S), 'versions': list(VARIANTS),
            'S_states': [deepcopy(per_species_states[s]) for s in fixed_S],
            'unestimable_remains_NA': True, 'fit_calls': 0}


def restore_exact_bytes(value):
    """Reverse only json_safe's exact singleton byte tag; no ID/coercion rules."""
    if type(value) is dict:
        if set(value) == {'__bytes_hex__'}:
            h = value['__bytes_hex__']
            need(type(h) is str and len(h) % 2 == 0 and
                 all(c in '0123456789abcdef' for c in h), 'exact lowercase byte hex')
            return bytes.fromhex(h)
        return {k: restore_exact_bytes(v) for k, v in value.items()}
    if type(value) is list:
        return [restore_exact_bytes(v) for v in value]
    return deepcopy(value)

def normalize_accepted_panel(panel_object, membership, signatures, source_refs):
    """Normalize existing accepted preoutcome only, no repeat responses.

    Event-ref bundle is a conservative whole-target provenance identity, NOT
    separate coordinate/exon/mask content checksums. New mapping must explicitly
    carry the same bundle identity for reuse. Missing fields remain missing.
    source_refs supplies accepted panel, member evidence, full signature pins.
    """
    import hashlib
    import json
    need(source_refs['accepted_complete_inputs'] is True, 'accepted local source gate')
    for name in ('panel', 'membership', 'signatures'):
        r = source_refs[name]
        need(type(r) is dict and type(r.get('path')) is str and
             type(r.get('sha256')) is str and len(r['sha256']) == 64, 'source pin')
    p = restore_exact_bytes(panel_object)
    need(len(signatures) == source_refs['signature_groups'], 'complete frozen signatures')
    need(p['assembly'] == source_refs['assembly'] and
         p['native_species_name'] == source_refs['native_species_name'], 'species identity')
    groups = {}; proteins = set()
    for m in membership:
        g = text(m['old_group']); pr = text(m['protein'])
        need(g not in groups and pr not in proteins, 'unique old group/protein membership')
        need(g in signatures, 'group has frozen signature')
        groups[g] = pr; proteins.add(pr)
    need(len(groups) == source_refs['species_target_count'], 'complete species membership')
    rows = {}
    for row in p['literal_rows50']:
        need(row['id'] not in rows and row['id'] in groups, 'unique literal primary row')
        rows[row['id']] = row
    events = {}
    for e in p['source_event_ledger']:
        need(e.get('schema') == 'ROOT_LOSSLESS_TARGET_EVENT_REFERENCE_V1' and
             e.get('complete_event_retained') is True, 'lossless event refs required')
        g = e['old_group_from_cache']; need(g in groups, 'event group membership')
        need(e['radius'] in (None, 50000, 25000), 'original radius only')
        need(type(e.get('rawline_sha256')) is str and len(e['rawline_sha256']) == 64 and
             type(e.get('physical_record_index')) is int and e['physical_record_index'] > 0 and
             type(e.get('source_ref')) is dict, 'exact raw physical source binding')
        group_events = events.setdefault(g, {})
        need(e['radius'] not in group_events, 'unique radius event')
        group_events[e['radius']] = e
    ids = p['panel']['ids']; weights = p['panel']['weights']
    need(len(ids) == len(weights) and len(ids) == len(set(ids)) and set(ids) <= set(groups),
         'accepted J ids/weights closure')
    weights = dict(zip(ids, weights)); targets = []; baseline = []; missing = []
    for g, protein in groups.items():
        ev = events.get(g, {}); reasons = []
        if set(ev) not in ({None}, {50000, 25000}):
            missing.append({'old_group': g, 'reason': 'INCOMPLETE_TARGET_EVENT_BUNDLE'})
            continue
        keys = {text(e['target_key']) for e in ev.values()}
        need(len(keys) == 1, 'same target across both radii')
        # Hash exact accepted reference objects, not reserialized raw event bodies.
        refs = [ev[k] for k in sorted(ev, key=lambda x: -1 if x is None else x)]
        token = 'WHOLE_ACCEPTED_TARGET_EVENT_REF:' + hashlib.sha256(
            json.dumps(refs, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        lineage = {k: token for k in LINEAGE}
        for e in refs:
            if e['row_preparation_status'] == 'TYPED_NA':
                reasons.append({'status': 'TYPED_NA', 'exact_reasons': 'NOT_IN_REFERENCE_ONLY_LEDGER',
                                'original_event_ref': deepcopy(e)})
        narrow = p['literal_narrow_records'].get(g)
        if narrow is None:
            reasons.append('NARROW_RECORD_MISSING')
        row = rows.get(g)
        if row is not None:
            sig = signatures[g]
            need(row['xA'] == sig['xA'] and row['xK'] == sig['xK'], 'frozen signature equality')
        t = {'old_group': g, 'assembly': p['assembly'], 'protein': protein,
             'target_key': next(iter(keys)), 'lineage': lineage,
             'lineage_basis': 'WHOLE_ACCEPTED_TARGET_EVENT_REF_NOT_SEPARATE_MASK_CHECKSUMS',
             'typed_NA': reasons, 'backend_row50': deepcopy(row),
             'narrow_record': deepcopy(narrow),
             'measurement_ref': {'status': 'PREOUTCOME_EVENT_REF_ONLY_RM_REFERENCE_NOT_SUPPLIED',
                                 'events': deepcopy(refs), 'panel_ref': deepcopy(source_refs['panel'])}}
        targets.append(t)
        if g in weights:
            need(row is not None, 'accepted J must have original primary row')
            need(finite(weights[g]) > 0, 'accepted baseline weight')
            baseline.append({**{k: deepcopy(t[k]) for k in ('old_group','assembly','protein','lineage')},
                             'weight': weights[g]})
    return {'status': 'NORMALIZED_PREOUTCOME_DRAFT_WITH_EXPLICIT_GAPS',
            'targets': targets, 'baseline_members': baseline, 'missing': missing,
            'RM_response_reads': 0, 'backend_calls': 0,
            'release_gaps': ['Root bind exact RM measurement refs and preserve typed RM NA',
                             'Root recover exact typed-NA reason bodies from accepted rows only if needed',
                             'Root mapping must carry same whole accepted bundle token'],
            'source_refs': deepcopy(source_refs)}
