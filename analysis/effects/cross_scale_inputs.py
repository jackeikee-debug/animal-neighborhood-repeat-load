"""Pure phase2 input preparation. No IO, fits, dispatch, source or cohort choice.

Root calls the unchanged backend.spatial_jackknife only after independent
acceptance, mathematical deduplication and exact once/resource/source release.
"""
from copy import deepcopy
import hashlib
import json
import math
import struct

BACKEND_SHA = '36a901dce504df25c8e2d8cc0f9df2810b41ceb004d4d64b2d53e149699f0c26'
RUNTIME_ID = 'Py3.12.14/NumPy2.3.5/SciPy1.16.3/one_BLAS_thread'
CACHE_SCHEMA = 'ROOT_COMPACT_PHASE2_ACCEPTED_BASELINE_AND_PAIRED_OBSERVATIONS_V1'
READY_A25 = 'PREPARED_FROM_EXACT_ACCEPTED_HELPERS'


def need(ok, reason):
    if not ok: raise ValueError('HOLD: '+reason)


def number(x):
    need(type(x) in (int, float) and math.isfinite(x), 'finite literal numeric value')
    return x


def _wire(value):
    """Lossless local signature codec; byte identity never converted to text."""
    if value is None: return ['none']
    if type(value) is bool: return ['bool', value]
    if type(value) is int: return ['int', str(value)]
    if type(value) is float:
        number(value); return ['float64', struct.pack('>d', value).hex()]
    if type(value) is str: return ['str', value]
    if type(value) is bytes: return ['bytes', value.hex()]
    if type(value) in (list, tuple): return ['sequence', [_wire(v) for v in value]]
    if type(value) is dict:
        pairs = [[_wire(k), _wire(v)] for k, v in value.items()]
        pairs.sort(key=lambda p: json.dumps(p[0], ensure_ascii=True, separators=(',', ':')))
        return ['mapping', pairs]
    raise ValueError('HOLD: unsupported signature value; no stringify')


def _f64(x):
    return struct.pack('>d', float(number(x))).hex()


def consumer_signature(effective, *, backend_sha=BACKEND_SHA, runtime_id=RUNTIME_ID):
    """Same routine accepts NEW and reconstructed old/PC effective dictionaries.

    Requires rows/weights/x/observations + effective A50/K50/A25. Logical
    discovery/target labels, original schema labels and unused proxy x fields do
    not distinguish a numerical call. Preserve ordered ID/location, consumed
    structure50, D, weights, actual x and actual float64 target fractions.
    Exact species/assembly identity is bound; reuse is within that unit only.
    D inclusion is conservative: equal signatures suffice for reuse; differing
    signatures do not prove that final ranked/regression inputs differ.
    """
    variant = effective['variant']; need(variant in ('A50', 'K50', 'A25'), 'consumer variant')
    rows, weights, x, observations = (effective[k] for k in ('rows','weights','x','observations'))
    need(len(rows) == len(weights) == len(x), 'effective array lengths')
    need(len({r['id'] for r in rows}) == len(rows), 'effective unique IDs')
    need(set(observations) == {r['id'] for r in rows}, 'effective observation closure')
    arrays = []
    for r, w, actual_x in zip(rows, weights, x):
        o = observations[r['id']]
        need(o.get('status') in (None, 'MEASURED'), 'signature has no typed NA')
        A, U, D = o['A'], o['U'], r['nonexonic_acgt']
        need(all(type(v) is int for v in (A,U,D)) and 0 <= U <= A <= D and D > 0,
             'integer outcome denominator/bounds')
        need(len(r['structure']) == 5 and number(w) > 0, 'effective structure/weights')
        numerator = A-U if variant == 'K50' else A
        arrays.append([r['id'], r['location'], [_f64(v) for v in r['structure']],
                       D, _f64(w), _f64(actual_x), _f64(numerator/D)])
    identity = ['CONSUMER_SUFFICIENT_SPATIAL_INPUT_V1', backend_sha, runtime_id,
                effective['species'], effective['assembly'], arrays]
    payload = json.dumps(_wire(identity), ensure_ascii=True, separators=(',', ':')).encode('ascii')
    return {'schema': 'CONSUMER_SUFFICIENT_SPATIAL_INPUT_V1',
            'sha256': hashlib.sha256(payload).hexdigest(), 'backend_sha256':backend_sha,
            'runtime_identity':runtime_id, 'ordered_rows':len(rows),
            'logical_labels_excluded':True, 'species_assembly_identity_bound':True,
            'same_hash_sufficient_not_necessary':True}


def _base(cache):
    need(cache['schema'] == CACHE_SCHEMA and cache['unchanged_J_weights_S'] is True,
         'accepted compact phase2 schema/frozen baseline')
    base = cache['base']; p = base['panel']; ids = p['ids']; rows = p['rows']; weights = p['weights']
    need(type(ids) in (list, tuple) and len(set(ids)) == len(ids), 'unique original J')
    need([r['id'] for r in rows] == list(ids) and len(weights) == len(ids), 'original J arrays/order')
    need(all(number(w) > 0 for w in weights), 'original positive weights')
    need((cache['native_species_name'],cache['assembly']) == (base['native_species_name'],base['assembly']),
         'exact species/assembly')
    return base, p


def _result(cache, panel, observations, labels, exclusions):
    ids = panel['ids']; rows = panel['rows']; variant = labels['backend_variant']
    need(set(observations) == set(ids), 'exact selected observations including NA')
    bad = []
    for r in rows:
        o = observations[r['id']]
        if o.get('status') == 'TYPED_NA':
            need(bool(o.get('reason')), 'typed NA reason preserved'); bad.append(r['id'])
        else:
            need(o.get('status') in (None,'MEASURED'), 'known observation state')
            A,U,D = o['A'],o['U'],r['nonexonic_acgt']
            need(all(type(v) is int for v in (A,U,D)) and 0 <= U <= A <= D and D > 0,
                 'actual selected response denominator')
    x = [r['xK' if variant == 'K50' else 'xA'] for r in rows]
    need(all(math.isfinite(number(v)) for v in x), 'unchanged finite discovery vector')
    effective = {'species':cache['native_species_name'], 'assembly':cache['assembly'],
                 'variant':variant, 'rows':deepcopy(rows), 'weights':list(panel['weights']),
                 'x':list(x), 'observations':deepcopy(observations)}
    ready = bool(ids) and not bad
    return {'schema':'PHASE2_PURE_SPATIAL_ADAPTER_INPUT_V1',
            'status':'PREPARED_NUMERIC_INPUT' if ready else 'TYPED_NA_NO_NUMERIC_CALL',
            'reason':None if ready else 'EMPTY_ACTIVE_SET' if not ids else 'FROZEN_OBSERVATIONS_HAVE_TYPED_NA',
            'minimal_panel':panel, 'observations':deepcopy(observations),
            'effective':effective, 'logical_labels':labels, 'typed_NA_ids':bad,
            'excluded':deepcopy(exclusions), 'selection_excluded':deepcopy(exclusions),
            'backend_expected_excluded':[], 'source_lineage':deepcopy(cache['source_lineage']),
            'consumer_signature':consumer_signature(effective) if ready else None,
            'backend_call':{'function':'spatial_jackknife','variant':variant,'narrow_rows':None},
            'new_numeric_calls':0, 'primary_reweighting_calls':0,
            'A25_requalification_calls':0, 'diagonal_calls':0,
            'Root_exact_release_and_global_dedup_required':True}


def prepare_offdiagonal(cache, *, discovery, target):
    """Only xA50->K50 and xK50->A50, same full B0 J/weights/structure.

    Returns input only; Root supplies minimal_panel/observations and recorded
    backend variant to unchanged spatial_jackknife. Never execute diagonals.
    """
    need(discovery in ('xA50','xK50') and target in ('A50','K50'), 'finite 2x2 labels')
    need((discovery,target) in (('xA50','K50'),('xK50','A50')), 'diagonal reuse only; no diagonal adapter call')
    base,p = _base(cache)
    rows = deepcopy(p['rows']); frozen_field = 'xA' if discovery == 'xA50' else 'xK'
    consumed_field = 'xK' if target == 'K50' else 'xA'
    # Preserve the other frozen field; only the actual consumer field is proxied.
    for source,row in zip(p['rows'],rows): row[consumed_field] = source[frozen_field]
    panel = {'ids':list(p['ids']), 'rows':rows, 'weights':list(p['weights'])}
    observations = cache['observations']['50000']
    need(set(observations) == set(p['ids']), 'full B0 original observation closure')
    labels = {'analysis':'POSTHOC_DISCOVERY_TARGET_2X2_OFFDIAGONAL',
              'discovery':discovery, 'target':target, 'response_radius':50000,
              'backend_variant':target, 'discovery_source_field':frozen_field,
              'proxied_consumer_field':consumed_field, 'J':'ORIGINAL_B0_FULL_J',
              'weights':'ORIGINAL_B0_UNCHANGED', 'structure':'ORIGINAL_50',
              'preprocessing':'ORIGINAL_BACKEND_RECOMPUTES_SAME_ROWS_RANKS_PER_FULL_DELETION',
              'original_diagonal_reused_not_called':True}
    return _result(cache,panel,observations,labels,p.get('excluded',[]))


def prepare_fixed_common_scale(cache, *, response_radius):
    """A25-qualified original IDs, old weights; xA50 and structure50 both scales.

    This is NOT original A25 (which uses narrow structure). Use backend A50 for
    both proxies. Only D changes at radius25; no A25 _variant/requalification,
    overlap recomputation or signature restandardization. No K25.
    """
    need(response_radius in (50000,25000), 'finite common-scale response radius')
    base,p = _base(cache); a25 = base['A25']
    if a25['status'] != READY_A25:
        return {'schema':'PHASE2_PURE_SPATIAL_ADAPTER_INPUT_V1',
                'status':'TYPED_NA_NO_NUMERIC_CALL','reason':'FROZEN_A25_TECHNICAL_METADATA_UNAVAILABLE',
                'accepted_A25_status':a25['status'],'minimal_panel':None,'observations':None,
                'effective':None,'typed_NA_ids':[], 'excluded':deepcopy(a25.get('excluded',[])),
                'selection_excluded':deepcopy(a25.get('excluded',[])), 'backend_expected_excluded':[],
                'logical_labels':{'analysis':'POSTHOC_FIXED_COMMON_MEMBER_SCALE','response_radius':response_radius,
                                  'discovery':'xA50','structure':'ORIGINAL_50','backend_variant':'A50'},
                'consumer_signature':None,'source_lineage':deepcopy(cache['source_lineage']),
                'new_numeric_calls':0,'primary_reweighting_calls':0,'A25_requalification_calls':0,
                'Root_exact_release_and_global_dedup_required':True}
    active = [r['id'] for r in a25['rows']]; active_set = set(active)
    need(len(active)==len(active_set) and active==[g for g in p['ids'] if g in active_set],
         'accepted A25 active IDs in original J order')
    positions={g:i for i,g in enumerate(p['ids'])}
    weights=[p['weights'][positions[g]] for g in active]
    need(list(a25['weights']) == weights, 'accepted corresponding original A25 weights')
    rows=[deepcopy(p['rows'][positions[g]]) for g in active]
    need(list(a25['xA']) == [r['xA'] for r in rows], 'accepted A25 uses frozen xA50')
    narrow_byid={r['id']:r for r in a25['rows']}
    for row in rows:
        g=row['id']; geometry=cache['accepted_nested_geometry'][g]
        need(geometry['same_mask_clip_exact'] is True and geometry['id']==g,
             'accepted nested source/mask geometry identity')
        need(geometry['D50']==row['nonexonic_acgt'] and
             geometry['D25']==narrow_byid[g]['nonexonic_acgt'], 'actual accepted radius denominators')
        D=geometry['D50' if response_radius==50000 else 'D25']
        need(type(D) is int and D>0, 'positive actual selected radius D')
        count=cache['coverage'][str(response_radius)][g]
        obs=cache['observations'][str(response_radius)][g]
        if count is None:
            need(obs.get('status')=='TYPED_NA', 'missing coverage preserves source NA')
        else:
            need(count['D']==D, 'saved count denominator equals actual radius')
            if obs.get('status')!='TYPED_NA':
                need((count['A'],count['U'])==(obs['A'],obs['U']), 'same accepted observations/counts')
        row['nonexonic_acgt']=D  # _design never consumes this coverage denominator
    panel={'ids':active,'rows':rows,'weights':weights}
    source_obs=cache['observations'][str(response_radius)]
    need(set(p['ids'])==set(source_obs), 'complete original radius observation ledger')
    observations={g:source_obs[g] for g in active}
    labels={'analysis':'POSTHOC_FIXED_COMMON_MEMBER_SCALE', 'discovery':'xA50',
            'target':'A'+str(response_radius//1000), 'response_radius':response_radius,
            'backend_variant':'A50', 'J':'ORIGINAL_A25_ACCEPTED_ACTIVE_SUBSET_IN_B0_ORDER',
            'weights':'CORRESPONDING_ORIGINAL_B0_UNRENORMALIZED', 'structure':'ORIGINAL_50',
            'denominator':'ACTUAL_RADIUS_D', 'not_original_A25':True,
            'preprocessing':'ORIGINAL_BACKEND_RECOMPUTES_SAME_50_STRUCTURE_ON_FIXED_COMMON_ROWS_PER_DELETION'}
    return _result(cache,panel,observations,labels,a25['excluded'])
