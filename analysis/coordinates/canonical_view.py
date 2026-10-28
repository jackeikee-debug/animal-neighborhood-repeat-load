"""Pure literal accepted-canonical view. No IO, source classification or acceptance.

All input is already accepted by Root. Caps are caller-supplied and bounded by
CAPS. A rejected input produces ValueError, never a truncated or partial view.
The downstream frozen binding_builder owns coordinate equality/geometry checks.
"""
import json
import math
from copy import deepcopy

MARKER = 'COPY_ACCEPTED_CANONICAL_STATUS_NO_RECLASSIFICATION'
CAPS = {'max_record_bytes':1048576, 'max_view_bytes':4194304,
        'max_nodes':50000, 'max_depth':32, 'max_container_items':10000,
        'max_string_chars':1048576, 'max_authority_refs':32}
# Whole present subtrees, no flag-name heuristic or leaf explosion. Unknown keys
# inside these subtrees survive literally; unknown top-level fields survive in
# the full unchanged record. Presence is not truth and values are not classified.
SOURCE_FIELDS = ('coordinate_status', 'NA_reason', 'input_statuses',
    'candidate_evidence', 'source_candidates', 'source_coordinates',
    'source_ancestry_issues', 'original_gene_placement_count',
    'five_prime_uncertainty_status', 'original_transcript_identity_ambiguity',
    'same_coordinate_original_RNA_IDs', 'raw_nomination_attributes',
    'raw_RNA_flags', 'raw_gene_attributes', 'CDS_relation')


def _need(ok, why):
    if not ok:
        raise ValueError(why)


def _bounded_json(value, lim, byte_cap):
    """Validate bounded JSON shape before serialization or deep-copy allocation."""
    state = {'nodes':0, 'text':0}
    active = set()
    def walk(v, depth):
        state['nodes'] += 1
        _need(state['nodes'] <= lim['max_nodes'], 'JSON node cap')
        _need(depth <= lim['max_depth'], 'JSON depth cap')
        t = type(v)
        if t is str:
            _need(len(v) <= lim['max_string_chars'], 'JSON string cap')
            state['text'] += len(v)
            _need(state['text'] <= byte_cap, 'JSON aggregate text cap')
        elif v is None or t is bool:
            pass
        elif t is int:
            _need(v.bit_length() <= 4096, 'JSON integer size cap')
        elif t is float:
            _need(math.isfinite(v), 'JSON number must be finite')
        elif t is dict or t is list:
            _need(len(v) <= lim['max_container_items'], 'JSON container cap')
            _need(id(v) not in active, 'JSON cycles forbidden')
            active.add(id(v))
            if t is dict:
                for k, child in v.items():
                    _need(type(k) is str, 'JSON object keys must be strings')
                    walk(k, depth+1)
                    walk(child, depth+1)
            else:
                for child in v:
                    walk(child, depth+1)
            active.remove(id(v))
        else:
            raise ValueError('JSON-only input; no bytes decoding or object coercion')
    walk(value, 0)
    # ASCII JSON is a conservative byte bound. This serialization is only for
    # size admission; it never replaces or normalizes any returned source value.
    raw = json.dumps(value, ensure_ascii=True, allow_nan=False,
                     separators=(',', ':'))
    _need(len(raw) <= byte_cap, 'serialized JSON byte cap')


def _ref(r):
    _need(type(r) is dict, 'explicit source authority reference')
    path = r.get('path')
    sha = r.get('sha256')
    _need(type(path) is str and bool(path) and path.startswith('/'),
          'explicit absolute authority path; no path normalization')
    _need(type(sha) is str and len(sha)==64 and
          all(c in '0123456789abcdef' for c in sha), 'literal authority sha256')


def prepare_canonical_view(record, source_authority_refs, policy, limits):
    """Return the frozen binding_builder canonical_views four-field input.

    policy must explicitly carry Root's MARKER and its own supplied authority
    ref. Empty blocking_states means no ADDITIONAL classification; an accepted
    canonical NA remains NA through its unchanged coordinate_status/NA_reason.
    This function validates shape/bounds, not source truth or authority contents.
    """
    _need(type(limits) is dict and set(limits)==set(CAPS) and
          all(type(limits[k]) is int and 0 < limits[k] <= CAPS[k] for k in CAPS),
          'explicit bounded literal-view limits')
    lim = limits
    _bounded_json(record, lim, lim['max_record_bytes'])
    _bounded_json(source_authority_refs, lim, lim['max_record_bytes'])
    _bounded_json(policy, lim, lim['max_record_bytes'])
    _need(type(record) is dict, 'canonical record object')
    required = {'schema','key','native_species_name','coordinate_status',
                'NA_reason','source_coordinates','primary_callable_or_model_eligible'}
    _need(required <= set(record), 'missing mandatory canonical field')
    _need(record['schema']=='canonical-source-coordinate-v1', 'canonical schema')
    _need(type(record['key']) is dict and
          all(type(record['key'].get(k)) is str and bool(record['key'][k])
              for k in ('assembly','protein_id')), 'exact assembly/protein identity')
    _need(type(record['native_species_name']) is str and bool(record['native_species_name']),
          'exact nonempty native species identity')
    _need(type(record['coordinate_status']) is str and bool(record['coordinate_status']),
          'accepted coordinate disposition must be present')
    _need(record['NA_reason'] is None or type(record['NA_reason']) is str,
          'accepted NA reason must remain literal string or null')
    _need(record['source_coordinates'] is None or type(record['source_coordinates']) is dict,
          'source coordinate container shape; no coordinate derivation')
    _need(record['primary_callable_or_model_eligible'] is False,
          'source-only canonical record; no invented eligibility')
    _need(type(policy) is dict and
          policy.get('schema')=='ROOT_LITERAL_CANONICAL_VIEW_POLICY_V1' and
          policy.get('approved_by')=='current-conversation-root' and
          policy.get('marker')==MARKER and 'ref' in policy,
          'explicit Root copy-accepted-status/no-reclassification policy')
    _ref(policy['ref'])
    _need(type(source_authority_refs) is list and bool(source_authority_refs),
          'nonempty actual Root source authority refs; no fabricated authority')
    for r in source_authority_refs:
        _ref(r)
    refs = deepcopy(source_authority_refs)
    if policy['ref'] not in refs:
        refs.append(deepcopy(policy['ref']))
    _need(len(refs) <= lim['max_authority_refs'], 'authority reference count cap')
    view = {'record':deepcopy(record),
            'source_flags':[{'literal_path':'/'+name, 'value':deepcopy(record[name])}
                            for name in SOURCE_FIELDS if name in record],
            'blocking_states':[], 'source_authority_refs':refs}
    _bounded_json(view, lim, lim['max_view_bytes'])
    return view
