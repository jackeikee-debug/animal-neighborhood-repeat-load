"""Pure PREPARATION ONLY adapters. No paths opened, SQL, network, QA or masks.
Root exclusively admits immutable input bytes/code, their roles and provenance,
resource/once/exact releases. Returned drafts prove no acceptance or eligibility.
"""
from copy import deepcopy
import hashlib
import json
import re

SUPPORTED = 'SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY'
MODERN = 'ACCEPTED_COMPACT_GFF3_GENOME_QA_JSON'
SQLITE = 'ACCEPTED_READONLY_SQLITE'
LEGACY = 'LEGACY32_ACCEPTED_COMPACT_RECEIVER_AND_NAMESPACE_JSON'
NA_ONLY = 'TYPED_NA_ONLY_NO_FACTS_OR_FAKE_QA'
ASSIGNMENTS = ('primary_callable_assigned', 'real_chromosome_assigned', 'model_eligible_assigned')
COORD_INTS = ('RNA_start1', 'RNA_end1', 'full_original_gene_start1', 'full_original_gene_end1')
RAW_FIELDS = ('canonical_seqid_raw', 'source_gene_id_raw', 'strand_raw', 'actual_token_raw')
CARRY = ('target', 'window_index', 'radius', 'canonical_seqid', 'requested_start0',
         'requested_end0', 'gene_count_including_target', 'full_original_gene_length',
         'structural_status', 'structural_exon_bp')

class AdapterHold(ValueError):
    pass

def need(ok, why):
    if not ok:
        raise AdapterHold(why)

def ref(r):
    need(type(r) is dict and type(r.get('path')) is str and r['path'].startswith('/')
         and re.fullmatch('[0-9a-f]{64}', r.get('sha256', '')) is not None,
         'explicit immutable Root-admitted reference')
    return r

def provenance(p, kind):
    """Declarations are not an authorization mechanism; Root verifies them outside."""
    need(p['schema'] == 'ROOT_FACT_VIEW_PROVENANCE_V1' and p['kind'] == kind,
         'explicit format-specific provenance')
    keys = ['catalogue_ref', 'catalogue_acceptance_ref', 'canonical_acceptance_ref',
            'source_disposition_ref', 'code_freeze_ref']
    if kind != NA_ONLY:
        keys += ['Genome_QA_ref', 'Genome_acceptance_ref', 'scientific_artifact_ref']
    if kind == MODERN:
        keys += ['original_Genome_manifest_ref', 'original_manifest_lineage_ref']
    for k in keys:
        ref(p[k])
    need(p['actual_projection_released_by_this_adapter'] is False,
         'adapter cannot release actual projection')
    return deepcopy(p)

def literal_ascii(v, name, *, empty=False):
    """Receiver JSON has ALREADY undone producer %25/%09 transport escaping.
    Encode the resulting literal ASCII; never URL-decode GFF %xx or strip/version
    normalize/unicode-normalize. SQLite __bytes_hex__ is decoded only by its
    original genome_protocol.unpack; this adapter never converts it again.
    """
    need(type(v) is str and v.isascii() and (empty or bool(v)), name + ' literal ASCII')
    return v.encode('ascii')

def targets_index(targets):
    need(type(targets) is list and len(targets) <= 20000, 'finite ordered target list')
    seen = set()
    for t in targets:
        need(type(t['target_key']) is str and t['target_key'] and t['target_key'] not in seen,
             'exact unique target order')
        need(type(t['coordinate_status']) is str and t['coordinate_status'], 'typed source status')
        seen.add(t['target_key'])
    return [t for t in targets if t['coordinate_status'] == SUPPORTED], [t for t in targets if t['coordinate_status'] != SUPPORTED]

def typed_na(record, original, *, compact):
    key = 'target' if compact else 'target_key'
    status = 'input_status' if compact else 'coordinate_status'
    reason = 'reason' if compact else 'NA_reason'
    need(record[key] == original['target_key'] and record[status] == original['coordinate_status'],
         'exact typed-NA identity/status/order')
    need(record.get(reason) == (original.get('NA_reason') or '') if compact else
         record.get(reason) == original.get('NA_reason'), 'literal typed-NA reason')
    return {'target_key': record[key], 'coordinate_status': record[status],
            'NA_reason': record.get(reason), 'adapter_original_typed_NA': deepcopy(record)}

def window_expectation(t, x, index):
    """Only equality against independent frozen manifest fields. No reconstruction
    of annotation, masks, source sequences, counts or scientific source acceptance.
    structural_exon_bp is the sum of ALREADY accepted manifest exon lengths.
    """
    need(t['coordinate_status'] == SUPPORTED, 'only defined original windows')
    for k in COORD_INTS + ('window_center0',):
        need(type(t[k]) is int, 'manifest exact integer ' + k)
    a, z, rs, re_ = [t[k] for k in ('full_original_gene_start1', 'full_original_gene_end1', 'RNA_start1', 'RNA_end1')]
    need(1 <= a <= rs <= re_ <= z and t['strand'] in ('+', '-'), 'manifest RNA/fullgene containment')
    need(t['window_center0'] == (rs - 1 if t['strand'] == '+' else re_), 'unchanged RNA center')
    radius = x['radius']; need(type(radius) is int and radius in (50000, 25000), 'frozen radius')
    need(x['full_original_gene_length'] == z-a+1, 'independent manifest fullgene length')
    ex = x['exon_intervals']
    need(type(ex) is list and all(type(q) is list and len(q) == 2 and all(type(n) is int for n in q) and q[0] < q[1] for q in ex), 'accepted manifest exon-pair schema')
    return dict(target=t['target_key'], window_index=index, radius=radius,
                canonical_seqid=t['canonical_seqid'], requested_start0=t['window_center0']-radius,
                requested_end0=t['window_center0']+radius,
                gene_count_including_target=x['gene_count_including_target'],
                full_original_gene_length=x['full_original_gene_length'], structural_status=x['structural_status'],
                structural_exon_bp=sum(q[1]-q[0] for q in ex))

def from_compact_audit(a, t, window, index):
    """Independent original manifest target ONLY; never the row binding under test."""
    expected = window_expectation(t, window, index)
    for k in CARRY:
        need(type(a[k]) is type(expected[k]) and a[k] == expected[k], 'full original target/window carry: ' + k)
    f = deepcopy(a)
    for k in COORD_INTS:
        need(k not in a or (type(a[k]) is int and a[k] == t[k]), 'no contradictory audit coordinate: ' + k)
        f[k] = t[k]
    fields = {'source_gene_id_raw': literal_ascii(t['source_gene_ID'], 'original source gene ID'),
              'strand_raw': literal_ascii(t['strand'], 'original strand'),
              'canonical_seqid_raw': literal_ascii(t['canonical_seqid'], 'original seqid'),
              'actual_token_raw': literal_ascii(a['actual_token'], 'actual token', empty=True),
              'actual_header_raw': literal_ascii(a['actual_header'], 'actual header', empty=True)}
    for k, v in fields.items():
        need(k not in a or a[k] == v, 'no contradictory audit byte field: ' + k)
    f.update(fields)
    for k in ASSIGNMENTS:
        need(type(f[k]) is int and f[k] == 0, 'preserve unassigned source fact')
    return f

class DraftView:
    """Driver-compatible .facts()/.typed_na() context; each iterator is single-use.
    Full physical input EOF/hash/stability are caller duties, not inferred by this
    in-memory closure. Early exit or skipped/extra/reordered row stays HOLD.
    No masks() method exists: segment counts cannot stand in for raw intervals.
    """
    def __init__(self, facts, nas, expected_facts, expected_nas, p, check=lambda: None):
        self._facts, self._nas = facts, nas
        self._nf, self._nn = expected_facts, expected_nas
        self.provenance = p; self._check = check
        self._started = set(); self._done = set(); self.counts = {}
    def __enter__(self):
        self._check(); return self
    def __exit__(self, exc_type, *args):
        if exc_type is None:
            need(self._done == {'facts', 'typed_na'}, 'both complete adapter iterators/EOF required')
        self._check()
    def _consume(self, name, rows, expected):
        need(name not in self._started, 'single-use iterator, no replay')
        self._started.add(name); n = 0
        for row in rows:
            self._check(); n += 1
            need(n <= expected, 'extra adapter rows')
            yield deepcopy(row)
        need(n == expected, 'missing adapter rows')
        self.counts[name] = n; self._done.add(name); self._check()
    def facts(self):
        return self._consume('facts', self._facts, self._nf)
    def typed_na(self):
        return self._consume('typed_na', self._nas, self._nn)
    def closure(self):
        need(self._done == {'facts', 'typed_na'}, 'draft EOF not reached')
        return {'status': 'DRAFT_ADAPTER_ITERATORS_COMPLETE_NOT_ROOT_ACCEPTANCE',
                'counts': deepcopy(self.counts), 'mask_intervals_read_or_reconstructed': False,
                'source_EOF_or_file_hash_verified_by_adapter': False,
                'eligibility_J_weights_models_assigned': False}

def compact_view(qa, original_manifest, p, *, check=lambda: None):
    """Root provides pinned accepted QA + independently pinned ORIGINAL manifest.
    No receiver/materializer rerun. Exact manifest semantic hash is the original
    producer convention (json sort_keys/compact/default ensure_ascii=True).
    """
    provenance(p, MODERN)
    m = original_manifest
    need(m['schema'] == 'GENOME_TARGET_WINDOW_INPUT_V1', 'original manifest family')
    semantic = hashlib.sha256(json.dumps(m, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    need(semantic == qa['manifest_sha256'] == p['original_manifest_semantic_sha256'], 'independent original manifest hash binding')
    need(qa['schema'] == 'GENOME_WINDOW_RECEIVER_V1' and qa['status'] == 'SOURCE_GEOMETRY_AND_COUNT_ACCOUNTING_ACCEPTED', 'existing accepted receiver family')
    for k in ('primary_callable_or_model_eligible', 'real_chromosome_assigned', 'repeat_outcomes_computed', 'weights_or_model_computed'):
        need(qa[k] is False, 'no scientific assignment in accepted QA')
    supported, nas = targets_index(m['targets']); expected = 2*len(supported)
    need(qa['terminal']['targets'] == len(m['targets']) and qa['terminal']['windows'] == expected
         and len(qa['window_audits']) == expected and len(qa['typed_NA']) == len(nas), 'all-target/window/typed-NA count closure')
    need(qa['terminal']['input_bytes'] == m['expected_source']['bytes'], 'original source accounting preserved')
    need(qa['output_bytes'] == p['scientific_artifact_ref']['bytes'] and qa['output_sha256'] == p['scientific_artifact_ref']['sha256'], 'accepted raw-derived descriptor binding')
    def facts():
        i = 0
        for t in supported:
            need([w['radius'] for w in t['windows']] == [50000, 25000], 'original full-target two-radius order')
            for window in t['windows']:
                check(); a = qa['window_audits'][i]; i += 1
                yield from_compact_audit(a, t, window, i)
    def nrows():
        for r, t in zip(qa['typed_NA'], nas):
            check(); yield typed_na(r, t, compact=True)
    return DraftView(facts(), nrows(), expected, len(nas), deepcopy(p), check)

def sqlite_view(existing_accepted_view, ordered_targets, p, *, check=lambda: None):
    """No DB open/hash here. Root must first use frozen open_view(path, proof),
    which rehashes the exact closed accepted DB. Pass original unpacked rows;
    no enrichment from the binding being checked and no QA or mask rerun.
    ordered_targets is a Root-pinned ordered identity/status expectation only.
    Caller owns/enters/closes the original DB context, enclosing this DraftView.
    """
    provenance(p, SQLITE); supported, nas = targets_index(ordered_targets)
    def facts():
        it = iter(existing_accepted_view.facts()); i = 0
        for t in supported:
            for radius in (50000, 25000):
                check(); f = next(it, None); i += 1
                need(f is not None and f['target'] == t['target_key'] and type(f['radius']) is int and f['radius'] == radius and type(f['window_index']) is int and f['window_index'] == i, 'exact SQLite fact target/window/order')
                for k in RAW_FIELDS:
                    need(type(f[k]) is bytes, 'original unpack literal bytes: ' + k)
                for k in ASSIGNMENTS:
                    need(type(f[k]) is int and f[k] == 0, 'unassigned SQLite fact')
                yield f
        need(next(it, None) is None, 'SQLite facts physical cursor EOF/no extras')
    def nrows():
        it = iter(existing_accepted_view.typed_na())
        for t in nas:
            check(); r = next(it, None); need(r is not None, 'missing SQLite typed NA')
            yield typed_na(r, t, compact=False)
        need(next(it, None) is None, 'SQLite typed NA physical cursor EOF/no extras')
    return DraftView(facts(), nrows(), 2*len(supported), len(nas), deepcopy(p), check)

def na_only_view(original_targets, p, *, check=lambda: None):
    provenance(p, NA_ONLY); supported, nas = targets_index(original_targets)
    need(not supported and p.get('Genome_QA_ref') is None and p.get('scientific_artifact_ref') is None,
         'honest no-defined-window/no-fake-QA route')
    rows = [typed_na(t, t, compact=False) for t in nas]
    return DraftView(iter(()), iter(rows), 0, len(nas), deepcopy(p), check)

def legacy32_view(*args, **kwargs):
    raise AdapterHold('LEGACY32_NEW_ADAPTER_ONLY_UNRESOLVED_NAMESPACE_SCHEMA: accepted source remains accepted; Root must observe exact pinned namespace_QC schema and bind original receiver/manifest/namespace authorities before a distinct adapter successor')

def na_only_driver_view(original_targets, p, *, check=lambda: None):
    """Exact existing skeleton empty_read_view route. Its groups/bindings emit
    source typed NA; it deliberately never consumes Genome typed_na on this path.
    Validate/preserve all original source NA separately, emit no fake Genome rows.
    Root's final target/NA accounting must compare source_typed_na with bindings.
    """
    source = na_only_view(original_targets, p, check=check)
    with source:
        need(list(source.facts()) == [], 'NA-only has no facts')
        records = list(source.typed_na())
    result = DraftView(iter(()), iter(()), 0, 0, deepcopy(p), check)
    result.source_typed_na = records
    return result
