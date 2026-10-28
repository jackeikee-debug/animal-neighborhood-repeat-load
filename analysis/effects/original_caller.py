"""Pure finite frozen-backend caller; no IO, imports of backend, RM or CLI.
Root alone verifies actual source hashes/EOF/roles/runtime/resources, claims and
accepts each saved result. No function grants execution or retry authority.
"""
from copy import deepcopy
from dataclasses import dataclass
from types import MappingProxyType
import json
import math
import re
import struct

BACKEND_SHA='36a901dce504df25c8e2d8cc0f9df2810b41ceb004d4d64b2d53e149699f0c26'
SUPPLEMENT_SHA='1eb89c955bb6420bf33435850d479cfa3989e645b21862b4d6871311e42d4ac1'
VARIANTS=('A50','K50','A25')
MODEL_LABELS=('A50','K50','A25','C_top')
LAYERS=('non-Bilateria','Spiralia','Ecdysozoa','Deuterostomia','Xenacoelomorpha_separate')
READY_A25='PREPARED_FROM_EXACT_ACCEPTED_HELPERS'


def need(ok,msg):
    if not ok:raise ValueError('HOLD: '+msg)


def _ref(r):
    need(type(r) is dict and type(r.get('path')) is str and r['path'].startswith('/')
         and type(r.get('sha256')) is str and re.fullmatch('[0-9a-f]{64}',r['sha256']) is not None,
         'explicit immutable Root source/acceptance ref')


def _finite(v):
    return type(v) in (int,float) and math.isfinite(v)


def _freeze(v):
    if type(v) is dict:
        need(all(type(k) is str for k in v),'literal string mapping keys')
        return MappingProxyType({k:_freeze(x) for k,x in v.items()})
    if type(v) in (tuple,list):return tuple(_freeze(x) for x in v)
    need(v is None or type(v) in (str,bytes,bool,int,float),'literal serialized value')
    if type(v) is float:need(math.isfinite(v),'prepared input must be finite')
    return v


def _thaw(v):
    if isinstance(v,MappingProxyType):return {k:_thaw(x) for k,x in v.items()}
    if type(v) is tuple:return [_thaw(x) for x in v]
    return v


def _authority(a,schema):
    need(type(a) is dict and a.get('schema')==schema and a.get('approved_by')=='current-conversation-root'
         and a.get('complete_source_hash_EOF_verified') is True,'Root admitted complete immutable inputs')


def _model_authority(a):
    _authority(a,'ROOT_ACCEPTED_FINITE_MODEL_INPUT_BINDING_V1')
    need(a['backend_sha256']==BACKEND_SHA and a['numerical_supplement_sha256']==SUPPLEMENT_SHA,
         'exact already accepted backend and numerical supplement')
    _ref(a['preoutcome_acceptance_ref'])


@dataclass(frozen=True)
class PreparedSpecies:
    native_species_name:str
    assembly:str
    fixed_study_role:str
    fixed_major_layer:str
    panel:object
    A25:object
    authority:object
    def materialize_panel(self):return _thaw(self.panel)
    def materialize_A25(self):return _thaw(self.A25)
    @property
    def ids(self):return self.panel['ids']


@dataclass(frozen=True)
class FrozenObservations:
    native_species_name:str
    assembly:str
    radius:int
    records:object
    authority:object
    def materialize(self):return _thaw(self.records)


@dataclass(frozen=True)
class FixedStudy:
    S_order:tuple
    records:object
    covariance:object
    authority:object


def restore_species_panel(record,binding):
    """Immutable constructor from Root-loaded accepted JSON, NO numerical calls.
    All panel/A25 keys and values survive; list representation is restored when
    materialized for the existing backend, which already accepts list weights.
    Full source-event/NA lineage remains in its pinned original panel artifact.
    """
    _model_authority(binding)
    for k in ('panel_ref','panel_acceptance_ref'):_ref(binding[k])
    need(record['native_species_name']==binding['native_species_name'] and record['assembly']==binding['assembly'],
         'exact native and assembly')
    role=record['fixed_study_role']
    need(role in ('new_normal_candidate','outgroup_background','Hemicordylus_method_exception','gene_tree_only'),
         'new51 roles only; old3 reused outside caller')
    p=record['panel'];ids=p['ids'];rows=p['rows'];w=p['weights']
    required={'status','ids','rows','weights','components','components_per_location','information','excluded','local_descriptive_ids'}
    need(required<=set(p) and type(ids) is list and len(ids)<=2000 and len(set(ids))==len(ids)
         and all(type(x) is str and x for x in ids),'complete accepted panel shape/IDs')
    need(binding['expected_main_J_count']==len(ids) and type(rows) is list and [r['id'] for r in rows]==ids,
         'exact frozenJ and row order, no relabeling')
    need(type(w) is list and len(w)==len(ids) and all(_finite(x) and x>0 for x in w)
         and len(p['components'])==len(ids),'literal saved weights/components, never recomputed')
    for row in rows:
        need(all(k in row for k in ('location','real_chromosome','start0','end0','complete','acgt_fraction','nonexonic_acgt','structure','xA','xK')),
             'full frozen backend row fields')
        need(type(row['structure']) is list and len(row['structure'])==5 and
             all(_finite(x) for x in row['structure']) and _finite(row['xA']) and _finite(row['xK']),
             'literal finite saved structure/signatures, no rescaling')
    a=record['A25']
    need(type(a) is dict and type(a.get('status')) is str,'literal accepted A25 state')
    if a['status']==READY_A25:
        need(a['main_J_unchanged'] is True and a['main_weights_unchanged'] is True
             and set(a['narrow_rows'])==set(ids),'A25 exact mainJ metadata')
        active=[r['id'] for r in a['rows']]
        need(len(set(active))==len(active) and set(active)<=set(ids)
             and active==[g for g in ids if g in set(active)],'saved preset A25 subset/order only')
        positions={g:i for i,g in enumerate(ids)}
        need(a['weights']==[w[positions[g]] for g in active],'exact corresponding main weights; no renormalization')
    return PreparedSpecies(record['native_species_name'],record['assembly'],role,record['fixed_major_layer'],
                           _freeze(p),_freeze(a),_freeze(binding))


def restore_observations(records,species,radius,binding):
    _authority(binding,'ROOT_ACCEPTED_FINITE_OBSERVATION_BINDING_V1')
    for k in ('source_ref','acceptance_ref'):_ref(binding[k])
    need(radius in (50000,25000) and type(radius) is int and binding['radius']==radius
         and binding['native_species_name']==species.native_species_name and binding['assembly']==species.assembly
         and binding['panel_sha256']==species.authority['panel_ref']['sha256'],'exact frozen panel/source/radius binding')
    need(type(records) is dict and set(records)==set(species.ids),'observations exact mainJ keys, including preset excluded A25 slots')
    for r in records.values():
        need(type(r) is dict,'literal observation object')
        if r.get('status')=='TYPED_NA':
            need(type(r.get('reason')) is str and r['reason'],'explicit source NA reason; no zero imputation')
        else:
            need(r.get('status') in (None,'MEASURED') and type(r.get('A')) is int and type(r.get('U')) is int
                 and 0<=r['U']<=r['A'],'integer measured A/U including truezero, or explicit typedNA')
    return FrozenObservations(species.native_species_name,species.assembly,radius,_freeze(records),_freeze(binding))


def restore_fixed_study(fixed_s,covariance,binding):
    _model_authority(binding)
    for k in ('fixed_S_ref','covariance_ref','covariance_acceptance_ref'):_ref(binding[k])
    order=fixed_s['S_order'];records=fixed_s['working_record_identity_template'];n=len(order)
    need(type(order) is list and 0<n<=45 and len(set(order))==n and n==binding['expected_S_count']
         and order==binding['expected_S_order']==covariance['S_order'],'exact accepted S order, no reselection')
    need(fixed_s['working_cohort_eligible'] is True and [r['species'] for r in records]==order
         and all(r['role']=='new_normal_metazoan' and r['spatial_eligible'] is True and r['layer'] in LAYERS for r in records),
         'already accepted fixed normal-new S/roles/layers')
    for name in ('C','C_top'):
        matrix=covariance[name]
        need(type(matrix) is list and len(matrix)==n and all(type(r) is list and len(r)==n
             and all(_finite(x) for x in r) for r in matrix),'saved finite same-order covariance; no matrix derivation')
    return FixedStudy(tuple(order),_freeze(records),_freeze(covariance),_freeze(binding))


def _once(claim,operation,inputs):
    r=claim(operation,deepcopy(inputs))
    need(type(r) is dict and r.get('schema')=='ROOT_ONCE_FINITE_CALL_CLAIM_V1'
         and r.get('approved_by')=='current-conversation-root' and r.get('operation_id')==operation
         and r.get('status')=='CLAIMED_ONCE' and r.get('inputs_immutable_bound') is True,
         'Root durable exact-once claim before backend call')
    return deepcopy(r)


class FiniteBackendHold(RuntimeError):
    def __init__(self,operation,receipt,error):
        super().__init__('HOLD: consumed/unknown finite call '+operation+': '+type(error).__name__+': '+str(error))
        self.operation_id=operation;self.claim_receipt=receipt;self.automatic_retry=False


def _invoke(operation,inputs,function,check,claim):
    check();receipt=_once(claim,operation,inputs);check()
    try:result=function()
    except Exception as error:raise FiniteBackendHold(operation,receipt,error) from error
    check();return result,receipt


def iter_species_versions(species,observations50,observations25,backend,check,claim):
    """Yield/persist each version before advancing; no duplicate fit_species call.
    Root must save each yielded result durably, preserve partial/UNKNOWN state on
    interruption and never resume/replay a consumed operation automatically.
    """
    for variant in VARIANTS:
        check();base={'schema':'FINITE_SPECIES_VERSION_RESULT_V1','native_species_name':species.native_species_name,
            'assembly':species.assembly,'variant':variant,'fixed_study_role':species.fixed_study_role,
            'fixed_J_ids':list(species.ids),'panel_ref':_thaw(species.authority['panel_ref']),
            'preoutcome_acceptance_ref':_thaw(species.authority['preoutcome_acceptance_ref']),
            'spatial_jackknife_calls':0,'backend_result':None,'automatic_retry':False}
        if species.fixed_study_role=='gene_tree_only':
            need(observations50 is None and observations25 is None,'gene-only has no fabricated repeat response')
            yield dict(base,status='NA',reason='FIXED_GENE_ONLY_NO_USABLE_REPEAT_RESPONSE');continue
        obs=observations25 if variant=='A25' else observations50
        need(isinstance(obs,FrozenObservations) and obs.native_species_name==species.native_species_name
             and obs.assembly==species.assembly and obs.radius==(25000 if variant=='A25' else 50000)
             and obs.authority['panel_sha256']==species.authority['panel_ref']['sha256'],'same accepted observations/panel/radius')
        base['observation_binding']=_thaw(obs.authority)
        a=species.materialize_A25()
        if variant=='A25' and a['status']!=READY_A25:
            yield dict(base,status='NA',reason='FROZEN_A25_TECHNICAL_METADATA_UNAVAILABLE',A25_metadata=a);continue
        panel=species.materialize_panel();measurements=obs.materialize()
        active=[r['id'] for r in a['rows']] if variant=='A25' else list(species.ids)
        byid={r['id']:r for r in panel['rows']}
        unavailable=[g for g in active if measurements[g].get('status')=='TYPED_NA']
        if unavailable:
            yield dict(base,status='NA',reason='FROZEN_OBSERVATIONS_HAVE_TYPED_NA',unavailable_ids=unavailable,
                       unavailable_records={g:measurements[g] for g in unavailable});continue
        for g in active:
            D=a['narrow_rows'][g]['nonexonic_acgt'] if variant=='A25' else byid[g]['nonexonic_acgt']
            need(0<=measurements[g]['U']<=measurements[g]['A']<=D,'frozen measured outcome within accepted D')
        op='SPATIAL:'+species.authority['panel_ref']['sha256']+':'+variant
        inputs={'panel_binding':_thaw(species.authority),'observation_binding':_thaw(obs.authority),'variant':variant}
        kwargs={'narrow_rows':a['narrow_rows']} if variant=='A25' else {}
        result,receipt=_invoke(op,inputs,lambda:backend.spatial_jackknife(panel,measurements,variant=variant,**kwargs),check,claim)
        # Entire full/deletion/NA/constant result is retained without filtering.
        yield dict(base,status='BACKEND_RETURNED',spatial_jackknife_calls=1,backend_result=deepcopy(result),claim_receipt=receipt)


def call_working_version(study,species_results,model_label,input_binding,backend,check,claim):
    check();need(model_label in MODEL_LABELS,'finite model labels only')
    variant='A50' if model_label=='C_top' else model_label
    _authority(input_binding,'ROOT_ACCEPTED_WORKING_RESULT_INPUT_BINDING_V1')
    need(input_binding['source_variant']==variant and set(input_binding['result_refs'])==set(study.S_order),
         'all exact fixed-S saved result sources, Ctop reuses A50')
    for refs in input_binding['result_refs'].values():
        _ref(refs['source_ref']);_ref(refs['acceptance_ref'])
    records=[];unestimable=[];diagnostics={}
    for template in study.records:
        name=template['species'];r=species_results.get(name,{}).get(variant)
        if r is None:
            unestimable.append(name);diagnostics[name]={'reason':'SAVED_VERSION_RESULT_UNOBSERVED'};continue
        need(r['schema']=='FINITE_SPECIES_VERSION_RESULT_V1' and r['native_species_name']==name
             and r['variant']==variant and r['preoutcome_acceptance_ref']==_thaw(study.authority['preoutcome_acceptance_ref']),
             'exact saved version/native/preoutcome lineage')
        result=r.get('backend_result');full=result.get('full') if isinstance(result,dict) else None
        beta=full.get('beta') if isinstance(full,dict) else None;v=result.get('vJ') if isinstance(result,dict) else None
        if not (isinstance(result,dict) and result.get('status')=='OK' and isinstance(full,dict)
                and full.get('status')=='OK' and _finite(beta) and _finite(v) and v>=0):
            unestimable.append(name);diagnostics[name]=deepcopy(r);continue
        records.append({**_thaw(template),'beta':beta,'vJ':v})
    base={'schema':'FINITE_WORKING_VERSION_RESULT_V1','model_label':model_label,'source_variant':variant,
        'fixed_S_order':list(study.S_order),'preoutcome_acceptance_ref':_thaw(study.authority['preoutcome_acceptance_ref']),
        'input_binding':deepcopy(input_binding),'working_backend_calls':0,'backend_result':None,'automatic_retry':False}
    if unestimable:
        return dict(base,status='NA',reason='FIXED_S_MEMBER_UNESTIMABLE_NO_SPECIES_DROPPED',
                    unestimable_species=unestimable,member_diagnostics=diagnostics)
    matrix=_thaw(study.covariance['C_top' if model_label=='C_top' else 'C'])
    op='WORKING:'+study.authority['preoutcome_acceptance_ref']['sha256']+':'+model_label
    result,receipt=_invoke(op,{'study_binding':_thaw(study.authority),'results_binding':input_binding,'model_label':model_label},
                           lambda:backend.fit_working_summary(records,matrix),check,claim)
    need(not isinstance(result,dict) or result.get('optimizer_calls',0)<=1,'accepted backend finite optimizer call count')
    return dict(base,status='BACKEND_RETURNED',working_backend_calls=1,backend_result=deepcopy(result),claim_receipt=receipt,
        species_records=deepcopy(records),beta_ordered=[r['beta'] for r in records],vJ_ordered=[r['vJ'] for r in records],
        layers_ordered=[r['layer'] for r in records],covariance_ref=_thaw(study.authority['covariance_ref']))


def call_primary_fixed_v_influence(study,saved_A50,input_binding,backend,check,claim):
    check();_authority(input_binding,'ROOT_ACCEPTED_FIXED_V_INPUT_BINDING_V1')
    _ref(input_binding['source_ref']);_ref(input_binding['acceptance_ref'])
    need(saved_A50['schema']=='FINITE_WORKING_VERSION_RESULT_V1' and saved_A50['model_label']=='A50'
         and saved_A50['fixed_S_order']==list(study.S_order)
         and saved_A50['preoutcome_acceptance_ref']==_thaw(study.authority['preoutcome_acceptance_ref']),
         'only independently accepted saved PRIMARY A50 V, exact frozenS')
    r=saved_A50.get('backend_result')
    base={'schema':'FINITE_PRIMARY_FIXED_V_INFLUENCE_V1','model_label':'A50','fixed_S_order':list(study.S_order),
          'saved_result_binding':deepcopy(input_binding),'variance_refits':0,'backend_calls':0,'automatic_retry':False}
    if not isinstance(r,dict) or r.get('status')!='OK' or 'V' not in r:
        return dict(base,status='NA',reason='ACCEPTED_PRIMARY_A50_V_UNAVAILABLE',backend_result=None)
    beta=deepcopy(saved_A50['beta_ordered']);layers=deepcopy(saved_A50['layers_ordered']);V=deepcopy(r['V'])
    need(layers==[t['layer'] for t in study.records] and len(beta)==len(study.S_order),'same saved beta/V/layer order')
    op='FIXED_V_PRIMARY_A50:'+study.authority['preoutcome_acceptance_ref']['sha256']
    result,receipt=_invoke(op,{'study_binding':_thaw(study.authority),'saved_A50_binding':input_binding},
                           lambda:backend.fixed_v_influence(beta,V,layers),check,claim)
    need(type(result) is list and len(result)<=5 and len({r['deleted_layer'] for r in result})==len(result)
         and all(r['deleted_layer'] in LAYERS for r in result),'at most5 predefined fixed-V layer outputs total')
    return dict(base,status='BACKEND_RETURNED',backend_calls=1,backend_result=deepcopy(result),claim_receipt=receipt)


_BYTES='__bytes_hex__'
_FLOAT='__float64_bits_hex__'

def _wire(value):
    active=set();nodes=[0]
    def cv(v,depth):
        nodes[0]+=1;need(nodes[0]<=2000000 and depth<=64,'result serialization bounds')
        t=type(v)
        if v is None or t in (str,bool,int):return v
        if t is float:
            return v if math.isfinite(v) else {_FLOAT:struct.pack('>d',v).hex()}
        if t is bytes:return {_BYTES:v.hex()}
        if t in (dict,list,tuple) or isinstance(v,MappingProxyType):
            need(id(v) not in active,'cyclic result');active.add(id(v))
            try:
                if t is dict or isinstance(v,MappingProxyType):
                    need(all(type(k) is str for k in v) and set(v) not in ({_BYTES},{_FLOAT}),
                         'literal result dict collides with reserved reversible tag')
                    return {k:cv(x,depth+1) for k,x in v.items()}
                return [cv(x,depth+1) for x in v]
            finally:active.remove(id(v))
        if t.__module__.split('.')[0]=='numpy':
            if callable(getattr(v,'tolist',None)):return cv(v.tolist(),depth+1)
            if callable(getattr(v,'item',None)):return cv(v.item(),depth+1)
        raise ValueError('HOLD: unsupported result value, no stringify/coercion')
    return cv(value,0)


def encode_result(result):
    """Lossless diagnostic codec; nonfinite tags never become model inputs."""
    raw=json.dumps(_wire(result),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    need(len(raw)<=64*2**20,'one finite-result packed byte cap; Root accounts whole output')
    return raw


def decode_result(raw):
    need(type(raw) is bytes and len(raw)<=64*2**20,'bounded serialized finite result')
    def pairs(items):
        d={}
        for k,v in items:need(k not in d,'duplicate result JSON key');d[k]=v
        if set(d)=={_BYTES}:
            h=d[_BYTES];need(type(h) is str and len(h)%2==0 and re.fullmatch('[0-9a-f]*',h) is not None,'literal byte hex tag')
            return bytes.fromhex(h)
        if set(d)=={_FLOAT}:
            h=d[_FLOAT];need(type(h) is str and re.fullmatch('[0-9a-f]{16}',h) is not None,'literal IEEE754 diagnostic tag')
            return struct.unpack('>d',bytes.fromhex(h))[0]
        return d
    def bad_constant(x):raise ValueError('HOLD: untagged nonfinite JSON constant')
    return json.loads(raw,object_pairs_hook=pairs,parse_constant=bad_constant)
