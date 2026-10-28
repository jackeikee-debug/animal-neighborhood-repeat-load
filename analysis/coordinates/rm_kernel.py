"""Pure bounded OLD-LITERAL RepeatMasker table/window union kernel, PREP ONLY.
No IO/network/DB/DNA/QA/model/backend. Root owns sources, original masks,
immutable raw NA references, format/exception authority, once and actual release.
Unsorted source rows supported; compact endpoint arrays bounded before append.
"""
from array import array
from copy import deepcopy
import bisect,hashlib,json,re,time

CLASSES=(b'DNA',b'LINE',b'LTR',b'SINE',b'RC',b'Retroposon',b'PLE',b'SVA',b'Unknown',b'Other',b'Excluded')
EXCLUDED={b'Low_complexity',b'Simple_repeat',b'Satellite',b'RNA',b'tRNA',b'rRNA',b'scRNA',b'snRNA',b'srpRNA'}
HEADER=b'SW score\t% div.\t% del.\t% ins.\tquery sequence\tpos in  query: begin\tend\t(left)\trepeat\tclass/family\tpos in repeat: begin\tend\t(left)\tID'
FORMAT='OLD_ACCEPTED_TABULAR_2LINE_REPEAT_OUT_V1'
MASK_BASIS='ALL_original_exons_OR_nonACGT'
NUM=re.compile(br'[0-9]+(?:[.][0-9]+)?\Z')
INTEGER=re.compile(br'[0-9]+\Z')
REP_POS=re.compile(br'[(]?[0-9,]+[)]?\Z')
REP_LEFT_NEG=re.compile(br'[(]-[0-9,]+[)]\Z')
QUERY_LEFT=re.compile(br'[(][0-9,]+[)]\Z')

class InterfaceHold(ValueError):pass
def poison_on_error(method):
    def guarded(self,*args,**kwargs):
        need(not self.failed,'kernel already HOLD; no automatic continuation/retry')
        try:return method(self,*args,**kwargs)
        except BaseException as error:
            self.failed=True;self.failure=type(error).__name__+':'+str(error);raise
    return guarded
def need(ok,why):
    if not ok:raise InterfaceHold(why)
def ref(r):
    need(type(r) is dict and type(r.get('path')) is str and r['path'] and
         re.fullmatch('[0-9a-f]{64}',r.get('sha256','')) is not None and len(r)<=16 and all(type(k) is str and type(v) in (str,int,bool,type(None)) and (type(v) is not str or len(v)<=2048) for k,v in r.items()),'explicit compact scalar Root pinned reference, no body')
def integer(n,name):need(type(n) is int,name+' exact integer');return n
def event_ref(r):
    need(type(r) is dict and set(r)=={'schema','complete_event_retained','source_ref','physical_record_index','rawline_sha256','event_kind','target_key','radius','row_preparation_status','old_group_from_cache','source_status','protection_status'} and all(type(r[k]) in (str,int,bool,type(None)) for k in r if k!='source_ref') and r.get('schema')=='ROOT_LOSSLESS_TARGET_EVENT_REFERENCE_V1' and r.get('complete_event_retained') is True,'full immutable original target event reference')
    ref(r['source_ref']);need(type(r.get('physical_record_index')) is int and r['physical_record_index']>0 and re.fullmatch('[0-9a-f]{64}',r.get('rawline_sha256','')) is not None,'exact raw event record/hash, not reconstructed NA')

def reasons(xs):
    need(type(xs) is list and len(xs)<=64 and all(type(x) is str and x and len(x)<=1024 for x in xs) and sum(map(len,xs))<=4096,'literal bounded NA reasons')

def category(raw):
    """Exact repeat_statistics.py category(): slash PREFIX, case-sensitive."""
    prefix=raw.split(b'/',1)[0]
    if prefix in CLASSES[:9]:return CLASSES.index(prefix)
    return 10 if prefix in EXCLUDED else 9

def parse_data(line,nr,aliases,known):
    """Literal old TAB14forward/TAB15C(+optional final TAB*) parser.
    Unlocalizable query/field count -> source HOLD. Other malformed values are
    localized bad rows; frozen exact known source_NR/raw literals must match.
    No native whitespace conversion, query alias guessing or new class mapping.
    """
    f=line.split(b'\t');n=len(f)-(f[-1]==b'*')
    need(14<=n<=15,'unlocalizable_repeat_fields')
    need(f[4] in aliases and INTEGER.fullmatch(f[5]) and INTEGER.fullmatch(f[6]),'unlocalizable_repeat_query')
    start,end=int(f[5]),int(f[6]);length=aliases[f[4]][1]
    need(1<=start<=end<=length,'unlocalizable_repeat_query')
    bad=not(all(NUM.fullmatch(x) for x in f[:4]) and QUERY_LEFT.fullmatch(f[7]))
    if f[8]==b'C':
        repeat,rawclass,rl,re_,rb,rid=(f[i] if i<len(f) else b'' for i in (9,10,11,12,13,14));bad=bad or n!=15
    else:repeat,rawclass,rb,re_,rl,rid=f[8:14];bad=bad or n!=14
    bad=bool(bad or not repeat or not rawclass or not INTEGER.fullmatch(rid)
             or not REP_POS.fullmatch(rb) or not REP_POS.fullmatch(re_)
             or not(REP_POS.fullmatch(rl) or REP_LEFT_NEG.fullmatch(rl)))
    if nr in known:need(line==known[nr] and bad,'known_exception_identity_changed')
    return {'canonical':aliases[f[4]][0],'query':aliases[f[4]][2],'start0':start-1,'end0':end,
            'bad':bad,'known':nr in known,'rawclass':rawclass,'category':None if bad else category(rawclass)}

def limits(l):
    keys={'max_targets','max_line_bytes','max_source_bytes','max_source_rows','max_events',
          'max_events_per_target','max_mask_segments','max_raw_classes','max_class_bytes',
          'max_anomalies','max_anomaly_bytes','max_aliases','max_alias_bytes','max_known_literal_bytes','max_target_metadata_bytes','max_seconds'}
    need(set(l)==keys and all(type(v) is int and v>0 for v in l.values()),'explicit positive finite limits')
    need(l['max_targets']<=2000 and l['max_line_bytes']<=2**20 and l['max_events']<=2000000
         and l['max_events_per_target']<=250000 and l['max_mask_segments']<=1000000
         and l['max_raw_classes']<=20000 and l['max_class_bytes']<=4*2**20
         and l['max_anomalies']<=100000 and l['max_anomaly_bytes']<=8*2**20 and l['max_aliases']<=100000 and l['max_alias_bytes']<=8*2**20 and l['max_known_literal_bytes']<=2**20
         and l['max_target_metadata_bytes']<=16384 and l['max_seconds']<=180,
         'finite metadata/event/mask/180s ceilings; Root256MiB tree guard still required')

def _bytes_json_size(x):
    def cv(v):
        if type(v) is bytes:return {'__bytes_hex__':v.hex()}
        if type(v) is dict:return {k:cv(z) for k,z in v.items()}
        if type(v) in (list,tuple):return [cv(z) for z in v]
        return v
    return len(json.dumps(cv(x),sort_keys=True,separators=(',',':'),allow_nan=False).encode())

class Kernel:
    """Root future use: construct metadata/masks -> feed each physical table line
    once -> close_source(independent Root receipt) -> iter_measurements() once.
    Each returned record is DRAFT; only Root complete provider/output acceptance
    permits reuse. No source read/mask reader is built into this class.
    """
    def __init__(self,contract,targets,mask_provider,hooks,l,check):
        limits(l);need(array('I').itemsize==4,'exact packed32-bit event storage');self.l=l;self.check_hook=check;self.started=time.monotonic();self.hooks=hooks;self.failed=False;self.failure=None
        self.contract=deepcopy(contract);self.closed=False;self.used=False
        self.stats={'physical_lines':0,'data_rows':0,'parsed_rows':0,'bad_rows':0,'known_rows':0,'table_bytes':0}
        self.digest=hashlib.sha256();self.no_lf=False;self.event_count=0;self.mask_count=0
        self.rawclasses={};self.classbytes=0;self.anomalies=[];self.anomaly_bytes=0;self.targets=[];self.byseq={};self.known_seen=set()
        self.check();c=self.contract
        need(type(c) is dict,'explicit source policy object')
        need(set(c)=={'schema','format','classification_rule','mask_basis','localized_anomaly_policy','count_mode','expected_data_rows','known_exceptions','aliases','source_identity_ref','format_authority_ref','exception_authority_ref','source_alias_authority_ref','frozen_preoutcome_acceptance_ref','frozen_panel_ref','protected_policy_ref','count_policy_ref'},'exact compact source-policy interface')
        need(c['schema']=='ROOT_RM_SOURCE_AND_MEASUREMENT_CONTRACT_V1' and c['format']==FORMAT,
             'HOLD_SOURCE_FORMAT_INTERFACE: exact frozen TAB table contract required; nativewhitespace unsupported')
        for k in ('source_identity_ref','format_authority_ref','exception_authority_ref','source_alias_authority_ref',
                  'frozen_preoutcome_acceptance_ref','frozen_panel_ref','protected_policy_ref','count_policy_ref'):ref(c[k])
        need(c['classification_rule']=='OLD_LITERAL_PREFIX_11_CLASSES' and c['mask_basis']==MASK_BASIS
             and c['localized_anomaly_policy']=='OUTER_A50_OVERLAP_NA_BOTH_SCALES',
             'literal accepted class/mask/anomaly rules, no redesign')
        need(c['count_mode'] in ('EXACT_EXPECTED_DATA_ROWS','ROOT_APPROVED_FIRST_FULL_EOF_ROWS_ACCOUNTING'),
             'explicit source-specific count authority')
        need((c['count_mode']=='EXACT_EXPECTED_DATA_ROWS' and type(c['expected_data_rows']) is int and c['expected_data_rows']>=0)
             or (c['count_mode']=='ROOT_APPROVED_FIRST_FULL_EOF_ROWS_ACCOUNTING' and c['expected_data_rows'] is None),
             'missing row count cannot silently become zero')
        self.known={};knownbytes=0
        need(type(c['known_exceptions']) is list,'explicit source-specific literal exceptions, including accepted emptylist')
        for k in c['known_exceptions']:
            nr=integer(k['source_NR'],'source_NR');raw=k['verbatim_raw_row']
            need(nr>=3 and nr not in self.known and type(raw) is bytes and b'\n' not in raw,'literal unique exception row')
            knownbytes+=len(raw);need(knownbytes<=l['max_known_literal_bytes'],'knownliteral cap');self.known[nr]=raw
        self.aliases={};canonical_lengths={};aliasbytes=0
        for a in c['aliases']:
            q,canon,n=a['query_token_raw'],a['canonical_seqid_raw'],a['length'];ref(a['owner_ref'])
            need(type(q) is bytes and q and b'\t' not in q and b'\n' not in q and q not in self.aliases
                 and type(canon) is bytes and canon and type(n) is int and n>0,'exact unique source-declared alias')
            need(canon not in canonical_lengths or canonical_lengths[canon]==n,'same canonical source length')
            aliasbytes+=len(q)+len(canon);need(len(self.aliases)<l['max_aliases'] and aliasbytes<=l['max_alias_bytes'],'explicit alias count/byte caps');canonical_lengths[canon]=n;self.aliases[q]=(canon,n,q)
        need(self.aliases and type(targets) is list and len(targets)<=l['max_targets'],'finite complete Root targets/alias table')
        need(hooks.verify_inputs(c,targets,l) is True,'Root exact source/target/mask/panel/protected authority before future body')
        seen=set()
        for original in targets:
            self.check();need(_bytes_json_size(original)<=l['max_target_metadata_bytes'],'compact target/ref metadata cap')
            t=deepcopy(original);need(set(t)=={'target_key','old_group','source_status','source_NA_reasons','frozen_J_member','frozen_target_binding_ref','frozen_J_binding_ref','raw_source_event_ref','windows'} | ({'canonical_seqid_raw'} if t['windows'] else set()),'exact compact target interface; full originals retained by event ref');key=t['target_key'];need(type(key) is str and re.fullmatch('t[0-9a-f]{64}',key) is not None and key not in seen,'unique exact already frozen target key');seen.add(key)
            need(type(t['frozen_J_member']) is bool and type(t['source_status']) is str and type(t['source_NA_reasons']) is list,'literal frozen target/status/J binding')
            ref(t['frozen_target_binding_ref']);ref(t['frozen_J_binding_ref']);event_ref(t['raw_source_event_ref']);reasons(t['source_NA_reasons']);need(type(t['old_group']) is str and t['old_group'] and t['raw_source_event_ref']['target_key']==key and t['raw_source_event_ref']['old_group_from_cache']==t['old_group'],'exact original target/group/event provenance')
            need(type(t['windows']) is list and len(t['windows']) in (0,2),'sourceNA or both fixed windows preserved')
            st={'binding':t,'events':array('I'),'localized_observable':False,'anomalies':0,'hits':{50000:0,25000:0},'active':{}}
            if not t['windows']:
                need(t['source_NA_reasons'] and t['raw_source_event_ref']['radius'] is None,'undefined source windows need explicit original typedNA reasons/ref')
                self.targets.append(st);continue
            need([w['radius'] for w in t['windows']]==[50000,25000],'exact frozen two-window order')
            big,small=t['windows'];lo=integer(big['start0'],'start0');hi=integer(big['end0'],'end0')
            need(type(small['start0']) is int and type(small['end0']) is int and hi-lo==100000 and (small['start0'],small['end0'])==(lo+25000,hi-25000),'same frozen center/boundaries')
            canon=t['canonical_seqid_raw'];need(type(canon) is bytes and canon,'literal canonical seqid');need(canon in canonical_lengths or all(w['measurement_status']=='ROOT_TYPED_NA' for w in t['windows']),'active target requires Root source-declared sequence alias');st['localized_observable']=canon in canonical_lengths
            st.update(canonical=canon,start=lo,end=hi)
            mask_intervals={}
            for w in t['windows']:
                need(set(w)=={'radius','start0','end0','measurement_status','NA_reasons','raw_event_ref','mask_acceptance_ref','protected_disposition_ref','nonexonic_ACGT_bp','observed_bp','combined_mask_bp','combined_mask_segments','mask_basis','geometry_status'},'exact compact accepted window metadata');r=w['radius'];event_ref(w['raw_event_ref']);ref(w['mask_acceptance_ref']);ref(w['protected_disposition_ref']);reasons(w['NA_reasons']);need(w['raw_event_ref']['target_key']==key and w['raw_event_ref']['radius']==r and w['raw_event_ref']['old_group_from_cache']==t['old_group'],'exact original window event binding')
                need(type(w['NA_reasons']) is list and w['measurement_status'] in ('ROOT_ADMITTED_MEASURABLE_WINDOW','ROOT_TYPED_NA'), 'Root explicit per-window measurement/NA disposition')
                if w['measurement_status']=='ROOT_TYPED_NA':need(w['NA_reasons'],'typedNA cannot be empty/default');continue
                d=integer(w['nonexonic_ACGT_bp'],'D');obs=integer(w['observed_bp'],'observed_bp')
                need(w['mask_basis']==MASK_BASIS and w['geometry_status']=='GEOMETRICALLY_COMPLETE'
                     and obs==2*r and 0<=w['start0']<w['end0']<=canonical_lengths[canon] and 0<d<=obs,
                     'accepted complete supported nonexonic geometry only; kernel does not select J')
                pairs=[];previous=None;bp=0;mh=hashlib.sha256()
                for pair in mask_provider(key,r):
                    self.check();need(type(pair) in (tuple,list) and len(pair)==2,'accepted mask interval pair')
                    a,z=pair;need(type(a) is int and type(z) is int and w['start0']<=a<z<=w['end0']
                                 and (previous is None or a>previous),'sorted merged accepted ALLexon_OR_nonACGT masks')
                    previous=z;bp+=z-a;self.mask_count+=1;need(self.mask_count<=l['max_mask_segments'],'masksegment cap')
                    pairs.append((a,z));mh.update(f'{a}:{z}\n'.encode())
                need(bp==w['combined_mask_bp']==obs-d and len(pairs)==w['combined_mask_segments'],'accepted original mask/D/segment closure')
                need(hooks.verify_mask_eof(t,w,len(pairs),bp,mh.hexdigest()) is True,'Root original mask cursor/fullEOF/hash authority')
                mask_intervals[r]=pairs;st['active'][r]=w
            if 50000 in mask_intervals and 25000 in mask_intervals:
                clipped=[(max(a,lo+25000),min(z,hi-25000)) for a,z in mask_intervals[50000] if max(a,lo+25000)<min(z,hi-25000)]
                need(clipped==mask_intervals[25000],'accepted A25 mask equals A50 clipped; no mask reconstruction')
            chosen=mask_intervals.get(50000,mask_intervals.get(25000,[]))
            for a,z in chosen:self.add(st,a,z,11)
            self.targets.append(st)
            if st['localized_observable']:self.byseq.setdefault(canon,[]).append(st)
        for canon,ss in self.byseq.items():
            ss.sort(key=lambda x:x['start']);self.byseq[canon]=(ss,[s['end'] for s in ss])
        self.check()

    def check(self):
        need(not self.failed,'kernel HOLD; no continuation')
        need(time.monotonic()-self.started<=self.l['max_seconds'],'one cumulative<=180s kernel deadline')
        self.check_hook()

    def add(self,st,a,z,cat):
        need(st['start']<=a<z<=st['end'],'clipped target endpoint interval')
        need(self.event_count+2<=self.l['max_events'] and len(st['events'])+2<=self.l['max_events_per_target'],'bounded compact endpoint capacity; no automatic partition/replay')
        st['events'].append(((a-st['start'])<<5)|(cat<<1));st['events'].append(((z-st['start'])<<5)|(cat<<1)|1);self.event_count+=2

    @poison_on_error
    def feed_line(self,raw):
        """Root passes each raw physical table line once; optional final noLF is
        allowed exactly as original AWK NR. No subsequent line after noLF.
        Root must separately prove physical EOF/stat/identity without source retry.
        """
        self.check();need(not self.closed and not self.no_lf and type(raw) is bytes and 0<len(raw)<=self.l['max_line_bytes'],'bounded singlepass physical line')
        self.stats['table_bytes']+=len(raw);need(self.stats['table_bytes']<=self.l['max_source_bytes'],'table byte cap');self.digest.update(raw)
        self.stats['physical_lines']+=1;nr=self.stats['physical_lines'];self.no_lf=not raw.endswith(b'\n');line=raw[:-1] if raw.endswith(b'\n') else raw
        if nr==1:need(line==HEADER,'repeat_header/sourceformat HOLD');return
        if nr==2:need(line==b'','repeat_second_line');return
        self.stats['data_rows']+=1;need(self.stats['data_rows']<=self.l['max_source_rows'],'source data rowcap')
        row=parse_data(line,nr,self.aliases,self.known)
        if row['known']:self.known_seen.add(nr);self.stats['known_rows']+=1
        if row['bad']:
            self.stats['bad_rows']+=1;need(len(self.anomalies)<self.l['max_anomalies'],'localized anomaly ledger cap')
            self.anomaly_bytes+=160+len(row['query']);need(self.anomaly_bytes<=self.l['max_anomaly_bytes'],'compact anomaly ledger bytecap');self.anomalies.append({'source_NR':nr,'query_token_raw':row['query'],'start0':row['start0'],'end0':row['end0'],'known':row['known'],'rawline_sha256':hashlib.sha256(raw).hexdigest()})
        else:
            self.stats['parsed_rows']+=1;c=row['rawclass']
            if c not in self.rawclasses:self.classbytes+=len(c);need(len(self.rawclasses)<self.l['max_raw_classes'] and self.classbytes<=self.l['max_class_bytes'],'literal rawclass ledger cap');self.rawclasses[c]=0
            self.rawclasses[c]+=1
        ss,ends=self.byseq.get(row['canonical'],([],[]));i=bisect.bisect_right(ends,row['start0'])
        while i<len(ss) and ss[i]['start']<row['end0']:
            self.check();st=ss[i];i+=1
            if row['bad']:st['anomalies']+=1;continue # outerA50 overlap invalidates BOTH scales
            if not st['active']:continue
            for r,w in st['active'].items():
                if row['start0']<w['end0'] and row['end0']>w['start0']:st['hits'][r]+=1
            a=max(row['start0'],st['start']);z=min(row['end0'],st['end']);self.add(st,a,z,row['category'])

    @poison_on_error
    def close_source(self,receipt):
        self.check();need(not self.closed and self.stats['physical_lines']==self.stats['data_rows']+2,'complete two-header/data physical accounting')
        need(self.stats['data_rows']==self.stats['parsed_rows']+self.stats['bad_rows'] and self.known_seen==set(self.known),'all source rows/known literal exception EOF closure')
        if self.contract['expected_data_rows'] is not None:need(self.stats['data_rows']==self.contract['expected_data_rows'],'exact independently supplied source row count')
        summary={'schema':'RM_KERNEL_FULL_SOURCE_ACCOUNTING_DRAFT_V1',**self.stats,'table_sha256':self.digest.hexdigest(),'rawclass_counts':[{'rawclass':k,'rows':v} for k,v in sorted(self.rawclasses.items())],'anomalies':deepcopy(self.anomalies),'endpoint_records':self.event_count,'mask_segments_read':self.mask_count}
        need(type(receipt) is dict and receipt.get('schema')=='ROOT_REPEAT_SOURCE_FULL_EOF_ONCE_RECEIPT_V1' and receipt.get('source_identity_ref')==self.contract['source_identity_ref'] and receipt.get('complete_physical_EOF') is True and receipt.get('source_stat_closed') is True and receipt.get('table_bytes')==self.stats['table_bytes'] and receipt.get('table_sha256')==self.digest.hexdigest(),'exact independent Root source identity/physicalEOF/stat/tableSHA receipt')
        need(self.hooks.verify_source_eof(self.contract,summary,receipt) is True,'Root independent physical EOF/bytes/SHA/stat/once proof before measurements; no missing=>zero')
        need(_bytes_json_size(receipt)<=16384,'compact Root EOF receipt, not a disguised source body');self.source_receipt=deepcopy(receipt);self.source_summary=summary;self.closed=True;return deepcopy(summary)

    def _unions(self,st):
        self.check();events=sorted(st['events']);opened=[0]*12;prev=0;totals={r:[0]*17 for r in st['active']};i=0
        while i<len(events):
            if i%1024==0:self.check()
            point=events[i]>>5
            if point>prev and opened[11]==0:
                classified=any(opened[c]>0 for c in range(8));un=opened[8]>0;types=sum(opened[c]>0 for c in range(9))
                for r,w in st['active'].items():
                    a=max(prev+st['start'],w['start0']);z=min(point+st['start'],w['end0'])
                    if a>=z:continue
                    bp=z-a;v=totals[r];v[0]+=bp*(classified or un);v[1]+=bp*classified;v[2]+=bp*un;v[3]+=bp*(types>1);v[4]+=bp*(opened[9]>0);v[5]+=bp*(opened[10]>0)
                    for c in range(11):v[6+c]+=bp*(opened[c]>0)
            while i<len(events) and events[i]>>5==point:
                code=events[i];cat=(code>>1)&15;opened[cat]+=(-1 if code&1 else 1);i+=1
            need(all(n>=0 for n in opened),'nonnegative repeat/mask union');prev=point
        need(not any(opened),'closed repeat/mask union EOF');return totals

    def iter_measurements(self):
        need(not self.failed,'kernel HOLD; no output continuation')
        try:yield from self._measurement_rows()
        except BaseException as error:
            self.failed=True;self.failure=type(error).__name__+':'+str(error);raise

    def _measurement_rows(self):
        need(self.closed and not self.used,'complete source accepted hook before zero/one-use draft outputs');self.used=True
        for st in self.targets:
            self.check();t=st['binding'];totals=self._unions(st) if st['active'] and not st['anomalies'] else {}
            for w in (t['windows'] or [None]):
                r=w['radius'] if w else None;reasons=list(w['NA_reasons'] if w else t['source_NA_reasons'])
                if st['anomalies']:reasons.append('localized_repeat_source_anomaly_outer_A50_both_scales_NA')
                value=None
                if r in totals:
                    A,C,U,conflict,other,excluded,*classes=totals[r];D=w['nonexonic_ACGT_bp'];I=C+U-A;K=A-U
                    need(0<=C<=A<=D and 0<=U<=A and 0<=I<=min(C,U) and K>=0 and K+U==A and all(0<=x<=D for x in classes),'old A/C/U/I/K exact union identities')
                    need(U==classes[8] and other==classes[9] and excluded==classes[10] and conflict<=A and I<=conflict,'old category/conflict closure; conflict is NOT I')
                    value={'A':A,'C':C,'U':U,'I':I,'K':K,'D':D,'broad_class_conflict_bp':conflict,'other_class_bp':other,'excluded_class_bp':excluded,'class_coverage':[{'class_raw':k,'bp':v} for k,v in zip(CLASSES,classes)]}
                else:need(reasons,'unmeasured/missing inputs cannot become a measured zero')
                yield {'schema':'RM_TARGET_WINDOW_MEASUREMENT_DRAFT_V1','target_key':t['target_key'],'old_group':t['old_group'],'radius':r,'status':'MEASURED_DRAFT' if value is not None else 'TYPED_NA','coverage':value,'NA_reasons':reasons,'valid_repeat_hits_intersecting_window':st['hits'].get(r,0) if r in st['active'] else None,'localized_anomaly_outer_overlaps':st['anomalies'] if st['localized_observable'] else None,'frozen_J_member':t['frozen_J_member'],'frozen_target_binding_ref':deepcopy(t['frozen_target_binding_ref']),'frozen_J_binding_ref':deepcopy(t['frozen_J_binding_ref']),'raw_source_event_ref':deepcopy(t['raw_source_event_ref']),'window_raw_event_ref':deepcopy(w['raw_event_ref']) if w else None,'protected_disposition_ref':deepcopy(w['protected_disposition_ref']) if w else None,'source_receipt':deepcopy(self.source_receipt),'J_weights_S_or_models_changed':False,'Root_independent_acceptance_required':True}
            st['events']=array('I') # free completed target endpoint storage
        self.check()
