"""Independent structural-mask v1 receiver: pure binary-stream validation only.

No file opens, network, subprocesses, source extraction, or eligibility decisions.
Caller supplies the independently pinned producer manifest and immutable stream.
"""
import hashlib
import json
import re
import time

SUPPORTED = 'SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY'
SOURCE_ONLY = 'STRUCTURAL_SOURCE_ONLY_COMPLETENESS_PENDING'
HOLD = 'TARGET_GENE_NOT_ONCE_OR_COORDINATE_CONFLICT_HOLD'
CAPS = {'max_targets','max_input_bytes','max_physical_lines','max_output_bytes',
        'max_index_bins','max_index_memberships','max_feature_window_pairs',
        'max_event_keys','max_gene_window_keys','max_retained_id_bytes','max_bin_probes'}
PROJECTION = {'bytes','feature_records','source_lines','sequence_region_directives',
              'distinct_feature_types','physical_lines'}
AUDIT = {'target','window_index','radius','seqid','requested_start0','requested_end0',
         'status','gene_count_including_target','full_original_gene_length',
         'exon_union_bp','union_segments','target_gene_records','inner_half_gene_count',
         'complete_window','primary_callable_assigned'}
TERMINAL = {'version','input_bytes','physical_lines','feature_records',
            'sequence_region_directives','comment_lines','blank_lines','targets','windows',
            'bin_index_keys','bin_probes','feature_window_pairs','event_keys',
            'gene_window_keys','primary_callable_assigned'}

def require(ok, why):
    if not ok:
        raise ValueError(why)

def number(s, signed=False):
    require(isinstance(s,str) and re.fullmatch(r'-?(0|[1-9][0-9]*)' if signed else r'0|[1-9][0-9]*',s) is not None,'integer syntax')
    require(len(s)<=20,'integer magnitude')
    return int(s)

def fields(line, prefix, exact=None):
    parts=line.split('\t')
    require(parts[0]==prefix,'record tag')
    result={}
    for item in parts[1:]:
        k,sep,v=item.partition('=')
        require(sep and k and k not in result,'malformed or duplicate field')
        result[k]=v
    if exact is not None:
        require(set(result)==exact,'exact record fields')
    return result

def manifest_plan(m):
    """Validate receiver-relevant frozen semantics, never a source body."""
    encoded=json.dumps(m,sort_keys=True,separators=(',',':')).encode()
    require(len(encoded)<=8*2**20,'receiver manifest cap')
    require(m.get('source_kind')=='ACCEPTED_FULL_ORIGINAL_STRUCTURAL_PROJECTION_V2','full original projection binding')
    p=m['expected_projection'];c=m['caps'];targets=m['targets']
    require(set(p)==PROJECTION and set(c)==CAPS,'manifest count/cap schema')
    require(all(type(v) is int and v>0 for v in c.values()),'positive integer caps')
    require(all(type(v) is int and v>=0 for v in p.values()),'projection counters')
    require(p['feature_records']>0 and p['bytes']>0,'nonempty accepted projection')
    require(p['physical_lines']==p['feature_records']+p['sequence_region_directives']+p['distinct_feature_types']+6,'F+Q+D+6')
    require(c['max_input_bytes']>=p['bytes'] and c['max_physical_lines']>=p['physical_lines'],'source cap binding')
    require(isinstance(targets,list) and 0<len(targets)<=min(c['max_targets'],20000),'receiver target cap')
    windows=[];nas=[];keys=set()
    for t in targets:
        for k in ['target_key','coordinate_status']:
            require(isinstance(t.get(k),str) and t[k] and t[k].isascii() and all(32<=ord(x)<127 for x in t[k]),'ASCII target token')
        key=t['target_key'];require(key not in keys,'unique target key');keys.add(key)
        if t['coordinate_status']!=SUPPORTED:
            reason=t.get('NA_reason') or ''
            require(isinstance(reason,str) and reason.isascii() and all(32<=ord(x)<127 for x in reason),'NA reason token')
            nas.append(t);continue
        a=t['full_original_gene_start1'];b=t['full_original_gene_end1'];center=t['window_center0']
        rna_start=t['RNA_start1'];rna_end=t['RNA_end1']
        require(all(type(v) is int for v in [a,b,rna_start,rna_end,center]) and 1<=a<=rna_start<=rna_end<=b and center>=0,'supported canonical RNA/full-gene containment')
        require(t['strand'] in ['+','-'] and center==(rna_start-1 if t['strand']=='+' else rna_end),'fixed canonical RNA strand-specific center')
        for k in ['seqid','source_gene_ID']:
            require(isinstance(t[k],str) and t[k] and t[k].isascii() and all(32<=ord(x)<127 for x in t[k]),'ASCII source token')
        require(re.fullmatch(r'(?:[A-Za-z0-9.:^*$@!+_?|\-]|%[0-9A-Fa-f]{2})+',t['seqid']) is not None,'GFF3 seqid requires reserved/space characters escaped')
        for radius in [50000,25000]:
            windows.append({'target':key,'window_index':len(windows)+1,'radius':radius,
                            'seqid':t['seqid'],'requested_start0':center-radius,
                            'requested_end0':center+radius,'full_original_gene_length':b-a+1})
    bins=set();memberships=0
    for w in windows:
        lo=w['requested_start0']//65536;hi=(w['requested_end0']-1)//65536
        memberships+=hi-lo+1
        for b in range(lo,hi+1):bins.add((w['seqid'],b))
        require(len(bins)<=c['max_index_bins'] and memberships<=c['max_index_memberships'],'frozen bin/index caps')
    return p,c,windows,nas,hashlib.sha256(encoded).hexdigest(),memberships

def validate(stream, manifest, *, expected_bytes=None, expected_sha256=None,
             max_bytes=None, max_line_bytes=2**20, max_seconds=180):
    """Consume binary readline(size) through physical EOF; return source-only QA.

    <=20,000 targets, <=8MiB manifest, <=1MiB line, <=180s sampled per line.
    Segment intervals are checked online and never retained. No gzip detection:
    pass a caller-owned decompressor stream and bind its transport independently.
    """
    started=time.monotonic()
    p,c,windows,nas,mhash,index_memberships=manifest_plan(manifest)
    require(type(max_line_bytes) is int and 0<max_line_bytes<=2**20,'line cap')
    require(isinstance(max_seconds,(int,float)) and 0<max_seconds<=180,'time cap')
    budget=c['max_output_bytes'] if max_bytes is None else max_bytes
    require(type(budget) is int and 0<budget<=c['max_output_bytes'],'output cap binding')
    require(expected_bytes is None or (type(expected_bytes) is int and 0<=expected_bytes<=budget),'expected output bytes')
    require(expected_sha256 is None or re.fullmatch('[0-9a-f]{64}',expected_sha256) is not None,'expected output hash')
    max_lines=5+len(nas)+len(windows)+c['max_event_keys']//2
    nbytes=0;nlines=0;digest=hashlib.sha256();state='header';ni=0;wi=0
    segment_count=0;segment_bp=0;previous_end=None;audits=[];upstream={};terminal=None
    clipped50=hashlib.sha256();actual25=hashlib.sha256();native_empty_zero=[]
    def next_state():
        return 'na' if ni<len(nas) else ('windows' if wi<len(windows) else 'source')
    while True:
        require(time.monotonic()-started<=max_seconds,'receiver time cap')
        raw=stream.readline(min(max_line_bytes+1,budget-nbytes+1))
        require(isinstance(raw,bytes),'binary stream required')
        if not raw:break
        nbytes+=len(raw);nlines+=1
        require(nbytes<=budget and len(raw)<=max_line_bytes and nlines<=max_lines,'output byte/line/count cap')
        require(raw.endswith(b'\n') and b'\r' not in raw and b'\x00' not in raw,'LF-only complete line')
        require(all(x==9 or x==10 or 32<=x<127 for x in raw),'ASCII producer transport')
        digest.update(raw);line=raw[:-1].decode('ascii')
        require(state!='done','content after terminal EOF')
        if state=='header':
            require(line=='##gff-version 3','single initial GFF3 header');state='schema';continue
        if state=='schema':
            f=fields(line,'#STRUCTURAL_MASK_SCHEMA',{'version','semantics','target_manifest_sha256'})
            require(f=={'version':'1','semantics':'all_original_gene_pseudogene_and_ALL_exon_union','target_manifest_sha256':mhash},'schema/manifest hash')
            state=next_state();continue
        if state=='na':
            f=fields(line,'#TARGET_NA',{'target','input_status','reason','primary_callable_assigned'})
            t=nas[ni]
            require(f=={'target':t['target_key'],'input_status':t['coordinate_status'],
                        'reason':t.get('NA_reason') or '','primary_callable_assigned':'0'},'exact unsupported typed NA')
            ni+=1;state=next_state();continue
        if state=='windows':
            w=windows[wi]
            if not line.startswith('#'):
                f=line.split('\t');require(len(f)==9,'nine GFF3 columns')
                require(f[0]==w['seqid'] and f[1:3]==['derived_all_original_exons','region'] and f[5:8]==['.','.','.'],'neutral region mask syntax')
                a=number(f[3])-1;b=number(f[4])
                require(0<=a<b and w['requested_start0']<=a and b<=w['requested_end0'],'positive GFF3 clipped interval')
                require(previous_end is None or a>previous_end,'sorted, disjoint, non-touching merged union')
                segment_count+=1;segment_bp+=b-a;previous_end=b
                require(f[8]==f"ID=mask_w{wi+1}_{segment_count};codex_window={wi+1};codex_mask_basis=ALL_original_exons",'exact unique mask/window attributes')
                if wi%2==0:
                    clip_a=max(a,w['requested_start0']+25000);clip_b=min(b,w['requested_end0']-25000)
                    if clip_a<clip_b:clipped50.update(f'{clip_a}:{clip_b}\n'.encode())
                else:actual25.update(f'{a}:{b}\n'.encode())
                continue
            f=fields(line,'#TARGET_WINDOW_AUDIT',AUDIT)
            audit={k:(number(v,signed=k=='requested_start0') if k not in {'target','seqid','status','complete_window'} else v) for k,v in f.items()}
            for k,v in w.items():require(audit[k]==v,'frozen target/window/center/full-gene binding: '+k)
            require(audit['status'] in [SOURCE_ONLY,HOLD] and audit['complete_window']=='NOT_ESTABLISHED' and audit['primary_callable_assigned']==0,'source-only status/eligibility')
            require(audit['union_segments']==segment_count and audit['exon_union_bp']==segment_bp,'segment/union conservation')
            require(segment_bp<=2*w['radius'],'window union bound')
            require(audit['gene_count_including_target']<=p['feature_records'] and audit['target_gene_records']<=p['feature_records'],'source gene count bound')
            require(audit['inner_half_gene_count']<=audit['gene_count_including_target'],'inner count subset')
            if audit['status']==SOURCE_ONLY:
                require(audit['target_gene_records']==1 and audit['gene_count_including_target']>=1,'source-only target once/including-target')
            if audit['target_gene_records']!=1:require(audit['status']==HOLD,'target multiplicity must remain HOLD')
            if wi%2==1:
                a50=audits[-1]
                require(a50['inner_half_gene_count']==audit['gene_count_including_target'],'INNERGENES50 = GENES25')
                require(audit['inner_half_gene_count']<=a50['inner_half_gene_count'],'inner25 subset')
                require(audit['target_gene_records']==a50['target_gene_records'] and audit['status']==a50['status'],'same source target audit both radii')
                require(audit['exon_union_bp']<=a50['exon_union_bp'],'nested window union bound')
                require(clipped50.digest()==actual25.digest(),'A25 union must equal A50 union clipped at the same center')
                clipped50=hashlib.sha256();actual25=hashlib.sha256()
            audits.append(audit);wi+=1;segment_count=0;segment_bp=0;previous_end=None;state=next_state();continue
        if state=='source':
            needed={'source_lines','feature_records','emitted_feature_records','comment_lines','blank_lines',
                    'sequence_region_directives','fasta_directives','fasta_tail_lines'}
            f=fields(line,'#UPSTREAM_SOURCE_COUNTS',needed)
            nums={k:number(f[k]) for k in needed}
            require(nums['source_lines']==p['source_lines'] and nums['feature_records']==nums['emitted_feature_records']==p['feature_records'],'upstream source identity/count')
            require(nums['source_lines']==sum(nums[k] for k in ['feature_records','comment_lines','blank_lines','fasta_tail_lines']),'upstream full-original physical accounting')
            require(nums['sequence_region_directives']==p['sequence_region_directives'],'original/projection sequence-region conservation')
            require(nums['fasta_directives'] in [0,1] and nums['comment_lines']>=nums['sequence_region_directives']+nums['fasta_directives'],'original directive/comment accounting')
            require(nums['fasta_tail_lines']==0 or nums['fasta_directives']==1,'original FASTA tail requires directive')
            upstream['source_counts']=f;state='upstream_eof';continue
        if state=='upstream_eof':
            needed={'version','source_lines','feature_records','emitted_feature_records'}
            f=fields(line,'#UPSTREAM_EOF_OK',needed)
            require(f['version']=='2','upstream EOF version/fields')
            require(number(f['source_lines'])==p['source_lines'] and number(f['feature_records'])==number(f['emitted_feature_records'])==p['feature_records'],'upstream EOF binding')
            for k in set(f)&set(upstream['source_counts']):require(f[k]==upstream['source_counts'][k],'upstream footer agreement: '+k)
            upstream['EOF']=f;state='terminal';continue
        require(state=='terminal','internal parser state')
        f=fields(line,'#STRUCTURAL_MASK_EOF_OK',TERMINAL)
        # These exact AWK counters are uninitialized empty strings until first use.
        # Preserve their provenance; an empty field is numeric zero only here.
        empty_allowed={'windows','bin_index_keys','bin_probes','feature_window_pairs','event_keys','gene_window_keys'}
        native_empty_zero=[k for k,v in f.items() if v=='' and k in empty_allowed]
        terminal={k:(0 if k in native_empty_zero else number(v)) for k,v in f.items()}
        for outkey,pk in [('input_bytes','bytes'),('physical_lines','physical_lines'),('feature_records','feature_records'),('sequence_region_directives','sequence_region_directives')]:
            require(terminal[outkey]==p[pk],'terminal source binding: '+outkey)
        require(terminal['version']==1 and terminal['targets']==len(manifest['targets']) and terminal['windows']==len(windows) and terminal['primary_callable_assigned']==0,'terminal target/window/schema/eligibility')
        require(terminal['comment_lines']==p['sequence_region_directives']+p['distinct_feature_types']+6 and terminal['blank_lines']==0,'accepted projection exact comment/blank accounting')
        require(terminal['physical_lines']==terminal['feature_records']+terminal['comment_lines']+terminal['blank_lines'],'terminal physical closure')
        for key,cap in [('bin_index_keys','max_index_memberships'),('bin_probes','max_bin_probes'),('feature_window_pairs','max_feature_window_pairs'),('event_keys','max_event_keys'),('gene_window_keys','max_gene_window_keys')]:
            require(terminal[key]<=c[cap],'terminal work cap: '+key)
        genes=sum(a['gene_count_including_target'] for a in audits);segments=sum(a['union_segments'] for a in audits)
        require(terminal['gene_window_keys']==genes,'gene/window dedup conservation')
        require(2*segments<=terminal['event_keys']<=2*terminal['feature_window_pairs'],'event/segment bound')
        require((segments==0)==(terminal['event_keys']==0),'zero union/event conservation')
        require(terminal['feature_window_pairs']>=genes+segments,'feature/window work bound')
        require(terminal['feature_window_pairs']<=p['feature_records']*len(windows),'feature/window source upper bound')
        require(terminal['bin_index_keys']==index_memberships,'exact deterministic window/bin memberships')
        require(terminal['feature_window_pairs']==0 or terminal['bin_probes']>0,'feature hits require bin probes')
        state='done'
    require(state=='done' and terminal is not None,'missing/truncated unique terminal EOF')
    require(expected_bytes is None or nbytes==expected_bytes,'output exact byte mismatch')
    actual=digest.hexdigest();require(expected_sha256 is None or actual==expected_sha256,'output SHA mismatch')
    return {'schema':'INDEPENDENT_STRUCTURAL_MASK_RECEIVER_V1','status':'SOURCE_OUTPUT_FORMAT_AND_ACCOUNTING_ACCEPTED',
            'output_bytes':nbytes,'output_lines':nlines,'output_sha256':actual,'target_manifest_sha256':mhash,
            'targets':len(manifest['targets']),'windows':len(windows),'typed_NA':len(nas),
            'source_HOLD_windows':sum(a['status']==HOLD for a in audits),'window_audits':audits,
            'typed_NA_targets':[{'target':t['target_key'],'input_status':t['coordinate_status'],'reason':t.get('NA_reason') or ''} for t in nas],
            'terminal':terminal,'terminal_native_empty_zero_fields':native_empty_zero,'upstream':upstream,'window_completeness_assigned':False,
            'callability_assigned':False,'primary_callable_or_model_eligible':False,
            'source_or_outcome_computed':False}
