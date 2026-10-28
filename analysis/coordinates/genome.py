"""Pure materializer for one existing Genome FASTA infile. No source/API IO."""
import hashlib
import json
import re

SUP='SOURCE_COORDINATE_SUPPORTED_NOT_CALLABILITY'
STRUCT='STRUCTURAL_SOURCE_ONLY_COMPLETENESS_PENDING'
STRUCT_HOLD='TARGET_GENE_NOT_ONCE_OR_COORDINATE_CONFLICT_HOLD'
CAPS={'max_targets','max_windows','max_input_bytes','max_physical_lines','max_line_bytes',
      'max_output_bytes','max_sequences','max_total_bases','max_retained_header_bytes',
      'max_exon_intervals','max_index_bins','max_index_memberships','max_bin_probes',
      'max_window_line_hits','max_scanned_window_bases','max_exon_intersections',
      'max_nonacgt_runs','max_mask_segments','max_program_bytes'}

def require(ok,why):
    if not ok:raise ValueError(why)

def token(v, empty=False):
    require(isinstance(v,str) and (v or empty) and v.isascii() and all(32<=ord(x)<127 for x in v),'printable ASCII token')
    return v

def plan(m):
    """Shared frozen INPUT schema only; no FASTA counting or mask implementation."""
    encoded=json.dumps(m,sort_keys=True,separators=(',',':')).encode()
    require(len(encoded)<=16*2**20,'manifest cap')
    require(m['schema']=='GENOME_TARGET_WINDOW_INPUT_V1' and m['source_transport']=='FASTA_ASCII_LF_FINAL_LF','exact schema and byte-counted source transport')
    token(m['assembly']);token(m['genome_source_hda'])
    for k in ['canonical_binding_sha256','structural_output_sha256']:
        require(re.fullmatch('[0-9a-f]{64}',m[k]) is not None,'external canonical/structural digest required')
    c=m['caps'];p=m['expected_source']
    require(set(c)==CAPS and all(type(x) is int and 0<x<2**53 for x in c.values()),'exact explicit integer caps')
    require(set(p)=={'bytes','physical_lines','sequences','total_bases'},'expected-source schema')
    require(type(p['bytes']) is int and 0<p['bytes']<=c['max_input_bytes']<2**53,'exact source byte binding')
    for k,cap in [('physical_lines','max_physical_lines'),('sequences','max_sequences'),('total_bases','max_total_bases')]:
        require(p[k] is None or (type(p[k]) is int and 0<p[k]<=c[cap]),'optional exact source counters')
    require(c['max_output_bytes']>=16384 and c['max_line_bytes']<=2**20,'bounded output margin and source line cap')
    require(c['max_sequences']<=100000 and c['max_retained_header_bytes']<=16*2**20,'shared producer/receiver source-receipt ceilings')
    ts=m['targets'];require(isinstance(ts,list) and 0<len(ts)<=min(c['max_targets'],20000),'target cap')
    windows=[];nas=[];keys=set();needed=set();exon_total=0
    for t in ts:
        key=token(t['target_key']);require(key not in keys and t['assembly']==m['assembly'],'target uniqueness/exact assembly');keys.add(key)
        token(t['coordinate_status'])
        if t['coordinate_status']!=SUP:
            token(t.get('NA_reason') or '',True);nas.append(t);continue
        s=token(t['canonical_seqid']);require(re.fullmatch(r'(?:[A-Za-z0-9.:^*$@!+_?|\-]|%[0-9A-Fa-f]{2})+',s) is not None,'valid canonical GFF3 seqid');needed.add(s)
        a=t['full_original_gene_start1'];b=t['full_original_gene_end1'];z=t['window_center0']
        require(all(type(x) is int for x in [a,b,z]) and 1<=a<=b<2**53-100001 and z>=0,'full gene coordinates')
        require(t['strand'] in ['+','-'] and z==(a-1 if t['strand']=='+' else b),'fixed canonical strand center')
        require([x['radius'] for x in t['windows']]==[50000,25000],'exact two frozen scales/order')
        outer=None
        for x in t['windows']:
            require(set(x)=={'radius','gene_count_including_target','full_original_gene_length','structural_status','exon_intervals'},'structural carry schema')
            require(type(x['gene_count_including_target']) is int and x['gene_count_including_target']>=0 and x['full_original_gene_length']==b-a+1,'accepted gene count/fullgene length')
            require(x['structural_status'] in [STRUCT,STRUCT_HOLD],'accepted structural status')
            r=x['radius'];lo=z-r;hi=z+r;previous=-1;bp=0
            require(isinstance(x['exon_intervals'],list),'accepted ALL-original-exon intervals')
            for pair in x['exon_intervals']:
                require(isinstance(pair,list) and len(pair)==2 and all(type(v) is int for v in pair),'interval pair')
                u,v=pair;require(0<=u<v and lo<=u and v<=hi and u>previous,'clipped sorted merged exon intervals');previous=v;bp+=v-u;exon_total+=1
            if outer is not None:
                clipped=[[max(u,lo),min(v,hi)] for u,v in outer if max(u,lo)<min(v,hi)]
                require(clipped==x['exon_intervals'],'A25 ALL-exon mask equals A50 clipped')
            else:outer=x['exon_intervals']
            windows.append(dict(x,target=key,window_index=len(windows)+1,canonical_seqid=s,
                                requested_start0=lo,requested_end0=hi,structural_exon_bp=bp))
    require(len(windows)<=c['max_windows'] and exon_total<=c['max_exon_intervals'],'window/exon caps')
    require(set(m['sequence_aliases'])==needed,'aliases only for exact needed canonical sequences')
    owner={}
    for canonical,aliases in m['sequence_aliases'].items():
        require(isinstance(aliases,list) and aliases and len(set(aliases))==len(aliases),'explicit unique source-declared aliases')
        for alias in aliases:
            token(alias);require(not any(ch.isspace() for ch in alias) and not alias.startswith('>') and alias not in owner,'global exact alias uniqueness');owner[alias]=canonical
    bins=set();memberships=0
    for w in windows:
        lo=w['requested_start0']//65536;hi=(w['requested_end0']-1)//65536
        memberships+=hi-lo+1
        for z in range(lo,hi+1):bins.add((w['canonical_seqid'],z))
        require(memberships<=c['max_index_memberships'] and len(bins)<=c['max_index_bins'],'fixed window index caps')
    return p,c,windows,nas,owner,hashlib.sha256(encoded).hexdigest(),memberships

def materialize(m):
    p,c,ws,nas,owner,mhash,_=plan(m)
    q=lambda s:json.dumps(s,ensure_ascii=True)
    init=['BEGIN { FS=OFS="\\t"; BIN=65536; OFMT=CONVFMT="%.17g"',
          'if (!("version" in PROCINFO)) die("GNU awk PROCINFO/RT capability required")',
          'if (length("\\303\\251")!=2) die("verified byte locale required")']
    init += [k+'='+str(v) for k,v in c.items()]
    init += ['expected_'+k+'='+str(v) for k,v in p.items() if v is not None]
    init += ['manifest_hash='+q(mhash),'NT='+str(len(m['targets'])),'NW='+str(len(ws)),
             'emit("##gff-version 3")','emit("#GENOME_WINDOW_SCHEMA\\tversion=1\\tsemantics=ALL_original_exons_OR_nonACGT\\tmanifest_sha256=" manifest_hash)']
    for t in nas:
        init.append('emit("#TARGET_NA\\ttarget=" '+q(t['target_key'])+' "\\tinput_status=" '+q(t['coordinate_status'])+' "\\treason=" '+q(t.get('NA_reason') or '')+' "\\tprimary_callable_assigned=0")')
    for alias,canonical in owner.items():init.append('OWNER['+q(alias)+']='+q(canonical))
    for w in ws:
        i=str(w['window_index'])
        for k,v in [('TK',w['target']),('WC',w['canonical_seqid']),('SSTATUS',w['structural_status'])]:init.append(k+'['+i+']='+q(v))
        for k,v in [('WR',w['radius']),('WS',w['requested_start0']),('WE',w['requested_end0']),('GCOUNT',w['gene_count_including_target']),('GLEN',w['full_original_gene_length']),('EXBP',w['structural_exon_bp'])]:init.append(k+'['+i+']='+str(v))
        for a,b in w['exon_intervals']:init.append('XS['+i+',++NX['+i+']]='+str(a)+';XE['+i+',NX['+i+']]='+str(b))
        init.append('CW[WC['+i+'],++NCW[WC['+i+']]]='+i+';addindex('+i+')')
    code='\n'.join(init)+'\n}\n'+AWK
    require(len(code.encode())<=c['max_program_bytes'],'materialized program cap')
    return code

AWK=r'''
function die(s){BAD=1;print "#GENOME_ERROR", "physical_line=" NR,s;exit 2}
function emit(s){if(outbytes+length(s)+1>max_output_bytes-2048)die("output margin cap");print s;outbytes+=length(s)+1;outlines++}
function esc(s){gsub(/%/,"%25",s);gsub(/\t/,"%09",s);return s}
function addindex(w,   lo,hi,z,k){
 lo=int(WS[w]/BIN);if(WS[w]<0 && WS[w]!=lo*BIN)lo--;hi=int((WE[w]-1)/BIN)
 for(z=lo;z<=hi;z++){
  k=WC[w] SUBSEP z;if(!(k in BN)){AB[WC[w],++NAB[WC[w]]]=z;activebins++;if(activebins>max_index_bins)die("active bin cap")}
  BN[k]++;BW[k SUBSEP BN[k]]=w;indexkeys++;if(indexkeys>max_index_memberships)die("bin membership cap")
 }
}
function flushmask(w,   s){
 if(!MOPEN[w])return
 MSEG[w]++;masksegments++;if(masksegments>max_mask_segments)die("merged mask segment cap")
 MTOT[w]+=ME[w]-MA[w]
 s=WC[w] "\tderived_exon_or_nonACGT\tregion\t" sprintf("%.0f",MA[w]+1) "\t" sprintf("%.0f",ME[w]) "\t.\t.\t.\tID=genome_mask_w" w "_" MSEG[w] ";codex_window=" w ";codex_mask_basis=ALL_original_exons_OR_nonACGT"
 emit(s);MOPEN[w]=0
}
function maskadd(w,a,b){
 if(a>=b)die("empty union candidate")
 if(MOPEN[w] && a<=ME[w]){if(b>ME[w])ME[w]=b}
 else{flushmask(w);MOPEN[w]=1;MA[w]=a;ME[w]=b}
}
function countwindow(w,a,b,text,    l,r,z,u,ac,gc,k,j,q,xa,xb,nrun,off,run,exa,exb,hasex){
 if(LAST[w]==NR)return;LAST[w]=NR
 l=(a>WS[w]?a:WS[w]);r=(b<WE[w]?b:WE[w]);if(l>=r)return
 windowhits++;scanned+=r-l;if(windowhits>max_window_line_hits || scanned>max_scanned_window_bases)die("target-window work cap")
 z=substr(text,l-a+1,r-l);u=z;ac=gsub(/[ACGT]/,"",u);u=z;gc=gsub(/[GC]/,"",u)
 OBS[w]+=r-l;AC[w]+=ac;GC[w]+=gc
 k=XP[w]+1;while(k<=NX[w] && XE[w,k]<=l)k++;XP[w]=k-1
 for(q=k;q<=NX[w] && XS[w,q]<r;q++){
  xa=(XS[w,q]>l?XS[w,q]:l);xb=(XE[w,q]<r?XE[w,q]:r)
  exonhits++;if(exonhits>max_exon_intersections)die("exon intersection work cap")
  EOBS[w]+=xb-xa;u=substr(text,xa-a+1,xb-xa);EAC[w]+=gsub(/[ACGT]/,"",u)
  u=substr(text,xa-a+1,xb-xa);EGC[w]+=gsub(/[GC]/,"",u)
 }
 for(j in RA)delete RA[j];for(j in RB)delete RB[j]
 nrun=0;off=l;run=z
 while(match(run,/[^ACGT]+/)){
  nonruns++;if(nonruns>max_nonacgt_runs)die("nonACGT run work cap")
  RA[++nrun]=off+RSTART-1;RB[nrun]=RA[nrun]+RLENGTH
  off+=RSTART-1+RLENGTH;run=substr(run,RSTART+RLENGTH)
 }
 q=k;j=1
 while((q<=NX[w] && XS[w,q]<r) || j<=nrun){
  hasex=(q<=NX[w] && XS[w,q]<r);exa=(hasex?(XS[w,q]>l?XS[w,q]:l):0)
  if(hasex && (j>nrun || exa<=RA[j])){exb=(XE[w,q]<r?XE[w,q]:r);maskadd(w,exa,exb);q++}
  else{maskadd(w,RA[j],RB[j]);j++}
 }
}
function visit(chr,z,a,b,text,    key,j,w){
 probes++;if(probes>max_bin_probes)die("bin probe work cap")
 key=chr SUBSEP z
 for(j=1;j<=BN[key];j++){w=BW[key SUBSEP j];countwindow(w,a,b,text)}
}
function finishseq(    j,w){
 if(seq=="")return
 if(bp==0)die("empty FASTA sequence")
 if(canon!=""){
  CLEN[canon]=bp
  for(j=1;j<=NCW[canon];j++){w=CW[canon,j];flushmask(w)}
 }
 emit("#FASTA_SEQUENCE_AUDIT\tindex=" nseq "\ttoken=" esc(seq) "\tactual_header=" esc(header) "\tlength=" bp "\tsequence_lines=" seq_lines "\tmatched=" (canon!=""?1:0) "\tcanonical_seqid=" esc(canon))
}
{
 if(RT!="\n")die("complete LF record separator required; truncated final LF")
 inbytes+=length($0)+length(RT)
 if(inbytes>max_input_bytes || NR>max_physical_lines || length($0)>max_line_bytes)die("source byte/physical/line cap")
 if(index($0,"\r") || $0~/[^\t -~]/)die("ASCII LF source transport")
}
/^>/{
 finishseq();header=substr($0,2);seq=header;sub(/[ \t].*$/,"",seq)
 if(seq=="" || SEEN[seq]++)die("empty or duplicate FASTA header token")
 nseq++;headerbytes+=length(seq);if(nseq>max_sequences)die("source sequence cap")
 canon=(seq in OWNER?OWNER[seq]:"")
 if(canon!=""){
  if(MATCHED[canon]++)die("multiple declared source aliases present; no first selection")
  CTOKEN[canon]=seq;CHEADER[canon]=header;matchedheaderbytes+=length(header)
 }
 if(headerbytes+matchedheaderbytes>max_retained_header_bytes)die("retained header memory cap")
 bp=0;seq_lines=0;next
}
{
 if($0==""){blanks++;next}
 if(seq=="" || $0!~/^[A-Za-z]+$/)die("FASTA letter-only body format")
 text=toupper($0);a=bp;bp+=length(text);total+=length(text);lines++;seq_lines++
 if(total>max_total_bases)die("source base cap")
 if(canon=="")next
 lo=int(a/BIN);hi=int((bp-1)/BIN)
 if(hi-lo+1<=NAB[canon])for(z=lo;z<=hi;z++)visit(canon,z,a,bp,text)
 else for(jj=1;jj<=NAB[canon];jj++){z=AB[canon,jj];if(z>=lo && z<=hi)visit(canon,z,a,bp,text)}
}
END{
 if(BAD)exit 2
 finishseq()
 if(!nseq || inbytes!=expected_bytes || NR!=nseq+lines+blanks)die("complete FASTA byte/physical accounting")
 if(expected_physical_lines && NR!=expected_physical_lines)die("expected physical line mismatch")
 if(expected_sequences && nseq!=expected_sequences)die("expected sequence count mismatch")
 if(expected_total_bases && total!=expected_total_bases)die("expected full-source length mismatch")
 for(w=1;w<=NW;w++){
  chr=WC[w];present=(chr in CLEN);complete=(present && WS[w]>=0 && WE[w]<=CLEN[chr] && OBS[w]==WE[w]-WS[w])
  status=(!present?"SEQID_UNOBSERVED":(complete?"GEOMETRICALLY_COMPLETE":"BOUNDARY_INCOMPLETE"))
  d=(AC[w]+0)-(EAC[w]+0);ng=(GC[w]+0)-(EGC[w]+0)
  if(d<0 || ng<0 || ng>d || (OBS[w]+0)-(MTOT[w]+0)!=d)die("mask/nonexonic denominator closure")
  fraction=(d?sprintf("%.17g",ng/d):"NA");af=(OBS[w]?sprintf("%.17g",AC[w]/OBS[w]):"NA")
  emit("#GENOME_WINDOW_AUDIT\ttarget=" TK[w] "\twindow_index=" w "\tradius=" WR[w] "\tcanonical_seqid=" esc(chr) "\trequested_start0=" sprintf("%.0f",WS[w]) "\trequested_end0=" sprintf("%.0f",WE[w]) "\tactual_token=" esc(CTOKEN[chr]) "\tactual_header=" esc(CHEADER[chr]) "\tsequence_length=" (present?sprintf("%.0f",CLEN[chr]):"NA") "\tobserved_bp=" (OBS[w]+0) "\tacgt_bp=" (AC[w]+0) "\tgc_bp=" (GC[w]+0) "\tstructural_exon_bp=" EXBP[w] "\tobserved_exon_bp=" (EOBS[w]+0) "\texon_acgt_bp=" (EAC[w]+0) "\texon_gc_bp=" (EGC[w]+0) "\tnonexonic_acgt_bp=" d "\tnonexonic_gc_bp=" ng "\tnonexonic_gc_fraction=" fraction "\tobserved_ACGT_fraction=" af "\tcombined_mask_bp=" (MTOT[w]+0) "\tcombined_mask_segments=" (MSEG[w]+0) "\tgeometric_complete=" complete "\tgeometry_status=" status "\tD_status=" (d?"POSITIVE":"D_ZERO_NA") "\tgene_count_including_target=" GCOUNT[w] "\tfull_original_gene_length=" GLEN[w] "\tstructural_status=" SSTATUS[w] "\tprimary_callable_assigned=0\treal_chromosome_assigned=0\tmodel_eligible_assigned=0")
 }
 emit("#GENOME_EOF_OK\tversion=1\tinput_bytes=" inbytes "\tphysical_lines=" NR "\tsequences=" nseq "\tsequence_lines=" (lines+0) "\tblank_lines=" (blanks+0) "\ttotal_bases=" (total+0) "\ttargets=" NT "\twindows=" NW "\tmask_segments=" (masksegments+0) "\twindow_line_hits=" (windowhits+0) "\tscanned_window_bases=" (scanned+0) "\texon_intersections=" (exonhits+0) "\tnonacgt_runs=" (nonruns+0) "\tbin_probes=" (probes+0) "\tbin_index_keys=" (indexkeys+0) "\theader_id_bytes=" (headerbytes+0) "\tmatched_header_bytes=" (matchedheaderbytes+0) "\tpayload_bytes_before_terminal=" outbytes "\toutput_lines_before_terminal=" outlines "\tprimary_callable_assigned=0\treal_chromosome_assigned=0\tmodel_eligible_assigned=0")
}
'''
