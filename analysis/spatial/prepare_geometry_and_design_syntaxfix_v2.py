"""One new source-free geometry/mask/design projection; no real shift responses."""
from pathlib import Path
import csv,gzip,hashlib,json,re,sys
import numpy as np
from scipy.stats import rankdata
P=Path(__import__('os').environ['ANRL_PROJECT_ROOT'])
O=P/'results/migration/2026-10-04-mac-handoff/new-mac-root-v1/three-species-neighborhood-input-statistics-20261013-v1'
ACCS={'G':'GCF_964340395.1','S':'GCF_978046115.1'}
PREFIX={'G':'Glandiceps_talaboti','S':'Saccoglossus_kowalevskii'}
csv.field_size_limit(2**20)
def rows(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='')as f:yield from csv.DictReader(f,delimiter='\t')
def merge(z):
    out=[]
    for a,b in sorted(z):
        assert 0<=a<b
        if out and a<=out[-1][1]:out[-1][1]=max(b,out[-1][1])
        else:out.append([a,b])
    return out
cohort=P/'results/p2/2026-09-24/gene-centered-extraction/cds-transcript-cohort.tsv'
prior={x['orthogroup']for x in rows(cohort)};assert len(prior)==207
forbidden={'G':{},'S':{}};labels=[];seen=set();matched=set();frame_hash=hashlib.sha256()
frame=P/'results/p2/2026-09-24/gene-centered-frame/analysis/anchor-frame-9214.tsv'
def hashed_lines():
    with frame.open('rb')as f:
        for raw in f:frame_hash.update(raw);yield raw.decode()
for x in csv.DictReader(hashed_lines(),delimiter='\t'):
    og=x['orthogroup'];assert og not in seen;seen.add(og)
    if og not in prior:continue
    matched.add(og);ww=json.loads(x['window_fields_by_species'])
    for sp in ACCS:
        y=ww[PREFIX[sp]];a=int(y['window_start0']);b=int(y['window_end0']);ch=y['seqid'];assert b-a==100000
        assert x[PREFIX[sp]+'_accession']==ACCS[sp]
        forbidden[sp].setdefault(ch,[]).append([a,b]);labels.append(dict(kind='prior207',species=sp,orthogroup=og,seqid=ch,start0=a,end0=b))
assert len(seen)==9214 and matched==prior
bed=P/'results/p2/2026-09-24/gene-centered-extraction/window-bed-mapping.tsv'
for x in rows(bed):
    sp=x['species_code']
    if sp not in ACCS or x['split']!='validation':continue
    a,b=int(x['start0']),int(x['end0']);assert b-a==100000
    forbidden[sp].setdefault(x['seqid'],[]).append([a,b]);labels.append(dict(kind='reserved_validation',species=sp,orthogroup=x['orthogroup'],seqid=x['seqid'],start0=a,end0=b))
neighbors=P/'results/p2/2026-09-24/gene-centered-frame/recruitment-ledger/astra-two-reserved-neighbor-overlap.tsv'
for x in rows(neighbors):
    sp=x['species']
    if sp not in ACCS:continue
    a,b=int(x['candidate_start0']),int(x['candidate_end0']);assert b-a==100000
    forbidden[sp].setdefault(x['seqid'],[]).append([a,b]);labels.append(dict(kind='reserved_neighbor',species=sp,orthogroup=x['candidate_OG'],seqid=x['seqid'],start0=a,end0=b))
exception=P/'results/p2/2026-09-24/gene-centered-repeat-exposure/analysis-composite/exception-ledger.tsv'
known={sp:{}for sp in ACCS};all_ex=list(rows(exception));assert len(all_ex)==13
for x in all_ex:
    sp=x['species']
    if sp not in ACCS:continue
    a,b=int(x['query_start0']),int(x['query_end0']);known[sp][x['source_NR']]=x['verbatim_raw_row']
    forbidden[sp].setdefault(x['query_seqid'],[]).append([a,b]);labels.append(dict(kind='unresolved_query',species=sp,seqid=x['query_seqid'],start0=a,end0=b,source_NR=x['source_NR']))
for sp in forbidden:
    for ch in forbidden[sp]:forbidden[sp][ch]=merge(forbidden[sp][ch])
saved=list(rows(O/'Stage-B-main-comparison-v1/main_50kb-saved-values.tsv.gz'))
ids=[x['orthogroup']for x in saved];assert len(ids)==1896 and len(set(ids))==1896 and ids==sorted(ids)
ix={og:i for i,og in enumerate(ids)}
table={(x['species'],x['orthogroup']):x for x in rows(O/'all-1896-triplets-two-scales-input-statistics.tsv.gz')if x['half_window_bp']=='50000'}
species={};arrays={};COVS=['nonexonic_GC_fraction','overlapping_gene_count','target_gene_length_bp','exon_fraction','nonACGT_fraction']
for sp in ACCS:
    j=json.loads((P/f'results/p2/2026-09-24/inputs/ncbi_{ACCS[sp]}_sequence_reports.page-001.json').read_text());assert len(j['reports'])==int(j['total_count'])
    lengths={x['refseq_accession']:int(x['length'])for x in j['reports']}
    code=gzip.decompress((O/f'{sp}-repeat.awk.gz').read_bytes()).decode()
    ws={int(w):(int(a),int(b),ch)for w,a,b,ch in re.findall(r'WS\[(\d+)\]=(\d+); WE\[\1\]=(\d+); WC\[\1\]="([^"]+)";',code)}
    assert len(ws)==1896
    events={}
    for w,a,b in re.findall(r'ev\((\d+),(\d+),"MASK",1\);ev\(\1,(\d+),"MASK",-1\);',code):events.setdefault(int(w),[]).append([int(a),int(b)])
    fixed=sorted(rows(O/f'{sp}-fixed-windows.tsv.gz'),key=lambda x:(x['seqid'],int(x['start0'])))
    fixed=[x for x in fixed if x['half_window_bp']=='50000'];assert len(fixed)==1896
    targets=[None]*1896
    for w,x in enumerate(fixed,1):
        a,b,ch=ws[w];assert(a,b,ch)==(int(x['start0']),int(x['end0']),x['seqid'])
        og=x['orthogroup'];i=ix[og];t=table[sp,og];assert(t['seqid'],int(t['window_start0']),int(t['window_end0']))==(ch,a,b)
        masks=merge(events.get(w,[]));parts=[];pos=a
        for l,r in masks:
            assert a<=l<r<=b
            if pos<l:parts.append([pos,l])
            pos=r
        if pos<b:parts.append([pos,b])
        d=int(t['nonexonic_ACGT_bp']);assert sum(r-l for l,r in parts)==d
        assert all(b<=l or r<=a for l,r in forbidden[sp].get(ch,[])),('zero_forbidden',sp,og,ch,a,b)
        targets[i]=dict(orthogroup=og,index=i,seqid=ch,start0=a,end0=b,callable=parts,D=d,A=int(t['annotated_interspersed_bp']),observed_rank=float(saved[i][sp+'_response_average_rank']))
    aligned=[table[sp,og]for og in ids];chrom=[x['seqid']for x in aligned];levels=sorted(set(chrom))
    assert len(levels)==(23 if sp=='G' else 22)
    values=dict(nonexonic_GC_fraction=[int(x['nonexonic_GC_bp'])/int(x['nonexonic_ACGT_bp'])for x in aligned],overlapping_gene_count=[int(x['overlapping_gene_count'])for x in aligned],target_gene_length_bp=[int(x['target_gene_length_bp'])for x in aligned],exon_fraction=[int(x['all_isoform_exon_union_bp'])/int(x['window_bp'])for x in aligned],nonACGT_fraction=[int(x['nonACGT_bp'])/int(x['window_bp'])for x in aligned])
    fit=json.loads((O/f'Stage-B-main-comparison-v1/main_50kb-{sp}-fixed-adjustment.json').read_text());columns=[np.ones(1896)];names=['intercept']
    for name in COVS:
        z=rankdata(values[name],method='average');sd=float(np.std(z))
        if sd==0:assert name in fit['dropped_constant_covariates'];continue
        assert abs(float(z.mean())-fit['covariate_transforms'][name]['rank_mean'])<1e-10 and abs(sd-fit['covariate_transforms'][name]['rank_sd_ddof0'])<1e-10
        columns.append((z-z.mean())/sd);names.append('rank_z:'+name)
    for ch in levels[1:]:columns.append(np.array([c==ch for c in chrom],float));names.append('chromosome:'+ch)
    X=np.column_stack(columns);assert names==list(fit['coefficients']) and levels==fit['chromosome_levels']
    arrays[sp+'_X']=X;arrays[sp+'_old_residual']=np.array([float(x[sp+'_residual'])for x in saved]);arrays[sp+'_old_rank']=np.array([float(x[sp+'_response_average_rank'])for x in saved])
    species[sp]=dict(assembly=ACCS[sp],lengths=lengths,chromosomes=levels,targets=targets,forbidden=forbidden[sp],known_exceptions=known[sp],X_columns=names,expected_source_rows=1367131 if sp=='G' else 1407175)
arrays['P_fixed_residual']=np.array([float(x['P_residual'])for x in saved])
Path('geometry-and-mask-design.json.gz').write_bytes(gzip.compress(json.dumps(dict(ids=ids,species=species,forbidden_labels=labels),separators=(',',':')).encode(),mtime=0))
np.savez_compressed('original-fixed-design-and-saved-residuals.npz',**arrays)
Path('return.json').write_text(json.dumps(dict(status='SOURCE_FREE_GEOMETRY_MASK_DESIGN_PREPARED_NO_REAL_NULL',targets=1896,G_chromosomes=23,S_chromosomes=22,prior207=len(prior),forbidden_label_rows=len(labels),frame_rows=len(seen),frame_SHA256=frame_hash.hexdigest(),source_reads={'Genome':0,'GFF':0,'RM':0},model_fits=0,real_null_draws=0,protected_outcomes_opened=0),indent=2)+'\n')
