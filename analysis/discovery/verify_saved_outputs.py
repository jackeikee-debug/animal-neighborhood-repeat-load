"""Check saved results and reporting; no refit, new random draws or source scan."""
import hashlib, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from pypdf import PdfReader

O=Path(__file__).resolve().parent
root=Path.cwd()
scope=json.loads((O/'execution-scope.json').read_text())
S=json.loads((O/'comparison-summary.json').read_text())
for pin in scope['inputs']:
    p=root/pin['path']
    assert p.stat().st_size==pin['bytes']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==pin['sha256']
assert hashlib.sha256((O/'run_comparison.py').read_bytes()).hexdigest()==scope['code_sha256']
assert S['actual_counts']=={'fits':6,'shifts':1998}
checks=[]; rows=[]
for s in S['versions']:
    v=s['version']; n=s['n_final']
    panel=json.loads((O/(v+'-panel.json')).read_text())
    a=pd.read_csv(O/(v+'-saved-values.tsv.gz'),sep='\t')
    q=pd.read_csv(O/(v+'-null.tsv.gz'),sep='\t')
    shifts=pd.read_csv(O/(v+'-cyclic-shifts.tsv.gz'),sep='\t')
    influence=pd.read_csv(O/(v+'-fixed-residual-chromosome-influence.tsv'),sep='\t')
    assert len(a)==n and a.orthogroup.is_unique
    assert a.orthogroup.tolist()==panel['final_groups']
    assert len(q)==999 and q.draw.tolist()==list(range(1,1000))
    assert len(shifts)==999*45 and shifts.groupby('draw').size().eq(45).all()
    assert ((shifts.offset>=0)&(shifts.offset<shifts.n)).all()
    rs=[]
    for x,y in [('G','P'),('G','S'),('P','S')]:
        key=x+y
        r=float(np.corrcoef(a[x+'_residual'],a[y+'_residual'])[0,1])
        raw=float(np.corrcoef(a[x+'_response_average_rank'],a[y+'_response_average_rank'])[0,1])
        assert abs(r-s['residual_Pearson'][key])<1e-12
        assert abs(raw-s['raw_Spearman'][key])<1e-12
        rs.append(r)
        rows.append(dict(version=v,n=n,pair=key,raw_Spearman=raw,residual_Pearson=r))
    assert abs(np.mean(rs)-s['T'])<1e-12
    assert np.allclose(q[['GP','GS','PS']].mean(axis=1),q['T'],rtol=0,atol=1e-12)
    exceed=int((q['T']>=s['T']).sum())
    assert (1+exceed)/1000==s['empirical_upper_tail_p']
    assert np.allclose(np.quantile(q['T'],[.025,.975]),[s['null_T_q025'],s['null_T_q975']],rtol=0,atol=1e-12)
    assert len(influence)==67 and influence.all_pairs_positive.all()
    assert np.allclose(influence[['GP','GS','PS']].mean(axis=1),influence['T'],rtol=0,atol=1e-12)
    assert abs(influence['T'].min()-s['fixed_residual_chromosome_influence']['T_min'])<1e-12
    assert abs(influence['T'].max()-s['fixed_residual_chromosome_influence']['T_max'])<1e-12
    fit_checks=[]
    for species in 'GPS':
        d=json.loads((O/f'{v}-{species}-fixed-adjustment.json').read_text())
        assert d['rank']==d['columns'] and d['checks']['finite'] and d['checks']['residual_positive_variance']
        assert d['checks']['relative_normal_equation_error']<1e-12
        fit_checks.append({'species':species,'rank':d['rank'],'columns':d['columns'],'condition_number':d['condition_number'],'normal_equation_error':d['checks']['relative_normal_equation_error']})
    checks.append({'version':v,'n':n,'saved_correlations_and_T_match':True,'background_rows':len(q),'shift_log_rows':len(shifts),'null_exceedances':exceed,'influence_rows':len(influence),'saved_fit_checks':fit_checks})
assert S['S2']['additional_fits']==0 and S['S2']['additional_shifts']==0
pdf=O/'three-species-main-and-sensitivity-figures-v2.pdf'
reader=PdfReader(pdf)
assert len(reader.pages)==3
assert all('X Lable' not in p.extract_text() and 'Y Lable' not in p.extract_text() for p in reader.pages)
v6=root/'results/migration/2026-10-04-mac-handoff/new-mac-root-v1/manuscript-v6-RNA-integration-20261013-v1/editorial-receipt.json'
v6_receipt=json.loads(v6.read_text())
for p in v6_receipt['outputs']:
    f=root/p['path']
    assert f.stat().st_size==p['bytes'] and hashlib.sha256(f.read_bytes()).hexdigest()==p['sha256']
out={'time_utc':datetime.now(timezone.utc).isoformat(),'status':'SAVED_OUTPUTS_VERIFIED','reviewer':'current conversation /root; not an independent agent or new scientific review','inputs_and_frozen_code_unchanged':True,'checks':checks,'new_fits':0,'new_random_draws':0,'S2_extra_fits':0,'v6_source_receipt':str(v6.relative_to(root)),'v6_files_hash_matched':len(v6_receipt['outputs']),'figure_pdf':pdf.name,'figure_pages':3,'visual_review':'All three final pages viewed; readable axes, data points and summaries; no observed clipping or unwanted point labels.','presentation_history':'First figure draft retained. ReportLab ScatterPlot defaults added point labels and dummy axis titles. Presentation-only v2 disables those defaults, using identical saved results.'}
with (O/'saved-output-verification.json').open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
with (O/'pairwise-results.tsv').open('x') as f:pd.DataFrame(rows).to_csv(f,sep='\t',index=False)
print(json.dumps({'status':out['status'],'versions':len(checks),'v6_hash_matches':out['v6_files_hash_matched'],'new_fits':0,'new_draws':0}))
