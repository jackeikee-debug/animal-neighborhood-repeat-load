"""First actual SVG display of accepted compact summaries. No scientific fits."""
from pathlib import Path
import csv,importlib.util,json,sys,time
def load(n,p):
    s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def build(cfg,out,check,events):
    start=time.monotonic();source=load('v12plot_sources',cfg['sources_code']);src=source.Sources(cfg,check,events);out=Path(out);out.mkdir();sys.path.insert(0,cfg['code_dir'])
    from v12_summary_figures_compact_v3 import figure1,figure4
    from v12_case_figure_v2 import render_svg as case_svg,supplement_rows
    treeplot=load('v12tree_plot',cfg['treeplot_code']);read=lambda n:json.loads(''.join(src.lines(n)))
    f3=read('fig3');context=read('cases');clade=read('clade');contrasts={r['contrast']:r for r in csv.DictReader(src.lines('contrasts'),delimiter='\t')};proofs=[read(n)for n in ['effect_proof','R23_proof','score_proof']]
    refs=[{k:x[k]for k in ['path','sha256']}for x in cfg['sources']]
    base=dict(accepted_compact_summaries=True,source_refs=refs,title='Accepted comparative neighborhood summaries')
    d1=base|dict(selection=dict(total=9214,selected=1896,eligible_unselected=5599,excluded=1719,unit='original coordinate-frame orthogroups'),denominators=[dict(n=n,label=label,definition=definition,unit='species')for n,label,definition in [(54,'Original tree roster','All original tree tips, including discovery species and outgroups'),(51,'Added coordinate targets','All added species, retaining original reporting roles'),(45,'Main-role candidates','Normal-method candidates; seven without qualifying members and four source NA'),(38,'Admitted profile routes','Nonempty normal J with accepted response route, including four source NA'),(34,'Finite main effects','Finite original full A50/K50/A25 effects'),(35,'Fixed working cohort','Pre-outcome S35; four source NA remain, P31 is separate')]],source_NA=dict(n=9,unit='anchor-by-species paired-measurement units',definition='Nine source-NA pairs across four species, corresponding to 48 abnormal source rows; distinct units'),scope_note='OG0000725 is an unresolved orthogroup exception:154genes/41species, not absent species. R23 composition source support:32species, separate diagnosis.')
    (out/'Figure-1-coverage.svg').write_text(figure1(d1))
    source.write_json(out/'Figure-1-display-data.json',d1)
    # A less crowded common axis keeps every finite value and zero visible.
    f3['effect_axis']['ticks']=[-.1,0,.2]
    (out/'Figure-3-working-tree.svg').write_text(treeplot.render_svg(f3));source.write_json(out/'Figure-3-display-data.json',f3)
    def cr(key,label):
        r=contrasts[key];return dict(label=label,median_abs_paired_delta=float(r['median_abs_delta_beta']),finite_pairs=int(r['finite_pairs']),sign_changes=int(r['strict_sign_changes']))
    retention=[]
    for layer,display in [('Deuterostomia','Deuterostomia'),('Ecdysozoa','Ecdysozoa'),('Spiralia','Spiralia'),('Xenacoelomorpha_separate','Xenacoelomorpha'),('non-Bilateria','Non-bilaterian')]:
        z={r['rule']:r for r in clade['group_summaries']if r['layer']==layer};row=dict(label=display)
        for target,rule in [('all_common','ALL_NONEMPTY_NORMAL_J'),('representatives','TWO_LARGEST_FROZEN_J')]:
            r=z[rule];row[target]=dict(median_retention=r['median_retention'],n_species=r['species'],finite_species=r['finite_full_beta'],common_group_count=r['common_J'])
        retention.append(row)
    d4=base|dict(paired_2x2=dict(statistic='median_abs_paired_delta',axis=[0,.06],unit='absolute standardized coefficient change',rows=[cr('TARGET_A50_PREDICTOR_K_VS_A','Predictor K vs A; target A'),cr('TARGET_K50_PREDICTOR_A_VS_K','Predictor A vs K; target K'),cr('PREDICTOR_A_TARGET_K_VS_A','Target K vs A; predictor A'),cr('PREDICTOR_K_TARGET_K_VS_A','Target K vs A; predictor K')],interpretation='Original J, weights, structure50 fixed; four purposes retain 34finite/11NA in45main-role species.'),scale_change=dict(statistic='median_abs_paired_delta',axis=[0,.03],unit='absolute standardized coefficient change',rows=[cr('FIXED_COMMON_SCALE_Y25_VS_Y50','Target25 vs50; common J')],interpretation='Fixed common J, corresponding original weights, xA50 and structure50. Original A25 has other preprocessing.'),retention=dict(unit='median species-level retained fraction of original frozen J',rows=retention,interpretation='G common oldgroups; S species; finite models; median is across species, not G divided by one common J.'),dependency=dict(blocks=42,nodes_percent=100*1399/1441,score_effective=3.095305224831465,mean_score_dependence=.08195283178344641,working_range=None,NA_reason='P31mean max squared-score share .5561337885>.20; effective scoreblocks3.095<10. All31species ranges fail their own gates, finite betas retained.'),positive_result_note='P31 conditional mean beta .08195; five layer-omission means .07227-.09581; each positive.',limits_note='Conditional descriptive panel; no calibrated working interval; source-NA S35 remains NA; mapping identical-input aliases are not independent replicates.',scope_note='Mapping own A50 median absolute changes .001651-.005058 across five cases; common-map A50 vs baseline .005681. SeeFigure3 and supplementary exact tables.')
    svg=figure4(d4).replace('Overlap blocks','Global closure blocks').replace('Nodes in scope','Largest block nodes').replace('Score effective','Effective score blocks').replace('Mean score dependence','Conditional mean beta').replace('information units','effective blocks')
    svg=svg.replace('</svg>','<text x="374" y="196" font-family="Arial" font-size="12">Mapping own median |Δβ| across five cases</text><text x="374" y="216" font-family="Arial" font-size="12">0.001651–0.005058; common-map 0.005681</text><text x="374" y="241" font-family="Arial" font-size="12">All contrasts: 34 finite pairs; 11 NA retained.</text></svg>')
    (out/'Figure-4-sensitivity.svg').write_text(svg);source.write_json(out/'Figure-4-display-data.json',d4)
    cases=[]
    for c in context['cases']:
        c=dict(c);c['selection_is_outcome_blind']=True;c['selection_note']='Predeclared five cases; lexicographic first shared originalJ group, no locus response/effect selection';r23={};windows={}
        for w in c['windows']:
            w=dict(w);r=str(w['radius']);co=w.pop('composition');r23[r]=dict(status='AVAILABLE',**co)if co is not None else dict(status='TYPED_NA',reason=c['composition_NA_reason']);w.update(old_group=c['old_group'],target_key=c['target_key']);windows[r]=w
        c['windows']=windows;c['r23']=r23;c['focal_gene']=dict(zip(['start0','end0'],c['focal_gene']))if c['focal_gene']is not None else None
        c['all_J_effects']={k:v|dict(denominator='ORIGINAL_J')for k,v in c['all_J_effects'].items()};cases.append(c)
    (out/'Figure-5-fixed-cases.svg').write_text(case_svg(cases));source.write_json(out/'Figure-5-display-data.json',cases)
    table=supplement_rows(cases)
    with(out/'Figure-5-exact-locus-counts-and-model-context.tsv').open('x',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]),delimiter='\t');w.writeheader();w.writerows(table)
    assert set(src.done)==set(src.sources)
    source.write_json(out/'terminal.json',dict(status='V12_FIRST_ACTUAL_FIGURES_1_3_4_5_SVG_COMPLETE_PENDING_VISUAL_QA',source_receipts=src.done,source_profiles=54,main_finite=34,case_group=context['selection_old_group'],case_common_groups=context['common_members'],new_models_SQL_QC_tree=0,seconds=time.monotonic()-start))
    return dict(numeric_call_counts=dict(new_figures=4,old_science=0),source_receipts=src.done)
