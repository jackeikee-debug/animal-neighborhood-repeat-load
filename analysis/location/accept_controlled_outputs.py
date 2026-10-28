"""Complete saved-output QA and logical alias/NA expansion; no projections."""
from pathlib import Path
import csv,gzip,hashlib,json
import numpy as np
Q=Path(__import__('os').environ['ANRL_R4_ROOT']);R=Q.parent;D=Q/'controlled-distances-first'
freeze=json.loads((Q/'geometry-and-saved1998-first/fixed-distance-admission-freeze.json').read_text());run=json.loads((D/'return.json').read_text());receipt=json.loads((D/'receipt.json').read_text());assert receipt['status']=='NORMAL_RETURN_ALL_OBSERVED_EXACT_STOP'and receipt['all_owned_exactSTOP']
geo=json.loads(gzip.decompress((R/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz').read_bytes()));N=len(geo['ids']);assert N==1896
with(D/'unique-controlled-offset-results.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==run['unique_results']==len(freeze['unique_admitted_cases'])
keys=[x['canonical_key']for x in rows];assert len(set(keys))==len(keys);bykey={x['canonical_key']:x for x in rows}
with np.load(D/'unique-controlled-offset-counts.npz',allow_pickle=False)as z:counts=z['counts'];assert counts.shape==(len(rows),2,N)and counts.dtype==np.uint32
assert hashlib.sha256(counts.tobytes()).hexdigest()==run['count_array_sha256']
for spi,sp in enumerate(['G','S']):
 den=np.array([t['D']for t in geo['species'][sp]['targets']],dtype=np.uint32);assert np.all(counts[:,spi,:]<=den[None,:])
for index,row in enumerate(rows):
 assert int(row['count_index'])==index
 r=[float(row[k])for k in ['r_GP','r_GS','r_PS']];T=float(row['T']);assert np.isfinite(r+[T]).all()and abs(sum(r)/3-T)<1e-14
expected_new=2*sum(x['origin']=='NEW_CONTROLLED_OFFSET_FIRST_EVALUATION'for x in rows);assert expected_new==run['new_rank_QR_projections']
progress=[json.loads(s)for s in(D/'unique-scenario-progress.jsonl').read_text().splitlines()];assert len(progress)==2*len(rows)
for i,k in enumerate(keys):assert progress[2*i]['canonical_key']==progress[2*i+1]['canonical_key']==k and progress[2*i]['state']=='STARTED'and progress[2*i+1]['state']=='RETURNED'
oldT=.22424191842389465;logical=[]
for case in freeze['cases']:
 rec=dict(case=case['case'],mode=case['mode'],distance_bp=case['distance_bp'],G_delta=case['G_delta'],S_delta=case['S_delta'],status='ADMITTED_SAVED_UNIQUE_REUSE'if case['admissible']else'NA_WHOLE_PANEL_GEOMETRY',canonical_key=case['canonical_key']or'',T=bykey[case['canonical_key']]['T']if case['admissible']else'NA',blocked_chromosomes=len(case['blocked_chromosomes']))
 logical.append(rec)
with Path('all72-logical-controlled-scenarios.tsv').open('x')as f:
 w=csv.DictWriter(f,fieldnames=list(logical[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(logical)
T=np.array([float(x['T'])for x in rows]);delta=T-oldT
summary=dict(status='ROOT_ALL_UNIQUE_COUNTS_AND_LOGICAL72_OUTPUT_QA_ACCEPTED',logical_scenarios=72,admitted=run['logical_admitted'],unique_results=len(rows),equivalent_mode_aliases=run['logical_admitted']-len(rows),NA_geometry=sum(not x['admissible']for x in freeze['cases']),max_admitted_distance_bp=max(int(x['distance_bp'])for x in rows),T_observed_reused=oldT,T_controlled_min=float(T.min()),T_controlled_max=float(T.max()),delta_T_min=float(delta.min()),delta_T_max=float(delta.max()),max_abs_delta_T=float(np.max(abs(delta))),positive_controlled_T=int(sum(T>0)),new_rank_QR_projections=run['new_rank_QR_projections'],old_saved_vector_reuses=run['old_saved_vector_reuses'],original_models_refitted=0,source_reads=0,new_pvalues=0,new_CIs=0,undefined=0,numerical_projections_reexecuted_in_QA=0,union_queries_reexecuted_in_QA=0,scope='Admissiblefixedpoints only reach1kb, shared100kb windows remain highlyoverlapping; cannotestimate fullwindowrange locationdecay',counts_cells_checked=int(counts.size))
Path('accepted-controlled-distance-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
report=f'''# 本轮空间定位补充：已接受结果

预先固定9个距离×4种G/S方向×2模式，共72个逻辑情景。完整1896成员及保护规则下16个标签通过，等价模式合并为8个实际新评价；其余56保持几何NA。通过距离只到1000bp，2500–100000bp均不允许，未改成员或放宽边界。

原T={oldT:.10f}复用。8个近距离情景T范围{T.min():.10f}–{T.max():.10f}，最大绝对变化{np.max(abs(delta)):.10f}。方向和每个情景逐项提供；不把不同方向支持的距离均值解释为同一条剂量曲线。

既有1998真实bp背景仅连接保存的T与几何：环形原窗口自重叠比例中位0.58408，线性0.77923；该量是整个窗口足迹的几何重叠，不是重复覆盖比例或可解释效应比例。全体originalwindows/D/X固定，不新增随机样本。

新评价只回答0–1kb的近位置敏感性。以现有100kb整窗和保护条件，无法给出完整窗口尺度的定位衰减检验。结果没有补齐S35源NA或P31区间，也不提升条件nominal p的校准地位。

科研源读取0，旧model/null/RM/QC/tree/HOG重复0；新固定响应rank/QR投影{run['new_rank_QR_projections']}次，零点和旧T未重算。所有NA、OG725及旧6644UNKNOWN保留。
'''
Path('空间定位特异性-接受结果.md').write_text(report)
Path('return.json').write_text(json.dumps(dict(status='CONTROLLED_ENDPOINTS_AND_72_LOGICAL_STATES_COMPLETE_QA_NO_REPLAY',summary='accepted-controlled-distance-summary.json',new_models=0,new_source_reads=0),indent=2)+'\n')
