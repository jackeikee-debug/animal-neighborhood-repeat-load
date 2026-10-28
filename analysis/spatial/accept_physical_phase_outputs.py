"""New complete saved-output QA only; no annotation/projection/model replay."""
from pathlib import Path
import csv,gzip,hashlib,json
import numpy as np
R3=Path(__import__('os').environ['ANRL_R3_ROOT'])
geo=json.loads(gzip.decompress((R3/'geometry-design-preparation-syntaxfix-once/geometry-and-mask-design.json.gz').read_bytes()))
D=np.array([[t['D']for t in geo['species'][sp]['targets']]for sp in ['G','S']],dtype=np.uint32)
accepted=[];refs=[]
for mode in ['circ','linear']:
    d=R3/(mode+'-all999-first');summary=json.loads((d/'return.json').read_text());receipt=json.loads((d/'receipt.json').read_text())
    assert receipt['status']=='NORMAL_RETURN_ALL_OBSERVED_EXACT_STOP'and receipt['all_owned_exactSTOP']and summary['draws']==999 and summary['undefined_draws']==summary['redraws']==0
    with(d/'all999-conditional-background-statistics.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==999 and [int(r['draw'])for r in rows]==list(range(1,1000))
    values=[]
    for r in rows:
        pairs=[float(r[k])for k in ['r_GP','r_GS','r_PS']];T=float(r['T']);assert all(np.isfinite(pairs))and all(abs(x)<=1+1e-14 for x in pairs)
        assert abs(T-sum(pairs)/3)<1e-14 and int(r['exceeds_observed'])==int(T>=summary['observed_T_reused'])
        assert float(r['G_residual_norm'])>0 and float(r['S_residual_norm'])>0;values.append(T)
    vals=np.array(values);exceed=sum(x>=summary['observed_T_reused']for x in values)
    assert exceed==summary['exceedances'] and (1+exceed)/1000==summary['conditional_nominal_upper_tail']
    for key,q in [('T_q025',.025),('T_median',.5),('T_q975',.975)]:assert abs(float(np.quantile(vals,q))-summary[key])<1e-14
    with np.load(d/'all999-two-species-frozen-window-counts.npz',allow_pickle=False)as z:
        counts=z['counts'];assert counts.shape==(999,2,1896)and counts.dtype==np.dtype('uint32')and np.all(counts<=D[None,:,:])
        assert hashlib.sha256(counts.tobytes()).hexdigest()==summary['count_array_sha256']
    with(d/'draw-progress.jsonl').open()as f:progress=[json.loads(x)for x in f]
    assert len(progress)==1998 and all(progress[2*i]['state']=='STARTED'and progress[2*i+1]['state']=='RETURNED'and progress[2*i]['draw']==progress[2*i+1]['draw']==i+1 for i in range(999))
    accepted.append(summary)
    for name in ['return.json','receipt.json','all999-conditional-background-statistics.tsv','all999-two-species-frozen-window-counts.npz']:
        f=d/name;refs.append(dict(path=str(f),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
Path('accepted-two-mode-physical-backgrounds.json').write_text(json.dumps(dict(status='ALL1998_TRUE_BP_DRAW_OUTPUTS_AND_ALL7576416_INTEGER_COUNTS_ACCEPTED',modes=accepted,saved_count_cells=999*2*1896*2,rank_QR_projections=3996,old_models_refitted=0,source_body_reads=0,union_query_reexecutions=0,all_source_NA_unchanged=True,OG725_unresolved=True,old6644UNKNOWN1_untouched=True,outputs=refs),indent=2)+'\n')
report='''# 真实碱基位移背景的接受结果

原A50观察统计量T=0.2242419184，原1896组、三个物种和原观察模型完全复用。新分析仅将G/S既有A类注释轨迹按真实bp位移，固定目标ALLexon/nonACGT掩膜、分母和原X列空间；P保存残差固定。保护/旧207/异常来源的反向窗口位置在真实轨迹读取前用于冻结允许相位。没有Genome/GFF/P重读或RepeatMasker/原模型重算。

|预定模式|draws|背景中位T|背景2.5–97.5百分位|最大T|达到观察值|条件名义上尾|
|---|---:|---:|---|---:|---:|---:|
'''
for s in accepted:report+=f"|{s['mode']}|999|{s['T_median']:.8f}|{s['T_q025']:.8f}–{s['T_q975']:.8f}|{s['T_max']:.8f}|{s['exceedances']}/999|{s['conditional_nominal_upper_tail']:.3f}|\n"
report+='''
共3996次新的rank/QR投影，所有1998个draw均完整、无undefined、丢弃或重抽。两套结果均高于全部新背景；这增加了原关联相对于这两种条件空间放置模型的证据。线性不绕回背景接近原观察值，说明保留局部位置约束的空间背景可解释部分正向对应。

保护/既往使用排除使circ的20/45条G/S染色体、linear的33/45条在已固定999次offset中没有达到100kb的位移。该100kb只作原窗口长度尺度诊断，不据此重新抽样。两模式不是独立复验，也不构成可自动套用精确置换定理的群。0.001是999次比较的分辨率下限，保留条件名义空间背景解释；不能由此声称所有空间混杂已消除、得到动物总体CI、或校准原S35共享依赖区间。

两模式来源、offset、counter、完整整数count数组/CRC/SHA、统计量和原T身份均经一次新输出验收。旧v12和先前校准记录保留；本分析与更早CNE失败校准区分。所有源NA/OG0000725和旧6644UNKNOWN1不变。
'''
Path('真实碱基位移-两模式接受结果.md').write_text(report)
