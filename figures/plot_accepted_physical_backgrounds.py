"""One scientific display from accepted saved statistics, no fits or new draws."""
from pathlib import Path
import csv,gzip,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R3=Path(__import__('os').environ['ANRL_R3_ROOT'])
O=Path(__import__('os').environ['ANRL_ORIGINAL_STATS_ROOT'])
qa=json.loads((R3/'physical-output-complete-QA-first/accepted-two-mode-physical-backgrounds.json').read_text());assert len(qa['modes'])==2
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
with gzip.open(O/'Stage-B-main-comparison-v1/main_50kb-null.tsv.gz','rt')as f:old=[float(x['T'])for x in csv.DictReader(f,delimiter='\t')]
curves=[('Saved index rotations',old,'#7b8188','--')]
for mode,label,color in [('circ','True bp circular shifts','#267c9f'),('linear','True bp no wrap shifts','#ca7d25')]:
    with(R3/(mode+'-all999-first')/'all999-conditional-background-statistics.tsv').open()as f:vals=[float(x['T'])for x in csv.DictReader(f,delimiter='\t')]
    assert len(vals)==999;curves.append((label,vals,color,'-'))
T=qa['modes'][0]['observed_T_reused'];fig=plt.figure(figsize=(7.2,4.8));ax=fig.add_axes([.11,.25,.85,.56])
for label,values,color,linestyle in curves:
    y=np.arange(1,len(values)+1)/len(values);ax.step(sorted(values),y,where='post',lw=1.7,color=color,ls=linestyle,label=label)
ax.axvline(T,color='#a33434',lw=1.6);ax.text(T+.004,.43,'Observed T = 0.2242',rotation=90,va='center',ha='left',color='#a33434',fontsize=9)
allvalues=[x for _,v,_,_ in curves for x in v];ax.set_xlim(min(allvalues)-.01,T+.025);ax.set_ylim(0,1.025)
ax.set_xlabel('Mean pairwise adjusted residual correlation T',labelpad=7);ax.set_ylabel('Empirical cumulative probability');ax.grid(axis='y',color='#e8e8e8',lw=.6)
fig.suptitle('Repeat track positions under physical displacement backgrounds',y=.965,fontsize=12,fontweight='bold')
fig.legend(*ax.get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.905),ncol=3,frameon=False,fontsize=8.8,handlelength=2,columnspacing=1.1)
fig.text(.11,.125,'Fixed original 1,896 groups, target masks and nuisance design; P residuals held fixed.',fontsize=8.6)
fig.text(.11,.083,'Each new background has 0 / 999 exceedances; conditional nominal upper tail = 0.001.',fontsize=8.6)
fig.text(.11,.041,'Protection limited phase sets; the no wrap check is not an independent replication.',fontsize=8.6)
fig.savefig('Figure-6-physical-backgrounds.png',dpi=300);fig.savefig('Figure-6-physical-backgrounds.svg');plt.close(fig)
outputs=[dict(path=q.name,bytes=q.stat().st_size,sha256=hashlib.sha256(q.read_bytes()).hexdigest())for q in Path('.').glob('Figure-6-*')]
Path('return.json').write_text(json.dumps(dict(status='FIGURE6_FROM_ACCEPTED_SAVED_STATISTICS_PENDING_VISUAL_QA',outputs=outputs,old_index_rotations_reexecuted=0,new_null_draws=0,new_models=0),indent=2)+'\n')
