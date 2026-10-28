"""One new Figure7 from accepted saved values; zero new scientific evaluations."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
Q=Path(__import__('os').environ['ANRL_R4_ROOT'])
qa=json.loads((Q/'controlled-output-QA-first/accepted-controlled-distance-summary.json').read_text());assert qa['status'].endswith('ACCEPTED')
with(Q/'geometry-and-saved1998-first/saved1998-distance-and-footprint-T.tsv').open()as f:saved=list(csv.DictReader(f,delimiter='\t'))
with(Q/'controlled-distances-first/unique-controlled-offset-results.tsv').open()as f:points=list(csv.DictReader(f,delimiter='\t'))
freeze=json.loads((Q/'geometry-and-saved1998-first/fixed-distance-admission-freeze.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(7.2,5.7));a=fig.add_axes([.09,.48,.39,.36]);b=fig.add_axes([.60,.48,.37,.36]);c=fig.add_axes([.16,.14,.81,.20])
for mode,label,color in [('circ','Circular','#267c9f'),('linear','No wrap','#ca7d25')]:
 rows=[x for x in saved if x['mode']==mode];assert len(rows)==999
 a.scatter([float(x['whole_window_self_overlap_fraction'])*100 for x in rows],[float(x['T_reused'])for x in rows],s=4,color=color,alpha=.30,rasterized=True,label=label)
a.axhline(qa['T_observed_reused'],color='#a33434',lw=1,ls='--');a.set_xlabel('Original whole-window footprint retained (%)');a.set_ylabel('Saved phase-background T');a.set_title('A  Placement background and footprint overlap',loc='left',fontsize=8.5,fontweight='bold');a.legend(loc='lower right',frameon=False,markerscale=2,fontsize=7)
colors=['#247ba0','#9a6f3c','#4e8c55','#98558c'];patterns=[(1,1),(1,-1),(-1,1),(-1,-1)]
for (gs,ss),color in zip(patterns,colors):
 rows=sorted([x for x in points if int(x['G_delta'])*gs>0 and int(x['S_delta'])*ss>0],key=lambda x:int(x['distance_bp']))
 xs=[0]+[int(x['distance_bp'])for x in rows];ys=[qa['T_observed_reused']]+[float(x['T'])for x in rows]
 b.plot(xs,ys,marker='o',ms=3,lw=1,color=color,label=f'G{"+"if gs>0 else"−"}, S{"+"if ss>0 else"−"}')
b.axhline(qa['T_observed_reused'],color='#a33434',lw=.8,ls='--');b.set_xticks([0,250,500,1000]);b.set_xlabel('Absolute fixed genomic shift (bp)');b.set_ylabel('Controlled T');b.set_title('B  Admitted near-position profiles',loc='left',fontsize=8.5,fontweight='bold');b.legend(frameon=False,fontsize=6.7,ncol=2,loc='best');b.ticklabel_format(axis='y',style='plain',useOffset=False)
matrix=np.array([[next(int(x['admissible'])for x in freeze['cases']if x['mode']=='linear'and x['distance_bp']==d and np.sign(x['G_delta'])==gs and np.sign(x['S_delta'])==ss)for d in freeze['grid_bp']]for gs,ss in patterns]);assert np.array_equal(matrix,np.array([[next(int(x['admissible'])for x in freeze['cases']if x['mode']=='circ'and x['distance_bp']==d and np.sign(x['G_delta'])==gs and np.sign(x['S_delta'])==ss)for d in freeze['grid_bp']]for gs,ss in patterns]))
c.imshow(matrix,aspect='auto',cmap=ListedColormap(['#e7e7e7','#267c9f']),vmin=0,vmax=1,interpolation='nearest')
for i in range(4):
 for j in range(9):c.text(j,i,'PASS'if matrix[i,j]else'NA',ha='center',va='center',fontsize=6,color='white'if matrix[i,j]else'#666666')
c.set_yticks(range(4),[f'G{"+"if g>0 else"−"}, S{"+"if s>0 else"−"}'for g,s in patterns]);c.set_xticks(range(9),['0.25','0.5','1','2.5','5','10','25','50','100']);c.set_xlabel('Predeclared absolute genomic displacement (kb)');c.set_title('C  Whole-panel admission, both modes (all 1,896 groups retained)',loc='left',fontsize=8.5,fontweight='bold');c.tick_params(length=0)
fig.suptitle('Spatial placement constraints and near-position sensitivity',y=.97,fontsize=12,fontweight='bold');fig.text(.09,.9,'Saved backgrounds retain overlapping windows; controlled offsets reach only 1 kb.',fontsize=8.5)
fig.text(.09,.035,'16 admitted mode labels = 8 distinct evaluations; gray NA cells reflect protected geometry.',fontsize=8);fig.text(.09,.008,'Genomic signs denote coordinate directions. Zero T and all background draws are reused.',fontsize=8)
fig.savefig('Figure-7-location-specificity.png',dpi=300);fig.savefig('Figure-7-location-specificity.svg');plt.close(fig)
Path('return.json').write_text(json.dumps(dict(status='FIGURE7_FROM_ACCEPTED_SAVED_VALUES_PENDING_VISUAL_QA',new_evaluations=0,new_source_reads=0),indent=2)+'\n')
