"""ReportLab scientific charts from saved results only: no fitting or new draws."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.graphics.shapes import Drawing, String, Line
from reportlab.graphics.charts.lineplots import ScatterPlot, LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.graphics import renderPDF
from reportlab.lib import colors
import pypdfium2 as pdfium

O=Path(__file__).resolve().parent
S=json.loads((O/'comparison-summary.json').read_text())
BLUE=colors.HexColor('#286785'); RED=colors.HexColor('#b04435'); GRAY=colors.HexColor('#66727c')
W,H=842,595
pdf=O/'three-species-main-and-sensitivity-figures-v2.pdf'
if pdf.exists():raise FileExistsError(pdf)


def text(d,x,y,t,size=9,color=colors.black):
    d.add(String(x,y,t,fontName='Helvetica',fontSize=size,fillColor=color))


def page_header(d,title,subtitle):
    text(d,38,563,title,17)
    text(d,38,544,subtitle,9,GRAY)
    text(d,38,24,'G: G. talaboti   P: P. flava   S: S. kowalevskii | Fixed panel; annotated interspersed repeats include Unknown.',8,GRAY)


def scatter(d,x,y,a,b,label,stat,rank=False):
    p=ScatterPlot()
    p.lineLabelFormat=None
    p.xLabel='';p.yLabel=''
    p.x=x;p.y=y;p.width=208;p.height=160
    p.data=[list(zip(map(float,a),map(float,b)))]
    p.joinedLines=0
    marker=makeMarker('FilledCircle');marker.size=1.65;marker.fillColor=colors.HexColor('#749dae');marker.strokeColor=None
    p.lines[0].symbol=marker
    for axis in (p.xValueAxis,p.yValueAxis):
        axis.labels.fontName='Helvetica';axis.labels.fontSize=7
        axis.strokeColor=GRAY;axis.strokeWidth=.6
        axis.valueMin=0 if rank else -3;axis.valueMax=1 if rank else 3
        axis.valueStep=.25 if rank else 1
        axis.labelTextFormat='%.2f' if rank else '%.0f'
    d.add(p)
    text(d,x,y+176,label+'   '+stat,10)
    text(d,x+55,y-27,label[0]+' '+('rank fraction' if rank else 'residual / SD'),8)
    text(d,x-18,y+164,label[-1],8)


c=canvas.Canvas(str(pdf),pagesize=(W,H),pageCompression=1)
c.setTitle('Three-species fixed neighborhood comparison - main and sensitivity')
c.setAuthor('P. flava / hemichordate research project')
for summary in S['versions']:
    v=summary['version']; n=summary['n_final']
    df=pd.read_csv(O/(v+'-saved-values.tsv.gz'),sep='\t')
    d=Drawing(W,H)
    title='Figure 1A. Main comparison: TSS +/-50 kb' if v=='main_50kb' else 'Figure 1B. Fixed sensitivity: TSS +/-25 kb'
    page_header(d,title,f'{n:,} orthogroup triplets | Every point retained; no line fitted to these panels.')
    text(d,38,519,'Raw correspondence (Spearman)',11,BLUE)
    text(d,38,284,'After the fixed structural adjustment (Pearson of rank-regression residuals)',11,BLUE)
    for k,(a,b) in enumerate([('G','P'),('G','S'),('P','S')]):
        pair=a+b;x=61+268*k
        scatter(d,x,335,df[a+'_response_average_rank']/(n+1),df[b+'_response_average_rank']/(n+1),a+' / '+b,
                'rho = %.3f'%summary['raw_Spearman'][pair],rank=True)
        ra=df[a+'_residual'];rb=df[b+'_residual']
        scatter(d,x,98,ra/ra.std(ddof=0),rb/rb.std(ddof=0),a+' / '+b,
                'r = %.3f'%summary['residual_Pearson'][pair])
    text(d,38,50,'Adjustment: ranked GC, gene count, target gene length, exon fraction, non-ACGT fraction, and chromosome effects.',8,GRAY)
    renderPDF.draw(d,c,0,0);c.showPage()

d=Drawing(W,H)
page_header(d,'Figure 2. Spatial background and chromosome influence','Saved results only | 999 chromosome-wise cyclic shifts per scale; P held fixed; seed 20261013.')
for k,summary in enumerate(S['versions']):
    v=summary['version']; x=62+400*k; y=327
    null=pd.read_csv(O/(v+'-null.tsv.gz'),sep='\t')['T'].to_numpy()
    hist,edges=np.histogram(null,bins=30)
    pts=[(float(edges[0]),0.)]
    for a,b,h in zip(edges[:-1],edges[1:],hist):pts.extend([(float(a),float(h)),(float(b),float(h))])
    pts.append((float(edges[-1]),0.))
    p=LinePlot();p.x=x;p.y=y;p.width=302;p.height=158;p.data=[pts];p.joinedLines=1
    p.lines[0].strokeColor=BLUE;p.lines[0].strokeWidth=1.3
    p.xValueAxis.valueMin=-.08;p.xValueAxis.valueMax=.30;p.xValueAxis.valueStep=.1;p.xValueAxis.labelTextFormat='%.2f'
    p.yValueAxis.valueMin=0;p.yValueAxis.valueMax=float(np.ceil(hist.max()/20)*20);p.yValueAxis.valueStep=40
    for axis in (p.xValueAxis,p.yValueAxis):axis.labels.fontSize=8;axis.strokeColor=GRAY
    d.add(p)
    observed=summary['T']; xx=x+(observed+.08)/.38*302
    d.add(Line(xx,y,xx,y+158,strokeColor=RED,strokeWidth=1.7))
    text(d,x,y+180,('Main +/-50 kb' if k==0 else 'S1 +/-25 kb')+f"  (n={summary['n_final']:,})",11)
    text(d,xx-20,y+166,'Observed',8,RED)
    text(d,x+58,y-25,'T = mean of the three residual correlations',8)
    text(d,x,y-48,'T = %.3f; empirical upper-tail p = %.3f'%(observed,summary['empirical_upper_tail_p']),10,RED)
    text(d,x,y-64,'Null 95%% range: [%.3f, %.3f]'%(summary['null_T_q025'],summary['null_T_q975']),9,GRAY)

text(d,38,229,'Fixed-residual chromosome deletion diagnostic',12)
text(d,38,211,'67 deletions per scale; no refitting. Bars are ranges, not confidence intervals.',9,GRAY)
left,width,lo,hi=250,430,.20,.28
for tick in [.20,.22,.24,.26,.28]:
    xx=left+(tick-lo)/(hi-lo)*width
    d.add(Line(xx,96,xx,183,strokeColor=colors.HexColor('#e4e8eb'),strokeWidth=.6))
    text(d,xx-10,83,'%.2f'%tick,8,GRAY)
for k,summary in enumerate(S['versions']):
    q=summary['fixed_residual_chromosome_influence'];y=162-k*43
    xa=left+(q['T_min']-lo)/(hi-lo)*width;xb=left+(q['T_max']-lo)/(hi-lo)*width;xc=left+(summary['T']-lo)/(hi-lo)*width
    text(d,52,y-3,'Main +/-50 kb' if k==0 else 'S1 +/-25 kb',10)
    d.add(Line(xa,y,xb,y,strokeColor=BLUE,strokeWidth=4))
    d.add(Line(xc,y-7,xc,y+7,strokeColor=RED,strokeWidth=2))
    text(d,699,y-3,'%.3f - %.3f'%(q['T_min'],q['T_max']),9)
text(d,38,57,'All three pairwise residual correlations remain positive after every deletion at both scales.',9)
text(d,38,42,'S2: no fixed window intersects mapped Hox core clusters; same dataset, no extra fit, not independent validation.',8,GRAY)
renderPDF.draw(d,c,0,0);c.showPage();c.save()

doc=pdfium.PdfDocument(str(pdf))
for i in range(len(doc)):
    page=doc[i]
    bitmap=page.render(scale=1.5)
    bitmap.to_pil().save(O/f'figure-v2-page-{i+1}.png')
    page.close()
doc.close()
print(json.dumps({'pdf':str(pdf),'pages':3,'new_fits':0,'new_random_draws':0}))
