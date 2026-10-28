"""V2 pure SVG: retention is accepted median of species-level fractions.
No file/network/data analysis, fitting or scientific statistic derivation.
All scientific values, labels, denominators and acceptance refs are caller data.
"""
from html import escape
import math

INK='#182c3c';BLUE='#16729b';RUST='#ad523b';MUTED='#607480';PALE='#eaf1f4'

def need(x,why):
    if not x:raise ValueError('HOLD: '+why)

def finite(x):
    need(type(x) in (int,float) and math.isfinite(x),'finite admitted summary number')
    return x

def nonnegative_int(x):
    need(type(x) is int and x>=0,'exact nonnegative admitted count')
    return x

def label(x):
    need(type(x) is str and bool(x) and len(x)<=300,'bounded explicit label')
    return x

def fmt(x):
    return 'NA' if x is None else format(finite(x),'.5g')

class SVG:
    def __init__(self,w,h,title,description):
        self.parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}" role="img" aria-labelledby="title desc">',
                    '<title id="title">'+escape(label(title))+'</title>',
                    '<desc id="desc">'+escape(description)+'</desc>',
                    f'<rect width="{w}" height="{h}" fill="#fff"/>']
    def text(self,x,y,s,size=15,color=INK,weight='normal',anchor='start'):
        self.parts.append(f'<text x="{x}" y="{y}" fill="{color}" font-family="Arial,Helvetica,sans-serif" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}">{escape(str(s))}</text>')
    def rect(self,x,y,w,h,fill=PALE):
        self.parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" fill="{fill}"/>')
    def line(self,x1,y1,x2,y2,color=MUTED,width=1):
        self.parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
    def circle(self,x,y,r,fill):
        self.parts.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}"/>')
    def finish(self):return '\n'.join(self.parts+['</svg>'])

def admitted(d):
    need(d.get('accepted_compact_summaries') is True,'Root summary acceptance gate')
    need(type(d.get('source_refs')) is list and len(d['source_refs'])>0,'explicit accepted summary refs')
    for r in d['source_refs']:
        need(type(r) is dict and type(r.get('path')) is str and type(r.get('sha256')) is str and len(r['sha256'])==64,'exact source ref')

def figure1(d):
    """Selection partition + parallel denominator cards; NO false nesting."""
    admitted(d);s=d['selection'];total=nonnegative_int(s['total'])
    parts=[('Selected',s['selected'],BLUE),('Eligible, unselected',s['eligible_unselected'],'#75aab8'),('Excluded',s['excluded'],'#c8d3d9')]
    need(total>0 and sum(nonnegative_int(n) for _,n,_ in parts)==total,'exact selection partition')
    cards=d['denominators'];need(1<=len(cards)<=8,'bounded parallel denominator cards')
    v=SVG(1180,800,'Figure 1. Coverage and selection',
      'Selection is one conserved partition. Species counts below are distinct scope or eligibility definitions, displayed in parallel and not as nested attrition. Source NA is a separate observation status.')
    v.text(45,44,d['title'],25,weight='bold');v.text(45,76,s['unit'],15,MUTED)
    v.text(45,113,f"{total:,} = "+' + '.join(f'{n:,}' for _,n,_ in parts),22,weight='bold')
    x=45
    for title,n,color in parts:
        width=1090*n/total;v.rect(x,134,width,46,color);x+=width
    for i,(title,n,color) in enumerate(parts):
        x=45+i*365;v.circle(x+6,213,6,color);v.text(x+22,218,f'{title}: {n:,}',17,weight='bold')
    v.text(45,266,'Distinct denominators — parallel definitions, not a nested funnel',18,weight='bold')
    for i,c in enumerate(cards):
        x=45+(i%3)*367;y=290+(i//3)*132
        n=nonnegative_int(c['n']);v.rect(x,y,345,115)
        v.text(x+16,y+35,f'{n:,}',29,BLUE,'bold');v.text(x+85,y+34,label(c['label']),15,weight='bold')
        v.text(x+16,y+67,label(c['definition']),13,MUTED)
        v.text(x+16,y+92,label(c['unit']),12,MUTED)
    n=d['source_NA'];y=710
    v.rect(45,y,1090,55,'#fbf0e9')
    v.text(62,y+23,f"Source NA: {nonnegative_int(n['n'])} {label(n['unit'])}",16,RUST,'bold')
    v.text(62,y+44,label(n['definition']),13,RUST)
    v.text(45,789,label(d['scope_note']),12,MUTED)
    return v.finish()

def paired_panel(v,x,y,title,p):
    """Draw accepted median absolute paired change, never effect endpoints."""
    rows=p['rows'];need(1<=len(rows)<=4,'bounded paired-summary rows')
    lo,hi=p['axis'];finite(lo);finite(hi)
    need(lo==0 and hi>0,'explicit nonnegative absolute-change axis')
    need(p['statistic']=='median_abs_paired_delta','exact accepted summary estimand')
    v.text(x,y,title,18,weight='bold');v.text(x,y+23,label(p['unit']),12,MUTED)
    left=x+165;right=x+490;top=y+57
    v.line(left,top-10,right,top-10)
    v.text(left,top-15,'0',11,MUTED);v.text(right,top-15,fmt(hi),11,MUTED,anchor='end')
    for i,r in enumerate(rows):
        yy=top+i*49;v.text(x,yy+4,label(r['label']),13)
        z=r['median_abs_paired_delta'];n=nonnegative_int(r['finite_pairs']);flips=nonnegative_int(r['sign_changes'])
        need(flips<=n,'sign changes bounded by finite paired availability')
        if z is None:
            v.text(left,yy+4,'NA — '+label(r['NA_reason']),12,RUST);continue
        finite(z);need(n>0 and 0<=z<=hi,'admitted absolute-change median in axis')
        end=left+z/hi*(right-left)
        v.rect(left,yy-5,max(0,end-left),10,BLUE);v.circle(end,yy,4,BLUE)
        v.text(left,yy+20,f'median |Δ| {fmt(z)}; pairs {n}; sign changes {flips}',11,MUTED)
    v.text(x,top+len(rows)*49+8,label(p['interpretation']),12,MUTED)


def figure4(d):
    admitted(d)
    v=SVG(1180,1000,'Figure 4. Sensitivity and dependency',
      'Precomputed paired summaries, median species-level common-member retention and distinct dependency diagnostics. Absolute paired changes and retention medians are descriptive; dependency blocks, nodes and effective information are different units. Working range remains NA.')
    v.text(45,43,d['title'],25,weight='bold')
    v.text(45,70,'Blue = accepted median absolute paired change; sign changes are reported separately; NA remains explicit.',12,MUTED)
    paired_panel(v,45,115,'A  Predictor × response: median |Δ|',d['paired_2x2'])
    paired_panel(v,625,115,'B  Fixed common members: median |Δ|',d['scale_change'])
    v.line(45,424,1135,424,'#d8e2e7')
    r=d['retention'];v.text(45,467,'C  Clade member retention',18,weight='bold')
    v.text(45,491,label(r['unit']),12,MUTED);need(1<=len(r['rows'])<=5,'bounded clade rows')
    for i,row in enumerate(r['rows']):
        yy=530+i*57;v.text(45,yy,label(row['label']),13)
        for j,(key,color) in enumerate([('all_common',MUTED),('representatives',BLUE)]):
            a=row[key];median=finite(a['median_retention'])
            need(0<=median<=1,'accepted median species retention fraction')
            ns=nonnegative_int(a['n_species']);nf=nonnegative_int(a['finite_species'])
            ng=nonnegative_int(a['common_group_count'])
            need(0<=nf<=ns and ns>0,'explicit species/finite availability, including zero finite fits')
            y=yy-12+j*20;v.rect(185,y,190*median,12,color)
            v.text(385,y+11,f"{fmt(median)}; G={ng}; S={ns}; finite={nf}",10,color)
    v.text(45,844,label(r['interpretation']),12,MUTED)
    dep=d['dependency'];nonnegative_int(dep['blocks'])
    need(0<=finite(dep['nodes_percent'])<=100,'node percentage units')
    need(finite(dep['score_effective'])>=0,'nonnegative effective information')
    finite(dep['mean_score_dependence'])
    v.text(625,467,'D  Distinct dependency diagnostics',18,weight='bold')
    entries=[('Overlap blocks',dep['blocks'],'blocks'),
             ('Nodes within reported scope',dep['nodes_percent'],'percent of nodes'),
             ('Score effective information',dep['score_effective'],'effective information units'),
             ('Mean score dependence',dep['mean_score_dependence'],'unitless score diagnostic')]
    for i,(title,value,unit) in enumerate(entries):
        yy=500+i*69;v.rect(625,yy,510,58)
        v.text(639,yy+22,title,14,weight='bold');v.text(1119,yy+23,fmt(value),21,BLUE,'bold',anchor='end')
        v.text(639,yy+45,unit,12,MUTED)
    need(dep['working_range'] is None,'working range remains unavailable')
    v.text(625,806,'Working range: NA',18,RUST,'bold')
    v.text(625,834,label(dep['NA_reason']),12,RUST)
    v.text(45,908,label(d['positive_result_note']),15,BLUE,'bold')
    v.text(45,938,label(d['limits_note']),12,MUTED)
    v.text(45,969,label(d['scope_note']),12,MUTED)
    return v.finish()
