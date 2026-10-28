"""Compact720px pure SVG layouts. Same v2 scientific data schema.
No IO/research. Root executes fixtures/render/visual QA after review.
"""
import json
from html import escape
from v12_summary_figures_v2 import SVG, admitted, need, finite, nonnegative_int, label, fmt, BLUE, RUST, MUTED

def short_lines(v,x,y,text,width=25):
    text=label(text);words=text.split();lines=['']
    for word in words:
        if len(lines[-1])+len(word)+1>width and lines[-1]:lines.append(word)
        else:lines[-1]+=(' ' if lines[-1] else '')+word
    need(len(lines)<=2,'provide a short display label of at most2lines; details belong in caption')
    for i,line in enumerate(lines):v.text(x,y+14*i,line,12)

def figure1(d):
    admitted(d);s=d['selection'];total=nonnegative_int(s['total'])
    parts=[('Selected',s['selected'],BLUE),('Eligible, unselected',s['eligible_unselected'],'#75aab8'),('Excluded',s['excluded'],'#b7c6cf')]
    need(total>0 and sum(nonnegative_int(n) for _,n,_ in parts)==total,'exact partition')
    v=SVG(720,650,'Figure1. Coverage and selection','Parallel denominators are distinct definitions, not nested attrition. Full definitions and sourceNA provenance remain in formal caption and metadata.')
    v.text(24,31,'Coverage and selection',20,weight='bold')
    v.text(24,64,f"{total:,} orthogroups",17,weight='bold');x=24
    for title,n,color in parts:
        w=672*n/total;v.rect(x,82,w,32,color);x+=w
    for i,(title,n,color) in enumerate(parts):
        x=24+i*228;v.text(x,142,title,12,color,'bold');v.text(x,165,f'{n:,}',19,color,'bold')
    v.text(24,205,'Distinct species denominators',16,weight='bold')
    need(len(d['denominators'])<=6,'compact6denominator cards')
    for i,c in enumerate(d['denominators']):
        x=24+(i%2)*344;y=226+(i//2)*101
        v.rect(x,y,328,88);v.text(x+12,y+34,nonnegative_int(c['n']),25,BLUE,'bold')
        short_lines(v,x+74,y+30,c['label'],29)
        v.text(x+12,y+71,'Parallel scope / eligibility definition',12,MUTED)
    n=d['source_NA'];v.rect(24,548,672,64,'#fbf0e9')
    v.text(36,574,f"Source NA: {nonnegative_int(n['n'])} anchor-by-species units",15,RUST,'bold')
    v.text(36,597,'Retained as unavailable; separate from species counts.',12,RUST)
    v.text(24,636,'Scope definitions and source lineage: formal caption and metadata.',12,MUTED)
    v.parts.append('<metadata>'+escape(json.dumps(d,ensure_ascii=False,sort_keys=True,allow_nan=False))+'</metadata>')
    return v.finish()

def compact_changes(v,x,y,title,p):
    need(p['statistic']=='median_abs_paired_delta','accepted median abs change')
    lo,hi=p['axis'];need(lo==0 and finite(hi)>0,'nonnegative change axis')
    need(1<=len(p['rows'])<=4,'max4summaryrows')
    v.text(x,y,title,15,weight='bold');v.text(x,y+22,'Median |Δ|; pairs / sign changes',12,MUTED)
    for i,r in enumerate(p['rows']):
        yy=y+52+i*52;short_lines(v,x,yy,r['label'],19)
        a=r['median_abs_paired_delta'];n=nonnegative_int(r['finite_pairs']);f=nonnegative_int(r['sign_changes'])
        need(f<=n,'finite paired signchanges')
        if a is None:v.text(x+158,yy,'NA',13,RUST);continue
        need(n>0 and 0<=finite(a)<=hi,'accepted median withinaxis')
        v.rect(x+158,yy-9,148*a/hi,10,BLUE)
        v.text(x+158,yy+17,f'{fmt(a)}; {n} / {f}',12)
    v.text(x,y+267,'Independent axis: 0 to '+fmt(hi),12,MUTED)

def figure4(d):
    admitted(d)
    v=SVG(720,900,'Figure4. Sensitivity and dependency','Accepted median absolute paired changes and median species-level retention; distinct dependency metrics. Detailed units, unavailable-result reasons and source refs remain in formal caption and metadata.')
    v.text(24,29,'Sensitivity and dependency',20,weight='bold')
    compact_changes(v,24,62,'A  Predictor × response',d['paired_2x2'])
    compact_changes(v,374,62,'B  Common-member scale',d['scale_change'])
    v.line(24,345,696,345,'#d2dfe5')
    v.text(24,373,'C  Median species-level retention',15,weight='bold')
    v.text(24,396,'Gray: all-common; blue: representatives. Labels: median; G / S / finite.',12,MUTED)
    rows=d['retention']['rows'];need(1<=len(rows)<=5,'compact5clades')
    for i,r in enumerate(rows):
        yy=428+i*43;short_lines(v,24,yy,r['label'],16)
        for j,(k,color) in enumerate([('all_common',MUTED),('representatives',BLUE)]):
            a=r[k];m=finite(a['median_retention']);ns=nonnegative_int(a['n_species']);nf=nonnegative_int(a['finite_species']);g=nonnegative_int(a['common_group_count'])
            need(0<=m<=1 and 0<=nf<=ns and ns>0,'median retention, including zero finite species')
            y=yy-9+j*18;v.rect(163,y,220*m,10,color)
            v.text(399,y+9,f'{fmt(m)}; {g} / {ns} / {nf}',12,color)
    v.line(24,650,696,650,'#d2dfe5');v.text(24,678,'D  Distinct dependency diagnostics',15,weight='bold')
    dep=d['dependency'];nonnegative_int(dep['blocks']);need(0<=finite(dep['nodes_percent'])<=100,'nodepercent')
    entries=[('Overlap blocks',dep['blocks'],'blocks'),('Nodes in scope',dep['nodes_percent'],'percent nodes'),
             ('Score effective',dep['score_effective'],'information units'),('Mean score dependence',dep['mean_score_dependence'],'unitless')]
    for i,(title,a,unit) in enumerate(entries):
        x=24+(i%2)*344;y=694+(i//2)*66;finite(a)
        v.rect(x,y,328,57);v.text(x+11,y+20,title,12,weight='bold')
        v.text(x+11,y+43,fmt(a),18,BLUE,'bold');v.text(x+120,y+42,unit,12,MUTED)
    need(dep['working_range'] is None,'workingrange unavailable')
    v.text(24,850,'Working range: NA',15,RUST,'bold')
    v.text(24,877,'Full units, positive findings, limits and source lineage: formal caption / metadata.',12,MUTED)
    v.parts.append('<metadata>'+escape(json.dumps(d,ensure_ascii=False,sort_keys=True,allow_nan=False))+'</metadata>')
    return v.finish()
