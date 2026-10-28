"""Compact pure Figure5 SVG text, v2. V1 and its source validation unchanged.

No file IO, actual rendering, source acquisition, SQL or scientific computation.
Detailed source counts/identity are retained in metadata and pure supplement rows.
"""
import json
import xml.etree.ElementTree as ET
from v12_case_figure import FIXED_CASES,RADII,VERSIONS,COLORS,NS,need,validate_case

FONT=11.5


def checked(cases):
    need(isinstance(cases,(list,tuple))and len(cases)==5,'exact five fixed cases required')
    for case,(sid,name)in zip(cases,FIXED_CASES):validate_case(case,sid,name)
    return cases


def supplement_rows(cases):
    """Ten caller-writeable rows; model values explicitly separate from locus bp."""
    checked(cases);rows=[]
    for case in cases:
        for radius in RADII:
            w=case['windows'][str(radius)];c=case['r23'][str(radius)];f=case.get('focal_gene')
            row=dict(native_species_id=str(case['native_species_id']),species=case['species'],assembly=case['assembly'],old_group=case['old_group'],target_key=case['target_key'],center0=case['center0'],radius=radius,start0=w['start0'],end0=w['end0'],combined_mask_segments_json=json.dumps(w['masks'],separators=(',',':')),D_bp=w['D'],A_bp=w['A'],U_bp=w['U'],K_bp=w['K'],repeat_NA_reason=w.get('NA_reason'),focal_full_start0=f['start0']if f else None,focal_full_end0=f['end0']if f else None,R23_status=c['status'],R23_NA_reason=c.get('reason'),F_bp=c.get('F'),O_bp=c.get('O'),I_bp=c.get('I'),R23_D_bp=c.get('D'),selection_note=case.get('selection_note'),selection_is_outcome_blind=True)
            for version in VERSIONS:
                effect=case['all_J_effects'][version]
                for key,value in [('status',effect['status']),('beta',effect.get('beta')),('n',effect.get('n')),('NA_reason',effect.get('reason'))]:row['model_'+version+'_'+key]=value
            rows.append(row)
    return rows


def element(parent,tag,**attrs):
    return ET.SubElement(parent,'{'+NS+'}'+tag,{k.replace('_','-'):str(v)for k,v in attrs.items()})


def text(parent,x,y,value,size=FONT,fill='#283541',weight=None,anchor=None):
    need(size>=FONT,'print font floor11.5px')
    attrs=dict(x=x,y=y,font_size=size,fill=fill)
    if weight:attrs['font_weight']=weight
    if anchor:attrs['text_anchor']=anchor
    e=element(parent,'text',**attrs);e.text=str(value);return e


def rect(parent,x,y,w,h,color,**attrs):
    return element(parent,'rect',x=x,y=y,width=max(0,w),height=h,fill=color,**attrs)


def title(parent,value):element(parent,'title').text=str(value)


def geometry(parent,case,top):
    x,w=83,400;center=case['center0']
    def px(c):return x+(c-center+50000)/100000*w
    for kb in (-50,-25,0,25,50):
        xx=px(center+kb*1000)
        text(parent,xx,top+29,str(kb),anchor='middle',fill='#65717e')
    for index,radius in enumerate(RADII):
        c=case['windows'][str(radius)];y=top+36+index*16
        text(parent,53,y+7,'±'+str(radius//1000),anchor='end')
        rect(parent,px(c['start0']),y,px(c['end0'])-px(c['start0']),7,'#f0f3f6',stroke='#8896a5',stroke_width=.65,data_kind='frozen-window',data_radius=radius)
        for a,b in c['masks']:
            rect(parent,px(a),y,px(b)-px(a),7,COLORS['mask'],data_kind='combined-mask',data_radius=radius,data_start0=a,data_end0=b)
    y=top+68;text(parent,53,y+7,'gene',anchor='end')
    f=case.get('focal_gene')
    if f is None:text(parent,x,y+7,'NA',fill='#687685');return
    a=max(f['start0'],center-50000);b=min(f['end0'],center+50000)
    if a>=b:text(parent,x,y+7,'outside window',fill='#687685');return
    rect(parent,px(a),y,px(b)-px(a),7,COLORS['focal'],data_kind='focal-full-gene',data_start0=f['start0'],data_end0=f['end0'])
    if f['start0']<a:element(parent,'path',d='M %s %s l 4 -3.5 v 7 z'%(px(a),y+3.5),fill=COLORS['focal'])
    if f['end0']>b:element(parent,'path',d='M %s %s l -4 -3.5 v 7 z'%(px(b),y+3.5),fill=COLORS['focal'])


def strips(parent,case,top):
    xcomp,xrepeat,width=83,309,190
    text(parent,xcomp,top+94,'F / O / I',weight='600')
    text(parent,xrepeat,top+94,'U + K = A',weight='600')
    for index,radius in enumerate(RADII):
        y=top+101+index*18;text(parent,53,y+9,str(radius//1000),anchor='end')
        c=case['r23'][str(radius)]
        g=element(parent,'g',data_panel='composition',data_radius=radius,data_composition_status=c['status'])
        rect(g,xcomp,y,width,10,'#f0f2f5',stroke='#bbc4ce',stroke_width=.5)
        if c['status']=='TYPED_NA':
            text(g,xcomp+width/2,y+9,'NA',anchor='middle',fill='#687685');title(g,c['reason'])
        elif c['D']==0:
            text(g,xcomp+width/2,y+9,'NA (D=0)',anchor='middle',fill='#687685')
        else:
            cursor=xcomp
            for k in ('F','O','I'):
                piece=width*c[k]/c['D'];rect(g,cursor,y,piece,10,COLORS[k],data_component=k,data_bp=c[k]);cursor+=piece
            title(g,'Callable-bp proportions only; F=%s O=%s I=%s D=%s'%(c['F'],c['O'],c['I'],c['D']))
        w=case['windows'][str(radius)]
        g=element(parent,'g',data_panel='repeat',data_radius=radius,data_denominator_bp=w['D'])
        rect(g,xrepeat,y,width,10,'#eef1f4',stroke='#bbc4ce',stroke_width=.5)
        if w['A']is None or w['D']==0:
            text(g,xrepeat+width/2,y+9,'NA',anchor='middle',fill='#687685');title(g,w.get('NA_reason','D=0; fraction undefined'))
        else:
            uw=width*w['U']/w['D'];kw=width*w['K']/w['D']
            rect(g,xrepeat,y,uw,10,COLORS['U'],data_component='U',data_bp=w['U'])
            rect(g,xrepeat+uw,y,kw,10,COLORS['K'],data_component='K',data_bp=w['K'])
            rect(g,xrepeat,y,uw+kw,10,'none',stroke='#283541',stroke_width=.8,data_component='A-union-outline')
            title(g,'Locus bp/D; A=%s U=%s K=%s D=%s; A=U+K'%(w['A'],w['U'],w['K'],w['D']))


def effects(parent,case,top):
    x=532
    for i,version in enumerate(VERSIONS):
        e=case['all_J_effects'][version];y=top+45+i*23
        if e['status']=='OK':
            value=version+' '+format(e['beta'],'+.3g')+' / '+str(e['n'])
        else:value=version+' NA'
        t=text(parent,x,y,value,size=12)
        title(t,json.dumps(dict(e,version=version),ensure_ascii=False))


def render_svg(cases):
    """Pure string,720x930; no physical render, external fonts or file writes."""
    checked(cases)
    root=ET.Element('{'+NS+'}svg',{'width':'720','height':'930','viewBox':'0 0 720 930','role':'img','font-family':'Arial, Helvetica, sans-serif'})
    title(root,'Five fixed species: verified locus geometry and original all-J model context')
    element(root,'desc').text='F/O/I bars show proportions, never other-gene or intergenic locations. Repeat U+K bars and R23 composition have separate meanings. Single-locus counts and original all-J model effects use different denominators.'
    element(root,'metadata').text=json.dumps(dict(schema='FIGURE5_COMPACT_V2_SOURCE_VALUES_AND_SUPPLEMENT',cases=cases,supplement_rows=supplement_rows(cases),genomic_positions_drawn=['window','accepted_combined_mask','accepted_focal_full_gene_or_NA'],other_gene_or_intergenic_positions_inferred=False,repeat_counts_assigned_to_FOI=False),ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False)
    rect(root,0,0,720,930,'white')
    text(root,20,25,'Fixed-species locus geometry and model context',size=17,weight='600')
    text(root,20,45,'Geometry: kb from center. Bar rows: ±50 / ±25 kb; proportions only.',fill='#526171')
    text(root,532,64,'Original all-J β / n',weight='600')
    for index,(case,(sid,name))in enumerate(zip(cases,FIXED_CASES)):
        top=74+150*index
        g=element(root,'g',data_case_native_id=sid,data_row_top=top,data_row_height=150,data_old_group=case['old_group'])
        element(g,'line',x1=18,y1=top-3,x2=702,y2=top-3,stroke='#d4dce5',stroke_width=.7)
        text(g,20,top+14,name,size=12.5,weight='600')
        assembly=case['assembly'];shown=assembly if len(assembly)<=29 else assembly[:28]+'…'
        t=text(g,260,top+14,shown,fill='#526171');title(t,assembly)
        title(g,'Native '+sid+'; old group '+case['old_group']+'; target '+case['target_key']+'; '+str(case.get('selection_note','Root-supplied outcome-blind locus')))
        geometry(g,case,top);strips(g,case,top);effects(g,case,top)
    element(root,'line',x1=18,y1=830,x2=702,y2=830,stroke='#bfc9d4',stroke_width=.7)
    for label,x,color in [('F focal gene',20,COLORS['F']),('O other genes',160,COLORS['O']),('I intergenic',320,COLORS['I']),('U Unknown',470,COLORS['U']),('K strict',590,COLORS['K'])]:
        rect(root,x,842,10,10,color);text(root,x+15,852,label)
    text(root,20,874,'Dark segments: exon ∪ non-ACGT mask. F/O/I bars are not genomic positions.')
    text(root,20,894,'U and K are disjoint; outline A=U+K. Locus D and model n denominators differ.')
    text(root,20,914,'β: SD/SD. Danio composition NA. Full-gene nonexon ≠ proven functional intron.')
    return ET.tostring(root,encoding='unicode',xml_declaration=False)
