"""Pure SVG text for the five fixed Figure5 species. No file IO or plotting run.

Source bindings and outcome-blind locus selection are Root inputs, not inferred
here. Only supplied window/focal/mask positions are drawn. R23 F/O/I are counts,
never invented other-gene/intergenic locations or repeat allocations. No models.
"""
import math
import xml.etree.ElementTree as ET

FIXED_CASES=(('26','Lineus longissimus'),('17','Ephydatia muelleri'),
             ('44','Sycon ciliatum'),('8','Bolinopsis microptera'),('14','Danio rerio'))
RADII=(50000,25000)
VERSIONS=('A50','K50','A25')
COLORS=dict(F='#238b8b',O='#d69b37',I='#d9dfe5',U='#9a4f96',K='#416da7',mask='#333c49',focal='#52745c')
NS='http://www.w3.org/2000/svg'
ET.register_namespace('',NS)


def need(value,message):
    if not value:raise ValueError(message)


def whole(v):return type(v)is int


def clipped(masks,a,b):
    return [[max(x,a),min(y,b)]for x,y in masks if max(x,a)<min(y,b)]


def validate_case(case,sid,name):
    need(str(case.get('native_species_id'))==sid,'fixed nativeID order is required; never reorder by effects')
    need(isinstance(case.get('species'),str) and ' '.join(case['species'].replace('_',' ').split()).startswith(name),'fixed species identity')
    need(all(isinstance(case.get(k),str)and case[k]for k in ('assembly','old_group','target_key')),'assembly/oldgroup/target identity')
    need(case.get('selection_is_outcome_blind')is True,'Root must supply outcome-blind locus selection')
    center=case['center0'];need(whole(center)and center>=50000,'valid existing complete window center')
    windows=case['windows'];need(set(windows)==set(map(str,RADII)),'exact two existing radii')
    for radius in RADII:
        w=windows[str(radius)]
        need(w['old_group']==case['old_group']and w['target_key']==case['target_key'],'same target and oldgroup across both radii')
        need(w['radius']==radius and whole(w['start0'])and whole(w['end0'])and (w['start0'],w['end0'])==(center-radius,center+radius),'frozen relative window geometry')
        need(w['complete']is True,'only complete accepted frozen window geometry can be positioned')
        right=w['start0'];maskbp=0
        for segment in w['masks']:
            need(len(segment)==2 and all(whole(v)for v in segment),'literal mask integer segments')
            a,b=segment;need(right<=a<b<=w['end0'],'ordered disjoint within-window combined masks');right=b;maskbp+=b-a
        need(whole(w['D'])and w['D']==2*radius-maskbp,'accepted mask/D closure')
        vals=[w[k]for k in ('A','U','K')]
        if any(v is None for v in vals):
            need(all(v is None for v in vals)and isinstance(w.get('NA_reason'),str)and w['NA_reason'],'all repeat counts explicitly NA with reason')
        else:
            need(all(whole(v)and 0<=v<=w['D']for v in vals)and w['A']==w['U']+w['K'],'A=U+K integer coverage closure')
    big,small=windows['50000'],windows['25000']
    need(clipped(big['masks'],small['start0'],small['end0'])==small['masks'],'same accepted two-radius mask clip required')
    focal=case.get('focal_gene')
    if focal is not None:need(all(whole(focal[k])for k in ('start0','end0'))and 0<=focal['start0']<focal['end0'],'literal original full-gene bounds or None')
    comp=case['r23'];need(set(comp)==set(map(str,RADII)),'two R23 radius states required')
    for radius in RADII:
        c=comp[str(radius)]
        if sid=='14':need(c.get('status')=='TYPED_NA','Danio ANN ledger composition must remain NA')
        if c.get('status')=='AVAILABLE':
            need(all(whole(c[k])and c[k]>=0 for k in ('F','O','I','D'))and c['F']+c['O']+c['I']==c['D']==windows[str(radius)]['D'],'R23 F+O+I=D and same locus denominator')
        else:
            need(c.get('status')=='TYPED_NA'and isinstance(c.get('reason'),str)and c['reason'],'explicit R23 typedNA reason')
            need(all(c.get(k)is None for k in ('F','O','I','D')),'R23 NA must not include invented zero counts')
    effects=case['all_J_effects'];need(set(effects)==set(VERSIONS),'exact original three model versions')
    for version in VERSIONS:
        e=effects[version];need(e['denominator']=='ORIGINAL_J','single-locus and original model denominators must stay separate')
        if e['status']=='OK':
            need(type(e['beta'])in(int,float)and math.isfinite(e['beta']),'finite original model effect')
            need(whole(e['n'])and e['n']>0,'original accepted model n')
        else:
            need(e['status']=='TYPED_NA'and e.get('beta')is None and e.get('reason'),'original model typedNA preserved')
    return case


def node(parent,tag,**attrs):
    return ET.SubElement(parent,'{'+NS+'}'+tag,{k.replace('_','-'):str(v)for k,v in attrs.items()})


def text(parent,x,y,value,size=12,fill='#273442',weight=None,anchor=None,**attrs):
    attrs.update(x=x,y=y,font_size=size,fill=fill)
    if weight:attrs['font_weight']=weight
    if anchor:attrs['text_anchor']=anchor
    n=node(parent,'text',**attrs);n.text=str(value);return n


def rect(parent,x,y,w,h,fill,**attrs):
    return node(parent,'rect',x=x,y=y,width=max(0,w),height=h,fill=fill,**attrs)


def shorten(value,n=73):
    return value if len(value)<=n else value[:n-1]+'…'


def label_with_title(parent,x,y,value,**attrs):
    n=text(parent,x,y,shorten(value),**attrs);node(n,'title').text=value;return n


def percent(n,d):return ('%.1f%%'%(100*n/d))if d else 'NA (D=0)'


def geometry_panel(parent,case,x,y,w):
    center=case['center0']
    def px(coord):return x+(coord-center+50000)/100000*w
    text(parent,x,y,'Single locus: verified geometry',size=13,weight='600')
    for offset in (-50000,-25000,0,25000,50000):
        p=px(center+offset);node(parent,'line',x1=p,y1=y+19,x2=p,y2=y+96,stroke='#e4e8ed',stroke_width=1)
        text(parent,p,y+12,str(offset//1000),size=10,anchor='middle')
    for index,radius in enumerate(RADII):
        window=case['windows'][str(radius)];top=y+27+index*27
        a,b=window['start0'],window['end0']
        rect(parent,px(a),top,px(b)-px(a),12,'#fafbfc',stroke='#8795a5',stroke_width=1,data_kind='frozen-window',data_radius=radius)
        for m,n in window['masks']:
            rect(parent,px(m),top,px(n)-px(m),12,COLORS['mask'],data_kind='combined-mask',data_start0=m,data_end0=n,data_radius=radius)
        text(parent,x-6,top+10,'±'+str(radius//1000),size=10,anchor='end')
    top=y+83;f=case.get('focal_gene')
    if f is None:
        text(parent,x,top+9,'full-gene bounds: NA',size=11,fill='#6a7582')
    else:
        left=max(f['start0'],center-50000);right=min(f['end0'],center+50000)
        if left<right:
            rect(parent,px(left),top,px(right)-px(left),9,COLORS['focal'],data_kind='focal-full-gene',data_start0=f['start0'],data_end0=f['end0'])
            if f['start0']<left:node(parent,'path',d='M %s %s l 5 -4 l 0 8 z'%(px(left),top+4.5),fill=COLORS['focal'])
            if f['end0']>right:node(parent,'path',d='M %s %s l -5 -4 l 0 8 z'%(px(right),top+4.5),fill=COLORS['focal'])
        else:text(parent,x,top+9,'full gene outside displayed window',size=10)
    text(parent,x,y+112,'kb from frozen center; dark segments = exon ∪ non-ACGT mask',size=9,fill='#5b6876')


def composition_panel(parent,case,x,y,w):
    text(parent,x,y,'R23 callable-bp composition (proportions, not positions)',size=11,weight='600')
    for index,radius in enumerate(RADII):
        c=case['r23'][str(radius)];top=y+12+index*41
        group=node(parent,'g',data_panel='composition',data_radius=radius,data_composition_status=c['status'])
        text(group,x-6,top+11,'±'+str(radius//1000),size=10,anchor='end')
        if c['status']=='TYPED_NA':
            rect(group,x,top,w,14,'#f2f3f5',stroke='#c3cad3',stroke_dasharray='3 2')
            text(group,x+5,top+11,'NA: composition source unavailable',size=10,fill='#6a7582')
            node(group,'title').text=c['reason'];continue
        D=c['D'];rect(group,x,top,w,14,'#f2f3f5');cursor=x
        for component in ('F','O','I'):
            width=w*c[component]/D if D else 0
            rect(group,cursor,top,width,14,COLORS[component],data_component=component,data_bp=c[component]);cursor+=width
        text(group,x,top+28,'F=%s  O=%s  I=%s; D=%s bp'%(c['F'],c['O'],c['I'],D),size=9)
    text(parent,x,y+91,'F focal full-gene nonexon  ·  O other genes  ·  I intergenic',size=9,fill='#5b6876')


def repeat_panel(parent,case,x,y,w):
    text(parent,x,y,'Single-locus repeat fractions: A = U + K',size=12,weight='600')
    text(parent,x,y+17,'U and strict K are disjoint; denominator = this locus D',size=10,fill='#5b6876')
    for index,radius in enumerate(RADII):
        c=case['windows'][str(radius)];top=y+31+index*58
        if c['A']is None:
            g=node(parent,'g',data_panel='repeat-fractions-NA',data_radius=radius)
            text(g,x,top+12,'±%s kb — Repeat counts: NA'%(radius//1000),size=11)
            node(g,'title').text=c['NA_reason'];continue
        D=c['D'];g=node(parent,'g',data_panel='repeat-fractions',data_radius=radius,data_denominator_bp=D,data_a_bp=c['A'],data_u_bp=c['U'],data_k_bp=c['K'])
        text(g,x-6,top+11,'±'+str(radius//1000),size=10,anchor='end')
        rect(g,x,top,w,16,'#eef1f4',stroke='#c2cad4',stroke_width=.6)
        uw=w*c['U']/D if D else 0;kw=w*c['K']/D if D else 0
        rect(g,x,top,uw,16,COLORS['U'],data_component='U')
        rect(g,x+uw,top,kw,16,COLORS['K'],data_component='K')
        if D:rect(g,x,top,uw+kw,16,'none',stroke='#25303b',stroke_width=1,data_component='A-union-outline')
        text(g,x,top+31,'A %s  U %s  K %s'%(percent(c['A'],D),percent(c['U'],D),percent(c['K'],D)),size=10)
        text(g,x,top+45,'D=%s bp; A=%s, U=%s, K=%s bp'%(D,c['A'],c['U'],c['K']),size=9,fill='#5b6876')
    text(parent,x,y+160,'No repeat counts are allocated to F/O/I regions.',size=10,fill='#5b6876')


def effects_panel(parent,case,x,y):
    text(parent,x,y,'Original all-J models',size=13,weight='600')
    text(parent,x,y+19,'β (SD/SD), original version-specific n',size=10,fill='#5b6876')
    for i,version in enumerate(VERSIONS):
        effect=case['all_J_effects'][version];top=y+46+i*35
        if effect['status']=='OK':
            text(parent,x,top,version+': '+format(effect['beta'],'+.4g'),size=14,weight='600')
            text(parent,x+155,top,'n='+str(effect['n']),size=11)
        else:
            text(parent,x,top,version+': NA',size=13)
            node(parent,'title').text=version+': '+str(effect['reason'])
    text(parent,x,y+159,'Not an effect estimated from this locus.',size=10,fill='#5b6876')
    text(parent,x,y+175,'Intervals are not supplied or invented.',size=10,fill='#5b6876')


def render_svg(cases):
    """Return a standalone SVG string. Input list is never mutated or reordered.

    Exactly five Root-bound cases in FIXED_CASES order; this module has no file
    IO, selection routine, inference, SQL, rendering backend or output creation.
    See v12_case_figure_INTERFACE.md for strict payload and missingness rules.
    """
    need(isinstance(cases,(list,tuple))and len(cases)==5,'exact five fixed cases required')
    for case,(sid,name)in zip(cases,FIXED_CASES):validate_case(case,sid,name)
    width,height=1500,1590
    svg=ET.Element('{'+NS+'}svg',{'width':str(width),'height':str(height),'viewBox':'0 0 %s %s'%(width,height),'role':'img','font-family':'Arial, Helvetica, sans-serif'})
    # ElementTree does not execute text; source labels are escaped by serialization.
    node(svg,'title').text='Fixed five-species case geometry, callable composition and original all-J effects'
    node(svg,'desc').text='Only supplied focal-gene and combined-mask positions are shown. Other-gene/intergenic counts have no implied positions. Single-locus counts and original all-J model effects use different denominators.'
    rect(svg,0,0,width,height,'#ffffff')
    text(svg,38,34,'Five fixed illustrative species: locus geometry and original model context',size=21,weight='600')
    text(svg,38,57,'Outcome-blind locus inputs; single-locus bp and all-J model denominators differ. No new effects or interval-specific repeat loads.',size=12,fill='#526171')
    for label,color,x in [('F',COLORS['F'],38),('O',COLORS['O'],99),('I',COLORS['I'],160),('U',COLORS['U'],221),('K',COLORS['K'],282)]:
        rect(svg,x,70,13,10,color);text(svg,x+18,80,label,size=10)
    text(svg,365,80,'F/O/I composition and U/K repeat loads are separate summaries.',size=10,fill='#526171')
    for index,(case,(sid,name))in enumerate(zip(cases,FIXED_CASES)):
        top=100+index*285;g=node(svg,'g',data_case_native_id=sid,data_old_group=case['old_group'])
        rect(g,25,top,1450,274,'#ffffff',stroke='#d8e0e8',stroke_width=1,rx=5)
        text(g,43,top+24,name+'  [native '+sid+']',size=16,weight='600')
        label_with_title(g,420,top+24,'Assembly '+case['assembly']+' · old group '+case['old_group'],size=11)
        node(g,'title').text='Species input: '+case['species']+'; assembly: '+case['assembly']+'; target: '+case['target_key']+'; '+str(case.get('selection_note','Root-supplied outcome-blind locus'))
        geometry_panel(g,case,87,top+47,515)
        composition_panel(g,case,87,top+173,515)
        repeat_panel(g,case,684,top+50,380)
        effects_panel(g,case,1123,top+50)
    text(svg,38,1547,'Only verified window/focal/mask positions are plotted. R23 bars show counts only; no O/I positions are inferred.',size=11,fill='#526171')
    text(svg,38,1567,'Danio ANN-ledger composition remains NA. Full-gene nonexonic sequence is not asserted to be a functional intron.',size=11,fill='#526171')
    return ET.tostring(svg,encoding='unicode',xml_declaration=False)
