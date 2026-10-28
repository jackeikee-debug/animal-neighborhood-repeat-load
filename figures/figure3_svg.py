"""Pure Figure 3 SVG compositor. No source reads, tree inference or statistics.

Root supplies accepted values and frozen children order. Rendering is explicit.
"""
from html import escape
import math
import json
from pathlib import Path

SCHEMA = 'FIG3_FROZEN_WORKING_TREE_V1'
ROLES = {'DISCOVERY': 'D', 'MAIN': 'M', 'METHOD_EXCEPTION': 'X',
         'GENE_TREE_ONLY': 'G', 'OUTGROUP': 'O'}
AVAIL = {'FINITE': 'F', 'SOURCE_NA': 'S', 'NO_QUALIFIED_ANCHORS': 'Q',
         'METHOD_NA': 'X', 'NOT_APPLICABLE': '-'}
WORK = {'PASS': 'W', 'INFORMATION_GATE_NA': 'I', 'SOURCE_NA': 'S',
        'NOT_APPLICABLE': '-'}


def _number(v, name):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(name + ' must be a finite number')
    return v


def _text(v, name, allow_na=False):
    if allow_na and v is None:
        return 'NA'
    if not isinstance(v, str) or not v or not v.isascii() or any(ord(c) < 32 for c in v):
        raise ValueError(name + ' requires nonempty printable ASCII')
    return v


def render_svg(data):
    """Return SVG text; reject ambiguous identities, dropped tips and clipped values."""
    if data.get('schema') != SCHEMA:
        raise ValueError('schema mismatch')
    order = data['leaf_order']
    rows = data['tips']
    if len(order) != 54 or len(set(order)) != 54 or set(rows) != set(order):
        raise ValueError('exact 54-tip order/row closure required')
    scenarios = data['mapping_scenario_ids']
    if len(scenarios) != 5 or len(set(scenarios)) != 5:
        raise ValueError('five exact mapping scenario IDs required')
    for value in order + scenarios:
        _text(value, 'ID')
    lo, hi = data['effect_axis']['min'], data['effect_axis']['max']
    _number(lo, 'axis.min'); _number(hi, 'axis.max')
    if not lo < 0 < hi:
        raise ValueError('common effect axis must contain zero')
    ticks = data['effect_axis']['ticks']
    if not ticks or ticks != sorted(set(ticks)):
        raise ValueError('explicit ordered unique ticks required')
    for t in ticks:
        _number(t, 'tick')
        if not lo <= t <= hi:
            raise ValueError('tick outside axis')
    seen_nodes = set(); visited = []; nodes = []
    def visit(n, depth):
        if id(n) in seen_nodes:
            raise ValueError('cycle or shared node object')
        seen_nodes.add(id(n))
        if depth > 108:
            raise ValueError('tree too deep')
        children = n.get('children', [])
        if not isinstance(children, list):
            raise ValueError('children must be ordered list')
        if 'branch_length' in n and n['branch_length'] is not None:
            if _number(n['branch_length'], 'branch_length') < 0:
                raise ValueError('negative branch length')
        if children:
            xy = [visit(c, depth + 1) for c in children]
            y = (xy[0][1] + xy[-1][1]) / 2
        else:
            key = _text(n.get('tip_id'), 'leaf.tip_id')
            visited.append(key)
            y = 60 + (len(visited) - 1) * 7.65
            xy = []
        # Equal depth edges explicitly schematic; provided branch lengths are not used.
        item = {'depth': depth, 'y': y, 'children': xy}
        nodes.append(item)
        return item, y
    root, _ = visit(data['tree'], 0)
    if visited != order:
        raise ValueError('tree traversal must equal frozen leaf_order exactly')
    maxdepth = max(n['depth'] for n in nodes)
    if maxdepth < 1:
        raise ValueError('invalid tree')
    identities = []
    def effect(record, label):
        status = record['status']
        if status == 'FINITE':
            value = _number(record['value'], label)
            if not lo <= value <= hi:
                raise ValueError(label + ' exceeds axis; widen frozen axis, never clip')
            return value
        if status not in {'SOURCE_NA', 'NO_QUALIFIED_ANCHORS', 'METHOD_NA', 'NOT_APPLICABLE'}:
            raise ValueError('unknown effect state')
        if record.get('value') is not None:
            raise ValueError('NA effect has numeric value')
        _text(record.get('reason'), label + '.reason')
        return None
    values = {}
    for key in order:
        row = rows[key]
        species = _text(row['species'], 'species')
        assembly = _text(row['assembly'], 'assembly', allow_na=True)
        if len(species) > 48:
            raise ValueError('species label exceeds reserved width; Root must expand layout')
        identities.append((species, assembly))
        if row['role'] not in ROLES or row['availability'] not in AVAIL or row['workrange_status'] not in WORK:
            raise ValueError('unrecognized role/availability/workrange')
        j = row['J']
        if j is None:
            _text(row.get('J_reason'), 'J_reason')
        elif isinstance(j, bool) or not isinstance(j, int) or not 0 <= j <= 1896:
            raise ValueError('invalid J')
        effects = [effect(row[k], key + '.' + k) for k in ('A50', 'K50', 'A25')]
        if (row['availability'] == 'FINITE') != (effects[0] is not None):
            raise ValueError('A50 availability inconsistent')
        if row['workrange_status'] == 'SOURCE_NA' and row['availability'] != 'SOURCE_NA':
            raise ValueError('workrange source NA inconsistent')
        if row['workrange_status'] in {'PASS', 'INFORMATION_GATE_NA'} and effects[0] is None:
            raise ValueError('numeric working-range state requires finite A50')
        dr = row['A50_deletion_range']
        if dr['status'] == 'FINITE':
            if effects[0] is None:
                raise ValueError('deletion range requires finite A50')
            lower = _number(dr['min'], 'deletion min'); upper = _number(dr['max'], 'deletion max')
            if not lo <= lower <= upper <= hi:
                raise ValueError('invalid/uncontained deletion range')
        else:
            if dr['status'] not in {'UNAVAILABLE', 'NOT_APPLICABLE'}:
                raise ValueError('unknown deletion state')
            if dr.get('min') is not None or dr.get('max') is not None:
                raise ValueError('NA deletion range cannot have bounds')
            _text(dr.get('reason'), 'deletion reason')
        env = row['mapping_own_A50']
        if env['status'] == 'FINITE':
            if env['scenario_ids'] != scenarios:
                raise ValueError('mapping envelope must cover exact five scenarios in frozen order')
            lower = _number(env['min'], 'mapping min'); upper = _number(env['max'], 'mapping max')
            if not lo <= lower <= upper <= hi:
                raise ValueError('invalid mapping envelope')
        else:
            if env['status'] not in {'INCOMPLETE_NA', 'NOT_APPLICABLE'}:
                raise ValueError('unknown mapping envelope state')
            if env.get('min') is not None or env.get('max') is not None:
                raise ValueError('NA mapping envelope cannot have bounds')
            _text(env.get('reason'), 'mapping reason')
        values[key] = effects
    if len(set(identities)) != 54:
        raise ValueError('duplicate species+assembly identity')
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="10in" height="7in" viewBox="0 0 720 504">',
           '<title>Conditional profiles on the original working topology</title>',
           '<desc>54 frozen tips. Tree edges schematic. A50 chromosome deletion ranges are not confidence intervals. Mapping envelope is scenario sensitivity, not uncertainty.</desc>',
           '<style>text{font-family:Arial,Helvetica,sans-serif;font-size:7.5px;fill:#20333e} .head{font-weight:bold} .muted{fill:#61727b}</style>',
           '<rect width="720" height="504" fill="white"/>']
    def text(x, y, s, cls='', anchor='start', title=None):
        extra = '' if title is None else '<title>' + escape(title) + '</title>'
        out.append(f'<text x="{x:.3f}" y="{y:.3f}" class="{cls}" text-anchor="{anchor}" dominant-baseline="middle">{extra}{escape(str(s))}</text>')
    def line(x1, y1, x2, y2, color='#9aabb3', width=.5):
        out.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{color}" stroke-width="{width}"/>')
    def point(x, y, color):
        out.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="1.6" fill="{color}"/>')
    def xpos(v, start, width):
        return start + (v - lo) / (hi - lo) * width
    text(10, 12, 'Conditional repeat profiles on the original working topology', 'head')
    text(10, 25, 'Original working topology; schematic edges | Reference hypotheses T1-T4 are separate sensitivity inputs.', 'muted')
    for x, label, anchor in [(10,'Tree','start'),(93,'Species (assembly identity in SVG metadata)','start'),
                             (281,'Role/Av/Wr','start'),(344,'J','end'),(392,'A50 + chr range','middle'),
                             (475,'K50','middle'),(540,'Original A25','middle'),(641,'Own mapping (5)','middle')]:
        text(x, 39, label, 'head', anchor)
    for start, width in [(350,84),(447,56),(512,56),(586,111)]:
        for tick in ticks:
            x = xpos(tick,start,width)
            line(x,56,x,470,'#b2bdc2' if tick==0 else '#e2e7ea',.6 if tick==0 else .3)
            text(x,49,f'{tick:g}','muted','middle')
    def draw_tree(n):
        x = 10 + 70*n['depth']/maxdepth
        if n['children']:
            line(x,n['children'][0][1],x,n['children'][-1][1])
            for child, y in n['children']:
                cx = 10 + 70*child['depth']/maxdepth
                line(x,y,cx,y); draw_tree(child)
        else:
            line(x,n['y'],87,n['y'],'#c7d1d6',.35)
    draw_tree(root)
    for i, key in enumerate(order):
        row = rows[key]; y = 60+i*7.65
        # Bands begin after tree, preserving branches and uninterrupted row alignment.
        if i % 2 == 0:
            out.append(f'<rect x="89" y="{y-3.75:.3f}" width="610" height="7.5" fill="#f3f6f7" fill-opacity="0.55"/>')
        identity = row['species'] + ' | assembly=' + ('NA' if row['assembly'] is None else row['assembly']) + ' | tip=' + key
        out.append('<g data-tip-id="'+escape(key,quote=True)+'" data-species="'+escape(row['species'],quote=True)+'" data-assembly="'+escape(row['assembly'] or 'NA',quote=True)+'">')
        out.append('<title>'+escape(identity)+'</title>')
        out.append('<metadata>'+escape(json.dumps(row,sort_keys=True,ensure_ascii=True))+'</metadata>')
        text(93,y,row['species'],title=identity)
        codes=ROLES[row['role']]+'/'+AVAIL[row['availability']]+'/'+WORK[row['workrange_status']]
        text(281,y,codes,title=row['role']+' / '+row['availability']+' / '+row['workrange_status'])
        text(344,y,'NA' if row['J'] is None else row['J'],anchor='end')
        for k, start, width, v in zip(('A50','K50','A25'),(350,447,512),(84,56,56),values[key]):
            if v is None:
                text(start+width/2,y,AVAIL[row[k]['status']],'muted','middle',row[k]['reason'])
            else:
                color = '#187a80' if v>0 else '#b24c37' if v<0 else '#20333e'
                if k=='A50' and row['A50_deletion_range']['status']=='FINITE':
                    dr=row['A50_deletion_range']
                    line(xpos(dr['min'],start,width),y,xpos(dr['max'],start,width),y,color,1.1)
                point(xpos(v,start,width),y,color)
        env=row['mapping_own_A50']
        if env['status']=='FINITE':
            line(xpos(env['min'],586,111),y,xpos(env['max'],586,111),y,'#7b5aa0',1.7)
            for v in (env['min'],env['max']):
                line(xpos(v,586,111),y-1.7,xpos(v,586,111),y+1.7,'#7b5aa0',.6)
        else:
            text(641,y,'NA' if env['status']=='INCOMPLETE_NA' else '-', 'muted','middle',env['reason'])
        out.append('</g>')
    text(10,478,'Role D discovery; M main; X method exception; G gene-tree only; O outgroup. Availability F finite; S source NA; Q no anchors; X method NA; - N/A.')
    text(10,488,'Wr: I information-gate NA (finite beta retained); S source NA; W passing working range; - N/A. Chr ranges and mapping envelopes are NOT CI.')
    text(10,498,'OG0000725: unresolved 154 genes / 41 species; not tip absence. Fixed-common scale paired changes: Figure 4. Exact values and identities: supplement.')
    out.append('</svg>')
    return '\n'.join(out)


def write_svg(data, output_path):
    """Create one new file only. Never overwrite v11 or any existing artifact."""
    svg = render_svg(data)
    path = Path(output_path)
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(svg)
    return path
