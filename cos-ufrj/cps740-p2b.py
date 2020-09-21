from graph import Graph, LaTeX

## Questão 1.

CACHE = {}

FNAME = 'cps740-p2-1.tikz'
SCALE = 1.4

w = 2.0
h = 1.5

_vmap = {
    'A': {
    'x': 0.5*w,
    'y': 0,
    'label': 'A',
    },
    'B': {
    'x': w,
    'y':-h,
    'label': 'B',
    },
    'C': {
    'x': 2*w,
    'y':-h,
    'label': 'C',
    },
    'D': {
    'x': 2.5*w,
    'y': 0,
    'label': 'D',
    },
    'E': {
    'x': w,
    'y': h,
    'label': 'E',
    },
    'F': {
    'x': 2*w,
    'y': h,
    'label': 'F',
    }
}

EXTRA = {
    'tx' : 4 * w,
    'ty' : 0,
}

_emap = {
    ('A', 'B') : {'value': 10, 'label': '$10$', 'swap': True},
    ('A', 'E') : {'value': 4, 'label': '$4$'},
    ('B', 'C') : {'value': 1, 'label': '$1$', 'swap': True},
    ('B', 'D') : {'value': 2, 'label': '$2$'},
    ('B', 'E') : {'value': 7, 'label': '$7$', 'swap': True},
    ('C', 'D') : {'value': 9, 'label': '$9$', 'swap': True},
    ('D', 'F') : {'value': 11, 'label': '$11$', 'swap': True},
    ('E', 'F') : {'value': 12, 'label': '$12$'}
    }

vmap = _vmap.copy()
emap = _emap.copy()

## for u,v in list(emap): emap[v,u] = emap[u,v]

kwargs = {
    'scale' : SCALE,
}

tikz = LaTeX.tikz(vmap, emap, **kwargs)

LaTeX.dump(FNAME, tikz)


## 1.1


def graph_table(_vmap: dict, _emap: dict, **extra):
    global CACHE
    vmap = _vmap.copy()
    emap = _emap.copy()
    
    header = [r"$v$", r"$d(u, v)$", r"$\mathbf{r}[v]$"]

    rows = []
    for v in vmap:
        ## rows
        d = extra['d'][v]
        d = d if d != float('inf') else r'\infty'

        r = extra['r'][v]
        r = r if r is not None else r'\square'

        rows.append([str(v), str(d), str(r)])

        ## clear
        if 'color' in vmap[v]: del vmap[v]['color']
        if 'fill' in vmap[v]: del vmap[v]['fill']

    for e in emap:
        if 'color' in emap[e]: del emap[e]['color']
        

    if 'vink' in extra:
        for c, v in extra['vink']:
            vmap[v]['color'] = c
            vmap[v]['fill'] = c

    if 'eink' in extra:
        for c, e in extra['eink']:
            if e in emap:
                emap[e]['color'] = c
            
    kwargs = {
        'scale' : SCALE,
    }

    minitable = LaTeX.table(header, rows)

    tx = extra['tx']
    ty = extra['ty']

    code = (
        f'\\node [] at ({tx}, {ty}) {{{minitable}}};',
    )

    minigraph = LaTeX.tikz(vmap, emap, *code, **kwargs)

    CACHE['vmap'] = vmap
    CACHE['emap'] = emap

    tikz = f"{minigraph}{LaTeX.ENDL}"

    fig = LaTeX.fig(tikz)

    LaTeX.dump(FNAME, fig, 'a')
    

def djkstra(_emap: dict, _vmap: dict, s: object, func: callable):
    global CACHE
    vmap = _vmap.copy()
    emap = _emap.copy()

    ##
    for v in vmap:
        vmap[v]['label'] = ''

    for e in emap:
        emap[e]['label'] = ''
    ##
    
    Q = set()

    d = {}
    r = {}

    for v in vmap:
        d[v] = float('inf')
        r[v] = None
        Q.add(v)

    d[s] = 0
    r[s] = s

    extra = {
            'd': d,
            'r': r,
            'vink' : [('blue!30', s)],
            'eink' : [],
            **EXTRA
        }

    visited = []
    

    func(vmap, emap, **extra)

    while Q:
        u = min(Q, key=lambda w: d[w])
        Q.remove(u)

        visited.append(('red!30', u))

        neighbours = set()

        extra['vink'] = [*visited.copy()]
        extra['eink'] = []

        for v in [w for w in Q if (u, w) in emap]: 
            x = d[u] + emap[u, v]['value']
            if x < d[v]:
                d[v] = x
                r[v] = u
                
            ## neighborhood
            extra['vink'].append(('blue!30', v))
            extra['eink'].append(('blue!40', (v, u)))
            extra['eink'].append(('blue!40', (u, v)))
            neighbours.add((v, u))
            neighbours.add((u, v))
            ##
        else:
            ## spanning tree
            for v in r:
                if v == r[v]: continue
                if (v, r[v]) in neighbours: continue
                if (r[v], v) in neighbours: continue
                extra['eink'].append(('red!40', (v, r[v])))
                extra['eink'].append(('red!40', (r[v], v)))
            
            func(vmap, emap, **extra)

        CACHE['r'] = r

FNAME = 'cps740-p2-1a.tikz'
SCALE = 1.0
LaTeX.dump(FNAME, '') ## clear file

vmap = _vmap.copy()
emap = _emap.copy()

for u,v in list(emap): emap[v,u] = emap[u,v]

djkstra(emap, vmap, 'A', graph_table)

##
FNAME = 'cps740-p2-1b.tikz'
vmap = CACHE['vmap']
emap = CACHE['emap']

r = CACHE['r']

E = set(emap.keys())
A = E - ({(r[v], v) for v in vmap if r[v] != v} | {(v, r[v]) for v in vmap if r[v] != v})

for e in E:
    if e in A: ## to remove
        del emap[e]
    else:
        emap[e]['label'] = ''
        emap[e]['color'] = 'violet!75'

for v in set(vmap.keys()):
    vmap[v]['label'] = v
    vmap[v]['color'] = 'violet!75'
    vmap[v]['fill'] = 'white'

SCALE = 1.4

kwargs = {
    'scale' : SCALE
}

tikz = LaTeX.tikz(vmap, emap, **kwargs)
fig = LaTeX.fig(tikz)
LaTeX.dump(FNAME, fig)

def prim(_emap: dict, _vmap: dict, s: object, func: callable):
    global CACHE
    vmap = _vmap.copy()
    emap = _emap.copy()

    ##
    for v in vmap:
        vmap[v]['label'] = ''

    for e in emap:
        emap[e]['label'] = ''
    ##
        
    P = set()
    Q = set()

    d = {}
    r = {}

    for v in vmap:
        d[v] = float('inf')
        r[v] = None
        Q.add(v)

    d[s] = 0
    r[s] = s

    P.add(s) ## omega
    Q.remove(s) ## resto

    extra = {
            'd': d,
            'r': r,
            'vink' : [('blue!30', s)],
            'eink' : [],
            **EXTRA
        }

    visited = []
    
    func(vmap, emap, **extra)

    visited.append(('red!30', s))

    while Q:
        neighbours = set()

        extra['vink'] = [*visited.copy()]
        extra['eink'] = []

        ## cut
        min_edge = (None, float('inf'))
        for u, v in [(p, q) for p in P for q in Q if (p, q) in emap]: ## cut
            if emap[u, v]['value'] < min_edge[1]: ## new minimum found
                min_edge = ((u, v), emap[u, v]['value'])
                
            ## neighborhood
            extra['vink'].append(('blue!30', v))
            extra['eink'].append(('blue!40', (v, u)))
            extra['eink'].append(('blue!40', (u, v)))
            neighbours.add((v, u))
            neighbours.add((u, v))
            ##
        else:
            uv, duv = min_edge
            if uv is not None:
                u, v = uv
                d[v] = d[u] + duv
                r[v] = u

                P.add(v)
                Q.remove(v)

                visited.append(('red!30', v))
   
            ## spanning tree
            for v in r:
                if v == r[v]: continue
                if (v, r[v]) in neighbours: continue
                if (r[v], v) in neighbours: continue
                extra['eink'].append(('red!40', (v, r[v])))
                extra['eink'].append(('red!40', (r[v], v)))
            
            func(vmap, emap, **extra)
    else:
        neighbours = set()
        extra['eink'] = []
        extra['vink'] = [*visited.copy()]
        ## spanning tree
        for v in r:
            if v == r[v]: continue
            if (v, r[v]) in neighbours: continue
            if (r[v], v) in neighbours: continue
            extra['eink'].append(('red!40', (v, r[v])))
            extra['eink'].append(('red!40', (r[v], v)))
        func(vmap, emap, **extra)

        CACHE['r'] = r
        
FNAME = 'cps740-p2-1d.tikz'
SCALE = 1.0
LaTeX.dump(FNAME, '') ## clear file

vmap = _vmap.copy()
emap = _emap.copy()

for u,v in list(emap): emap[v,u] = emap[u,v]

prim(emap, vmap, 'A', graph_table)

##
FNAME = 'cps740-p2-1e.tikz'
vmap = CACHE['vmap']
emap = CACHE['emap']

r = CACHE['r']

E = set(emap.keys())
A = E - ({(r[v], v) for v in vmap if r[v] != v} | {(v, r[v]) for v in vmap if r[v] != v})

for e in E:
    if e in A: ## to remove
        del emap[e]
    else:
        emap[e]['label'] = ''
        emap[e]['color'] = 'violet!75'

for v in set(vmap.keys()):
    vmap[v]['label'] = v
    vmap[v]['color'] = 'violet!75'
    vmap[v]['fill'] = 'white'

SCALE = 1.4

kwargs = {
    'scale' : SCALE
}

tikz = LaTeX.tikz(vmap, emap, **kwargs)
fig = LaTeX.fig(tikz)
LaTeX.dump(FNAME, fig)
