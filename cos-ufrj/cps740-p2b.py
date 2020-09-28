from graph import Graph, LaTeX
from copy import deepcopy as dc
## Questão 1.

CACHE = {}

FNAME = 'cps740-p2-2.tikz'
SCALE = 1.0

w = 2.0
h = 1.5

_vmap = {
    's': {
    'x': -2*w,
    'y': 0,
    'label': r'$s$',
    },
    'v1': {
    'x':-1*w,
    'y': h,
    'label': r'$v_1$',
    },
    'v2': {
    'x':-1*w,
    'y':-h,
    'label': r'$v_2$',
    },
    'v3': {
    'x': 1*w,
    'y': h,
    'label': r'$v_3$',
    },
    'v4': {
    'x': 1*w,
    'y':-h,
    'label': r'$v_4$',
    },
    't': {
    'x': 2*w,
    'y': 0,
    'label': '$t$',
    }
}

EXTRA = {
    'tx' : 4 * w,
    'ty' : 0,
}

_emap = {
    ('s', 'v1') : {'value': (16, 9), 'label': '$16/9$'},
    ('s', 'v2') : {'value': (13, 5), 'label': '$13/5$', 'swap': True},
    ('v2', 'v1') : {'value': (4, 3), 'label': '$4/3$'},
    ('v1', 'v3') : {'value': (12, 12), 'label': '$12/12$'},
    ('v3', 'v2') : {'value': (9, 7), 'label': '$9/7$', 'swap': True},
    ('v2', 'v4') : {'value': (14, 9), 'label': '$14/9$', 'swap': True},
    ('v4', 'v3') : {'value': (7, 7), 'label': '$7/7$', 'swap': True},
    ('v3', 't') : {'value': (20, 12), 'label': '$20/12$'},
    ('v4', 't') : {'value': (4, 2), 'label': '$4/2$', 'swap': True}
    }

vmap = dc(_vmap)
emap = dc(_emap)

kwargs = {
    'scale' : SCALE,
}

tikz = LaTeX.graph(vmap, emap, **kwargs)

LaTeX.dump(FNAME, tikz)
## -------------- ##
def graphres(_vmap, _emap, _path=None):
    vmap = dc(_vmap)

    ## capacity left
    cemap = dc(_emap)

    for k in list(cemap.keys()):
        x, y = cemap[k]['value']
        if x == y:
            del cemap[k]
        else:
            cemap[k]['value'] = (x - y)
            cemap[k]['label'] = f"${cemap[k]['value']}$"

    ## residual
    remap = dc(_emap)

    ## merge
    for c in list(remap.keys()):
        a, b = c
        k = b, a
        remap[k] = remap[c]
        del remap[c]
        x, y = remap[k]['value']
        if y == 0:
            del remap[k]
        else:
            remap[k]['color'] = 'violet!60'
            remap[k]['value'] = y
            remap[k]['label'] = f"${y}$"
            remap[k]['bend'] = 1
            if 'swap' in remap[k] and remap[k]['swap']:
                remap[k]['bend'] = - remap[k]['bend']
                
    emap = {**cemap, **remap}

    ## path augment
    if _path is not None:
        path = dc(_path)

        for u in path:
            vmap[u]['color'] = 'blue!60'

        bottleneck = float('inf')

        for u, v in zip(path[:-1], path[1:]):
            emap[u, v]['color'] = 'blue!60'
            if emap[u, v]['value'] < bottleneck:
                bottleneck = emap[u, v]['value']


        x = 0
        y = min([vmap[u]['y'] for u in vmap])

        pathtext = ", ".join([vmap[u]['label'][1:-1] for u in path])

        
        text_a = r"\color{blue!60} Caminho aumentante: " + f"$({pathtext})$"
        text_b = r"\color{blue!60} Gargalo: " + f"${bottleneck}$"

        node_a = LaTeX.node(x, y - h/2, text_a)
        node_b = LaTeX.node(x, y - h/1, text_b)
    
        graphtex = LaTeX.graph(vmap, emap, node_a, node_b, **kwargs)

    else:
        graphtex = LaTeX.graph(vmap, emap, **kwargs)

    return LaTeX.fig(graphtex)

    

FNAME = 'cps740-p2-2a.tikz'

_path = ('s', 'v2', 'v3', 't')

fig = graphres(_vmap, _emap, _path)

LaTeX.dump(FNAME, fig)

##
_emap = {
    ('s', 'v1') : {'value': (16, 9), 'label': '$16/9$'},
    ('s', 'v2') : {'value': (13, 12), 'label': '$13/15$', 'swap': True},
    ('v2', 'v1') : {'value': (4, 3), 'label': '$4/3$'},
    ('v1', 'v3') : {'value': (12, 12), 'label': '$12/12$'},
    ('v3', 'v2') : {'value': (9, 0), 'label': '$9/7$', 'swap': True},
    ('v2', 'v4') : {'value': (14, 9), 'label': '$14/9$', 'swap': True},
    ('v4', 'v3') : {'value': (7, 7), 'label': '$7/7$', 'swap': True},
    ('v3', 't') : {'value': (20, 19), 'label': '$20/12$'},
    ('v4', 't') : {'value': (4, 2), 'label': '$4/2$', 'swap': True}
    }

_path = ('s', 'v1', 'v2', 'v4', 't')

fig = graphres(_vmap, _emap, _path)

LaTeX.dump(FNAME, fig, mode='a')

##
_emap = {
    ('s', 'v1') : {'value': (16, 11), 'label': '$16/11$'},
    ('s', 'v2') : {'value': (13, 12), 'label': '$13/12$', 'swap': True},
    ('v2', 'v1') : {'value': (4, 1), 'label': '$4/1$'},
    ('v1', 'v3') : {'value': (12, 12), 'label': '$12/12$'},
    ('v3', 'v2') : {'value': (9, 0), 'label': '$9/0$', 'swap': True},
    ('v2', 'v4') : {'value': (14, 11), 'label': '$14/11$', 'swap': True},
    ('v4', 'v3') : {'value': (7, 7), 'label': '$7/7$', 'swap': True},
    ('v3', 't') : {'value': (20, 19), 'label': '$20/19$'},
    ('v4', 't') : {'value': (4, 4), 'label': '$4/4$', 'swap': True}
    }

_path = None #('s', 'v1', 'v2', 'v4', 't')

fig = graphres(_vmap, _emap, _path)

LaTeX.dump(FNAME, fig, mode='a')

FNAME = 'cps740-p2-2b.tikz'

_path = None #('s', 'v1', 'v2', 'v4', 't')

tikz = LaTeX.graph(_vmap, _emap, **kwargs)

fig = LaTeX.fig(tikz)

LaTeX.dump(FNAME, fig, mode='w')

##
FNAME = 'cps740-p2-2c.tikz'

_emap = {
    ('s', 'v1') : {'value': (16, 8), 'label': '$16/8$'},
    ('s', 'v2') : {'value': (13, 6), 'label': '$13/6$', 'swap': True},
    ('v2', 'v1') : {'value': (4, 4), 'label': '$4/4$'},
    ('v1', 'v3') : {'value': (12, 12), 'label': '$12/12$'},
    ('v3', 'v2') : {'value': (9, 9), 'label': '$9/9$', 'swap': True},
    ('v2', 'v4') : {'value': (14, 11), 'label': '$14/11$', 'swap': True},
    ('v4', 'v3') : {'value': (7, 7), 'label': '$7/7$', 'swap': True},
    ('v3', 't') : {'value': (20, 10), 'label': '$20/10$'},
    ('v4', 't') : {'value': (4, 4), 'label': '$4/4$', 'swap': True}
    }


_path = None #('s', 'v1', 'v2', 'v4', 't')

tikz = LaTeX.graph(_vmap, _emap, **kwargs)

fig = LaTeX.fig(tikz, title="Uma rede de fluxo maximal mas que não é máximo.")

LaTeX.dump(FNAME, fig, mode='w')
