""" Graph
    =====

    Um vértice é simplesmente um inteiro `int`:
    Ex. 1, 2, 9, -10, 0

    Uma aresta é uma tupla `tuple` com dois vértices:
    Ex. (1, 5), (3, 9), (0, 0), (12, 42), (-1, 3)
"""
from functools import wraps
from collections import defaultdict
import numpy as np
import _thread as thread

def cache(callback: callable):
    """
    """
    
    @wraps(callback)
    def new_callback(*args, **kwargs):
        if not hasattr(callback, '__cache__'):
            setattr(callback, '__cache__', callback(*args, **kwargs))
        return getattr(callback, '__cache__')

    def trigger(trigger_callback):
        @wraps(trigger_callback)
        def new_callback(*args, **kwargs):
            if hasattr(callback, '__cache__'):
                delattr(callback, '__cache__')
            return trigger_callback(*args, **kwargs)
        return new_callback

    setattr(new_callback, 'trigger', trigger)
    return new_callback

class cache_property(object):
    "Emulate PyProperty_Type() in Objects/descrobject.c"

    def __init__(self, fget=None, fset=None, fdel=None, doc=None):
        self.fget = cache(fget)
        self.fset = fset
        self.fdel = fdel
        if doc is None and fget is not None:
            doc = fget.__doc__
        self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        if self.fget is None:
            raise AttributeError("unreadable attribute")
        return self.fget(obj)

    def __set__(self, obj, value):
        if self.fset is None:
            raise AttributeError("can't set attribute")
        self.fset(obj, value)

    def __delete__(self, obj):
        if self.fdel is None:
            raise AttributeError("can't delete attribute")
        self.fdel(obj)

    def getter(self, fget):
        return type(self)(fget, self.fset, self.fdel, self.__doc__)

    def setter(self, fset):
        return type(self)(self.fget, fset, self.fdel, self.__doc__)

    def deleter(self, fdel):
        return type(self)(self.fget, self.fset, fdel, self.__doc__)

    def trigger(self, func):
        return self.fget.trigger(func)

def threaded(callback: callable):
    @wraps(callback)
    def new_callback(*args, **kwargs):
        return thread.start_new_thread(callback, args, kwargs)
    return new_callback

inf = float('inf')
nan = float('nan')

class BaseGraph(object):
    """
    """
    @cache_property
    def ord_vertex(self):
        return sorted(self._V)

    @cache_property
    def adj_matrix(self):
        """
        """
        return np.array([[((v, w) in self) for w in self] for v in self], dtype=np.int64)
        
    @cache_property
    def adj_struct(self):
        """
        """
        return dict(self._oN)

    def __init__(self, V: set=None, E: set=None):
        ## Basic Graph information
        self._V = set() ## Vertices Set
        self._E = set() ## Edges Set
        
        self._oN = defaultdict(set) ## Vertex Neighbourhood Dictionary (out)
        ## self._iN = defaultdict(set)

        ## distance between vertices
        self._D = defaultdict(lambda: float('inf'))

        if V is not None:
            assert type(V) in {set, list}
            self.add_vertices(*V)

        if E is not None:
            assert type(E) in {set, list}
            self.add_edges(*E)

    def __len__(self):
        return len(self.ord_vertex)

    def __iter__(self):
        return iter(self.ord_vertex)

    def __contains__(self, x: object):
        if self.is_vertex(x):
            return x in self._V
        elif self.is_edge(x):
            return x in self._E
        else:
            return False

    def __getitem__(self, x: object):
        if self.is_vertex(x):
            return self._oN[x]
        elif self.is_edge(x):
            return (x in self._E)
        else:
            raise ValueError()

    def __invert__(self):
        return self._complement()

    def __add__(self, x: object):
        if self.is_vertex(x):
            copy = self.copy
            copy.add_vertex(x)
            return copy
        elif self.is_edge(x):
            copy = self.copy
            copy.add_edge(*x)
            return copy
        elif self.are_vertices(*x):
            copy = self.copy
            copy.add_vertices(*x)
            return copy
        elif self.are_edges(*x):
            copy = self.copy
            copy.add_edges(*x)
            return copy
        else:
            raise TypeError(f'Invalid operand of type `{type(x)}`')

    def __sub__(self, v: object):
        assert type(v) is int
        if type(v) is int:
            copy = self.copy
            copy.rmv_vertex(v)
            return copy
        else:
            raise TypeError('Operand must be of type `int` or `set`.')

    @classmethod
    def is_vertex(cls, x: object):
        return type(x) is int

    @classmethod
    def is_edge(cls, x: object):
        return (type(x) is tuple) and (len(x) == 2) and cls.are_vertices(*x)

    @classmethod
    def are_vertices(cls, *x: tuple):
        return all(cls.is_vertex(y) for y in x)

    @classmethod
    def are_edges(cls, *x: object):
        return all(cls.is_edge(y) for y in x)

    def _add_vertex(self, v: int):
        raise NotImplementedError

    def _rmv_vertex(self, v: int):
        raise NotImplementedError

    def _add_edge(self, v: int, w: int):
        raise NotImplementedError

    def _rmv_edge(self, v: int, w: int):
        raise NotImplementedError

    @ord_vertex.trigger
    @adj_matrix.trigger
    @adj_struct.trigger
    def add_vertex(self, v: int):
        return self._add_vertex(v)

    @ord_vertex.trigger
    @adj_matrix.trigger
    @adj_struct.trigger
    def rmv_vertex(self, v: int):
        return self._rmv_vertex(v)

    @adj_matrix.trigger
    @adj_struct.trigger
    def add_edge(self, v: int, w: int):
        self._add_edge(v, w)

    @adj_matrix.trigger
    @adj_struct.trigger
    def rmv_edge(self, v: int, w: int):
        self._rmv_edge(v, w)
        
    def add_vertices(self, *vertices: tuple):
        for v in vertices: self.add_vertex(v)

    def rmv_vertices(self, *vertices: tuple):
        for v in vertices: self.rmv_vertex(v)

    def add_edges(self, *edges: tuple):
        for v, w in edges: self.add_edge(v, w)

    def rmv_edges(self, *edges: tuple):
        for v, w in edges: self.rmv_edge(v, w)

    def d(self, v: int, w: int):
        """ Distância entre o vértice `v` e o vértice `w`.
        """
        if (v, w) in self._D:
            return self._D[v, w]
        else:
            if w in self._oN[v]:
                self._D[v, w] = 1
                return self._D[v, w]
            else:
                return min(self.d(u, w) for u in self._oN[v]) + 1

    @property
    def copy(self):
        copy = type(self)()
        copy.__dict__.update({k:v.copy for k,v in self.__dict__.items()})
        return copy

    @property
    def E(self):
        return self._E.copy()

    @property
    def V(self):
        return self._V.copy()

    @property
    def m(self):
        return len(self._E)

    @property
    def n(self):
        return len(self._V)

class DiGraph(BaseGraph):
    """ >>> G = DiGraph()
    """

    def __init__(self, vertices: set=None, edges: set=None):

        self._iN = defaultdict(set) ## Vertex Neighbourhood Dictionary (in)

        BaseGraph.__init__(self, vertices=vertices, edges=edges)
        
    def _add_vertex(self, v: int):
        assert self.is_edge(v)
        if v not in self._V:
            self._V.add(v)
        else:
            raise ValueError(f'Vertex `{v}` is already in the graph.')

    def _rmv_vertex(self, v: int):
        assert self.is_edge(v)
        if v in self._V:
            self._V.remove(v)
            while len(self._oN[v]) > 0: ## clears all links in v's neighbourhood
                w = self._oN[v].pop()
                self._iN[w].remove(v)
                self._E.remove((v, w))
            else:
                del self._oN[v]
        else:
            raise ValueError(f'Vertex `{v}` is not in the graph.')

    def _add_edge(self, v: int, w: int):
        assert self.are_edges(v, w)
        if (v, w) not in self._E:
            if v not in self._V: self._V.add(v)
            if w not in self._V: self._V.add(w)
            self._oN[v].add(w)
            self._iN[w].add(v)
            self._E.add((v, w))
            for u in self._oN[w]: self._D[v, u] = min(self._D[w, u] + 1, self._D[v, u])
        else:
            raise ValueError(f'Edge `({v}, {w})` is already in the graph.')

    def _rmv_edge(self, v: int, w: int):
        assert self.are_edges(v, w)
        if (v, w) in self._E:
            self._oN[v].remove(w)
            self._iN[w].remove(v)
            self._E.remove((v, w))
        else:
            raise ValueError(f'Edge `({v}, {w})` is already in the graph.')

    def ideg(self, v: int):
        """ in-degree(v)
        """
        assert self.is_vertex(v)
        return len(self._iN[v])

    def odeg(self, v: int):
        """ out-degree(v)
        """
        assert self.is_vertex(v)
        return len(self._oN[v])

class Graph(BaseGraph):
    
    def __init__(self, *args, **kwargs):
        BaseGraph.__init__(self, *args, **kwargs)

    def _complement(self):
        V = self.V
        E = {(v, w) for v in V for w in V if w > v and not self[v, w]}
        return type(self)(V, E)

    def _add_vertex(self, v: int):
        assert self.is_vertex(v)
        if v not in self._V:
            self._V.add(v)
        else:
            raise ValueError(f'Vertex `{v}` is already in the graph.')

    def _rmv_vertex(self, v: int):
        assert type(v) is int
        if v in self._V:
            self._V.remove(v)
            while len(self._oN[v]) > 0: ## clears all links in v's neighbourhood
                w = self._oN[v].pop()
                self._oN[w].remove(v)
                self._E.remove((w, v))
                self._E.remove((v, w))
            else:
                del self._oN[v]
        else:
            raise ValueError(f'Vertex `{v}` is already in the graph.')

    def _add_edge(self, v: int, w: int):
        assert type(v) is int and type(w) is int
        if (v, w) not in self._E:
            if v not in self._V:
                self._V.add(v)
            if w not in self._V:
                self._V.add(w)
                
            self._oN[v].add(w)
            self._oN[w].add(v)
            
            self._E.add((v, w))
            self._E.add((w, v))
            
            self._D[v, w] = self._D[w, v] = 1
        else:
            raise ValueError(f'Edge `({v}, {w})` is already in the graph.')

    def _rmv_edge(self, v: int, w: int):
        assert type(v) is int and type(w) is int
        if (v, w) in self._E:
            self._oN[v].remove(w)
            self._oN[w].remove(v)
            self._E.remove((v, w))
            self._E.remove((w, v))
        else:
            raise ValueError(f'Edge `({v}, {w})` is already in the graph.')

    def deg(self, v: int):
        assert self.is_vertex(v)
        return len(self._oN[v])

class BaseTree(object):

    def __init__(self, *subtrees):
        assert self.are_trees(*subtrees)
        self.__subtrees = subtrees

    @classmethod
    def is_tree(self, tree: object):
        return isinstance(tree, type(self))

    @classmethod
    def are_trees(cls, *trees):
        return all(cls.is_tree(t) for t in trees)

class Tree(BaseTree):

    def __init__(self, *subtrees):
        BaseTree.__init__(self, *subtrees)

class LaTeX:

    ENDL = r"\\"
    LCUR = r"{"
    RCUR = r"}"

    def __init__(self, env: str, *code: str):
        self.env = env
        self.code = "\n".join(code)

    def __enter__(self, *args, **kwargs):
        ...

    def __exit__(self, *args, **kwargs):
        ...

    @classmethod
    def kwget(cls, key: object, kwargs: dict, default=None):
        try:
            return kwargs[key]
        except KeyError:
            return default

    @classmethod
    def _table(cls, header, rows):
        yield r"\begin{tabular}" + f"{{|{'|'.join(['c' for x in header])}|}}"

        yield r"\hline"

        yield "\t" + r" & ".join([x for x in header]) + r"\\" + "\n" + r"\hline" + "\n"
        
        yield ("\n").join([("\t" + r" & ".join([x for x in row]) + r"\\") for row in rows])

        yield r"\hline"

        yield r"\end{tabular}"


    @classmethod
    def table(cls, header, rows):
        return "\n".join([x for x in cls._table(header, rows)])

    @classmethod
    def vertex(cls, v, x, y, **kwargs):
        color = cls.kwget('color', kwargs, 'black')
        label = cls.kwget('label', kwargs, '')
        fill = f"fill={cls.kwget('fill', kwargs, 'white')}"
        return f"\\node ({v}) [v, {color}, {fill}] at ({x}, {y}) {{{label}}};"

    @classmethod
    def edge(cls, v, w, **kwargs):
        color = cls.kwget('color', kwargs, None)
        label = cls.kwget('label', kwargs, '')
        thick = cls.kwget('thick', kwargs, None)

        in_ = cls.kwget('in', kwargs, None)
        out = cls.kwget('out', kwargs, None)

        
        swap = 'swap' if cls.kwget('swap', kwargs, False) else None
        args = cls.kwget('args', kwargs, ())

        options = ", ".join([x for x in ('edge', color, swap, *args) if x is not None])

        inout = f'[{", ".join([x for x in (in_, out) if x is not None])}]'
        
        return f"\\draw [{options}] ({v}) edge node {{{label}}} ({w});"

    @classmethod
    def _fig(cls, *code, title=None):
        yield r"\begin{fig}" + (f"[{title}]" if title is not None else "")

        yield from code

        yield r"\end{fig}"

    @classmethod
    def fig(cls, *code, title=None):
        return "\n".join(cls._fig(*code, title=title))

    @classmethod
    def _tikz(cls, vmap: dict, emap: dict, *code, **kwargs):
        """
        """

        scale = round(float(cls.kwget('scale', kwargs, 1.0)), 1)

        options = {
            r'>= stealth' : None,
            r'auto' : None,
            r'every node/.style' : f"{cls.LCUR}scale={scale}{cls.RCUR}",
            r'v/.style' : r'{draw, circle}',
            r'edge/.style' : r'{draw, -}',
            **kwargs
        }

        options = ", ".join([((f"{k} = {v}") if (v is not None) else (k)) for (k, v) in options.items()])
        
        yield r"\begin{tikzpicture}" + f"[{options}]"

        yield r"%% Vertices"

        for v in vmap:
            yield f"\t{cls.vertex(v, **vmap[v])}"

        yield ""

        yield r"%% Edges"
        
        for v, w in emap:
            yield f"\t{cls.edge(v, w, **emap[v, w])}"

        yield from code

        yield r"\end{tikzpicture}"

    @classmethod
    def tikz(cls, vmap: dict, emap: dict, *code, **kwargs):
        return "\n".join([line for line in cls._tikz(vmap, emap, *code, **kwargs)])

    @classmethod
    def dump(cls, fname: str, string: str, mode: str='w'):
        with open(fname, mode) as file:
            file.write(string)
