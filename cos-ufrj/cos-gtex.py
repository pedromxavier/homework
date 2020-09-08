# tex.py

import math

def circle(r, x, y, n, ry=None):
    if ry is None:
        rx = ry = r
    else:
        rx = r
    return [(rx * math.cos(2*math.pi*k/n) + x,
             ry * math.sin(2*math.pi*k/n) + y) for k in range(n)]

OPEN = r"""\begin{fig}\begin{tikzpicture}[
>=stealth,
shorten >=1pt,
auto,
thick,
every node/.style={minimum size=0pt, minimum width = 0em, minimum height = 0em, font={\bfseries}},
vertex/.style={circle,fill=white,draw, minimum width = 3em, minimum height = 3em,
font={\normalsize\bfseries}},
edge/.style={-}
]
"""
EXIT = r"""\end{tikzpicture}\end{fig}"""

def vertex(i, x, y, **kwargs):
    lbl = kwargs['lbl'] if 'lbl' in kwargs else '';

    color = (', fill=' + kwargs['color']) if ('color' in kwargs) else (', fill=white');
    
    print(f"""%% {i}
\\node [vertex{color}] ({i}) [] at ({x:0.2f}, {y:0.2f}) [minimum size=10pt] {{}};""")

def edge(i, j, **kwargs):

    lbl = kwargs['lbl'] if 'lbl' in kwargs else '';
    xs = kwargs['xs'] if 'xs' in kwargs else 0;
    ys = kwargs['ys'] if 'ys' in kwargs else 0;

    swap = (', swap') if (kwargs['swap'] if 'swap' in kwargs else False) else ('')

    color = (', ' + kwargs['color']) if ('color' in kwargs) else (', black');
    
    print(f"""%% {i} -> {j}
\draw [edge{color}{swap}] ({i}) -- node[xshift={xs:0.2f}, yshift={ys:0.2f}] {{}} ({j});""")



def graph(V, ink_r=None,ink_b=None, ink_rr=None,ink_bb=None ):
    if ink_r is None: ink_r = set()
    if ink_b is None: ink_b = set()
    if ink_rr is None: ink_rr = set()
    if ink_bb is None: ink_bb = set()
    
    V = sorted(V)

    RX = 1.5
    RY = 1.5

    c = circle(RX, 0, 0, len(V), RY)

    for k in range(len(V)):
        i, j = ij = V[k]
        if ij == (2,3):
            x1, y1 = c[k+1%len(c)]
            x2, y2 = c[k+2%len(c)]
            x = (x1+x2)/2
            y = (y1+y2)/2
        elif ij == (7,5):
            x1, y1 = c[k-1%len(c)]
            x2, y2 = c[k-2%len(c)]
            x = (x1+x2)/2
            y = (y1+y2)/2
        else:
            x, y = c[k]
            
        if ij in {(2,3), (3,2), (3,5), (5,3), (5,7), (7,5)}:
            x *= 1.6
            y *= 1.6
        else:
            x *= 0.9
            y *= 0.9
            
        if (i, j) in ink_b:
            vertex(f"{i}v{j}", x, y, lbl = ij, color = 'blue!30')
        elif (i, j) in ink_r:
            vertex(f"{i}v{j}", x, y, lbl = ij, color = 'red!30')
        else:
            vertex(f"{i}v{j}", x, y, lbl = ij)
        


    E1 = [((i, j), (j, i), i + j) for (i, j) in V if i < j]
    E2 = [((a, b), (c, d), b + d) for a,b in V for c,d in V if (a == c and b < d)]


    for k in range(len(E1)):
        a, b, c = ((i1, j1), (i2, j2), lbl) = E1[k]
        i = f"{i1}v{j1}"
        j = f"{i2}v{j2}"
        swap = False
        M = {(2,3), (3,2), (3,5), (5,3), (5,7), (7,5), (2, 5), (5, 2), (3, 7), (7,3)}
        if (i1, i1) in M or (i2, j2) in M:
            swap = True
            
        if (a, b) in ink_rr:
            edge(i, j, lbl=lbl, swap=swap, color='red!40')
        elif (a, b) in ink_bb:
            edge(i, j, lbl=lbl, swap=swap, color='blue!40')
        else:
            edge(i, j, lbl=lbl, swap=swap, color='black')

    for k in range(len(E2)):
        a, b, c = ((i1, j1), (i2, j2), lbl) = E2[k]
        i = f"{i1}v{j1}"
        j = f"{i2}v{j2}"
        swap = False
        if (a, b) in ink_rr:
            edge(i, j, lbl=lbl, swap=swap, color='red!40')
        elif (a, b) in ink_bb:
            edge(i, j, lbl=lbl, swap=swap, color='blue!40')
        else:
            edge(i, j, lbl=lbl, swap=swap, color='black')
            
inf = float('inf')

def matrix(X: dict, K: list, x=0, y=0):
    print(r"\node [] (mtx) at (" + str(x) + "," + str(y) + r") {", end="")
    print(r'\begin{tabular}{|c|c|}')
    print(r'\hline')
    print(r'$v$ & $d(u, v)$\\')
    print(r'\hline')
    for i in K:
        x = X[i]
        print(r'$' + str(i) + r'$ & ', end='')
        if x == inf:
            print(r'$\infty$', end='')
        else:
            print(r'$' + str(x) + r'$', end='')
        print(r'\\')
    else:
        print(r'\hline')
    print(r'\end{tabular}', end='')
    print(r'};')

def block(d, V, ink_r, ink_b, ink_rr, ink_bb):
    print(OPEN)

    graph(V, ink_r=ink_r, ink_b=ink_b, ink_rr=ink_rr, ink_bb=ink_bb)
    
    matrix(d, V, 4)

    print(EXIT)
    
    
def main(s = (2, 3)):
    N = [2, 3, 5, 7]
    V = [(m, n) for m in N for n in N if m != n]
    d = {v: (0 if s == v else float('inf')) for v in V}

    E1 = [((i, j), (j, i), i + j) for (i, j) in V if i < j]
    E2 = [((a, b), (c, d), b + d) for a,b in V for c,d in V if (a == c and b < d)]

    E = {(a,b) for (a,b,c) in (E1 + E2)}

    P = {(a,b): c for (a,b,c) in (E1 + E2)}

    Q = set(V)

    p = {} # path

    ink_r = set()
    ink_rr = set()
    while Q:
        ink_b = set()
        ink_bb = set()

        u = min(Q, key=lambda v: d[v])
        Q.remove(u)
        
        ink_r.add(u)
        y = None
        for v in [w for w in Q if (u, w) in E]:
            x = d[u] + P[u, v]
            if x < d[v]:
                d[v] = x
                p[v] = u
            ink_b.add(v)
            ink_bb.add((u, v))
            ink_bb.add((v, u))
        else:
            
            ink_rr = set()
            for w in ink_r:
                subpath = set()
                while w != s and w in p:
                    z = p[w]
                    subpath.add((w, z))
                    subpath.add((z, w))
                    w = z
                else:
                    if w == s:
                        ink_rr.update(subpath)
            block(d, V, ink_r, ink_b, ink_rr, ink_bb)
            
            if Q: print(r'~\\')

if __name__ == '__main__':
    main()
