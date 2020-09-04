# tex.py

import math

def circle(r, x, y, n, ry=None):
    if ry is None:
        rx = ry = r
    else:
        rx = r
    return [(rx * math.cos(2*math.pi*k/n) + x,
             ry * math.sin(2*math.pi*k/n) + y) for k in range(n)]

OPEN = r"""\begin{tikzpicture}[
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
EXIT = r"""\end{tikzpicture}
"""

def vertex(i, x, y, **kwargs):
    lbl = kwargs['lbl'] if 'lbl' in kwargs else '';
    
    print(f"""%% {i}
\\node [vertex] ({i}) [] at ({x:0.2f}, {y:0.2f}) [minimum size=20pt] {{{lbl}}};""")

def edge(i, j, **kwargs):

    lbl = kwargs['lbl'] if 'lbl' in kwargs else '';
    xs = kwargs['xs'] if 'xs' in kwargs else 0;
    ys = kwargs['ys'] if 'ys' in kwargs else 0;

    swap = ',swap' if (kwargs['swap'] if 'swap' in kwargs else False) else ''

    color = ',' + kwargs['color'] if 'color' in kwargs else ',black';
    
    print(f"""%% {i} -> {j}
\draw [edge {color} {swap}] ({i}) -- node[xshift={xs:0.2f}, yshift={ys:0.2f}] {{{lbl}}} ({j});""")


print(OPEN)

N = [2, 3, 5, 7]

V = [(m, n) for m in N for n in N if m != n]
V = sorted(V)

RX = 4.5
RY = 3.5

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
    vertex(f"{i}v{j}", x, y, lbl = ij)


E1 = [((i, j), (j, i)) for (i, j) in V if i < j]
E2 = [((a, b), (c, d)) for a,b in V for c,d in V if (a == c and b < d)]

for k in range(len(E1)):
    ((i1, j1), (i2, j2)) = E1[k]
    i = f"{i1}v{j1}"
    j = f"{i2}v{j2}"
    swap = False
    M = {(2,3), (3,2), (3,5), (5,3), (5,7), (7,5), (2, 5), (5, 2), (3, 7), (7,3)}
    if (i1, i1) in M or (i2, j2) in M:
        swap = True
    edge(i, j, lbl=(i1 + j1), swap=swap, color='blue!40')

for k in range(len(E2)):
    ((i1, j1), (i2, j2)) = E2[k]
    i = f"{i1}v{j1}"
    j = f"{i2}v{j2}"
    swap = False
    edge(i, j, lbl=(i1 + j1), swap=swap, color='red!40')

print(EXIT)
