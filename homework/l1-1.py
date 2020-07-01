from textypes import tex, texlist, boldmath, ENDL
from texlib import urange

L = texlist([2,7,5,6,9,0,1,4,8,5,3])

def sort(L: list, n: int):
    k = 0
    steps = 0
    total = 0
    for i in urange(n-1, 1):
        substeps = 0
        m = float('-inf')
        for j in urange(i, 0):
            substeps += 1
            if L[j] > m:
                m = L[j]
                k = j
        steps += 1
        print(f'&{tex(texlist([L[x] if x not in {i,k} else boldmath(L[x]) for x in range(n)]))}{steps:2d} \\times {substeps:2d} = {steps*substeps:02d}{ENDL}')
        L[i], L[k] = L[k], L[i]
        total += steps * substeps
        
    print(f'&{tex(texlist([x for x in L]))} total = {total}')
    return L

if __name__ == '__main__':
    print(r"\begin{align*}")
    sort(L, len(L))
    print(r"\end{align*}")