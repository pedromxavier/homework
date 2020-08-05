from textypes import tex, texlist, boldmath, ENDL
from texlib import urange

L = texlist([2,7,5,6,9,0,1,4,8,5,3])

def sort(L: list, n: int):
    steps = 0
    total = 0
    i = 1
    print(f'L = &{tex(texlist([x for x in L]))}{ENDL}')
    while i < n:
        steps += 1
        if L[i-1] > L[i]:
            L[i-1], L[i] = L[i], L[i-1]
            print(f'&{tex(texlist([L[x] if x not in {i,i-1} else boldmath(L[x]) for x in range(n)]))}{steps}{ENDL}')
            total += steps
            steps = 0
            if i > 1:
                i = i - 1
            else:
                i = i + 1
        else:
            i = i + 1
    total += steps
    print(f'&{tex(texlist([x for x in L]))} \\text{{Total de  ${total}$ passos}}')
    return L

if __name__ == '__main__':
    print(r"\begin{align*}")
    sort(L, len(L))
    print(r"\end{align*}")
