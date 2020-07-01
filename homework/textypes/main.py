LCUR = r'\{'
RCUR = r'\}'

## escaping
LESC = r'\left'
RESC = r'\right'

ENDL = r'\\'

def tex(x: object):
    if hasattr(x, '__tex__'):
        return x.__tex__()
    else:
        return x.__str__()


def boldmath(x: object) -> str:
    cmd = r'\mathbf'
    return f'{cmd}{{{tex(x)}}}'

class texlist(list):
    
    LD = f'{LESC}{LCUR}' ## left delimiter
    RD = f'{RESC}{RCUR}' ## right delimiter

    def __tex__(self):
        return f'{self.LD}{",".join([tex(x) for x in self])}{self.RD}'