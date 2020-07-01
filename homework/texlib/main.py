def urange(start: int, stop: int=None, step: int=None):
    if step is None:
        if start <= stop:
            step = 1
        else:
            step = -1
    
    if stop is None:
        start = 1
        stop = start

    if start > stop:
        while start >= stop:
            yield start
            start += step
    else:
        while start <= stop:
            yield start
            start += step