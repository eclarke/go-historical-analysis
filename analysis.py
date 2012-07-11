def tofloat(x):
    try:
        return float(x)
    except ValueError:
        return None


def mean(seq):
    if seq:
        return sum(seq) / len(seq)
    else:
        return 0


def convergence(resultset):
    rs = []
    for row in resultset:
        rs.append([tofloat(x) for x in row if tofloat(x) and tofloat(x) < 1])
    t = [x[-1] - x[0] for x in rs if x]
    print t
    return mean(t)
