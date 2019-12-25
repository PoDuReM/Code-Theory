from lib.util import choose


def hamming_border(d, n, k):
    m = 2 ** k
    upper = 2 ** n
    bottom = 0
    for i in range((d - 1) / 2):
        bottom = bottom + choose(n, i)
    return m <= upper / bottom


def vg_border(d, n, k):
    qn = 2 ** (n - k)
    right = 0
    for i in range(d - 2):
        right = right + choose(n - 1, i)
    return qn > right
