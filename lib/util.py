from expressions import *


def tobin(x):
    return bin(x)[2:]


def numToAlpha(x, field):
    if x == 0:
        return Num(0)
    else:
        p = field.powerOf(x)
        return Alpha(p)
