import itertools

from lib.expressions import Num, Alpha


def tobin(x):
    return bin(x)[2:]


def numToAlpha(x, field):
    if x == 0:
        return Num(0)
    else:
        p = field.powerOf(x)
        return Alpha(p)


def all_combinations(list):
    for l in range(0, len(list) + 1):
        for subset in itertools.combinations(list, l):
            yield subset


def isDepended(lists):
    numOfOne = 0
    for i in range(len(lists[0])):
        count = 0
        for j in range(len(lists)):
            count ^= lists[j][i]
        numOfOne += count
    return numOfOne == 0
