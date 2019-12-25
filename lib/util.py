from itertools import combinations
from lib.expressions import Num, Alpha


def tobin(x):
    return bin(x)[2:]


def numToAlpha(x, field):
    if x == 0:
        return Num(0)
    else:
        p = field.powerOf(x)
        return Alpha(p)


def all_combinations(lst):
    list_with_tags = list(enumerate(lst))
    for l in range(0, len(list_with_tags) + 1):
        for subset in combinations(list_with_tags, l):
            yield subset


def xor_vectors(a, b):
    return [p1 ^ p2 for (p1, p2) in zip(a, b)]


def isDepended(lists):
    numOfOne = 0
    for i in range(len(lists[0])):
        count = 0
        for j in range(len(lists)):
            count ^= lists[j][i]
        numOfOne += count
    return numOfOne == 0


def cyclic_shift_vector(lst, shift):
    n = len(lst)
    new_list = [None] * n
    for i in range(n):
        new_list[(n + i + shift) % n] = lst[i]
    return new_list
