from itertools import combinations


def tobin(x):
    return bin(x)[2:]


def all_combinations(lst):
    list_with_tags = list(enumerate(lst))
    for l in range(0, len(list_with_tags) + 1):
        for subset in combinations(list_with_tags, l):
            yield subset


def xor_vectors(a, b):
    return [p1 ^ p2 for (p1, p2) in zip(a, b)]


def is_depended(lists):
    sum = 0
    for i in range(len(lists[0])):
        count = 0
        for j in range(len(lists)):
            count ^= lists[j][i]
        sum += count
    return sum == 0


def cyclic_shift_vector(lst, shift):
    n = len(lst)
    new_list = [None] * n
    for i in range(n):
        new_list[(n + i + shift) % n] = lst[i]
    return new_list


def choose(n, k):
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


def get_cyclic_classes(modulus):
    result = [[0]]
    unused_count = modulus - 1
    used = {0: True}
    for i in range(1, modulus):
        new_class = []
        cur = i
        while not (cur in used or cur in new_class):
            new_class.append(cur)
            cur = (cur * 2) % modulus
        if not (cur in used) and cur in new_class:
            for it in new_class:
                used[it] = True
            unused_count -= len(new_class)
            result.append(new_class)
    return result
