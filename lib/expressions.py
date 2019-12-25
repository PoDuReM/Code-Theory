from functools import reduce
from lib.util import *

class Expression(object):
    def calc(self, vars, field):
        raise NotImplementedError("calc is not implemented for %s" % str(self))

    def mult(self, other, field):
        raise TypeError("can't multiply %s and %s" % (str(self), str(other)))

    def add(self, other, field):
        raise TypeError("can't multiply %s and %s" % (str(self), str(other)))

    def is_null(self):
        return False


class Num(Expression):
    def __init__(self, a):
        self.num = int(a)

    def calc(self, vars, field):
        return self.num

    def __str__(self):
        return str(self.num)

    def mult(self, other, field):
        if self.num == 1:
            return other
        if self.num == 0:
            return self
        if not isinstance(other, Num):
            raise TypeError("can't multiply %s and %s" % (str(self), str(other)))
        if other.num == 1:
            return self
        if other.num == 0:
            return other
        raise TypeError("can't multiply %s and %s" % (str(self), str(other)))

    def add(self, other, field):
        if self.num == 0:
            return other
        if self.num == 1:
            return Alpha(0).add(other, field)
        raise TypeError("can't add %s and %s" % (str(self), str(other)))

    def is_null(self):
        return self.num == 0


class Var(Expression):
    def __init__(self, power):
        self.power = power

    def calc(self, vars, field):
        if not ("x" in vars):
            raise LookupError("can't find variable x in %s" % vars)
        return field.power(vars["x"], self.power)

    def __str__(self):
        return "x^%d" % self.power


class Alpha(Expression):
    def __init__(self, power):
        self.power = power

    def calc(self, vars, field):
        return field.get(self.power)

    def mult(self, other, field):
        if isinstance(other, Alpha):
            return Alpha((self.power + other.power) % field.maxPower())
        if isinstance(other, Num):
            return other.mult(self, field)
        raise TypeError("can't multiply %s and %s" % (str(self), str(other)))

    def add(self, other, field):
        if isinstance(other, Alpha):
            result = field.add(field.get(self.power), field.get(other.power))
            return Alpha(field.powerOf(result)) if result != 0 else Num(0)
        if isinstance(other, Num):
            return other.add(self, field)
        raise TypeError("can't add %s and %s" % (str(self), str(other)))

    def __str__(self):
        return "a^%d" % self.power


class Binary(Expression):
    def __init__(self, a, b):
        self.left = a if isinstance(a, Expression) else Num(a)
        self.right = b if isinstance(b, Expression) else Num(b)

    def _op(self, field):
        return field.add

    def _op_str(self):
        return "+"

    def calc(self, vars, field):
        return (self._op(field))(self.left.calc(vars, field), self.right.calc(vars, field))

    def __str__(self):
        return "(%s %s %s)" % (str(self.left), self._op_str(), str(self.right))


class Sum(Binary):
    pass


class Mult(Binary):
    def _op(self, field):
        return field.multiply

    def _op_str(self):
        return "*"


class MultList(Expression):
    def __init__(self, *args):
        self.args = map(lambda a: a if isinstance(a, Expression) else Num(a), args[0])

    def calc(self, vars, field):
        results = map(lambda it: it.calc(vars, field), self.args)
        return reduce(lambda x, y: field.multiply(x, y), results)

    def __str__(self):
        return " * ".join(map(str, self.args))

    def toPolynomial(self, field):
        result = []
        cur = self.args[0]
        for i, next in enumerate(self.args[1:]):
            cur = cur.mult(next, field)
            result.append(MultList([cur] + self.args[i+2:]))
        return result


def _polynom_str((i, it), max_power):
    if isinstance(it, Num) and it.num == 1:
        if i == max_power:
            return str(it)
        elif i == max_power - 1:
            return "x"
        else:
            return "x^%d" % (max_power - i)
    if i == max_power:
        return str(it)
    elif i == max_power - 1:
        return "%s * x" % str(it)
    else:
        return "%s * x^%d" % (str(it), max_power - i)

def binaryToPolynomial(s):
    last_index = s.rindex("1")
    return Polynomial([Num(0) if c == '0' else Num(1) for c in s[0:last_index+1]])

class Polynomial(Expression):
    # constant goes first
    def __init__(self, lst):
        self.list = map(lambda a: a if isinstance(a, Expression) else Num(a), lst)

    def calc(self, vars, field):
        results = map(lambda (i, it): Mult(it, Var(i)).calc(vars, field), list(enumerate(self.list)))
        return reduce(lambda x, y: field.add(x, y), results)

    def mult(self, other, field):
        if not isinstance(other, Polynomial):
            raise TypeError("argument is not Polynomial")
        temp = {}
        for i, a in enumerate(self.list):
            for j, b in enumerate(other.list):
                k = i + j
                c = a.mult(b, field)
                if not (k in temp):
                    temp[k] = c
                else:
                    temp[k] = temp[k].add(c, field)
        max_pow = 0
        for (p, val) in temp.items():
            if not val.is_null():
                max_pow = max(max_pow, p)
        result = [Num(0)] * (max_pow + 1)
        for (p, val) in temp.items():
            if not val.is_null():
                result[p] = val
        return Polynomial(result)

    def __str__(self):
        indexed = list(enumerate(self.list[::-1]))

        return "( %s )" % " + ".join(
            map(lambda it: _polynom_str(it, len(self.list) - 1), filter(lambda (i, it): not it.is_null(), indexed)))
