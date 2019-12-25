class Field:
    def __init__(self, mod):
        if mod <= 1:
            raise ValueError("Invalid modulus")

        # x^5 + x + 1 == 0b100011 or 35
        self.modulus = mod

        self._all_elems = None

        # The number of (unique) elements in this field
        self.size = 1 << (mod.bit_length() - 1)

    def equals(self, x, y):
        return self._check(x) == self._check(y)

    def negate(self, x):
        return self._check(x)

    def add(self, x, y):
        return self._check(x) ^ self._check(y)

    def subtract(self, x, y):
        return self.add(x, y)

    def multiply(self, x, y):
        self._check(x)
        self._check(y)
        result = 0
        while y != 0:
            if y & 1 != 0:
                result ^= x
            x <<= 1
            if x >= self.size:
                x ^= self.modulus
            y >>= 1
        return result

    def power(self, x, a):
        self._check(x)
        if a == 0:
            return 1
        result = x
        while a > 1:
            result = self.multiply(result, x)
            a -= 1
        return result

    def get(self, i):
        if i < 0 or i > self.maxPower():
            raise ValueError("can't get element of field #%d" % i)
        return self.get_all_elems()[i]

    def maxPower(self):
        return len(self.get_all_elems()) - 1

    def powerOf(self, n):
        self._check(n)
        try:
            return self.get_all_elems().index(n)
        except ValueError:
            raise ValueError("field doesn't have element %d: %s" % (n, str(self.get_all_elems())))

    def reciprocal(self, w):
        # Extended Euclidean GCD algorithm
        x = self.modulus
        y = self._check(w)
        if y == 0:
            raise ValueError("Division by zero")
        a = 0
        b = 1
        while y != 0:
            q, r = self._divide_and_remainder(x, y)
            if q == self.modulus:
                q = 0
            x, y = y, r
            a, b = b, (a ^ self.multiply(q, b))
        if x == 1:
            return a
        else:  # All non-zero values must have a reciprocal
            raise AssertionError("Field modulus is not irreducible")

    # Returns a new tuple containing the pair of values (x div y, x mod y).
    def _divide_and_remainder(self, x, y):
        quotient = 0
        ylen = y.bit_length()
        for i in reversed(range(x.bit_length() - ylen + 1)):
            if x.bit_length() == ylen + i:
                x ^= y << i
                quotient |= 1 << i
        return (quotient, x)

    def _check(self, x):
        if not isinstance(x, int):
            raise TypeError()
        if not (0 <= x < self.size):
            raise ValueError("Not an element of this field: " + str(x))
        return x

    def get_all_elems(self):
        if self._all_elems is None:
            if self.modulus > 2:
                cur = 1
                output = []
                while not (cur in output):
                    output.append(cur)
                    cur = self.multiply(cur, 2)
                self._all_elems = output
            else:
                self._all_elems = [1]
        return self._all_elems
