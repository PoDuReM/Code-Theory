from lib.field import Field
from lib.matrix import Matrix
from lib.util import *
from lib.expressions import *
from functools import reduce

rtl = True
task1_input = 57
task1_errors = 4
task1_codename = 'RC'
task2_input = 23
task2_code = "111111110101011"
task2_code_is_binary = True
task2_errors = 2
task2_codename = 'RC'


def task1():
    print("=== TASK 1 ===")
    print("")

    oct_parsed = int(str(task1_input), 8)

    # wheter to interpret 1011 as 1 + x2 + x^3 or 1 + x + x^3
    # else clause: takes parsed number, converts to binary, adds zeroes to the left, reverses, parses as binary
    # (it's the same as reversing number's binary form)
    # e.g. 101001 -> 100101
    modulus = oct_parsed if rtl else int(tobin(oct_parsed).zfill(len(str(task1_input)) * 3)[::-1], 2)

    binary_modulus = tobin(modulus)
    mod_length = len(binary_modulus)
    print("modulus = %s" % binary_modulus)
    # print("GF(2^%s)" % (mod_length - 1))
    print("GF(2^???)")
    print("")

    d = task1_errors * 2 + 1
    n = 2 ** (mod_length - 1) - 1
    k = n - d + 1

    deg = d - 1
    print("floor((d-1)/2) = errors_num => d = (errors_num * 2 + 1) = %d" % d)
    print("n = 2^%d - 1 = %d" % (mod_length - 1, n))
    print("k <= n - d + 1 =(RS-cod udovletvoryaet granitse Singletona)> k = %d" % k)
    print("deg(g(x)) = d - 1 = %d" % deg)
    print("")

    field = Field(modulus)
    elements = field.get_all_elems()
    bin_elements = list(map(lambda it: tobin(it).zfill(mod_length - 1), elements))
    for i, item in enumerate(bin_elements):
        print("%s: %s" % (str(i).zfill(2), item))

    print("")

    def gen_brackets():
        for i in range(1, deg + 1):
            yield Polynomial([Alpha(i), Num(1)])

    brackets = MultList(list(gen_brackets()))
    print("g(x) = %s" % str(brackets))
    # temp = MultList([Polynomial([Alpha(1), Num(1)]), Polynomial([Alpha(2), Num(1)])])
    # print(str(temp))
    # print(str(temp.toPolynomial(field)[0]))
    print("")

    steps = brackets.toPolynomial(field)
    for step in steps:
        print(str(step))
    print("")

    print("g(x) = %s" % str(steps[-1]))
    print("")


def task2():
    print("=== TASK 2 ===")
    print("")

    oct_parsed = int(str(task2_input), 8)

    # wheter to interpret 1011 as 1 + x2 + x^3 or 1 + x + x^3
    # else clause: takes parsed number, converts to binary, adds zeroes to the left, reverses, parses as binary
    # (it's the same as reversing number's binary form)
    # e.g. 101001 -> 100101
    modulus = oct_parsed if rtl else int(tobin(oct_parsed).zfill(len(str(task2_input)) * 3)[::-1], 2)

    binary_modulus = tobin(modulus)
    mod_length = len(binary_modulus)
    print("modulus = %s" % binary_modulus)
    # print("GF(2^%s)" % (mod_length - 1))
    print("GF(2^???)")
    print("")

    field = Field(modulus)
    elements = field.get_all_elems()
    bin_elements = map(lambda it: tobin(it).zfill(mod_length - 1), elements)
    for i, item in enumerate(bin_elements):
        print("%s: %s" % (str(i).zfill(2), item))

    print("")

    code_binary = task2_code if task2_code_is_binary else tobin(int(task2_code, 8)).zfill(len(str(task2_code)) * 3)
    code_polynomial = binaryToPolynomial(code_binary)
    print("g(x) = %s" % str(code_polynomial))

    syndrome_len = task2_errors * 2
    syndrome_values_show = [(i, code_polynomial.calc({"x": field.get(i)}, field))for i in range(1, syndrome_len+1)]
    for (i, val) in syndrome_values_show:
        if val == 0:
            print("g(a^%d) = 0" % i)
        else:
            p = field.powerOf(val)
            print("g(a^%d) = a^%d" % (i, p))
    syndrome_values = map(lambda tuple: numToAlpha(tuple[1], field), syndrome_values_show)
    syndrome = Polynomial(syndrome_values)
    print("syndrome(x) = %s" % str(syndrome))

task1()
task2()
