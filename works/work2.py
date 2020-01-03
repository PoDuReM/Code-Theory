from lib.field import Field
from lib.expressions import *
from lib.matrix import get_m_matrix, Matrix

rtl = True
task1_input = 23
task1_errors = 2
task1_is_rc_code = True
task2_input = 23
task2_code = "6660440"
task2_code_is_binary = False
task2_errors = 2
# DOESN'T MATTER
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
    M = mod_length - 1
    print("modulus = %s" % binary_modulus)
    print("GF(2^%s)" % M)
    print("")

    d = task1_errors * 2 + 1

    field = Field(modulus)
    elements = field.get_all_elems()
    bin_elements = list(map(lambda it: tobin(it).zfill(mod_length - 1), elements))
    for i, item in enumerate(bin_elements):
        print("%s: %s" % (str(i).zfill(2), item))

    print("")

    if task1_is_rc_code:
        def gen_brackets():
            for i in range(1, d):
                yield Polynomial([Alpha(i), Num(1)])

        brackets = MultList(list(gen_brackets()))
        print("g(x) = %s" % str(brackets))
        print("")

        steps = brackets.toPolynomialCalc(field)
        for step in steps:
            print(str(step))
        print("")

        print("g(x) = %s" % str(steps[-1]))
        print("")
        n = 2 ** M - 1
        k = n - d + 1
        print("d = %d, n = %d, k = %d" % (d, n, k))
    else:
        print('classes: ')
        classes = get_cyclic_classes(2 ** M - 1)
        for c in classes:
            print("C_%d = %s" % (c[0], c))
        include_classes_nums = {}
        include_classes = []
        for i in range(1, d):
            for c in classes:
                if not (c[0] in include_classes_nums) and i in c:
                    include_classes_nums[c[0]] = True
                    include_classes.append(c)
        print('classes that include numbers 1..%d: %s' % (d - 1, list(include_classes_nums.keys())))
        print('list of M_i:')

        def gen_brackets(c):
            for i in c:
                yield Polynomial([Alpha(i), Num(1)])

        m_array = [MultList(list(gen_brackets(c))) for c in include_classes]
        m_array_reduced = list(map(lambda x: x.toPolynomial(field), m_array))
        for i in range(len(m_array)):
            print("M_%d = %s = %s" % (include_classes[i][0], str(m_array[i]), str(m_array_reduced[i])))
        gx = reduce(lambda x, y: x.mult(y, field), m_array)
        print("g(x) = %s" % str(gx))
        poly = gx.toPolynomial(field)
        print("g(x) = %s" % poly)
        n = 2 ** M - 1
        h = poly.max_power()
        k = n - h
        print("d = %d, n = %d, k = %d" % (d, n, k))


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
    print("GF(2^%s)" % (mod_length - 1))
    print("")

    field = Field(modulus)
    elements = field.get_all_elems()
    bin_elements = map(lambda it: tobin(it).zfill(mod_length - 1), elements)
    for i, item in enumerate(bin_elements):
        print("%s: %s" % (str(i).zfill(2), item))

    print("")

    code_polynomial = binaryToPolynomial(task2_code) if task2_code_is_binary else octToPolynomial(task2_code, field)
    print("g(x) = %s" % str(code_polynomial))

    syndrome_len = task2_errors * 2
    syndrome_values_show = [(i, code_polynomial.calc({"x": field.get(i)}, field)) for i in range(1, syndrome_len + 1)]
    for (i, val) in syndrome_values_show:
        if val == 0:
            print("g(a^%d) = 0" % i)
        else:
            p = field.powerOf(val)
            print("g(a^%d) = a^%d" % (i, p))
    syndrome_values = list(map(lambda tuple: numToAlpha(tuple[1], field), syndrome_values_show))
    syndrome_values_raw = list(map(lambda a: a.calc({}, field), syndrome_values))
    syndrome = Polynomial(syndrome_values)
    print("syndrome(x) = %s" % str(syndrome))
    v = task2_errors + 1
    D = 0
    M = None
    while D == 0:
        v = v - 1
        print(syndrome_values_raw[0:v+1])
        M = get_m_matrix(syndrome_values_raw[0:v+1], field)
        print("with M = ")
        print(M.str_field())
        D = M.determinant_and_ref()
        print("v = %d, D = %d" % (v, D))
    M = get_m_matrix(syndrome_values_raw[0:v + 1], field)
    MSolve = M.append_col(syndrome_values_raw[len(syndrome_values_raw)-v:])
    print(MSolve.str_field())
    lambdas = MSolve.solve()
    print(MSolve.str_field())
    lambdas_coef = list(map(lambda it: numToAlpha(it, field), lambdas))
    lambda_poly = Polynomial(lambdas_coef + [1])
    print(lambda_poly)
    roots = lambda_poly.find_roots(field)
    print("roots: %s" % list(map(lambda it: str(numToAlpha(it, field)), roots)))
    inverse_roots = list(map(lambda it: field.reciprocal(it), roots))
    print("inverse_roots (locators): %s" % list(map(lambda it: str(numToAlpha(it, field)), inverse_roots)))
    y_left = Matrix(len(inverse_roots), len(inverse_roots), field)
    y_left.values = [list(map(lambda it: field.power(it, i), inverse_roots)) for i in range(1, len(inverse_roots)+1)]
    system = y_left.append_col(syndrome_values_raw[:len(inverse_roots)])
    print("values system: ")
    print(system.str_field())
    y_arr = system.solve()
    print("solved: ")
    print(system.str_field())

task1()
task2()
