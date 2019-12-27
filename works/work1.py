from lib.expressions import Mult, Var, Num, Sum
from lib.field import Field
from lib.matrix import Matrix

# input_matrix_str = """
#     0 0 1 0 1 1 1 1 1
#     0 1 1 0 0 1 0 0 1
#     1 0 0 1 1 1 0 0 1
#     1 0 1 0 1 1 1 1 1
#     0 0 1 1 0 0 1 0 1
# """

# pls no whitespaces in-between bits!
from lib.util import tobin

input_matrix_str = """
    1 1 0 1 1 1 0 0 1 1
    1 0 0 1 1 0 1 1 0 1
    0 0 1 0 1 1 1 0 0 1
    0 0 0 1 1 0 0 1 1 0
"""

input_len_1 = 9
input_dim_1 = 3
input_len_2 = 9
input_dist_2 = 4
rtl = True
input_polynom = 15
input_a = 3
input_b = 5
input_c = 5
input_d = 3
input_x0 = 3

H = Matrix.fromString(input_matrix_str)
print('input (check matrix aka proverochnaya): ')
print(str(H))

# prints both steps and itself
(sys, ii, jj) = H.systematic_form()

trimmed = sys.rtrim()
print('after right-trimming:')
print(str(trimmed))
trans = trimmed.transpose()
print('after transposing:')
print(str(trans))
G = trans.lextend()
G.swap_cols(ii, jj)
print('after left-extending (generator matrix aka porojdayuschaya):')
print(str(G))

print('minimal distance of straight code:')
print(H.find_distance())

print('minimal distance of dual code:')
print(G.find_distance())


print("")

print('=== generating syndrome decoding table ===')
print('after input matrix transposing:')
t = H.transpose()
print(str(t))
table = t.get_syndrome_table()
print("syndrome decoding table:")
print(table.to_str_with_delimiter(t.column_count()))

print("")

print('=== generating minimal span form for straight code ===')
span_straight = G.to_minimal_span_form()
print(str(span_straight))
print("profile: " + str(tuple(map(lambda num: "2^%d" % num, span_straight.span_profile()))))

print("")

print('=== generating minimal span form for dual code ===')
span_dual = H.to_minimal_span_form()
print(str(span_dual))
print("profile: " + str(tuple(map(lambda num: "2^%d" % num, span_straight.span_profile()))))

print("")

print('=== computing field ===')
oct_parsed = int(str(input_polynom), 8)
modulus = oct_parsed if rtl else int(tobin(oct_parsed).zfill(len(str(input_polynom)) * 3)[::-1], 2)
binary_modulus = tobin(modulus)
print("modulus = %s" % binary_modulus)
mod_length = len(binary_modulus)
field = Field(modulus)
elements = field.get_all_elems()
bin_elements = list(map(lambda it: tobin(it).zfill(mod_length - 1), elements))
for i, item in enumerate(bin_elements):
    print("%s: %s" % (str(i).zfill(2), item))
print("")
upper = Sum(Mult(Num(input_a), Var(1)), Num(input_b)).calc({"x": input_x0}, field)
bottom = Sum(Mult(Num(input_c), Var(1)), Num(input_d)).calc({"x": input_x0}, field)
bottom_inverted = field.reciprocal(bottom)
upper_div_bottom = field.multiply(upper, bottom_inverted)
print('f(x) = ( %d * x + %d ) / ( %d * x + %d )' % (input_a, input_b, input_c, input_d))
print('f(%d) = %d / %d = %d * %d = %d' % (input_x0, upper, bottom, upper, bottom_inverted, upper_div_bottom))

