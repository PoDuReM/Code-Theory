from lib.field import Field
from lib.matrix import Matrix
from lib.util import *
from lib.expressions import *
from functools import reduce

# input_matrix_str = """
#     0 0 1 0 1 1 1 1 1
#     0 1 1 0 0 1 0 0 1
#     1 0 0 1 1 1 0 0 1
#     1 0 1 0 1 1 1 1 1
#     0 0 1 1 0 0 1 0 1
# """

# pls no whitespaces in-between bits!
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

matrix = Matrix.fromString(input_matrix_str)
print('input (check matrix aka proverochnaya): ')
print(str(matrix))

# prints both steps and itself
sys = matrix.systematic_form()

trimmed = sys.rtrim()
print('after right-trimming:')
print(str(trimmed))
trans = trimmed.transpose()
print('after transposing:')
print(str(trans))
gen = trans.lextend()
print('after left-extending (generator matrix aka porojdayuschaya):')
print(str(gen))

print("")

print('=== generating syndrome decoding table ===')
print('after input matrix transposing:')
t = matrix.transpose()
print(str(t))
