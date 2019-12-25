from lib.matrix import Matrix

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

H = Matrix.fromString(input_matrix_str)
print('input (check matrix aka proverochnaya): ')
print(str(H))

# prints both steps and itself
sys = H.systematic_form()

trimmed = sys.rtrim()
print('after right-trimming:')
print(str(trimmed))
trans = trimmed.transpose()
print('after transposing:')
print(str(trans))
G = trans.lextend()
print('after left-extending (generator matrix aka porojdayuschaya):')
print(str(G))

print('minamal distance of straight code:')
print(H.findDistance())

print('minamal distance of dual code:')
print(G.findDistance())


print("")

print('=== generating syndrome decoding table ===')
print('after input matrix transposing:')
t = H.transpose()
print(str(t))
