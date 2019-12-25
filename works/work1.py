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

print('minimal distance of straight code:')
print(H.findDistance())

print('minimal distance of dual code:')
print(G.findDistance())


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
print("profile: " + str(list(map(lambda num: "2^%d" % num, span_straight.span_profile()))))

print("")

print('=== generating minimal span form for dual code ===')
span_dual = H.to_minimal_span_form()
print(str(span_dual))
print("profile: " + str(list(map(lambda num: "2^%d" % num, span_straight.span_profile()))))
