from lib.matrix import Matrix


matrix = Matrix.fromString("""
    1 0 0 0 0 0 1 0 1 1
    0 1 0 0 0 0 0 1 1 0
    0 0 1 0 0 0 1 1 0 0
    0 0 0 1 0 0 1 0 1 0
    0 0 0 0 1 0 0 0 1 1
    0 0 0 0 0 1 0 1 0 1
""")
print(matrix)
span = matrix.to_minimal_span_form()
print(span)
print("profile: " + str(list(map(lambda num: "2^%d" % num, span.span_profile()))))