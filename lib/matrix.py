from functools import reduce
from itertools import combinations
from lib.field import Field
from lib.util import is_depended, tobin, all_combinations, xor_vectors, cyclic_shift_vector

def get_m_matrix(lst, field):
    size = (len(lst) + 1) // 2
    m = Matrix(size, size, field)
    for (i, v) in enumerate(lst):
        x, y = 0, i
        while y >= size:
            y = y - 1
            x = x + 1
        while x < size and y >= 0:
            m.set(y, x, v)
            y = y - 1
            x = x + 1
    return m

class Matrix(object):
    def __init__(self, rows, cols, field=Field(2), fill=None):
        """Constructs a blank matrix with the given number of rows and columns,
        with operations from the given field. All the elements are initially None."""
        if rows <= 0 or cols <= 0:
            raise ValueError("Invalid number of rows or columns")
        if not isinstance(field, Field):
            raise TypeError()

        # The field used to operate on the values in the matrix.
        self.f = field
        # The values of the matrix stored in row-major order, with each element initially None.
        self.values = [[fill] * cols for _ in range(rows)]

    @staticmethod
    def fromString(s, field=Field(2)):
        lines = list(map(lambda line: list(map(lambda n: int(n), line.strip().split())), s.strip().splitlines()))
        output = Matrix(len(lines), len(lines[0]), field)
        for (i, line) in enumerate(lines):
            for (j, v) in enumerate(line):
                output.set(i, j, v)
        return output

    # -- Basic matrix methods --

    def row_count(self):
        """Returns the number of rows in this matrix, which is a positive integer."""
        return len(self.values)

    def column_count(self):
        """Returns the number of columns in this matrix, which is a positive integer."""
        return len(self.values[0])

    def get(self, row, col):
        """Returns the element at the given location in this matrix. The result may be None."""
        if not (0 <= row < len(self.values) and 0 <= col < len(self.values[row])):
            raise IndexError("Row or column index out of bounds")
        return self.values[row][col]

    def set(self, row, col, val):
        """Stores the given element at the given location in this matrix. The value to store can be None."""
        if not (0 <= row < len(self.values) and 0 <= col < len(self.values[row])):
            raise IndexError("Row or column index out of bounds")
        self.values[row][col] = val

    def clone(self):
        """Returns a clone of this matrix. The field and elements are shallow-copied because they are
        assumed to be immutable. Any matrix element can be None when performing this operation."""
        result = Matrix(self.row_count(), self.column_count(), self.f)
        result.values = [list(row) for row in self.values]
        return result

    def transpose(self):
        """Returns a new matrix that is equal to the transpose of this matrix. The field and elements are shallow-copied
        because they are assumed to be immutable. Any matrix element can be None when performing this operation."""
        rows = self.row_count()
        cols = self.column_count()
        result = Matrix(cols, rows, self.f)
        for i in range(rows):
            for j in range(cols):
                result.values[j][i] = self.values[i][j]
        return result

    def __str__(self):
        """Returns a string representation of this matrix. The format is subject to change."""
        result = "[ (%d x %d matrix)\n" % (self.row_count(), self.column_count())
        for (i, row) in enumerate(self.values):
            if i > 0:
                result += "\n"
            result += "\t" + " ".join(str(val) for val in row)
        return result + "\n]"

    def to_str_with_delimiter(self, where):
        result = "[ (%d x %d matrix)\n" % (self.row_count(), self.column_count())
        for (i, row) in enumerate(self.values):
            if i > 0:
                result += "\n"
            result += "\t" + " ".join(str(val) for val in row[0:where]) + " | " + " ".join(
                str(val) for val in row[where:])
        return result + "\n]"

    # -- Simple matrix row operations --

    def swap_rows(self, row0, row1):
        """Swaps the two given rows of this matrix. If the two row indices are the same, the swap is a no-op.
        Any matrix element can be None when performing this operation."""
        if not (0 <= row0 < len(self.values) and 0 <= row1 < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[row0], self.values[row1] = self.values[row1], self.values[row0]

    def multiply_row(self, row, factor):
        """Multiplies the given row in this matrix by the given factor. In other words, row *= factor.
        The elements of the given row should all be non-None when performing this operation."""
        if not (0 <= row < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[row] = [self.f.multiply(val, factor) for val in self.values[row]]

    def add_rows(self, srcrow, destrow, factor):
        """Adds the first given row in this matrix multiplied by the given factor to the second given row.
        In other words, destdow += srcrow * factor. The elements of the given two rows
        should all be non-None when performing this operation."""
        if not (0 <= srcrow < len(self.values) and 0 <= destrow < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[destrow] = [self.f.add(destval, self.f.multiply(srcval, factor))
                                for (srcval, destval) in zip(self.values[srcrow], self.values[destrow])]

    def add_vector(self, destrow, vector):
        if not (len(vector) == len(self.values[0]) and 0 <= destrow < len(self.values)):
            raise IndexError("Row index out of bounds")
        self.values[destrow] = [self.f.add(s1, s2) for (s1, s2) in zip(self.values[destrow], vector)]

    def multiply(self, other):
        """Returns a new matrix representing this matrix multiplied by the given matrix. Requires the given matrix to have
        the same number of rows as this matrix's number of columns. Remember that matrix multiplication is not commutative.
        All elements of both matrices should be non-None when performing this operation.
        The time complexity of this operation is O(self.rows * self.cols * other.cols)."""
        rows = self.row_count()
        cols = other.column_count()
        cells = self.column_count()
        result = Matrix(rows, cols, self.f)
        for i in range(rows):
            for j in range(cols):
                sum = self.f.zero()
                for k in range(cells):
                    sum = self.f.add(self.f.multiply(self.get(i, k), other.get(k, j)), sum)
                result.set(i, j, sum)
        return result

    # -- Advanced matrix operations --

    def reduced_row_echelon_form(self, start=0):
        """Converts this matrix to reduced row echelon form (RREF) using Gauss-Jordan elimination.
        All elements of this matrix should be non-None when performing this operation.
        Always succeeds, as long as the field follows the mathematical rules and does not raise an exception.
        The time complexity of this operation is O(rows * cols * min(rows, cols))."""
        rows = self.row_count()
        cols = self.column_count()

        # Compute row echelon form (REF)
        numpivots = 0
        for j in range(start, cols):  # For each column
            if numpivots >= rows:
                break
            pivotrow = numpivots
            while pivotrow < rows and self.f.equals(self.get(pivotrow, j), 0):
                pivotrow += 1
            if pivotrow == rows:
                continue  # Cannot eliminate on this column
            self.swap_rows(numpivots, pivotrow)
            pivotrow = numpivots
            numpivots += 1

            # Simplify the pivot row
            self.multiply_row(pivotrow, self.f.reciprocal(self.get(pivotrow, j)))

            # Eliminate rows below
            for i in range(pivotrow + 1, rows):
                self.add_rows(pivotrow, i, self.f.negate(self.get(i, j)))

        # Compute reduced row echelon form (RREF)
        for i in reversed(range(numpivots)):
            # Find pivot
            pivotcol = 0
            while pivotcol < cols and self.f.equals(self.get(i, pivotcol), 0):
                pivotcol += 1
            if pivotcol == cols:
                continue  # Skip this all-zero row

            # Eliminate rows above
            for j in range(i):
                self.add_rows(i, j, self.f.negate(self.get(j, pivotcol)))
        return self

    def systematic_form(self):
        print('converting to systematic form...')
        for i in range(self.column_count()):
            for j in range(i, self.column_count()):
                m = self.clone()
                m.swap_cols(i, j)
                m._try_systematic_form()
                if m._is_in_non_reduced_systematic_form():
                    if i != j:
                        print('swapping colomns %d and %d worked: ' % (i + 1, j + 1))
                    print(str(m))
                    print('success! after reduction: ')
                    after = m._reduce_presystematic_form()
                    print(str(after))
                    return after, i, j
                elif i != j:
                    print('swapping columns %d and %d and transforming is not systematic' % (i + 1, j + 1))
        raise ValueError("Failed to transform matrix to systematic form: %s" % str(self))

    def swap_cols(self, col0, col1):
        if col0 == col1:
            return
        for row in self.values:
            row[col0], row[col1] = row[col1], row[col0]

    def _reduce_presystematic_form(self):
        y, x = self.row_count() - 1, self.column_count() - 1
        while x > 0 and y > 0:
            for j in range(y):
                if self.get(j, x) > 0:
                    self.add_rows(y, j, 1)
            x = x - 1
            y = y - 1
        return self

    def _try_systematic_form(self):
        """Tries to convert matrix to systematic_form"""
        rows = self.row_count()
        cols = self.column_count()
        return self.reduced_row_echelon_form(cols - min(rows, cols))

    def rtrim(self):
        rows = self.row_count()
        cols = self.column_count()
        rc = cols - min(rows, cols)
        result = Matrix(self.row_count(), self.column_count(), self.f)
        result.values = [list(row[0:rc]) for row in self.values]
        return result

    def lextend(self):
        result = Matrix(self.row_count(), self.column_count() + self.row_count(), self.f)
        result.values = [[0] * self.row_count() + list(row) for row in self.values]
        x, y = 0, 0
        while x < result.row_count() and y < result.column_count():
            result.values[y][x] = 1
            x = x + 1
            y = y + 1
        return result

    def _is_in_non_reduced_systematic_form(self):
        y, x = self.row_count() - 1, self.column_count() - 1
        while x > 0 and y > 0:
            if self.get(y, x) != 1:
                return False
            x = x - 1
            y = y - 1
        return True

    def invert(self):
        """Replaces the values of this matrix with the inverse of this matrix. Requires the matrix to be square.
        All elements of this matrix should be non-None when performing this operation.
        Raises an exception if the matrix is singular (not invertible). If an exception is raised, this matrix is unchanged.
        The time complexity of this operation is O(rows^3)."""
        rows = self.row_count()
        cols = self.column_count()
        if rows != cols:
            raise RuntimeError("Matrix dimensions are not square")

        # Build augmented matrix: [this | identity]
        temp = Matrix(rows, cols * 2, self.f)
        for i in range(rows):
            for j in range(cols):
                temp.set(i, j, self.get(i, j))
                temp.set(i, j + cols, (self.f.one() if i == j else self.f.zero()))

        # Do the main calculation
        temp.reduced_row_echelon_form()

        # Check that the left half is the identity matrix
        for i in range(rows):
            for j in range(cols):
                if not self.f.equals(temp.get(i, j), (self.f.one() if i == j else self.f.zero())):
                    raise RuntimeError("Matrix is not invertible")

        # Extract inverse matrix from: [identity | inverse]
        for i in range(rows):
            for j in range(cols):
                self.set(i, j, temp.get(i, j + cols))

    def determinant_and_ref(self):
        """Returns the determinant of this matrix, and as a side effect converts the matrix to row echelon form (REF).
        Requires the matrix to be square. The leading coefficient of each row is not guaranteed to be one.
        All elements of this matrix should be non-None when performing this operation.
        Always succeeds, as long as the field follows the mathematical rules and does not raise an exception.
        The time complexity of this operation is O(rows^3)."""
        rows = self.row_count()
        cols = self.column_count()
        if rows != cols:
            raise RuntimeError("Matrix dimensions are not square")
        det = 1

        # Compute row echelon form (REF)
        numpivots = 0
        for j in range(cols):  # For each column
            # Find a pivot row for this column
            pivotrow = numpivots
            while pivotrow < rows and self.f.equals(self.get(pivotrow, j), 0):
                pivotrow += 1

            if pivotrow < rows:
                # This column has a nonzero pivot
                if numpivots != pivotrow:
                    self.swap_rows(numpivots, pivotrow)
                    det = self.f.negate(det)
                pivotrow = numpivots
                numpivots += 1

                # Simplify the pivot row
                temp = self.get(pivotrow, j)
                self.multiply_row(pivotrow, self.f.reciprocal(temp))
                det = self.f.multiply(temp, det)

                # Eliminate rows below
                for i in range(pivotrow + 1, rows):
                    self.add_rows(pivotrow, i, self.f.negate(self.get(i, j)))

            # Update determinant
            det = self.f.multiply(self.get(j, j), det)
        return det

    def find_distance(self):
        for l in range(1, self.column_count() + 1):
            for subset in combinations(self.transpose().values, l):
                # print(isDepended(subset))
                # print(L)
                if (is_depended(subset)):
                    print("solving details: " + str(subset))
                    return l
        return 999

    def previous_set_bit(self, row, index):
        for i in range(index, -1, -1):
            if self.values[row][i] != 0:
                return i
        return -1

    def next_set_bit(self, row, index):
        for i in range(index, self.column_count()):
            if self.values[row][i] != 0:
                return i
        return -1

    def get_syndrome_table(self):
        p = 2 ** self.column_count()
        table = Matrix(p, self.column_count() + self.row_count(), Field(2), 0)
        for i in range(p):
            bin_form = tobin(i).zfill(self.column_count())
            for (j, c) in enumerate(bin_form):
                table.set(i, j, c)
            for subset in all_combinations(self.values):
                sum = reduce(lambda x, y: xor_vectors(x, y[1]), subset, [0] * self.column_count())
                num = int("".join(list(map(str, sum))), 2)
                if num == i:
                    for v in subset:
                        table.set(i, self.column_count() + v[0], 1)
                    break
        return table

    def to_minimal_span_form(self):
        m = self.clone()
        span_heads = [0] * m.row_count()
        span_tails = [0] * m.row_count()
        head_unique_indices = [-1] * m.column_count()
        tail_unique_indices = [-1] * m.column_count()
        for i in range(m.row_count()):
            span_tails[i] = m.previous_set_bit(i, m.column_count() - 1)
            span_heads[i] = m.next_set_bit(i, 0)
            if span_heads[i] == -1:
                raise ValueError("rows are linear-dependent")
            while True:
                cur_head = span_heads[i]
                cur_unique_row = head_unique_indices[cur_head % m.column_count()]
                if cur_unique_row != -1:
                    shift = cur_head - span_heads[cur_unique_row]
                    vector_to_add = cyclic_shift_vector(m.values[cur_unique_row], shift)
                    m.add_vector(i, vector_to_add)
                    if span_tails[cur_unique_row] == span_tails[i]:
                        span_tails[i] = m.previous_set_bit(i, m.column_count() - 1)
                        if span_tails[i] == -1:
                            raise ValueError("rows are linear-dependent")
                    elif span_tails[cur_unique_row] > span_tails[i]:
                        span_tails[i] = span_tails[cur_unique_row]

                    span_heads[i] = m.next_set_bit(i, 0)

                    if span_heads[i] == -1:
                        raise ValueError("rows are linear-dependent")
                else:
                    head_unique_indices[cur_head] = i
                    break
        for i in range(m.row_count()):
            index_to_verify = i
            while True:
                cur_tail = span_tails[index_to_verify]
                cur_unique_row = tail_unique_indices[cur_tail]
                if cur_unique_row != -1:
                    if span_heads[cur_unique_row] > span_heads[index_to_verify]:
                        shift = cur_tail - span_tails[cur_unique_row]
                        shifted_row = cyclic_shift_vector(m.values[cur_unique_row], shift)
                        m.add_vector(index_to_verify, shifted_row)
                        span_tails[index_to_verify] = m.previous_set_bit(index_to_verify, m.column_count() - 1)
                        if span_tails[index_to_verify] == -1:
                            raise ValueError("rows are linear-dependent")
                    else:
                        shift = span_tails[cur_unique_row] - cur_tail
                        shifted_row = cyclic_shift_vector(m.values[index_to_verify], shift)
                        tail_unique_indices[cur_tail] = index_to_verify
                        m.add_vector(cur_unique_row, shifted_row)
                        span_tails[cur_unique_row] = m.previous_set_bit(cur_unique_row, m.column_count() - 1)
                        if span_tails[cur_unique_row] == -1:
                            raise ValueError("rows are linear-dependent")
                        index_to_verify = cur_unique_row
                else:
                    tail_unique_indices[cur_tail] = index_to_verify
                    break
        return m

    def span_profile(self):
        ranges = [(self.next_set_bit(i, 0), self.previous_set_bit(i, self.column_count() - 1))
                  for i in range(self.row_count())]
        result = [0] * self.column_count()
        for i in range(self.column_count()):
            for (start, end) in ranges:
                if start <= i < end:
                    result[i] = result[i] + 1
        return [0] + result + [0]
