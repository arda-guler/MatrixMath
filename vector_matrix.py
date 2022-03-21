class vector:
    def __init__(self, data=[]):
        self.data = data
        self.dimensions = len(data)

    def __repr__(self):
        return str(self.data)

    def get(self, i=None):
        if not i:
            return self.data
        else:
            return self.data[i]

    def add(self, vec):
        if not self.dimensions == vec.dimensions:
            print("Cannot add vectors of different dimensions!")
            return

        for i in range(self.dimensions):
            self.data[i] += vec.data[i]

    def sub(self, vec):
        if not self.dimensions == vec.dimensions:
            print("Cannot subtract vectors of different dimensions!")
            return

        self.add(vec.scale(-1))

    def scale(self, scalar):
        for i in range(self.dimensions):
            self.data[i] *= scalar

    def get_length(self):
        s_sum = 0
        for element in self.data:
            s_sum += element**2

        return s_sum**(0.5)

    def get_p_norm(self, p):
        setsum = 0
        for element in self.data:
            setsum += abs(element)**p

        return setsum**(1/p)

    def dot(self, vec):
        if not self.dimensions == vec.dimensions:
            print("Cannot compute dot product of vectors of different dimensions!")
            return

        result = 0
        for i in range(self.dimensions):
            result += self.data[i] * vec.data[i]

        return result

    def cross(self, vec):
        if not self.dimensions == 3 or not vec.dimensions == 3:
            print("Cannot compute cross product for non-3D vectors!")
            return

        result = [0,0,0]
        result[0] = self.data[1] * vec.data[2] - self.data[2] * vec.data[1]
        result[1] = - (self.data[0] * vec.data[2] - self.data[2] * vec.data[0])
        result[2] = self.data[0] * vec.data[1] - self.data[1] * vec.data[0]

        return result

class matrix:
    def __init__(self, data=[]):
        self.data = data

        try:
            self.dimensions = [len(data), len(data[0])]
        except:
            self.dimensions = [1, len(data)]

    def __repr__(self):
        if self.dimensions[0] == 1:
            return str(self.data) + "\n" + str(self.dimensions[0]) + "x" + str(self.dimensions[1])
        
        return_str = ""
        for row in self.data:
            return_str += str(row) + "\n"

        return_str += str(self.dimensions[0]) + "x" + str(self.dimensions[1])

        return return_str

    def get(self, i=None, j=None):
        if not i and not j:
            return self.data
        else:
            return self.data[i][j]

    def safe_copy(self):
        copy_data = []
        for i in range(self.dimensions[0]):
            copy_data.append([])
            for j in range(self.dimensions[1]):
                copy_data[i].append(self.data[i][j])

        return matrix(copy_data)

    def add(self, mat):
        if not self.dimensions == mat.dimensions:
            print("Cannot add matrices of different dimensions!")
            return

        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                self.data[i][j] += mat.data[i][j]

    def sub(self, mat):
        if not self.dimensions == mat.dimensions:
            print("Cannot subtract matrices of different dimensions!")
            return

        self.add(mat.scale(-1))

    def scale(self, scalar):
        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                self.data[i][j] *= scalar

        return self

    def mult(self, mat):
        if not self.dimensions[1] == mat.dimensions[0]:
            print("Incompatible matrix dimensions for multiplication!")
            return

        A = self.data
        B = mat.data

        # row, column
        result_dimensions = [self.dimensions[0], mat.dimensions[1]]
        result_data = [[sum(a * b for a, b in zip(A_row, B_col))
                        for B_col in zip(*B)]
                       for A_row in A]

        result = matrix(result_data)
        return result

    def transpose(self):
        transposed_matrix = []
        for i in range(self.dimensions[0]):
            transposed_matrix.append([])
            for j in range(self.dimensions[1]):
                transposed_matrix[i].append(self.data[j][i])

        result = matrix(transposed_matrix)
        return result

    def submatrix(self, si, sj, ni=None, nj=None):
        # this one is used for easily removing a row and column
        if (not ni) and (not nj):
            si -= 1
            sj -= 1

            submatrix_data = []
            for i in range(self.dimensions[0]):
                submatrix_data.append([])
                for j in range(self.dimensions[1]):
                    submatrix_data[i].append(self.data[i][j])

            for row in submatrix_data:
                del row[sj]

            del submatrix_data[si]

            submatrix = matrix(submatrix_data)
            return submatrix
        
        # for everything else in general
        else:
            si -= 1
            sj -= 1
            submatrix_data = []
            rows = 0
            for i in range(si, si + ni):
                submatrix_data.append([])
                for j in range(sj, sj + nj):
                    submatrix_data[rows].append(self.data[i][j])
                rows += 1
            return matrix(submatrix_data)

    def cofactor(self, ci, cj):
        smx = self.submatrix(ci+1, cj+1)
        return (-1)**(ci+cj) * smx.determinant()

    def determinant(self):
        if not self.dimensions[0] == self.dimensions[1]:
            return None

        n = self.dimensions[0]
        if n > 2:
            summation = 0
            for j in range(n):
                summation += self.data[0][j] * self.cofactor(0, j)
        elif n == 2:
            summation = self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0]
        elif n == 1:
            if type(self.data[0]==list):
                summation = self.data[0][0]
            else:
                summation = self.data[0]

        return summation

    def cofactor_matrix(self):
        cofactor_matrix = []
        for i in range(self.dimensions[0]):
            cofactor_matrix.append([])
            for j in range(self.dimensions[1]):
                cofactor_matrix[i].append(0)

        for i in range(self.dimensions[0]):
            for j in range(self.dimensions[1]):
                cofactor_matrix[i][j] = self.cofactor(i, j)

        return matrix(cofactor_matrix)

    def adjugate(self):
        cm = self.cofactor_matrix()
        cmt = cm.transpose()
        return cmt
    
    def inverse(self):
        if not self.dimensions[0] == self.dimensions[1]:
            return None

        if self.determinant() == 0:
            return None

        recidet = 1/self.determinant()
        cmt = self.adjugate()

        return cmt.scale(recidet)

    def rank(self):
        #if self.determinant():
        #    return self.dimensions[0]
        
        N = min(self.dimensions) + 1
        nn = N

        # loop to determine size
        while nn > 0:
            i_max = self.dimensions[0]+1 - nn
            j_max = self.dimensions[1]+1 - nn

            # loop to determine starting point
            for i_start in range(i_max):
                for j_start in range(j_max):
                    cmatrix = []
                    # loop to create matrix with given starting point and size
                    for i in range(nn):
                        cmatrix.append([])
                        for j in range(nn):
                            cmatrix[i].append(self.data[i+i_start][j + j_start])

                    cmatrix = matrix(cmatrix)

                    if not cmatrix.determinant() == 0:
                        return nn

            nn -= 1

        return None

    def augment(self, m2):
        if not self.dimensions[0] == m2.dimensions[0]:
            return

        s = self.safe_copy()
        
        for i in range(self.dimensions[0]):
            for j in range(m2.dimensions[1]):
                s.data[i].append(m2.data[i][j])

        return s
