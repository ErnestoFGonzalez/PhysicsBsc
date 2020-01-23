# Author: Ernesto González
# Date: 31/10/2019


import numpy as np
import pylab



class SystemOfEquationsSolver():
    def __init__(self, matrix, independent_terms):
        self.matrix = matrix
        self.independent_terms = independent_terms


    def backsubstitution_method(self):
        """Solves system of equations with Backsubstitution Method"""
        n = len(self.matrix)
        x = [0 for i in range(n)]
        x[n-1] = self.independent_terms[n-1]/self.matrix[n-1][n-1]

        for i in range(n-2,-1,-1):
            sum = 0
            for j in range(i,n):
                sum += self.matrix[i][j]*x[j]
            x[i] = (self.independent_terms[i]-sum)/self.matrix[i][i]

        return x


    def gauss_elimination_method(self):
        """Solves system of equations with Gauss Elimination Method"""
        n = len(self.matrix)

        for k in range(0,n-1):
            for i in range(k+1,n):
                m = self.matrix[i][k]/self.matrix[k][k]
                for j in range(k,n):
                    self.matrix[i][j] -= m*self.matrix[k][j]
                self.independent_terms[i] -= m*self.independent_terms[k]

        x = [0 for i in range(n)]
        x[n-1] = self.independent_terms[n-1]/self.matrix[n-1][n-1]
        for i in range(n-2,-1,-1):
            sum = 0
            for j in range(i,n):
                sum += self.matrix[i][j]*x[j]
            x[i] = (self.independent_terms[i]-sum)/self.matrix[i][i]
        return x


    def gauss_elimination_method_with_partial_pivot_choice(self):
        """Solves system of equations with Gauss Elimination Method
           With Partial Pivot Choice"""
        n = len(self.matrix)

        for k in range(0,n-1):
            max = k
            for i in range(k+1,n):
                if abs(self.matrix[i][k])>abs(self.matrix[max][k]):
                    max = i
            if max != k:
                for i in range(k,n):
                    temp = self.matrix[k][i]
                    self.matrix[k][i] = self.matrix[max][i]
                    self.matrix[max][i] = temp
                temp = self.independent_terms[k]
                self.independent_terms[k] = self.independent_terms[max]
                self.independent_terms[max] = temp
            for i in range(k+1,n):
                m = self.matrix[i][k]/self.matrix[k][k]
                for j in range(k,n):
                    self.matrix[i][j] = self.matrix[i][j]-m*self.matrix[k][j]
                self.independent_terms[i] = self.independent_terms[i] - m*self.independent_terms[k]

        x = [0 for i in range(n)]
        x[n-1] = self.independent_terms[n-1]/self.matrix[n-1][n-1]

        for i in range(n-2,-1,-1):
            sum = 0
            for j in range(i,n):
                sum += self.matrix[i][j]*x[j]
            x[i] = (self.independent_terms[i]-sum)/self.matrix[i][i]

        return x


    def lu_decomposition_method(self):
        """Solves system of equations with LU Decomposition Method"""
        n = len(self.matrix)

        P = [[float(i==j) for i in range(n)] for j in range(n)]
        L = [[float(i==j) for i in range(n)] for j in range(n)]
        U = [[0.0] * n for i in range(n)]

        for i in range(0,n):
            sum = 0
            for k in range(0,i-1):
                sum += U[k][i]*L[i][k]
            max = 0
            k_of_max = 0
            for k in range(i+1,n):
                if abs(self.matrix[k][i]-sum) > max:
                    max = abs(self.matrix[k][i]-sum)
                    k_of_max = k
            if k_of_max != i:
                temp = P[i]
                P[i] = P[k_of_max]
                P[k_of_max] = temp
            for j in range(0,n):
                if j>=i:
                    sum = 0
                    for k in range(0,i):
                        sum += U[k][j]*L[i][k]
                    U[i][j] = self.matrix[i][j] - sum
                else:
                    sum = 0
                    for k in range(0,j):
                        sum += U[k][j]*L[i][k]
                    L[i][j] = (self.matrix[i][j]-sum)/ U[j][j]

        y = [0 for i in range(n)]
        for i in P[0]:
            if i != 0:
                p_0 = P[0].index(i)
                break
        y[0] = self.independent_terms[p_0]/L[0][0]

        for i in range(1,n):
            sum = 0
            for j in range(0,i):
                sum += L[i][j]*y[j]
            for k in P[i]:
                if k != 0:
                    p_i = P[i].index(k)
                    break
            # y[i] = (self.independent_terms[p_i]-sum)/L[i][i]
            y[i] = (self.independent_terms[i]-sum)/L[i][i]

        x = [0 for i in range(n)]
        x[n-1] = y[n-1]/U[n-1][n-1]
        for i in range(n-2,-1,-1):
            sum = 0
            for j in range(i,n):
                sum += U[i][j]*x[j]
            x[i] = (y[i]-sum)/U[i][i]

        return x, L, U



############################# APLICAÇÕES #######################################
print("#################################################################")
print("a.")
############# ALÍNEA a.
# R_3=0
first_matrix = SystemOfEquationsSolver(
                    matrix=[[1,-1,-1],
                            [0,2,-1],
                            [0,0,-1]],
                    independent_terms=[0,-2,-7]
                    )
first_matrix_inverse_substitution = first_matrix.backsubstitution_method()
print(first_matrix_inverse_substitution)

# R_3=2
second_matrix = SystemOfEquationsSolver(
                    matrix=[[1,-1,-1],
                            [0,2,-1],
                            [-2,0,-1]],
                    independent_terms=[0,-2,-7]
                    )
second_matrix_inverse_substitution = second_matrix.backsubstitution_method()
print(first_matrix_inverse_substitution) # resultado incorrecto: assume que a
                                         # matriz já está em forma de escada

print("#################################################################")
print("b.")
######## ALÍNEA b.
first_matrix_gauss_elimination = first_matrix.gauss_elimination_method()
print(first_matrix_gauss_elimination)

second_matrix_gauss_elimination = second_matrix.gauss_elimination_method()
print(second_matrix_gauss_elimination)


print("#################################################################")
print("c.")
######### ALÍNEA c.
#R_3=2; R_2=0
third_matrix = SystemOfEquationsSolver(
                    matrix=[[1,-1,-1],
                            [0,0,-1],
                            [-2,0,-1]],
                    independent_terms=[0,-2,-7]
                    )
# third_matrix_gauss_elimination = third_matrix.gauss_elimination_method()
# print(third_matrix_gauss_elimination) #ZeroDivisionError

third_matrix_gauss_elim_wpc = third_matrix.gauss_elimination_method_with_partial_pivot_choice()
print(third_matrix_gauss_elim_wpc)

print("#################################################################")
print("e.")
print("loading plot...")
###### ALÍNEA e.
v_2 = []
i_1 = []
i_2 = []
i_3 = []

counter = -10
while counter<10.1:
    v_2.append(counter)
    temp_matrix = SystemOfEquationsSolver(
                        matrix=[[1,-1,-1],
                                [0,2,-1],
                                [-2,0,-1]],
                        independent_terms=[0,-counter,-(counter+5)]
                        )
    temp_matrix_gauss_el_wpc = temp_matrix.gauss_elimination_method_with_partial_pivot_choice()
    i_1.append(temp_matrix_gauss_el_wpc[0])
    i_2.append(temp_matrix_gauss_el_wpc[1])
    i_3.append(temp_matrix_gauss_el_wpc[2])
    counter +=0.1

pylab.figure(figsize=(10,10))
pylab.plot(v_2, i_1, 'ro', label=r"$I_1$")
pylab.plot(v_2, i_2, 'ko', label=r"$I_2$")
pylab.plot(v_2, i_3, 'bo', label=r"$I_3$")
pylab.legend()
pylab.xlabel(r"$V_2$")
pylab.ylabel(r"$I$")
pylab.show()


print("#################################################################")
print("f.")
##### ALÍNEA f.
test = SystemOfEquationsSolver(matrix=[[1,-1,-1],
                                       [0,2,-1],
                                       [-2,0,-1]],
                               independent_terms=[0,-2,-7])
test_lu_result = test.lu_decomposition_method()
print("X:",test_lu_result[0])
print("L:",test_lu_result[1])
print("U:",test_lu_result[2])
