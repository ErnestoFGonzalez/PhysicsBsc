# Author: Ernesto González
# Date: 07/11/2019


import numpy as np
import pylab
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
import math


TINY = 10**(-8)

class LinearSystemOfEquationsSolver():
    def __init__(self, matrix, independent_terms):
        self.matrix = matrix
        self.independent_terms = independent_terms


    def gauss_seidel_method(self, convergence_criteria, k_max):
        """Solves a system of equations using Gauss-Seidel Method without relaxation
           Requires a stop-order based on precision (convergence_criteria)
           and maximum of iterations (k_max)."""
        n = len(self.matrix)
        x = [0,0,0,0,0]
        max_error_per_iter = [0]
        k=1
        max_error = convergence_criteria + TINY

        while (max_error>convergence_criteria) and (k<k_max+1):
            max_error=0
            for i in range(n):
                sum_pre = 0
                sum_pos = 0
                for j in range(i):
                    sum_pre += self.matrix[i][j]*x[j]
                for j in range(i+1,n):
                    sum_pos += self.matrix[i][j]*x[j]
                previous_x_i = x[i]
                x[i] = (self.independent_terms[i]-sum_pre-sum_pos)/self.matrix[i][i]
                error = abs((x[i]-previous_x_i)/x[i])
                if error > max_error:
                    max_error = error
            max_error_per_iter.append(max_error)
            k += 1

        return x, max_error_per_iter


    def gauss_seidel_method_with_relaxation(self, lambd, convergence_criteria, k_max):
        """Solves a system of equations using Gauss-Seidel Method without relaxation
           Requires an initial guess x_0, a stop-order based on precision
           (convergence_criteria) and maximum of iterations (k_max). The relaxation
           factor is lambd"""
        n = len(self.matrix)
        x = [0,0,0,0,0]
        max_error_per_iter = [0]
        k=1
        max_error = convergence_criteria + TINY

        while (max_error>convergence_criteria) and (k<k_max+1):
            max_error=0
            for i in range(n):
                sum_pre = 0
                sum_pos = 0
                for j in range(i):
                    sum_pre += self.matrix[i][j]*x[j]
                for j in range(i+1,n):
                    sum_pos += self.matrix[i][j]*x[j]
                previous_x_i = x[i]
                x[i] = (lambd*(self.independent_terms[i]-sum_pre-sum_pos)/self.matrix[i][i]
                        + (1-lambd)*x[i] )
                error = abs((x[i]-previous_x_i)/x[i])
                if error > max_error:
                    max_error = error
            max_error_per_iter.append(max_error)
            k += 1

        return x, max_error_per_iter


######################### APLICAÇÕES ###########################################

############################ PARTE 1 ###########################################
#### ALÍNEA a.
blocos_e_molas = LinearSystemOfEquationsSolver(
    matrix = [[-5,3,0,0,0],
             [3,-6,3,0,0],
             [0,3,-6,3,0],
             [0,0,3,-6,3],
             [0,0,0,3,-5]],
    independent_terms = [-80,0,0,60,0]
)

bloco_e_molas_gauss_seidel = blocos_e_molas.gauss_seidel_method(
    convergence_criteria=10**(-4),
    k_max=100
)
print("Solução do sistema com método de Gauss-Seidel:", bloco_e_molas_gauss_seidel[0])

bloco_e_molas_gauss_seidel_relax_first_lambd = blocos_e_molas.gauss_seidel_method_with_relaxation(
    lambd=0.5,
    convergence_criteria=10**(-4),
    k_max=200
)
print(r"Solução do sistema com método de Gauss-Seidel com relaxação, λ=0.5:",
    bloco_e_molas_gauss_seidel_relax_first_lambd[0])

bloco_e_molas_gauss_seidel_relax_second_lambd = blocos_e_molas.gauss_seidel_method_with_relaxation(
    lambd=1,
    convergence_criteria=10**(-4),
    k_max=200
)
print(r"Solução do sistema com método de Gauss-Seidel com relaxação, λ=1:",
    bloco_e_molas_gauss_seidel_relax_second_lambd[0])

bloco_e_molas_gauss_seidel_relax_third_lambd = blocos_e_molas.gauss_seidel_method_with_relaxation(
    lambd=1.2,
    convergence_criteria=10**(-4),
    k_max=200
)
print(r"Solução do sistema com método de Gauss-Seidel com relaxação, λ=1.2:",
    bloco_e_molas_gauss_seidel_relax_third_lambd[0])

bloco_e_molas_gauss_seidel_relax_fourth_lambd = blocos_e_molas.gauss_seidel_method_with_relaxation(
    lambd=2,
    convergence_criteria=10**(-4),
    k_max=200
)
print(r"Solução do sistema com método de Gauss-Seidel com relaxação, λ=2:",
    bloco_e_molas_gauss_seidel_relax_fourth_lambd[0])

#### ALÍNEA b.
pylab.figure(figsize=(12,10))
pylab.loglog(range(len(bloco_e_molas_gauss_seidel_relax_first_lambd[1])),
    bloco_e_molas_gauss_seidel_relax_first_lambd[1], 'r', label=r"$\lambda=0.5$")
pylab.loglog(range(len(bloco_e_molas_gauss_seidel_relax_second_lambd[1])),
    bloco_e_molas_gauss_seidel_relax_second_lambd[1], 'k', label=r"$\lambda=1$")
pylab.loglog(range(len(bloco_e_molas_gauss_seidel_relax_third_lambd[1])),
    bloco_e_molas_gauss_seidel_relax_third_lambd[1], 'b', label=r"$\lambda=1.2$")
pylab.loglog(range(len(bloco_e_molas_gauss_seidel_relax_fourth_lambd[1])),
    bloco_e_molas_gauss_seidel_relax_fourth_lambd[1], 'y', label=r"$\lambda=2$")
pylab.legend(loc='lower left', prop={'size': 24})
pylab.xlabel("Número de iterações", fontsize=27)
pylab.ylabel("Precisão", fontsize=27)
pylab.xticks(fontsize=26)
pylab.yticks(fontsize=26)
pylab.savefig("erroPorIteracao.png")
pylab.show()


############################### PARTE 2 ########################################
#### ALÍNEA a.
def first_curve(x):
    return math.sqrt(5-(x**2))

def second_curve(x):
    return (x**2)-1

first_curve_x_values = []
first_curve_y_values = []
first_curve_x_values_neg = []
first_curve_y_values_neg = []
i = -math.sqrt(5)+0.0001
while i<math.sqrt(5):
    first_curve_x_values.append(i)
    first_curve_y_values.append(first_curve(i))
    first_curve_x_values_neg.append(i)
    first_curve_y_values_neg.append(-first_curve(i))
    i += 0.001

second_curve_x_values = []
second_curve_y_values = []
i = -5
while i<5:
    second_curve_x_values.append(i)
    second_curve_y_values.append(second_curve(i))
    i += 0.1

pylab.figure(figsize=(12,10))
pylab.plot(first_curve_x_values, first_curve_y_values, 'y-', label=r"$x^2=5-y^2$")
pylab.plot(second_curve_x_values, second_curve_y_values, 'g-', label=r"$y+1=x^2$")
pylab.plot(first_curve_x_values_neg, first_curve_y_values_neg, 'y-')
pylab.plot([1.56155,-1.60049],[1.56155,1.60049], 'ko')
pylab.legend(loc='upper center',prop={'size': 24})
pylab.xlabel(r"$x$", fontsize=27)
pylab.ylabel(r"$y$", fontsize=27)
pylab.xticks(fontsize=26)
pylab.yticks(fontsize=26)
pylab.savefig("graficosistema2.png")
pylab.show()

### ALÍNEA b.
def f_1(x,y):
    return x**2 + y**2 - 5

def f_2(x,y):
    return -x**2 + y + 1

def df_1_dx(x,y):
    return 2*x

def df_1_dy(x,y):
    return 2*y

def df_2_dx(x,y):
    return -2*x

def df_2_dy(x,y):
    return 1

def function1(x,y):
    return [f_1(x,y),f_2(x,y)]

def jacobian1(x,y):
    return [[df_1_dx(x,y),df_1_dy(x,y)],
            [df_2_dx(x,y),df_2_dy(x,y)]]


def gauss_elimination_method_with_partial_pivot_choice(a,b):
    """Solves system of equations with Gauss Elimination Method
       With Partial Pivot Choice"""
    n = len(a)

    for k in range(0,n-1):
        max = k
        for i in range(k+1,n):
            if abs(a[i][k])>abs(a[max][k]):
                max = i
        if max != k:
            for i in range(k,n):
                temp = a[k][i]
                a[k][i] = a[max][i]
                a[max][i] = temp
            temp = b[k]
            b[k] = b[max]
            b[max] = temp
        for i in range(k+1,n):
            m = a[i][k]/a[k][k]
            for j in range(k,n):
                a[i][j] = a[i][j]-m*a[k][j]
            b[i] = b[i] - m*b[k]

    x = [0 for i in range(n)]
    x[n-1] = b[n-1]/a[n-1][n-1]

    for i in range(n-2,-1,-1):
        sum = 0
        for j in range(i,n):
            sum += a[i][j]*x[j]
        x[i] = (b[i]-sum)/a[i][i]

    return x


class NonLinearSystemOfEquationsSolver():

    def __init__(self, function, jacobian):
        self.function = function
        self.jacobian = jacobian


    def newton_method(self, x_0, convergence_criteria, k_max):
        """Solves system of equations using Newton's Method"""
        n = 2
        x = x_0
        k=1
        max_error = convergence_criteria + TINY

        x_values = [x[0]]
        y_values = [x[1]]

        while (abs(max_error)>convergence_criteria) and (k<k_max+1):
            max_error = 0
            J = self.jacobian(x[0],x[1])
            F = self.function(x[0],x[1])
            d = gauss_elimination_method_with_partial_pivot_choice(J,F)
            for i in range(n):
                x[i] -= d[i]
                if abs(d[i]/x[i])>max_error:
                    max_error = d[i]/x[i]
            x_values.append(x[0])
            y_values.append(x[1])
            k += 1
        J = self.jacobian(x[0],x[1])

        return x, J, x_values, y_values

non_linear_system = NonLinearSystemOfEquationsSolver(
    function=function1,
    jacobian=jacobian1,
)

non_linear_system_newton_first = non_linear_system.newton_method(
    x_0=[3,-3],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_second = non_linear_system.newton_method(
    x_0=[12,-21],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_third = non_linear_system.newton_method(
    x_0=[5,7.32],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_fourth = non_linear_system.newton_method(
    x_0=[0.213,0.12],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_fifth = non_linear_system.newton_method(
    x_0=[4,-1],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_sixth = non_linear_system.newton_method(
    x_0=[-60,-15],
    convergence_criteria=10**(-6),
    k_max=100)
non_linear_system_newton_seventh = non_linear_system.newton_method(
    x_0=[-60,5],
    convergence_criteria=10**(-6),
    k_max=100)

#### ALÍNEA c.
pylab.figure(figsize=(10,10))
pylab.plot(non_linear_system_newton_first[2], non_linear_system_newton_first[3], 'ro-', label=r"$x_0=(3,-3)$")
pylab.plot(non_linear_system_newton_second[2], non_linear_system_newton_second[3], 'bo-', label=r"$x_0=(12,-21)$")
pylab.plot(non_linear_system_newton_third[2], non_linear_system_newton_third[3], 'go-', label=r"$x_0=(5,7.32)$")
pylab.plot(non_linear_system_newton_fourth[2], non_linear_system_newton_fourth[3], 'yo-', label=r"$x_0=(0.213,0.12)$")
pylab.plot(non_linear_system_newton_fifth[2], non_linear_system_newton_fifth[3], 'ko-', label=r"$x_0=(4,-1)$")
# pylab.plot(non_linear_system_newton_sixth[2], non_linear_system_newton_sixth[3], 'co-', label=r"$x_0=(-60,-15)$")
pylab.plot(non_linear_system_newton_seventh[2], non_linear_system_newton_seventh[3], 'mo-', label=r"$x_0=(-60,5)$")
pylab.legend(prop={'size': 24})
pylab.xlabel(r"$x$", fontsize=27)
pylab.ylabel(r"$y$", fontsize=27)
pylab.xticks(fontsize=26)
pylab.yticks(fontsize=26)
pylab.savefig("nonlinearsystemtrajectorie.png")
pylab.show()

print("Solução do sistema (2) pelo método de Newton:")
print(r"\t - para $(x_0,y_0)=(3,-3)$ é ",non_linear_system_newton_first[0])
print(r"\t - para $(x_0,y_0)=(12,-21)$ é ",non_linear_system_newton_second[0])
print(r"\t - para $(x_0,y_0)=(5,7.32)$ é ",non_linear_system_newton_third[0])
print(r"\t - para $(x_0,y_0)=(0.213,0.12)$ é ",non_linear_system_newton_fourth[0])
print(r"\t - para $(x_0,y_0)=(4,-1)$ é ",non_linear_system_newton_fifth[0])
print(r"\t - para $(x_0,y_0)=(-60,-15)$ é ",non_linear_system_newton_sixth[0])
print(r"\t - para $(x_0,y_0)=(-60,5)$ é ",non_linear_system_newton_seventh[0])
########################### PARTE 3 ############################################

class FunctionRoots():
    def __init__(self, function, function_derivative):
        self.function = function
        self.function_derivative = function_derivative


    def newton_method(self, x_0, convergence_criteria, k_max):
        """Finds a root of function using Newton's Method """
        f_0 = self.function(x_0)
        k = 1
        d = f_0/self.function_derivative(x_0)

        steps = []
        estimated_x_per_step = []
        error_per_step = []

        while (abs(d)>convergence_criteria) and (k < k_max+1):
            x_0 = x_0 - d
            f_0 = self.function(x_0)
            steps.append(k)
            estimated_x_per_step.append(x_0)
            error_per_step.append(abs(d))

            d = f_0/self.function_derivative(x_0)
            k += 1

        return x_0, steps, estimated_x_per_step, error_per_step


def func(x):
    return -4*x*np.exp(-2*(x**2)-50)

def func_dx(x):
    return -4*np.exp(-2*(x**2)-50) + 16*(x**2)*np.exp(-2*(x**2)-50)


potencial_min = FunctionRoots(
    function=func,
    function_derivative=func_dx
)

potencial_min_newt = potencial_min.newton_method(
    x_0=0.2,
    convergence_criteria=10**(-6),
    k_max=150,
)
print("Raíz de dU(x)/dx",potencial_min_newt[0])
