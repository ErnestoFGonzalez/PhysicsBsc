# Author: Ernesto González
# Date: 24/10/2019

import pylab
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)
import math
import numpy as np

import time



class OptimizationManager1D():
    def __init__(self, function, function_derivative):
        self.function = function
        self.function_derivative = function_derivative


    def gold_number_method(self, initial_interval, convergence_criteria):
        a = initial_interval[0]
        b = initial_interval[1]

        phi = (1+math.sqrt(5))/2
        x_0 = b - (b-a)/phi
        x_1 = a + (b-a)/phi

        k = 2
        x_per_iteration = [x_0,x_1]
        iterations = [0,1]

        while abs(b-a)/(abs(x_0)+abs(x_1))>convergence_criteria:
            if self.function(x_0)<self.function(x_1):
                b = x_1
            else:
                a = x_0
            x_0 = b - (b-a)/phi
            x_1 = a + (b-a)/phi

            x_per_iteration.append((b+a)/2)
            iterations.append(k)
            k += 1

        return (b+a)/2, x_per_iteration, iterations


    def ratio_method(self, ratio, initial_interval, convergence_criteria):
        a = initial_interval[0]
        b = initial_interval[1]

        x_0 = b - (b-a)/ratio
        x_1 = a + (b-a)/ratio

        k = 2
        x_per_iteration = [x_0,x_1]
        iterations = [0,1]

        while abs(b-a)/(abs(x_0)+abs(x_1))>convergence_criteria:
            if self.function(x_0)<self.function(x_1):
                b = x_1
            else:
                a = x_0
            x_0 = b - (b-a)/ratio
            x_1 = a + (b-a)/ratio

            x_per_iteration.append((b+a)/2)
            iterations.append(k)
            k += 1

        return (b+a)/2, x_per_iteration, iterations


    def gradient_method(self, x_0, lambd, k_max, convergence_criteria):
        t_0 = time.time_ns()
        d = self.function_derivative(x_0)
        x_1 = x_0 - lambd*d

        k = 1
        x_per_iteration = [x_0]
        iterations = [0]

        while (abs(lambd*d) > convergence_criteria) and (k < k_max+1):
            iterations.append(k)
            x_per_iteration.append(x_1)
            x_0 = x_1
            d = self.function_derivative(x_0)
            x_1 = x_0 - lambd*d
            k += 1
        t_f = time.time_ns()
        time_for_finding_x = t_f-t_0

        return x_1, x_per_iteration, iterations, time_for_finding_x


    def gradient_method_barzilai_bowen(self, x_0, lambd, k_max, convergence_criteria):
        d = self.function_derivative(x_0)
        x_1 = x_0 - lambd*d
        k = 1
        x_per_iteration = [x_0]
        iterations = [0]
        while (abs(lambd*d) > convergence_criteria) and (k < k_max+1):
            iterations.append(k)
            x_per_iteration.append(x_1)
            x_0 = x_1
            d = self.function_derivative(x_0)
            x_1 = x_0 - lambd*d
            d_1 = self.function_derivative(x_1)
            if d_1 == 0:
                iterations.append(k+1)
                x_per_iteration.append(x_1)
                return x_1, x_per_iteration, iterations
            lambd = ( abs((x_1-x_0)*(d_1-d))
                     / ((d_1-d)**2) )
            k += 1

        return x_1, x_per_iteration, iterations


class OptimizationManager2D():
    def __init__(self, function, function_d_dx, function_d_dy):
        self.function = function
        self.function_d_dx = function_d_dx
        self.function_d_dy = function_d_dy


    def gradient_method_2D(self, x_0, y_0, lambd, k_max, convergence_criteria):
        dx = self.function_d_dx(x_0,y_0)
        dy = self.function_d_dy(x_0,y_0)
        x_1 = x_0 - lambd*dx
        y_1 = y_0 - lambd*dy
        k = 1
        x_list = [x_0]
        y_list = [y_0]
        iterations = [0]

        while (abs(lambd*math.sqrt(dx*dx+dy*dy)) > convergence_criteria) and (k < k_max+1):
            x_list.append(x_1)
            y_list.append(y_1)
            iterations.append(k)
            x_0 = x_1
            y_0 = y_1
            dx = self.function_d_dx(x_0,y_0)
            dy = self.function_d_dy(x_0,y_0)
            x_1 = x_0 - lambd*dx
            y_1 = y_0 - lambd*dy
            k +=1
        return x_1, y_1, x_list, y_list, iterations


class OptimizationManager3D():
    def __init__(self, function, function_d_dx, function_d_dy, function_d_dz):
        self.function = function
        self.function_d_dx = function_d_dx
        self.function_d_dy = function_d_dy
        self.function_d_dz = function_d_dz


    def gradient_method_3D(self, x_0, y_0, z_0, lambd, k_max, convergence_criteria):
        dx = self.function_d_dx(x_0,y_0,z_0)
        dy = self.function_d_dy(x_0,y_0,z_0)
        dz = self.function_d_dz(x_0,y_0,z_0)
        x_1 = x_0 - lambd*dx
        y_1 = y_0 - lambd*dy
        z_1 = z_0 - lambd*dz
        k = 1
        x_list = [x_0]
        y_list = [y_0]
        z_list = [z_0]
        iterations = [0]

        while (abs(lambd*math.sqrt(dx*dx+dy*dy+dz*dz)) > convergence_criteria) and (k < k_max+1):
            x_list.append(x_1)
            y_list.append(y_1)
            z_list.append(z_1)
            iterations.append(k)
            x_0 = x_1
            y_0 = y_1
            z_0 = z_1
            dx = self.function_d_dx(x_0,y_0,z_0)
            dy = self.function_d_dy(x_0,y_0,z_0)
            dz = self.function_d_dz(x_0,y_0,z_0)
            x_1 = x_0 - lambd*dx
            y_1 = y_0 - lambd*dy
            z_1 = z_0 - lambd*dz
            k +=1
        return x_1, y_1, z_1, x_list, y_list, z_list, iterations


class FunctionRoots():
    def __init__(self, function, function_derivative):
        self.function = function
        self.function_derivative = function_derivative


    def newton_method(self, x_0, x_1, convergence_criteria, k_max):
        #remove x_1 from parameters
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

##########################APLICAÇÃO#############################################

############################# PARTE I ##########################################
def mola(x):
    return 0.5*((x-2)**2)

def d_mola_dx(x):
    return x-2


mola_curva = OptimizationManager1D(function=mola,
                                   function_derivative=d_mola_dx)

first_interval_golden = mola_curva.gold_number_method(initial_interval=[-0.7,2.6],
                              convergence_criteria=0.00001)
second_interval_golden = mola_curva.gold_number_method(initial_interval=[0.4,1.7],
                              convergence_criteria=0.00001)
first_interval_ratio_e = mola_curva.ratio_method(ratio=math.e,
                              initial_interval=[-0.7,2.6],
                              convergence_criteria=0.00001)
first_interval_ratio_pi = mola_curva.ratio_method(ratio=math.pi,
                              initial_interval=[-0.7,2.6],
                              convergence_criteria=0.00001)
first_interval_ratio_one_point_five = mola_curva.ratio_method(ratio=1.5,
                              initial_interval=[-0.7,2.6],
                              convergence_criteria=0.00001)
print("Posição de equilíbrio usando método do nº de ouro para intervalo [-0.7,2.6]:",first_interval_golden[0])
print("Posição de equilíbrio usando método do nº de ouro para intervalo [0.4,1.7]:",second_interval_golden[0])
print("Posição de equilíbrio usando método do razão para ratio=e intervalo [-0.7,2.6]:",first_interval_ratio_e[0])
print("Posição de equilíbrio usando método do razão para ratio=π intervalo [-0.7,2.6]:",first_interval_ratio_pi[0])

# plots x per iteration for the various ratio methods (including gold_number_method)
pylab.plot(first_interval_golden[2], first_interval_golden[1], '--go', label=r"$r=\frac{1+\sqrt{5}}{2}$")
pylab.plot(first_interval_ratio_e[2], first_interval_ratio_e[1], '--ro', label=r"$r=e$")
pylab.plot(first_interval_ratio_pi[2], first_interval_ratio_pi[1], '--bo', label=r"$r=\pi$")
pylab.plot(first_interval_ratio_one_point_five[2], first_interval_ratio_one_point_five[1], '--ko', label=r"$r=1.5$")
pylab.legend(prop={'size': 15})
pylab.xlabel("Número de iterações", fontsize=15)
pylab.ylabel("Posição de equilíbrio", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("ratiomethods.png")
pylab.show()

first_lambd_gradient = mola_curva.gradient_method(x_0=0,
                           lambd=0.1,
                           k_max=10,
                           convergence_criteria=10**(-5))
second_lambd_gradient = mola_curva.gradient_method(x_0=0,
                           lambd=0.5,
                           k_max=10,
                           convergence_criteria=10**(-5))
third_lambd_gradient = mola_curva.gradient_method(x_0=0,
                           lambd=1,
                           k_max=10,
                           convergence_criteria=10**(-5))
fourth_lambd_gradient = mola_curva.gradient_method(x_0=0,
                           lambd=2,
                           k_max=10,
                           convergence_criteria=10**(-5))
fifth_lambd_gradient = mola_curva.gradient_method(x_0=0,
                           lambd=2.1,
                           k_max=10,
                           convergence_criteria=10**(-5))
sixth_lambd_gradient = mola_curva.gradient_method_barzilai_bowen(x_0=0,
                           lambd=2.1,
                           k_max=10,
                           convergence_criteria=10**(-5))

print("Posição de equílibrio do corpo encontrada método do gradiente para λ=0.1:",first_lambd_gradient[0])
print("Posição de equílibrio do corpo encontrada método do gradiente para λ=0.5:",second_lambd_gradient[0])
print("Posição de equílibrio do corpo encontrada método do gradiente para λ=1:",third_lambd_gradient[0])
print("Posição de equílibrio do corpo encontrada método do gradiente para λ=2:",fourth_lambd_gradient[0])
print("Posição de equílibrio do corpo encontrada método do gradiente para λ=2.1:",fifth_lambd_gradient[0])
print("Posição de equílibrio do corpo encontrada método do gradiente(Barzilai-Borwein):",sixth_lambd_gradient[0])

pylab.figure(figsize=(10,10))
pylab.plot(first_lambd_gradient[2], first_lambd_gradient[1], '--ro', label="λ=0.1")
pylab.plot(second_lambd_gradient[2], second_lambd_gradient[1], '--ko', label="λ=0.5")
pylab.plot(third_lambd_gradient[2], third_lambd_gradient[1], '--go', label="λ=1")
pylab.plot(fourth_lambd_gradient[2], fourth_lambd_gradient[1], '--yo', label="λ=2")
pylab.plot(fifth_lambd_gradient[2], fifth_lambd_gradient[1], '--bo', label="λ=2.1")
pylab.legend()
pylab.xlabel("Número de iterações")
pylab.ylabel("Posição de equilíbrio")
pylab.show()

pylab.plot(fifth_lambd_gradient[2], fifth_lambd_gradient[1], '--bo', lw=2, label="λ=2.1")
pylab.plot(sixth_lambd_gradient[2], sixth_lambd_gradient[1], '--ro',
    label=r"$λ_0$=2.1, $λ_n=\frac{|(x_n-x_{n-1})(f'(x_n)-f'(x_{n-1}))|}{(f'(x_n)-f'(x_{n-1}))^2}$")
pylab.legend(loc=3, prop={'size': 15})
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.xlabel("Número de iterações", fontsize=15)
pylab.ylabel("Valor mínimo", fontsize=15)
pylab.savefig("otimizacaometodogradiente.png")
pylab.show()

########################### FIM PARTE I ########################################


############################# PARTE II #########################################
def potencial(r, A=80, C=10, B=2):
    return A*np.exp(-B*r)- C/r

def d_potencial_dr(r, A=80, C=10, B=2):
    return -B*A*np.exp(-B*r)+ C/(r**2)

def d2_potencial_dr2(r, A=80, C=10, B=2):
    return (B**2)*A*np.exp(-B*r)- C/(r**3)

r=1
r_list = []
potencial_per_r_list = []
while r < 5:
    r_list.append(r)
    potencial_per_r_list.append(potencial(r))
    r += 0.1


pylab.plot(r_list, potencial_per_r_list, '--ro')
pylab.xlabel("Raio r")
pylab.ylabel("Potencial")
pylab.show()

potencial_minimo = OptimizationManager1D(function=potencial,
                                         function_derivative=d_potencial_dr)
potencial_grad_meth = potencial_minimo.gradient_method(x_0=2,
                                 lambd=0.7,
                                 k_max=20,
                                 convergence_criteria=10**(-5))
print("Distância de equilíbrio encontrada pelo método do gradiente:",potencial_grad_meth[0])

potencial_minimo_newt = FunctionRoots(function=d_potencial_dr,
                                      function_derivative=d2_potencial_dr2)

potencial_newt_meth = potencial_minimo_newt.newton_method(x_0=2,
                            x_1=2.5,
                            convergence_criteria=10**(-5),
                            k_max=20)
print("Distância de equilíbrio encontrada pelo método de Newton:",potencial_newt_meth[0])


# plots equilibrium position per number of iterations
# for both Newton's and gradient methods
pylab.figure(figsize=(10,10))
pylab.plot(potencial_grad_meth[2], potencial_grad_meth[1], '--ro', label="Método do gradiente")
pylab.plot(potencial_newt_meth[1], potencial_newt_meth[2], '--bo', label="Método de Newton")
pylab.legend()
pylab.xlabel("Número de iterações")
pylab.ylabel("Posição de equilíbrio")
pylab.show()


##################### CONSIDER U=U(x,y) ########################################
def potential2D(x, y, A=80, C=10, B=2):
    return A*math.exp(-B*math.sqrt((x**2)+(y**2))) - C/math.sqrt((x**2)+(y**2))

def d_potential2D_dx(x, y, A=80, C=10, B=2):
    return ( ((-B*A*x)/math.sqrt((x**2)+(y**2)))*math.exp(-B*math.sqrt((x**2)+(y**2)))
            + (C*x)/(math.sqrt((x**2)+(y**2))**3) )

def d_potential2D_dy(x, y, A=80, C=10, B=2):
    return ( ((-B*A*y)/math.sqrt((x**2)+(y**2)))*math.exp(-B*math.sqrt((x**2)+(y**2)))
            + (C*y)/(math.sqrt((x**2)+(y**2))**3) )

potencial2D_minimo = OptimizationManager2D(function=potential2D,
                                           function_d_dx=d_potential2D_dx,
                                           function_d_dy=d_potential2D_dy)
potencial_grad_meth_first= potencial2D_minimo.gradient_method_2D(
                                            x_0=5,
                                            y_0=-5,
                                            lambd=0.5,
                                            k_max=100,
                                            convergence_criteria=10**(-7))

# prints (x,y) minimum found by 2D gradient method
print("Mínimo de U(x,y) pelo método do gradiente 2D: (x,y)=({},{})".format(
    potencial_grad_meth_first[0],potencial_grad_meth_first[1]))

# plots y and x per iteration
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(potencial_grad_meth_first[2],potencial_grad_meth_first[3],
           potencial_grad_meth_first[4], c='r', marker='o')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Iterações')
pylab.show()

# plots:
#       - minimum x per iteration
#       - minimum y per iteration
#       - y = y(x)
fig, ax1 = pylab.subplots()
ax1.plot(potencial_grad_meth_first[4], potencial_grad_meth_first[2],
         'o', c='b', lw=1.8, alpha=0.8, label='x por iteração')
ax1.plot(potencial_grad_meth_first[4], potencial_grad_meth_first[3],
         'o', c='r', lw=1.8, alpha=0.8, label='y por iteração')
ax1.set_xlabel("Iterações", fontsize=15)
ax1.set_ylabel("Coordenada", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
ax1.legend(loc=0)
ax2 = pylab.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.15,0.37,0.4,0.4])
ax2.set_axes_locator(ip)
ax2.plot(potencial_grad_meth_first[2], potencial_grad_meth_first[3],
         'x', c='g', mew=2, alpha=0.8,
         label='y=y(x)')
ax2.legend(loc=0)
ax2.set_xlabel("x", fontsize=15)
ax2.set_ylabel("y", fontsize=15)
ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
ax2.tick_params(axis='x', which='major', pad=8)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("potencialGradiente2D.png")
pylab.show()


def potential3D(x, y, z, A=80, C=10, B=2):
    return A*math.exp(-B*math.sqrt((x**2)+(y**2)+(z**2))) - C/math.sqrt((x**2)+(y**2)+(z**2))

def d_potential3D_dx(x, y, z, A=80, C=10, B=2):
    return ( ((-B*A*x)/math.sqrt((x**2)+(y**2)))*math.exp(-B*math.sqrt((x**2)+(y**2)+(z**2)))
            + (C*x)/(math.sqrt((x**2)+(y**2)+(z**2))**3) )

def d_potential3D_dy(x, y, z, A=80, C=10, B=2):
    return ( ((-B*A*y)/math.sqrt((x**2)+(y**2)+(z**2)))*math.exp(-B*math.sqrt((x**2)+(y**2)+(z**2)))
            + (C*y)/(math.sqrt((x**2)+(y**2)+(z**2))**3) )

def d_potential3D_dz(x, y, z, A=80, C=10, B=2):
    return ( ((-B*A*z)/math.sqrt((x**2)+(y**2)+(z**2)))*math.exp(-B*math.sqrt((x**2)+(y**2)+(z**2)))
            + (C*z)/(math.sqrt((x**2)+(y**2)+(z**2))**3) )


potencial3D_minimo = OptimizationManager3D(function=potential3D,
                                           function_d_dx=d_potential3D_dx,
                                           function_d_dy=d_potential3D_dy,
                                           function_d_dz=d_potential3D_dz)
potencial_grad_meth_3D= potencial3D_minimo.gradient_method_3D(
                                            x_0=5,
                                            y_0=-5,
                                            z_0=5,
                                            lambd=0.3,
                                            k_max=100,
                                            convergence_criteria=10**(-7))
print("Mínimo de U(x,y,z) pelo método do gradiente 3D: (x,y,z)=({},{},{})".format(
    potencial_grad_meth_3D[0],potencial_grad_meth_3D[1],potencial_grad_meth_3D[2]))

print("r=",math.sqrt(1.9031056416116172**2 + (-0.783340438394618)**2 + 0.783340438394618**2))
