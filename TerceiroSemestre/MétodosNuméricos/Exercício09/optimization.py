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


    def gradient_method_2D(self, x_0, y_0, lambd, k_max, convergence_criteria, X, Y):
        dx = self.function_d_dx(x_0,y_0,X,Y)
        dy = self.function_d_dy(x_0,y_0,X,Y)
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
            dx = self.function_d_dx(x_0,y_0,X,Y)
            dy = self.function_d_dy(x_0,y_0,X,Y)
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
