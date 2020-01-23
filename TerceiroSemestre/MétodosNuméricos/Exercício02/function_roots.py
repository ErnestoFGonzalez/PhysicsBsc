import math
import numpy as np
import pylab



class FunctionRoots():
    def __init__(self, function, function_derivative):
        self.function = function
        self.function_derivative = function_derivative


    def bissection_method(self, initial_interval, convergence_criteria):
        """Finds a root of function with Bissection Method
         with a max error of convergence_criteria """
        a_0 = initial_interval[0]
        b_0 = initial_interval[1]

        a = a_0
        b = b_0

        f_a = self.function(a)
        f_b = self.function(b)

        d = (b-a)/2
        x = a+d

        k = 1

        steps = []
        estimated_x_per_step = []
        error_per_step = []

        while d > convergence_criteria:
            f_x = self.function(x)
            if f_x * f_a >= 0:
                a = x
                f_a = f_x
            else:
                b = x
                f_b = f_x

            steps.append(k)
            estimated_x_per_step.append(x)
            error_per_step.append(abs(d))

            d = (b-a)/2
            x = a + d

            k += 1

        return x, steps, estimated_x_per_step, error_per_step


    def newton_method(self, x_0, x_1, convergence_criteria, k_max):
        #remove x_1 from parameters
        """Finds a root of function using Newton's Method """
        f_0 = self.function(x_0)
        k = 1
        d = f_0/self.function_derivative(x_0)

        steps = []
        estimated_x_per_step = []
        error_per_step = []

        while (abs(d)>convergence_criteria) and k < k_max+1:
            x_0 = x_0 - d
            f_0 = self.function(x_0)
            steps.append(k)
            estimated_x_per_step.append(x_0)
            print("iteração {}; raíz {}".format(k,x_0))
            error_per_step.append(abs(d))

            d = f_0/self.function_derivative(x_0)
            k += 1

        return x_0, steps, estimated_x_per_step, error_per_step


    def secant_method(self, x_0, x_1, convergence_criteria, k_max):
        """Finds a root of function using Secant Method """
        f_0 = self.function(x_0)
        f_1 = self.function(x_1)

        k = 1
        d = f_1*(x_1-x_0)/(f_1-f_0)

        steps = []
        estimated_x_per_step = []
        error_per_step = []

        while (abs(d)> convergence_criteria) and (k < k_max+1):
            x_2 = x_1 - d
            x_0 = x_1
            x_1 = x_2
            f_0 = f_1
            f_1 = self.function(x_1)

            steps.append(k)
            estimated_x_per_step.append(x_2)
            error_per_step.append(abs(d))

            d = f_1*(x_1-x_0)/(f_1-f_0)
            k += 1

        return x_2, steps, estimated_x_per_step, error_per_step


class FindRoot3D():
    def __init__(self, function, gradient, hessian):
        self.function = function
        self.gradient = gradient
        self.hessian = hessian


    # def newthon_method_3d(self, x0, x1, convergence_criteria, k_max):
    #     # Use jacobian matrix vs Use Hessian matrix ?
    #     f_0 = self.function(x_0)
    #     k = 1
    #     d = -self.hessian * self.gradient(x_0)
    #
    #     while (d > convergence_criteria) and (k < k_max+1):

########### APLICATIONS ################
def curva(x):
    return x**2 - 4

def velocidade(x):
    return 2*x

curva_em_estudo = FunctionRoots(function=curva,function_derivative=velocidade)

first_interval_bissection = curva_em_estudo.bissection_method(initial_interval=[0.7, 2.6],
                                        convergence_criteria=10**(-5))
second_interval_bissection = curva_em_estudo.bissection_method(initial_interval=[0.4, 1.7],
                                        convergence_criteria=10**(-5))
third_interval_bissection = curva_em_estudo.bissection_method(initial_interval=[-3, 0.6],
                                        convergence_criteria=10**(-5))

# Plot estimated x value per
# number of iterations for the three different intervals
pylab.figure(figsize=(10,10))
# pylab.grid(color='k', linestyle='--', linewidth=0.5)
pylab.plot(first_interval_bissection[1], first_interval_bissection[2], 'ro', label="x=[0.7, 2.6]")
pylab.plot(second_interval_bissection[1], second_interval_bissection[2], 'ko', label="x=[0.4, 1.7]")
pylab.plot(third_interval_bissection[1], third_interval_bissection[2], 'bo', label="x=[-3, 0.6]")
pylab.legend()
pylab.xlabel("Número de iterações")
pylab.ylabel("Valor estimado")
pylab.show()

################################################################################

first_interval_newton = curva_em_estudo.newton_method(x_0=0.7, x_1=2.6,
                                                      convergence_criteria=10**(-5),
                                                      k_max=40)
print("O método de newton devolve:",first_interval_newton[0])

first_interval_secant = curva_em_estudo.secant_method(x_0=0.7, x_1=2.6,
                                                      convergence_criteria=10**(-5),
                                                      k_max=40)
print("O método da secante devolve:",first_interval_secant[0])
################################################################################

########################################################################################
pylab.figure(figsize=(10,10))
pylab.grid(color='k', linestyle='--', linewidth=0.5)
pylab.plot(first_interval_bissection[1], [math.log(i) for i in first_interval_bissection[3]], 'ro', label = 'Método da Bisseção')
pylab.plot(first_interval_newton[1], [math.log(i) for i in first_interval_newton[3]], 'ko', label = 'Método da Newton')
pylab.plot(first_interval_secant[1], [math.log(i) for i in first_interval_secant[3]], 'bo', label = 'Método da Secante')
pylab.legend()
pylab.xlabel("Número de iterações")
pylab.ylabel("Logaritmo do critério de convergência")
pylab.show()
################################################################################

#########################PARTE2#################################################
def corrente_oscilante(t):
    return 9*(np.exp(-t))*math.sin(2*math.pi*t) - 1.5

# Plots corrente_oscilante from in t=[0,30]
t_values = []
current_value = []
t = -10
while t<30:
    t_values.append(t)
    current_value.append(corrente_oscilante(t))
    t += 0.01
pylab.figure(figsize=(10,10))
pylab.plot(t_values, current_value, 'bo')
pylab.title("Corrente em função do tempo")
pylab.xlabel("Tempo em segundos")
pylab.ylabel("Corrente em miliAmperes")
pylab.show()

def derivada_corrente_oscilante(t):
    return -9*(np.exp(-t))*math.sin(2*math.pi*t) + 18*(np.exp(-t))*math.pi*math.cos(2*math.pi*t)

curva_corrente_oscilante = FunctionRoots(function=corrente_oscilante,
                                         function_derivative=derivada_corrente_oscilante)

# Defines the x0 values to use in Newton's Method for curva_corrente_oscilante
list_x_0s = [0.6, 0.7, 0.75, 0.8, 0.9]
# Newton's Method results for curva_corrente_oscilante
newton_method_results_cco = []
for x_0 in list_x_0s:
    newton_method_results_cco.append([x_0,curva_corrente_oscilante.newton_method(
                                           x_0=x_0,
                                           x_1=0,
                                           convergence_criteria=10**(-6),
                                           k_max=40)])

table_vals = []
col_labels = ["x0", "steps", "x"]
for i in range(len(newton_method_results_cco)):
    table_vals.append([str(newton_method_results_cco[i][0]),
                       str(newton_method_results_cco[i][1][1][-1]),
                       str(newton_method_results_cco[i][1][0])])

fig, ax = pylab.subplots()

# hide axes
fig.patch.set_visible(False)

ax.axis('off')
ax.axis('tight')
ax.table(cellText=table_vals,
         colLabels=col_labels, loc='center')
fig.tight_layout()
pylab.show()
################################################################################

#############################PARTE3#############################################
# y=1.80m; x=90m
def projetil(w):
    return 90*math.tan(w) - (44.145/(math.cos(w)**2)) - 0.8

def projetil_dt(w):
    return 90*(math.sec(w)**2) - 44.145*(math.sec(w)**2)*math.tan(w)


# Plotting of projetil to locate its roots, so we can later pick good initial values
# for the secant method
w_values = []
w_result = []
w = 0
while w<1.4:
    w_values.append(w)
    w_result.append(projetil(w))
    w += 0.0001
pylab.figure(figsize=(10,10))
pylab.plot(w_values, w_result, 'bo')
pylab.title("Equação da trajetória do projétil para x=90 e y=1.8")
pylab.xlabel("Ângulo de lançamento em radianos")
pylab.ylabel("Resultado da equação")
pylab.show()


projetil_x90 = FunctionRoots(function=projetil,
                             function_derivative=projetil_dt)

# Returns a table with initial guess (x0), number of steps and found root
proj_initial_values = [[0.8,0.7],[0.8,0.85],[0.7,0.6],[2,2.3],[4,5]]
proj_table_vals = []
proj_col_labels = ["theta", "theta0", "theta1", "steps", "precision"]
for i_value in proj_initial_values:
    results = projetil_x90.secant_method(x_0=i_value[0],x_1=i_value[1],
                                         convergence_criteria=(10**-6), k_max=40)
    proj_table_vals.append([str(results[0]),
                            str(i_value[0]),
                            str(i_value[1]),
                            str(results[1][-1]),
                            str(10**-6)])

fig, ax = pylab.subplots()

# hide axes
fig.patch.set_visible(False)

ax.axis('off')
ax.axis('tight')
ax.table(cellText=proj_table_vals,
         colLabels=proj_col_labels, loc='center')
fig.tight_layout()
pylab.savefig("parte3AlineaA.png")
pylab.show()

################################################################################

################################################################################
# y=1.80m
# Can I use Newton's method for 3 dimensions?
def projetil_w_x(w,x):
    return x*math.tan(2) - (9.81*(x**2))/(2*(30**2)*(math.cos(w)**2)) - 0.8

def projetil_w_x_gradient(w,x):
    return [x*(math.sec(w)**2)- (9.81*(x**2)*(math.sec(w)**2)*math.tan(w))/(2*(30**2)),
            math.tan(w) - (9.81*x)/((30**2)*(math.cos(w)**2))]
