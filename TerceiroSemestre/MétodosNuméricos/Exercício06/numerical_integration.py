# Author: Ernesto González
# Date: 14/11/2019


import numpy as np
import pylab
from scipy.interpolate import interp1d


TINY = 10**(-7)

class IntegrationManager():
    def __init__(self, function):
        self.function = function

    def trapezoidal_rule(self, a, b, n):

        h = (b-a)/n
        soma_h = self.function(a)+self.function(b)
        for i in range(n):
            soma_h = soma_h+2*self.function(a+i*h)
        I = (h/2)*soma_h

        return I


    def simpson_rule(self, a, b, n):

        h = (b-a)/n
        S_h = self.function(a)+self.function(b)
        for i in range(1,n+1,2):
            S_h += 4*self.function(a+i*h)
        for i in range(2,n,2):
            S_h += 2*self.function(a+i*h)
        S_h *= (h/3)

        S_2h = self.function(a)+self.function(b)
        for i in range(1,int((n+1)/2),2):
            S_2h += 4*self.function(a+i*2*h)
        for i in range(2,int(n/2),2):
            S_2h += 2*self.function(a+i*2*h)
        S_2h *= (2*h/3)
        erro = abs(S_h-S_2h)/15
        return S_h, erro


    def romberg_method(self, a, b, k_max, convergence_criteria):

        I = [[0 for i in range(k_max+1)] for j in range(k_max+1)]

        k = 0
        n = 0
        I[0][0] = self.trapezoidal_rule(a=a,b=b,n=(2**n))

        error_per_iter = []

        error = convergence_criteria + TINY

        while (error > convergence_criteria) and (k<k_max):
            n += 1
            I[n][0] = self.trapezoidal_rule(a=a,b=b,n=(2**n))
            k += 1
            for j in range(1,k+1):
                # i = k
                for i in range(1,n+1):
                    I[i][j] = ((2**(2*k))*I[i][j-1]-I[i-1][j-1])/(2**(2*k)-1)
            error = abs(I[k][k]-I[k][k-1]/I[k][k])
            error_per_iter.append(error)

        return I, error_per_iter


###################################### APLICAÇÕES ##############################

#################################### PARTE 1 ###################################
print("PARTE 1")
# # teste a.
# def test_func(x):
#     return x**3 + x**2
#
# integral = IntegrationManager(
#     function=test_func,
# )
# print("Resultado Real:",((5**4)/4)+((5**3)/3)-((2**4)/4)-((2**3)/3))
# resultado_trapezio = integral.trapezoidal_rule(a=2,b=5,n=100)
# print("Resultado Trapezio:",resultado_trapezio)
#
# resultado_simpson = integral.simpson_rule(a=2,b=5,n=100)
# print("Resultado Simpson:",resultado_simpson[0])
#
# resultado_romberg = integral.romberg_method(a=2,b=5,k_max=25,convergence_criteria=10**(-6))
# print("Resultado Romberg:",resultado_romberg[25][25])

#### ALÍNEA a.
print("ALÍENA a.")
def volume_imerso(h):
    """Volume imerso de uma esfera de raio r em que h é a altura da secção
    à superfície"""
    return (4*np.pi*(2**3))/3 - (np.pi/3)*(h**2)*(3*2-h)

def impulsao(h,densidade=1000):
    return densidade*volume_imerso(h)*9.81

v = 0
v_lista = []
impulsao_lista = []
while v<4.1:
    v_lista.append(v)
    impulsao_lista.append(impulsao(v))
    v +=0.1

pylab.figure(figsize=(10,10))
pylab.plot(v_lista,impulsao_lista, 'ro')
pylab.xlabel(r"$v$", fontsize=16)
pylab.ylabel(r"$I$", fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.show()



trabalho_impulsao = IntegrationManager(
    function=impulsao,
)
# Trabalho da impulsão calculado analiticamente: 123276.0957
trabalho_impulsao_trapezio = trabalho_impulsao.trapezoidal_rule(a=2,b=4,n=100)
print("Trabalho da impulsão, pela Regra do Trapezio:",trabalho_impulsao_trapezio)

trabalho_impulsao_simpson = trabalho_impulsao.simpson_rule(a=2,b=4,n=100)
print("Trabalho da impulsão, pela Regra de Simpson:",trabalho_impulsao_simpson[0])

trabalho_impulsao_romberg = trabalho_impulsao.romberg_method(a=2,b=4,k_max=25,convergence_criteria=10**(-6))
print("Trabalho da impulsão, pelo Método de Romberg:",trabalho_impulsao_romberg[0][25][25])

#### ALÍNEA b.
print("ALÍENA b.")
print("... a carregar gráfico")
valor_real = 123276.0957
numero_divisoes = [i for i in range(1,101)]
desvio_regra_trapezio = []
desvio_regra_simpson = []
for n in range(1,101):
    desvio_regra_trapezio.append(trabalho_impulsao.trapezoidal_rule(a=2,b=4,n=n)-valor_real)
    desvio_regra_simpson.append(trabalho_impulsao.simpson_rule(a=2,b=4,n=n)[0]-valor_real)

pylab.plot(numero_divisoes, desvio_regra_trapezio, 'bo', label="Regra do Trapézio")
pylab.plot(numero_divisoes, desvio_regra_simpson, 'ro', label="Regra de Simpson")
pylab.legend(prop={'size': 30})
pylab.xlabel("$n$", fontsize=30)
pylab.ylabel("$d$", fontsize=30)
pylab.xticks(fontsize=30)
pylab.yticks(fontsize=30)
pylab.savefig("desviointegral.png")
pylab.show()

#### ALÍNEA c.
print("ALÍENA c.")
print("... a carregar gráfico")
pylab.figure(figsize=(10,10))
pylab.plot([i for i in range(len(trabalho_impulsao_romberg[1]))],
           trabalho_impulsao_romberg[1],
           'ro', label="Método de Romberg")
pylab.yscale('log')
pylab.legend(prop={'size': 16})
pylab.xlabel("Número de iterações", fontsize=16)
pylab.ylabel("Desvio do valor real", fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.show()

#### ALÍNEA d.
print("ALÍENA d.")
print("... a carregar gráfico")
def inverse_square(x):
    return 1/(x**2)

int_inverse_square = IntegrationManager(
    function=inverse_square,
)
# Trabalho do integral do inverse square calculado no Integrate[]: 9.98004
# Trabalho do integral do inverse square calculado analiticamente: 9.98003992
valor_real = 9.98003992
numero_divisoes = [i for i in range(1,25)]
desvio_regra_trapezio = []
desvio_regra_simpson = []
desvio_metodo_romberg = []

int_inverse_square_romberg = int_inverse_square.romberg_method(a=0.1,b=50.1,k_max=25,convergence_criteria=10**(-6))

for n in range(1,25):
    desvio_regra_trapezio.append(int_inverse_square.trapezoidal_rule(a=0.1,b=50.1,n=n)-valor_real)
    desvio_regra_simpson.append(int_inverse_square.simpson_rule(a=0.1,b=50.1,n=n)[0]-valor_real)
    desvio_metodo_romberg.append(int_inverse_square_romberg[0][n][25]-valor_real)


pylab.loglog(numero_divisoes, desvio_regra_trapezio, 'ro', label="Regra do Trapézio")
pylab.loglog(numero_divisoes, desvio_regra_simpson, 'ko', label="Regra de Simpson")
pylab.loglog(numero_divisoes, desvio_metodo_romberg, 'bo', label="Método de Romberg")
pylab.legend(prop={'size': 20})
pylab.xlabel(r"$n$", fontsize=20)
pylab.ylabel(r"$d$", fontsize=20)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.savefig("devioinversesquare.png")
pylab.show()


#################################### PARTE 2 ###################################
print("PARTE 2")
print("ALÍENA a.")
x_values = [0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]
f_values = [0,37,71,104,134,161,185,207,225,239,250]


pylab.figure(figsize=(10,10))
pylab.plot(x_values,f_values, 'ro')
pylab.xlabel(r"x", fontsize=16)
pylab.ylabel(r"F", fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.show()


def trapezoidal_rule_discrete_data(x_data,y_data):
    n = len(y_data)-1
    h = (x_data[n]-x_data[0])/n

    sum_h = y_data[0]+y_data[n]
    for i in range(1,n):
        sum_h += 2*y_data[i]

    I = (h/2)*sum_h

    return I


print("O trabalho da força é",trapezoidal_rule_discrete_data(x_values,f_values))
velocidade_saida = np.sqrt((2*trapezoidal_rule_discrete_data(x_values,f_values))/0.075)
print("A velocidade de saída da flecha é:", velocidade_saida)


#### ALÍNEA c.
print("ALÍENA c.")
def force_interp_mathematica(x):
    # desenvolvimento de taylor em 0.25 de ordem 3 (source Mathematica Series[])
    return 162 + 510*(x-0.25)-600*(x-0.25)**2-(2.6148*(10**-11))*((x-0.25)**2)

trabalho_interp_mathematica = IntegrationManager(
    function=force_interp_mathematica,
)
trabalho_interp_mathematica_simpson = trabalho_interp_mathematica.simpson_rule(a=0,b=0.50,n=100)
velocidade_saida_interp_mathematica_simpson = np.sqrt((2*trabalho_interp_mathematica_simpson[0])/0.075)
print("""Trabalho da força interpolada pelo Mathematica, pela Regra de Simpson:{} ± {};
    traduz-se numa velocidade de saída da flecha de v={}""".format(trabalho_interp_mathematica_simpson[0],
    trabalho_interp_mathematica_simpson[1], velocidade_saida_interp_mathematica_simpson))


x_array = np.array([0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50])
y_array = np.array([0,37,71,104,134,161,185,207,225,239,250])
force_interp = interp1d(x_array,y_array, kind='cubic')

trabalho_interp = IntegrationManager(
    function=force_interp,
)

trabalho_interp_simpson = trabalho_interp.simpson_rule(a=0,b=0.50,n=100)
velocidade_saida_interp_simpson = np.sqrt((2*trabalho_interp_simpson[0])/0.075)
print("""Trabalho da força interpolada, pela Regra de Simpson:{} ± {};
    traduz-se numa velocidade de saída da flecha de v={}""".format(trabalho_interp_simpson[0],
    trabalho_interp_simpson[1], velocidade_saida_interp_simpson))


def simpson_rule_discrete_data(x_data, y_data):
    n = len(y_data)-1
    h = (x_data[n]-x_data[0])/n

    S_h = y_data[0]+y_data[n]
    for i in range(1,n,2):
        S_h += 4*y_data[i]
    for i in range(2,n,2):
        S_h += 2*y_data[i]
    S_h *= (h/3)

    S_2h = y_data[0]+y_data[n]
    for i in range(1,int((n+1)/2),2):
        S_2h += 4*y_data[2*i]
    for i in range(2,int(n/2),2):
        S_2h += 2*y_data[2*i]
    S_2h *= (2*h/3)
    erro = abs(S_h-S_2h)/15

    return S_h, erro

x_lista = [0.00,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50]
y_lista = [0,37,71,104,134,161,185,207,225,239,250]
simpson_discrete_result = simpson_rule_discrete_data(x_lista, y_lista)
velocidade_saida_simpson_discrete = np.sqrt((2*simpson_discrete_result[0])/0.075)
print("""Trabalho da força pela Regra de Simpson para dados discretos:{} ± {};
    traduz-se numa velocidade de saída da flecha de v={}""".format(simpson_discrete_result[0],
     simpson_discrete_result[1], velocidade_saida_simpson_discrete))

print("------------------------------------------------------------------")
print("\t\tTABELA DO CÁLCULO DO INTEGRAL DE SIMPSON")
print("\tmétodo\t\t\ttrabalho\t\tvelocidade")
print("scipy.interpolate.interp1d\t"
    +str(trabalho_interp_simpson[0])+"\t"+str(velocidade_saida_interp_simpson))
print("Series[Interpolate[]]\t"+"\t"
    +str(trabalho_interp_mathematica_simpson[0])+"\t"+str(velocidade_saida_interp_mathematica_simpson))
print("Simpson Discreto"+"\t\t"
    +str(simpson_discrete_result[0])+"\t"+str(velocidade_saida_simpson_discrete))
