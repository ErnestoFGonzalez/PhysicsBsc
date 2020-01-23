#### Author: Ernesto González 52857
#### Date: 12/10/2019


import pylab


class DerivativeApproxManager():
    def __init__(self, function, function_derivative, x_point, h_values):
        self.function = function
        self.function_derivative = function_derivative
        self.x = x_point
        self.h_values = h_values

    def approx_function_derivative(self, h):
        return (self.function(self.x+h)-self.function(self.x))/h

    def approx_abs_error(self, h):
        return abs(self.function_derivative(self.x)-self.approx_function_derivative(h))

    def absolute_error_per_h(self, plot_title, x_label, y_label, plot_filename):
        absolute_error_per_h = []
        for h in self.h_values:
            absolute_error_per_h.append(self.approx_abs_error(h))

        # Plot absolute_error_per_h per h in log scale axis
        pylab.figure(figsize=(10,10))
        pylab.loglog(self.h_values, absolute_error_per_h, 'ro')
        pylab.grid(color='k', linestyle='--', linewidth=0.5)
        # pylab.title(plot_title, fontsize=20, fontweight='bold')
        pylab.xticks(fontsize=22)
        pylab.yticks(fontsize=22)
        pylab.xlabel(x_label, fontsize=23)
        pylab.ylabel(y_label, fontsize=23)
        pylab.savefig(plot_filename)
        pylab.show()


#####################APLICAÇÃO######################
# Define a função a ser usada neste exemplo
def f(x):
    return x**2
def f_linha(x):
    return 2*x

# Primeiro caso: h com valores em decimais
valores_h10 = [10**(-i) for i in range(21)]
erros_base10 = DerivativeApproxManager(x_point=1,
                                       function=f,
                                       function_derivative=f_linha,
                                       h_values=valores_h10)
erros_base10.absolute_error_per_h(
    plot_title="Erros Absolutos de Aproximações à primeira derivada\n da função quadrática usando h de base decimal",
    x_label="h",
    y_label="Erro Absoluto",
    plot_filename="erroh10.png")

# Segundo caso: h com valores em base binária
valores_h2 = [2**(-i) for i in range(61)]
erros_base2 = DerivativeApproxManager(x_point=1,
                                       function=f,
                                       function_derivative=f_linha,
                                       h_values=valores_h2)
erros_base2.absolute_error_per_h(
    plot_title="Erros Absolutos de Aproximações à primeira derivada\n da função quadrática usando h de base binária",
    x_label="h",
    y_label="Erro Absoluto",
    plot_filename="erroh2.png")
