### Author: Ernesto González
### Date: 21/11/2019


import numpy as np
import pylab



def euler_method(func, t_0, t_max, y_0, n):

    h = (t_max-t_0)/n
    y = [0 for j in range(n+1)]
    y[0] = y_0
    t = [t_0]

    for i in range(n):
        y[i+1] = y[i] + func(y[i])*h
        t_0 += h
        t.append(t_0)

    return t, y


def function(y, k=-2.3):
    return y*k



#### ALÍNEA a.
first_h_results = euler_method(func=function, t_0=0, t_max=5, y_0=1, n=10)
second_h_results = euler_method(func=function, t_0=0, t_max=5, y_0=1, n=5)
third_h_results = euler_method(func=function, t_0=0, t_max=5, y_0=1, n=7)

t = 0
t_values = []
y_values = []
while t<5.01:
    t_values.append(t)
    y_values.append(np.exp(-2.3*t))
    t += 0.01

pylab.plot(first_h_results[0], first_h_results[1], '--go', label=r"$h=0.5$")
pylab.plot(second_h_results[0], second_h_results[1], '--ro', label=r"$h=1$")
pylab.plot(third_h_results[0], third_h_results[1], '--bo', label=r"$h=0.7$")
pylab.plot(t_values, y_values, '-k', label=r"$N(t)=-2.3 e^t$")
pylab.legend(prop={'size': 15})
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$N(t)$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

#### ALÍNEA b.

def fourth_order_runge_kutta_method(func, t_0, t_max, y_0, n):

    h = (t_max-t_0)/n
    y = [0 for j in range(n+1)]
    y[0] = y_0
    t = [t_0]

    for i in range(n):
        k_1 = func(y[i])*h
        k_2 = func(y[i] + k_1/2)*h
        k_3 = func(y[i] + k_2/2)*h
        k_4 = func(y[i] + k_3)*h
        y[i+1] = y[i] + (k_1/6) + (k_2/3) + (k_3/3) + (k_4/4)
        t_0 += h
        t.append(t_0)

    return t, y

first_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=5)
second_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=10)
third_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=20)
fourth_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=40)
fifth_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=80)
sixth_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=160)
seventh_euler_result = euler_method(func=function,t_0=0,t_max=5,y_0=1,n=320)

first_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=5)
second_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=10)
third_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=20)
fourth_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=40)
fifth_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=80)
sixth_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=160)
seventh_runge_kutta_result = fourth_order_runge_kutta_method(func=function,t_0=0,t_max=5,y_0=1,n=320)

euler_results = [
    first_euler_result,
    second_euler_result,
    third_euler_result,
    fourth_euler_result,
    fifth_euler_result,
    sixth_euler_result,
    seventh_euler_result
]
runge_kutta_results = [
    first_runge_kutta_result,
    second_runge_kutta_result,
    third_runge_kutta_result,
    fourth_runge_kutta_result,
    fifth_runge_kutta_result,
    sixth_runge_kutta_result,
    seventh_runge_kutta_result
]

euler_results_final_deviation_per_h = []
runge_kutta_results_final_deviation_per_h = []

for i in range(7):
    euler_results_final_deviation_per_h.append(
        abs(euler_results[i][1][-1] - np.exp(-2.3*euler_results[i][0][-1])))
    runge_kutta_results_final_deviation_per_h.append(
        abs(runge_kutta_results[i][1][-1] - np.exp(-2.3*runge_kutta_results[i][0][-1])))

h_values = [1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64]

pylab.loglog(h_values,
    euler_results_final_deviation_per_h, 'go', label=r"Euler Method")
pylab.loglog(h_values,
    runge_kutta_results_final_deviation_per_h, 'ro', label=r"Runge-Kutta Method")
pylab.legend(prop={'size': 15})
pylab.xlabel("$h$", fontsize=15)
pylab.ylabel("$d$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

print("Orden de convergência:")
print("\tEuler:")
print("\t\t",(np.log(1.0124327094724461e-05)-np.log(1.9329826037131946e-06))/(np.log(0.5)-np.log(0.015625)))
print("\t\t",(np.log(3.712940130093597)-np.log(1.0124327094724461e-05))/(np.log(1)-np.log(0.5)))
print("\tRunge-Kutta:")
print("\t\t",(np.log(5.488213686921787e-06)-np.log(1.6685272683902158e-12))/(np.log(0.5)-np.log(0.015625)))
print("\t\t",(np.log(0.026323064644376577)-np.log(5.488213686921787e-06))/(np.log(1)-np.log(0.5)))

#### ALÍNEA c.

def zombieland_simulation_euler(t_0, t_max, H_0, Z_0, n, a, b, c):

    h = (t_max-t_0)/n
    H = [0 for j in range(n)]
    Z = [0 for j in range(n)]
    H[0], Z[0] = H_0, Z_0
    t = [t_0]

    for i in range(n-1):
        H[i+1] = H[i] + (- b*H[i]*Z[i] - c*H[i]*Z[i])*h
        Z[i+1] = Z[i] + (c*H[i]*Z[i] - a*H[i]*Z[i])*h
        t_0 += h
        t.append(t_0)

    return t, H, Z


simulation_1_euler = zombieland_simulation_euler(
    t_0=0,t_max=10,
    H_0=90,Z_0=30,
    n=1000,
    a=0.05,b=0.06,c=0.02)
simulation_2_euler = zombieland_simulation_euler(
    t_0=0,t_max=10,
    H_0=55,Z_0=5,
    n=1000,
    a=0.05,b=0.06,c=0.02)
simulation_3_euler = zombieland_simulation_euler(
    t_0=0,t_max=10,
    H_0=70,Z_0=20,
    n=1000,
    a=0.05,b=0.06,c=0.02)

pylab.plot(simulation_1_euler[1], simulation_1_euler[2],
    'ro', label=r"$H_0=90, Z_0=30$")
pylab.plot(simulation_2_euler[1], simulation_2_euler[2],
    'bo', label=r"$H_0=55, Z_0=5$")
pylab.plot(simulation_3_euler[1], simulation_3_euler[2],
    'go', label=r"$H_0=70, Z_0=20$")
pylab.legend(prop={'size': 15})
pylab.xlabel("$H$", fontsize=15)
pylab.ylabel("$Z$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("epidemiazombieeuler.png")
pylab.show()


def zombieland_simulation_runge_kutta(t_0, t_max, H_0, Z_0, n, a, b, c):

    h = (t_max-t_0)/n
    H = [0 for j in range(n)]
    Z = [0 for j in range(n)]
    H[0], Z[0] = H_0, Z_0
    t = [t_0]

    def H_jump(H,Z):
        return -b*H*Z - c*H*Z

    def Z_jump(H,Z):
        return c*H*Z - a*H*Z

    for i in range(n-1):
        k_1 = H_jump(H[i],Z[i])*h
        k_2 = H_jump(H[i]+k_1/2,Z[i]+k_1/2)*h
        k_3 = H_jump(H[i]+k_2/2,Z[i]+k_2/2)*h
        k_4 = H_jump(H[i]+k_3,Z[i]+k_3)*h
        H[i+1] = H[i] + (k_1/6) + (k_2/3) + (k_3/3) + (k_4/4)

        k_1 = Z_jump(H[i],Z[i])*h
        k_2 = Z_jump(H[i]+k_1/2,Z[i]+k_1/2)*h
        k_3 = Z_jump(H[i]+k_2/2,Z[i]+k_2/2)*h
        k_4 = Z_jump(H[i]+k_3,Z[i]+k_3)*h
        Z[i+1] = Z[i] + (k_1/6) + (k_2/3) + (k_3/3) + (k_4/4)
        t_0 += h
        t.append(t_0)

    return t, H, Z

simulation_1 = zombieland_simulation_runge_kutta(
    t_0=0,t_max=10,
    H_0=120,Z_0=80,
    n=1000,
    a=0.01,b=0.03,c=0.04)
simulation_2 = zombieland_simulation_runge_kutta(
    t_0=0,t_max=10,
    H_0=140,Z_0=60,
    n=1000,
    a=0.01,b=0.03,c=0.04)
simulation_3 = zombieland_simulation_runge_kutta(
    t_0=0,t_max=10,
    H_0=90,Z_0=100,
    n=1000,
    a=0.01,b=0.03,c=0.04)


pylab.plot(simulation_1[1], simulation_1[2],
    'ro', label=r"$H_0=120, Z_0=80$")
pylab.plot(simulation_2[1], simulation_2[2],
    'bo', label=r"$H_0=140, Z_0=60$")
pylab.plot(simulation_3[1], simulation_3[2],
    'go', label=r"$H_0=90, Z_0=100$")
pylab.legend(prop={'size': 15})
pylab.xlabel("$H$", fontsize=15)
pylab.ylabel("$Z$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("epidemiazombierk.png")
pylab.show()



####################### PARTE 2 ################################################
#### ALÍNEA a.

def dvx_dt(v_x, x, t):
    return 0


def dx_dt(v_x, x, t):
    return v_x


def dvy_dt(v_y, y, t):
    return -9.81


def dy_dt(v_y, y, t):
    return v_y


def euler_method_for_system(func_x, func_y, d_func_x, d_func_y, t_0, t_max, x_0, y_0, v_x_0, v_y_0, n):

    h = (t_max-t_0)/n
    x = [0 for j in range(n+1)]
    y = [0 for j in range(n+1)]
    v_x = [0 for j in range(n+1)]
    v_y = [0 for j in range(n+1)]
    x[0], v_x[0], y[0], v_y[0] = x_0, v_x_0, y_0, v_y_0
    t = [t_0]

    for i in range(n):
        x[i+1] = x[i] + func_x(v_x[i], x[i],t[i])*h
        v_x[i+1] = v_x[i] + d_func_x(v_x[i], x[i],t[i])*h

        y[i+1] = y[i] + func_y(v_y[i], y[i],t[i])*h
        v_y[i+1] = v_y[i] + d_func_y(v_y[i], y[i],t[i])*h

        t_0 += h
        t.append(t_0)

    return t, x, y


def fourth_order_runge_kutta_method_for_system(func_x, func_y, d_func_x, d_func_y, t_0, t_max, x_0, y_0, v_x_0, v_y_0, n):

    h = (t_max-t_0)/n
    x = [0 for j in range(n+1)]
    y = [0 for j in range(n+1)]
    v_x = [0 for j in range(n+1)]
    v_y = [0 for j in range(n+1)]
    x[0], v_x[0], y[0], v_y[0] = x_0, v_x_0, y_0, v_y_0
    t = [t_0]

    for i in range(n):
        k_11 = func_x(v_x[i], x[i],t[i])
        k_21 = d_func_x(v_x[i], x[i],t[i])
        k_12 = func_x(v_x[i]+(h/2)*k_21, x[i]+(h/2)*k_11,t[i])
        k_22 = d_func_x(v_x[i]+(h/2)*k_21, x[i]+(h/2)*k_11,t[i])
        k_13 = func_x(v_x[i]+(h/2)*k_22, x[i]+(h/2)*k_12,t[i])
        k_23 = d_func_x(v_x[i]+(h/2)*k_22, x[i]+(h/2)*k_12,t[i])
        k_14 = func_x(v_x[i]+(h)*k_23, x[i]+(h)*k_13,t[i])
        k_24 = d_func_x(v_x[i]+(h)*k_23, x[i]+(h)*k_13,t[i])

        x[i+1] = x[i] + (h/6)*(k_11+2*k_12+2*k_13+k_14)
        v_x[i+1] = v_x[i] + (h/6)*(k_21+2*k_22+2*k_23+k_24)

        k_11 = func_y(v_y[i], y[i],t[i])
        k_21 = d_func_y(v_y[i], y[i],t[i])
        k_12 = func_y(v_y[i]+(h/2)*k_21, y[i]+(h/2)*k_11,t[i])
        k_22 = d_func_y(v_y[i]+(h/2)*k_21, y[i]+(h/2)*k_11,t[i])
        k_13 = func_y(v_y[i]+(h/2)*k_22, y[i]+(h/2)*k_12,t[i])
        k_23 = d_func_y(v_y[i]+(h/2)*k_22, y[i]+(h/2)*k_12,t[i])
        k_14 = func_y(v_y[i]+(h)*k_23, y[i]+(h)*k_13,t[i])
        k_24 = d_func_y(v_y[i]+(h)*k_23, y[i]+(h)*k_13,t[i])

        y[i+1] = y[i] + (h/6)*(k_11+2*k_12+2*k_13+k_14)
        v_y[i+1] = v_y[i] + (h/6)*(k_21+2*k_22+2*k_23+k_24)

        t_0 += h
        t.append(t_0)

    return t, x, y


system_result_euler = euler_method_for_system(
    func_x=dx_dt, func_y=dy_dt,
    d_func_x=dvx_dt, d_func_y=dvy_dt,
    t_0=0, t_max=5,
    x_0=0, y_0=0,
    v_x_0=20*np.cos(np.pi/4), v_y_0=20*np.sin(np.pi/4),
    n=10)
system_result_rk = fourth_order_runge_kutta_method_for_system(
    func_x=dx_dt, func_y=dy_dt,
    d_func_x=dvx_dt, d_func_y=dvy_dt,
    t_0=0, t_max=5,
    x_0=0, y_0=0,
    v_x_0=20*np.cos(np.pi/4), v_y_0=20*np.sin(np.pi/4),
    n=10)


x_values = []
y_values = []
t = 0
while t < 5.001:
    x_values.append(20*np.cos(np.pi/4)*t)
    y_values.append(20*np.sin(np.pi/4)*t + 0.5*(-9.81)*(t**2))
    t += 0.001

pylab.plot(system_result_euler[1], system_result_euler[2], 'go', label=r"Euler")
pylab.plot(system_result_rk[1], system_result_rk[2], 'ro', label=r"Runge-Kutta")
pylab.plot(x_values, y_values, '-k', label=r"Analítica")
pylab.legend(prop={'size': 15})
pylab.xlabel("$x$", fontsize=15)
pylab.ylabel("$y$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("trajetoriaprojetil.png")
pylab.show()
print("Euler", system_result_euler[2][-1])
print("Runge-kutta", system_result_rk[2][-1])
print("Analítica", y_values[-1])
#### ALÍNEA b.
gamma=0.1

def dvx_dt_dissipativo(v_x, x, t):
    return -gamma*(v_x**2)


def dx_dt_dissipativo(v_x, x, t):
    return v_x


def dvy_dt_dissipativo(v_y, y, t):
    if v_y >= 0:
        return -9.81 - gamma*(v_y**2)
    if v_y < 0:
        return -9.81 + gamma*(v_y**2)


def dy_dt_dissipativo(v_y, y, t):
    return v_y


system_result_euler_dissipativo = euler_method_for_system(
    func_x=dx_dt_dissipativo, func_y=dy_dt_dissipativo,
    d_func_x=dvx_dt_dissipativo, d_func_y=dvy_dt_dissipativo,
    t_0=0, t_max=5,
    x_0=0, y_0=0,
    v_x_0=20*np.cos(np.pi/4), v_y_0=20*np.sin(np.pi/4),
    n=10)
system_result_rk_dissipativo = fourth_order_runge_kutta_method_for_system(
    func_x=dx_dt_dissipativo, func_y=dy_dt_dissipativo,
    d_func_x=dvx_dt_dissipativo, d_func_y=dvy_dt_dissipativo,
    t_0=0, t_max=5,
    x_0=0, y_0=0,
    v_x_0=20*np.cos(np.pi/4), v_y_0=20*np.sin(np.pi/4),
    n=10)


pylab.plot(system_result_euler_dissipativo[1], system_result_euler_dissipativo[2], 'go', label=r"Euler")
pylab.plot(system_result_rk_dissipativo[1], system_result_rk_dissipativo[2], 'ro', label=r"Runge Kutta")
pylab.legend(prop={'size': 15})
pylab.xlabel("$x$", fontsize=15)
pylab.ylabel("$y$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.savefig("trajetoriaprojetildissipativo.png")
pylab.show()
