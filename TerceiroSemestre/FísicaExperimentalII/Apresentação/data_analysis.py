import numpy as np
from sklearn.linear_model import LinearRegression
import pylab


def linear_regression_line_coordinates(x,y,min_x,max_x):
    x = np.array(x).reshape(-1,1)
    y = np.array(y)
    model = LinearRegression().fit(x,y)
    r_sq = model.score(x,y)
    slope = model.coef_
    intercept = model.intercept_
    x_line = np.linspace(min_x,max_x, 100)
    y_line = slope*x_line + intercept
    return x_line,y_line


with open('rt_metal.csv', 'r') as file:
    readlines = file.readlines()

    df_metal = []

    for line in readlines:
        temp = line.split('\t')
        temp[0] = float(temp[0].replace(',','.'))
        temp[1] = float(temp[1].replace('\n','').replace(',','.'))

        if temp[1] != 0:
            df_metal.append([temp[0], temp[1]])


with open('rt_semicondutor.csv', 'r') as file:
    readlines = file.readlines()

    df_semicondutor = []

    for line in readlines:
        temp = line.split('\t')
        temp[0] = float(temp[0].replace(',','.'))
        temp[1] = float(temp[1].replace('\n','').replace(',','.'))

        if temp[1] != 0:
            df_semicondutor.append([temp[0], temp[1]])


pylab.plot([line[0] for line in df_metal], [line[1] for line in df_metal], 'bo', label='Metal')
lin_reg_metal = linear_regression_line_coordinates([line[0] for line in df_metal], [line[1] for line in df_metal], 278,365)
pylab.plot(lin_reg_metal[0], lin_reg_metal[1], 'b')
pylab.plot([line[0] for line in df_semicondutor], [line[1] for line in df_semicondutor], 'ro', label='Semicondutor')
# pylab.text(330,93,'$R(T)=0.3869T-17.604$', fontsize=24)
# pylab.text(330,80,'$R^2=0.9998$', fontsize=24)
# pylab.text(330,69,r'$\sigma=4.3261 \times 10^{-4}$', fontsize=24)
pylab.legend(fontsize=24)
pylab.xlabel("$T\,(K)$", fontsize=24)
pylab.ylabel("$R\,(\Omega)$", fontsize=24)
pylab.xticks(fontsize=24)
pylab.yticks(fontsize=24)
pylab.show()


with open('rt_semicondutor_logado.csv', 'r') as file:
    readlines = file.readlines()

    df_semicondutor_logado = []

    for line in readlines:
        temp = line.split('\t')

        if len(temp) != 1:
            temp[0] = float(temp[0].replace(',','.'))
            temp[1] = float(temp[1].replace('\n','').replace(',','.'))

            df_semicondutor_logado.append([temp[0], temp[1]])


pylab.plot([line[0] for line in df_semicondutor_logado],
           [line[1] for line in df_semicondutor_logado], 'ro', label='Semicondutor')
lin_reg_semicondutor = linear_regression_line_coordinates(
    [line[0] for line in df_semicondutor_logado],
    [line[1] for line in df_semicondutor_logado], 0.0027,0.0036)
pylab.plot(lin_reg_semicondutor[0], lin_reg_semicondutor[1], 'r')
pylab.text(0.003,2.2,r'$ln(R)=1461.7\frac{1}{T}-2.8686$', fontsize=24)
pylab.text(0.003,2.1,r'$R^2=0.9999$', fontsize=24)
pylab.text(0.003,2.0,r'$\sigma=1.1556$', fontsize=24)
pylab.legend(fontsize=24)
pylab.xlabel("$1/T\,(K^{-1})$", fontsize=24)
pylab.ylabel("$ln(R)\,(ln(\Omega))$", fontsize=24)
pylab.xticks(fontsize=24)
pylab.yticks(fontsize=24)
pylab.show()


with open('termopar.csv', 'r') as file:
    readlines = file.readlines()

    # df = [[Temperature, Emf]]
    df = []

    for line in readlines:
        temp = line.split('\t')

        temp[0] = temp[0].replace(',','.')
        temp[1] = temp[1].replace(',','.')

        temp[0] = float(temp[0])
        temp[1] = float(temp[1].replace('\n',''))

        df.append([temp[0], temp[1]])


pylab.plot([line[0] for line in df], [line[1] for line in df], 'bo')
lin_reg = linear_regression_line_coordinates([line[0] for line in df],[line[1] for line in df],14,70)
pylab.plot(lin_reg[0], lin_reg[1], 'b')
pylab.text(20,2,'$\epsilon=0.0374\Delta T+0.0218$', fontsize=24)
pylab.text(20,1.8,'$R^2=0.9892$', fontsize=24)
pylab.text(20,1.6,r'$\sigma=1.0444 \times 10^{-3}$', fontsize=24)
pylab.xlabel("$\Delta T \, (K)$", fontsize=24)
pylab.ylabel("$\epsilon \,(V)$", fontsize=24)
pylab.xticks(fontsize=24)
pylab.yticks(fontsize=24)
pylab.show()
