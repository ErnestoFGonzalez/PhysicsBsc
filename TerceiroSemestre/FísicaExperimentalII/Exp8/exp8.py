import csv
import pylab


data_volume = []

with open("adiabaticaroldana.txt", "r") as file:
    lis = [line.split() for line in file]        # create a list of lists
    for i, x in enumerate(lis):
        x[1] = x[1].replace(',', '.')
        data_volume.append(float(x[1]))


data_pressure = []

with open("adiabaticapressao.txt", "r") as file:
    lis = [line.split() for line in file]        # create a list of lists
    for i, x in enumerate(lis):
        x[1] = x[1].replace(',', '.')
        data_pressure.append(float(x[1]))


pylab.plot(data_volume, data_pressure, '--go')
pylab.xlabel("$V$", fontsize=15)
pylab.ylabel("$p$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()
