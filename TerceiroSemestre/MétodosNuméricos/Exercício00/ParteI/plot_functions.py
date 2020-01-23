import math
import pylab

# Change font size in plots
pylab.rcParams.update({'font.size': 22})

# Y=X**2 FUNCTION
x_points = []
y_points = []
for i in range(10):
    x_points.append(i)
    y_points.append(i*i)
pylab.figure(figsize=(10,10))
pylab.plot(x_points,y_points)
pylab.title("Parabolic Function")
pylab.xlabel("i")
pylab.ylabel("i**2")
pylab.savefig("parabola.pdf")
# pylab.savefig("parabola.png")

# SINE FUNCTION
x_points = []
y_points = []
x = 0.0
while x<10:
    x_points.append(x)
    y_points.append(math.sin(x))
    x += 0.1
pylab.figure(figsize=(13,10))
pylab.plot(x_points,y_points)
pylab.title("Sine Function")
pylab.xlabel("x")
pylab.ylabel("sin(x)")
pylab.savefig("sin.pdf")
# pylab.savefig("sin.png")

# LOG FUNCTION
x_points = []
y_points = []
x = 0.1
while x<20:
    x_points.append(x)
    y_points.append(math.log(x))
    x += 0.1
pylab.figure(figsize=(10,10))
pylab.plot(x_points,y_points)
pylab.title("Logarithmic Function")
pylab.xlabel("x")
pylab.ylabel("log(x)")
pylab.savefig("log.pdf")
# pylab.savefig("log.png")

# EXP FUNCTION
x_points = []
y_points = []
x = -10
while x<20:
    x_points.append(x)
    y_points.append(math.exp(x))
    x += 0.1
pylab.figure(figsize=(10,10))
pylab.plot(x_points,y_points)
pylab.title("Exponential Function")
pylab.xlabel("x")
pylab.ylabel("exp(x)")
pylab.savefig("exp.pdf")
# pylab.savefig("exp.png")
