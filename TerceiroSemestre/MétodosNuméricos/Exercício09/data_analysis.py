# Author: Ernesto González
# Date: 05/12/2019


import numpy as np
import pylab
from PIL import Image
import optimization


def average_value(array):
    sum = 0
    for value in array:
        sum += value
    average = sum/len(array)
    return average


def median_value(array):
    sorted_list = sorted(array)
    list_length = len(array)
    index = (list_length - 1) // 2

    if (list_length % 2):
        return sorted_list[index]
    else:
        return (sorted_list[index] + sorted_list[index + 1])/2.0


def variance(array):
    average = average_value(array)
    desvio_sum = 0
    for value in array:
        desvio_sum += (value - average)**2
    variancia = desvio_sum/len(array)
    return variancia


def covariance(X, Y):
    X_deviated = []
    x_mean = average_value(X)
    for x in X:
        X_deviated.append(x-x_mean)
    Y_deviated = []
    y_mean = average_value(Y)
    for y in Y:
        Y_deviated.append(y-y_mean)
    return average_value([X_deviated[i]*Y_deviated[i] for i in range(len(X_deviated))])


def linear_regression(X, Y):
    n = len(X)
    # determine m
    m = covariance(X,Y)/variance(X)
    # determine b
    b = average_value(Y) - m*average_value(X)
    return m, b


############################### PARTE 1 ########################################
with open('metabol.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]

    data_frame_metabol = []

    for line in readlines:
        temp = line.split("\t")
        data_frame_metabol.append([float(temp[0]), float(temp[1])])


m, b = linear_regression(
        X=[np.log(data_frame_metabol[i][0]) for i in range(len(data_frame_metabol))],
        Y=[np.log(data_frame_metabol[i][1]) for i in range(len(data_frame_metabol))]
        )

log_X_fit = np.linspace(-2,6,100)
log_Y_fit = m*log_X_fit + b

pylab.plot([np.log(data_frame_metabol[i][0]) for i in range(len(data_frame_metabol))],
            [np.log(data_frame_metabol[i][1]) for i in range(len(data_frame_metabol))],
            'ko')
pylab.plot(log_X_fit, log_Y_fit, '-b')
pylab.text(-1,0,r"$log(W) = {}\,log(m) + {}$".format(m,b), fontsize=16)
pylab.xlabel("$log(m)$", fontsize=18)
pylab.ylabel("$log(W)$", fontsize=18)
pylab.xticks(fontsize=18)
pylab.yticks(fontsize=18)
pylab.show()

real_X_fit = np.linspace(0,400,1000)
real_Y_fit = np.exp(b)*(real_X_fit**m)

pylab.plot([data_frame_metabol[i][0] for i in range(len(data_frame_metabol))],
            [data_frame_metabol[i][1] for i in range(len(data_frame_metabol))],
            'ko')
pylab.plot(real_X_fit, real_Y_fit, '-b')
pylab.text(50,10,r"$m = e^{"+str(b)+"}\cdot W^{"+str(m)+"}$", fontsize=20)
pylab.xlabel("$m$", fontsize=18)
pylab.ylabel("$W$", fontsize=18)
pylab.xticks(fontsize=18)
pylab.yticks(fontsize=18)
pylab.show()


def S(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (Y[i]-m*X[i]-b)**2
    return sum


def dS_dm(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (Y[i]-m*X[i]-b)*X[i]
    sum *= -2
    return sum


def dS_db(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += Y[i]-m*X[i]-b
    sum *= -2
    return sum


regression_manager = optimization.OptimizationManager2D(
    function=S,
    function_d_dx=dS_dm,
    function_d_dy=dS_db
)

m_guess = 0.5
b_guess = 4

regression_gradient_result = regression_manager.gradient_method_2D(
    x_0=m_guess, y_0=b_guess, lambd=0.01, k_max=100, convergence_criteria=10**(-6),
    X=[np.log(data_frame_metabol[i][0]) for i in range(len(data_frame_metabol))],
    Y=[np.log(data_frame_metabol[i][1]) for i in range(len(data_frame_metabol))])

lambda_values = [0.01,0.05,0.1]
gradient_method_results_per_lambda = []

for lambd in lambda_values:
    gradient_results = regression_manager.gradient_method_2D(
        x_0=m_guess, y_0=b_guess, lambd=lambd, k_max=100, convergence_criteria=10**(-6),
        X=[np.log(data_frame_metabol[i][0]) for i in range(len(data_frame_metabol))],
        Y=[np.log(data_frame_metabol[i][1]) for i in range(len(data_frame_metabol))])
    gradient_method_results_per_lambda.append([lambd, gradient_results[0], gradient_results[1]])

for line in gradient_method_results_per_lambda:
    print(line)


def S2(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (Y[i] - b*(X[i]**m))**2
    return sum


def dS2_dm(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (Y[i] - b*(X[i]**m))*(X[i]**m)*np.log(X[i])
    sum = -2*b*sum
    return sum


def dS2_db(m, b, X, Y):
    n = len(X)
    sum = 0
    for i in range(n):
        sum += (Y[i] - b*(X[i]**m))*(X[i]**m)
    sum *= -2
    return sum


regression_manager2 = optimization.OptimizationManager2D(
    function=S2,
    function_d_dx=dS2_dm,
    function_d_dy=dS2_db
)

regression_gradient_result2 = regression_manager2.gradient_method_2D(
    x_0=3, y_0=1, lambd=0.001, k_max=200, convergence_criteria=10**(-7),
    X=[data_frame_metabol[i][0] for i in range(len(data_frame_metabol))],
    Y=[data_frame_metabol[i][1] for i in range(len(data_frame_metabol))])

print("$W=Cm^{\alpha}$")
print(r"$\alpha$ é",regression_gradient_result2[2])
print(r"$C$ é",regression_gradient_result2[3])
################################################################################

################################ PARTE 2 #######################################
img = Image.open('rocks.jpg')
rgb_matrix = np.array(img)

gray_scale_matrix = []

for line in rgb_matrix:
    line_gray_scale = []
    for pixel in line:
        R = int(pixel[0])
        G = int(pixel[1])
        B = int(pixel[2])
        I = (R+G+B)/(3*255)
        line_gray_scale.append(I)
    gray_scale_matrix.append(line_gray_scale)

gray_scale_array = np.asarray(gray_scale_matrix)

gray_image = Image.fromarray((gray_scale_array*255).astype('uint8'), mode='L')
gray_image.save('rocksgrayscale.png')
gray_image.show()

gray_image_pixels = []
for line in gray_scale_matrix:
    for pixel in line:
        gray_image_pixels.append(pixel)


n, bins, patches = pylab.hist(gray_image_pixels, 20, facecolor='b', edgecolor='black',
                              linewidth=2, align='mid', range=[0,1])
pylab.xlabel(r'$I$', fontsize=15)
pylab.ylabel(r'$N$', fontsize=15)
pylab.axis([0, 1, 0, 11000])
pylab.xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
pylab.show()


threshold_point = 0.25
thresholded_gray_scale_matrix = [[] for line in gray_scale_matrix]
for line in gray_scale_matrix:
    index = gray_scale_matrix.index(line)
    for pixel in line:
        if pixel < threshold_point:
            thresholded_gray_scale_matrix[index].append(0)
        if pixel >= threshold_point:
            thresholded_gray_scale_matrix[index].append(1)

thresholded_gray_scale_array = np.asarray(thresholded_gray_scale_matrix)
thresholded_gray_image = Image.fromarray((thresholded_gray_scale_array*255).astype('uint8'), mode='L')
thresholded_gray_image.show()

# empty space is given by black pixels in thresholded_gray_image
counter_of_empty_space = 0
counter_total_space = 0
for line in thresholded_gray_scale_matrix:
    for pixel in line:
        if pixel == 1:
            counter_of_empty_space += 1
        counter_total_space +=1

porosity_of_material = (counter_of_empty_space/counter_total_space)*100
print("A porosidade do material é",porosity_of_material)


# Edge detection
def L_x(a11,a12,a13,a21,a22,a23,a31,a32,a33):
    return a11-a13+2*a21-2*a23+a31-a33


def L_y(a11,a12,a13,a21,a22,a23,a31,a32,a33):
    return a11+2*a12+a13-a31-2*a32-a33


A = gray_scale_array
edge_detected_matrix = []
for i in range(len(gray_scale_array)):
    line = []
    if i+1 == len(gray_scale_array):
        i_plus_one = 0
    else:
        i_plus_one = i+1

    for j in range(len(gray_scale_array[i])):
        if j+1 == len(gray_scale_array[i]):
            j_plus_one = 0
        else:
            j_plus_one = j+1

        Lx = L_x(A[i-1][j-1],A[i-1][j],A[i-1][j_plus_one],
                A[i][j-1],A[i][j],A[i][j_plus_one],
                A[i_plus_one][j-1],A[i_plus_one][j],A[i_plus_one][j_plus_one],)
        Ly = L_y(A[i-1][j-1],A[i-1][j],A[i-1][j_plus_one],
                A[i][j-1],A[i][j],A[i][j_plus_one],
                A[i_plus_one][j-1],A[i_plus_one][j],A[i_plus_one][j_plus_one],)
        L = np.sqrt(Lx**2+Ly**2)
        line.append(L)
    edge_detected_matrix.append(line)



edge_detected_array = np.asarray(edge_detected_matrix)
edge_detected_image= Image.fromarray((edge_detected_array*255).astype('uint8'), mode='L')
# edge_detected_image.save('edgedetected.png')
edge_detected_image.show()
