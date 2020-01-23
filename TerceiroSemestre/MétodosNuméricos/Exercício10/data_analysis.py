# Author: Ernesto González
# Date: 12/12/2019


import numpy as np
import random
import pylab
from mpl_toolkits.mplot3d import Axes3D
import csv
from PIL import Image



################################################################################
############################## PARTE 1 #########################################

def kMeans(data, k_max, max_iterations):
    N = len(data)

    # algorithm works only for multidimensional datapoints
    # this blocks converts unidimensional datapoints to multidimensional
    try:
        unidimensional = False
        len(data[0])
    except:
        unidimensional = True
        for i in range(len(data)):
            data[i] = [data[i], 0]

    dim = len(data[0])

    centroids = [data[random.randint(0,N-1)] for k in range(k_max)]
    group = [0 for i in range(N)]

    for n in range(max_iterations+1):
        count = [0 for k in range(k_max)]

        # for each datapoint find closest centroid and assign datapoint
        # to relevant cluster
        for i in range(N):
            k_min = 0
            d_min = 1000000
            for k in range(k_max):
                dist = 0
                for d in range(dim):
                    dist += (data[i][d]-centroids[k][d])**2
                dist **= 0.5
                if dist < d_min:
                    d_min = dist
                    k_min = k
            group[i] = k_min
            count[k_min] += 1

        # reposition centroids to center of mass of found cluster
        centroids = [[0 for i in range(dim)] for k in range(k_max)]
        for i in range(N):
            for d in range(dim):
                centroids[group[i]][d] += data[i][d]/count[group[i]]

    # sum of the distance of each point of the cluster to respective centroid
    sum_distance_to_centroids = []
    for k in range(k_max):
        sum_distances = 0
        for i in range(N):
            if group[i] == k:
                dist = 0
                for d in range(dim):
                    dist += (data[i][d]-centroids[k][d])**2
                # dist **= 0.5
                sum_distances += dist
        sum_distance_to_centroids.append(sum_distances)

    if unidimensional:
        for i in range(len(centroids)):
            centroids[i] = centroids[i][0]
        for i in range(len(data)):
            data[i] = data[i][0]

    return centroids, group, sum_distance_to_centroids


with open('finland.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]

    players_coordinates = []

    for line in readlines:
        temp = line.split('\t')
        temp[0] = int(temp[0])
        temp[1] = int(temp[1])
        players_coordinates.append([temp[1], temp[0]])


kmeans_2clusters = kMeans(players_coordinates, k_max=2, max_iterations=10)
kmeans_3clusters = kMeans(players_coordinates, k_max=3, max_iterations=10)

kmeans_2clusters_cluster1_points = []
kmeans_2clusters_cluster2_points = []
for i in range(len(players_coordinates)):
    if kmeans_2clusters[1][i] == 0:
        kmeans_2clusters_cluster1_points.append(players_coordinates[i])
    if kmeans_2clusters[1][i] == 1:
        kmeans_2clusters_cluster2_points.append(players_coordinates[i])

kmeans_3clusters_cluster1_points = []
kmeans_3clusters_cluster2_points = []
kmeans_3clusters_cluster3_points = []
for i in range(len(players_coordinates)):
    if kmeans_3clusters[1][i] == 0:
        kmeans_3clusters_cluster1_points.append(players_coordinates[i])
    if kmeans_3clusters[1][i] == 1:
        kmeans_3clusters_cluster2_points.append(players_coordinates[i])
    if kmeans_3clusters[1][i] == 2:
        kmeans_3clusters_cluster3_points.append(players_coordinates[i])

# plots centroids and clusters found for Kmeans (2 and 3 centroids) for the
# Finnish players' coordinates.
fig, (plt1, plt2) = pylab.subplots(1,2, sharey=True)
plt1.plot([point[0] for point in kmeans_2clusters_cluster1_points],
          [point[1] for point in kmeans_2clusters_cluster1_points], 'bo',
          label='Cluster 1')
plt1.plot([point[0] for point in kmeans_2clusters_cluster2_points],
          [point[1] for point in kmeans_2clusters_cluster2_points], 'ro',
          label='Cluster 2')
plt1.plot([centroid[0] for centroid in kmeans_2clusters[0]],
          [centroid[1] for centroid in kmeans_2clusters[0]],'ko',
          label='Centróides')
plt1.legend(fontsize=15)
plt1.set_xlabel(r'$Longitude$', fontsize=15)
plt1.set_ylabel(r'$Latitude$', fontsize=15)
plt2.plot([point[0] for point in kmeans_3clusters_cluster1_points],
          [point[1] for point in kmeans_3clusters_cluster1_points], 'bo',
          label='Cluster 1')
plt2.plot([point[0] for point in kmeans_3clusters_cluster2_points],
          [point[1] for point in kmeans_3clusters_cluster2_points], 'ro',
          label='Cluster 2')
plt2.plot([point[0] for point in kmeans_3clusters_cluster3_points],
          [point[1] for point in kmeans_3clusters_cluster3_points], 'go',
          label='Cluster 3')
plt2.plot([centroid[0] for centroid in kmeans_3clusters[0]],
          [centroid[1] for centroid in kmeans_3clusters[0]],'ko', label='Centróides')
plt2.legend(fontsize=15)
plt2.set_xlabel(r'$Longitude$', fontsize=15)
pylab.show()


################################################################################
############################## PARTE 2 #########################################

img = Image.open('Imagem1.jpg')
rgb_matrix = np.array(img)

img_pixels = []

for line in rgb_matrix:
    for pixel in line:
        img_pixels.append(pixel)

img_kmeans_result = kMeans(data=img_pixels,
                           k_max=3,
                           max_iterations=15)

pixel_group_image = [img_kmeans_result[1][x:x+164] for x in range(0, len(img_kmeans_result[1]), 164)]

black_cluster_image = [[0 for j in range(164)] for i in range(108)]
for i in range(108):
    for j in range(164):
        if pixel_group_image[i][j] == 0:
            black_cluster_image[i][j] = img_kmeans_result[0][0]
        else:
            black_cluster_image[i][j] = [255, 255, 255]

for i in range(len(black_cluster_image)):
    black_cluster_image[i] = np.array(black_cluster_image[i])

black_cluster_image = np.asarray(black_cluster_image)
black_cluster_image = Image.fromarray((black_cluster_image).astype('uint8'), mode='RGB')
black_cluster_image.show()

blue_cluster_image = [[0 for j in range(164)] for i in range(108)]
for i in range(108):
    for j in range(164):
        if pixel_group_image[i][j] == 1:
            blue_cluster_image[i][j] = img_kmeans_result[0][1]
        else:
            blue_cluster_image[i][j] = [255, 255, 255]

for i in range(len(blue_cluster_image)):
    blue_cluster_image[i] = np.array(blue_cluster_image[i])


blue_cluster_image = np.asarray(blue_cluster_image)
blue_cluster_image = Image.fromarray((blue_cluster_image).astype('uint8'), mode='RGB')
blue_cluster_image.show()

red_cluster_image = [[0 for j in range(164)] for i in range(108)]
for i in range(108):
    for j in range(164):
        if pixel_group_image[i][j] == 2:
            red_cluster_image[i][j] = img_kmeans_result[0][2]
        else:
            red_cluster_image[i][j] = [255, 255, 255]

for i in range(len(red_cluster_image)):
    red_cluster_image[i] = np.array(red_cluster_image[i])

red_cluster_image = np.asarray(red_cluster_image)
red_cluster_image = Image.fromarray((red_cluster_image).astype('uint8'), mode='RGB')
red_cluster_image.show()


# convert image to gray scale
gray_scale_matrix = []
gray_img_pixels = []
for line in rgb_matrix:
    line_gray_scale = []
    for pixel in line:
        R = int(pixel[0])
        G = int(pixel[1])
        B = int(pixel[2])
        I = (R+G+B)/(3*255)
        line_gray_scale.append(I)
        gray_img_pixels.append(I)
    gray_scale_matrix.append(line_gray_scale)

gray_scale_array = np.asarray(gray_scale_matrix)

gray_image = Image.fromarray((gray_scale_array*255).astype('uint8'), mode='L')
gray_image.show()


num_clusters = 3
gray_img_kmeans_result = kMeans(data=gray_img_pixels,
                           k_max=num_clusters,
                           max_iterations=15)

print(gray_img_kmeans_result[0])

# datapoints_per_cluster = [...,[COORDINATES OF POINTS OF CLUSTER_k],...]
datapoints_per_cluster = [[] for i in range(num_clusters)]
markers = ['r.', 'bx', 'k1', 'g^']

for k in range(num_clusters):
    for j in range(len(gray_img_pixels)):
        if gray_img_kmeans_result[1][j] == k:
            datapoints_per_cluster[k].append(gray_img_pixels[j])
    pylab.plot(datapoints_per_cluster[k], markers[k])
pylab.show()

pixel_group_image = [gray_img_kmeans_result[1][x:x+164] for x in range(0, len(gray_img_kmeans_result[1]), 164)]

image_names = ['graykmeans1.png', 'graykmeans2.png', 'graykmeans3.png']

for k in range(num_clusters):
    k_cluster_image = [[0 for j in range(164)] for i in range(108)]
    for i in range(108):
        for j in range(164):
            if pixel_group_image[i][j] == k:
                k_cluster_image[i][j] = np.array((0,0,0))
            else:
                k_cluster_image[i][j] = np.array((255, 255, 255))

    for i in range(len(k_cluster_image)):
        k_cluster_image[i] = np.array(k_cluster_image[i])

    k_cluster_image = np.asarray(k_cluster_image)
    k_cluster_image = Image.fromarray((k_cluster_image).astype('uint8'), mode='RGB')
    k_cluster_image.save(image_names[k])
    k_cluster_image.show()


red_cluster = []
blue_cluster = []
black_cluster = []

for i in range(len(img_pixels)):
    if img_kmeans_result[1][i] == 0:
        red_cluster.append(img_pixels[i])
    elif img_kmeans_result[1][i] == 1:
        blue_cluster.append(img_pixels[i])
    elif img_kmeans_result[1][i] == 2:
        black_cluster.append(img_pixels[i])

# plot clusters found by K-Means in RGB space
fig = pylab.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot([point[0] for point in red_cluster],
           [point[1] for point in red_cluster],
           [point[2] for point in red_cluster], 'r.')
ax.plot([point[0] for point in blue_cluster],
           [point[1] for point in blue_cluster],
           [point[2] for point in blue_cluster], 'b1')
ax.plot([point[0] for point in black_cluster],
           [point[1] for point in black_cluster],
           [point[2] for point in black_cluster], 'kx')
ax.set_xlabel('R')
ax.set_ylabel('G')
ax.set_zlabel('B')
pylab.show()



################################################################################
############################## PARTE 3 #########################################
def average_value(array):
    sum = 0
    for value in array:
        sum += value
    average = sum/len(array)
    return average


def standard_deviation(array):
    average = average_value(array)
    desvio_sum = 0
    for value in array:
        desvio_sum += (value - average)**2
    variancia = desvio_sum/len(array)
    return variancia**0.5


def normalize(data):
    N = len(data)
    dim = len(data[0])

    normalized_data = [[0 for d in range(dim)] for n in range(N)]

    for d in range(dim):
        mean = average_value([point[d] for point in data])
        sd = standard_deviation([point[d] for point in data])

        for i in range(N):
            normalized_coord = (data[i][d] - mean) / sd
            normalized_data[i][d] = normalized_coord

    return normalized_data


def dot_product(v1, v2):
    D = len(v1)
    sum = 0
    for d in range(D):
        sum += v1[d]*v2[d]
    return sum


def PCA(data, k_max, tolerance):
    N = len(data)
    dim = len(data[0])

    normalized_data = normalize(data)

    pc = []

    # find principal component 1 and principal component 2
    for j in range(2):

        p = [random.randint(-3,3) for d in range(dim)]
        r = [x/(dot_product(p,p)**0.5) for x in p]

        error = 0

        k = 0
        # calculates principal component
        while (k < k_max) and (error < tolerance):
            s = [0 for d in range(dim)]

            for i in range(1,N):
                dot = dot_product(normalized_data[i], r)
                temp = [x*dot for x in normalized_data[i]]
                s = [s[i]+temp[i] for i in range(len(s))]
            prior_r = r
            dot = dot_product(s,s)
            r = [x/(dot**0.5) for x in s]
            error = 0

            for i in range(len(r)):
                if (r[i] - prior_r[i]) > error:
                    error = r[i]-prior_r[i]

            k += 1

        pc.append(r)

        for i in range(len(normalized_data)):
            temp = dot_product(normalized_data[i], r)/dot_product(r, r)
            proj = [x*temp for x in r]

            # remove projection of datapoint to found principal component
            for x in range(len(normalized_data[i])):
                normalized_data[i][x] -= proj[x]

    return pc




with open('semiconductors.csv', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]

    semiconductors = []
    semiconductors_parameters = []

    for line in readlines:
        if readlines.index(line) != 0:
            temp = line.split(',')

            for i in range(1,len(temp)):
                temp[i] = float(temp[i])

            semiconductors.append(temp)
            semiconductors_parameters.append(temp[1:])

semiconductors_kmeans_result = kMeans(data=semiconductors_parameters, k_max=3, max_iterations=15)
for i in range(len(semiconductors_kmeans_result[1])):
    print("{}\t{}".format(i+2, semiconductors_kmeans_result[1][i]))


pca_result = PCA(data=semiconductors_parameters,
                 k_max=15,
                 tolerance=4)

print("Vector próprio PC1:", pca_result[0])
print("Vector próprio PC2:", pca_result[1])


markers = ['ro', 'bo', 'ko', 'go', 'co', 'yo']
labels = [r'Atomic no',r'Melting point',r'VE',r'radii',r'EN',r'lattice const. (ang)']
for i in range(6):
    pylab.plot(pca_result[0][i], pca_result[1][i], markers[i], label=labels[i])

pylab.ylabel("Principal Component 2", fontsize=20)
pylab.xlabel("Principal Component 1", fontsize=20)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.legend(fontsize=20)
pylab.show()


norm1 = 0
for x in pca_result[0]:
    norm1 += x**2
norm2 = 0
for x in pca_result[1]:
    norm2 += x**2

projection_of_data_in_pc = []
for line in semiconductors_parameters:
    dot1 = dot_product(line, pca_result[0])
    dot2 = dot_product(line, pca_result[1])
    proj1_scale = (dot1/norm1)
    proj2_scale = (dot2/norm2)
    projection_of_data_in_pc.append([proj1_scale, proj2_scale])


pylab.plot([line[0] for line in projection_of_data_in_pc],
           [line[1] for line in projection_of_data_in_pc], 'b.')
pylab.ylabel(r"$\frac{u \cdot pc_2}{u \cdot u}$", fontsize=20)
pylab.xlabel(r"$\frac{u \cdot pc_1}{u \cdot u}$", fontsize=20)
pylab.xticks(fontsize=20)
pylab.yticks(fontsize=20)
pylab.show()


# appends new result of pca method to dataframe
with open('resultats_pca.csv', 'a') as file:
    file.write("{}\t{}\n".format(pca_result[0],pca_result[1]))
