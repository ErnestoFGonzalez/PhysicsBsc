import math



## {10.5, 9.3, 11.4, 10.9, 13.0, 8.4, 9.2, 8.9, 10.3, 11.2, 12.1, 8.4, 9.2, 9.9, 10.1}
array_ = [ 10.5, 9.3, 11.4, 10.9, 13.0, 8.4, 9.2, 8.9, 10.3, 11.2, 12.1, 8.4, 9.2, 9.9, 10.1]

# if array_ has even number of values
if len(array_) % 2 == 0:
    for i in range(0,len(array_),2):
        print(str(array_[i])+'\t'+str(array_[i+1]))
# if array has odd number of values
else:
    for i in range(0,(len(array_)-1),2):
        print(str(array_[i])+'\t'+str(array_[i+1]))
    print(str(array_[-1])+'\t'+' -')


def average_value(array):
    sum = 0
    for value in array:
        sum += value
    average_value = sum/len(array)
    return average_value


def variancia(array):
    average = average_value(array)
    desvio_sum = 0
    for value in array:
        desvio_sum += (value - average)**2
    variancia = desvio_sum/len(array)
    return variancia


def desvio_padrao(array):
    variancia_ = variancia(array)
    desvio_padrao = math.sqrt(variancia_)
    return desvio_padrao



print("Average:",average_value(array_))
print("Variance:",variancia(array_))
print("Mean Deviation:",desvio_padrao(array_))


with open('array.txt', 'w') as file:
    # if array_ has even number of values
    if len(array_) % 2 == 0:
        for i in range(0,len(array_),2):
            print("{:f}\t\t{:f}".format(array_[i],array_[i+1]), file=file)
    # if array has odd number of values
    else:
        for i in range(0,(len(array_)-1),2):
            print("{:.2f}\t\t{:.2f}".format(array_[i],array_[i+1]), file=file)
        print("{:.2f}\t\t{}".format(array_[-1],' -'), file=file)


with open('array.txt', 'r') as file:
    value_columns = []
    new_array = []
    lines = file.readlines()
    for line in lines:
        line = line.strip('\n')
        line = line.split('\t\t')
        value_columns.append(line)
    for line in value_columns:
        for value in line:
            if value != ' -':
                new_array.append(float(value))
    print(new_array)
    file.close()
