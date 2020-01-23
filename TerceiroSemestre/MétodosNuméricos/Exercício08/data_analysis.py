# Author: Ernesto González
# Date: 28/11/2019

# source elections2016.txt: https://data.world/garyhoov/2016-pres-election-by-county/workspace/file?filename=2016+Presidential+Election+Analysis.xlsx
# source averageglobaltemperature.txt: http://climexp.knmi.nl/getindices.cgi?STATION=Berkeley_land_temperature&TYPE=i&WMO=BerkeleyData/t2m_land_best&id=$id

import numpy as np
import math
import pylab
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
import csv
import datetime
import re


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


with open('allwithnamestate.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]
    data_frame_2008 = [[0,0,0,0] for i in range(len(readlines))]

    # separates first line in a list of the column values
    readlines[0] = readlines[0].split('\t')
    readlines[0] = ' '.join(readlines[0])
    readlines[0] = readlines[0].split(' ')

    for line in readlines:
        if readlines.index(line) != 0:
            # separates all other lines in a list of the column values
            readlines[readlines.index(line)] = line.split(' ')

    for i in range(len(readlines)):
        r_line = readlines[i]
        df_line = data_frame_2008[i]
        df_line[0], df_line[1] = int(r_line[0]), int(r_line[1])
        df_line[3] = r_line[-1]

        if len(r_line) > 4:
            if '' in r_line:
                r_line.remove('')
            r_line.remove(r_line[0])
            r_line.remove(r_line[0])
            r_line.remove(r_line[-1])
            county = ' '.join(r_line)
            df_line[2] = county
        else:
            df_line[2] = r_line[2]


def county_results(df, first_candidate, second_candidate):
    """ Studies results of votes for the American presidential elections in each county.
        Arguments:
                - df is a list of of the results in each county. The result of each
                county is a list [VOTES_FIRST_CANDIDATE, VOTES_SECOND_CANDIDATE, COUNTY, STATE];
                - first_candidate is the name of the first candidate;
                - second_candidate is the name of the second candidate.
        Returns:
                - List with mean, median and variance of votes per state for each candidate
                in the form [STATE_NAME,
                            AVERAGE_FIRST_CANDIDATE, MEDIAN_FIRST_CANDIDATE, VARIANCE_FIRST_CANDIDATE,
                            AVERAGE_SECOND_CANDIDATE, MEDIAN_SECOND_CANDIDATE, VARIANCE_SECOND_CANDIDATE];
                - Plots histograms of relative difference of votes in each county, considering all counties
                in one histogram and only major counties (plus 20000 voters) in other histogram
    """
    states = []
    for line in df:
        if line[3] not in states:
            states.append(line[3])

    # vote_data_per_state = [...,[STATE_NAME,
    #                       AVERAGE_FIRST_CANDIDATE, MEDIAN_FIRST_CANDIDATE, VARIANCE_FIRST_CANDIDATE,
    #                       AVERAGE_SECOND_CANDIDATE, MEDIAN_SECOND_CANDIDATE, VARIANCE_SECOND_CANDIDATE],...]
    # taking counties as unit event
    vote_data_per_state = [[0,0,0,0,0,0,0] for i in range(len(states))]
    # vote_data_national = [NAT_AVERAGE_FIRST_CANDIDATE, NAT_MEDIAN_FIRST_CANDIDATE, NAT_VARIANCE_FIRST_CANDIDATE,
    #                       NAT_AVERAGE_SECOND_CANDIDATE, NAT_MEDIAN_SECOND_CANDIDATE, NAT_VARIANCE_SECOND_CANDIDATE]

    # taking states as unit event
    vote_data_national_1 = []

    # taking counties as unit event
    vote_data_national_2 = []

    # state_total_votes = [...,[STATE_NAME, VOTES_FIRST_CANDIDATE, VOTES_SECOND_CANDIDATE],...]
    state_total_votes = []

    vote_data_national_2.append(average_value([df[i][0] for i in range(len(df))]))
    vote_data_national_2.append(median_value([df[i][0] for i in range(len(df))]))
    vote_data_national_2.append(variance([df[i][0] for i in range(len(df))]))
    vote_data_national_2.append(average_value([df[i][1] for i in range(len(df))]))
    vote_data_national_2.append(median_value([df[i][1] for i in range(len(df))]))
    vote_data_national_2.append(variance([df[i][1] for i in range(len(df))]))

    for state in states:
        state_first_candidate_votes = []
        state_second_candidate_votes = []
        vote_data_per_state[states.index(state)][0] = state
        for line in df:
            county_first_candidate_votes = line[0]
            county_second_candidate_votes = line[1]

            # if county in state
            if line[-1] == state:
                state_first_candidate_votes.append(county_first_candidate_votes)
                state_second_candidate_votes.append(county_second_candidate_votes)

        state_first_candidate_total_votes = 0
        for county_votes in state_first_candidate_votes:
            state_first_candidate_total_votes += county_votes

        state_second_candidate_total_votes = 0
        for county_votes in state_second_candidate_votes:
            state_second_candidate_total_votes += county_votes

        state_total_votes.append([state, state_first_candidate_total_votes, state_second_candidate_total_votes])

        vote_data_per_state[states.index(state)][1] = average_value(state_first_candidate_votes)
        vote_data_per_state[states.index(state)][2] = median_value(state_first_candidate_votes)
        vote_data_per_state[states.index(state)][3] = variance(state_first_candidate_votes)
        vote_data_per_state[states.index(state)][4] = average_value(state_second_candidate_votes)
        vote_data_per_state[states.index(state)][5] = median_value(state_second_candidate_votes)
        vote_data_per_state[states.index(state)][6] = variance(state_second_candidate_votes)

    vote_data_national_1.append(average_value([state_total_votes[i][1] for i in range(len(state_total_votes))]))
    vote_data_national_1.append(median_value([state_total_votes[i][1] for i in range(len(state_total_votes))]))
    vote_data_national_1.append(variance([state_total_votes[i][1] for i in range(len(state_total_votes))]))
    vote_data_national_1.append(average_value([state_total_votes[i][2] for i in range(len(state_total_votes))]))
    vote_data_national_1.append(median_value([state_total_votes[i][2] for i in range(len(state_total_votes))]))
    vote_data_national_1.append(variance([state_total_votes[i][2] for i in range(len(state_total_votes))]))


    relative_difference_of_votes_all_counties = []
    relative_difference_of_votes_major_counties = []

    absolute_difference_of_votes_all_counties = []
    absolute_difference_of_votes_major_counties = []

    for county in df:
        # consider only counties with more than 20000 electors
        if county[0]+county[1] >=20000:
            relative_difference_of_votes_major_counties.append((county[0]-county[1])/(county[0]+county[1]))
            absolute_difference_of_votes_major_counties.append((county[0]-county[1]))
        # consider all counties
        relative_difference_of_votes_all_counties.append((county[0]-county[1])/(county[0]+county[1]))
        absolute_difference_of_votes_all_counties.append((county[0]-county[1]))

    relative_difference_of_votes_major_counties.sort()
    relative_difference_of_votes_all_counties.sort()

    # list of frequency of values in intervals from -1 to 1 with amplitude 0.1
    # first element in list is frequency at [-1,-0.9[, second is frequency at [-0.9, -0.8[ and so on
    rel_diff_intervals_major_counties = [0 for i in range(20)]
    rel_diff_intervals_all_counties = [0 for i in range(20)]

    for rel_diff in relative_difference_of_votes_major_counties:
            index_of_interval = int(math.floor((rel_diff+1)*10))
            rel_diff_intervals_major_counties[index_of_interval] +=1

    for rel_diff in relative_difference_of_votes_all_counties:
            index_of_interval = int(math.floor((rel_diff+1)*10))
            rel_diff_intervals_all_counties[index_of_interval] +=1

    hist_intervals_min = [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    # plots histogram
    fig, (plt1, plt2) = pylab.subplots(1,2)
    plt1.bar(hist_intervals_min, rel_diff_intervals_all_counties, width=0.1, align='edge',
              color=['blue','blue','blue','blue','blue','blue','blue','blue','blue',
              'blue','green','green','green','green','green','green','green','green','green','green'],
              edgecolor='black', linewidth=2)
    blue_patch = mpatches.Patch(color='blue', label=second_candidate)
    green_patch = mpatches.Patch(color='green', label=first_candidate)
    plt1.text(0.5,300, "(a)", fontsize=13)
    plt1.legend(handles=[blue_patch, green_patch], fontsize=12)
    plt1.set_xticks([-1,-0.5,0,0.5,1])
    plt1.set_xlabel(r'$v_{}-v_{}$'.format(first_candidate[0], second_candidate[0]), fontsize=15)
    plt1.set_ylabel(r'$N_c$', fontsize=15)
    plt1.xaxis.set_tick_params(labelsize=15)
    plt1.yaxis.set_tick_params(labelsize=15)
    plt2.bar(hist_intervals_min, rel_diff_intervals_major_counties, width=0.1, align='edge',
              color=['blue','blue','blue','blue','blue','blue','blue','blue','blue',
              'blue','green','green','green','green','green','green','green','green','green','green'],
              edgecolor='black', linewidth=2)
    blue_patch = mpatches.Patch(color='blue')
    green_patch = mpatches.Patch(color='green')
    plt2.text(-0.8,130, "(b)", fontsize=13)
    plt2.set_xticks([-1,-0.5,0,0.5,1])
    plt2.set_xlabel(r'$v_{}-v_{}$'.format(first_candidate[0], second_candidate[0]), fontsize=15)
    plt2.xaxis.set_tick_params(labelsize=15)
    plt2.yaxis.set_tick_params(labelsize=15)
    pylab.show()

    # histogram of absolute county votes difference.
    n, bins, patches = pylab.hist(absolute_difference_of_votes_major_counties, 20,
            align='right',range=[-133000,130000],facecolor='blue',edgecolor='black', linewidth=2)
    pylab.xlabel(r'$v_{}-v_{}$'.format(first_candidate[0], second_candidate[0]), fontsize=15)
    pylab.ylabel(r'$N_c$', fontsize=15)
    pylab.axis([-133000, 130000, 0, 550])
    pylab.show()

    return vote_data_per_state, vote_data_national_1, vote_data_national_2


with open('elections2016.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]

    # removes header line
    readlines.remove(readlines[0])

    data_frame_2016 = [[0,0,0,0] for i in range(len(readlines))]

    # separates all lines in a list of the column values
    for line in readlines:
        readlines[readlines.index(line)] = line.split('\t')

    for i in range(len(readlines)):
        r_line = readlines[i]
        df_line = data_frame_2016[i]

        df_line[2], df_line[3] = r_line[1], r_line[0]

        # converts values to integers
        r_line[4] = r_line[4].strip(' ')
        r_line[4] = r_line[4].split(' ')
        df_line[0] = int(''.join(r_line[4]))
        r_line[5] = r_line[5].strip(' ')
        r_line[5] = r_line[5].split(' ')
        df_line[1] = int(''.join(r_line[5]))



county_results_2008 = county_results(data_frame_2008, 'Obama', 'McCain')
county_results_2016 = county_results(data_frame_2016, 'Clinton', 'Trump')

print(county_results_2008[2])
print(county_results_2016[2])


# ############################### PARTE 2 ########################################

with open('spots.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]
    data_frame_spots = []

    for line in readlines:
        temp = line.split()
        temp[0], temp[1], temp[2] = int(temp[0]), int(temp[1]), float(temp[2])
        data_frame_spots.append([datetime.date(temp[0], temp[1], 1), temp[2]])

# plots solar spots monthly Wolf's number time series
fig, ax = pylab.subplots()
ax.plot([data_frame_spots[i][0] for i in range(len(data_frame_spots))],
        [data_frame_spots[i][1] for i in range(len(data_frame_spots))])
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$W$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()


def autocorrelation(array):
    n = len(array)
    autocorrelation = []
    y = array
    y_mean = average_value(array)

    for k in range(n):
        sum = 0
        for t in range(n-k):
            sum += (y[t]-y_mean)*(y[t+k]-y_mean)
        numerator = sum/(n-k)

        sum = 0
        for t in range(n):
            sum += (y[t]-y_mean)**2
        denominator = sum/n

        r_k = numerator / denominator
        autocorrelation.append(r_k)

    return autocorrelation


r = autocorrelation([line[1] for line in data_frame_spots])

# plots solar spots monthly Wolf's number autocorrelation time series
fig, ax = pylab.subplots()
ax.plot([i for i in range(len(data_frame_spots))],
        r)
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$k$", fontsize=15)
pylab.ylabel("$r_k$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

# seasonal_difference list with values of w(t) = f(t) - f(t-T), where f(t) is
# the value of wolf's number at time t and T is the period (given in months)
# of the sasonality observed
seasonal_difference = []
for i in range(128):
    seasonal_difference.append(data_frame_spots[i][1])
for i in range(128,len(data_frame_spots)):
    seasonal_difference.append(data_frame_spots[i][1]-data_frame_spots[i-128][1])

fig, ax = pylab.subplots()
ax.plot([data_frame_spots[i][0] for i in range(len(data_frame_spots))], seasonal_difference)
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$W(t)- W(t-T)$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()


with open('averageglobaltemperature.txt', 'r') as infile:
    readlines = [line.rstrip('\n') for line in infile]

    data_frame_avg_temperature = readlines[7:]

    for line in data_frame_avg_temperature:
        temp = line.split(' ')
        temp = [x for x in temp if x != '']
        for elem in temp:
            temp[temp.index(elem)] = float(elem)
        data_frame_avg_temperature[data_frame_avg_temperature.index(line)] = temp

# converts two first columns to a single datetime.date format column
for line in data_frame_avg_temperature:
    temp = [datetime.date(int(line[0]),int(line[1]),1), line[2]]
    data_frame_avg_temperature[data_frame_avg_temperature.index(line)] = temp

# plots average monthly temperature time series
fig, ax = pylab.subplots()
ax.plot([data_frame_avg_temperature[i][0] for i in range(len(data_frame_avg_temperature))],
        [data_frame_avg_temperature[i][1] for i in range(len(data_frame_avg_temperature))])
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$T(\circ C)$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()


def pearson_correlation_coeff(X, Y):
    return covariance(X,Y)/((variance(X)**0.5)*(variance(Y)**0.5))

# reduce data_frame_avg_temperature and data_frame_spots to same date interval
data_frame_spots = data_frame_spots[12:]
data_frame_avg_temperature = data_frame_avg_temperature[:2807]

# data_frame_avg_temperature lacking values for 1/12/1751 so split data frames into
# two continuous time intervals (there will be a discontinuity point in correlation plot)
first_interval_data_frame_spots = data_frame_spots[:22]
first_interval_data_frame_avg_temperature = data_frame_avg_temperature[:22]
second_interval_data_frame_spots = data_frame_spots[24:]
second_interval_data_frame_avg_temperature = data_frame_avg_temperature[23:]


# plots wolf's number variation with temperature each month
first_int_months = len(first_interval_data_frame_spots)
second_int_months = len(second_interval_data_frame_spots)
pylab.plot([first_interval_data_frame_spots[i][1] for i in range(first_int_months)],
    [first_interval_data_frame_avg_temperature[i][1] for i in range(first_int_months)], 'bo')
pylab.plot([second_interval_data_frame_spots[i][1] for i in range(second_int_months)],
    [second_interval_data_frame_avg_temperature[i][1]for i in range(second_int_months)], 'bo')
pylab.xlabel("$T$", fontsize=15)
pylab.ylabel("$W$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

# join the two intervals' data into only one list, for each variable
first_interval_data_frame_spots.extend(second_interval_data_frame_spots)
first_interval_data_frame_avg_temperature.extend(second_interval_data_frame_avg_temperature)

print("Coeficiente de correlação de Pearson entre a medição mensal do número de Wolf e a temperatura média mensal:",
    pearson_correlation_coeff(
    [first_interval_data_frame_spots[i][1] for i in range(2806)],
    [first_interval_data_frame_avg_temperature[i][1] for i in range(2806)]
    ))


############################# PARTE 3 ##########################################
with open('earthquakedata.txt', 'r') as infile:
    readlines = infile.read()

    readlines = readlines.split("<")

    data_frame_earthquake = []

    for line in readlines:
        if line.startswith("|Entity"):
            readlines.remove(line)

    for line in readlines:
        date = re.search(r"DateObject\[\{(.*?)},",line)
        if date != None:
            date = date.group(1).strip(".").split(",")
            for elem in date:
                date[date.index(elem)] = int(elem.strip(" "))
            date = datetime.datetime(date[0], date[1], date[2], date[3], date[4], date[5])

            magnitude =  re.search(r"\"Magnitude\" -> (.*?),",line)
            if magnitude != None:
                magnitude = float(magnitude.group(1))

                data_frame_earthquake.append([date, magnitude])

# plots earthquake intensity's time series
fig, ax = pylab.subplots()
ax.plot([data_frame_earthquake[i][0] for i in range(len(data_frame_earthquake))],
        [data_frame_earthquake[i][1] for i in range(len(data_frame_earthquake))], '--bo')
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$I$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

# plots earthquake intensity's time series, without first earthquake
fig, ax = pylab.subplots()
ax.plot([data_frame_earthquake[i][0] for i in range(1,len(data_frame_earthquake))],
        [data_frame_earthquake[i][1] for i in range(1,len(data_frame_earthquake))], '--bo')
ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
fig.autofmt_xdate()
pylab.xlabel("$t$", fontsize=15)
pylab.ylabel("$I$", fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.show()

max_intensity = 0
index_max_intensity = 0
for earthquake in data_frame_earthquake:
    if earthquake[1] > max_intensity:
        max_intensity = earthquake[1]
        index_max_intensity = data_frame_earthquake.index(earthquake)

print(""""O terramoto com maior amplitude ocorreu em {}
    e teve amplitude {} na escala de Richter.""".format(
        data_frame_earthquake[index_max_intensity][0], data_frame_earthquake[index_max_intensity][1]))


# histogram of earthquakes
n, bins, patches = pylab.hist([data_frame_earthquake[i][1] for i in range(1,len(data_frame_earthquake))],
7, align='mid',range=[0,7],facecolor='blue',edgecolor='black', linewidth=2)
pylab.xlabel(r'$I$', fontsize=15)
pylab.ylabel(r'$F_I$', fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.axis([0, 7, 0, 70])
pylab.show()

a_1 = 3
b_1 = (np.log10(115)-3)/(-2.5)
M_1 = np.linspace(0, 7, 100)
N_1 = 10**(a_1-b_1*M_1)

a_2 = 3
b_2 = (np.log10(120)-3)/(-1.8)
M_2 = np.linspace(0, 7, 100)
N_2 = 10**(a_2-b_2*M_2)

a_3 = 3
b_3 = (np.log10(200)-3)/(-1)
M_3 = np.linspace(0, 7, 100)
N_3 = 10**(a_3-b_3*M_3)

# histogram of earthquakes N that have at least intensity I
n, bins, patches = pylab.hist([data_frame_earthquake[i][1] for i in range(1,len(data_frame_earthquake))],
7, align='mid',range=[0,7], cumulative=-1, facecolor='blue',edgecolor='black', linewidth=2)
pylab.plot(M_1, N_1, label="$a=3,\,b=0.38$")
pylab.plot(M_2, N_2, label="$a=3,\,b=0.51$")
pylab.plot(M_3, N_3, label="$a=3,\,b=0.70$")
pylab.legend(fontsize=15)
pylab.xlabel(r'$I$', fontsize=15)
pylab.ylabel(r'$N$', fontsize=15)
pylab.xticks(fontsize=15)
pylab.yticks(fontsize=15)
pylab.axis([0, 7, 0, 180])
pylab.show()
