import time
import gurobipy as gp
from gurobipy import GRB
import argparse

# Define INSTANCE
parser = argparse.ArgumentParser(description="Code to solve the ILP model")
parser.add_argument('instance_name', type=str, help="Name of the instance")
args = parser.parse_args()
name = args.instance_name

name = 'M_10_268_4'

name_char = name.split('_')

# Parameters
G = name_char[1]    # granularity
L = name_char[2]    # length
NC = name_char[3]   # number of quality classes
gra = {'01': 1000, '02': 2000, '03': 3000, '05': 5000, '10': 10000}
nO = 113651         # number of observations

# Default settings
N = 1                   # number of observation per day
V = [0, 0.25, 0.5, 1]   # class quality weights
W = [1, 1]              # [1/beta, 1], score for AOI

#  Information on PIs
cPIAV = []
wPI = []
with open(f"{name}/{name}(1).txt", 'r') as file:
    lines = file.readlines()
    nPI = len(lines)        # number of PI
    for line in lines:
        line = list(line.strip().split(' '))
        cPIAV.append(V[int(line[5])])
        wPI.append(W[int(line[6])])

#  PI in observations
with open(f"{name}/{name}(2).txt", "r") as file:
    PIinO = [[int(i) for i in line.strip().split("\t")]
            if line != "" and line != "\n" else [] for line in file.readlines()]

# Class of the quality of PIs in observations, weighted through set V
with open(f"{name}/{name}(4).txt", "r") as file:
    cPIinO = [[V[int(i)] for i in line.strip().split("\t")]
                    if line != "" and line != "\n" else [] for line in file.readlines()]

# Day of the observations
with open("day.txt", "r") as file:
    day = [line for line in file.readlines()]

# Heuristic solution used to initialize the x variables
with open(f"H2_{name}.txt", "r") as file:
    obj = float(file.readline().split("=")[1])
    file.readline()
    file.readline()
    file.readline()
    selected_h2 = [int(file.readline()) for i in range(nO)]

# Ranges of observations belonging to the same day
ranges = {}
beg = 0
while beg < nO:
    end = beg
    while end+1 < nO and day[end+1] == day[beg]:
        end = end+1
    ranges[day[beg]] = range(beg, end+1)
    beg = end+1

def retrieve_solution(xo, bestobj, bestsol, bestvis):
    qualitytemp = cPIAV[:]
    visited = [0] * nPI

    for i in range(nO):
        if xo[i] == 1:
            for j in range(len(PIinO[i])):
                visited[PIinO[i][j]] = 1
                qualitytemp[PIinO[i][j]] = max(qualitytemp[PIinO[i][j]], cPIinO[i][j])
    ob = 0
    for j in range(nPI):
        ob += (qualitytemp[j] - cPIAV[j]) * wPI[j]
    if ob > bestobj:
        bestobj = ob
        bestsol = list(xo.values())
        bestvis = sum(visited)

    return bestobj, bestsol, bestvis

def solution(lambdaM, bestobj, bestsol, bestvis):
    qualitytemp = cPIAV[:]
    visited = [0] * nPI
    xo = {o: 0 for o in range(nO)}

    for _ in range(N):
        for i in ranges.keys():
            idx = -1
            maxVal = -1.0
            for j in ranges[i]:
                curVal = 0.0
                for k in range(len(PIinO[j])):
                    curVal += lambdaM[PIinO[j][k], cPIinO[j][k]]
                if curVal > maxVal:
                    maxVal = curVal
                    idx = j
            if idx >= 0:
                xo[idx] = 1
                for k in range(len(PIinO[idx])):
                    visited[PIinO[idx][k]] = 1
                    qualitytemp[PIinO[idx][k]] = max(qualitytemp[PIinO[idx][k]], cPIinO[idx][k])

    ob = 0
    for j in range(nPI):
        ob += (qualitytemp[j] - cPIAV[j]) * wPI[j]

    if ob > bestobj:
        bestobj = ob
        bestsol = list(xo.values())
        bestvis = sum(visited)

    return bestobj, bestsol, bestvis

##################################### MAIN ###################################################
# parameters
PC = {(PIinO[o][i], cPIinO[o][i]): 1 for o in range(nO) for i in range(len(PIinO[o]))}  # feasible (p, c) couples, c = class of quality q
lambdaM = {(i, j): 1 for (i, j) in PC.keys()}
delta = 2.0
bestbound = 999999999.0
LB = obj
counter = 0
bestsol = [0]*nO
bestobj = 0
bestvis = 0

t1 = time.time()  # start time
while delta > 0.5:

    bound = 0.0
    y = {(p, c): 0 for (p, c) in PC.keys()}  # 0-1
    NumPC = {(p, c): 0 for (p, c) in PC.keys()}  # number of time PI is covered with class c
    xo = {o: 0 for o in range(nO)}  # 0-1

    # lagrangian solution construction
    for _ in range(N):
        for i in ranges.keys():  # i=day
            idx = -1
            maxVal = -1.0
            for j in ranges[i]:  # for each observation in the same day
                curVal = 0.0
                for k in range(len(PIinO[j])):  # compute the value of the observation based on lambda
                    curVal += lambdaM[PIinO[j][k], cPIinO[j][k]]
                if curVal > maxVal:  # find the maximum value
                    maxVal = curVal
                    idx = j
            if idx >= 0:
                xo[idx] = 1  # select the observation with maximum value
                for k in range(len(PIinO[idx])):
                    y[PIinO[idx][k], cPIinO[idx][k]] += 1  # PI covered with quality c

    bestobj, bestsol, bestvis = retrieve_solution(xo, bestobj, bestsol, bestvis)  # retrieve real solution value (PI covered only once)

    # update lambda with subgradient algorithm
    denominator = 0.0
    for (p, c) in PC.keys():
        if lambdaM[p, c] <= (c - cPIAV[p]) * wPI[p]:
            y[p, c] = 1
        bound += y[p, c] * (c - cPIAV[p]) * wPI[p] + lambdaM[p, c]*(NumPC[p, c] - y[p, c])
        denominator += (y[p, c] - NumPC[p, c]) * (y[p, c] - NumPC[p, c])

    if bestbound > bound:
        bestbound = bound
        counter = 0
    else:
        counter += 1

    theta = delta * (bound - LB) / denominator
    for (p, c) in PC.keys():
        lambdaM[p, c] = max(0.0, lambdaM[p, c] + theta * (y[p, c] - NumPC[p, c]))
    if counter > 20:
        delta /= 2
        counter = 0

bestobj, bestsol, bestvis = solution(lambdaM, bestobj, bestsol, bestvis) # build the last solution and compute the real value

t2 = time.time()  # end time

f = open(f"LagrangianSolution_{name}.txt", "w+")
f.write(f"objective function = {bestobj}\nvisited PI = {bestvis}\ntotal time = {t2-t1}\nUB = {bestbound}\n")
for i in range(nO):
    f.write("%d\n" % bestsol[i])
f.close()


