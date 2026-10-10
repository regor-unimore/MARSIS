import time
import gurobipy as gp
from gurobipy import GRB
import argparse

# Define INSTANCE
#parser = argparse.ArgumentParser(description="Code to solve the ILP model")
#parser.add_argument('instance_name', type=str, help="Name of the instance")
#args = parser.parse_args()
#name = args.instance_name

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
with open(f"Opt-instance/{name}/{name}(1).txt", 'r') as file:
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

############################# MODEL ################################################

t1 = time.time()  # start construction time
model = gp.Model()
model.setParam("TimeLimit", 3600)
model.setParam(GRB.Param.Threads, 1)
model.setParam(GRB.Param.Method, 2)
model.setParam("MIPGap", 1e-6)

z = 0
PC = {(PIinO[o][i], cPIinO[o][i]): 1 for o in range(nO) for i in range(len(PIinO[o]))}  # feasible (p, c) couples, c = class of quality q
NumPC = {(p, c): 0 for (p, c) in PC.keys()}  # number of time PI is covered with class c
NumObsPerDay = {day[o]: 0 for o in range(nO)} # number of observations selected for each day

# Variables
xo = model.addVars(range(nO), lb=0, ub=1, vtype=GRB.INTEGER, name='X')
y = model.addVars(PC.keys(), lb=0, ub=1, vtype=GRB.INTEGER, name='y')

model.update()

for o in range(nO):
    NumObsPerDay[day[o]] += xo[o]
    for i in range(len(PIinO[o])):
        NumPC[(PIinO[o][i], cPIinO[o][i])] += xo[o]

# Objective function
for (p, c) in PC.keys():
    if c > cPIAV[p]:
        z += y[(p, c)] * (c - cPIAV[p]) * wPI[p]

model.setObjective(z, sense=gp.GRB.MAXIMIZE)

# Constraints
model.addConstrs(y[(p, c)] <= NumPC[(p, c)] for (p, c) in PC.keys())
model.addConstrs(y.sum(p, '*') <= 1 for p in range(nPI))
model.addConstrs(NumObsPerDay[d] <= N for d in NumObsPerDay.keys())

t2 = time.time()  # finish construction time
model.optimize()
t3 = time.time()  # finish execution time

# Save solution and compute coverage
visited = [0]*nPI
selected = [0]*nO
for o in range(nO):
    selected[o] = xo[o].x
    if selected[o] == 1:
        for p in PIinO[o]:
            visited[p] = 1
coverage = sum(visited)

# Print output
print(f'Is the solution optimal? {model.status == gp.GRB.OPTIMAL}\n')
print(f'objective function = {z.getValue()}')
print(f'number of PI visited = {coverage} \n')
print(f'total time = {t3 - t1}')
print(f'construction time = {t2 - t1}')
print(f'execution time = {t3 - t2}')
print(f'UB = {model.ObjBound}')
print(f'Gap = {model.MIPGap}')

#  Save output
f = open(f"ModelSolution_{name}.txt", "w")
f.write(f'Gra: {G}, L: {L}, NC: {NC} \n')
f.write(f'objective function = {z.getValue()} \n')
f.write(f'number of PI visited = {coverage} \n')
f.write(f'total time = {t3 - t1} \n')
f.write(f'construction time = {t2 - t1} \n')
f.write(f'execution time = {t3 - t2} \n')
f.write(f'UB = {model.ObjBound} \n')
f.write(f'Gap = {model.MIPGap} \n')
f.write(f'explored nodes = {model.NodeCount} \n')
for i in selected:
    f.write("%d\n" % i)
f.close()




