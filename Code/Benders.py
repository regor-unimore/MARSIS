import time
import argparse
import gurobipy as gp
import numpy as np
from gurobipy import GRB

# Define INSTANCE
#parser = argparse.ArgumentParser(description="Code to solve the Benders decomposition")
#parser.add_argument('instance_name', type=str, help="Name of the instance")
#args = parser.parse_args()
#name = args.instance_name

name = 'M_10_134_4'

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

# Day of the observations
with open("day.txt", "r") as file:
    day = [line for line in file.readlines()]

# Class of the quality of PIs in observations for each scenario s, weighted through set V
with open(f"{name}/{name}(4).txt", "r") as file:
    cPIinO= [[V[int(i)] for i in line.strip().split("\t")]
             if line != "" and line != "\n" else [] for line in file.readlines()]

def real_obj_f(selected):
    q = np.array(cPIAV)
    for o in range(nO):
        if selected[o] == 1:
            for j in range(len(PIinO[o])):
                q[PIinO[o][j]] = max(q[PIinO[o][j]], cPIinO[o][j])
    real_obj = sum((q[j] - cPIAV[j]) for j in range(nPI))
    return real_obj, q

def lazy(model, where):
    if where == GRB.Callback.MIPSOL:
        # retrieves the solution
        x = model.cbGetSolution(model._x)
        ts = model.cbGetSolution(model._t)
        real_obj, q = real_obj_f(x)  # compute the real objective function and the maximum quality with which each point is covered
        if ts != real_obj:
            adv = [sum((cPIinO[i][j] - q[PIinO[i][j]]) for j in range(len(PIinO[i])) if cPIinO[i][j] > q[PIinO[i][j]])
                   if x[i] == 0 else 0 for i in range(nO)]  # advantage I would get from each new selected observation
            model.cbLazy(model._t <= real_obj + sum(adv[i] * model._x[i] for i in range(nO)))  # Cut B0f

############################## MASTER PROBLEM - MODEL #######################################
t1 = time.time()  # start time
model = gp.Model()
model.setParam("TimeLimit", 3600)
model.setParam(GRB.Param.Threads, 1)
model.setParam(GRB.Param.Method, 2)
model.Params.lazyConstraints = 1

# variables
xo = model.addVars(range(nO), lb=0, ub=1, vtype=GRB.INTEGER, name='X')
theta = model.addVar(lb=0, vtype=GRB.CONTINUOUS, name='T')
model.update()

# objective function
model.setObjective(theta, sense=gp.GRB.MAXIMIZE)

# constraints
nbPicPerDay = {day[o]: 0 for o in range(nO)}
for o in range(nO): nbPicPerDay[day[o]] += xo[o]
model.addConstrs(nbPicPerDay[i] <= N for i in nbPicPerDay.keys())

ub = [sum((cPIinO[o][j] - cPIAV[PIinO[o][j]]) for j in range(len(PIinO[o])) if cPIinO[o][j] > cPIAV[PIinO[o][j]])
      for o in range(nO)]
UB = sum(xo[o] * ub[o] for o in range(nO))
model.addConstr(theta <= UB)

t2 = time.time()  # finish construction time
model._t = theta
model._x = xo
model.optimize(lazy)
t3 = time.time()  # finish execution time

visited = [0]*nPI
selected = [0]*nO
for o in range(nO):
    selected[o] = xo[o].x
    if selected[o] == 1:
        for p in PIinO[o]:
            visited[p] = 1
coverture = sum(visited)
obj, q = real_obj_f(selected)

print(f'Is the solution optimal? {model.status == gp.GRB.OPTIMAL}\n')
print(f'model objective function = {theta.x}')
print(f'real objective function = {obj}')
print(f'number of PI visited = {coverture}')

f = open(f"Benders_{name}.txt", "w")
f.write(f'objective function = {obj}\n')
f.write(f'number of PI visited = {coverture}\n')
f.write(f'total time = {t3-t1}\n')
f.write(f'construction time = {t2-t1}\n')
f.write(f'execution time = {t3-t2}\n')
f.write(f'UB = {model.ObjBound}\n')
f.write(f'Gap = {model.MIPGap}\n')
f.write(f'explored nodes= {model.NodeCount}\n')
for i in selected:
    f.write("%d\n" % i)
f.close()
