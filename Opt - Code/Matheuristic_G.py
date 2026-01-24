import math
import numpy as np
from statistics import mean
import time
import gurobipy as gp
from gurobipy import GRB
from shapely.geometry import Point, Polygon, MultiPolygon

# Define INSTANCE
#parser = argparse.ArgumentParser(description="Code to solve the ILP model")
#parser.add_argument('instance_name', type=str, help="Name of the instance")
#parser.add_argument('variant_name', type=str, help="R1, R2 or R3")
#args = parser.parse_args()
#name = args.instance_name
#variant = args.variant_name

name = 'M_10_134_4'
variant = 'G1'

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

# test parameters
seeds_G1 = [5, 85, 582, 133, 3371, 8759, 17803, 24093, 69350, 81604, 59, 855, 5482, 1233, 39371, 859, 1803, 94093, 6930, 8104]
seeds_G2 = [9, 54, 782, 533, 8371, 2759, 27803, 44093, 89350, 11604, 97, 504, 7282, 33, 80371, 24759, 7803, 4093, 8950, 1604]
seeds_G3 = [5, 85, 582, 133, 3371, 8759, 17803, 24093, 69350, 81604, 509, 835, 5782, 1833, 371, 80759, 86803, 2403, 6950, 8604]
ttype = {'01': 'VLR', '02': 'LR', '03': 'MR', '05': 'HR', '10': 'VHR'}
limits = {'VHR': [160, 80, 40, 20], 'HR': [80, 40, 20, 10], 'MR': [40, 20, 10, 5], 'LR': [20, 10, 5, 2], 'VLR': [10, 5, 2, 1]}

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
    PIinO = [[int(i) for i in line.strip().split("\t")] if line != "" and line != "\n" else []
             for line in file.readlines()]

# Class of the quality of PIs in observations, weighted through set V
with open(f"{name}/{name}(4).txt", "r") as file:
    cPIinO = [[V[int(i)] for i in line.strip().split("\t")] if line != "" and line != "\n" else []
              for line in file.readlines()]

# Day of the observations
with open(f"day.txt", "r") as file:
    day = [line for line in file.readlines()]

# Orbit of the observations
with open(f"orbit.txt", "r") as file:
    orbit = [line for line in file.readlines()]

# Heuristic solution used to initialize the x variables
with open(f"H2_{name}.txt", "r") as file:
    obj = float(file.readline().split("=")[1])
    file.readline()
    file.readline()
    file.readline()
    selected_h2 = [int(file.readline()) for i in range(nO)]

with open("x.dat", "r") as file:
    x = [float(line) for line in file.readlines()]

with open("y.dat", "r") as file:
    y = [float(line) for line in file.readlines()]

# Associate each observation starting point to an angle
round_angle = 2 * math.pi
angle = [0] * nO
for o in range(nO):
    angle[o] = math.atan2(y[o], x[o])
    if angle[o] < 0:
        angle[o] = round_angle + angle[o]

# ######################## MODEL ####################################################
model = gp.Model()
model.setParam("TimeLimit", 360)
model.setParam(GRB.Param.Threads, 1)
model.setParam(GRB.Param.Method, 2)

z = 0
PC = {(PIinO[o][i], cPIinO[o][i]): 1 for o in range(nO) for i in
      range(len(PIinO[o]))}  # feasible (p, c) couples, c = class of quality q
NumPC = {(p, c): 0 for (p, c) in PC.keys()}  # number of time PI is covered with class c
NumObsPerDay = {day[o]: 0 for o in range(nO)}  # number of observations selected for each day

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

# ######################## MODEL EXECUTION ###############################

# parameters
n_iter = 20
step = 0.10
coeff_0 = 0.15
n_slices = 5

obj_values = []
coverture_values = []
time_values = []
iter_values = []

for k in range(n_iter):

    t0 = time.time()  # starting time
    tf = 0  # total execution time
    tmodel = 0  # total model time

    iter_w_imp = 0
    iter_count = 0
    ilim = 0
    coeff = coeff_0
    if variant == 'G3':  # width of the geographical slice to fix
        delta = coeff * round_angle / n_slices
    else:
        delta = coeff * round_angle
    alpha = [np.random.random() * (round_angle - delta)]  # random rotation for variant G2

    iter_results = []
    selected = selected_h2  # initialize the solution
    if variant == 'G1': np.random.seed(seeds_G1[k])
    if variant == 'G2': np.random.seed(seeds_G2[k])
    if variant == 'G3': np.random.seed(seeds_G3[k])

    while tf < 3600:

        iter_count += 1
        if iter_w_imp == limits[ttype[G]][ilim]:
            if ilim == 3:
                break
            ilim += 1
            coeff += step
            iter_w_imp = 0
            alpha = [np.random.random() * (round_angle - delta)]  # new random rotation for variant G2
            if variant == 'G3':
                delta = coeff * round_angle / n_slices
            else:
                delta = coeff * round_angle
        print(f'Iteration: {iter_count}, Coeff: {coeff}')

        # Selection of random observations to unfix #
        if variant == 'G1': alpha = [np.random.random() * (round_angle - delta)] # new random rotation for variant G1
        if variant == 'G3': alpha = [np.random.random() * (round_angle - delta) for _ in range(n_slices)]  # new random rotations for variant G3

        for o in range(nO):
            if any(a <= angle[o] < a + delta for a in alpha):  # variables to let free
                xo[o].lb = 0
                xo[o].ub = 1
            else:
                xo[o].lb = selected[o]
                xo[o].ub = selected[o]
            xo[o].setAttr("start", selected[o])

        tm1 = time.time() - t0  # start model
        model.optimize()
        tm2 = time.time() - t0  # end model
        tmodel += tm2 - tm1  # model time

        #print(f'Is the solution optimal? {model.status == gp.GRB.OPTIMAL}')
        #print(f'Solution: {z.getValue()}\n')

        obj_prec = obj
        obj = z.getValue()
        selected = [0] * nO
        for o in range(nO):
            selected[o] = xo[o].x
        if obj_prec == obj:
            iter_w_imp += 1
        else:
            iter_w_imp = 0
            iter_results.append([iter_count, coeff, obj, tm1, tm2, time.time() - t0])

        # rolling slice for variant G3
        if variant == 'G2':
            if alpha[0] + 2 * delta <= round_angle:
                alpha[0] = alpha[0] + delta
            else:
                alpha = [np.random.random() * (round_angle - delta)]  # random rotation
        tf = time.time() - t0

    # end while
    visited = [0] * nPI
    for o in range(nO):
        if selected[o] == 1:
            for p in PIinO[o]:
                visited[p] = 1
    coverture = sum(visited)

    tf = time.time() - t0
    time_values.append(tf)
    iter_values.append(iter_count)
    obj_values.append(obj)
    coverture_values.append(coverture)

    # Save partial solutions
    g = open(f"{ttype[G]}_{coeff_0}_{step}_{variant}_{name}.txt", "a")
    g.write(f'Objective Value = {obj} \n')
    g.write(f'Number of PI covered = {coverture} \n')
    g.write(f'Total execution time = {tf} \n')
    g.write(f'Total model time = {tmodel} \n')
    g.write(f'Number of iterations = {iter_count} \n')
    g.write('Iter Coeff Obj StartM EndM Tend\n')
    for i in iter_results:
        g.write(f'{i[0]} {i[1]} {i[2]} {i[3]} {i[4]} {i[5]} \n')
    g.write('\n')
    g.close()

    # Save final solution
    f = open(f"{ttype[G]}_{coeff_0}_{step}_{variant}(sol)_{name}.txt", "a")
    f.write(f'Objective Value = {obj} \n')
    for i in selected:
        f.write("%d\n" % i)
    f.close()

f = open(f"{ttype[G]}_{coeff_0}_{step}_{variant}_{name}.txt", "a")
f.write('Max objective function = ' + str(max(obj_values)) + '\n')
f.write('Max number of PI covered = ' + str(max(coverture_values)) + '\n')
f.write('Max number of iterations = ' + str(max(iter_values)) + '\n')
f.write('Max total time = ' + str(max(time_values)) + '\n\n')
f.write('Average objective function = ' + str(mean(obj_values)) + '\n')
f.write('Average number of PI covered = ' + str(mean(coverture_values)) + '\n')
f.write('Average number of iterations = ' + str(mean(iter_values)) + '\n')
f.write('Average total time = ' + str(mean(time_values)) + '\n\n')
f.close()

