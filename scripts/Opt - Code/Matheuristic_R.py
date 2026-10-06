import random
from statistics import mean
import numpy as np
import time
import gurobipy as gp
from gurobipy import GRB
from opt_cli import make_parser, parse_args, add_gurobi_args, add_initial_solution_arg

# Define INSTANCE
parser = make_parser("Code to solve the R matheuristics")
parser.add_argument('variant_name', type=str, choices=['R1', 'R2', 'R3'], help="R1, R2 or R3")
add_gurobi_args(parser, time_limit=360)
add_initial_solution_arg(parser)
parser.add_argument('--n-iter', type=int, default=20, help="Number of repetitions (one per random seed)")
parser.add_argument('--first-iter', type=int, default=1, help="Index of the first repetition, i.e. of the first seed used")
parser.add_argument('--max-time', type=float, default=3600, help="Time limit of each repetition (s)")
parser.add_argument('--coeff0', type=float, default=0.85, help="Initial fixing coefficient")
parser.add_argument('--step', type=float, default=0.10, help="Update step of the fixing coefficient")
parser.add_argument('--policy', type=int, choices=[1, 2, 3, 4, 5], default=None,
                    help="Updating policy (default: best policy for the granularity G, "
                         "i.e. 1 for G=10, 2 for G=05, 3 for G=03, 4 for G=02 and 5 for G=01)")
args = parse_args(parser)
name = args.instance_name
variant = args.variant_name

name_char = name.split('_')

# Parameters
G = name_char[1]    # granularity
L = name_char[2]    # length
NC = name_char[3]   # number of quality classes
gra = {'01': 1000, '02': 2000, '03': 3000, '05': 5000, '10': 10000}
nO = 113651         # number of observations

# Default settings
N = args.N              # number of observation per day
V = args.V              # class quality weights
W = args.W              # [1/beta, 1], score for AOI

#test parameters
seeds_R1 = [3, 91, 572, 103, 4246, 9999, 19873, 14593, 75380, 91704, 2, 78, 412, 36, 97654, 2378, 456, 5, 10980, 3780]
seeds_R2 = [9, 71, 672, 193, 6846, 9099, 19073, 14093, 70380, 91004, 14, 34, 150, 249, 533, 1033, 2345, 3572, 8473, 10002]
seeds_R3 = [8, 51, 739, 291, 846, 5099, 12073, 1093, 30380, 81004, 88, 531, 7739, 8291, 1846, 4099, 22073, 13093, 70380, 51004]
policies = {1: [160, 80, 40, 20], 2: [80, 40, 20, 10], 3: [40, 20, 10, 5], 4: [20, 10, 5, 2], 5: [10, 5, 2, 1]}  # max n. of iterations without improvement for each coefficient value
best_policy = {'10': 1, '05': 2, '03': 3, '02': 4, '01': 5}  # default updating policy for each granularity
policy = args.policy if args.policy is not None else best_policy[G]
limits = policies[policy]
if not (args.first_iter >= 1 and args.n_iter >= 1 and args.first_iter + args.n_iter - 1 <= len(seeds_R1)):
    parser.error(f"--first-iter and --n-iter must select repetitions between 1 and {len(seeds_R1)} (one per available seed)")

#  Information on PIs
cPIAV = []
wPI = []
with open(args.instance_path / f"{name}(1).txt", 'r') as file:
    lines = file.readlines()
    nPI = len(lines)        # number of PI
    for line in lines:
        line = list(line.strip().split(' '))
        cPIAV.append(V[int(line[5])])
        wPI.append(W[int(line[6])])

#  PI in observations
with open(args.instance_path / f"{name}(2).txt", "r") as file:
    PIinO = [[int(i) for i in line.strip().split("\t")] if line != "" and line != "\n" else []
             for line in file.readlines()]

# Class of the quality of PIs in observations, weighted through set V
with open(args.instance_path / f"{name}(4).txt", "r") as file:
    cPIinO = [[V[int(i)] for i in line.strip().split("\t")] if line != "" and line != "\n" else []
              for line in file.readlines()]

# Day of the observations
with open(args.general_dir / "day.txt", "r") as file:
    day = [line for line in file.readlines()]

# Orbit of the observations
with open(args.general_dir / "orbit.txt", "r") as file:
    orbit = [line for line in file.readlines()]

# Heuristic solution used to initialize the x variables
with open(args.initial_solution, "r") as file:
    obj = float(file.readline().split("=")[1])
    file.readline()
    file.readline()
    selected_h2 = [int(file.readline()) for i in range(nO)]

def variable_fixing_R1():
    var_fix = np.random.random(nO)  # Selection of random observations to fix

    for o in range(nO):
        if var_fix[o] < coeff:  # fix to the values of the previous solution
            xo[o].lb = selected[o]
            xo[o].ub = selected[o]
        else:
            xo[o].lb = 0
            xo[o].ub = 1
        xo[o].setAttr("start", selected[o])

def variable_fixing_R2():
    day_fix = {day[i]: 0 for i in range(nO)}  # 1 if day d is fixed, 0 otherwise
    day_fix_r = np.random.random(len(day_fix))
    idx = 0
    for d in day_fix.keys():
        if day_fix_r[idx] < coeff:
            day_fix[d] = 1
        idx += 1

    # Fix the selected variables
    for o in range(nO):
        if day_fix[day[o]] == 1:  # fix to the values of the previous solution
            xo[o].lb = selected[o]
            xo[o].ub = selected[o]
        else:
            xo[o].lb = 0
            xo[o].ub = 1
        xo[o].setAttr("start", selected[o])

def variable_fixing_R3():
    orbit_fix = {orbit[i]: 0 for i in range(nO)}
    orbit_fix_r = np.random.random(len(orbit_fix))
    idx = 0
    for o in orbit_fix.keys():
        if orbit_fix_r[idx] < coeff:
            orbit_fix[o] = 1
        idx += 1

    for o in range(nO):
        if orbit_fix[orbit[o]] == 1:
            xo[o].lb = selected[o]
            xo[o].ub = selected[o]
        else:
            xo[o].lb = 0
            xo[o].ub = 1
        xo[o].setAttr("start", selected[o])

######################### MODEL ####################################################
model = gp.Model()
model.setParam("TimeLimit", args.time_limit)
model.setParam(GRB.Param.Threads, args.threads)
model.setParam(GRB.Param.Method, 2)

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

######################### MODEL EXECUTION ###############################

# parameters
n_iter = args.n_iter
reps = "" if (args.first_iter, n_iter) == (1, len(seeds_R1)) else f"_rep{args.first_iter}-{args.first_iter + n_iter - 1}"  # repetitions in the output file name
step = args.step
coeff_0 = args.coeff0

obj_values = []
coverture_values = []
time_values = []
iter_values = []

for k in range(args.first_iter - 1, args.first_iter - 1 + n_iter):

    t0 = time.time()  # starting time
    tf = 0  # total execution time
    tmodel = 0  # total model time

    iter_w_imp = 0
    iter_count = 0
    ilim = 0
    coeff = coeff_0

    iter_results = []
    selected = selected_h2  # initialize the solution
    if variant == 'R1': np.random.seed(seeds_R1[k])
    if variant == 'R2': np.random.seed(seeds_R2[k])
    if variant == 'R3': np.random.seed(seeds_R3[k])

    while tf < args.max_time:

        # Update fixing coefficient
        iter_count += 1
        if iter_w_imp == limits[ilim]:
            if ilim == 3:
                break
            ilim += 1
            coeff -= step
            iter_w_imp = 0
        print(f'Iteration: {iter_count}, Coeff: {coeff}')

        if variant == 'R1': variable_fixing_R1()
        if variant == 'R2': variable_fixing_R2()
        if variant == 'R3': variable_fixing_R3()

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
            iter_results.append([iter_count, coeff, obj, tm1, tm2, time.time()-t0])
        tf = time.time() - t0  # T END

    # end while
    visited = [0] * nPI
    for o in range(nO):
        if selected[o] == 1:
            for p in PIinO[o]:
                visited[p] = 1
    coverture = sum(visited)
    tf = time.time()-t0
    time_values.append(tf)
    iter_values.append(iter_count)
    obj_values.append(obj)
    coverture_values.append(coverture)

    # Save partial solutions
    g = open(args.output_dir / f"P{policy}_{coeff_0}_{step}_{variant}_{name}{reps}.txt", "a")
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
    f = open(args.output_dir / f"P{policy}_{coeff_0}_{step}_{variant}(sol)_{name}{reps}.txt", "a")
    f.write(f'Objective Value = {obj} \n')
    for i in selected:
        f.write("%d\n" % i)
    f.close()

f = open(args.output_dir / f"P{policy}_{coeff_0}_{step}_{variant}_{name}{reps}.txt", "a")
f.write('Max objective function = ' + str(max(obj_values)) + '\n')
f.write('Max number of PI covered = ' + str(max(coverture_values)) + '\n')
f.write('Max number of iterations = ' + str(max(iter_values)) + '\n')
f.write('Max total time = ' + str(max(time_values)) + '\n\n')
f.write('Average objective function = ' + str(mean(obj_values)) + '\n')
f.write('Average number of PI covered = ' + str(mean(coverture_values)) + '\n')
f.write('Average number of iterations = ' + str(mean(iter_values)) + '\n')
f.write('Average total time = ' + str(mean(time_values)) + '\n\n')
f.close()
