import time
from opt_cli import make_parser, parse_args

# Define INSTANCE
parser = make_parser("Code to solve the H2 heuristic")
args = parse_args(parser)
name = args.instance_name

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

# PI in observations
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

# Ranges of observations belonging to the same day
ranges = {}
beg = 0
while beg < nO:
    end = beg
    while end+1 < nO and day[end+1] == day[beg]:
        end = end+1
    ranges[day[beg]] = range(beg, end+1)
    beg = end+1

############################# MAIN #####################################

selected = [0] * nO  # 1 if observation o is selected, 0 otherwise
visited = [0] * nPI  # 1 if point p is covered by at least one observation, 0 otherwise
w_quality = cPIAV[:]  # max quality with which each point is covered, initialized to the past quality
days = {i: 0 for i in set(day)}  # counter for the number of observations selected in each day

t1 = time.time()  # start solution time

z_o = [sum(max(0, (cPIinO[i][j] - w_quality[PIinO[i][j]]) * wPI[PIinO[i][j]]) for j in range(len(PIinO[i])))
       for i in range(nO)]  # compute the expected profit for each observation

while sum(selected) < len(days)*N:

    curS = sorted(range(nO), key=lambda i: z_o[i])
    maxobs = curS[-1]  # observation with the highest profit
    if z_o[maxobs] == -1: break
    selected[maxobs] = 1  # select the observation
    for l, j in enumerate(PIinO[maxobs]):
        w_quality[j] = max(w_quality[j], cPIinO[maxobs][l])  # update the maximum quality for each p in the selected observation
        visited[j] = 1  # points within the selected observation are covered

    days[day[maxobs]] += 1  # update day counter
    if days[day[maxobs]] == N or day[maxobs] == '2023NOV17':
        for i in ranges[day[maxobs]]:
            z_o[i] = -1

    for i in range(nO):  # update the expected profit for each observation
        if z_o[i] != -1:
            z_o[i] = sum(max(0, (cPIinO[i][j] - w_quality[PIinO[i][j]]) * wPI[PIinO[i][j]]) for j in range(len(PIinO[i])))

# compute objective function
ob = 0
for p in range(nPI):
    ob += (w_quality[p] - cPIAV[p]) * wPI[p]

t2 = time.time()  # end solution time

# print output
print("objective function = ", ob)
print("visited PI = ", sum(visited))
print("execution time = ", t2-t1)

# save output
f = open(args.output_dir / f"H2_{name}.txt", "w+")
f.write(f"objective function = {ob}\nvisited PI = {sum(visited)}\nexecution time = {t2 - t1}\n")
for i in selected:
    f.write("%d\n" % i)
f.close()