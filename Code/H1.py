import time
import argparse

# Define INSTANCE
#parser = argparse.ArgumentParser(description="Code to solve the H1 heuristic")
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

# Class of the quality of PIs in observations, weighted through set V
with open(f"{name}/{name}(4).txt", "r") as file:
    cPIinO = [[V[int(i)] for i in line.strip().split("\t")]
                    if line != "" and line != "\n" else [] for line in file.readlines()]

# Day of the observations
with open(f"day.txt", "r") as file:
    day = [line for line in file.readlines()]

############################# MAIN #####################################

selected = [0]*nO  # 1 if observation o is selected, 0 otherwise
visited = [0]*nPI  # 1 if point p is covered by at least one observation, 0 otherwise
w_quality = cPIAV[:]  # max quality with which each point is covered, initialized to the past quality

t1 = time.time()  # start solution time

beg = 0
while beg < nO:
    end = beg
    while end+1 < nO and day[end+1] == day[beg]:
        end = end+1
    cur = PIinO[beg:(end+1)]  # set of observations belonging to the same day
    curI = range(beg, end+1)  # indices of observations belonging to the same day

    for k in range(N):  # N = number of observation to select per day
        curS = sorted(zip(cur, curI), key=lambda i: sum(max(0, (cPIinO[i[1]][j] - w_quality[i[0][j]]) * wPI[i[0][j]])
                                        for j in range(len(i[0]))))  # compute the expected profit for each observation
        if len(curS) >= 1:
            maxobs = curS[-1]  # observation with the highest profit
            selected[maxobs[1]] = 1  # select the observation
            for l, j in enumerate(maxobs[0]):
                w_quality[j] = max(w_quality[j], cPIinO[maxobs[1]][l])  # update the maximum quality for each p in the selected observation
                visited[j] = 1  # points within the selected observation are covered
    beg = end+1  # go to the next day

# compute objective function
ob = 0
for p in range(nPI):
    ob += (w_quality[p] - cPIAV[p]) * wPI[p]

t2 = time.time()  # end solution time

# print output
print("objective function = ", ob)
print("visited PI = ", sum(visited))
print("execution time = ", t2 - t1)

# save output
f = open(f"H1_{name}.txt", "w+")
f.write(f"objective function = {ob}\nvisited PI = {sum(visited)}\nexecution time = {t2 - t1}\n")
for i in selected:
    f.write("%d\n" % i)
f.close()

