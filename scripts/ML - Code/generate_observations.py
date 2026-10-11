from math import floor, ceil, sqrt
import argparse
import re
from pathlib import Path
import pandas as pd

"""
    Output qualità dei punti coperti per ogni osservazione programmata:
    file (3) con la qualità continua e file (4) con la classe di qualità di ogni PI coperto
"""

REPO_ROOT = Path(__file__).resolve().parents[2]

# Define INSTANCE
parser = argparse.ArgumentParser(description="Code to compute the quality of the PIs covered by each observation "
                                             "(files (3) and (4) of an instance) from the ML predictions",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('instance_name', type=str, help="Name of the instance (e.g. M_10_134_4)")
parser.add_argument('--dataset', type=Path, default=REPO_ROOT / "data" / "ML - Data" / "MARSIS_dataset_2023_2025.zip",
                    help="Dataset for the future planning")
parser.add_argument('--predictions', type=Path, default=REPO_ROOT / "data" / "Opt - General Data" / "predictions.txt",
                    help="ML predictions for each sampled time")
parser.add_argument('--instance-dir', type=Path, default=REPO_ROOT / "data" / "Opt - Instances",
                    help="Directory containing the extracted instance folders (files (1) and (2) are read from here)")
parser.add_argument('--output-dir', type=Path, default=REPO_ROOT / "output" / "Instances",
                    help="Directory where files (3) and (4) are saved")
args = parser.parse_args()
name = args.instance_name
if not re.fullmatch(r"M_(01|02|03|05|10)_\d+_[345]", name):
    parser.error(f"invalid instance name '{name}', expected M_<G>_<L>_<NC> (e.g. M_10_134_4)")
for path in (args.dataset, args.predictions):
    if not path.is_file():
        parser.error(f"file {path} not found")

name_char = name.split('_')

# Parameters
G = name_char[1]    # granularity
L = int(name_char[2])    # length
NC = name_char[3]   # number of quality classes
gra = {'01': 1000, '02': 2000, '03': 3000, '05': 5000, '10': 10000}
nO = 113651         # number of observations
thresholds = {'3': [9, 15], '4': [9, 12, 15], '5': [9, 12, 13.5, 15]}  # quality thresholds (dB) between classes

instance_path = args.instance_dir / name
if not (instance_path / f"{name}(1).txt").is_file():
    parser.error(f"instance files not found in {instance_path}. Extract them first with: python scripts/setup_data.py {name}")

df = pd.read_csv(args.dataset, sep=";")
time = [float(t) for t in list(df.ephemeris_time_2023_2025)]
altitude = [float(a) for a in list(df.altitude_2023_2025)]
xo = [float(x) for x in list(df.x_2023_2025)]
yo = [float(y) for y in list(df.y_2023_2025)]

minAlt = min(altitude)
maxAlt = max(altitude)

with open(args.predictions, 'r') as file:
    quality = [float(i) for i in file.readlines()]

if len(time) != nO or len(quality) != nO:
    parser.error(f"the dataset and the predictions must contain {nO} sampled times "
                 f"(found {len(time)} and {len(quality)})")

#  PI in observations
with open(instance_path / f"{name}(2).txt", "r") as file:
    PIinO = [[int(i) for i in line.strip().split("\t")]
            if line.strip() != "" else [] for line in file.readlines()]

#  Information on PIs
PI = {}  # id of the PI in each grid point
with open(instance_path / f"{name}(1).txt", 'r') as file:
    for line in file.readlines():
        line = list(line.strip().split(' '))
        PI[(int(line[1]), int(line[2]))] = int(line[0])

args.output_dir.mkdir(parents=True, exist_ok=True)
g = open(args.output_dir / f"{name}(3).txt", "w")
h = open(args.output_dir / f"{name}(4).txt", "w")

step = gra[G]  # distance between two consecutive grid points
for i in range(nO):

    obs = [[xo[i], yo[i], altitude[i], quality[i]]]  # observation starting in i

    end = i
    while end + 1 < nO and time[end + 1] - time[i] <= L:
        end = end + 1  # alla fine end sarà l'ultimo istante incluso

    for j in range(i, end):  # tutti questi indici di istanti racchiudono punti a una distanza di max L secondi da i
        for k in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            if (time[j] + k * (time[j + 1] - time[j])) - time[i] <= L:  # prendo degli istanti di tempo intermedi
                obs.append([xo[j] + k * (xo[j + 1] - xo[j]), yo[j] + k * (yo[j + 1] - yo[j]),
                            altitude[j] + k * (altitude[j + 1] - altitude[j]),
                            quality[j] + k * (quality[j + 1] - quality[j])])
            else:
                break

    q_cov = {p: 0 for p in PIinO[i]}

    for j in obs:  # obs[i] = x,y,a,q
        radius = (1500 + 2000 * (j[2] - minAlt) / (maxAlt - minAlt))
        for x in range(floor((j[0] - radius) / step) * step, ceil((j[0] + radius) / step) * step + 1, step):
            for y in range(floor((j[1] - radius) / step) * step, ceil((j[1] + radius) / step) * step + 1, step):
                if (x, y) in PI and PI[(x, y)] in q_cov:
                    if sqrt((x - j[0]) * (x - j[0]) + (y - j[1]) * (y - j[1])) <= radius:
                        q_cov[PI[(x, y)]] = max(q_cov[PI[(x, y)]], j[3])

    # output
    for k in q_cov.keys():
        g.write("%f\t" % q_cov[k])
        h.write("%d\t" % sum(q_cov[k] >= t for t in thresholds[NC]))  # quality class
    g.write("\n")
    h.write("\n")

g.close()
h.close()
