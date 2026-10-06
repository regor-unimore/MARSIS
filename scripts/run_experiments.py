"""Run the experiments listed in results/experiments.csv.

Each row of the manifest corresponds to one command and to the rows of results/Computational Results.xlsx
that it reproduces. The missing instances are extracted automatically, and heuristic H2 is run before the
Lagrangian relaxation and the matheuristics when its solution (used as initial solution) is not available.

Usage examples:
    python scripts/run_experiments.py --dry-run                              # print all the commands
    python scripts/run_experiments.py --sheet Heuristics --instance M_10_134_4
    python scripts/run_experiments.py --algorithm ILP Lagrangian --skip-existing
    python scripts/run_experiments.py --id 001 002
"""
import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "results" / "experiments.csv"
INSTANCE_DIR = REPO_ROOT / "data" / "Opt - Instances"

parser = argparse.ArgumentParser(description="Run the experiments of the results workbook",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--manifest', type=Path, default=DEFAULT_MANIFEST, help="Experiment manifest")
parser.add_argument('--id', nargs='+', help="Run only the experiments with these ids")
parser.add_argument('--sheet', nargs='+', help="Run only the experiments of these workbook sheets "
                                               "(ILP, Heuristics, Lagrangian, 'Matheuristic all')")
parser.add_argument('--algorithm', nargs='+', help="Run only these algorithms (ILP, H1, H2, Lagrangian, Matheuristic)")
parser.add_argument('--instance', nargs='+', help="Run only the experiments on these instances")
parser.add_argument('--variant', nargs='+', help="Run only these matheuristic variants (R1, ..., G3)")
parser.add_argument('--N', nargs='+', help="Run only the experiments with these numbers of observations per day")
parser.add_argument('--skip-existing', action='store_true', help="Skip the experiments whose output files already exist")
parser.add_argument('--dry-run', action='store_true', help="Print the commands without running them")
args = parser.parse_args()

with open(args.manifest, newline='', encoding='utf-8') as file:
    manifest = list(csv.DictReader(file))

filters = {'id': args.id, 'sheet': args.sheet, 'algorithm': args.algorithm, 'instance': args.instance,
           'variant': args.variant, 'N': args.N}
selected = [e for e in manifest if all(v is None or e[k] in v for k, v in filters.items())]
if not selected:
    sys.exit("No experiment matches the selection")


def outputs_exist(experiment):
    return all((REPO_ROOT / f).is_file() for f in experiment['output_files'].split(';'))


def run(command):
    print(f"\n>>> {command}", flush=True)
    if args.dry_run:
        return
    argv = shlex.split(command)
    if argv[0] == 'python':
        argv[0] = sys.executable
    result = subprocess.run(argv, cwd=REPO_ROOT)
    if result.returncode != 0:
        sys.exit(f"Command failed with exit code {result.returncode}: {command}")


# H2 runs needed as initial solution, taken from the manifest
h2 = {(e['instance'], e['N']): e for e in manifest if e['algorithm'] == 'H2'}
scheduled = set()  # commands already run (or printed, with --dry-run)

print(f"{len(selected)} experiment(s) selected")
for experiment in selected:
    if args.skip_existing and outputs_exist(experiment):
        print(f"\n--- {experiment['id']}: output already available, skipping")
        continue
    instance = experiment['instance']
    setup = f'python scripts/setup_data.py {instance}'
    if not (INSTANCE_DIR / instance / f"{instance}(1).txt").is_file() and setup not in scheduled:
        run(setup)
        scheduled.add(setup)
    if experiment['depends_on']:
        dependency = h2[(instance, experiment['N'])]
        if not outputs_exist(dependency) and dependency['command'] not in scheduled:
            print(f"\n--- {dependency['id']}: required as initial solution")
            run(dependency['command'])
            scheduled.add(dependency['command'])
    scheduled.add(experiment['command'])
    print(f"\n--- {experiment['id']}: {experiment['sheet']} (workbook rows {experiment['workbook_rows']})")
    run(experiment['command'])
