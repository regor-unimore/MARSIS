"""Extract the compressed optimization instances into the layout expected by the scripts.

Each archive data/Opt - Instances/<instance>.zip is extracted into
data/Opt - Instances/<instance>/, which then contains <instance>(1).txt ... <instance>(4).txt.

Usage:
    python scripts/setup_data.py M_10_134_4 M_10_268_4   # extract selected instances
    python scripts/setup_data.py --all                   # extract all 25 instances
"""
import argparse
import sys
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INSTANCE_DIR = REPO_ROOT / "data" / "Opt - Instances"

parser = argparse.ArgumentParser(description="Extract the optimization instances")
parser.add_argument('instances', nargs='*', help="Names of the instances to extract (e.g. M_10_134_4)")
parser.add_argument('--all', action='store_true', help="Extract all the available instances")
parser.add_argument('--instance-dir', type=Path, default=DEFAULT_INSTANCE_DIR,
                    help="Directory containing the instance archives (default: %(default)s)")
parser.add_argument('--force', action='store_true', help="Re-extract instances that are already extracted")
args = parser.parse_args()
if args.all and args.instances:
    parser.error("specify either instance names or --all, not both")

available = sorted(p.stem for p in args.instance_dir.glob("M_*.zip"))
if not available:
    parser.error(f"no instance archives (M_*.zip) found in {args.instance_dir}")
names = available if args.all else args.instances
if not names:
    parser.error(f"specify at least one instance or --all. Available instances: {', '.join(available)}")

missing = [n for n in names if n not in available]
if missing:
    sys.exit(f"Unknown instance(s): {', '.join(missing)}. Available instances: {', '.join(available)}")

for name in names:
    target = args.instance_dir / name
    expected = [target / f"{name}({i}).txt" for i in range(1, 5)]
    if not args.force and all(p.is_file() for p in expected):
        print(f"{name}: already extracted, skipping")
        continue
    with zipfile.ZipFile(args.instance_dir / f"{name}.zip") as archive:
        archive.extractall(target)
    absent = [p.name for p in expected if not p.is_file()]
    if absent:
        sys.exit(f"{name}: archive does not contain {', '.join(absent)}")
    print(f"{name}: extracted to {target}")
