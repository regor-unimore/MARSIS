"""Command-line interface and default file locations shared by the optimization scripts."""
import argparse
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INSTANCE_DIR = REPO_ROOT / "data" / "Opt_instances"
DEFAULT_GENERAL_DIR = REPO_ROOT / "data" / "Opt_general_data"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "output"

# Class quality weights (set V) used in the paper, depending on the number of quality classes NC
CLASS_WEIGHTS = {'3': [0, 0.5, 1], '4': [0, 0.25, 0.5, 1], '5': [0, 0.25, 0.5, 0.75, 1]}

# Instance names: M_<G>_<L>_<NC>, with granularity G in km and length L in s
INSTANCE_NAME = re.compile(r"M_(01|02|03|05|10)_\d+_[345]")


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show the default values in --help, except for the options whose default is described in their help text."""
    def _get_help_string(self, action):
        if action.default is None:
            return action.help
        return super()._get_help_string(action)


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"must be a positive integer, got {value}")
    return number


def positive_float(value):
    number = float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError(f"must be positive, got {value}")
    return number


def make_parser(description):
    parser = argparse.ArgumentParser(description=description, formatter_class=HelpFormatter)
    parser.add_argument('instance_name', type=str, help="Name of the instance (e.g. M_10_134_4)")
    parser.add_argument('--N', type=positive_int, default=1, help="Number of observations per day")
    parser.add_argument('--W', type=float, nargs=2, default=[1, 1], metavar=('W_OUT', 'W_AOI'),
                        help="Score of the PIs outside and inside an AOI ([1/beta, 1])")
    parser.add_argument('--V', type=float, nargs='+', default=None,
                        help="Class quality weights (default: [0, 0.5, 1] for NC=3, [0, 0.25, 0.5, 1] for NC=4, "
                             "[0, 0.25, 0.5, 0.75, 1] for NC=5)")
    parser.add_argument('--instance-dir', type=Path, default=DEFAULT_INSTANCE_DIR,
                        help="Directory containing the extracted instance folders")
    parser.add_argument('--general-dir', type=Path, default=DEFAULT_GENERAL_DIR,
                        help="Directory containing the general data files (day.txt, orbit.txt, ...)")
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_DIR, help="Directory where results are saved")
    return parser


def add_gurobi_args(parser, time_limit):
    parser.add_argument('--time-limit', type=positive_float, default=time_limit, help="Gurobi time limit (s)")
    parser.add_argument('--threads', type=positive_int, default=1, help="Number of Gurobi threads")


def add_initial_solution_arg(parser):
    parser.add_argument('--initial-solution', type=Path, default=None,
                        help="Solution file of heuristic H2 used as initial solution "
                             "(default: <output-dir>/H2_<instance>.txt)")


def parse_args(parser, check=None):
    """Parse and check the arguments; check(args) runs further script-specific checks."""
    args = parser.parse_args()
    name = args.instance_name
    if not INSTANCE_NAME.fullmatch(name):
        parser.error(f"invalid instance name '{name}', expected M_<G>_<L>_<NC> (e.g. M_10_134_4)")
    NC = name.split('_')[3]
    if args.V is None:
        args.V = list(CLASS_WEIGHTS[NC])
    if len(args.V) != int(NC):
        parser.error(f"--V must contain {NC} weights for instance {name}")
    args.instance_path = args.instance_dir / name
    if not (args.instance_path / f"{name}(1).txt").is_file():
        option = "" if args.instance_dir == DEFAULT_INSTANCE_DIR else f' --instance-dir "{args.instance_dir}"'
        parser.error(f"instance files not found in {args.instance_path}. "
                     f"Extract them first with: python scripts/setup_data.py {name}{option}")
    if not (args.general_dir / "day.txt").is_file():
        parser.error(f"general data files not found in {args.general_dir} (e.g. day.txt)")
    if hasattr(args, 'initial_solution'):
        if args.initial_solution is None:
            args.initial_solution = args.output_dir / f"H2_{name}.txt"
        if not args.initial_solution.is_file():
            parser.error(f"initial solution {args.initial_solution} not found. Run heuristic H2 first with: "
                         f"python scripts/Opt_code/h2.py {name} --N {args.N} --output-dir \"{args.output_dir}\"")
    if check is not None:
        check(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    return args
