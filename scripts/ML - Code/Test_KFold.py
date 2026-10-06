import time
from statistics import mean
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.base import clone
from joblib import dump

REPO_ROOT = Path(__file__).resolve().parents[2]
ML_DATA_DIR = REPO_ROOT / "data" / "ML - Data"

# Define MODEL
parser = argparse.ArgumentParser(description="Code to train the regression models",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('model_name', type=str, choices=['linear_regression', 'knn', 'decision_tree', 'gr_boosting', 'random_forest'],
                    help="Regression model to train")
parser.add_argument('--dataset', type=Path, default=ML_DATA_DIR / "MARSIS_historical_dataset.csv",
                    help="Historical dataset (to be downloaded from the link in data/ML - Data/MARSIS_historical_dataset)")
parser.add_argument('--orbits-to-remove', type=Path, default=ML_DATA_DIR / "orbit_to_remove",
                    help="List of the orbits to exclude from the training")
parser.add_argument('--output-dir', type=Path, default=REPO_ROOT / "output" / "ML", help="Directory where results are saved")
args = parser.parse_args()
if not args.dataset.is_file():
    parser.error(f"dataset {args.dataset} not found. Download it from the link in data/ML - Data/MARSIS_historical_dataset")
args.output_dir.mkdir(parents=True, exist_ok=True)
model_name = args.model_name

orbit_to_remove = []
with open(args.orbits_to_remove) as file:
    for line in file:
        orbit_to_remove.append(float(line))
df = pd.read_csv(args.dataset, sep=";")

frequency_to_keep = 4000000.0
df = df[df['FM_data_frequency'] == frequency_to_keep]
df = df[~df.FM_data_orbit_number.isin(orbit_to_remove)]
df['FM_data_solar_longitude_cos'] = np.cos(df['FM_data_solar_longitude'])
df['FM_data_solar_longitude_sin'] = np.sin(df['FM_data_solar_longitude'])
X = df.drop(columns=['FM_data_ephemeris_time', 'FM_data_F10_7_index', 'FM_data_frequency',
                     'FM_data_median_corrected_echo_power', 'FM_data_orbit_number', 'FM_data_peak_corrected_echo_power',
                     'FM_data_peak_distorted_echo_power', 'FM_data_peak_simulated_echo_power', 'FM_data_solar_longitude'])
X = X.to_numpy()
y = df['FM_data_peak_distorted_echo_power'].to_numpy()

#--- Regression model parameters ---#
algorithm = {
    'knn': KNeighborsRegressor(),
    'linear_regression': LinearRegression(),
    'decision_tree': DecisionTreeRegressor(random_state=0, criterion='squared_error', max_depth=20),
    'gr_boosting': GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=25, random_state=0, criterion='squared_error'),
    'random_forest': RandomForestRegressor(n_estimators=100, criterion='squared_error', random_state=0, max_depth=25, n_jobs=16)
}

kf = KFold(n_splits=10, shuffle=False)
kf.get_n_splits(X)

f = open(args.output_dir / f"MAE_{model_name}.txt", "w")
fold = -1
y_all_pred = np.zeros(y.shape)

t0 = time.time()  # start of the cross-validation
for train_index, test_index in kf.split(X, y):
    fold = fold + 1
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    f.write('TRAIN: ' + str(train_index)+'\n')
    f.write('TEST: ' + str(test_index)+'\n')

    t1 = time.time()
    model = clone(algorithm[model_name])
    if model_name in ["knn", "linear_regression"]:
        model = make_pipeline(StandardScaler(), model)  # the fitted scaler is saved together with the model
    model.fit(X_train, y_train)
    dump(model, args.output_dir / f'{model_name}_fold{fold}.joblib')
    y_pred = model.predict(X_test)
    t2 = time.time()

    f.write(f'  MAE = {mean_absolute_error(y_test, y_pred)}\n')
    f.write(f'  MAPE = {mean_absolute_percentage_error(y_test, y_pred)}\n')
    f.write(f'  MSE = {mean_squared_error(y_test, y_pred)}\n')
    f.write(f'  Execution time = {t2 - t1}\n')
    y_all_pred[test_index] = y_pred

t3 = time.time()
f.write('\nGlobal MAE = ' + str(mean_absolute_error(y, y_all_pred)))
f.write('\nGlobal MAPE = ' + str(mean_absolute_percentage_error(y, y_all_pred)))
f.write('\nGlobal MSE = ' + str(mean_squared_error(y, y_all_pred)))
f.write(f'\nTime = {t3-t0}')
f.close()
np.savetxt(args.output_dir / f"y_all_pred_noantennas_{model_name}.txt", y_all_pred)

