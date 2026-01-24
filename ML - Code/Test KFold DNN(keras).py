import time
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error, mean_squared_error
from keras.models import Sequential, save_model
from keras.layers import Dense, Dropout
from keras.callbacks import EarlyStopping
from sklearn.preprocessing import StandardScaler
from joblib import dump

orbit_to_remove = []
with open('orbit_to_remove') as file:
    for line in file:
        orbit_to_remove.append(float(line))
df = pd.read_csv("all_past.csv", sep=";")

frequency_to_keep = 4000000.0
df = df[df['FM_data_frequency'] == frequency_to_keep]
df = df[~df.FM_data_orbit_number.isin(orbit_to_remove)]
df['FM_data_solar_longitude_cos'] = np.cos(df['FM_data_solar_longitude'])
df['FM_data_solar_longitude_sin'] = np.sin(df['FM_data_solar_longitude'])

X = df.drop(columns=['FM_data_ephemeris_time', 'FM_data_F10_7_index', 'FM_data_frequency',
                     'FM_data_median_corrected_echo_power', 'FM_data_orbit_number', 'FM_data_peak_corrected_echo_power',
                     'FM_data_peak_distorted_echo_power', 'FM_data_peak_simulated_echo_power', 'FM_data_solar_longitude'])
col_names = X.columns.tolist()
X = X.to_numpy()

y = df['FM_data_peak_distorted_echo_power'].to_numpy()

kf = KFold(n_splits=10, shuffle=False)
kf.get_n_splits(X)

fold = -1
y_all_pred = np.zeros(y.shape)
f = open("MAE_nn_chrono.txt", "w")
t1 = time.time()
for train_index, test_index in kf.split(X):
    fold = fold + 1
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    f.write('TRAIN: ' + str(train_index)+'\n')
    f.write('TEST: ' + str(test_index)+'\n')
    dump(scaler, f'scaler_chrono_{fold}.save')

    # Neural Network
    es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=10)
    model = Sequential()
    model.add(Dense(800, input_dim=X.shape[1], activation='relu'))  # 300
    model.add(Dropout(0.5))
    model.add(Dense(400, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(200, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(100, activation='relu'))
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='relu'))
    model.compile(loss='mse', optimizer='adam', metrics=['mean_absolute_error'])
    model.fit(X_train, y_train, validation_split=0.1, epochs=1000, callbacks=[es])
    save_model(model, f'nn_chrono_model_{fold}.keras')
    # evaluate the model
    _, train_mse = model.evaluate(X_train, y_train, verbose=0)
    _, test_mse = model.evaluate(X_test, y_test, verbose=0)
    print('Train: %.3f, Test: %.3f' % (train_mse, test_mse))
    t2 = time.time()

    y_pred = model.predict(X_test)
    f.write(f'MAE = {mean_absolute_error(y_test,y_pred)}\n')
    f.write(f'MAPE = {mean_absolute_percentage_error(y_test, y_pred)}\n')
    f.write(f'MSE = {mean_squared_error(y_test, y_pred)}\n')
    f.write(f'Execution time = {t2 - t1}\n')
    for i in range(len(test_index)):
        y_all_pred[test_index[i]] = y_pred[i]

t3 = time.time()
f.write(f'\nGlobal MAE = {mean_absolute_error(y, y_all_pred)}')
f.write(f'\nGlobal MAPE = {mean_absolute_percentage_error(y, y_all_pred)}')
f.write(f'\nGlobal MSE = {mean_squared_error(y, y_all_pred)}')
f.write(f'\nTime = {t3-t1}')
f.close()
np.savetxt("nn_y_all_pred_chrono.txt", y_all_pred)