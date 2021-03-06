import numpy as np
import pandas as pd
import xgboost as xgb
import scipy.stats as ss
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split


def to_categorical(df, cols):
    for col in cols:
        df[col] = pd.Categorical(df[col])
        df[col] = df[col].cat.codes
    return df


def to_numeric(df, cols):
    for col in cols:
        df[col] = pd.to_numeric(df[col])
    return df


def split_data(df):
    # x = df.drop(columns=['ma_binding_score'])
    # y = df['ma_binding_score']
    train = df['days_since_g0'] < 20200401
    train_data = df[train]
    test_data = df[~train]
    x_train = train_data.drop(columns=['ma_binding_score'])
    y_train = train_data['ma_binding_score']
    x_test = test_data.drop(columns=['ma_binding_score'])
    y_test = test_data['ma_binding_score']
    # x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.1, random_state=101)
    return x_train, x_test, y_train, y_test


def load_data(filename='dataframe.pickle'):
    df = pd.read_pickle(filename)
    df = to_categorical(df, ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'hla', 'location'])
    df = to_numeric(df, ['pos', 'edit_dist', 'days_since_g0'])
    x_train, x_test, y_train, y_test = split_data(df)
    return x_train, x_test, y_train, y_test


def train(x_train, y_train):
    reg = xgb.XGBRegressor()
    reg.fit(x_train, y_train)
    return reg


def evaluate(reg, x_test, y_test):
    preds = pd.Series(reg.predict(x_test))
    preds = np.nan_to_num(preds.fillna(0))
    y_test = np.nan_to_num(y_test.fillna(0))
    rmse = np.sqrt(mean_squared_error(y_test, preds))
    range_val = np.ptp(y_test)
    nrmse = rmse/range_val
    print('NRMSE:', nrmse)


def distribution_test(dist):
    chisquare_val, p_val = ss.chisquare(dist)
    print('Chisquare:', chisquare_val)
    print('P:', p_val)


if __name__ == '__main__':
    x_train, x_test, y_train, y_test = load_data('dataframe_347k.pickle')
    reg = train(x_train, y_train)
    evaluate(reg, x_test, y_test)
    distribution_test(y_test)
