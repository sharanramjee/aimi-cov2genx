import numpy as np
import pandas as pd
import xgboost as xgb
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
    x = df.drop(columns=['ma_binding_score'])
    y = df['ma_binding_score']
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.1, random_state=101)
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
    preds = reg.predict(x_test)
    rmse = np.sqrt(mean_squared_error(y_test, preds))
    print('RMSE:', rmse)


if __name__ == '__main__':
    x_train, x_test, y_train, y_test = load_data('dataframe.pickle')
    reg = train(x_train, y_train)
    evaluate(reg, x_test, y_test)
