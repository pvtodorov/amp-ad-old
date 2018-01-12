import pandas as pd
import random
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import uuid
from tqdm import tqdm
import csv
import os
from collections import defaultdict
import numpy as np
from itertools import combinations


def recursivedict():
    return defaultdict(recursivedict)


def check_or_create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)


def load_msbb_data(path):
    df = pd.read_table(path)
    df = df.dropna(axis=0, how='any')
    return df


def get_train_test_df(df, column, subset):
    df_train = df[~df[column].isin([subset])]
    df_test = df[df[column].isin([subset])]
    return df_train, df_test


def get_data_cols(df, meta_cols):
    """ Given a dataframe and its meta columns, get back a list of the data
    columns from the dataframe.
    """
    cols = df.columns.tolist()
    data_cols = [x for x in cols if x not in meta_cols]
    return data_cols


def sample_data_cols(data_cols, k):
    """ Select a random sample of k-items from a list of columns.
    """
    sampled_cols = random.sample(data_cols, k)
    return [x for x in data_cols if x in sampled_cols]


def get_X_y(df, target, data_cols):
    """ Given a dataframe, a target col, and a sample of data column names,
    generate X, an numpy array of the data, and y, a list of target values
    corresponding to each row.
    """
    X = df.as_matrix(columns=data_cols)
    y = df[target].tolist()
    return X, y


def binarize_Braak_scores(y):
    """ Take Braak scores and binarize them such that:
    [1, 2, 3] -> 0
    [4, 5] -> 1
    """
    y_bin = [1 if x >= 4 else 0 for x in y]
    return y_bin


def fit_RF_regressor(X, y, rf_params):
    """ Fit and return a RandomForestRegressor object to the provided data """
    rf = RandomForestRegressor(**rf_params)
    rf.fit(X, y)
    return rf


def fit_RF_classifier(X, y, rf_params):
    """ Fit and return a RandomForestClassifier object to the provided data """
    rf = RandomForestClassifier(**rf_params)
    rf.fit(X, y)
    return rf


def fit_RF(X, y, task_type, rf_params=None):
    if rf_params is None:
        rf_params = {"n_estimators": 100,
                     "max_features": 'auto',
                     "n_jobs": 1,
                     "oob_score": True}
    if task_type == 'classification':
        rf = fit_RF_classifier(X, y, rf_params)
    elif task_type == 'regression':
        rf = fit_RF_regressor(X, y, rf_params)
    return rf


def predict_RF(rf, X):
    predictions = rf.predict_proba(X)
    predictions = [x[1] for x in predictions]
    return predictions


def gen_background_predictions(df, target, data_cols,
                               subset_col, subsets,
                               interval=10, max_cols=1000):
    bcg_predictions = recursivedict()
    for subset in subsets:
        df_train, df_test = get_train_test_df(df, subset_col, subset)
        for k in tqdm(range(10, max_cols + interval, interval)):
            selected_cols = sample_data_cols(data_cols, k)
            X_train, y_train = get_X_y(df_train, target, selected_cols)
            y_train = binarize_Braak_scores(y_train)
            X_test, y_test = get_X_y(df_test, target, selected_cols)
            rf = fit_RF(X_train, y_train, 'classification')
            predictions = predict_RF(rf, X_test)
            sub_ids = df_test['ID'].tolist()
            for id, p in zip(sub_ids, predictions):
                bcg_predictions[subset][id][k] = p
    return bcg_predictions


def gen_bcg_predictions(estimator, df, target, data_cols,
                        subset_col, subsets,
                        interval=10, max_cols=1000):
    bcg_predictions = recursivedict()
    for subset in subsets:
        df_train, df_test = get_train_test_df(df, subset_col, subset)
        for k in tqdm(range(10, max_cols + interval, interval)):
            selected_cols = sample_data_cols(data_cols, k)
            X_train, y_train = get_X_y(df_train, target, selected_cols)
            y_train = [int(x) for x in y_train]
            y_train = np.array(y_train)
            X_test, y_test = get_X_y(df_test, target, selected_cols)
            e = estimator()
            e = e.fit(X_train, y_train)
            predictions = e.predict(X_test)
            sub_ids = df_test['ID'].tolist()
            for id, p in zip(sub_ids, predictions):
                bcg_predictions[subset][id][k] = p
    return bcg_predictions


def save_predictions(predictions, folder):
    for subset in predictions:
        df = pd.DataFrame(predictions[subset])
        df = df.transpose().reset_index().rename(columns={'index': 'ID'})
        outfile = str(uuid.uuid4())
        outfolder = folder + '/' + subset
        check_or_create_dir(outfolder)
        df.to_csv(outfolder + '/' + outfile + '.csv', index=False)


def gen_background_performance(df, target, data_cols,
                               interval=10,
                               max_cols=1000):
    """ Given a dataframe, a target variable, data columns, a max number of
    columns to fit, and an interval to increase the sample of columns. Fits
    the data over the intervals and returns the oob R2 scores for each fit
    as a dict of {n_features: R2, ...}
    """
    oob_scores = {}
    for k in tqdm(range(10, max_cols + interval, interval)):
        selected_cols = sample_data_cols(data_cols, k)
        X, y = get_X_y(df, target, selected_cols)
        rf = fit_RF(X, y)
        oob_scores[k] = rf.oob_score_
    return oob_scores


def save_oob_scores(oob_scores, folder):
    """ Given a dict of oob_scores and a folder path, dump a csv of the
    oob_scores into the folder """
    d = oob_scores
    scores = pd.DataFrame(list(d.items()), columns=['n_features', 'R2'])
    scores = scores.sort_values('n_features').reset_index(drop=True)
    outfile = str(uuid.uuid4())
    scores.to_csv(folder + '/' + outfile + '.csv', index=False)


def get_gene_list_intersect(gene_list, data_cols):
    """ return the intersection between the current gene list HGNC symbols and
    the columns in the dataset. return a second list, `missing` for any genes
    that are missing.
    """
    intersect = [x for x in gene_list if x in data_cols]
    missing = [x for x in gene_list if x not in data_cols]
    return intersect, missing


def standardize_gmt(gmt):
    """ Takes a loaded list from a .gmt file and reformats it, if necessary,
    so that the html id is always at index 0 and the description is at index 1
    """
    if 'http' in gmt[0][1]:
        gmt_standard = [[x[1]] + [x[0]] + x[2:] for x in gmt]
    else:
        gmt_standard = gmt
    return gmt_standard


def read_gmt(path):
    """ given a filepath, reads the gmt or txt file at that location, returning
    a list that can be used in the scripts
    """
    if os.path.isfile(path):
        gmt = []
        with open(path) as f:
            rd = csv.reader(f, delimiter="\t", quotechar='"')
            for row in rd:
                gmt.append(row)
        gmt = standardize_gmt(gmt)
        gmt_suffix = path.split('/')[-1][:-4]
        return gmt, gmt_suffix
    else:
        files = os.listdir(path)
        gmt = []
        gmt_suffix = path.split('/')[-2]
        for f in files:
            print(f)
            with open(path + f) as fd:
                rd = csv.reader(fd, delimiter="\t", quotechar='"')
                gene_list = []
                for row in rd:
                    gene_list.append(row[0])
                gmt.append([f, 'user defined'] + gene_list)
        return gmt, gmt_suffix


def build_pairs_list(data, filter_column, value, split):
    data_f = data[data[filter_column] == value]
    samples = [(x, y) for x, y in
               zip(data_f['ID'].tolist(), data_f[split].tolist())]
    samples = sorted(samples, key=lambda x: x[1])
    pairs_list = [x for x in combinations(samples, 2)]
    return pairs_list
