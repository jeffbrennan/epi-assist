import scipy.stats as st
import pandas as pd
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

def prob_calc_normal(query):
    mean = float(query[2])
    sd = float(query[3])
    option = query[4]
    obs = float(query[5])

    z_result = (obs - mean) / sd
    if option == 'atleast':
        area = 1 - (st.norm.cdf(z_result))
    elif option == 'atmost':
        area = st.norm.cdf(z_result)

    text = ('cdf Method: P(X>=' + str(obs) + ')=' + str(area))
    return area, text

def prob_calc(query):
    if query[1] == 'n':
        sub_result = prob_calc_normal(query)

    return sub_result

def regress(args, options):

    regr = linear_model.LinearRegression()

    if len(args) > 2:  # multiple regression
        y_var = args[0]
        x_vars = args[1:]
    elif len(args) == 2:  # simple linear regression
        y_var = args[0]
        y_mid = len(y_var) // 2
        y_train = y_var[:y_mid]
        y_test = y_var[y_mid:]

        x_var = args[1]
        x_mid = len(x_var) // 2
        x_train = x_var[:x_mid]
        x_test = x_var[x_mid:]

        regr.fit(x_train, y_train)
        y_predict = regr.predict(x_test)

    elif len(args) < 2:
        print('Regression requires two variables')

    coef = regr.coef_
    MSE = mean_squared_error(y_test, y_predict)
    r2 = r2_score(y_test, y_predict)

    sub_result = [coef, MSE, r2]

    return sub_result

def dataset_import(args, options):
    data_path = args[0]
    dataset = pd.read_csv(data_path)

# TODO: add ability to complete several functions before exit (eg import dataset then do stuff with it)
def calc_selector(command, args, options):
    if command[0] == 'sum':
        user_result = sum(args)
        print(user_result)
    elif command[0] == 'import':
        dataset_import(args, options)
    elif command[0] == 'regress':
        user_result = regress(args, options)

    elif command[0] == 'probcalc':
        result, output_text = prob_calc(command)
        print(output_text)

def main():
    args = []
    options = []

    user_query = [x for x in input().split()]
    command = user_query[0]

    try:
        comma_loc = [i for i, s in enumerate(user_query) if ',' in s][0]
        args.extend(user_query[1:comma_loc])
        args.extend(user_query[comma_loc][:-1])
        options.extend(user_query[comma_loc + 1:])
    except IndexError:
        args.extend(user_query[1:])
        options = None

    calc_selector(command, args, options)

if __name__ == "__main__":
    main()
