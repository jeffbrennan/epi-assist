import scipy.stats as st
import pandas as pd
import numpy as np

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
    y_var = args[0]
    x_vars = args[1:]

    return sub_result

# TODO: add ability to import datasets and refer to them
def calc_selector(command, args, options):
    if command[0] == 'sum':
        user_result = sum(args)
        print(user_result)

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
