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
def calc_selector(user_query):
    if user_query[0] == 'sum':
        nums = [float(i) for i in user_query[1:]]
        user_result = sum(nums)
        print(user_result)

    elif user_query[0] == 'probcalc':
        result, output_text = prob_calc(user_query)
        print(output_text)

def main():
    user_query = [x for x in input().split()]
    calc_selector(user_query)

if __name__ == "__main__":
    main()
