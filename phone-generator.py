from random import randint
import pandas as pd

case_phones = pd.read_csv('phone-list.csv').squeeze()
control_phones = []

for number in case_phones:
    control_phones.extend([number[:8] + str(randint(0000, 9999))])

phone_df = pd.DataFrame(control_phones, columns=["Number"])
phone_df.to_csv('control-numbers.csv')
