import math
import scipy.stats as st
import numpy as np
import pandas as pd


def directadjustment(round_val):

    print('======== DIRECT ADJUSTMENT =========')
    bins = int(input('Enter number of groups: '))
    popFactor = int(input('Enter population factor eg 1000, 10000...'))

    loc1Rates = [float(x) for x in input('Enter location 1 rates separated by a space:\n').split()]
    loc2Rates = [float(x) for x in input('Enter location 2 rates separated by a space:\n').split()]
    standardPop = [float(x) for x in input('Enter standard populations separated by a space:\n').split()]

    totalStandardPop = sum(standardPop)

    # multiplies the condition rates by the population in the standard group
    # saves as list, then sums and creates rate by comparing expected to total
    expectedPop1 = [(loc1Rates[i] * standardPop[i]) for i in range(bins)]
    totalExpected1 = sum(expectedPop1)
    adj_rate_1 = (totalExpected1 / totalStandardPop) * popFactor

    expectedpop2 = [(loc2Rates[i] * standardPop[i]) for i in range(bins)]
    totalexpected2 = sum(expectedpop2)
    adj_rate_2 = (totalexpected2 / totalStandardPop) * popFactor

    popFactorText = ' per ' + str(popFactor) + ' persons'

    results = [adj_rate_1, adj_rate_2]
    output = [str(round(result, round_val)) for result in results]

    print('Location 1: ' + output[0] + popFactorText +
          ' | Location 2: ' + output[1] + popFactorText)

    return results


# Zscore calculations - 1 #
def zscore_calcs(round_val):

    print('================== Z SCORE FINDER ================== ')
    print('1 - Z Score | 2 - Z Score -> Area | 3 - Area -> Z Score | ' +
          '4 - Obs given percentile')

    typeChoice = str(input())
    zscore_chooser(round_val, typeChoice)


# Function navigates to the appropriate zscore related calculator based on user choice
def zscore_chooser(round_val, choice):

    if choice == '1':
        zscore_value(round_val)
    elif choice == '2':
        zscore_to_area(round_val, None)
    elif choice == '3':
        area_to_zscore(round_val, None)
    elif choice == '4':
        zscore_observation(round_val)


# Returns z score given observed value, pop mean, and standard deviation
def zscore_value(round_val):

    Obs = float(input('Enter the observation: '))
    Mean = float(input('Enter the mean: '))
    SD = float(input('Enter SD: '))

    zResult = (Obs - Mean) / SD
    zResult = float(round(zResult, round_val))

    print('Z Score: ' + str(zResult))

    zscore_to_area(round_val, zResult)

    return zResult


# Converts a given z score to its corresponding percentile, relates to zscore_value function
def zscore_to_area(round_val, calc_zscore):

    # Uses z score from zscore_value, otherwise prompts user input
    if calc_zscore:
        zscore = calc_zscore
    else:
        zscore = round(float(input('Enter your desired Z score (returns area): ')), 2)

    areacalc = st.norm.cdf(zscore)
    inversearea = 1 - areacalc
    areapercent = areacalc * 100

    area_result = round(areacalc, round_val)
    inverse_result = round(inversearea, round_val)
    area_pct_result = round(areapercent, round_val)

    print('Area: ' + str(area_result) + ' (' + str(area_pct_result) + '%)' +
          ' | Remaining area: ' + str(inverse_result))

    return areacalc


def area_to_zscore(round_val, percentile):

    if percentile:
        area = percentile
    else:
        area = float(input('Enter your desired area: '))

    # have to include output in if statement otherwise returns unboundlocalerror
    if area > 1 and area <= 100:
        zResult = st.norm.ppf(area / 100)
        z_output = round(zResult, round_val)
        print('Z score: ' + str(z_output))
        return zResult

    elif area < 1:
        zResult = st.norm.ppf(area)
        z_output = round(zResult, round_val)
        print('Z score: ' + str(z_output))
        return zResult

    elif area > 100:
        print('Area must be between 0 and 1 - restarting...')
        area_to_zscore(round_val, None)


def zscore_observation(round_val):  # returns observation at a given percentile

    percentile = float(input('Enter percentile: '))
    SD = float(input('Enter SD: '))
    Mean = float(input('Enter mean: '))

    zResult = area_to_zscore(round_val, percentile)

    observation = (zResult * SD) + Mean
    observation_output = round(observation, round_val)

    print(str(percentile) + ' percentile observation: ' + str(observation_output))
    return observation


# Bayes - 2 #
def bayes(round_val):
    print('====== BAYES CALC ========')
    print('Enter values as a proportion')

    prevalence = float(input('Enter disease prevalence: '))
    disease_acc = float(input('Enter accuracy for disease test: '))
    complement_acc = float(input('Enter accuracy for complement test: '))

    numerator = (prevalence * disease_acc)
    denominator = (numerator + ((1 - prevalence) * complement_acc))

    result = numerator / denominator
    result_output = round(result, round_val)

    print(result_output)
    return result


# TODO: add life table adjustments etc
# def disease_occurence(round_val):


def twobytwo(round_val):

    print('====== 2x2 TABLE SOLVER ======')

    print('1 - table empty | 2 - table full')
    typeChoice = str(input())

    if typeChoice == '1':
        print('Enter sensitivity and specificity or PPV and NPV')
        sensitivity = input('Enter sensitivity: ')
        specificity = input('Enter specificity: ')

        PPV = input('Enter PPV: ')
        NPV = input('Enter NPV: ')

        if sensitivity and specificity:

            sensitivity = float(sensitivity)
            specificity = float(specificity)

            popprev = float(input('Enter prevalence of population: '))
            popsize = float(input('Enter size of population: '))

            AC = popprev * popsize
            A = sensitivity * AC
            C = AC - A

            BD = popsize - AC
            D = BD * specificity
            B = BD - D

        elif PPV and NPV:

            PPV = float(PPV)
            NPV = float(NPV)

            popsize = int(input('Enter size of population: '))
            totalposresults = int(input('Enter number of positive test results: '))

            AB = totalposresults
            A = totalposresults * PPV
            B = AB - A

            CD = popsize - totalposresults
            D = CD * NPV
            C = CD - D

    elif typeChoice == '2':

        A = int(input('Enter value in cell A: '))
        B = int(input('Enter value in cell B: '))
        C = int(input('Enter value in cell C: '))
        D = int(input('Enter value in cell D: '))

    # calls helper functions to calculate basic values
    sensitivity, specificity = sens_spec(A, B, C, D)
    PPV, NPV = PPV_NPV(A, B, C, D)
    RR, OR = RR_OR(A, B, C, D)

    # Exposed incidence, nonexposed incidence, population incidence, adjusted rate, percent adjusted rate
    eIncidence, nonEIncidence, popIncidence, AR, PAR = incidence_2x2(A, B, C, D)

    errorChoice = str(input('Does the table have information bias (1: yes |2: no) '))
    if errorChoice == '1':
        errorRR, errorOR = error_2x2(A, B, C, D)

        errorRR = round(errorRR, round_val)
        errorOR = round(errorOR, round_val)
    else:
        errorRR = RR
        errorOR = OR

    col1 = [A, C]
    col2 = [B, D]

    tableData = {'+': col1, '-': col2}
    twoTable = pd.DataFrame(tableData, index=['Exposed', 'Nonexposed'])

    results = [PPV, NPV, sensitivity, specificity, OR, errorOR, RR, errorRR, eIncidence,
               nonEIncidence, popIncidence, AR, PAR, A, B, C, D]

    output = [str(round(result, round_val)) for result in results]

    print('======================== RESULTS ========================')

    print(twoTable)

    print('\n' + 'Positive predictive value (PPV): ' + output[0])
    print('Negative predictive value (NPV): ' + output[1])
    print('Sensitivity: ' + output[2])
    print('Specificity: ' + output[3])
    print('Odds Ratio (case control): ' + output[4])
    print('Error adjusted OR: ' + output[5])
    print('Risk Ratio (cohort): ' + output[6])
    print('Error adjusted RR: ' + output[7])
    print('Exposed Incidence: ' + output[8])
    print('Unexposed Incidence: ' + output[9])
    print('Total population incidence: ' + output[10])
    print('Attributable risk: ' + output[11])
    print('Percent attributable risk: ' + output[12])

    return results


def RR_OR(A, B, C, D):

    RR = (A / (A + B)) / (C / (C + D))
    OR = (A * D) / (B * C)

    return RR, OR


def sens_spec(A, B, C, D):
    sensitivity = A / (A + C)
    specificity = D / (B + D)

    return sensitivity, specificity


def PPV_NPV(A, B, C, D):
    PPV = A / (A + B)
    NPV = D / (C + D)

    return PPV, NPV


def incidence_2x2(A, B, C, D):
    totalPop = A + B + C + D

    eIncidence = (A / (A + B))
    nonEIncidence = (C / (C + D))
    popIncidence = (A + C) / totalPop

    AR = eIncidence - nonEIncidence
    PAR = ((popIncidence - nonEIncidence) / popIncidence) * 100

    return eIncidence, nonEIncidence, popIncidence, AR, PAR


def error_2x2(A, B, C, D):

    errorType = str(input('Enter error type (1: nondifferential | 2: differential) '))
    errorRate = float(input('Enter misclassification rate: '))
    errorDirection = str(input('Enter error direction (1: Error->No error |' +
                               '2: No error->Error): '))

    if errorType == '1' and errorDirection == '1':
        A -= (A * errorRate)
        B -= (B * errorRate)
        C += (A * errorRate)
        D += (B * errorRate)

    if errorType == '1' and errorDirection == '2':
        A += (C * errorRate)
        B += (D * errorRate)
        C -= (C * errorRate)
        D -= (D * errorRate)

    if errorType == '2':
        groupSelector = str(input('Enter differential group (1: cases | 2: controls) '))

        if groupSelector == '1' and errorDirection == '1':
            A -= (A * errorRate)
            C += (A * errorRate)

        elif groupSelector == '1' and errorDirection == '2':
            A += (C * errorRate)
            C -= (C * errorRate)

        elif groupSelector == '2' and errorDirection == '1':
            B -= (B * errorRate)
            D += (B * errorRate)

        elif groupSelector == '2' and errorDirection == '2':
            B += (D * errorRate)
            D -= (D * errorRate)

    errorRR, errorOR = RR_OR(A, B, C, D)
    return errorRR, errorOR


# TODO: add dataset import functionality (unrealistic to manually enter data for histogram etc)
def histogram_feat(round_val):
    print('Enter histogram data set as a list separated by spaces')
    data = [float(x) for x in input().split()]

    kurt_calc = st.kurtosis(data)
    skew_calc = st.skew(data)
    shapiro_calc = st.shapiro(data)
    mean = np.mean(data)
    median = np.median(data)

# This gross
    if mean > median:
        skew_text = 'Histogram is right-skewed'
    elif mean < median:
        skew_text = 'Histogram is left-skewed'
    elif mean == median:
        skew_text = 'Histogram is symmetric'

    if kurt_calc > 0:
        kurt_text = 'Data is leptokurtic'
    elif kurt_calc < 0:
        kurt_text = 'Data is platykurtic'

    if shapiro_calc[1] < 0.05:
        shapiro_text = 'Data is not normal'
    elif shapiro_calc[1] > 0.05:
        shapiro_text = 'Data is normal'

    skew_out = round(skew_calc, round_val)
    kurt_out = round(kurt_calc, round_val)
    shapiro_out = round(shapiro_calc[1], round_val)

    print('Skewness: ' + str(skew_out) + ' | ' + skew_text)
    print('Kurtosis: ' + str(kurt_out) + ' | ' + kurt_text)
    print('Shapiro-Wilk Test: ' + str(shapiro_out) + ' | ' + shapiro_text)

    return skew_calc, kurt_calc, shapiro_calc


def estimation_calcs(round_val):
    print('Enter desired estimation: (1: Population mean | 2: Sample mean | 3: Variation) ')
    typeChoice = str(input())

    if typeChoice == '1':
        pop_calc(round_val)
    elif typeChoice == '2':
        sample_calc(round_val)
    elif typeChoice == '3':
        var_calc(round_val)


def pop_calc(round_val):
    sample_size = int(input('Enter number in sample: '))

    print('Enter observations for all N in sample separated by a space')
    observations = [float(x) for x in input().split()]

    if len(observations) == sample_size:

        mean, variance, SD, SE = list_calcs(observations, sample_size)

        results = [mean, variance, SD, SE]
        output = [str(round(result, round_val)) for result in results]
        print('Population mean: ' + output[0])
        print('Population variance: ' + output[1])
        print('Population SD: ' + output[2])

        return results

    else:
        print('Entered incorrect number of observations - restarting...')
        estimation_calcs(round_val)


def sample_calc(round_val):
    sample_size = int(input('Enter number in sample: '))

    print('Enter observations for all N in sample')
    observations = [float(x) for x in input().split()]

    if len(observations) == sample_size:
        mean, variance, SD, SE = list_calcs(observations, sample_size)

        results = [mean, variance, SD, SE]
        output = [str(round(result, round_val)) for result in results]

        print('Sample mean: ' + output[0])
        print('Sample variance: ' + output[1])
        print('Sample SD: ' + output[2])
        print('Sample SE: ' + output[3])

        return results

    else:
        print('Entered incorrect number of observations - restarting...')
        estimation_calcs(round_val)


def var_calc(round_val):
    s_var = float(input('Enter sample variance: '))
    sample_size = int(input('Enter the sample size: '))
    CI = float(input('Enter desired confidence (as decimal): '))

    CI_tails = (1 - CI) / 2
    chi_low = st.chi2.ppf(1 - CI_tails, sample_size - 1)
    chi_high = st.chi2.ppf(CI_tails, sample_size - 1)

    var_low = ((sample_size - 1) * s_var) / chi_low
    var_high = ((sample_size - 1) * s_var) / chi_high

    CI_PCT = CI * 100

    results = [CI_PCT, var_low, var_high]
    output = [str(round(result, round_val)) for result in results]
    print('We are ' + output[0] + '% confident that the variance is between ' + output[1] + ' and ' + output[2])


def list_calcs(numberList, sample_size):

    mean = sum(numberList) / sample_size

    squareDif = [(i - mean)**2 for i in numberList]
    variance = sum(squareDif) / sample_size

    SD = math.sqrt(variance)
    SE = variance / (math.sqrt(sample_size))

    return mean, variance, SD, SE


def discrete_calcs(round_val):
    discrete_choice = str(input('Enter desired discrete calc (1: Binomial | 2: Poisson) '))

    if discrete_choice == '1':
        binomial_calc(round_val)
    elif discrete_choice == '2':
        poisson_calc(round_val)


def binomial_calc(round_val):

    n = int(input('Enter n: '))
    p = float(input('Enter p: '))

    num_events = int(input('Enter desired number of events [int]: '))

    prob_sum = 0
    for i in range(num_events):
        choose = math.factorial(n) / (math.factorial(i) * math.factorial(n - i))
        prob = (choose * (p ** i) * (1 - p) ** (n - i))
        print('P(X = ' + str(i) + ') = ' + str(prob))

        prob_sum += prob

    inverse_prob = 1 - prob_sum
    results = [prob_sum, inverse_prob]

    output = [str(round(result, round_val)) for result in results]

    print('Probability: ' + output[0] + ' | ' + output[1])

    return results


def poisson_calc(round_val):
    poisson_rate = float(input('Enter event rate: '))
    poisson_time = float(input('Enter desired event time: '))
    poisson_pop = int(input('Enter size of population: '))
    num_events = int(input('Enter desired number of events [int]: '))

    expected_num = poisson_rate * poisson_time * poisson_pop

    prob_sum = 0
    for i in range(num_events):
        prob = (((math.e) ** -expected_num) * (expected_num ** i)) / math.factorial(i)
        print('P(X = ' + str(i) + ') = ' + str(prob))

        prob_sum += prob

    final_prob = 1 - prob_sum

    results = [expected_num, num_events, final_prob]
    output = [str(round(result, round_val)) for result in results]

    print('Expected number of events: ' + output[0])
    print('Probability of ' + output[1] + ' or more events: ' + output[2])

    return results


def hypothesis_calcs(round_val):
    print('================== HYPOTHESIS TESTER ================== ')

    data_type = str(input('Enter the data type: (1: Numerical | 2: Categorical) '))
    # Numerical Tests (t test etc)
    if data_type == '1':
        num_test_choice = str(input('1: Z-test | 2: T-Test| 3: Correlation: '))

        tailChoice = str(input('1: One-tailed test | 2: Two-tailed test: '))

        if tailChoice == '1':
            print('1: HA > H0 | 2: HA < H0')
            tailSide = str(input())

            if num_test_choice == '1':
                hypothesis_zTest(round_val, tailChoice, tailSide)
            elif num_test_choice == '2':
                hypothesis_tTest(round_val, tailChoice, tailSide)
            elif num_test_choice == '3':
                hypothesis_chisquare(round_val)

        elif tailChoice == '2':
            if num_test_choice == '1':
                hypothesis_zTest(round_val, tailChoice, None)
            elif num_test_choice == '2':
                hypothesis_tTest(round_val, tailChoice, None)
            elif num_test_choice == '3':
                hypothesis_correlation(round_val)

    # Categorical tests (proportion, difference of proportion etc)
    elif data_type == '2':
        cat_test_choice = str(input('Enter desired test: (1: Chi-Squared | 2: Fishers exact) '))
        if cat_test_choice == '1':
            hypothesis_chisquare(round_val)
        elif cat_test_choice == '2':
            hypothesis_fishers(round_val)


def hypothesis_zTest(round_val, tailChoice, tailSide):

    nullValue = float(input('Enter null hypothesis value: '))
    sample_size = float(input('Enter sample size: '))
    sampleMean = float(input('Enter sample mean: '))
    sampleSD = float(input('Enter sample SD: '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

    popSD = sampleSD / math.sqrt(sample_size)
    zValue = (sampleMean - nullValue) / popSD
    zArea = zscore_to_area(round_val, zValue)

    if tailChoice == '1':
        pValue = round(1 - float(zArea), int(round_val))
    elif tailChoice == '2':
        pValue = round(2 * (1 - float(zArea)), int(round_val))

    print('P-value: ' + str(pValue))

    return zValue, zArea, pValue


def hypothesis_tTest(round_val, tailChoice, tailSide):

    tType = str(input('Enter desired type of t-test: (1: One Sample | 2: 2-Sample | 3: Indepdendent 2-Sample) '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

# TODO: replace with functions
    if tType == '1' or tType == '2':
        nullValue = float(input('Enter null hypothesis value: '))
        sample_size = float(input('Enter sample size: '))
        sampleMean = float(input('Enter sample mean: '))
        sampleSD = float(input('Enter sample SD: '))

        df = sample_size - 1

        tScore = (sampleMean - nullValue) / (sampleSD / math.sqrt(sample_size))
        criticalT = math.fabs(st.t.ppf(alphaValue, df))

        print('T-Score: ' + str(tScore))
        print('Critical T-Score: ', criticalT)

    elif tType == '3':
        xbar1 = float(input('Enter sample mean 1: '))
        s1 = float(input('Enter sample SD 1: '))
        n1 = int(input('Enter sample size 1: '))

        xbar2 = float(input('Enter sample mean 2: '))
        s2 = float(input('Enter sample SD 2: '))
        n2 = int(input('Enter sample size 2: '))

        tScore_numerator = xbar1 - xbar2
        tScore_denominator = math.sqrt((s1**2 / n1) + (s2**2 / n2))
        tScore = tScore_numerator / tScore_denominator

        df = hypothesis_df(s1, n1, s2, n2)
        criticalT = math.fabs(st.t.ppf(alphaValue, df))

        print(round(tScore, round_val))
        print(round(criticalT, round_val))

    if tScore < criticalT:
        print('FTR Null')
    elif tScore > criticalT:
        print('Reject Null')


def table_maker(table, rows, i):
    while True:
        input_text = 'Row #' + str(i + 1) + ' | Enter row values separated by a space: '
        new_row = [int(x) for x in input(input_text).split()]

        if i > 0:
            if len(new_row) == len(table[i - 1]):
                return new_row
                break
            else:
                print('\n Number of row items entered incorrectly. Retry input \n')
        else:
            return new_row
            break


def hypothesis_chisquare(round_val):
    table = []
    rows = int(input('Enter number of rows: '))

    for i in range(rows):
        table.append(table_maker(table, rows, i))

    chi, p, dof, expected = st.chi2_contingency(table)
    results = [chi, p, dof]
    output = [str(round(result, round_val)) for result in results]

    print('Chi Square value: ' + output[0])
    print('P-value: ' + (output[1]))
    print('Degrees of freedom: ' + output[2])

    return results


def hypothesis_fishers(round_val):

    table = []
    rows = int(input('Enter number of rows: '))

    for i in range(rows):
        table.append(table_maker(table, rows, i))

    odds, p = st.fisher_exact(table)
    results = [odds, p]
    output = [str(round(result, round_val)) for result in results]

    print('Odds ratio: ' + output[0])
    print('P-value: ' + output[1])

    return results


def hypothesis_correlation(round_val):

    print('Enter group 1 numbers (continuous) separated by a space:')
    group_1_nums = [float(x) for x in input().split()]

    print('Enter continuous group 2 numbers (continuous) separated by a space:')
    group_2_nums = [float(x) for x in input().split()]

    data = [group_1_nums, group_2_nums]

    normal_1 = normality_checker(data[0])
    normal_2 = normality_checker(data[1])

    if not normal_1 or not normal_2:
        print('Data is not normal')

    elif normal_1 and normal_2:
        calc_choice = str(input('|1: Spearman|2: Pearson| '))

        if calc_choice == '1':
            r, p = st.spearmanr(data[0], data[1])
        elif calc_choice == '2':
            r, p = st.pearsonr(data[0], data[1])

        results = [r, p]
        output = [str(round(result, round_val)) for result in results]

        print('R-coefficient: ' + output[0])
        print('P-value: ' + output[1])

        return results


def hypothesis_df(s1, n1, s2, n2):

    # independent two sample test
    if s1 and s2 and n1 and n2:
        numerator = ((s1**2 / n1) + (s2**2 / n2))**2
        denominator = ((s1**2 / n1)**2 / (n1 - 1) + (s2**2 / n2)**2 / (n2 - 1))

        df = int(math.floor(numerator / denominator))

    print('Degrees of Freedom: ' + str(df))
    return df


def normality_checker(data):
    t, p = st.shapiro(data)

    # TODO: test this

    # if p > 0.05:
    #     return True
    return p > 0.05


def confidence_interval(round_val, calc_SE):

    pointEstimate = float(input('Enter the point estimate: '))

    if calc_SE:
        SE = calc_SE
    else:
        SE_Choice = str(input('Enter SE type ( 1: Z Score | 2: Proportions | 3: Diff Proportions ) '))
        if SE_Choice == '1':
            SD = float(input('Enter the standard deviation: '))
            sample_size = float(input('Enter the sample size: '))
            SE = SD / math.sqrt(sample_size)
        elif SE_Choice == '2':
            phat = float(input('Enter phat value: '))
            sample_size = float(input('Enter the sample size: '))
            SE = math.sqrt((phat * (1 - phat) / sample_size))
        elif SE_Choice == '3':
            phat_1 = float(input('Enter first phat value: '))
            sample_size_1 = float(input('Enter first sample size: '))

            phat_2 = float(input('Enter second phat value: '))
            sample_size_2 = float(input('Enter first sample size: '))

            SE = math.sqrt((phat_1 * (1 - phat_1) / sample_size_1) +
                           (phat_2 * (1 - phat_2) / sample_size_2))

    print('Enter desired % confidence (1: 90 | 2: 95 | 3: 99 | 4: Custom ) ')

    # TODO: replace this with a dictionary or automated function call

    CI_Choice = str(input())

    if CI_Choice == '1':
        z = 1.645
    elif CI_Choice == '2':
        z = 1.96
    elif CI_Choice == '3':
        z = 2.58
    elif CI_Choice == '4':
        z = zscore_to_area(round_val, None)

    CIlow = pointEstimate - (z * SE)
    CIhigh = pointEstimate + (z * SE)
    marginError = (CIhigh - CIlow) / 2

    results = [CIlow, CIhigh, marginError]
    output = [str(round(result, round_val)) for result in results]

    print('There is a 95% chance that the true population [parameter] ' +
          'lies between the interval: (' + output[0] + ',' +
          output[1] + ') +/- ' + output[2])

    return results


def nonparametric_calcs(round_val):
    print('Enter group 1 numbers (continuous) separated by a space:')
    group_1_nums = [float(x) for x in input().split()]

    print('Enter continuous group 2 numbers (continuous) separated by a space:')
    group_2_nums = [float(x) for x in input().split()]

    data = [group_1_nums, group_2_nums]

    normal_1 = normality_checker(data[0])
    normal_2 = normality_checker(data[1])

    if not normal_1 or not normal_2:
        print('Data is not normal')
    elif normal_1 and normal_2:
        test_choice = str(input('Are the observations dependent? (1: Yes | 2: No) '))
        if test_choice == '1':
            statistic, p = st.wilcoxon(data[0], data[1])
        elif test_choice == '2':
            statistic, p = st.ranksums(data[0], data[1])

    results = [statistic, p]
    output = [str(round(result, round_val)) for result in results]

    print('Statistic: ' + output[0])
    print('P-value: ' + output[1])

    return results


def calcselection(start_calc, i=4):
    chooser = {'1': zscore_calcs, '2': bayes, '3': directadjustment, '4': twobytwo,
               '5': histogram_feat, '6': estimation_calcs, '7': discrete_calcs,
               '8': hypothesis_calcs, '9': confidence_interval, '10': nonparametric_calcs}
    chooser[start_calc](i)


def main():
    print('Enter the calc you would like to use')
    print('1 - zscore')
    print('2 - bayes')
    print('3 - direct adjustment')
    print('4 - 2X2 Solver')
    print('5 - Histogram')
    print('6 - Estimations')
    print('7 - Discrete distributions')
    print('8 - Hypothesis Testing')
    print('9 - Confidence Interval')
    print('10 - Nonparametric Tests')

    start_calc = str(input())
    round_val = int(input('Round result to how many decimal places? '))
    calcselection(start_calc, round_val)


if __name__ == "__main__":
    main()
