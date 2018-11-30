import math
import scipy.stats as st
import numpy as np
import pandas as pd

def directadjustment(roundingValue):

    print('======== DIRECT ADJUSTMENT =========')
    bins = int(input('Enter number of groups: '))
    popFactor = int(input('Enter population factor eg 1000, 10000...'))

    print('Enter location 1 rates separated by a space:')
    loc1Rates = [float(x) for x in input().split()]

    print('Enter location 2 rates separated by a space:')
    loc2Rates = [float(x) for x in input().split()]

    print('Enter Standard populations separated by a space:')  # might help to have explanation on this point
    standardPop = [float(x) for x in input().split()]

    totalStandardPop = sum(standardPop)
    # multiples the condition rates by the population in the standard group
    # saves as list, then sums and creates rate by comparing expected to total
    expectedPop1 = [(loc1Rates[i] * standardPop[i]) for i in range(bins)]
    totalExpected1 = sum(expectedPop1)
    adjustedRate1 = (totalExpected1 / totalStandardPop) * popFactor

    expectedpop2 = [(loc2Rates[i] * standardPop[i]) for i in range(bins)]
    totalexpected2 = sum(expectedpop2)
    adjustedrate2 = (totalexpected2 / totalStandardPop) * popFactor

    popFactorText = ' per ' + str(popFactor) + ' persons'

    print('Location 1: ' + str(adjustedRate1) + popFactorText +
          ' | Location 2: ' + str(adjustedrate2) + popFactorText)

    return adjustedRate1, adjustedrate2

# Zscore calculations - 1 #
def zscore_calcs(roundingValue):

    print('================== Z SCORE FINDER ================== ')
    print('1 - Z Score | 2 - Z Score -> Area | 3 - Area -> Z Score | ' +
          '4 - Obs given percentile')

    typeChoice = str(input())
    zscore_chooser(roundingValue, typeChoice)

# Function navigates to the appropriate zscore related calculator based on user choice
def zscore_chooser(roundingValue, choice):

    if choice == '1':
        zscore_value(roundingValue)
    elif choice == '2':
        zscore_toarea(roundingValue, None)
    elif choice == '3':
        area_tozscore(roundingValue, None)
    elif choice == '4':
        zscore_observation(roundingValue)

# Returns z score given observed value, pop mean, and standard deviation
def zscore_value(roundingValue):

    Obs = float(input('Enter the observation: '))
    Mean = float(input('Enter the mean: '))
    SD = float(input('Enter SD: '))

    zResult = (Obs - Mean) / SD
    zResult = float(round(zResult, roundingValue))

    print('Z Score: ' + str(zResult))

    zscore_toarea(roundingValue, zResult)

    return zResult

# Converts a given z score to its corresponding percentile, relates to zscore_value function
def zscore_toarea(roundingValue, calc_zscore):
    if calc_zscore:
        zscore = calc_zscore
    else:
        zscore = round(float(input('Enter your desired Z score (returns area): ')), 2)

    areacalc = st.norm.cdf(zscore)
    inversearea = 1 - areacalc
    areapercent = areacalc * 100
    areacalc = round(areacalc, roundingValue)
    inversearea = round(inversearea, roundingValue)
    areapercent = round(areapercent, roundingValue)

    print('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' +
          '|Remaining area: ' + str(inversearea))

    return areacalc

# What does this do
def area_tozscore(roundingValue, percentile):
    area = input('Enter your desired area, blank if auto (returns zscore): ')
    if area:
        area = float(area)
    else:
        area = percentile

    if area > 1:
        zResult = st.norm.ppf(percentile / 100)
    else:
        zResult = st.norm.ppf(percentile)

    print(zResult)
    return zResult

def zscore_observation(roundingValue):  # returns observation at a given percentile

    percentile = float(input('Enter percentile: '))
    SD = float(input('Enter SD: '))
    Mean = float(input('Enter mean: '))

    zResult = area_tozscore(roundingValue, percentile)

    observation = (zResult * SD) + Mean
    observation = round(observation, roundingValue)

    print(str(percentile) + ' percentile observation: ' + str(observation))
    return observation

# Bayes - 2 #
def bayes(roundingValue):
    print('====== BAYES CALC ========')
    print('Enter values as a proportion')

    prevalence = float(input('Enter disease prevalence: '))
    diseaseAccuracy = float(input('Enter accuracy for disease test: '))
    complementAccuracy = float(input('Enter accuracy for complement test: '))

    numerator = (prevalence * diseaseAccuracy)
    denominator = (numerator + ((1 - prevalence) * complementAccuracy))

    result = numerator / denominator
    result = round(result, roundingValue)

    print(result)

    return result

def twobytwo(roundingValue):

    print('====== 2x2 TABLE SOLVER ======')

    print('1 - table empty | 2 - table full')
    typeChoice = str(input())

    if typeChoice == '1':
        print('Enter sensitivity and specificity or PVP and PVN')
        sensitivity = input('Enter sensitivity: ')
        specificity = input('Enter specificity: ')

        PVP = input('Enter PVP: ')
        PVN = input('Enter PVN: ')

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

        elif PVP and PVN:

            PVP = float(PVP)
            PVN = float(PVN)

            popsize = int(input('Enter size of population: '))
            totalposresults = int(input('Enter number of positive test results: '))

            AB = totalposresults
            A = totalposresults * PVP
            B = AB - A

            CD = popsize - totalposresults
            D = CD * PVN
            C = CD - D

    elif typeChoice == '2':

        A = int(input('Enter value in cell A: '))
        B = int(input('Enter value in cell B: '))
        C = int(input('Enter value in cell C: '))
        D = int(input('Enter value in cell D: '))

    # calls helper functions to calculate basic values
    sensitivity, specificity = sens_spec(A, B, C, D)
    PVP, PVN = PVP_PVN(A, B, C, D)
    RR, OR = RR_OR(A, B, C, D)
    RR, OR = RR_OR(A, B, C, D)

    # Exposed incidence, nonexposed incidence, population incidence, adjusted rate, percent adjusted rate
    eIncidence, nonEIncidence, popIncidence, AR, PAR = incidence_2x2(A, B, C, D)

    results = [PVP, PVN, sensitivity, specificity, OR, RR, eIncidence,
               nonEIncidence, popIncidence, AR, PAR, A, B, C, D]

    output = [round(i, roundingValue) for i in results]

    errorChoice = str(input('Does the table have information bias (1 - y|2 - n): '))
    if errorChoice == '1':
        errorRR, errorOR = error_2x2(A, B, C, D)

        errorRR = round(errorRR, roundingValue)
        errorOR = round(errorOR, roundingValue)
    else:
        errorRR = 'N/A: No information bias declared'

    col1 = [A, C]
    col2 = [B, D]

    tableData = {'+': col1, '-': col2}
    twoTable = pd.DataFrame(tableData, index=['Exposed', 'Nonexposed'])

    print('======================== RESULTS ========================')

    print(twoTable)

    # consider converting to dictionary
    print('\n' + 'Predictive positive value (PVP): ' + str(output[0]))
    print('Predictive negative value (PVN): ' + str(output[1]))
    print('Sensitivity: ' + str(output[2]))
    print('Specificity: ' + str(output[3]))
    print('Odds Ratio (case control): ' + str(output[4]))
    print('Error adjusted OR: ' + str(errorOR))
    print('Risk Ratio (cohort): ' + str(output[5]))
    print('Error adjusted RR: ' + str(errorRR))
    print('Exposed Incidence: ' + str(output[6]))
    print('Unexposed Incidence: ' + str(output[7]))
    print('Total population incidence: ' + str(output[8]))
    print('Attributable risk: ' + str(output[9]))
    print('Percent attributable risk: ' + str(output[10]))

    return results

def RR_OR(A, B, C, D):

    RR = (A / (A + B)) / (C / (C + D))
    OR = (A * D) / (B * C)

    return RR, OR

def sens_spec(A, B, C, D):
    sensitivity = A / (A + C)
    specificity = D / (B + D)

    return sensitivity, specificity

def PVP_PVN(A, B, C, D):
    PVP = A / (A + B)
    PVN = D / (C + D)

    return PVP, PVN

def incidence_2x2(A, B, C, D):
    totalPop = A + B + C + D

    eIncidence = (A / (A + B))
    nonEIncidence = (C / (C + D))
    popIncidence = (A + C) / totalPop
    AR = eIncidence - nonEIncidence
    PAR = ((popIncidence - nonEIncidence) / popIncidence) * 100
    return eIncidence, nonEIncidence, popIncidence, AR, PAR

def error_2x2(A, B, C, D):

    errorType = str(input('Enter error type (1 - nondifferential |2 - differential): '))
    errorRate = float(input('Enter misclassification rate'))
    errorDirection = str(input('Enter error direction (1 - E->NE |2 - NE->E): '))

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
        groupSelector = str(input('Enter differential group (1 - cases|2- controls): '))

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

def histogram_feat(roundingValue):
    print('Enter histogram data set as a list')
    data = [float(x) for x in input().split()]

    kurtosisCalc = st.kurtosis(data)
    skewCalc = st.skew(data)
    shapiroCalc = st.shapiro(data)
    mean = np.mean(data)
    median = np.median(data)

    print('=============RESULT===========')

    if mean > median:
        print('Histogram is right-skewed')
    elif mean < median:
        print('Histogram is left-skewed')
    elif mean == median:
        print('Histogram is symmetric')

    # !! Need to check decisions
    if kurtosisCalc > 3:
        print('Data is leptokurtic')
    elif kurtosisCalc < 3:
        print('Data is platykurtic ')

    skewCalc = round(skewCalc, roundingValue)
    kurtosisCalc = round(kurtosisCalc, roundingValue)

    print('Skewness: ' + str(skewCalc))
    print('Kurtosis: ' + str(kurtosisCalc))
    print('Shapiro-Wilk Test: ' + str(shapiroCalc))

    return skewCalc, kurtosisCalc, shapiroCalc

def estimation(roundingValue):
    print('1 - Population parameters | 2 - Sample parameters')
    typeChoice = str(input())

    if typeChoice == '1':
        sampleSize = int(input('Enter number in sample: '))

        print('Enter observations for all N in sample separated by a space')
        sampleValues = [float(x) for x in input().split()]

        if len(sampleValues) == sampleSize:

            mean, variance, SD, SE = list_calcs(sampleValues, sampleSize)

            results = [mean, variance, SD, SE]
            output = [round(i, roundingValue) for i in results]
            print('Population mean: ' + str(output[0]))
            print('Population variance: ' + str(output[1]))
            print('Population SD: ' + str(output[2]))

            return results

        else:
            print('Entered incorrect number of observations - restarting...')
            estimation(roundingValue)

    elif typeChoice == '2':

        sampleSize = int(input('Enter number in sample: '))

        print('Enter observations for all N in sample')
        sampleValues = [float(x) for x in input().split()]

        if len(sampleValues) == sampleSize:

            mean, variance, SD, SE = list_calcs(sampleValues, sampleSize)

            results = [mean, variance, SD, SE]
            output = [round(i, roundingValue) for i in results]

            print('Sample mean: ' + str(output[0]))
            print('Sample variance: ' + str(output[1]))
            print('Sample SD: ' + str(output[2]))
            print('Sample SE: ' + str(output[3]))

            return results

        else:
            print('Entered incorrect number of observations - restarting...')
            estimation(roundingValue)

def list_calcs(numberList, sampleSize):

    mean = sum(numberList) / sampleSize

    squareDif = [(i - mean)**2 for i in numberList]
    variance = sum(squareDif) / sampleSize

    SD = math.sqrt(variance)
    SE = variance / (math.sqrt(sampleSize))

    return mean, variance, SD, SE

def binomial(roundingValue):

    n = int(input('Enter n: '))
    x = int(input('Enter x: '))
    p = float(input('Enter p: '))

    choose = math.factorial(n) / (math.factorial(x) * math.factorial(n - x))
    result = (choose * (p ** x) * (1 - p) ** (n - x))
    resultpercent = result * 100

    result = round(result, roundingValue)
    resultpercent = round(resultpercent, str(int(roundingValue) - 2))

    print('Probability: ' + result + ' | ' + resultpercent + '%')

    return result

def hypothesis_calcs(roundingValue):
    print('================== HYPOTHESIS TESTER ================== ')

    data_type = str(input('1: Numerical | 2: Categorical'))
    # Numerical Tests (t test etc)
    if data_type == '1':
        num_test_choice = str(input('|1: Z-test|2: T-Test|3: Correlation|'))

        tailChoice = str(input('1: One-tailed test | 2: Two-tailed test: '))

        if tailChoice == '1':
            print('Enter greater or less than null: 1 - greater (upper) | 2 - less than (lower)')
            tailSide = str(input())

            if num_test_choice == '1':
                hypothesis_zTest(roundingValue, tailChoice, tailSide)
            elif num_test_choice == '2':
                hypothesis_tTest(roundingValue, tailChoice, tailSide)
            elif num_test_choice == '3':
                hypothesis_chisquare(roundingValue)

        elif tailChoice == '2':
            if num_test_choice == '1':
                hypothesis_zTest(roundingValue, tailChoice, None)
            elif num_test_choice == '2':
                hypothesis_tTest(roundingValue, tailChoice, None)
            elif num_test_choice == '3':
                hypothesis_correlation(roundingValue)

    # Categorical tests (proportion, difference of proportion etc)
    elif data_type == '2':
        cat_test_choice = str(input('1: Chi-Squared | 2: Fishers exact: '))
        if cat_test_choice == '1':
            hypothesis_chisquare(roundingValue)
        elif cat_test_choice == '2':
            hypothesis_fishers(roundingValue)

def hypothesis_zTest(roundingValue, tailChoice, tailSide):

    nullValue = float(input('Enter null hypothesis value: '))
    sampleSize = float(input('Enter sample size: '))
    sampleMean = float(input('Enter sample mean: '))
    sampleSD = float(input('Enter sample SD: '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

    popSD = sampleSD / math.sqrt(sampleSize)

    zValue = (sampleMean - nullValue) / popSD

    zArea = zscore_toarea(roundingValue, zValue)

    if tailChoice == '1':
        pValue = round(1 - float(zArea), int(roundingValue))
    elif tailChoice == '2':
        pValue = round(2 * (1 - float(zArea)), int(roundingValue))
        print('P-value: ' + str(pValue))

    print(hypothesis_decision(pValue, alphaValue, tailChoice, tailSide))

def hypothesis_tTest(roundingValue, tailChoice, tailSide):

    tType = str(input('1: One Sample | 2: 2-Sample | 3: Indepdendent 2-Sample: '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

    if tType == '1' or tType == '2':

        nullValue = float(input('Enter null hypothesis value: '))
        sampleSize = float(input('Enter sample size: '))
        sampleMean = float(input('Enter sample mean: '))
        sampleSD = float(input('Enter sample SD: '))

        df = sampleSize - 1

        tScore = (sampleMean - nullValue) / (sampleSD / math.sqrt(sampleSize))
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

        print(round(tScore, roundingValue))
        print(round(criticalT, roundingValue))

    if tScore < criticalT:
        print('FTR Null')
    elif tScore > criticalT:
        print('Reject Null')

def dynamic_table_maker(table, rows, i):
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

def hypothesis_chisquare(roundingValue):
    table = []
    rows = int(input('Enter number of rows: '))

    for i in range(rows):
        table.append(dynamic_table_maker(table, rows, i))

    chi, p, dof, expected = st.chi2_contingency(table)
    results = [chi, p, dof]
    output = [round(i, roundingValue) for i in results]

    print('Chi Square value: ' + str(output[0]))
    print('P-value: ' + str(output[1]))
    print('Degrees of freedom: ' + str(output[2]))

    return results

def hypothesis_fishers(roundingValue):

    table = []
    rows = int(input('Enter number of rows: '))

    for i in range(rows):
        table.append(dynamic_table_maker(table, rows, i))

    odds, p = st.fisher_exact(table)
    results = [odds, p]
    output = [round(i, roundingValue) for i in results]

    print('Odds ratio: ' + str(output[0]))
    print('P-value: ' + str(output[1]))

    return results

def hypothesis_correlation(roundingValue):

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
        calc_choice = str(input('1: Spearman|2: Pearson|'))

        if calc_choice == '1':
            r, p = st.spearmanr(data[0], data[1])
        elif calc_choice == '2':
            r, p = st.pearsonr(data[0], data[1])

        results = [r, p]
        output = [round(i, roundingValue) for i in results]

        print('R-coefficient: ' + str(output[0]))
        print('P-value: ' + str(output[1]))

        return results
# Consider deleting and handle each decision individually (becomes difficult w/ proportion based tests)
def hypothesis_decision(testValue, alphaValue, tailChoice, tailSide):

    if tailChoice == '1':
        if tailSide == '1':
            if testValue > alphaValue:
                return(str(testValue) + ' > ' + str(alphaValue) + ': Reject null hypothesis')
            else:
                return(str(testValue) + ' < ' + str(alphaValue) + ': FTR null hypothesis')
        elif tailSide == '2':
            if testValue < alphaValue:
                return(str(testValue) + ' < ' + str(alphaValue) + ': Reject null hypothesis')
            else:
                return(str(testValue) + ' > ' + str(alphaValue) + ': FTR null hypothesis')
    elif tailChoice == '2':
        if abs(testValue) < abs(alphaValue):
            return(str(testValue) + ' < ' + str(alphaValue) + ': Reject null hypothesis')
        else:
            return(str(-testValue) + ' <= ' + str(alphaValue) + ' <= ' + str(testValue) + ': FTR null hypothesis')

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

    if p > 0.05:
        return True

def confidence_interval(roundingValue, calc_SE):

    pointEstimate = float(input('Enter the point estimate: '))

    if calc_SE:
        SE = calc_SE
    else:
        SE_Choice = str(input('Enter SE type|1: Z Score|2: Proportions|3: Diff Proportions: '))
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

    print('Enter desired % confidence: |1: 90|2: 95|3: 99|4: Custom|')
    CI_Choice = str(input())

    if CI_Choice == '1':
        z = 1.645
    elif CI_Choice == '2':
        z = 1.96
    elif CI_Choice == '3':
        z = 2.58
    elif CI_Choice == '4':
        z = zscore_toarea(roundingValue, None)

    CIlow = pointEstimate - (z * SE)
    CIhigh = pointEstimate + (z * SE)
    marginError = (CIhigh - CIlow) / 2

    results = [CIlow, CIhigh, marginError]
    output = [round(i, roundingValue) for i in results]

    print('There is a 95% chance that the true population [parameter] ' +
          'lies between the interval: (' + str(output[0]) + ',' +
          str(output[1]) + ') ± ' + str(output[2]))

    return results

def nonparametric_calcs(roundingValue):
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
        test_choice = str(input('|1: Wilcoxon Signed (dependent)|2: Wilcoxon Rank Sum (independent)| '))
        if test_choice == '1':
            statistic, p = st.wilcoxon(data[0], data[1])
        elif test_choice == '2':
            statistic, p = st.ranksums(data[0], data[1])

    results = [statistic, p]
    output = [round(i, roundingValue) for i in results]

    print('Statistic: ' + str(output[0]))
    print('P-value: ' + str(output[1]))

    return results

def calcselection(calcChoice, roundingValue):
    if calcChoice == 1:
        zscore_calcs(roundingValue)
    elif calcChoice == 2:
        bayes(roundingValue)
    elif calcChoice == 3:
        directadjustment(roundingValue)
    elif calcChoice == 4:
        twobytwo(roundingValue)
    elif calcChoice == 5:
        histogram_feat(roundingValue)
    elif calcChoice == 6:
        estimation(roundingValue)
    elif calcChoice == 7:
        binomial(roundingValue)
    elif calcChoice == 8:
        hypothesis_calcs(roundingValue)
    elif calcChoice == 9:
        confidence_interval(roundingValue)
    elif calcChoice == 10:
        nonparametric_calcs(roundingValue)
def chooser():
    print('Enter the calc you would like to use')
    print('1 - zscore')
    print('2 - bayes')
    print('3 - direct adjustment')
    print('4 - 2X2 Solver')
    print('5 - Histogram')
    print('6 - Estimations')
    print('7 - Binomial')
    print('8 - Hypothesis Testing')
    print('9 - Confidence Interval')
    print('10 - Nonparametric Tests')
    calcChoice = int(input())

    roundingValue = int(input('Round result to how many decimal places? '))

    if roundingValue:
        calcselection(calcChoice, roundingValue)
    else:
        calcselection(calcChoice, 4)

chooser()
