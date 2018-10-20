import math
import scipy.stats as st
import numpy as np
import pandas as pd

def directadjustment (roundingValue): 
    
    print ('======== DIRECT ADJUSTMENT =========')
    bins = int(input('Enter number of groups: '))
    popFactor = int(input('Enter population factor eg 1000, 10000...'))

    print ('Enter location 1 rates separated by a space:')
    loc1Rates = [float(x) for x in input().split()]

    print ('Enter location 2 rates separated by a space:')
    loc2Rates = [float(x) for x in input().split()]

    print ('Enter Standard populations separated by a space:') ## might help to have explanation on this point
    standardPop = [float(x) for x in input().split()]

    totalStandardPop = sum(standardPop)
    
    ##multiples the condition rates by the population in the standard group
    ##saves as list, then sums and creates rate by comparing expected to total
    expectedPop1 = [(loc1Rates[i] * standardPop[i]) for i in range(bins)] 
    totalExpected1 = sum(expectedPop1)
    adjustedRate1 = (totalExpected1 / totalStandardPop) * popFactor
    
    expectedpop2 = [(loc2Rates[i] * standardPop[i]) for i in range(bins)]
    totalexpected2 = sum(expectedpop2)
    adjustedrate2 = (totalexpected2 / totalStandardPop) * popFactor
   
    popFactorText = ' per ' + str(popFactor) + ' persons'

    print (resultsDivider)
    print('Location 1: ' + str(adjustedRate1) + popFactorText + 
        ' | Location 2: ' + str(adjustedrate2) + popFactorText)

    return adjustedRate1, adjustedrate2

### Zscore calculations - 1 ###
def zscore_calcs(roundingValue):
    
    print ('================== Z SCORE FINDER ================== ')
    print ('1 - Z Score | 2 - Z Score -> Area | 3 - Area -> Z Score | ' +
            '4 - Obs given percentile | 5 - Confidence Interval')

    typeChoice = str(input())
    zscore_chooser(roundingValue, typeChoice)

def zscore_chooser(roundingValue, choice):
    
    if choice == '1':
        zscore_value(roundingValue)
    elif choice == '2':
        zscore_toarea(roundingValue, None)
    elif choice == '3':
        area_tozscore(roundingValue, None)
    elif choice == '4':
        zscore_observation(roundingValue)
    elif choice == '5':
        confidence_interval(roundingValue, None)

def zscore_value(roundingValue):
    
    Obs = float(input('Enter the observation: '))
    Mean = float(input('Enter the mean: '))
    SD = float(input('Enter SD: '))

    zResult = (Obs - Mean) / SD
    zResult = float(round(zResult, roundingValue))

    print (resultsDivider)
    print('Z Score: ' + str(zResult))

    zscore_toarea(roundingValue, zResult)

    return zResult

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
    
    print(resultsDivider)
    print('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' +
                '|Remaining area: ' + str(inversearea))

    return areacalc

def area_tozscore(roundingValue, percentile):
    area = input('Enter your desired area, blank if auto (returns zscore): ')
    if area:
        area = float(area)
    else: 
        area = percentile
    
    if area > 1:
        zResult = st.norm.ppf(percentile/100)
    else:
        zResult = st.norm.ppf(percentile)

    print (zResult)
    return zResult

def zscore_observation(roundingValue): ## returns observation at a given percentile
    
    percentile = float(input('Enter percentile: '))
    SD = float(input('Enter SD: '))
    Mean = float(input('Enter mean: '))

    zResult = area_tozscore(roundingValue, percentile)

    observation = (zResult * SD) + Mean
    observation = round(observation, roundingValue)

    print (str(percentile) + ' percentile observation: ' + str(observation))
    return observation

def confidence_interval(roundingValue, calc_SE):
    
    pointEstimate = float(input('Enter the point estimate: '))
    
    if calc_SE:
        SE = calc_SE
    else:
        SD = float(input('Enter the standard deviation: '))
        sampleSize = float(input('Enter the sample size: '))

        SE = SD / math.sqrt(sampleSize)

    print ('Enter desired % confidence - 1: 90 | 2: 95 | 3: 99 | 4: Custom')
    choice = str(input())

    if choice == '1':
        z = 1.645
    elif choice == '2':
        z = 1.96
    elif choice =='3':
        z = 2.58
    elif choice == '4':
        z = zscore_toarea(roundingValue, None)

    CIlow = pointEstimate - (z * SE)
    CIhigh = pointEstimate + (z * SE)
    marginError = (CIhigh - CIlow) / 2

    results = [CIlow, CIhigh, marginError]
    output = [round(i, roundingValue) for i in results]

    print ('There is a 95% chance that the true population [parameter] ' +  
    'lies between the interval: (' + str(output[0]) + ',' + str(output[1]) +
     ') ± ' + str(output[2]))

    return results

### Bayes - 2 ###      
def bayes(roundingValue):
    print ('====== BAYES CALC ========')
    print('Enter values as a proportion')

    prevalence = float(input('Enter disease prevalence: '))
    diseaseAccuracy = float(input('Enter accuracy for disease test: '))
    complementAccuracy = float(input('Enter accuracy for complement test: '))

    numerator = (prevalence * diseaseAccuracy)
    denominator = (numerator + ((1-prevalence) * complementAccuracy))

    result = numerator / denominator
    result = round(result, roundingValue)

    print (resultsDivider)
    print (result)

    return result

def twobytwo(roundingValue):

    print('====== 2x2 TABLE SOLVER ======')

    print ('1 - table empty | 2 - table full')
    typeChoice = str(input())

    if typeChoice == '1':
        print ('Enter sensitivity and specificity or PVP and PVN')
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
    
    ## calls helper functions to calculate basic values
    sensitivity, specificity = sens_spec(A, B, C, D)
    PVP, PVN = PVP_PVN(A, B, C, D)
    RR, OR = RR_OR(A, B, C, D)
    RR, OR = RR_OR(A, B, C, D)

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

    col1 = [A,C]
    col2 = [B,D]

    tableData = {'+':col1,'-':col2}
    twoTable = pd.DataFrame(tableData, index=['Exposed', 'Nonexposed'])

    print ('======================== RESULTS ========================')

    print (twoTable)

        ## consider converting to dictionary
    print ('\n' + 'Predictive positive value (PVP): ' + str(output[0]))
    print ('Predictive negative value (PVN): ' + str(output[1]))
    print ('Sensitivity: ' + str(output[2]))
    print ('Specificity: ' + str(output[3]))
    print ('Odds Ratio (case control): ' + str(output[4]))
    print ('Error adjusted OR: ' + str(errorOR))
    print ('Risk Ratio (cohort): ' + str(output[5]))
    print ('Error adjusted RR: ' + str(errorRR))
    print ('Exposed Incidence: ' + str(output[6]))
    print ('Unexposed Incidence: ' + str(output[7]))
    print ('Total population incidence: ' + str(output[8]))
    print ('Attributable risk: ' + str(output[9]))
    print ('Percent attributable risk: ' + str(output[10]))

    return results

def RR_OR(A, B, C, D):

    RR = (A/(A+B)) / (C/(C+D))
    OR = (A*D) / (B*C)

    return RR, OR

def sens_spec(A, B, C, D):
    sensitivity = A / (A+C)
    specificity = D / (B+D)

    return sensitivity, specificity

def PVP_PVN(A, B, C, D):
    PVP = A / (A+B)
    PVN = D / (C+D)

    return PVP, PVN

def incidence_2x2(A, B, C, D):
    totalPop = A+B+C+D
   
    eIncidence = (A/(A+B))
    nonEIncidence = (C/(C+D))
    popIncidence = (A+C) / totalPop
    AR = eIncidence - nonEIncidence
    PAR = ((popIncidence - nonEIncidence) / popIncidence) * 100
    return eIncidence, nonEIncidence, popIncidence, AR, PAR

def error_2x2(A, B, C, D):

    errorType = str(input('Enter error type (1 - nondifferential |2 - differential): '))
    errorRate = float(input('Enter misclassification rate'))
    errorDirection = str(input('Enter error direction (1 - E->NE |2 - NE->E): '))
   
    if errorType == '1' and errorDirection == '1':
        A -= (A*errorRate)
        B -= (B*errorRate)
        C += (A*errorRate)
        D += (B*errorRate)

    if errorType == '1' and errorDirection == '2':
        
        A += (C*errorRate)
        B += (D*errorRate)
        C -= (C*errorRate)
        D -= (D*errorRate)

    if errorType == '2':
        groupSelector = str(input('Enter differential group (1 - cases|2- controls): '))

        if groupSelector == '1' and errorDirection == '1':
            A -= (A*errorRate)
            C += (A*errorRate)

        elif groupSelector == '1' and errorDirection == '2':
            A+= (C*errorRate)
            C-= (C*errorRate)

        elif groupSelector == '2' and errorDirection == '1':
            B -= (B*errorRate)
            D += (B*errorRate)

        elif groupSelector == '2' and errorDirection == '2':
            B += (D*errorRate)
            D -= (D*errorRate)
            
    errorRR, errorOR = RR_OR(A, B, C, D)
    return errorRR, errorOR

def histogram_feat(roundingValue):
    print ('Enter histogram data set as a list')
    data = [float(x) for x in input().split()]
    
    kurtosisCalc = st.kurtosis(data)
    skewCalc = st.skew(data)
    shapiroCalc = st.shapiro(data)
    mean = np.mean(data)
    median = np.median(data)

    print ('=============RESULT===========')
   
    if mean > median:
        print('Histogram is right-skewed')
    elif mean < median:
        print('Histogram is left-skewed')
    elif mean == median:
        print('Histogram is symmetric')
   
    ## !! Need to check decisions
    if kurtosisCalc > 3:
        print ('Data is leptokurtic')
    elif kurtosisCalc < 3:
        print ('Data is platykurtic ')

    skewCalc = round(skewCalc, roundingValue)
    kurtosisCalc = round(kurtosisCalc, roundingValue)

    print ('Skewness: ' + str(skewCalc))
    print ('Kurtosis: ' + str(kurtosisCalc))
    print ('Shapiro-Wilk Test: ' + str(shapiroCalc))

    return skewCalc, kurtosisCalc, shapiroCalc

def estimation(roundingValue):
    print ('1 - Population parameters | 2 - Sample parameters')
    typeChoice = str(input())

    if typeChoice == '1':
        sampleSize = int(input('Enter number in sample: '))
        
        print('Enter observations for all N in sample separated by a space')
        sampleValues = [float(x) for x in input().split()]

        if len(sampleValues) == sampleSize:
     
            mean, variance, SD, SE = listCalculations(sampleValues, sampleSize)

            results = [mean, variance, SD, SE]
            output = [round(i, roundingValue) for i in results]
            print ('Population mean: ' + str(output[0]))
            print ('Population variance: ' + str(output[1]))
            print ('Population SD: ' + str(output[2]))

            return results

        else:
            print ('Entered incorrect number of observations - restarting...')
            estimation(roundingValue)

    elif typeChoice == '2':
        
        sampleSize = int(input('Enter number in sample: '))

        print('Enter observations for all N in sample')
        sampleValues = [float(x) for x in input().split()]
        
        if len(sampleValues) == sampleSize:

            mean, variance, SD, SE = listCalculations(sampleValues, sampleSize)

            results = [mean, variance, SD, SE]
            output = [round(i, roundingValue) for i in results]

            print ('Sample mean: ' + str(output[0]))
            print ('Sample variance: ' + str(output[1]))
            print ('Sample SD: ' + str(output[2]))
            print ('Sample SE: ' + str(output[3]))

            return results

        else:
            print ('Entered incorrect number of observations - restarting...')
            estimation(roundingValue)  

def listCalculations (numberList, sampleSize):

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
    
    choose = math.factorial(n) / (math.factorial(x) * math.factorial(n-x))
    result = (choose * (p ** x) * (1 - p)  ** (n - x))
    resultpercent = result * 100

    result = round(result, roundingValue)
    resultpercent = round(resultpercent, str(int(roundingValue) - 2))

    print ('Probability: ' + result + ' | ' + resultpercent + '%')

    return result

def hypothesis_calcs(roundingValue):
    print ('================== HYPOTHESIS TESTER ================== ')
    
    testChoice = str(input('1: Z-test | 2: T-Test: '))
    tailChoice = str(input('1: One-tailed test | 2: Two-tailed test: '))
    
    if tailChoice == '1':
        print ('Enter greater or less than null: 1 - greater (upper) | 2 - less than (lower)')
        tailSide = str(input())
        
        if testChoice == '1':        
            hypothesis_zTest(roundingValue, tailChoice, tailSide)
        elif testChoice == '2':
            hypothesis_tTest(roundingValue, tailChoice, tailSide)

    elif tailChoice == '2':
        
        if testChoice == '1':
            hypothesis_zTest(roundingValue, tailChoice, None)
        elif testChoice == '2':
            hypothesis_tTest(roundingValue, tailChoice, None)

def hypothesis_zTest(roundingValue, tailChoice, tailSide):
    
    nullValue = float(input('Enter null hypothesis value: '))
    sampleSize = float(input('Enter sample size: '))
    sampleMean = float(input('Enter sample mean: '))
    sampleSD = float(input('Enter sample SD: '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

    popSD =  sampleSD / math.sqrt(sampleSize)

    zValue = (sampleMean - nullValue) / popSD

    zArea = zscore_toarea(roundingValue, zValue)
    
    if tailChoice == '1':
        pValue = round(1 - float(zArea), int(roundingValue))
    elif tailChoice == '2':
        pValue = round(2 * (1 - float(zArea)), int(roundingValue))
        print ('P-value: ' + str(pValue))

    print(hypothesis_decision(pValue, alphaValue, tailChoice, tailSide))

def hypothesis_tTest(roundingValue, tailChoice, tailSide):

    tType = str(input('1: One Sample | 2: 2-Sample | 3: Indepdendent 2-Sample: '))
    alphaValue = float(input('Enter desired significance (0.10|0.05|0.01): '))

    if tType == '1' or tType == '2':

        nullValue = float(input('Enter null hypothesis value: '))
        sampleSize = float(input('Enter sample size: '))
        sampleMean = float(input('Enter sample mean: '))
        sampleSD = float(input('Enter sample SD: '))

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

        df = hypothesis_degreesFreedom(s1, n1, s2, n2)
        criticalT = math.fabs(st.t.ppf(alphaValue, df))
        
        print (tScore)
        print (criticalT)

        if tScore < criticalT:
            print ('FTR Null')
        elif tScore > criticalT:
            print ('Reject Null')
            
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
    
def hypothesis_degreesFreedom(s1, n1, s2, n2):

    ## independent two sample test
    if s1 and s2 and n1 and n2:
        numerator = ((s1**2 / n1) + (s2**2 / n2))**2
        denominator = ((s1**2/n1)**2 / (n1-1) + (s2**2/n2)**2 / (n2-1))

        df = int(math.floor(numerator / denominator))

    print ('Degrees of Freedom: ' + str(df))
    return df

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

def chooser():
    print ('Enter the calc you would like to use')
    print ('1 - zscore')
    print ('2 - bayes')
    print ('3 - direct adjustment')
    print ('4 - 2X2 Solver')
    print ('5 - Histogram')
    print ('6 - Estimations')
    print ('7 - Binomial')
    print ('8 - Hypothesis Testing')
    calcChoice = int(input())  
    
    roundingValue = int(input('Round result to how many decimal places? '))

    if roundingValue:
        calcselection(calcChoice, roundingValue)
    else:
        calcselection(calcChoice, 4)

resultsDivider = ('======================== RESULTS ========================')
chooser()