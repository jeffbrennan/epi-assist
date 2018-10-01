import math
import scipy.stats as st
import numpy as np

def directadjustment (roundingValue): 
    
    print ('======== DIRECT ADJUSTMENT =========')
    bins = int(input('Enter number of groups: '))
    popFactor = int(input('Enter population factor eg 1000, 10000...'))

    print ('Enter location 1 rates separated by a space')
    loc1Rates = [float(x) for x in input().split()]

    print ('Enter location 2 rates separated by a space')
    loc2Rates = [float(x) for x in input().split()]

    print ('Enter Standard populations separated by a space')
    standardPop = [float(x) for x in input().split()]

    totalStandardPop = sum(standardPop)
    
    expectedPop1 = [(loc1Rates[i] * standardPop[i]) for i in range(bins)]
    totalExpected1 = sum(expectedPop1)
    adjustedRate1 = (totalExpected1 / totalStandardPop) * popFactor
    
    expectedpop2 = [(loc2Rates[i] * standardPop[i]) for i in range(bins)]
    totalexpected2 = sum(expectedpop2)
    adjustedrate2 = (totalexpected2 / totalStandardPop) * popFactor
   
    popFactorText = ' per ' + str(popFactor) + ' persons'

    print (resultsDivider)
    print('Location 1: ' + str(adjustedRate1) + popFactorText + ' | Location 2: ' + str(adjustedrate2) + popFactorText)

    return adjustedRate1, adjustedrate2


##Zscore calculations##
def zScoreCalcs(roundingValue):
    
    print ('================== Z SCORE FINDER ================== ')
    print ('1 - Z Score | 2 - Z Score -> Area | 3 - Area -> Z Score | ' +
            '4 - Obs given percentile | 5 - Confidence Interval')

    typechoice = str(input())
    zchooser(roundingValue, typechoice)
def zchooser(roundingValue, choice):
    
    if choice == '1':
        zscore_value(roundingValue)
    elif choice == '2':
        zscore_toarea(roundingValue, None)
    elif choice == '3':
        area_tozscore(roundingValue)
    elif choice == '4':
        zscore_observation(roundingValue)
    elif choice == '5':
        confidence_interval(roundingValue, None)
def zscore_value(roundingValue):
    
    Obs = float(input('Enter the observation: '))
    Mean = float(input('Enter the mean: '))
    SD = float(input('Enter SD: '))

    zresult = (Obs - Mean) / SD
    zresult = float(round(zresult, roundingValue))

    print (resultsDivider)
    print('Z Score: ' + str(zresult))

    zscore_toarea(roundingValue, zresult)

    return zresult

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

def area_tozscore(roundingValue):
    area = input('Enter your desired area (returns zscore): ')
    area = float(area)
    zresult = st.norm.ppf(area)
    zresult = round(zresult, roundingValue)
    print (zresult)
    
    return zresult

def zscore_observation(roundingValue): ## returns observation at a given percentile
    
    percentile = float(input('Enter percentile: '))
    SD = float(input('Enter SD: '))
    Mean = float(input('Enter mean: '))

    if percentile > 1:
        zresult = st.norm.ppf(percentile/100)
    else:
        zresult = st.norm.ppf(percentile)

    observation = (zresult * SD) + Mean
    observation = round(observation, roundingValue)

    print (str(percentile) + ' percentile observation: ' + str(observation))

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
    typechoice = str(input())

    if typechoice == '1':
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

            PVP = A / (A + B)
            PVN = D / (C + D)
           
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

            sensitivity = A / (A + C)
            specificity = B / (B + D)


        results = [PVP, PVN, sensitivity, specificity, A, B, C, D]
        output = [round(i, roundingValue) for i in results]
        dividerLength = len(str(output[4])) + len(str(output[5])) + 3

        print ('======================== RESULTS ========================')
        print ('_' * dividerLength)
        print ('|' + str(output[4]) + '|' + str(output[5]) + '|')
        print ('-' * dividerLength)
        print ('|' + str(output[6]) + '|' + str(output[7]) + '|')
        print ('_' * dividerLength + '\n')

        print ('Predictive positive value (PVP): ' + str(output[0]))
        print ('Predictive negative value (PVN): ' + str(output[1]))
        print ('Sensitivity: ' + str(output[2]))
        print ('Specificity: ' + str(output[3]))


    elif typechoice == '2':
        A = int(input('Enter value in cell A: '))
        B = int(input('Enter value in cell B: '))
        C = int(input('Enter value in cell C: '))
        D = int(input('Enter value in cell D: '))

        totalPop = A+B+C+D

        sensitivity = A / (A+C)
        specificity = D / (B+D)

        PVP = A / (A+B)
        PVN = D / (C+D)

        RR = (A/(A+B)) / (C/(C+D))
        OR = (A*D) / (B*C)

        eIncidence = (A/(A+B))
        nonEIncidence = (C/(C+D))
        popIncidence = (A+C) / totalPop

        PAR = (popIncidence - nonEIncidence) / popIncidence

        results = [PVP, PVN, sensitivity, specificity, OR, 
                    RR, eIncidence, nonEIncidence, popIncidence, PAR]

        output = [round(i, roundingValue) for i in results]

        print ('Predictive positive value (PVP): ' + str(output[0]))
        print ('Predictive negative value (PVN): ' + str(output[1]))
        print ('Sensitivity: ' + str(output[2]))
        print ('Specificity: ' + str(output[3]))
        print ('Odds Ratio (case control): ' + str(output[4]))
        print ('Risk Ratio (cohort): ' + str(output[5]))
        print ('Exposed Incidence: ' + str(output[6]))
        print ('Unexposed Incidence: ' + str(output[7]))
        print ('Total population incidence: ' + str(output[8]))
        print ('Population attributable risk: ' + str(output[9]))

        return results

def histogramfeat(roundingValue):
    print ('Enter histogram data set as a list')
    data = [float(x) for x in input().split()]
    
    kurtosisCalc = st.kurtosis(data)
    skewCalc = st.skew(data)
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

def estimation(roundingValue):
    print ('1 - Population parameters | 2 - Sample parameters')
    typechoice = str(input())

    if typechoice == '1':
        sampleN = int(input('Enter number in sample: '))
        
        print('Enter observations for all N in sample')
        samplevalues = [float(x) for x in input().split()]

        if len(samplevalues) == sampleN:
     
            mean = sum(samplevalues) / sampleN
            
            squareDif = [(i - mean)**2 for i in samplevalues]
            variance = sum(squareDif) / sampleN

            SD = math.sqrt(variance)
            SE = variance / (math.sqrt(sampleN))

            results = [mean, variance, SD, SE]
            output = [round(i, roundingValue) for i in results]
            print ('Population mean: ' + str(output[0]))
            print ('Population variance: ' + str(output[1]))
            print ('Population SD: ' + str(output[2]))

            return results

        else:
            print ('Entered incorrect number of observations - restarting...')
            estimation(roundingValue)

    elif typechoice == '2':
        
        sampleN = int(input('Enter number in sample: '))

        print('Enter observations for all N in sample')
        samplevalues = [float(x) for x in input().split()]
        
        if len(samplevalues) == sampleN:
            mean = sum(samplevalues) / sampleN
            
            squareDif = [(i - mean)**2 for i in samplevalues]
            variance = sum(squareDif) / (sampleN - 1)

            SD = math.sqrt(variance)
            
            SE = variance / (math.sqrt(sampleN))

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

def hypothesis(roundingValue):
    print ('================== HYPOTHESIS TESTER ================== ')

    def tailTest(tailChoice, tailSide):
        popSD =  sampleSD / math.sqrt(sampleSize)
    
        zValue = (sampleMean - nullHypo) / popSD
        zArea = zscore_toarea(roundingValue, zValue)
        
        if tailChoice == '1':
            pValue = round(1 - float(zArea), int(roundingValue))
        elif tailChoice == '2':
            pValue = round(2 * (1 - float(zArea)), int(roundingValue))
            print ('P-value: ' + str(pValue))

        print(hypoDecision(pValue, alphaValue, tailChoice, tailSide))

    def hypoDecision(testValue, alphaValue, tailChoice, tailSide):
        
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
    
    tailChoice = str(input('1: One-tailed test | 2: Two-tailed test: '))
    nullHypo = float(input('Enter null hypothesis value: '))
    sampleSize = float(input('Enter sample size: '))
    sampleMean = float(input('Enter sample mean: '))
    sampleSD = float(input('Enter sample SD: '))
    alphaValue = float(input('Enter desired significance: '))

    if tailChoice == '1':
        print ('Enter greater or less than null: 1 - greater (upper) | 2 - less than (lower)')
        tailSide = str(input())
        tailTest(tailChoice, tailSide)
    elif tailChoice == '2':
        tailTest(tailChoice, None)

def calcselection(choice, roundingValue):
    if choice == 1:
        zScoreCalcs(roundingValue)
    elif choice == 2:
        bayes(roundingValue)
    elif choice == 3: 
        directadjustment(roundingValue)
    elif choice == 4:
        twobytwo(roundingValue)
    elif choice == 5:
        histogramfeat(roundingValue)
    elif choice == 6:
        estimation(roundingValue)
    elif choice == 7:
        binomial(roundingValue)
    elif choice == 8:
        hypothesis(roundingValue)

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
    choice = int(input())  
    
    roundingValue = int(input('Round result to how many decimal places? '))

    if roundingValue:
        calcselection(choice, roundingValue)
    else:
        calcselection(choice, '4')
resultsDivider = ('======================== RESULTS ========================')
chooser()