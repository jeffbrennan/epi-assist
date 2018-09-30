import scipy.stats as st
import math
import numpy as np

def directadjustment (decimalchoice): 
    
    loc1Rates = []
    loc2Rates = []
    standardpop = []
    expectedpop1 = []
    expectedpop2 = []

    print ('======== DIRECT ADJUSTMENT =========')
    print ('Enter number of groups')
    bins = int(input())

    print ('Enter population factor eg 1000, 10000...')
    popfactor = int(input())

    print ('Enter location 1 rates in order')
    for _ in range(bins):
        loc1Rates.append(float(input()) / popfactor)

    print ('Enter location 2 rates in order')
    for _ in range(bins):
        loc2Rates.append(float(input()) / popfactor)

    print ('Enter Standard populations in order')
    for _ in range(bins):
        standardpop.append(int(input()))

    for i in range(bins):
        expectedpop1.append(loc1Rates[i] * standardpop[i])
        expectedpop2.append(loc2Rates[i] * standardpop[i])

    totalexpected1 = sum(expectedpop1)
    totalexpected2 = sum(expectedpop2) 
    totalstandardpop = sum(standardpop)

    adjustedrate1 = (totalexpected1 / totalstandardpop) * popfactor
    adjustedrate2 = (totalexpected2 / totalstandardpop) * popfactor

    adjustedrate1 = formatter(adjustedrate1, decimalchoice)
    adjustedrate2 = formatter(adjustedrate2, decimalchoice)

    popfactortext = ' per ' + str(popfactor) + ' persons'

    outputter(('Location 1: ' + str(adjustedrate1) + popfactortext + ' | Location 2: ' + str(adjustedrate2) + popfactortext))
    
def zScoreCalcs(decimalchoice):
    def zscore_value():
        print ('Enter the observation')
        Obs = float(input())
        print ('Enter the mean')
        Mean = float(input())
        print ('Enter SD')
        SD = float(input())

        zresult = (Obs - Mean) / SD

        zresult = float(formatter(zresult, decimalchoice))
        outputter('Z Score: ' + str(zresult))

        zscore_area(zresult)
    
    def zscore_area(calc_zscore): 
        if calc_zscore:
            zscore = calc_zscore    
        else:
            print ('Enter your desired Z score')
            zscore = input()

            print ('Enter your desired area')
            area = input()

        if zscore:
            zscore = round(float(zscore), 2) ## biostats class uses table with 2 decimal place z scores

            areacalc = st.norm.cdf(zscore)
            inversearea = 1 - areacalc
            areapercent = areacalc * 100
            areacalc = formatter(areacalc, decimalchoice)
            inversearea = formatter(inversearea, decimalchoice)
            areapercent = formatter(areapercent, decimalchoice)    
            outputter(('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' +
                        '|Remaining area: ' + str(inversearea)))

            return areacalc

        elif area:
            area = float(area)
            zresult = st.norm.ppf(area)
            zresult = formatter(zresult, decimalchoice)
            print (zresult)
            
            return zresult

    def zscore_observation():
       
        print ('Enter percentile')
        percentile = float(input())
        print ('Enter SD')
        SD = float(input())
        print ('Enter mean')
        Mean = float(input())

        zresult = st.norm.ppf(percentile)

        observation = (zresult * SD) + Mean
        observation = formatter(observation, decimalchoice)

        print ('Observation at the desired percentile: ' + observation)

    def confidence_interval(calc_SE):
        print ('Enter the point estimate')
        pointEstimate = float(input())
        
        if calc_SE:
            SE = calc_SE 
        else:
            
            print ('Enter the standard deviation')
            SD = float(input())
            print ('Enter the sample size')
            sampleSize = float(input())

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
            z = zscore_area(None)

        CIlow = pointEstimate - (z * SE)
        CIhigh = pointEstimate + (z * SE)
        marginError = (CIhigh - CIlow) / 2

        results = [CIlow, CIhigh, marginError]
        output = [formatter(i, decimalchoice) for i in results]

        print ('There is a 95% chance that the true population [parameter]' +  
        'lies between the interval of ' + output[0] + ' and' + output[1])
    
    def zchooser(choice):
        if choice == '1':
            zscore_value()
        elif choice == '2':
            zscore_area(None)
        elif choice == '3':
            zscore_observation()
        elif choice == '4':
            confidence_interval(None)
   
    print ('================== Z SCORE FINDER ================== ')

    print ('1 - Z Score | 2 - Area given Z Score | 3 - Obs given percentile | 4 - Confidence Interval')

    typechoice = str(input())
    zchooser(typechoice)

def zscore_area(decimalchoice, calc_zscore): 
    if calc_zscore:
        zscore = calc_zscore    
    else:
        print ('Enter your desired Z score')
        zscore = input()

        print ('Enter your desired area')
        area = input()

    if zscore:
        zscore = round(float(zscore), 2) ## biostats class uses table with 2 decimal place z scores

        areacalc = st.norm.cdf(zscore)
        inversearea = 1 - areacalc
        areapercent = areacalc * 100
        areacalc = formatter(areacalc, decimalchoice)
        inversearea = formatter(inversearea, decimalchoice)
        areapercent = formatter(areapercent, decimalchoice)    
        outputter(('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' +
                    '|Remaining area: ' + str(inversearea)))

        return areacalc

    elif area:
        area = float(area)
        zresult = st.norm.ppf(area)
        zresult = formatter(zresult, decimalchoice)
        print (zresult)
        
        return zresult

def bayes(decimalchoice):
    print ('====== BAYES CALC ========')

    print ('Enter disease prevalence (as a proportion)')
    prevalence = float(input())
    print ('Enter accuracy for disease test (as a proportion)')
    diseaseAccuracy = float(input())
    print ('Enter accuracy for complement test (as a proportion)')
    complementAccuracy = float(input())

    numerator = (prevalence * diseaseAccuracy)
    denominator = (numerator + ((1-prevalence) * complementAccuracy))

    result = numerator / denominator
    result = formatter(result, decimalchoice)
    
    outputter(result)

def twobytwo(decimalchoice):

    print('====== 2x2 TABLE SOLVER ======')

    print ('1 - table empty | 2 - table full')
    typechoice = str(input())

    if typechoice == '1':

        print ('Enter sensitivity')
        sensitivity = input()
        print ('Enter specificity')
        specificity = input()

        print ('Enter PVP')
        PVP = input()
        print ('Enter PVN')
        PVN = input()

        if sensitivity and specificity:

            sensitivity = float(sensitivity)
            specificity = float(specificity)

            print ('Enter prevalence of population')
            popprev = float(input())
            print ('Enter size of population')
            popsize = float(input())

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
            
            print ('Enter size of population')
            popsize = int(input())

            print ('Enter number of positive test results')
            totalposresults = int(input())
               
            AB = totalposresults
            A = totalposresults * PVP
            B = AB - A

            CD = popsize - totalposresults
            D = CD * PVN
            C = CD - D

            sensitivity = A / (A + C)
            specificity = B / (B + D)

        results = [PVP, PVN, sensitivity, specificity, A, B, C, D]
        output = [formatter(i, decimalchoice) for i in results]

        print ('======================== RESULTS ========================')
        print (str(output[4]) + '|' + str(output[5]))
        print (str(output[6]) + '|' + str(output[7]) + '\n')
        print ('Predictive positive value (PVP): ' + str(output[0]))
        print ('Predictive negative value (PVN): ' + str(output[1]))
        print ('Sensitivity: ' + str(output[2]))
        print ('Specificity: ' + str(output[3]))

    elif typechoice == '2':
        print ('Enter value in cell A')
        A = int(input())
        print ('Enter value in cell B')
        B = int(input())
        print ('Enter value in cell C')
        C = int(input())
        print ('Enter value in cell D')
        D = int(input())

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

        output = [formatter(i, decimalchoice) for i in results]

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

def histogramfeat(decimalchoice):
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

    if kurtosisCalc > 3:
        print ('Data is leptokurtic')
    elif kurtosisCalc < 3:
        print ('Data is platykurtic ')

    skewCalc = formatter(skewCalc, decimalchoice)
    kurtosisCalc = formatter(kurtosisCalc, decimalchoice)

    print ('Skewness: ' + str(skewCalc))
    print ('Kurtosis: ' + str(kurtosisCalc))

def estimation(decimalchoice):
    print ('1 - Population parameters | 2 - Sample parameters | 3 - CI')
    typechoice = str(input())

    if typechoice == '1':
        print('Enter number in sample')
        sampleN = int(input())

        print('Enter observations for all N in sample')
        samplevalues = [float(x) for x in input().split()]

        if len(samplevalues) == sampleN:
     
            mean = sum(samplevalues) / sampleN
            
            squareDif = [(i - mean)**2 for i in samplevalues]
            variance = sum(squareDif) / sampleN

            SD = math.sqrt(variance)
            SE = variance / (math.sqrt(sampleN))
            results = [mean, variance, SD, SE]
            output = [formatter(i, decimalchoice) for i in results]

            print ('Population mean: ' + str(output[0]))
            print ('Population variance: ' + str(output[1]))
            print ('Population SD: ' + str(output[2]))
        else:
            print ('Entered incorrect number of observations - restarting...')
            estimation(decimalchoice)

    elif typechoice == '2':
        
        print('Enter number in sample')
        sampleN = int(input())

        print('Enter observations for all N in sample')
        samplevalues = [float(x) for x in input().split()]
        
        if len(samplevalues) == sampleN:
            mean = sum(samplevalues) / sampleN
            
            squareDif = [(i - mean)**2 for i in samplevalues]
            variance = sum(squareDif) / (sampleN - 1)

            SD = math.sqrt(variance)
            
            SE = variance / (math.sqrt(sampleN))

            results = [mean, variance, SD, SE]
            output = [formatter(i, decimalchoice) for i in results]

            print ('Sample mean: ' + str(output[0]))
            print ('Sample variance: ' + str(output[1]))
            print ('Sample SD: ' + str(output[2]))
            print ('Sample SE: ' + str(output[3]))

        else:
            print ('Entered incorrect number of observations - restarting...')
            estimation(decimalchoice)  

    elif typechoice =='3':
        
        print('Enter desired CI as a probability')
        CI = float(input())

        print ('Enter sample SD')
        sampleSD = input()

        print ('Enter sample size')
        sampleSize = input()

        print ('Enter sample point estimate (mean etc.)')
        sampleEstimate = input()

        if sampleSD and sampleSize and sampleEstimate:
            sampleSD = float(sampleSD)
            sampleSize = float(sampleSize)
            sampleEstimate = float(sampleEstimate)

            zPercentile = CI + ((1 - CI) / 2) ## outputs needed percentile to find upper z-score
            zresult = st.norm.ppf(zPercentile) ## gets zscore at calculated percentile
        
            SE = sampleSD / math.sqrt(sampleSize)
        
            ## CI = point estimate +/- (z-score * SE)
            lowCI = sampleEstimate - (zresult * SE)
            highCI = sampleEstimate + (zresult * SE)

            results = [-zresult, zresult, SE, lowCI, highCI]
            output = [formatter(i, decimalchoice) for i in results]

            print ('========================== RESULT ==========================')
            print ('Low Z-score: ' + str(output[0]))
            print ('High Z-score: ' + str(output[1]))
            print ('SE: ' + str(output[2]))
            print ('CI: (' + str(output[3]) + ', ' + str(output[4]) + ')')
            print ('We are 95% confident that the true population parameter is between ' + str(output[3]) + ' and ' + str(output[4]))

        else:
            
            zPercentile = CI + ((1 - CI) / 2) ## outputs needed percentile to find upper z-score
            zresult = st.norm.ppf(zPercentile) ## gets zscore at calculated percentile   results = [-zresult, zresult, SE, lowCI, highCI]
            
            results = [-zresult, zresult]
            output = [formatter(i, decimalchoice) for i in results]

            print ('========================== RESULT ==========================')
            print ('Low Z-score: ' + str(output[0]))
            print ('High Z-score: ' + str(output[1]))
            
def binomial(decimalchoice):
    
    print ('Enter n')
    n = int(input())
    
    print ('Enter x')
    x = int(input())

    print ('Enter p')
    p = float(input())
    
    choose = math.factorial(n) / (math.factorial(x) * math.factorial(n-x))
    result = (choose * (p ** x) * (1 - p)  ** (n - x))
    resultpercent = result * 100

    result = formatter(result, decimalchoice)
    resultpercent = formatter(resultpercent, str(int(decimalchoice) - 2))

    print ('Probability: ' + result + ' | ' + resultpercent + '%')

def hypothesis(decimalchoice):
    print ('================== HYPOTHESIS TESTER ================== ')

    def tailTest(tailChoice, tailSide):
        popSD =  sampleSD / math.sqrt(sampleSize)
    
        zValue = (sampleMean - nullHypo) / popSD
        zArea = zscore_area(decimalchoice, zValue)
        
        if tailChoice == '1':
            pValue = round(1 - float(zArea), int(decimalchoice))
        elif tailChoice == '2':
            pValue = round(2 * (1 - float(zArea)), int(decimalchoice))
            print (pValue)

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

def outputter(result):
    print ('=========================== RESULT ===========================')
    print (result)

def formatter(rawNum, decimalchoice):
        decimalformat = ('{:.' + decimalchoice + 'f}')
        rawNum = decimalformat.format(rawNum)
        return rawNum

def calcselection(choice, decimalchoice):
    if choice == 1:
        zScoreCalcs(decimalchoice)
    elif choice == 2:
        bayes(decimalchoice)
    elif choice == 3: 
        directadjustment(decimalchoice)
    elif choice == 4:
        twobytwo(decimalchoice)
    elif choice == 5:
        histogramfeat(decimalchoice)
    elif choice == 6:
        estimation(decimalchoice)
    elif choice == 7:
        binomial(decimalchoice)
    elif choice == 8:
        hypothesis(decimalchoice)

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
    
    print ('Enter desired number of decimal places')
    decimalchoice = str(input())

    if decimalchoice:
        calcselection(choice, decimalchoice)
    else:
        calcselection(choice, '4')

chooser()