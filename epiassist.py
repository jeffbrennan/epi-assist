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
    
def zscore(decimalchoice):
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
            zscore = float(input())
            zscore = float(formatter(zscore, '2')) ## biostats class uses table with 2 decimal place z scores

        areacalc = st.norm.cdf(zscore)
        areapercent = areacalc * 100

        areacalc = formatter(areacalc, decimalchoice)
        areapercent = formatter(areapercent, decimalchoice)    

        outputter(('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' ))

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

    def zchooser(choice):
        if choice == '1':
            zscore_value()
        elif choice == '2':
            zscore_area(None)
        elif choice == '3':
            zscore_observation()
   
    print ('================== Z SCORE FINDER ================== ')

    print ('1 - Z Score | 2 - Area given Z Score | 3 - Obs given percentile')

    typechoice = str(input())
    zchooser(typechoice)

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
           
            ## insert pandas df with 2x2 table

            output = [formatter(i, decimalchoice) for i in [PVP, PVN, sensitivity, specificity]]

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

            output = [formatter(i, decimalchoice) for i in [PVP, PVN, sensitivity, specificity]]

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
    
        sensitivity = A / (A+C)
        specificity = D / (B+D)
        PVP = A / (A+B)
        PVN = D / (C+D)

        output = [formatter(i, decimalchoice) for i in [PVP, PVN, sensitivity, specificity]]

        print ('Predictive positive value (PVP): ' + str(output[0]))
        print ('Predictive negative value (PVN): ' + str(output[1]))
        print ('Sensitivity: ' + str(output[2]))
        print ('Specificity: ' + str(output[3]))

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

    print ('Probability: ' + result + ' | ' + (resultpercent) + '%')

def outputter(result):
    print ('=========================== RESULT ===========================')
    print (result)

def formatter(rawNum, decimalchoice):
        decimalformat = ('{:.' + decimalchoice + 'f}')
        rawNum = decimalformat.format(rawNum)
        return rawNum

def calcselection(choice, decimalchoice):
    if choice == 1:
        zscore(decimalchoice)
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

def chooser():
    print ('Enter the calc you would like to use')
    print ('1 - zscore')
    print ('2 - bayes')
    print ('3 - direct adjustment')
    print ('4 - 2X2 Solver')
    print ('5 - Histogram')
    print ('6 - Estimations')
    print ('7 - Binomial')
    choice = int(input())  
    
    print ('Enter desired number of decimal places')
    decimalchoice = str(input())

    if decimalchoice:
        calcselection(choice, decimalchoice)
    else:
        calcselection(choice, '4')

chooser()