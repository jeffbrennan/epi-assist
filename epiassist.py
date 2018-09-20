import scipy.stats as st


def directadjustment (): 
    
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

    print ('Enter desired number of decimal places')
    decimalchoice = str(input())

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
    
def zscore():
    print ('======= Z SCORE FINDER ==========')

    print ('1 - Z Score | 2 - Area given Z Score | 3 - Obs given percentile')

    typechoice = str(input())

    if typechoice == '1':
        print ('Enter the observation')
        Obs = float(input())
        print ('Enter the mean')
        Mean = float(input())
        print ('Enter SD')
        SD = float(input())

        print ('Enter desired number of decimal places')
        decimalchoice = str(input())

        zresult = (Obs - Mean) / SD

        zresult = formatter(zresult, decimalchoice)

        outputter('Z Score: ' + str(zresult))

    elif typechoice == '2':

        print ('Enter your desired Z score')
        zscore = float(input())
        zscore = float(formatter(zscore, '2')) ## biostats class uses table with 2 decimal place z scores

        print ('Enter desired number of decimal places')
        decimalchoice = str(input())

        areacalc = st.norm.cdf(zscore)
        areapercent = areacalc * 100

        areacalc = formatter(areacalc, decimalchoice)
        areapercent = formatter(areapercent, decimalchoice)    

        outputter(('Area: ' + str(areacalc) + '|' + str(areapercent) + '%' ))

    elif typechoice == '3':
       
        print ('Enter percentile')
        percentile = float(input())
        print ('Enter SD')
        SD = float(input())
        print ('Enter mean')
        Mean = float(input())

        print ('Enter desired number of decimal places')
        decimalchoice = str(input())

        zresult = st.norm.ppf(percentile)

        observation = (zresult * SD) + Mean
        observation = formatter(observation, decimalchoice)

        print ('Observation at the desired percentile: ' + observation)

def bayes():
    print ('====== BAYES CALC ========')

    print ('Enter disease prevalence (as a proportion)')
    prevalence = float(input())
    print ('Enter accuracy for disease test (as a proportion)')
    diseaseAccuracy = float(input())
    print ('Enter accuracy for complement test (as a proportion)')
    complementAccuracy = float(input())

    print ('Enter desired number of decimal places')
    decimalchoice = str(input())

    numerator = (prevalence * diseaseAccuracy)
    denominator = (numerator + ((1-prevalence) * complementAccuracy))

    result = numerator / denominator
    result = formatter(result, decimalchoice)
    
    outputter(result)

def twobytwo():

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

        print ('Enter desired number of decimal places')
        decimalchoice = str(input())

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

        print ('Enter desired number of decimal places')
        decimalchoice = str(input())
    
        sensitivity = A / (A+C)
        specificity = D / (B+D)
        PVP = A / (A+B)
        PVN = D / (C+D)

        output = [formatter(i, decimalchoice) for i in [PVP, PVN, sensitivity, specificity]]

        print ('Predictive positive value (PVP): ' + str(output[0]))
        print ('Predictive negative value (PVN): ' + str(output[1]))
        print ('Sensitivity: ' + str(output[2]))
        print ('Specificity: ' + str(output[3]))

def skew():
    print ('Enter histogram data set as a list')
    data = [float(x) for x in input().split()]
    
    print ('Enter desired number of decimal places')
    decimalchoice = str(input())


    skewcalc = st.skew(data)
    print ('=============RESULT===========')
    if skewcalc > 1: 
        print ('Your data is right skewed')
    elif skewcalc < 1:
        print ('Your data is left skewed')
    elif skewcalc == 0: 
        print ('Your data is normal')

    skewcalc = formatter(skewcalc, decimalchoice)

    print ('Skewness: ' + str(skewcalc))

def outputter(result):
    print ('=========================== RESULT ===========================')
    print (result)

def formatter(rawNum, decimalchoice):
        decimalformat = ('{:.' + decimalchoice + 'f}')
        rawNum = decimalformat.format(rawNum)
        return rawNum

def calcselection(choice):
    if choice == 1:
        zscore()
    elif choice == 2:
        bayes()
    elif choice == 3: 
        directadjustment()
    elif choice == 4:
        twobytwo()
    elif choice == 5:
        skew()

def chooser():
    print ('Enter the calc you would like to use')
    print ('1 - zscore')
    print ('2 - bayes')
    print ('3 - direct adjustment')
    print ('4 - 2X2 Solver')
    print ('5 - Histogram skew')

    choice = int(input())  
    
    calcselection(choice)

chooser()