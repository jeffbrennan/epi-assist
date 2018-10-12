import random

## complete random assigner
def randomizer(patientIDs, p):

    treatmentGroup = []
    controlGroup = []

    while len(controlGroup) < p and len(treatmentGroup) < p:
        for patient in patientIDs:

            checker = random.uniform(0,1)
            if checker > 0.5:
                if len(controlGroup) < p:
                    controlGroup.append(patient)
                else:
                    treatmentGroup.append(patient)
            elif checker < 0.5:
                if len(treatmentGroup) < p:
                    treatmentGroup.append(patient)
                else:
                    controlGroup.append(patient)
    
    return treatmentGroup, controlGroup

# randomized block
def chunks(patientIDs, chunkNumber):
    output = []
    for i in range(0, len(patientIDs), chunkNumber):
         output.append(patientIDs[i:i + chunkNumber])
    return output

def chooser(inputList, choice):
    p = len(inputList) // 2
    
    if choice == '1':
        treatmentIDs, controlIDs = randomizer(inputList, p)
        print('Treatment IDs: ' + str(treatmentIDs))
        print('Control IDs: ' + str(controlIDs))

    elif choice =='2':
        chunkNumber = int(input('Enter number of blocks needed: '))
        random.shuffle(inputList)
        print(chunks(inputList, chunkNumber))

    # elif choice =='3':

print ('Enter subjects (IDs, names etc) separated by a space')
inputList = [x for x in input().split()]

try:
    inputList = [int(x) for x in inputList]
except ValueError:
    inputList = [str(x) for x in inputList]

print ('Number of subjects entered: ' + str(len(inputList)))

print ('Which method of randomization would you like to use? 1 - Complete |'+
        '2 - Blocks | 3 - Hierarchical')
choice = str(input())

chooser(inputList, choice)

## hiearchical 


