precision = float(1e-15)
printErrors = True

filePath1 = "/Users/zamanlantra/phd/repo/nbody/oppic/simpic_oppic_soa/files/"
filePath2 = "/Users/zamanlantra/phd/repo/nbody/OP-PIC/app_simpic_oppic/files/"

fileNames = [
    "F_0_E.dat",
    "F_0_density.dat",
    "F_0_phi.dat",
    "F_0_vxx.dat"]


import csv

simulationPass = True

def isWithinPrecision(a1, a2):
    if abs(a1-a2) > precision:
        return False  
    return True

def getFileArray(fileName):
    fileArray1 = []
    with open(fileName, 'r') as file1:
        csvreader1 = csv.reader(file1)
        for row in csvreader1:
            fileArray1.append(row) 
    return fileArray1

print(f'Path1: {filePath1}')
print(f'Path2: {filePath2}\n')

for fileName in fileNames:
    print(f'CHECKING FILE {fileName}')
    file1 = getFileArray(filePath1 + fileName)
    file2 = getFileArray(filePath2 + fileName)

    filePass = True
    fileErrors = 0

    for i in range(len(file1)):

        line1 = file1[i]
        line2 = file2[i]

        lineError = False

        for j in range(1, len(line1)):

            if not isWithinPrecision(float(line1[j]), float(line2[j])):
                simulationPass = False
                filePass = False
                lineError = True
                fileErrors+=1
                
        if printErrors and lineError:
            print(f'ERROR... fileName={fileName} line1={line1} line2={line2}')

    if filePass:
        print("\tALL VALUES OK IN FILE, PASS")
    else:
        print(f'\tFAILED, CHECK ERRORS {fileErrors}')

if simulationPass:
    print("\nALL VALUES OK, PASS...")
else:
    print("\nFAILED, CHECK ERRORS")
