precision = float(1e-15)
printErrors = True

filePath1 = "/ext-home/zl/phd/OP-PIC/app_fempic_oppic/files_genseq/"
filePath2 = "/ext-home/zl/phd/OP-PIC/app_fempic_oppic/files_cuda/"

# filePath2 = "/ext-home/zl/phd/OP-PIC/backup_fempic/files/"

# fileNames = [
#     "F_50_node_charge_density.dat",
#     "F_50_node_potential.dat",
#     "F_50_cell_electric_field.dat",
#     "F_50_part_position.dat",
#     "F_50_part_velocity.dat",
#     "F_50_part_cell_index.dat"]

# fileNames = [
#     "F_1_node_charge_density.dat",
#     "F_1_node_potential.dat",
#     "F_1_cell_electric_field.dat",
#     "K_1_cell_determinants.dat",
#     "K_1_part_weights.dat",
#     "F_1_part_position.dat",
#     "F_1_part_velocity.dat",
#     "F_1_part_cell_index.dat",
#     "F_1_part_weights.dat",
#     "K_1_cell_volume.dat",
#     "K_1_node_volume.dat",
#     "B_1_part_weights.dat",
#     "B_1_node_charge_density.dat",
#     "B_1_node_volume.dat",
#     "A_1_part_weights.dat",
#     "A_1_node_charge_density.dat",
#     "A_1_node_volume.dat"]

# fileNames = [
#     "B_1_part_weights.dat",
#     "B_1_node_charge_density.dat",
#     "B_1_node_volume.dat",
#     "A_1_part_weights.dat",
#     "A_1_node_charge_density.dat",
#     "A_1_node_volume.dat"]

fileNames = [
    "B_1_cell_determinants.dat",
    "B_1_part_position.dat",
    "B_1_cell_volume.dat",
    "B_1_part_weights.dat",
    "B_1_node_charge_density.dat",
    "B_1_node_volume.dat",
    "A_1_node_charge_density.dat"]

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

print(f'Precision - {precision}')

for fileName in fileNames:
    print(f'FILE PATH 1 - {filePath1}')
    print(f'FILE PATH 2 - {filePath2}')
    print(f'CHECKING FILE NAME {fileName}')

    file1 = getFileArray(filePath1 + fileName)
    file2 = getFileArray(filePath2 + fileName)

    filePass = True
    fileErrors = 0

    for i in range(len(file1)):

        line1 = file1[i]
        line2 = file2[i]

        lineError = False
        diff = 0.0

        for j in range(1, len(line1)):
            
            if not isWithinPrecision(float(line1[j]), float(line2[j])):
                simulationPass = False
                filePass = False
                lineError = True
                fileErrors+=1
                diff = max(diff, abs(float(line1[j])-float(line2[j])))
                
        if printErrors and lineError:
            print(f'ERROR... fileName={fileName} line1={line1} line2={line2} | max_diff={diff}' )

    if filePass:
        print("\tALL VALUES OK IN FILE, PASS")
    else:
        print(f'\tFAILED, CHECK ERRORS {fileErrors}')

if simulationPass:
    print("\nALL VALUES OK, PASS...")
else:
    print("\nFAILED, CHECK ERRORS")
