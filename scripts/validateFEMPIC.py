precision = float(1e-15)
printErrors = False

filePath1 = "/ext-home/zl/phd/OP-PIC/fempic_new/files/"
filePath2 = "/ext-home/zl/phd/OP-PIC/fempic_new/files/"

# F_1_c_cell_electric_field_acomp.dat  F_1_s_cell_electric_field_acomp.dat          1e-6 this should be the problem
# F_1_c_cell_electric_field.dat        F_1_s_cell_electric_field.dat                1e-25
# F_1_c_cell_shape_deriv.dat           F_1_s_cell_shape_deriv.dat                   1e-25
# F_1_c_node_charge_density_acomp.dat  F_1_s_node_charge_density_acomp.dat          1e+1 this should not be a problem
# F_1_c_node_charge_density_am.dat     F_1_s_node_charge_density_am.dat             1e-12
# F_1_c_node_potential_acomp.dat       F_1_s_node_potential_acomp.dat               1e-25 wow :)
# F_1_c_part_mesh_relation_am.dat      F_1_s_part_mesh_relation_am.dat              1e-25
# F_1_c_part_mesh_relation_bm.dat      F_1_s_part_mesh_relation_bm.dat              1e-25
# F_1_c_part_position_am.dat           F_1_s_part_position_am.dat                   1e-18
# F_1_c_part_position_bi.dat           F_1_s_part_position_bi.dat                   1e-25
# F_1_c_part_position_bm.dat           F_1_s_part_position_bm.dat                   1e-19          
# F_1_c_part_velocity_am.dat           F_1_s_part_velocity_am.dat                   1e-25
# F_1_c_part_velocity_bm.dat           F_1_s_part_velocity_bm.dat                   1e-25
# F_1_c_part_lc_am.dat                 F_1_s_part_lc_am.dat                         1e-14
# F_1_c_cell_volume_am.dat             F_1_s_cell_volume_am.dat                     1e-25    
# F_1_c_cell_determinants_am.dat       F_1_s_cell_determinants_am.dat               1e-25          

fileNames1 = [
    "F_23_c_cell_electric_field_bi.dat",
    "F_23_c_part_mesh_relation_bm.dat",  
    "F_23_c_part_position_bm.dat",
    "F_23_c_part_velocity_bm.dat",
    "F_20_c_node_charge_density_acomp.dat",
    "F_20_c_node_charge_density_am.dat",
    "F_20_c_node_potential_acomp.dat",
    "F_20_c_part_lc_am.dat",
    ]

fileNames2 = [
    "F_23_s_cell_electric_field_bi.dat",
    "F_23_s_part_mesh_relation_bm.dat",  
    "F_23_s_part_position_bm.dat",
    "F_23_s_part_velocity_bm.dat",
    "F_20_s_node_charge_density_acomp.dat",
    "F_20_s_node_charge_density_am.dat",
    "F_20_s_node_potential_acomp.dat",
    "F_20_s_part_lc_am.dat",
    ]

import csv

simulationPass = True

def isWithinPrecision(a1, a2):
    if abs(a1-a2) > precision:
        return False  
    return True

def getFileArray(fileName):
    fileArray1 = []
    with open(fileName, 'r') as file1:
        csvreader1 = csv.reader(file1, delimiter=' ')
        for row in csvreader1:
            fileArray1.append(row) 
    return fileArray1

print(f'Precision - {precision}')

for i in range(len(fileNames1)):

    fileName1 = fileNames1[i]
    fileName2 = fileNames2[i]

    print(f'FILE PATH 1 - {filePath1} -- {fileName1}')
    print(f'FILE PATH 2 - {filePath2} -- {fileName2}')

    file1 = getFileArray(filePath1 + fileName1)
    file2 = getFileArray(filePath2 + fileName2)

    filePass = True
    fileErrors = 0
    lineCounter = 0
    maxDiff = 0.0
    lineErrorCount = 0

    for i in range(len(file1)):
        lineCounter += 1
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
                maxDiff = max(maxDiff, diff)
        if lineError:
            lineErrorCount += 1

        if printErrors and lineError:
            print(f'ERROR... {fileName1}|{lineCounter} line1={line1} line2={line2} | max_diff={diff}' )

    if filePass:
        print("\t********* ALL VALUES OK IN FILE, PASS *********\n")
    else:
        print(f'\t********* FAILED, CHECK ERRORS Line Erros: {lineErrorCount} Total Errors:{fileErrors} | Max Difference: {maxDiff} *********\n')

if simulationPass:
    print("\nALL VALUES OK, PASS...")
else:
    print("\nFAILED, CHECK ERRORS")
