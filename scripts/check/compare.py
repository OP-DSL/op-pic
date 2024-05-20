printLog = True
errorCount = 0

def Compare(fileName1, fileName2):
    # Open File in Read Mode
    file_1 = open(fileName1, 'r')
    file_2 = open(fileName2, 'r')
    
    print("Comparing files ", " @ " + fileName1, " # " + fileName2, sep='\n')
    
    file_1_line = file_1.readline()
    file_2_line = file_2.readline()
    
    # Use as a Counter
    line_no = 1
    errorCount = 0

    print("Difference Lines in Both Files")
    while file_1_line != '' or file_2_line != '':
    
        # Removing whitespaces
        file_1_line = file_1_line.rstrip()
        file_2_line = file_2_line.rstrip()
    
        # Compare the lines from both file
        if file_1_line != file_2_line:
            
            if (printLog and errorCount < 80): #  
                # otherwise output the line on file1 and use @ sign
                if file_1_line == '':
                    print("@", "Line-%d" % line_no, file_1_line)
                else:
                    print("@-", "Line-%d" % line_no, file_1_line)
                    
                # otherwise output the line on file2 and use # sign
                if file_2_line == '':
                    print("#", "Line-%d" % line_no, file_2_line)
                else:
                    print("#+", "Line-%d" % line_no, file_2_line)
        
                # Print a empty line
                print()

            errorCount += 1

        # Read the next line from the file
        file_1_line = file_1.readline()
        file_2_line = file_2.readline()
    
        line_no += 1
    
    print(f'Total Error Count: {errorCount}')

    file_1.close()
    file_2.close()

# -----------------------------------------------------------------------------------

# fileName1 =  '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/a/struct_com96_r0_Edge_Bef_M'
# fileName2 = '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/b/struct_com96_r0_Edge_Bef_M'

fileName1 =  '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/struct_com96_r0_Edge_Aft_R'
fileName2 = '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/struct_com96_r48_Edge_Aft_R'


Compare(fileName1, fileName2)

# print()

# fileName1 = '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/struct_com96_r0_Edge_Aft_R'
# fileName2 = '/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/files/struct_com48_r0_Edge_Aft_R'

# Compare(fileName1, fileName2)
