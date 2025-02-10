
import numpy as np
import os 
import shutil
from ComplexDSbCounter import *


def RdChecker(file_path):
    try:
        with open(file_path, 'r') as file:
            # Count the number of lines in the file
            lines = 0
            for _ in file:
                lines += 1
                if lines >= 10:
                    return True
        return False
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return False




def read_after_header(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    start_processing = False
    processed_lines = []

    for line in lines:
        if start_processing:
            if line.strip():  # Check if the line is not empty
                processed_lines.append(line.strip())
        elif "***EndOfHeader***" in line:
            start_processing = True

    return processed_lines

eps = np.array([40,41,42,43,44,45,46,47,48,49])/1000
volume = "outer"


doses  = np.array([3,5,7,9,11,13,15,17,19])
for dose in doses:

    EndPoints = np.zeros(3+len(eps))

    NewFileIndex = 1

    HomeIndex = 1
    LeadingIndex = 1

    count = 0 

    denominator = 0

    NewFolder = "x250keV_"+volume+"_electron_"+str(dose)+"Gy"

    os.mkdir(NewFolder)
    shutil.copyfile("header.txt", NewFolder+"/"+"x250keV_"+volume+"_electron_"+str(dose)+"Gy_"+str(NewFileIndex)+"_SDDOutput.txt")

    while LeadingIndex < 1000:

        if RdChecker("x250keV_"+volume+"_electron/rd_x250keV_"+volume+"_electron_"+str(LeadingIndex)+".txt") == True:
            lines = read_after_header("x250keV_"+volume+"_electron/x250keV_"+volume+"_electron_"+str(LeadingIndex)+"_SDDOutput.txt")

            myfile =  open(NewFolder+"/"+"x250keV_"+volume+"_electron_"+str(dose)+"Gy_"+str(NewFileIndex)+"_SDDOutput.txt", "a")


            for line in lines:
                myfile.write("\n"+line)
            myfile.close()

            count += 1 

        LeadingIndex += 1 

        if count == dose:

            denominator += 1
            #print("count: "+str(count))
            #print(NewFolder+"/x250keV_"+volume+"_electron_"+str(NewFileIndex)+"_SDDOutput.txt"+"\n")
            EndPoints += clusterer(NewFolder+"/"+"x250keV_"+volume+"_electron_"+str(dose)+"Gy_"+str(NewFileIndex)+"_SDDOutput.txt", eps)

            count = 0 
            HomeIndex = LeadingIndex
            NewFileIndex += 1
            shutil.copyfile("header.txt", NewFolder+"/"+"x250keV_"+volume+"_electron_"+str(dose)+"Gy_"+str(NewFileIndex)+"_SDDOutput.txt")

    EndPoints /= denominator

    f = open(volume+"_EndPoints.txt", "a")
    f.write(str(dose)+","+",".join(f"{point:.3f}" for point in EndPoints[:3+len(eps)]) + "\n")
    f.close()






