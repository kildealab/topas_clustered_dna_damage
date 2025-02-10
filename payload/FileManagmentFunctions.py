import csv
import os
import shutil
import random
import re 
import glob
from ComplexDSbCounter import clusterer
#----------------- Defining some useful functions --------------------------

def ReInitiateFolder(folder_name):
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)

def check_for_results_files(directory):
    if os.path.exists(directory):
        for root, dirs, files in os.walk(directory):
            if len(files) > 0:
                print("Cannot initalize if there is results files")
                print("Exiting")
                exit()
       



def write_damage_variables(filename, SimNumber, JobName, nThreads):
    with open(filename, 'a') as damage_param_file:
        r1 = random.randint(10000,99999)
        damage_param_file.write('\n')
        damage_param_file.write('includeFile = supportFiles/spectra/spectrum_' + JobName + '.txt' + '\n')
        damage_param_file.write('includeFile = supportFiles/relative_doses/reldose_' + JobName + '.txt' + '\n')
        damage_param_file.write('s:Sc/ClusterScorer/FileAllEvents = "damage_results/' + JobName +'/' + JobName + "_" + str(SimNumber) + '_AllEvents.txt"' + '\n')
        damage_param_file.write('s:Sc/ClusterScorer/SDDOutputFile = "damage_results/' + JobName +'/' + JobName + "_" + str(SimNumber) + '_SDDOutput.txt"' + '\n')
        damage_param_file.write('s:Sc/ClusterScorer/RealDoseFile =  "damage_results/' + JobName +'/' + 'rd_' + JobName + "_" + str(SimNumber)  + '.txt"' +'\n')
        damage_param_file.write('i:Ts/Seed = '+ str(r1) + '\n')
        damage_param_file.write('i:Ts/NumberOfThreads = ' + str(nThreads))

# write the repair parameter files
def write_repair_variables(filename, SimNumber, JobName, NumRepairRepeat, includeHR):
    with open(filename, 'a') as repair_param_file:
    # write the second and third columns to the output file
        r2 = random.randint(10000,99999)
        repair_param_file.write('\n')
        repair_param_file.write('i:Ts/Seed = '+ str(r2) + '\n')
        repair_param_file.write('s:Ch/DaMaRiS/STDFormatDamageFileName = "/home/ndesja7/scratch/payload/damage_results/' + JobName +'/' + JobName + "_" + str(SimNumber) + '_SDDOutput.txt"' + '\n')
        repair_param_file.write('i:Ch/DaMaRiS/BiologyRepeatNumber = ' + str(NumRepairRepeat) + '\n')
        repair_param_file.write('i:Ch/DaMaRiS/AlternativeRunID = ' + str(SimNumber) + '\n')
        if includeHR == True:
            repair_param_file.write('includefile = /home/ndesja7/scratch/payload/supportFiles/damaris/pathwayHR.txt')
        else:
             repair_param_file.write('includefile = /home/ndesja7/scratch/payload/supportFiles/damaris/pathwayNHEJ.txt')


def write_sbatch_file(sbatch_filename, time, memory,JobName, nThreads, JobArray, SimType):
    with open(sbatch_filename, 'w') as sbatch_file:
            sbatch_file.write('#!/bin/bash' + '\n' + '#' '\n')
            sbatch_file.write('#SBATCH --time=' + time +'\n')
            sbatch_file.write('#SBATCH --mem='+ memory + 'G' + '\n')
            sbatch_file.write('#SBATCH --job-name='+ JobName + '\n')
            sbatch_file.write('#SBATCH --nodes=1' + '\n')
            sbatch_file.write('#SBATCH -n ' + nThreads + ' # number of cores' + '\n')
            sbatch_file.write('#SBATCH --array=' + JobArray + '\n')

            if SimType == "damage":
                sbatch_file.write('#SBATCH --output=' + SimType  +'_results/' + JobName + '/' + JobName + '-%a.out' + '\n \n')
            else:
                sbatch_file.write('#SBATCH --output=' + JobName + '-%a.out' + '\n \n')

            sbatch_file.write('echo "Starting run at: `date`"' + '\n')
            sbatch_file.write('~/topas/bin/topas ~/scratch/payload/' + SimType + '_params/' + JobName + '/' + JobName + '$SLURM_ARRAY_TASK_ID.txt' + '\n')
            sbatch_file.write('echo "Job finished with exit code $? at: `date`"' + '\n')

def write_sbatch_list(filename,JobName,SimType):
     with open(filename, 'a') as f:
        if SimType == "damage":
            f.write("sbatch " + SimType + "_sbatch/" + JobName + ".sh;\n")
        else:
            f.write("cd ~/scratch/payload/repair_results/" + JobName +";\n")
            f.write("sbatch ~/scratch/payload/" + SimType + "_sbatch/" + JobName + ".sh;\n")

def check_for_xray_damage(SimType, JobName):
    if SimType == "damage" and "x250keV" in JobName:
        return True
    else:
        return False


def Initial_files(SimType,ressource_table, NumDamageRepeat, NumRepairRepeat, includeHR = True):


    # Intiating some directories
    # The subdirectories of 'damage_results' and 'repair_results' will be fill up by the simulations
    folder_to_create = [SimType +'_params', SimType + '_sbatch', SimType + '_results']

    check_for_results_files(folder_to_create[2])

    for folder_name in folder_to_create:
        ReInitiateFolder(folder_name)

    #Clearing job list
    job_list = SimType + "_JobList"
    if os.path.exists(job_list):
        os.remove(job_list)

    with open(ressource_table, 'r') as file:
        reader = csv.reader(file)
        next(reader) # skip the header row
        for row in reader:
            JobName, nThreads, memory, time = row
            
            # Create new subdirectories for each particle 
            os.mkdir(SimType + '_params/' + JobName)
            os.mkdir(SimType + '_results/' + JobName)

            # Each of the simulation repeat get its own top parameter file 
            for SimNumber in range(1,NumDamageRepeat+1):
                filename = SimType + '_params/' + JobName + '/' + JobName + str(SimNumber) + '.txt'

                #Copies the parameters that needs to be in the top parameter file and which are common to all simulations for that SimType
                if check_for_xray_damage(SimType, JobName):
                    shutil.copy('damage_xray_static_top_params.txt', filename) 
                else:
                    shutil.copy(SimType + '_static_top_params.txt', filename) 
                #Write the simulations specific parameters 
                if SimType == "damage":
                    write_damage_variables(filename, SimNumber, JobName, nThreads)
                elif SimType == "repair":
                    write_repair_variables(filename, SimNumber, JobName, NumRepairRepeat, includeHR)
                else:
                    Print("ERROR: Invalid SimType")

            #Writing the sbatch file for a the particle specie/energy   
            sbatch_filename = SimType + '_sbatch/' + JobName + '.sh'
            JobArray = '1-' + str(NumDamageRepeat)
            write_sbatch_file(sbatch_filename, time, memory,JobName, nThreads, JobArray, SimType)

            #Writting the job list
            write_sbatch_list(job_list, JobName, SimType)


def count_checker(filename,expected, condition= "equal", HasHeader = False):

    if os.path.exists(filename):
        with open(filename, 'r') as file:
            count = sum(1 for line in file if line.strip())
            
        if HasHeader == False:
            if condition == "equal" and count == int(expected):
                return True
            if condition == "more" and count > int(expected):
                return True
            if condition == "less" and count < int(expected):
                return True
            else:
                return False

        if HasHeader == True:
            if condition == "equal" and count == int(expected)+1:
                return True
            if condition == "more" and count > int(expected)+1:
                return True
            if condition == "less" and count < int(expected)+1:
                return True
            else:
                return False
    else:
        return False

def add_RealDose(filename):
    total = 0 
    if os.path.exists(filename):
        with open(filename, "r") as f:
            for line in f:
                total += float(line.strip())
    return total



def run_checker(SimType,ressource_table, NumDamageRepeat, NumRepairRepeat):


    ReInitiateFolder('new_' + SimType + '_sbatch')

    with open(ressource_table, 'r') as file:
        reader = csv.reader(file)
        next(reader) # skip the header row
        for row in reader:

            JobName, nThreads, memory, time = row

            failed_sims = []
            for SimNumber in range(1,NumDamageRepeat+1):


                if SimType == "damage":
                    filename = 'damage_results/' + JobName + '/' + "rd_" + JobName + "_" + str(SimNumber) + '.txt'
                    success = count_checker(filename, nThreads)

                elif SimType == "repair":
                    filename = 'repair_results/' + JobName + '/' + "PerRepResults" + str(SimNumber) + '.out'
                    success = count_checker(filename, NumRepairRepeat, HasHeader = True)
                    
                else:
                    print("ERROR: Invalid SimType")

                if not success:
                    failed_sims.append(SimNumber)

            if failed_sims:       
                sbatch_filename = 'new_' + SimType + '_sbatch/' + JobName + '.sh'
                JobArray = str(failed_sims)[1:-1]
                write_sbatch_file(sbatch_filename, time, memory, JobName, nThreads, JobArray, SimType)


def run_aggregator(ressource_table, NumDamageRepeat, eps, NumRepairRepeat = 0):


    if NumRepairRepeat == 0: 

        ReInitiateFolder('final_damage_results')

        with open(ressource_table, 'r') as file:
            reader = csv.reader(file)
            next(reader) # skip the header row
            for row in reader:

                JobName, nThreads, memory, time = row

                for SimNumber in range(1,NumDamageRepeat+1):

                    rd = 'damage_results/' + JobName + '/' + "rd_" + JobName + "_" +str(SimNumber) + '.txt'
                    rd_success = count_checker(rd, nThreads)

                    dsb = 'damage_results/' + JobName + '/' + JobName + "_" + str(SimNumber) + '_SDDOutput.txt'
                    dsb_success = count_checker(dsb, 26, condition= "more")

                    FinalFile = "final_damage_results/" + JobName +"_DamageAllRuns.csv"

                    RealDose = add_RealDose(rd)

                    if dsb_success:
                        
                        EndPoints =  clusterer(dsb,eps)

                        with open(FinalFile, "a") as f2:

                            OuputLine = str(SimNumber)
                            for endpoint in EndPoints:
                                OuputLine += "," + str(int(endpoint))
                            OuputLine += "," +str(RealDose) + "\n"

                            f2.write(OuputLine)

                    if rd_success and not dsb_success:
                        with open(FinalFile, "a") as f2:

                            OuputLine = str(SimNumber)
                            for endpoint in EndPoints:
                                OuputLine += ",0"
                            OuputLine += ","+str(RealDose) + "\n"

                            f2.write(OuputLine)



    if NumRepairRepeat > 0: 

        ReInitiateFolder('final_results')

        with open(ressource_table, 'r') as file:
            reader = csv.reader(file)
            next(reader) # skip the header row
            for row in reader:

                JobName, nThreads, memory, time = row

                for SimNumber in range(1,NumDamageRepeat+1):

                    rd = 'damage_results/' + JobName + '/' + "rd_" + JobName + "_" +str(SimNumber) + '.txt'
                    rd_success = count_checker(rd, nThreads)

                    dsb = 'damage_results/' + JobName + '/' + JobName + "_" + str(SimNumber) + '_SDDOutput.txt'
                    dsb_success = count_checker(dsb, 26, condition= "more")

                    repair = 'repair_results/' + JobName + '/' + "PerRepResults" + str(SimNumber) + '.out'
                    repair_success = count_checker(repair, NumRepairRepeat, HasHeader = True)


                    FinalFile = "final_results/" + JobName +"_RepAllRuns.csv"
                    
                    if repair_success:

                        RealDose = add_RealDose(rd)
                        with open(repair, "r") as f1:
                            lines = f1.readlines()

                        with open(FinalFile, "a") as f2:
                            comma_separated_results = ",".join(lines[1].rstrip().split())
                            f2.write(str(SimNumber) + "," + comma_separated_results + "," + str(RealDose) + "\n")

                    if rd_success and not dsb_success:
                        with open(FinalFile, "a") as f2:
                            f2.write(str(SimNumber) + ",0,0,0,0,0,0,0," + str(RealDose) + "\n")

def DamarisEnabler(file_path):
    
    with open(file_path, "r") as f:
        content = f.readlines()
        
    with open(file_path, "w") as f:
        
        Header = True 
        
        for line in content:
            
            if Header == False:
                if re.compile(";").split(line)[5][7] != "0":
                    f.write(line)
                    
            else:
                if line == "Data entries, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0;\n":
                    f.write("Data entries, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;\n")
                else:
                    f.write(line)
                    
                if line == "***EndOfHeader***;\n":
                    Header = False


def DamarisSddMaker(src_dir, dest_dir, file_ending):
    
    if os.path.exists(dest_dir):
        print(f"Destination directory {dest_dir} already exists. Please remove it or choose another location.")
    else:
        # Copy the entire directory
        shutil.copytree(src_dir, dest_dir)

    # Find all files in the directory that match the pattern
    file_list = glob.glob(os.path.join(dest_dir, "**/*" + file_ending))

    # Pass each file to the function
    for file in file_list:
        DamarisEnabler(file)




def copy_james_file_structure():

    scoring_stings = ['inner','inter','outer']

    energy_strings = [
    'n1eV','n100eV','n1keV','n10keV','n50keV','n100keV','n500keV','n700keV','n800keV','n900keV', \
    'n1MeV','n1-1MeV','n1-2MeV','n1-3MeV','n1-5MeV','n2MeV','n5MeV','n10MeV','x250keV'
    ]

    particles_strings = ['proton', 'alpha', 'electron']

    os.mkdir("james_file_structure") 

    for scoring in scoring_stings:
        for energy in energy_strings:
            for particle in particles_strings:
                os.makedirs(os.path.join("james_file_structure", scoring, energy, particle))

    folder_path = "final_damage_results"  # Replace with the path to your folder

    # Get all files in folder
    files = os.listdir(folder_path)

    # Print name of each file
    for file in files:
        energy, scoring, particle  = file.split("_")[0:3]
        os.rename(os.path.join(folder_path,file),os.path.join("james_file_structure", scoring, energy, particle,file))


                    
        

           




