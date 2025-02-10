from FileManagmentFunctions import *
import numpy as np

# specify the input file name
damage_ressources = 'damage_ressources.csv'
repair_ressources = 'repair_ressources.csv'

# Specify the number of simulations of each energy
NumDamageRepeat = 1000

#Specify the number of repair simulations for each damage simulations
#This means number of repair sims = NumRepairRepeat*NumDamageRepeat
NumRepairRepeat = 0

# if set to True both NHEJ and HR will be simulated
# if set to False only NHEJ will be simulated 
includeHR = True

eps = np.array([1,2,3,4,5,30,70,100])/1000

#----- pre-simulations -------

#Initial_files("damage",damage_ressources, NumDamageRepeat, NumRepairRepeat)

#DamarisSddMaker("damage_results", "damage_results_damaris", "SDDOutput.txt")

#Initial_files("repair",repair_ressources, NumDamageRepeat, NumRepairRepeat,includeHR)

# ----- pos-simulations ------
#run_checker("damage",damage_ressources, NumDamageRepeat, NumRepairRepeat)
#run_checker("repair",repair_ressources, NumDamageRepeat, NumRepairRepeat)

run_aggregator(damage_ressources, NumDamageRepeat, eps, NumRepairRepeat)
#copy_james_file_structure()
