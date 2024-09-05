from FileManagmentFunctions import *

# specify the input file name
damage_ressources = 'damage_ressources.csv'
repair_ressources = 'repair_ressources.csv'

# Specify the number of simulations of each energy
NumDamageRepeat = 125:

#Specify the number of repair simulations for each damage simulations
#This means number of repair sims = NumRepairRepeat*NumDamageRepeat
NumRepairRepeat = 0

# if set to True both NHEJ and HR will be simulated
# if set to False only NHEJ will be simulated 
includeHR = True

#----- pre-simulations -------

#Initial_files("damage",damage_ressources, NumDamageRepeat, NumRepairRepeat)
#Initial_files("repair",repair_ressources, NumDamageRepeat, NumRepairRepeat,includeHR)

# ----- posy 
#run_checker("damage",damage_ressources, NumDamageRepeat, NumRepairRepeat)
#run_checker("repair",repair_ressources, NumDamageRepeat, NumRepairRepeat)

run_aggregator(damage_ressources, NumDamageRepeat, NumRepairRepeat)
#copy_james_file_structure()
