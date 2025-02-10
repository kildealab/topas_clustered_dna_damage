from RbeCalcFunctions import *


region_strings = ['inner', 'inter', 'outer']

for region in region_strings:

    rbe = np.zeros((12,18))
    rbe_err = np.zeros((12,18))

    for i in np.arange(1,13):
        settings = {
        # 'dataset': ['inner','inter','outer'],
        'dataset': ['outer','inter','inner'],
        # 'dataset': ['inter'],
        'dose_column': 12,
        'dose_scaling': True,
        'dose_target': 1, #Gy
        'input_path_prefix': '/Users/nicolas/Documents/NiceGit/topas_clustered_dna_damage/plotting/Processed_PreRep_Continuous/',
        'input_path_sufix': '_DamageAllRuns.csv',
        'num_yield_columns': 13, # combined
        # 'num_yield_columns': 14, # for direct or indirect only
        'num_jobs_n': 100,
        'num_jobs_x': 950,
        'yield_column': i, # complex DSB - combined action
        'energy_strings': ['n1eV','n100eV','n1keV','n10keV','n50keV','n100keV','n500keV','n700keV','n800keV','n900keV', \
        'n1MeV','n1-1MeV','n1-2MeV','n1-3MeV','n1-5MeV','n2MeV','n5MeV','n10MeV']
        }

        rbe[i-1], rbe_err[i-1] = get_yield_rbe_vs_energy(region, settings)

    np.savetxt("PreRep_RBE_"+region+".csv", rbe)
    np.savetxt("PreRep_RBE_ERR_"+region+".csv", rbe_err)
        
