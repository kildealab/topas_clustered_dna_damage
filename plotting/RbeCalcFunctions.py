#!/usr/bin/env python
# coding: utf-8

# In[137]:


import numpy as np
import os, copy
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import AutoMinorLocator
import matplotlib
#matplotlib.rcParams['text.usetex'] = True



#matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
from scipy.interpolate import interp1d
from dateutil import parser
# import plot_helpers
from scipy import stats
from cycler import cycler


relative_dose_map = {
    'proton': 0,
    'electron': 9,
    'alpha': 3,
}

energy_map = {
    'n1eV': 0,
    'n100eV': 1,
    'n1keV': 2,
    'n10keV': 3,
    'n50keV': 4,
    'n100keV': 5,
    'n500keV': 6,
    'n700keV': 7,
    'n800keV': 8,
    'n900keV': 9,
    'n1MeV': 10,
    'n1-1MeV': 11,
    'n1-2MeV': 12,
    'n1-3MeV': 13,
    'n1-5MeV': 14,
    'n2MeV': 15,
    'n5MeV': 16,
    'n10MeV': 17,
    'x250keV': 18,
}


def import_yields_species(path_prefix,path_sufix,dataset,energy,species,num_jobs,num_yield_columns):
    '''Import array of all types of DNA damage yields. Each row contains the yields for a distinct
    simulation job. The first 3 columns non-yield info. The remaining 5 columns contain yields for
    5 types of DNA damage.'''
    
    file_yields = path_prefix + energy +"_"+ dataset +"_"+ species + path_sufix

    if os.path.isfile(file_yields):
        # yields = np.genfromtxt(file_yields, delimiter=',')
        # yields = yields[0:num_jobs,:]
#         yields = np.array(pd.read_csv(file_yields, header=None).sample(num_jobs))
#         print('\t',pd.read_csv(file_yields, header=None).values.shape)
        yields = np.array(pd.read_csv(file_yields, header=None).iloc[:num_jobs])
    else:
        yields = np.zeros((num_jobs,num_yield_columns))

    return yields

def import_delivered_dose_species(path_prefix,path_sufix, dataset,energy,species,num_jobs,dose_column):
    '''Import vector of delivered dose values across multiple simulation jobs for the specified
    scoring dataset, incident particle energy, and secondary particle species.'''
    file_doses = path_prefix + energy +"_"+ dataset +"_"+ species + path_sufix
    
    if os.path.isfile(file_doses):
        doses = np.genfromtxt(file_doses, delimiter=',')
        # doses = doses[:,dose_column]
        doses = doses[0:num_jobs,dose_column]
    else:
        doses = np.zeros((num_jobs,1))

    return doses

def import_relative_doses(path_prefix,dataset):
    '''Import relative dose values for all incident particle energies and secondary particle species
    for a particular scoring dataset. Result is a 2D array wherein the rows span the secondary
    particle species and the columns span the incident particle energies.'''
    file_doses = path_prefix + 'relative_doses/relative_dose_' + dataset + '.csv'
    relative_doses = np.genfromtxt(file_doses, delimiter=',')

    # Normalize the relative dose matrix to account for only species of interest
    doses_normalized = relative_doses / np.sum(relative_doses[[relative_dose_map['proton'], \
        relative_dose_map['electron'],relative_dose_map['alpha']],:],axis=0)

    return doses_normalized


def scale_species_yields_by_dose(yields,dataset,energy,species,target_dose,path_prefix,path_sufix,num_jobs,dose_column):
    '''Scale (reduce) DNA damage yields for a secondary particle species using the amount of dose
    delivered relative to the target dose scaled by relative dose contribution for that species.'''
    
    # Vector of deliverd dose values
    delivered_doses = import_delivered_dose_species(path_prefix,path_sufix,dataset,energy,species,num_jobs,dose_column)
    
    # Matrix of relative dose values for all species, normalized to species of interest
    relative_doses = import_relative_doses(path_prefix,dataset)
    
    # Relative dose value of interest for current incident particle energy & secondary particle species
    relative_dose = relative_doses[relative_dose_map[species],energy_map[energy]]

    # Calculate dose scaling factor, i.e. the ratio of delivered dose to the target dose for this
    # secondary particle species. Handle case where the target dose is zero (& thus delivered doses
    # are zero b/c no simulations with that particle were run). In this case, the yields will also
    # be zero, so can set scaling factors to 1

    if relative_dose > 0.001:
        dose_scaling = delivered_doses / (target_dose*relative_dose)
    else:
        dose_scaling = np.ones(len(delivered_doses))

    # Apply scaling factor and return
    #print(dataset, ", ", energy, ", ", species)
    yields = yields / dose_scaling

    return yields


def import_yields(path_prefix,path_sufix,dataset,energy,scale_by_dose,target_dose,num_jobs,num_yield_columns,yield_column,dose_column):
    '''Import DNA damage yields for all species of interest, scale by dose delivered relative to
    target dose on a species-by-species if desired, and combine results for all species.'''
    # print('*******************************************************************')
    # print(energy)
    # print(yield_column)
    if (yield_column >= 0):
        yields_p = import_yields_species(path_prefix,path_sufix,dataset,energy,'proton',num_jobs,num_yield_columns)
        yields_e = import_yields_species(path_prefix,path_sufix,dataset,energy,'electron',num_jobs,num_yield_columns)
        yields_a = import_yields_species(path_prefix,path_sufix,dataset,energy,'alpha',num_jobs,num_yield_columns)

        # Keep only yield column of interest (e.g. SSB or DSB)
        yields_p = yields_p[:,yield_column]
        yields_e = yields_e[:,yield_column]        
        yields_a = yields_a[:,yield_column]

        # print(yields_p)
    elif (yield_column == -1):
        clusters_p = import_clusters_species(path_prefix,dataset,energy,'comp_dsb','proton')
        clusters_e = import_clusters_species(path_prefix,dataset,energy,'comp_dsb','electron')
        clusters_a = import_clusters_species(path_prefix,dataset,energy,'comp_dsb','alpha')

        if (clusters_p.size != 0):
            clusters_p = split_clusters_by_job(clusters_p)
        if (clusters_e.size != 0):
            clusters_e = split_clusters_by_job(clusters_e)
        if (clusters_a.size != 0):
            clusters_a = split_clusters_by_job(clusters_a)

        yields_p = np.zeros((num_jobs,1))[:,0]
        yields_e = np.zeros((num_jobs,1))[:,0]
        yields_a = np.zeros((num_jobs,1))[:,0]

        # print(clusters_p)
        # print(clusters_a[1].size)

        for job in clusters_p:
            job_id = int(job[0][0]) - 1
            job = job[:,-1]
            yields_p[job_id] = np.count_nonzero(job >= 2)

        for job in clusters_e:
            job_id = int(job[0][0]) - 1
            job = job[:,-1]
            yields_e[job_id] = np.count_nonzero(job >= 2)

        for job in clusters_a:
            job_id = int(job[0][0]) - 1
            job = job[:,-1]
            yields_a[job_id] = np.count_nonzero(job >= 2)

        # print(yields_p[:,0])

    # Scale by dose if requested
    if (scale_by_dose):
        yields_p = scale_species_yields_by_dose(yields_p,dataset,energy,'proton',target_dose,path_prefix,path_sufix,num_jobs,dose_column)
        yields_e = scale_species_yields_by_dose(yields_e,dataset,energy,'electron',target_dose,path_prefix,path_sufix,num_jobs,dose_column)
        yields_a = scale_species_yields_by_dose(yields_a,dataset,energy,'alpha',target_dose,path_prefix,path_sufix,num_jobs,dose_column)

    # Concatenate yields from all species together
    yields = yields_p+yields_e+yields_a

    return yields



def get_damage_yields_vs_energy(dataset, settings):
    energy_strings = settings['energy_strings']
    yields_per_energy = np.zeros((settings['num_jobs_n'],len(energy_strings)))

    for i in range(len(energy_strings)):        
        yields = import_yields( \
            settings['input_path_prefix'], \
            settings['input_path_sufix'], \
            dataset, \
            energy_strings[i], \
            settings['dose_scaling'], \
            settings['dose_target'], \
            settings['num_jobs_n'], \
            settings['num_yield_columns'], \
            settings['yield_column'], \
            settings['dose_column'])

        yields_per_energy[:,i] = yields #/6.32736 # normalize per bp

    mean_yield_per_energy = np.mean(yields_per_energy,axis=0)
    stderr_yield_per_energy = stats.sem(yields_per_energy,axis=0)

    # reference x-rays
    yields_ref = import_yields( \
        settings['input_path_prefix'], \
        settings['input_path_sufix'], \
        dataset, \
        'x250keV', \
        settings['dose_scaling'], \
        settings['dose_target'], \
        settings['num_jobs_x'], \
        settings['num_yield_columns'], \
        settings['yield_column'], \
        settings['dose_column'])
    
    yields_ref = yields_ref  #/6.32736 # normalize per bp
    mean_yield_ref = np.mean(yields_ref)
    stderr_yield_ref = stats.sem(yields_ref)
    
    return mean_yield_per_energy, stderr_yield_per_energy, mean_yield_ref, stderr_yield_ref

def get_yield_rbe_vs_energy(dataset, settings):
    mean_yield_per_energy, stderr_yield_per_energy, mean_yield_ref, stderr_yield_ref = get_damage_yields_vs_energy(dataset, settings)
    rbe_per_energy = mean_yield_per_energy / mean_yield_ref
    stderr_rbe_per_energy = get_propagated_error_div(mean_yield_per_energy, stderr_yield_per_energy, mean_yield_ref, stderr_yield_ref)
    return rbe_per_energy, stderr_rbe_per_energy

def get_propagated_error_mult(x_val, x_err, y_val, y_err):
    prod = x_val*y_val
    a = x_err/x_val
    b = y_err/y_val
    return np.abs(prod) * np.sqrt(np.power(a,2) + np.power(b,2))

def get_propagated_error_div(x_val, x_err, y_val, y_err):
    quot = x_val/y_val
    a = x_err/x_val
    b = y_err/y_val
    return np.abs(quot) * np.sqrt(np.power(a,2) + np.power(b,2))