#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import re 


# In[2]:


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


# In[3]:


def get_block_end(beginning, FullBreak):
    # Split the string by " / " to separate the groups
    groups = FullBreak.split(" / ")
    
    # Get the last group (will also work if there is only one group)
    last_group = groups[-1]
    
    # Split the last group by "," to get the numbers
    numbers_in_last_group = last_group.split(", ")
    
    # Get the first number of the last group
    first_number = numbers_in_last_group[0]
    
    return int(first_number) + beginning - 1


# In[4]:


def Count_ComplexClusters(table):
    ComplexDSBs = 0
    i = 0 
    Damages = table[i][3:6]
    while i < len(table)-1:
        ChromoID = table[i][0]
        BP_start = table[i][1]
        BP_end = table[i][2]

        NextChromoID = table[i+1][0]
        NextBP_start = table[i+1][1]
        NextBP_end = table[i+1][2]
        NextDamages = table[i+1][3:6]

        if ChromoID == NextChromoID and (NextBP_start-BP_end)<40:
            Damages += NextDamages
        else:
            if Damages[2] > 0 and np.sum(Damages) > 1:
                ComplexDSBs += 1
            Damages = NextDamages
        i += 1 
        
    return ComplexDSBs


# In[15]:


def Count_BaioccoCluster(table):
    Baiocco_cluster = 0
    i = 0
    while i < len(table)-1:
        ChromoID = table[i][0]
        BP_start = table[i][1]
        BP_end = table[i][2]

        NextChromoID = table[i+1][0]
        NextBP_start = table[i+1][1]
        NextBP_end = table[i+1][2]

        if ChromoID == NextChromoID and (NextBP_start-BP_end)< 25:
            Baiocco_cluster += 1 
            i += 2
        else:
            i += 1
    return Baiocco_cluster
        


def GeoCluster(SDDFilePath,eps):
    
    blocks = read_after_header(SDDFilePath)
    n = len(blocks)
    
    x, y, z = np.zeros((3,n))
    for i in range(n):
        if re.compile(";").split(blocks[i])[5][7] != "0":
            v = blocks[i].split(";")[1].split(" / ")[1].split(", ")
            x[i] = v[0]
            y[i] = v[1]
            z[i] = v[2]
        else:
            x[i] = np.nan
            y[i] = np.nan
            z[i] = np.nan

    x = x[np.logical_not(np.isnan(x))]
    y = y[np.logical_not(np.isnan(y))]
    z = z[np.logical_not(np.isnan(z))]

    nDSB = len(x)
    d, clusters_pos, SumVector = np.zeros((3,nDSB))
    
    for i in range(nDSB):
        clusters_pos = np.zeros(nDSB)
        
        d = np.sqrt((x-x[i])**2 + (y-y[i])**2 + (z-z[i])**2)
        clusters_pos[np.where(d<eps)] = 1
        
        SumVector += clusters_pos
        
    return (np.sum(SumVector)-nDSB)/2


def clusterer(SDD_file_path,eps):

    EndPoints = np.zeros(3+len(eps))
    blocks = read_after_header(SDD_file_path)
    Nblocks = len(blocks)
    block_table = np.zeros((Nblocks,6), dtype = "int")
    

    for i, block in enumerate(blocks):
        parts = block.split(";")
        block_table[i, 0] = int(parts[2].split(',')[1])

        first_bp = int(parts[3].split(',')[0])
        block_table[i, 1] = first_bp
        block_table[i, 2] = get_block_end(first_bp, parts[6])

        block_table[i, 3:6] = np.fromstring(parts[5], dtype=int, sep=',')
        
    sorted_indices = np.lexsort((block_table[:, 1], block_table[:, 0]))
    block_table = block_table[sorted_indices]
        
    
    ComplexDSB_table = block_table.copy()
    ComplexDSB_table[:,4] = ComplexDSB_table[:,4]-2*ComplexDSB_table[:,5]
    
    Baiocco_table = block_table.copy()
    Baiocco_table = np.array([arr for arr in Baiocco_table if arr[-1] != 0])

    EndPoints[0] = len(Baiocco_table)
        
    EndPoints[1] = Count_ComplexClusters(ComplexDSB_table)

    EndPoints[2] = Count_BaioccoCluster(Baiocco_table)
    
    k = 3 
    for d in eps:
        EndPoints[k] = GeoCluster(SDD_file_path,d)
        k += 1
        
    return  EndPoints
        
        
        

