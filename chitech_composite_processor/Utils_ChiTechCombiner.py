"""Combined Chi-tech data from multiple files to form composite cross sections"""

# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 11:38:48 2021

@author: Asus
"""

import numpy as np

def BuildCombinedChiTechData (chi_list, N_density):
    
    data = {}
    
    cs_OverallList = []
    neutron_gs = []
    gamma_gs = []
    for i in range (len(chi_list)):
        cf = open(chi_list[i], 'r')
        cross_section = []
        neutron_gs = []
        gamma_gs = []
        for line in cf:
            # Get the group structure
            if line.find('GROUPS') != -1:
                line_list  = line.split(' ')
                if i == 0:
                    num_group = int(line_list[1])
                else:
                    if (int(line_list[1]) != num_group):
                        err_msg = "The group_structure is not the same in " + str(chi_list[i]) + line_list[1] + str(num_group)
                        raise ValueError(err_msg)
            
            # Get the group structures
            elif line.find('GS_BEGIN') != -1:
                line_list = line.split('_')
                particle = line_list[0]
                line_string = cf.readline()
                while line_string.find('END') == -1:
                    gs_data = []
                    line_list = line_string.split()
                    if line_list[0].find('GS_END') == -1:
                        for i in range (len(line_list)):
                            gs_data.append(float(line_list[i]))
                        if particle == 'NEUTRON':
                            neutron_gs.append(gs_data)
                        elif particle == 'GAMMA':
                            gamma_gs.append(gs_data)
                    line_string = cf.readline()
                        
                
            # Get the moment
            elif ((line.find('MOMENTS') !=  -1) and (line.find('TRANSFER') == -1)):
                line_list = line.split(' ')
                if i == 0:
                    num_moment = int(line_list[1])
                else:
                    if (int(line_list[1]) != num_moment):
                        err_msg = "The group moment is not the same in " + str(chi_list[i]) + line_list[1] + str(num_moment)
                        raise ValueError(err_msg)
            
            # Get the cross_section
            elif ((line.find('TRANSFER') == -1) and (line.find('BEGIN') != -1)):
                line_list = line.split('_')
                cs_name = line_list[0] + '_' + line_list[1]
                cross_section.append(cs_name)
        cs_OverallList.append(cross_section)
        cf.close() 
    
    # Flip the neuton and gamma gs 
    # neutron_gs = np.flip(neutron_gs)
    # gamma_gs = np.flip(gamma_gs)

    data["G"] = num_group
    data["M"] = num_moment
    
    data["neutron_gs"] = neutron_gs
    data["gamma_gs"] = gamma_gs
    
    #==================================== Find all cross section labels
    union = set.union(*[set(x) for x in cs_OverallList])
   
    #===================================  Process Cross Sections
    cs_dict = {}
    for cs_name in union:  
        cs_value = [] 
        cs_value2 = np.zeros(num_group)
        for i in range (len(chi_list)):
            cf = open(chi_list[i], 'r')
            for line in cf:
                if ((line.find(cs_name) != -1) and (line.find('BEGIN') != -1)):
                    for g in range (0,num_group):
                        value_line = cf.readline().split()
                        if value_line != []:
                            cs_value2[g] += float(value_line[1])*N_density[i]
            cf.close()
        
        cs_dict[cs_name.lower()] = cs_value2
        
    for key in cs_dict:
        print(key, len(cs_dict[key]))
    data.update(cs_dict)
    
    #======================================= Process Transfer Matrix
    # Create an empty space to store the transfer_matrices
    transfer_mats = []
    for m in range (0, num_moment):
        transfer_mats.append(np.zeros((num_group,num_group)))
    
    for i in range (len(chi_list)):
        m = 0
        cf = open(chi_list[i], 'r')
        for line in cf:
            if line.find('MOMENTS_BEGIN') != -1:
                moment_line = cf.readline().split()
                while (moment_line[0] != 'TRANSFER_MOMENTS_END'):
                    moment_line = cf.readline().split()
                    if ((moment_line[0] == '#') or (moment_line[0] == 'TRANSFER_MOMENTS_END')):
                        m += 1
                    else:
                        # Replace values in transfer matrixes
                        XS_value = float(moment_line[-1])*N_density[i]
                        moment_index  = int(moment_line[1])
                        group_index   = int(moment_line[2])
                        nonzero_index = int(moment_line[3])
                        transfer_mats[moment_index][group_index,nonzero_index] += XS_value
        cf.close()
    
    data["transfer_matrices"] = transfer_mats
    
    return data


"""
    #======================================= Process Transfer Matrix
    # Create an empty space to store the transfer_matrices
    transfer_mats = []
    for m in range (0, num_moment):
        moment_list = []
        for g in range (0, num_group):
            null_array = np.zeros(num_group)
            moment_list.append(null_array)
        transfer_mats.append(moment_list)
    
    # Create an empty space to store the transfer_matrices
    transfer_mats_nonzeros = []
    for m in range (0, num_moment):
        moment_list = []
        for g in range (0, num_group):
            null_list = []
            moment_list.append(null_list)
        transfer_mats_nonzeros.append(moment_list)

    total_moment = []
    for i in range (len(chi_list)):
        element_moment = []
        moment_value = []
        group_moment = []
        m = 0
        g = 0
        cf = open(chi_list[i], 'r')
        for line in cf:
            if line.find('MOMENTS_BEGIN') != -1:
                moment_line = cf.readline().split()
                while (moment_line[0] != 'TRANSFER_MOMENTS_END'):
                    moment_line = cf.readline().split()
                    if ((moment_line[0] == '#') or (moment_line[0] == 'TRANSFER_MOMENTS_END')):
                        moment_value.append(group_moment)
                        group_moment = []
                        element_moment.append(moment_value)
                        moment_value = []
                        m += 1
                        g = 0
                    else:
                        if (int(moment_line[2]) != g):
                            moment_value.append(group_moment)
                            group_moment = []
                            g += 1
                        
                        # Replace values in transfer matrixes
                        XS_value = float(moment_line[-1])*N_density[i]
                        moment_index = int(moment_line[1])
                        group_index = int(moment_line[2])
                        nonzero_index = int(moment_line[3])
                        transfer_mats[moment_index][group_index][nonzero_index] += XS_value
                        
                        # Save the index for non_zeros values in transfer matrixes
                        transfer_mats_nonzeros[moment_index][group_index].append(nonzero_index)
                        group_moment.append(nonzero_index)
        total_moment.append(element_moment)
                    
                                          
        cf.close()
    
    # Eliminate duplicates in the non-zero transfer matrices
    for m in range (0, num_moment):
        for g in range (0, num_group):
            transfer_mats_nonzeros[m][g] = list(set(transfer_mats_nonzeros[m][g]))
    
    data["transfer_matrices"] = transfer_mats
    data["transfer_matrices_sparsity"] = transfer_mats_nonzeros
    
    return data
"""