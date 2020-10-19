'''Combines NJOY raw data into comprehensible transport data'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#====================================================================
def BuildCombinedData(raw_njoy_data, plot=False):
    ''' Combines enjoy raw data into a dictionary of vectors and matrices '''
    # ================================= Get dictionaries
    group_structures = raw_njoy_data["group_structures"]
    cross_sections = raw_njoy_data["cross_sections"]
    transfer_matrices = raw_njoy_data["transfer_matrices"]
    
    # ================================= Determine # of groups
    neutn_gs = group_structures["neutron"]
    gamma_gs = group_structures["gamma"] if "gamma" in group_structures else []
    
    G_n = len(neutn_gs)
    G_g = len(gamma_gs)
    G   = G_n + G_g 
      
    # ================================= Combine sig_t
    sig_t = np.zeros(G)
    
    if ("(n,total)" in cross_sections):
        sig_t_ndata = cross_sections["(n,total)"]
        for entry in sig_t_ndata:
            g = entry[0]
            v = entry[1]
            sig_t[G_n-g-1] += v
    
    if ("(g,total)" in cross_sections):
        sig_t_gdata = cross_sections["(g,total)"]
        for entry in sig_t_gdata:
            g = entry[0]
            v = entry[1]
            sig_t[G_n + G_g-g-1] += v
    
    # ================================= Inverse velocity term
    inv_v = np.zeros(G)
    
    if ("inv_velocity" in cross_sections):
        inv_v_ndata = cross_sections["inv_velocity"]
        for entry in inv_v_ndata:
            g = entry[0]
            v = entry[1]
            inv_v[G_n-g-1] += v
    
    # ================================= Fission data
    sig_f = np.zeros(G)
    if ("(n,fission)" in cross_sections):
        sig_f_data = cross_sections["(n,fission)"]
        for entry in sig_f_data:
            g = entry[0]
            v = entry[1]
            sig_f[G_n-g-1] += v
    
    nu_total = np.zeros(G)
    if ("total_nubar" in cross_sections):
        total_nubar_data = cross_sections["total_nubar"]
        for entry in total_nubar_data:
            g = entry[0]
            v = entry[1]
            nu_total[G_n-g-1] += v
    
    nu_prompt = np.zeros(G)
    if ("prompt_nubar" in cross_sections):
        prompt_nubar_data = cross_sections["prompt_nubar"]
        for entry in prompt_nubar_data:
            g = entry[0]
            v = entry[1]
            nu_prompt[G_n-g-1] += v
    
    chi_prompt = np.zeros(G)
    if ("prompt_chi" in cross_sections):
        prompt_chi_data = cross_sections["prompt_chi"]
        for entry in prompt_chi_data:
            g = entry[0]
            v = entry[1]
            chi_prompt[G_n-g-1] += v
    
    # ================================= Delayed neutron data
    decay_const=[]
    if ("decay_constants" in cross_sections):
        for entry in cross_sections["decay_constants"]:
            decay_const.append(entry[1])
        decay_const = np.asarray(decay_const)
    J = len(decay_const) # number of precursor groups
    
    
    nu_delayed = np.zeros(G)
    if ("delayed_nubar" in cross_sections):
        delayed_nubar_data = cross_sections["delayed_nubar"]
        for entry in delayed_nubar_data:
            g = entry[0]
            v = entry[1]
            nu_delayed[G_n-g-1] += v
    
    chi_delayed = np.zeros((G,J))
    if ("delayed_chi" in cross_sections):
        delayed_chi_data = cross_sections["delayed_chi"]
        for entry in delayed_chi_data:
            g = entry[0]
            v = entry[1:]
            chi_delayed[G_n-g-1] += v
    
    gamma = np.zeros(J)
    if (np.sum(nu_delayed)>0 and np.sum(chi_delayed)>0):
        gamma = np.sum(chi_delayed,axis=0)
    chi_delayed /= gamma
    
    # ================================= Combine transfer matrices
    # Keys available to all isotopes
    n_to_n_elastic_keys = ["(n,elastic)"]
    n_to_g_transfer_keys = [ "(n,g)", "(n,inel)", "(n,np)", 
                             "(n,nd)", "(n,p)", "(n,d)", 
                             "(n,t)", "(n,a)"]
    g_to_g_transfer_keys = ["(g,coherent)", "(g,incoherent)", \
                            "(g,pair_production)"]
    n_to_n_freegas_keys = ["mt221"]
    
    # Keys for S(alpha, beta) thermal corrections. Note that this 
    # will NOT be correct for inelastic scattering if data for 
    # H2O and ZrH exist in H1 files. As it stands, all reaction 
    # keys are additive. For isotopes with only one S(alpha, beta) 
    # material, all other reactions will be skipped over because 
    # those reactions do not live in the file. A flag could be used to 
    # choose the correct reaction numbers. 
    
    # graphite, H in ZrH, Zr in ZrH
    n_to_n_sab_elastic_keys = ["mt230", "mt226", "mt236"] 
    
    # H in H2O, graphite, H in ZrH, Zr in ZrH
    n_to_n_sab_inelastic_keys = ["mt222", "mt229", "mt226", "mt235"]
    
    # # graphite, H in ZrH, Zr in ZrH
    # n_to_n_sab_elastic_keys = ["mt250"] 
    # # H in H2O, graphite, H in ZrH, Zr in ZrH
    # n_to_n_sab_inelastic_keys = ["mt249"]
    

    # ===== Get the transfer matrices
    # Adding all the elastic scattering data
    nranges_to_nranges_elastic = []
    for rxn in n_to_n_elastic_keys:
        if rxn in transfer_matrices:
            nranges_to_nranges_elastic.append(transfer_matrices[rxn])
    
    # Adding all the (n,nxx) data, inelastic data
    nranges_to_nranges_inelastic = []
    for nn in range(1,24+1):
        rx_name = "(n,n{:02d})".format(nn)
        if rx_name in transfer_matrices:
            mat = transfer_matrices[rx_name]
            nranges_to_nranges_inelastic.append(mat)
    
    # Adding all the free gas elastic scattering data
    nranges_to_nranges_freegas = []
    for rxn in n_to_n_freegas_keys:
        if rxn in transfer_matrices:
            nranges_to_nranges_freegas.append(transfer_matrices[rxn])
    
    # Adding all the elastic scattering S(\alpha,\beta) data
    nranges_to_nranges_sab_el = []
    for rxn in n_to_n_sab_elastic_keys:
        if rxn in transfer_matrices:
            nranges_to_nranges_sab_el.append(transfer_matrices[rxn])
    
    # Adding all the inelastic scattering S(\alpha,\beta) data
    nranges_to_nranges_sab_inel = []
    for rxn in n_to_n_sab_inelastic_keys:
        if rxn in transfer_matrices:
            nranges_to_nranges_sab_inel.append(transfer_matrices[rxn])
    
    # Adding all the neutron to gamma data
    nranges_to_granges = []
    for rxn in n_to_g_transfer_keys:
        if rxn in transfer_matrices:
            nranges_to_granges.append(transfer_matrices[rxn])
    
    # Adding all the gamma to gamma data
    granges_to_granges = []
    for rxn in g_to_g_transfer_keys:
        if rxn in transfer_matrices:
            granges_to_granges.append(transfer_matrices[rxn])
    
    # ===== Computing the max number of moments
    max_num_moms = 0
    for range_data in nranges_to_nranges_elastic:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_nranges_inelastic:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_nranges_freegas:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_nranges_sab_el:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_nranges_sab_inel:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_granges:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in granges_to_granges:
        if range_data:
            max_num_moms = max(max_num_moms,len(range_data[0])-2)
    
    # ===== Initializing the transfer matrices
    transfer_mats = []
    for m in range(0,max_num_moms):
        transfer_mats.append(np.zeros((G,G)))
    transfer_mats_nonzeros = []
    
    #=======================================
    # Lambda-ish to add neutron data
    def AddTransferNeutron(data_vals,offset=0):
        for entry in data_vals:
            num_moms = len(entry)-2
            gprime = G_n - entry[0] - 1
            g      = G_n - entry[1] - 1 + offset
            
            for m in range(0,num_moms):
                v = entry[m+2]
                transfer_mats[m][gprime,g] += v
    
    #=======================================
    # Lambda-ish to add gamma data
    def AddTransferGamma(data_vals):
      # (g,coherent)
        for entry in data_vals:
            gprime = G_n + G_g - entry[0] - 1
            g      = G_n + G_g - entry[1] - 1
        
            # Determine number of moments
            num_moms = len(entry)-2
            for m in range(0,num_moms):
                v = entry[m+2]
                transfer_mats[m][gprime,g] += v
    
    # ===== Format the transfer matrices, store scattering
    transfer_mats_standard       = np.copy(transfer_mats)
    # transfer_mats_standard_el    = np.copy(transfer_mats)
    # transfer_mats_standard_inel  = np.copy(transfer_mats)
    # transfer_mats_standard_gamma = np.copy(transfer_mats)
    transfer_mats_freegas        = np.copy(transfer_mats)
    transfer_mats_sab_el         = np.copy(transfer_mats)
    transfer_mats_sab_inel       = np.copy(transfer_mats)
    
    # Regular elastic scatter
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_elastic:
        AddTransferNeutron(range_data)
    sig_s_el = np.sum(transfer_mats[0],axis=1)
    # transfer_mats_standard_el = np.copy(transfer_mats)
    transfer_mats_standard = np.copy(transfer_mats)
    
    # Regular inelastic scatter
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_inelastic:
        AddTransferNeutron(range_data)
    sig_s_inel = np.sum(transfer_mats[0],axis=1)
    # transfer_mats_standard_inel = np.copy(transfer_mats)
    for m in range(0,max_num_moms):
        transfer_mats_standard[m] += transfer_mats[m]
    
    # Regular n,\gamma transfer
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_granges:
        AddTransferNeutron(range_data,offset=G_g)
    for m in range(0,max_num_moms):
        transfer_mats_standard[m] += transfer_mats[m]
    
    # Regular \gamma,\gamma transfer
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in granges_to_granges:
        AddTransferGamma(range_data)
    for m in range(0,max_num_moms):
        transfer_mats_standard[m] += transfer_mats[m]
    sig_s_uncorr = np.sum(transfer_mats_standard[0],axis=1)    
    
    # Freegas elastic scatter
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_freegas:
        AddTransferNeutron(range_data)
    sig_s_freegas = np.sum(transfer_mats[0],axis=1)
    transfer_mats_freegas = np.copy(transfer_mats)
    
    # S(alpha,beta) elastic scatter
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_sab_el:
        AddTransferNeutron(range_data)
    sig_s_el_sab = np.sum(transfer_mats[0],axis=1)
    transfer_mats_sab_el = np.copy(transfer_mats)
    
    # S(alpha,beta) inelastic scatter
    transfer_mats = [0.0*mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_sab_inel:
        AddTransferNeutron(range_data)
    sig_s_inel_sab = np.sum(transfer_mats[0],axis=1)
    transfer_mats_sab_inel = np.copy(transfer_mats)
    
    # ===== Store stuff
    # Uncorrected quantities
    sig_t_uncorr = sig_t
    sig_a = sig_t_uncorr - sig_s_uncorr
    # Correction terms
    sig_s_sab = sig_s_el_sab + sig_s_inel_sab
    
    # ===== Make the termal scattering corrections
    # This is a bit complex. First, the freegas transfer
    # matrices are overlaid on the standard transfer matrices.
    # Then, the same process is carried out with the 
    # inelastic S(alpha, beta) terms. Finally, the elastic
    # S(alpha, beta) terms are added to the result.
    transfer_mats = np.copy(transfer_mats_standard)
    m0_transfer_mat_freegas = transfer_mats_freegas[0]
    m0_transfer_mat_sab_inel = transfer_mats_sab_inel[0]
    for m in range(0,max_num_moms):
        m_transfer_mat_freegas = transfer_mats_freegas[m]
        m_transfer_mat_sab_inel = transfer_mats_sab_inel[m]
        for gprime in range(0,G):
            for g in range(0,G):
                if (np.abs(m0_transfer_mat_freegas[gprime,gprime]) > 1.0e-18):
                    transfer_mats[m][gprime,g] = m_transfer_mat_freegas[gprime,g]
                if (np.abs(m0_transfer_mat_sab_inel[gprime,gprime]) > 1.0e-18):
                    transfer_mats[m][gprime,g] = m_transfer_mat_sab_inel[gprime,g]
        transfer_mats[m] += transfer_mats_sab_el[m]

    # ===== Update cross section vectors
    sig_s = np.sum(transfer_mats[0],axis=1)
    sig_t = sig_a + sig_s
        
    # ===== Determine sparsity of the transfer matrices
    for m in range(0,max_num_moms):
        mat_non_zeros = []
        for gprime in range(0,G):
            non_zeros = []
            for g in range(0,G):
                if (abs(transfer_mats[m][gprime,g]) > 1.0e-18):
                    non_zeros.append(g)
            mat_non_zeros.append(non_zeros)
        transfer_mats_nonzeros.append(mat_non_zeros)
    
    # ===== Plot the matrix
    if plot:
        Atest = np.copy(transfer_mats[0])
        nz = np.nonzero(Atest)
        Atest[nz] = np.log10(Atest[nz]) + 10.0
        
        plt.figure(figsize=(6,6))
        im = plt.imshow(Atest, cmap=cm.Greys)
        plt.xticks(np.arange(0,G,10), [str(g) for g in range(0,G,10)])
        plt.yticks(np.arange(0,G,10), [str(g) for g in range(0,G,10)])
        plt.xlabel('Destination energy group')
        plt.ylabel('Source energy group')
        plt.gca().xaxis.set_ticks_position('top')
        plt.gca().xaxis.set_label_position('top')
        plt.show()
        # plt.savefig("SERPENTTransferMatrix.png")

        #================================== Build group structures
        np_neutn_gs = np.matrix(neutn_gs)
        nbin_lo = np_neutn_gs[:,1]
        nbin_hi = np_neutn_gs[:,2]
        nbin_center = (0.5*(nbin_lo+nbin_hi))[::-1]
        
        #================================== Plot Crosssections
        plt.figure(figsize=(6,6))
        plt.semilogx(nbin_center,sig_t,label=r"$\sigma_t$")
        plt.semilogx(nbin_center,sig_a,"r",label=r"$\sigma_a$")
        
        plt.semilogx(nbin_center,sig_s_freegas,label=r"$\sigma_s$ freegas")
        plt.semilogx(nbin_center,sig_s_uncorr,label=r"$\sigma_s$ elastic unbound")
        plt.semilogx(nbin_center,sig_s,label=r"$\sigma_s$ total SAB")
        plt.semilogx(nbin_center,sig_s_el_sab,label=r"$\sigma_s$ elastic SAB")
        plt.semilogx(nbin_center,sig_s_inel_sab,label=r"$\sigma_s$ inelastic SAB")
        plt.legend()
        plt.show()


    
    #================================== Build return data
    return_data = {}
    return_data["neutron_gs"] = neutn_gs
    return_data["gamma_gs"] = gamma_gs
    return_data["sigma_t"] = sig_t
    return_data["sigma_t_uncorr"] = sig_t_uncorr
    return_data["sigma_a"] = sig_a
    return_data["sigma_s"] = sig_s
    return_data["sigma_s_uncorr"] = sig_s_uncorr
    return_data["sigma_s_freegas"] = sig_s_freegas
    return_data["sigma_s_sab"] = sig_s_sab
    return_data["sigma_s_el"] = sig_s_el
    return_data["sigma_s_inel"] = sig_s_inel
    return_data["sigma_s_sab_el"] = sig_s_el_sab
    return_data["sigma_s_sab_inel"] = sig_s_inel_sab
    return_data["sigma_f"] = sig_f
    return_data["nu_total"] = nu_total
    return_data["nu_prompt"] = nu_prompt
    return_data["nu_delayed"] = nu_delayed
    return_data["chi_prompt"] = chi_prompt
    return_data["chi_delayed"] = chi_delayed
    return_data["decay_constants"] = decay_const
    return_data["gamma"] = gamma
    return_data["inv_velocity"] = inv_v
    return_data["transfer_matrices"] = transfer_mats
    return_data["transfer_matrices_sparsity"] = transfer_mats_nonzeros
    return return_data 
