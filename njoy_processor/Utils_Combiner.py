"""Combines NJOY raw data into comprehensible transport data"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ====================================================================
def BuildCombinedData(raw_njoy_data, plot = False, verbose = False):
    """ Combines enjoy raw data into a dictionary of vectors and matrices """
    # ================================= Get dictionaries
    group_structures = raw_njoy_data["group_structures"]
    cross_sections = raw_njoy_data["cross_sections"]
    transfer_matrices = raw_njoy_data["transfer_matrices"]
    print(cross_sections.keys())
    
    # ================================= Determine # of groups
    neutn_gs = group_structures["neutron"] if "neutron" in group_structures else []
    gamma_gs = group_structures["gamma"] if "gamma" in group_structures else []

    G_n = len(neutn_gs)
    G_g = len(gamma_gs)
    G = G_n + G_g

    with_sab = False

    # ================ Get a description of the problem
    problem_description = {}
    if neutn_gs != []:
        if gamma_gs != []:
            problem_description['problem_type'] = "Neutron + Gamma"
        elif gamma_gs == []:
            problem_description['problem_type'] = "Neutron only"
    elif gamma_gs != []:
        problem_description['problem_type'] = "Gamma only"
    else:
        problem_description['problem_type'] ="Unknown"
        print(problem_description)
        raise Exception('unknown particle type')
    
    problem_description['G_n'] = G_n
    problem_description['G_g'] = G_g

    # ================================= Combine sig_t
    sig_t = np.zeros(G)
    if "(n,total)" in cross_sections:
        data = cross_sections["(n,total)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_t[G_n - g - 1] += v

    if "(g,total)" in cross_sections:
        data = cross_sections["(g,total)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_t[G_n + G_g - g - 1] += v
    
    # ================================= Combine Heating Cross-Sections
    sig_heat = np.zeros(G)

    if "(n,heat)" in cross_sections:
        data = cross_sections["(n,heat)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_heat[G_n - g - 1] += v

    if "(g,heat)" in cross_sections:
        data = cross_sections["(g,heat)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_heat[G_n + G_g - g - 1] += v

    # ================================= Scattering terms
    sig_el = np.zeros(G)
    if "(n,elastic)" in cross_sections:
        data = cross_sections["(n,elastic)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_el[G_n - g - 1] += v
    
    sig_inel = np.zeros(G)
    if "(n,inel)" in cross_sections:
        data = cross_sections["(n,inel)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_inel[G_n - g - 1] += v

    sig_freegas = np.zeros(G)
    if "free_gas" in cross_sections:
        data = cross_sections["free_gas"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_freegas[G_n - g - 1] += v

    sig_el_sab = np.zeros(G)
    if "elastic_s(a,b)" in cross_sections:
        with_sab = True
        data = cross_sections["elastic_s(a,b)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_el_sab[G_n - g - 1] += v

    sig_inel_sab = np.zeros(G)
    if "inelastic_s(a,b)" in cross_sections:
        with_sab = True
        data = cross_sections["inelastic_s(a,b)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_inel_sab[G_n - g - 1] += v

    sig_nxn = np.zeros(G)
    for nn in range(2, 4 + 1):
        rx = "(n,{:01d}n)".format(nn)
        if rx in cross_sections:
            data = cross_sections[rx]
            for entry in data:
                g = entry[0]
                v = entry[1]
                sig_nxn[G_n - g - 1] += v

    sig_coht = np.zeros(G_g)
    if "(g,coherent)" in cross_sections:
        data = cross_sections["(g,coherent)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_coht[G_g - g - 1] += v

    sig_incoht = np.zeros(G_g)
    if "(g,incoherent)" in cross_sections:
        data = cross_sections["(g,incoherent)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_incoht[G_g - g - 1] += v


    # ================================= Inverse velocity term
    inv_v = np.zeros(G)
    if "inv_velocity" in cross_sections:
        data = cross_sections["inv_velocity"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            inv_v[G_n - g - 1] += v

    # ================================= Fission data
    sig_f = np.zeros(G)
    if "(n,fission)" in cross_sections:
        data = cross_sections["(n,fission)"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            sig_f[G_n - g - 1] += v

    nu_total = np.zeros(G)
    if "total_nubar" in cross_sections:
        data = cross_sections["total_nubar"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            nu_total[G_n - g - 1] += v

    nu_prompt = np.zeros(G)
    if "prompt_nubar" in cross_sections:
        data = cross_sections["prompt_nubar"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            nu_prompt[G_n - g - 1] += v

    chi_prompt = np.zeros(G)
    if "prompt_chi" in cross_sections:
        data = cross_sections["prompt_chi"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            chi_prompt[G_n - g - 1] += v

    # ================================= Delayed neutron data
    decay_const = []
    if "decay_constants" in cross_sections:
        data = cross_sections["decay_constants"]
        for entry in data:
            decay_const.append(entry[1])
        decay_const = np.asarray(decay_const)
    J = len(decay_const)  # number of precursor groups

    nu_delayed = np.zeros(G)
    if "delayed_nubar" in cross_sections:
        data = cross_sections["delayed_nubar"]
        for entry in data:
            g = entry[0]
            v = entry[1]
            nu_delayed[G_n - g - 1] += v

    chi_delayed = np.zeros((G, J))
    if "delayed_chi" in cross_sections:
        data = cross_sections["delayed_chi"]
        for entry in data:
            g = entry[0]
            v = entry[1:]
            chi_delayed[G_n - g - 1] += v

    gamma = np.zeros(J)
    if np.sum(nu_delayed) > 0 and np.sum(chi_delayed) > 0:
        gamma = np.sum(chi_delayed, axis = 0)
    chi_delayed /= gamma

    # ================================= Combine transfer matrices
    # Normal elastic scattering
    n_to_n_el_keys = ["(n,elastic)"]
    # Normal inelastic scattering
    n_to_n_inel_keys = []
    for nn in range(1, 40 + 1):
        rx_name = "(n,n{:02d})".format(nn)
        n_to_n_inel_keys.append(rx_name)
    n_to_n_inel_keys.append("(n,nc)")
    # (n,xn) reactions
    n_to_n_nxn_keys = []
    for nn in range(2, 4 + 1):
        rx_name = "(n,{:01d}n)".format(nn)
        n_to_n_nxn_keys.append(rx_name)
    # Freegas thermal scattering
    n_to_n_freegas_keys = ["mt221"]
    # S(alpha,beta) thermal scattering keys
    n_to_n_sab_inel_keys = ["mt222", "mt229",
                            "mt225", "mt235"]
    n_to_n_sab_el_keys = ["mt230", "mt226", "mt236"]
    # Neutron to neutron other reactions
    n_to_n_other_keys = ["(n,np)", "(n,nd)", "(n,na)",
                         "(n,n2a)", "(n,n3a)"]
    # Neutron to gamma reactions
    n_to_g_transfer_keys = ["(n,g)", "(n,inel)",
                            "(n,np)", "(n,nd)",
                            "(n,p)", "(n,d)", "(n,t)", "(n,a)", "(n,2a)",  
                            "(n,na)", "(n,n2a)", "(n,n3a)",
                            "(n,nc)"]
    # Gamma to gamma reactions
    g_to_g_transfer_keys = ["(coherent)", "(incoherent)",
                            "(pair_production)"]

    # ================================= Get the transfer matrices
    # Adding all the elastic scattering data
    nranges_to_nranges_el = []
    for rxn in n_to_n_el_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_el.append(mat)

    # Adding all the inelastic data
    nranges_to_nranges_inel = []
    for rxn in n_to_n_inel_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_inel.append(mat)

    # Adding all the (n,xn) data
    nranges_to_nranges_nxn = []
    for rxn in n_to_n_nxn_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_nxn.append(mat)

    # Adding all the free gas elastic scattering data
    nranges_to_nranges_freegas = []
    for rxn in n_to_n_freegas_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_freegas.append(mat)

    # Adding all the elastic scattering S(\alpha,\beta) data
    nranges_to_nranges_sab_el = []
    for rxn in n_to_n_sab_el_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_sab_el.append(mat)
            with_sab = True

    # Adding all the inelastic scattering S(\alpha,\beta) data
    nranges_to_nranges_sab_inel = []
    for rxn in n_to_n_sab_inel_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_sab_inel.append(mat)
            with_sab = True

    # Adding all the other neutron to neutron data
    nranges_to_nranges_other = []
    for rxn in n_to_n_other_keys:
        if rxn in transfer_matrices['neutron']:
            mat = transfer_matrices['neutron'][rxn]
            nranges_to_nranges_other.append(mat)

    # Adding all the neutron to gamma data
    nranges_to_granges = []
    for rxn in n_to_g_transfer_keys:
        if rxn in transfer_matrices['gamma']:
            nranges_to_granges.append(transfer_matrices['gamma'][rxn])

    # Adding all the gamma to gamma data
    granges_to_granges = []
    for rxn in g_to_g_transfer_keys:
        if rxn in transfer_matrices['gamma']:
            granges_to_granges.append(transfer_matrices['gamma'][rxn])

    # ===== Computing the max number of moments
    max_num_moms = 0
    for range_data in nranges_to_nranges_el:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_nranges_inel:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_nranges_freegas:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_nranges_sab_el:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_nranges_sab_inel:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_nranges_other:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in nranges_to_granges:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)
    for range_data in granges_to_granges:
        if range_data:
            max_num_moms = max(max_num_moms, len(range_data[0]) - 2)

    # ===== Initializing the transfer matrices
    transfer_mats = []
    for m in range(0, max_num_moms):
        transfer_mats.append(np.zeros((G, G)))

    # =======================================
    # Lambda-ish to add neutron data
    # Always operates on transfer_mats
    def AddTransferNeutron(data_values, offset = 0):
        for data_value in data_values:
            gprime = G_n - data_value[0] - 1
            g = G_n - data_value[1] - 1 + offset

            num_moms = len(data_value) - 2
            for m in range(0, num_moms):
                v = data_value[m + 2]
                transfer_mats[m][gprime, g] += v

    # =======================================
    # Lambda-ish to add gamma data
    # Always operates on transfer_mats
    def AddTransferGamma(data_values):
        # (g,coherent)
        for data_value in data_values:
            gprime = G_n + G_g - data_value[0] - 1
            g = G_n + G_g - data_value[1] - 1

            # Determine number of moments
            num_moms = len(data_value) - 2
            for m in range(0, num_moms):
                v = data_value[m + 2]
                transfer_mats[m][gprime, g] += v

    # ===== Format the transfer matrices, store scattering
    transfer_el = np.copy(transfer_mats)
    transfer_inel = np.copy(transfer_mats)
    transfer_nxn = np.copy(transfer_mats)
    transfer_sab = np.copy(transfer_mats)
    transfer_sab_el = np.copy(transfer_mats)
    transfer_sab_inel = np.copy(transfer_mats)
    transfer_freegas = np.copy(transfer_mats)
    transfer_n_other = np.copy(transfer_mats)

    # Regular elastic scatter (MT-2)
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_el:
        AddTransferNeutron(range_data)
    transfer_el = np.copy(transfer_mats)

    # Regular inelastic scattering
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_inel:
        AddTransferNeutron(range_data)
    transfer_inel = np.copy(transfer_mats)

    # (n,xn) reactions
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_nxn:
        AddTransferNeutron(range_data)
    transfer_nxn = np.copy(transfer_mats)

    # Freegas thermal scattering
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_freegas:
        AddTransferNeutron(range_data)
    transfer_freegas = np.copy(transfer_mats)

    # Elastic S(alpha, beta)
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_sab_el:
        AddTransferNeutron(range_data)
    transfer_sab_el = np.copy(transfer_mats)

    # Inelastic S(alpha, beta)
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_sab_inel:
        AddTransferNeutron(range_data)
    transfer_sab_inel = np.copy(transfer_mats)

    # (n,n+other) reactions
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for range_data in nranges_to_nranges_other:
        AddTransferNeutron(range_data)
    transfer_n_other = np.copy(transfer_mats)

    # Build base transfer matrices from regular scattering processes
    transfer_mats = [0.0 * mat for mat in transfer_mats]
    for m in range(0, max_num_moms):
        transfer_mats[m] += transfer_el[m]
        transfer_mats[m] += transfer_inel[m]
        transfer_mats[m] += transfer_nxn[m]
        transfer_mats[m] += transfer_n_other[m]

    # We already have the freegas thermal scattering matrix
    # transfer_freegas

    # Build SAB matrix (if necessary)
    if with_sab:
        for m in range(0, max_num_moms):
            transfer_sab[m] = transfer_sab_el[m] + \
                              transfer_sab_inel[m]

    # Apply freegas thermal scattering
    for m in range(0, max_num_moms):
        mat = transfer_freegas[m]
        for gprime in range(0, G):
            for g in range(0, G):
                if np.abs(mat[gprime, g]) > 1.0e-18:
                    transfer_mats[m][gprime, :] = mat[gprime, :]
                    break  # from for g

    # Apply SAB thermal treatment
    if with_sab:
        for m in range(0, max_num_moms):
            mat = transfer_sab[m]
            for gprime in range(0, G):
                for g in range(0, G):
                    if np.abs(mat[gprime, g]) > 1.0e-18:
                        transfer_mats[m][gprime, :] = mat[gprime, :]
                        break  # from for g

    # Compute new sigma_t, sigma_a, sigma_s 
    # FIXME: Ask Dr. Ragusa why ?
    sig_s = np.sum(transfer_mats[0], axis = 1)
    #Correction for gamma scattering
    if G_g > 0:
        if (len(sig_coht) > 0) and (len(sig_incoht) > 0):
            for i in range (len(sig_coht)):
                sig_s[G_n + i] = sig_coht[i] + sig_incoht[i]

    #sig_t = sig_a + sig_s
    ### This is to correct the total due to across-particule-type transfers
    # aat this stage, sig_t, sig_el, and sig_inel only contains within-particle interactions
    # thus sig_a is only the within-particle absorption
    sig_a = sig_t - sig_el - sig_inel
    
    # Correction for gamma absorption
    if G_g > 0:
        sig_pairprod = np.zeros(G_g)
        if "(g,pair_production)" in cross_sections:
            data = cross_sections["(g,pair_production)"]
            for entry in data:
                g = entry[0]
                v = entry[1]
                sig_pairprod[G_g - g - 1] += v
        
        for i in range (len(sig_pairprod)):
            sig_a[G_n + i] = sig_t[G_n + i] - sig_s[G_n + i] - sig_pairprod[i]
            
    #
    sig_sab = sig_el_sab + sig_inel_sab

    #for range_data in nranges_to_nranges_other:
    #    AddTransferNeutron(range_data)

    # Regular n,\gamma transfer <- add to above
    for range_data in nranges_to_granges:
        AddTransferNeutron(range_data, offset = G_g)

    # Regular \gamma,\gamma transfer <- add to above
    for range_data in granges_to_granges:
        AddTransferGamma(range_data)
        
    # ===== Print outs
    if verbose:
        diff_el = np.sum(sig_el - np.sum(transfer_el[0], axis = 1))
        print('\nElastic:\n\t', diff_el)
        diff_inel = np.sum(sig_inel - np.sum(transfer_inel[0], axis = 1))
        print('\nInelastic:\n\t', diff_inel)
        if with_sab:
            diff_sab_el = np.sum(sig_el_sab - np.sum(transfer_sab_el[0], axis = 1))
            print('\nElastic S(alpha, beta):\n\t', diff_sab_el)
            diff_sab_inel = np.sum(sig_inel_sab - np.sum(transfer_sab_inel[0], axis = 1))
            print('\nInelastic S(alpha, beta):\n\t', diff_sab_inel)
        else:
            diff_freegas = np.sum(sig_freegas - np.sum(transfer_freegas[0], axis = 1))
            print('\nFreegas:\n\t', diff_freegas)
        print('\n')

    # ===== Determine sparsity of the transfer matrices
    transfer_mats_nonzeros = []
    for m in range(0, max_num_moms):
        mat_non_zeros = []
        for gprime in range(0, G):
            non_zeros = []
            for g in range(0, G):
                if np.abs(transfer_mats[m][gprime, g]) > 1.0e-16:
                    non_zeros.append(g)
            mat_non_zeros.append(non_zeros)
        transfer_mats_nonzeros.append(mat_non_zeros)

    # ===== Plot the matrix
    if plot:
        Atest = np.copy(transfer_mats[0])
        nz = np.nonzero(Atest)
        Atest[nz] = np.log10(Atest[nz]) + 10.0

        plt.figure(figsize = (6, 6))
        im = plt.imshow(Atest, cmap = cm.Greys)
        plt.xticks(np.arange(0, G, 10), [str(g) for g in range(0, G, 10)])
        plt.yticks(np.arange(0, G, 10), [str(g) for g in range(0, G, 10)])
        plt.xlabel('Destination energy group')
        plt.ylabel('Source energy group')
        plt.gca().xaxis.set_ticks_position('top')
        plt.gca().xaxis.set_label_position('top')
        plt.savefig("TransferMatrix_NJOY.png")

        if neutn_gs != []:
            # ================================== Build group structures
            np_neutn_gs = np.matrix(neutn_gs)
            nbin_lo = np_neutn_gs[:, 1]
            nbin_hi = np_neutn_gs[:, 2]
            nbin_center = (0.5 * (nbin_lo + nbin_hi))[::-1]

            # ================================== Plot Cross sections
            fig, ax = plt.subplots(ncols = 2, figsize = (6, 6))
            ax[0].semilogx(nbin_center, sig_t[:G_n], label = r"$\sigma_t$")
            ax[0].semilogx(nbin_center, sig_a[:G_n], label = r"$\sigma_a$")
            ax[0].semilogx(nbin_center, sig_s[:G_n], label = r"$\sigma_s$")
            ax[0].legend()
            ax[0].grid(True)

            # ================================== Plot scattering
            ax[1].semilogx(nbin_center, sig_s[:G_n], label = r"$\sigma_s$")
            ax[1].semilogx(nbin_center, sig_el[:G_n], label = r"$\sigma_s$ elastic")
            ax[1].semilogx(nbin_center, sig_inel[:G_n], label = r"$\sigma_s$ inelastic")
            if not with_sab:
                ax[1].semilogx(nbin_center, sig_freegas[:G_n], label = r"$\sigma_s$ freegas")
            else:
                ax[1].semilogx(nbin_center, sig_sab[:G_n], label = r"$\sigma_s$ total SAB")
                ax[1].semilogx(nbin_center, sig_el_sab[:G_n], label = r"$\sigma_s$ elastic SAB")
                ax[1].semilogx(nbin_center, sig_inel_sab[:G_n], label = r"$\sigma_s$ inelastic SAB")
            ax[1].legend()
            ax[1].grid(True)

            plt.show()
            plt.savefig("CrossSections_n_NJOY.png")

        if gamma_gs != []:
            # ================================== Build group structures
            np_gamma_gs = np.matrix(gamma_gs)
            nbin_lo = np_gamma_gs[:, 1]
            nbin_hi = np_gamma_gs[:, 2]
            nbin_center = (0.5 * (nbin_lo + nbin_hi))[::-1]

            # ================================== Plot Cross sections
            plt.figure(figsize = (6, 6))
            plt.semilogx(nbin_center, sig_t[G_n:], label = r"$\sigma_t$")
            plt.semilogx(nbin_center, sig_a[G_n:], label = r"$\sigma_a$")
            plt.semilogx(nbin_center, sig_s[G_n:], label = r"$\sigma_s$")
            plt.legend()
            plt.grid(True)

            plt.show()
            plt.savefig("CrossSections_g_NJOY.png")

    # ================================== Build return data
    return_data = {"neutron_gs": neutn_gs,
                   "gamma_gs": gamma_gs,
                   "sigma_t": sig_t,
                   "sigma_a": sig_a,
                   "sigma_s": sig_s,
                   "sigma_s_el": sig_el,
                   "sigma_s_inel": sig_inel,
                   "sigma_s_freegas": sig_freegas,
                   "sigma_s_sab": sig_sab,
                   "sigma_s_sab_el": sig_el_sab,
                   "sigma_s_sab_inel": sig_inel_sab,
                   "sigma_f": sig_f,
                   "sigma_heat": sig_heat,
                   "nu_total": nu_total,
                   "nu_prompt": nu_prompt,
                   "nu_delayed": nu_delayed,
                   "chi_prompt": chi_prompt,
                   "chi_delayed": chi_delayed,
                   "decay_constants": decay_const,
                   "gamma": gamma,
                   "inv_velocity": inv_v,
                   "transfer_matrices": transfer_mats,
                   "transfer_matrices_sparsity": transfer_mats_nonzeros}
    return return_data, problem_description
