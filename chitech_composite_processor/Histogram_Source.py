import numpy as np
import dtra_n_g_mcnp_post_process as mcnp_reader
import matplotlib.pyplot as plt
def Source_list(mcnp_histo_source, source_gs):
    # ================== Get the MCNP tally values
    tally_list = ['2']
    tally, tally_dict, atom_density, gram_density = mcnp_reader.ReadMcnpFile(mcnp_histo_source, tally_list)
    #Check if there's only one tally from MCNP
    tally_bins, tally_value = [], []
    if len(tally) == 1:
        for value in tally:
            tally_bins = value['values'][:,0]
            tally_value = value['values'][:,1]

    print(len(tally_bins), len(tally_value))
    tally_bins = np.insert(tally_bins,0, 0.001)
    tally_sum  = np.sum(tally_value)
    tally_value = [x/tally_sum for x in tally_value]
    print(len(tally_bins), len(tally_value))
    print(np.sum(tally_value))
    source_term = np.zeros(len(source_gs) -1)
    for i in range (len(source_gs) - 1):
        upper_bound = 0
        lower_bound = 0
        for j in range (len(tally_bins)):
            # Get the upper and lower bounds first
            if (source_gs[i] >= tally_bins[j]):
                lower_bound = j
                #source_prob += tally_value[j]
            elif (source_gs[i+1] <= tally_bins[j]) and (source_gs[i+1] > tally_bins[j-1]):
                upper_bound = j
        #print("For group %.d, from %.3f to %.3f, source bins from index %d to %d"%(i,source_gs[i], source_gs[i+1], lower_bound, upper_bound))
        
        if upper_bound != 0:
        #Start combining the tally value into the source term
            source_prob = 0
            for k in range (lower_bound, upper_bound):
                # If the source bin is inside a tally bin:
                if (upper_bound - lower_bound) == 1:
                    dE_source = source_gs[i+1] - source_gs[i]
                    dE_tally = tally_bins[k+1] - tally_bins[k]
                    ratio = dE_source / dE_tally
                    source_prob += tally_value[k] * ratio
                else:
                    #Case 1: Overlapping lower bounds
                    if tally_bins[k] < source_gs[i]:
                        dE_source = tally_bins[k+1] - source_gs[i]
                        dE_tally = tally_bins[k+1] - tally_bins[k]
                        ratio = dE_source / dE_tally
                        source_prob += tally_value[k] * ratio 
                    #Case 2: Overlapping upper bounds
                    elif tally_bins[k+1] > source_gs[i+1]:
                        dE_source = source_gs[i+1] - tally_bins[k]
                        dE_tally = tally_bins[k+1] - tally_bins[k]
                        ratio = dE_source / dE_tally
                        source_prob += tally_value[k] * ratio
                    else:
                        source_prob += tally_value[k]
            source_term[i] = source_prob

    source_term = np.flip(source_term)
    #Write it in a text file
    cf = open('48g_Gamma_HistogramS.txt','w')
    cf.write("Group | Value" + "\n")
    for n in range(0, len(source_term)):
        cf.write("{:<7d}".format(n+1) + " ")
        cf.write("{:<g}".format(source_term[n]) + " ")
        cf.write("\n")
    cf.close()

    print("Sum of 100g HistoS = %.4f"%(np.sum(tally_value)))
    print("Sum of 48g HistoS = %.4f"%(np.sum(source_term)))
    #Plot the 100g source vs the 48g collapsed source
    plt.figure(dpi = 750)
    
    #Plot the 100 Histo Source
    bndrys, vals = [], []
    for g in range (0, len(tally_bins) - 1):
        low_bound = tally_bins[g]
        high_bound = tally_bins[g+1] 
        bndrys += [low_bound, high_bound]
        vals += [tally_value[g], tally_value[g]]

    bndrys = np.array(bndrys)
    vals = np.array(vals)
    plt.plot(bndrys, vals, color = 'r', linewidth = 2.0, label='Uncollapsed 100g HistoS')

    #Plot the 100 Histo Source
    bndrys, vals = [], []
    source_term = np.flip(source_term)
    for g in range (0, len(source_gs) - 1):
        low_bound = source_gs[g]
        high_bound = source_gs[g+1] 
        bndrys += [low_bound, high_bound]
        vals += [source_term[g], source_term[g]]
    bndrys = np.array(bndrys)
    vals = np.array(vals)
    plt.plot(bndrys, vals, color = 'm', linewidth = 1.0, label='Collapsed 48g HistoS')
    plt.xlim([0,6])
    plt.grid()
    plt.xlabel('Energy [MeV]')
    plt.ylabel('Source Fraction')
    plt.savefig('100g_HistoS_vs_Collapsed_48g_HistoS.png')
    source_term = np.flip(source_term)


    return source_term
    