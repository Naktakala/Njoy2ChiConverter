''' Converts neutron gamma njoy output to chitech-format '''
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#====================================================================
def PrintLine(line_num,line):
    ''' Prints a file line without the line feed. '''
    print("Line: ",line_num,line,end='') 

#====================================================================
def ProcessGroupStructure(nL_i,lines):
    ''' Processes a group structure.
        params:
            nL_i  File line number before data starts.
            lines An array of files lines.
            
        return:
            returns the group_structure '''
    group_struct = []
    nL = nL_i + 1
    while (True):
        words = lines[nL].split()
        if not words: 
            break 

        group_struct.append( [int(words[0])-1,
                              float(words[1]),
                              float(words[3])])
        nL += 1
    
    return group_struct

#====================================================================
def ProcessCrossSection(nL_i, lines,header_size=5,line_incr=1):
    ''' Reads 1D neutron cross-sections from the lines. 
        params:
            nL_i    File line number before data starts.
            lines   An array of file lines.
        
        return:
            A table containing the group wise xs. '''
    xs=[]
    # Skip header
    nL = nL_i + header_size
    words = lines[nL].split()
    num_words = len(words)

    while (num_words>=2):
        number_word = words[1]
        number_word = number_word.replace("+","E+")
        number_word = number_word.replace("-","E-")
        xs.append([int(words[0])-1,float(number_word)])
        nL += line_incr

        words = lines[nL].split()
        num_words = len(words)

    return xs

#====================================================================
def ProcessTransferMatrix(nL_i, lines):    
    ''' Reads transfer matrices from the lines. 
        params:
            nL_i    File line number before data starts.
            lines   An array of file lines.
        
        return:
            A table containing the group wise transfer coefficients. '''    
    matrix = []
    # Skip 4 lines
    add_skip = 0
    if lines[nL_i+1].find("particle emission")>=0:
        add_skip = 1
    nL = nL_i + 5 + add_skip
    words = lines[nL].split()
    num_words = len(words)

    while (num_words>=2):
        # Transform words 0.000-00 to 0.000E+00
        for i in range(0,num_words):
            word = words[i]
            loc_of_sign = word.find("-",2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+",2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+","E+")
            secnd_part = secnd_part.replace("-","E-")

            words[i] = first_part+secnd_part

        # Develop table entry
        entry = []
        entry.append(int(words[0])-1)
        entry.append(int(words[1])-1)
        for i in range(2,num_words):
            entry.append(float(words[i]))

        matrix.append(entry)
        nL += 1

        words = lines[nL].split()
        num_words = len(words)

    return matrix  

#====================================================================
def ProcessTransferMatrixB(nL_i, lines, process_overflow=False):    
    ''' Reads transfer matrices from the lines. 
        params:
            nL_i    File line number before data starts.
            lines   An array of file lines.
            line_incr Line increment to use
        
        return:
            A table containing the group wise transfer coefficients. '''    
    matrix = []
    # Skip 3 lines
    nL = nL_i + 4
    words = lines[nL].split()
    num_words = len(words)

    while (num_words>=2):
        # Transform words 0.000-00 to 0.000E+00
        for i in range(0,num_words):
            word = words[i]
            loc_of_sign = word.find("-",2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+",2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+","E+")
            secnd_part = secnd_part.replace("-","E-")

            words[i] = first_part+secnd_part

        # Develop table entry
        entry = []
        entry.append(int(words[0])-1)
        entry.append(int(words[1])-1)
        for i in range(2,num_words):
            entry.append(float(words[i]))

        # Process overflow
        if process_overflow:
            words = lines[nL+1].split()

            word = words[0]
            loc_of_sign = word.find("-",2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+",2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+","E+")
            secnd_part = secnd_part.replace("-","E-")

            words[0] = first_part+secnd_part

            entry.append(float(words[0]))
            nL += 1

        # Complete matrix entry
        matrix.append(entry)

        # Move to next line
        nL += 1
        words = lines[nL].split()
        num_words = len(words)

        # Skip lines with xsec and heat
        while (num_words>=2 and 
            not words[1].isnumeric()):
            nL += 1
            words = lines[nL].split()
            num_words = len(words)
            

    return matrix  

#====================================================================
def ReadNJOYfile(njoy_filename="output"):
    ''' Reads an NJOY output file
        params:
            njoy_filename Name of the NJOY file to process. 
            
        return:
            Returns a complex dictionary of raw data. '''

    njoy_file = open(njoy_filename, 'r')
    file_lines = njoy_file.readlines()

    njoy_raw_data = {}

    group_structures = {}
    cross_sections = {}
    transfer_matrices = {}

    flag_run_processed             = False
    flag_gamma_structure_processed = False

    nL=-1
    while (nL<(len(file_lines)-1)):
        nL += 1
        line = file_lines[nL]
        words = file_lines[nL].split()
        num_words = len(words)

        if (line.find("neutron group structure") != -1 and 
                words[3] != "option"):
            PrintLine(nL,line)
            group_structures["neutron"] = ProcessGroupStructure(nL,file_lines)

        if (line.find("gamma group structure") != -1):
            if not flag_gamma_structure_processed:
                PrintLine(nL,line)
                group_structures["gamma"] = ProcessGroupStructure(nL,file_lines)
                flag_gamma_structure_processed = True 

        if (line.find("for mf") != -1 and
                line.find("mt") != -1 ):
            PrintLine(nL,line)  

            if (words[2] == "3" and words[5] == "1" ):
                cross_sections["(n,total)"] = \
                    ProcessCrossSection(nL,file_lines,line_incr=2)  

            if (words[2] == "3" and words[5] == "2" ):
                cross_sections["(n,elastic)"] = \
                    ProcessCrossSection(nL,file_lines)        

            if (words[2] == "3" and words[5] == "4" ):
                cross_sections["(n,inelastic)"] = \
                    ProcessCrossSection(nL,file_lines) 

            if (words[2] == "3" and words[5] == "16" ):
                cross_sections["(n,2n)"] = \
                    ProcessCrossSection(nL,file_lines) 

            if (words[num_words-1] == "matrix"):
                transfer_matrices[words[num_words-3]] = \
                    ProcessTransferMatrix(nL,file_lines)

            if (words[1] == "mf23" and words[3] == "mt501" ):
                cross_sections["(g,total)"] = \
                    ProcessCrossSection(nL,file_lines,header_size=4) 

            if (words[1] == "mf23" and words[3] == "mt502" ):
                cross_sections["(g,coherent)"] = \
                    ProcessCrossSection(nL,file_lines,header_size=4)   

            if (words[1] == "mf23" and words[3] == "mt504" ):
                cross_sections["(g,incoherent)"] = \
                    ProcessCrossSection(nL,file_lines,header_size=4)  

            if (words[1] == "mf23" and words[3] == "mt516" ):
                cross_sections["(g,pair_production)"] = \
                    ProcessCrossSection(nL,file_lines,header_size=4)      

            if (words[1] == "mf26" and words[3] == "mt502"):
                transfer_matrices["(g,coherent)"] = \
                    ProcessTransferMatrixB(nL,file_lines,process_overflow=True)

            if (words[1] == "mf26" and words[3] == "mt504"):
                transfer_matrices["(g,incoherent)"] = \
                    ProcessTransferMatrixB(nL,file_lines,process_overflow=True)

            if (words[1] == "mf26" and words[3] == "mt516"):
                transfer_matrices["pair_production"] = \
                    ProcessTransferMatrixB(nL,file_lines,process_overflow=False)         

    njoy_file.close()

    njoy_raw_data["group_structures"] = group_structures
    njoy_raw_data["cross_sections"] = cross_sections
    njoy_raw_data["transfer_matrices"] = transfer_matrices

    return njoy_raw_data

#====================================================================
def BuildCombinedData(raw_njoy_data):
    ''' Combines enjoy raw data into a dictionary of vectors and matrices '''
    # ================================= Determine # of groups
    neutn_gs = raw_njoy_data["group_structures"]["neutron"]
    gamma_gs = raw_njoy_data["group_structures"]["gamma"]

    G_n = len(neutn_gs)
    G_g = len(gamma_gs)
    G   = G_n + G_g 

    # ================================= Combine sig_t
    sig_t_ndata = raw_njoy_data["cross_sections"]["(n,total)"]
    sig_t_gdata = raw_njoy_data["cross_sections"]["(g,total)"]

    sig_t = np.zeros(G)

    for entry in sig_t_ndata:
        g = entry[0]
        v = entry[1]
        sig_t[G_n-g-1] += v

    for entry in sig_t_gdata:
        g = entry[0]
        v = entry[1]
        sig_t[G_n + G_g-g-1] += v

    # ================================= Combine sig_s
    # This cross-section still needs some work but
    # is ultimately not used by Chi-Tech
    nelastic_data = raw_njoy_data["cross_sections"]["(n,elastic)"]
    ninelstc_data = raw_njoy_data["cross_sections"]["(n,inelastic)"]
    n_to_n_data   = raw_njoy_data["cross_sections"]["(n,2n)"]

    coherent_data = raw_njoy_data["cross_sections"]["(g,coherent)"]
    incohrnt_data = raw_njoy_data["cross_sections"]["(g,incoherent)"]
    pp_xs_data    = raw_njoy_data["cross_sections"]["(g,pair_production)"]

    sig_s = np.zeros(G)

    def AddSigSNeutron(data_vals):
        for entry in data_vals:
            g = entry[0]
            v = entry[1]
            sig_s[G_n-g-1] += v

    def AddSigSGamma(data_vals):
        for entry in coherent_data:
            g = entry[0]
            v = entry[1]
            sig_s[G_n + G_g-g-1] += v 

    AddSigSNeutron(nelastic_data)  
    AddSigSNeutron(ninelstc_data) 
    AddSigSNeutron(n_to_n_data)       

    AddSigSGamma(coherent_data)
    AddSigSGamma(incohrnt_data)
    AddSigSGamma(pp_xs_data)   

    # ================================= Combine multiplication data
    sig_f = np.zeros(G)

    nu = np.zeros(G)

    chi = np.zeros(G)

    # ================================= Combine transfer matrices
    nranges_to_nranges = []
    nranges_to_nranges.append(raw_njoy_data["transfer_matrices"]["(n,elastic)"])
    nranges_to_nranges.append(raw_njoy_data["transfer_matrices"]["(n,2n)"])

    #Adding all the (n,nxx) data, inelastic data
    for nn in range(1,24+1):
        rx_name = "(n,n{:02d})".format(nn)
        if rx_name in raw_njoy_data["transfer_matrices"]:
            mat = raw_njoy_data["transfer_matrices"][rx_name]
            nranges_to_nranges.append(mat)

    nranges_to_granges = []
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,g)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,inel)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,np)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,nd)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,p)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,d)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,t)"])
    nranges_to_granges.append(raw_njoy_data["transfer_matrices"]["(n,a)"])

    granges_to_granges = []
    granges_to_granges.append(raw_njoy_data["transfer_matrices"]["(g,coherent)"])
    granges_to_granges.append(raw_njoy_data["transfer_matrices"]["(g,incoherent)"])
    granges_to_granges.append(raw_njoy_data["transfer_matrices"]["pair_production"])

    max_num_moms = 0
    for range_data in nranges_to_nranges:
        max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in nranges_to_granges:
        max_num_moms = max(max_num_moms,len(range_data[0])-2)
    for range_data in granges_to_granges:
        max_num_moms = max(max_num_moms,len(range_data[0])-2)

    transfer_mats = []
    for m in range(0,max_num_moms):
        transfer_mats.append(np.zeros((G,G)))

    transfer_mats_nonzeros = []

    #=======================================
    # Lambda-ish to add neutron data
    def AddTransferNeutron(data_vals,offset=0):
        # (n,elastic)
        for entry in data_vals:
            gprime = G_n - entry[0] - 1
            g      = G_n - entry[1] - 1 + offset

            # Determine number of moments
            num_moms = len(entry)-2
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

    for range_data in nranges_to_nranges:
        AddTransferNeutron(range_data)
    for range_data in nranges_to_granges:
        AddTransferNeutron(range_data,offset=G_g)
    for range_data in granges_to_granges:
        AddTransferGamma(range_data)

    # Determine sparsity of the transfer matrices
    for m in range(0,max_num_moms):
        mat_non_zeros = []
        for gprime in range(0,G):
            non_zeros = []
            for g in range(0,G):
                if (abs(transfer_mats[m][gprime,g]) > 1.0e-18):
                    non_zeros.append(g)
            mat_non_zeros.append(non_zeros)
        transfer_mats_nonzeros.append(mat_non_zeros)

    # Plot the matrix
    plt.figure(figsize=(6,6))
    Atest = transfer_mats[0]
    plt.imshow(np.log10(Atest)+10, cmap=cm.Greys)
    plt.xlabel('Destination energy group')
    plt.ylabel('Source energy group')
    # plt.savefig("SERPENTTransferMatrix.png")
    plt.show()

    #================================== Build return data
    return_data = {}
    return_data["neutron_gs"] = neutn_gs
    return_data["gamma_gs"] = gamma_gs
    return_data["sigma_t"] = sig_t
    return_data["sigma_s"] = sig_s
    return_data["sigma_f"] = sig_f
    return_data["nu"] = nu
    return_data["chi"] = chi
    return_data["transfer_matrices"] = transfer_mats
    return_data["transfer_matrices_sparsity"] = transfer_mats_nonzeros

    return return_data 

#====================================================================
def WriteChiTechFile(data,chi_filename="output.cxs",comment="# Output"):
    cf = open(chi_filename,'w')

    sig_t = data["sigma_t"]
    sig_s = data["sigma_s"]
    sig_f = data["sigma_f"]
    nu = data["nu"]
    chi = data["chi"]
    transfer_mats = data["transfer_matrices"]
    transfer_mats_nonzeros = data["transfer_matrices_sparsity"]

    G = np.size(sig_t)
    M = len(transfer_mats)

    cf.write(comment+"\n")
    cf.write("NUM_GROUPS "+str(G)+"\n")
    cf.write("NUM_MOMENTS "+str(M)+"\n")

    cf.write("SIGMA_T_BEGIN"+"\n")
    for g in range(0,G):
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(sig_t[g]))
        cf.write("\n")
    cf.write("SIGMA_T_END"+"\n")

    cf.write("SIGMA_S_BEGIN"+"\n")
    for g in range(0,G):
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(sig_s[g]))
        cf.write("\n")
    cf.write("SIGMA_S_END"+"\n")

    cf.write("SIGMA_F_BEGIN"+"\n")
    for g in range(0,G):
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(sig_f[g]))
        cf.write("\n")
    cf.write("SIGMA_F_END"+"\n")

    cf.write("NU_BEGIN"+"\n")
    for g in range(0,G):
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(nu[g]))
        cf.write("\n")
    cf.write("NU_END"+"\n")

    cf.write("CHI_BEGIN"+"\n")
    for g in range(0,G):
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(chi[g]))
        cf.write("\n")
    cf.write("CHI_END"+"\n")

    cf.write("TRANSFER_MOMENTS_BEGIN"+"\n")
    for m in range(0,M):
        cf.write("TRANSFER_MOMENT_BEGIN "+str(m)+"\n")
        for gprime in range(0,G):
            cf.write("GPRIME_TO_G " + str(gprime))
            num_limits = G 
            cf.write(" " + str(num_limits)+"\n")

            for g in transfer_mats_nonzeros[m][gprime]:
                cf.write("{:<4d}".format(gprime)+ " ")
                cf.write("{:<4d}".format(g)+ " ")
                cf.write("{:<g}".format(transfer_mats[m][gprime,g]))
                cf.write("\n")
        cf.write("TRANSFER_MOMENT_END"+"\n")
    cf.write("TRANSFER_MOMENTS_END"+"\n")

    cf.close()

#====================================================================
def InfiniteMediumSpectrum(data):
    neutron_gs = data["neutron_gs"]
    gamma_gs = data["gamma_gs"]
    sig_t = data["sigma_t"]
    transfer_mats = data["transfer_matrices"]
    transfer_mats_nonzeros = data["transfer_matrices_sparsity"]

    #======================================= Convert data to numpy data
    # for infinite medium
    G         = np.size(sig_t)
    G_neutron = len(neutron_gs)
    G_gamma   = len(gamma_gs)

    v_sig_t = np.array(sig_t)
    M_sig_gp_to_g = np.zeros([G,G])

    for gprime in range(0,G):
        for g in transfer_mats_nonzeros[0][gprime]:
            M_sig_gp_to_g[gprime,g] = transfer_mats[0][gprime,g]

    S = M_sig_gp_to_g.transpose()
    T = np.diag(v_sig_t)

    A = T - S
    A_inv = np.linalg.inv(A)

    v_src = np.zeros(G)
    v_src[0] = 1.0

    v_psi = np.matmul(A_inv,v_src)

    #======================================= Build data/energy
    neutron_group_bndries = []
    neutron_bndry_edge_values = []
    gamma_group_bndries = []
    gamma_bndry_edge_values = []
    for g in range(0,len(neutron_gs)):
        gprime = G_neutron - g -1

        lo_bound = neutron_gs[g][1]/1e6
        hi_bound = neutron_gs[g][2]/1e6

        bin_width = hi_bound-lo_bound

        neutron_group_bndries.append(lo_bound)
        neutron_group_bndries.append(hi_bound)

        neutron_bndry_edge_values.append(v_psi[gprime]/bin_width)
        neutron_bndry_edge_values.append(v_psi[gprime]/bin_width)

    for g in range(0,len(gamma_gs)):
        gprime = (G_gamma - g -1) + G_neutron

        lo_bound = gamma_gs[g][1]/1e6
        hi_bound = gamma_gs[g][2]/1e6

        bin_width = hi_bound-lo_bound

        gamma_group_bndries.append(lo_bound)
        gamma_group_bndries.append(hi_bound)

        gamma_bndry_edge_values.append(v_psi[gprime]/bin_width)
        gamma_bndry_edge_values.append(v_psi[gprime]/bin_width)

    last_nval = neutron_bndry_edge_values[len(neutron_bndry_edge_values)-1]
    last_gval = gamma_bndry_edge_values[len(gamma_bndry_edge_values)-1]

    for i in range(0,len(neutron_bndry_edge_values)):
        neutron_bndry_edge_values[i] /= last_nval

    for i in range(0,len(gamma_bndry_edge_values)):
        gamma_bndry_edge_values[i] /= last_gval

    plt.plot(neutron_group_bndries, neutron_bndry_edge_values)
    plt.yscale("log")
    plt.show()

    plt.plot(gamma_group_bndries, gamma_bndry_edge_values)
    plt.yscale("log")
    plt.show()


#####################################################################
# Stand-alone usage
if __name__ == "__main__":
    plt.close('all')
    
    raw_njoy_data = ReadNJOYfile()

    data = BuildCombinedData(raw_njoy_data)
    WriteChiTechFile(data)

    InfiniteMediumSpectrum(data)

    # print(raw_njoy_data["transfer_matrices"]["(n,elastic)"][0])
    # print(raw_njoy_data["transfer_matrices"]["coherent"][11])
    # print(raw_njoy_data["transfer_matrices"]["incoherent"][11])
    # print(raw_njoy_data["transfer_matrices"]["pair_production"][8])