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
    nL = nL_i + 5
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
    nelastic_data = raw_njoy_data["cross_sections"]["(n,elastic)"]
    ninelstc_data = raw_njoy_data["cross_sections"]["(n,inelastic)"]
    n_to_n_data   = raw_njoy_data["cross_sections"]["(n,2n)"]

    coherent_data = raw_njoy_data["cross_sections"]["(g,coherent)"]
    incohrnt_data = raw_njoy_data["cross_sections"]["(g,incoherent)"]
    pp_xs_data    = raw_njoy_data["cross_sections"]["(g,pair_production)"]

    sig_s = np.zeros(G)

    for entry in nelastic_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n-g-1] += v

    for entry in ninelstc_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n-g-1] += v   

    for entry in n_to_n_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n-g-1] += v   

    for entry in coherent_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n + G_g-g-1] += v 

    for entry in incohrnt_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n + G_g-g-1] += v 

    for entry in pp_xs_data:
        g = entry[0]
        v = entry[1]
        sig_s[G_n + G_g-g-1] += v     

    # ================================= Combine multiplication data
    sig_f = np.zeros(G)

    nu = np.zeros(G)

    chi = np.zeros(G)

    # ================================= Combine transfer matrices
    nelastic_data = raw_njoy_data["transfer_matrices"]["(n,elastic)"]
    ninelstc_data = raw_njoy_data["transfer_matrices"]["(n,inel)"]

    n_2n_mat_data = raw_njoy_data["transfer_matrices"]["(n,2n)"]

    ngamma_data = raw_njoy_data["transfer_matrices"]["(n,g)"]

    coherent_data = raw_njoy_data["transfer_matrices"]["(g,coherent)"]
    incohrnt_data = raw_njoy_data["transfer_matrices"]["(g,incoherent)"]
    pp_trnfr_data = raw_njoy_data["transfer_matrices"]["pair_production"]

    #Determine number of moments
    num_moms = len(nelastic_data[0])-2

    transfer_mats = []
    for m in range(0,num_moms):
        transfer_mats.append(np.zeros((G,G)))

    transfer_mats_nonzeros = []
    
    # (n,elastic)
    for entry in nelastic_data:
        gprime = G_n - entry[0] - 1
        g      = G_n - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    # (n,inelastic)
    for entry in ninelstc_data:
        gprime = G_n - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    # (n,2n)
    for entry in n_2n_mat_data:
        gprime = G_n - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    # (n,gamma)
    for entry in ngamma_data:
        gprime = G_n - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    # (g,coherent)
    for entry in coherent_data:
        gprime = G_n + G_g - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v
    
    # (g,incoherent)
    for entry in incohrnt_data:
        gprime = G_n + G_g - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    # (pair_production)
    for entry in pp_trnfr_data:
        gprime = G_n + G_g - entry[0] - 1
        g      = G_n + G_g - entry[1] - 1

        for m in range(0,num_moms):
            v = entry[m+2]
            transfer_mats[m][gprime,g] += v

    for m in range(0,num_moms):
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


#####################################################################
# Stand-alone usage
if __name__ == "__main__":
    raw_njoy_data = ReadNJOYfile()

    data = BuildCombinedData(raw_njoy_data)
    WriteChiTechFile(data)

    # print(raw_njoy_data["transfer_matrices"]["(n,elastic)"][0])
    # print(raw_njoy_data["transfer_matrices"]["coherent"][11])
    # print(raw_njoy_data["transfer_matrices"]["incoherent"][11])
    # print(raw_njoy_data["transfer_matrices"]["pair_production"][8])
