"""
@author: ragusa
"""
import numpy as np
import matplotlib.pyplot as plt
import copy

def convert_spectrum(spectrum):
    energy_bins = spectrum[1]
    values = spectrum[0][1:]
    bndrys = []
    vals = []
    spectrum = []
    dE = np.diff(energy_bins)
    for g in range (0, len(energy_bins) - 1):
        low_bound = energy_bins[g]
        high_bound = energy_bins[g+1]
        width = dE[g]
        spectrum_value = (values[g]/width) 
        
        bndrys += [low_bound, high_bound]
        vals += [spectrum_value, spectrum_value]


    bndrys = np.array(bndrys)
    vals = np.array(vals)
    return bndrys, vals

def ReadMcnpFile(output_filename):
    #tally_lists = ['4', '24', '6', '26']
    tally_lists = ['4', '24', '6', '26', '104', '204', '106', '206']

    with open(output_filename, "r") as fileHandler:
        tally = [] 
        for line in fileHandler:
            # print every single line
            # print(line.strip())
            #==================== Get the atom and gram density
            if line.count('density') > 0:
                next(fileHandler)
                next_line = fileHandler.readline().strip()
                number = next_line.split()
                atom_density = float(number[3])
                gram_density = float(number[4])
                print(atom_density, gram_density)
            # found a tally
            elif line.count('1tally')>0:
                # split that line and check tally ID
                my_list = line.split()
                
                # tally ID match
                tally_str = []
                tally_str = [match for match in tally_lists if my_list[1] == match]
                
                if tally_str:
                    print('Processing tally ',tally_str)
                    # skip some many lines before the good part
                    for _ in range(9):
                        next(fileHandler)
                    # read until blank line
                    tally2 = []
                    read_block = lambda: fileHandler.readline().strip()
                    for block in iter(read_block, ''):
                        tally2.append(block)
                
                    # remove last line of tally2 (total)
                    if tally2[-1][0:5]=='total':
                        tally2.pop()
                    # convert each entry (str) to list of floats
                    aux =[[float(v) for v in x.split()] for x in tally2]
                    # convert to np array and remove 3rd column
                    t2 = np.asarray(aux)
                    t2 = t2[:,0:2]
                    # check if user put a 0 as first entry of the En card
                    check = np.isclose(t2[0,0]/10,0, atol=1e-13)
                    if check:
                        t2 = t2[1:,:]

                    my_dict = { 'tallyID': tally_str[0], 'values': t2}
                    tally.append(copy.deepcopy(my_dict))

    #Putting tallies into catagories:
    tally_dict = {}
    neutron_list = []
    gamma_list = []
    for tally_name in tally_lists:
        # If deals with neutrons        
        if tally_name.count('2') == 0:
            neutron_list.append(tally_name)
        elif tally_name.count('2') > 0:
            gamma_list.append(tally_name)
    tally_dict["Neutron"] = neutron_list
    tally_dict["Gamma"] = gamma_list
    for particle_type in tally_dict.keys():
        flux = []
        heating = []
        for tally_name in tally_dict[particle_type]:
            if tally_name.count('4') > 0:
                flux.append(tally_name)
            elif tally_name.count('6') > 0:
                heating.append(tally_name)
        tally_dict[particle_type] = {}
        tally_dict[particle_type]["Flux"] = flux
        tally_dict[particle_type]["Heating"] = heating

    #Start storing output in .txt file
    MCNP_Writer(tally, atom_density, gram_density, output_filename)

    return tally, tally_dict, atom_density, gram_density

def MCNP_Writer (tally, atom_density, gram_density, output_filename):
    string = output_filename.split('.out')
    new_output = string[0] + '.txt'
    print("Begin Storing the MCNP output at" + new_output)
    
    cf = open(new_output, 'w')
    cf.write("# OUTPUT" + "\n")
    cf.write("ATOM_DENSITY = " + str(atom_density) + "\n")
    cf.write("MASS DENSITY = " + str(gram_density) + "\n")
    cf.write("\n")

    for tally_value in tally:
        E = tally_value['values'][:,0]
        Flx = tally_value['values'][:,1]
        tally_ID = tally_value['tallyID']
        cf.write("TALLY_" + tally_ID +"_BEGIN" + "\n")
        cf.write("|NUM |E" + 9*" " +"|Value" +"\n")
        for i in range(0, len(Flx)):
            cf.write("{:<4d}".format(i) + " ")
            cf.write("{:<10g}".format(E[i]) + " ")
            cf.write("{:<g}".format(Flx[i]) + " ")
            cf.write("\n")
        cf.write("TALLY_" + tally_ID +"_END" + "\n\n")
    
    cf.close()


    

