"""Writes cross-sections to file using the combined Chi-tech data"""
import numpy as np

def WriteChiTechFile(data, chi_filename, comment = "# Output"):
    #Start writing the file
    num_group = data['G']
    num_moment = data['M']
    
    cf = open(chi_filename, 'w')
    cf.write(comment + "\n")
    cf.write("NUM_GROUPS " + str(num_group) + "\n")
    cf.write("NUM_MOMENTS " + str(num_moment) + "\n")
    cf.write("\n")
    
    
    for key in data:
        #print(key)
        #print(bool(((key != "G") and (key != "M"))))
        if ((key != "G") and (key != "M")):
            if (key == 'neutron_gs'):
                neutron_gs = data[key]
                print(len(neutron_gs))
                if neutron_gs != [] :
                    cf.write("NEUTRON_GS_BEGIN" + "\n")
                    for n in range(0, len(neutron_gs)):
                        cf.write("{:<4d}".format(n) + " ")
                        cf.write("{:<8g}".format(neutron_gs[n][1]) + " ")
                        cf.write("{:<g}".format(neutron_gs[n][2]) + " ")
                        cf.write("\n")
                    cf.write("NEUTRON_GS_END" + "\n\n")
            
            elif (key == 'gamma_gs'):
                gamma_gs = data[key]
                print(len(gamma_gs))
                if gamma_gs != []:
                    cf.write("GAMMA_GS_BEGIN" + "\n")
                    for g in range(0, len(gamma_gs)):
                        cf.write("{:<4d}".format(g) + " ")
                        cf.write("{:<8g}".format(gamma_gs[g][1]) + " ")
                        cf.write("{:<g}".format(gamma_gs[g][2]) + " ")
                        cf.write("\n")
                    cf.write("GAMMA_GS_END" + "\n\n")
                
            elif (key == 'transfer_matrices'):
                # ============= Start writing the transfer matrix
                transfer_mats = data[key]
                transfer_mats_nonzeros = data["transfer_matrices_sparsity"]
                cf.write("TRANSFER_MOMENTS_BEGIN" + "\n")
                for m in range (num_moment):
                    cf.write("# l = " + str(m) + "\n")
                    for gprime in range (num_group):
                        for g in transfer_mats_nonzeros[m][gprime]:
                            cf.write("M_GPRIME_G_VAL" + " ")
                            cf.write("{:<4d}".format(m) + " ")
                            cf.write("{:<4d}".format(gprime) + " ")
                            cf.write("{:<4d}".format(g) + " ")
                            cf.write("{:<g}".format(transfer_mats[m][gprime][g]))
                            cf.write("\n")
                cf.write("TRANSFER_MOMENTS_END" + "\n\n")    
                
            elif (key != "transfer_matrices_sparsity"):
                #================  Start writing the cross_sections
                cf.write(key.upper()+ "_BEGIN" + "\n")
                for i in range(len(data[key])):
                    cf.write("{:<4d}".format(i) + " ")
                    cf.write("{:<g}".format(data[key][i]))
                    cf.write("\n")
                cf.write(key.upper() + "_END" + "\n\n")
        

    cf.close()