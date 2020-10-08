'''Prepares input/output files for NJOY2016 and executes NJOY2016'''
import sys 
import argparse
import os

#====================================== Check version
if sys.version_info[0] < 3:
    print("\n Error: This script requires python3 but was executed "
          "with version:\n\n"+sys.version+"\n")
    sys.exit(1)

#====================================== Setup arguments
ign = [*range(1,35)]
igg = [*range(0,11)]
niwt = [*range(1,13)]
giwt = [*range(1,4)]
argparser = argparse.ArgumentParser(
    description="Prepares input/output files "
                "for NJOY2016 and executes NJOY2016")

argparser.add_argument("--path_to_neutron_endf",
                       help="Path to the isotope's endf file for "
                       "incident neutrons",
                       default="",required=True)
argparser.add_argument("--temperature",
                       help="NJOY temperature output",type=float,
                       default=296.0)                       
argparser.add_argument("--output_directory",
                       help="Directory where output file needs to be stored."
                       " If not supplied, will default to current directory.",
                       default="")
argparser.add_argument("--output_filename",
                       help="Output file name. If not supplied then a default "
                       "file name will be constructed.",
                       default="")                                              
argparser.add_argument("--path_to_gamma_endf",
                       help="Path to the isotope's endf file for "
                       "incident gammas",
                       default="")
argparser.add_argument("--path_to_sab",
                       help="Path to the isotope's thermal scattering data",
                       default="")
argparser.add_argument("--neutron_group_structure",
                       help="ign in the NJOY2016 manual. 1 for user specified "
                       "(requires --custom_neutron_gs_file to be set)",
                       choices=ign,
                       default=9)
argparser.add_argument("--gamma_group_structure",
                       help="igg in the NJOY2016 manual. 1 for user specified "
                       "(requires --custom_gamma_gs_file to be set)",
                       choices=igg,
                       default=0)
argparser.add_argument("--neutron_weight_function",
                       help="Neutron weight function. 1 for user specified "
                       "(requires --custom_neutron_wt_file to be set)",
                       choices=niwt,
                       default=8)
argparser.add_argument("--gamma_weight_function",
                       help="Gamma weight function. 1 for user specified "
                       "(requires --custom_gamma_wt_file to be set)",
                       choices=giwt,
                       default=2)    

argparser.add_argument("--custom_neutron_gs_file",
                       help="Custom neutron group structure file path.",
                       default="")    
argparser.add_argument("--custom_gamma_gs_file",
                       help="Custom gamma group structure file path.",
                       default="")   
argparser.add_argument("--custom_neutron_wt_file",
                       help="Custom neutron weight function file path.",
                       default="")   
argparser.add_argument("--custom_gamma_wt_file",
                       help="Custom gamma weight function file path.",
                       default="")                                                                                                                                     

args = argparser.parse_args()

#====================================== Check argument rules
if (args.neutron_group_structure == 1 and 
        args.custom_neutron_gs_file == ""):
    print("Error: When specifying \"1\" for neutron group structure "
          "then argument \"custom_neutron_gs_file\" must be supplied.")

with_gamma = False
if (args.path_to_gamma_endf != ""):
    with_gamma = True

if (with_gamma and args.gamma_group_structure == 1 and 
        args.custom_gamma_gs_file == ""):
    print("Error: When specifying \"1\" for gamma group structure "
          "then argument \"custom_gamma_gs_file\" must be supplied.")

if (args.neutron_weight_function == 1 and 
        args.custom_neutron_wt_file == ""):
    print("Error: When specifying \"1\" for neutron weight function "
          "then argument \"custom_neutron_wt_file\" must be supplied.")

if (with_gamma and args.gamma_weight_function == 1 and 
        args.custom_gamma_wt_file == ""):
    print("Error: When specifying \"1\" for gamma weight function "
          "then argument \"custom_gamma_wt_file\" must be supplied.")

os.system("ln -fs "+args.path_to_neutron_endf+" tape20")

#====================================== Extract isotope information from
#                                       endf neutron file
endf_file = open("tape20","r")
endf_file.readline() #skip line
endf_file.readline() #skip line
endf_file.readline() #skip line
endf_file.readline() #skip line
endf_file.readline() #skip line

isotope_text = endf_file.readline(11); endf_file.readline()
isotope_atomic_number = int(isotope_text[0:3])
isotope_symbol        = isotope_text[4:6].strip()
isotope_mass_number   = int(isotope_text[7:10])
isotope_metastable    = isotope_text[10:].strip()

endf_file.readline()   #skip line
endf_file.readline(11) #skip field
endf_file.readline(11) #skip field
material_number = int(endf_file.readline(22)[9:])

output_filename = args.output_filename
if (args.output_filename == ""):
    output_filename = isotope_symbol + str(isotope_mass_number) + \
                      isotope_metastable+".txt"

endf_file.close()

#====================================== Extract information from S(a,b) file
with_sab = False 
material_number_sab = 0
if (args.path_to_sab != ""):
    with_sab = True

    sab_file = open(args.path_to_sab,"r")
    os.system("ln -fs "+args.path_to_sab+" tape51")

    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline() #skip line
    sab_file.readline(11) #skip field
    sab_file.readline(11) #skip field
    material_number_sab = int(sab_file.readline(22)[9:])

    print("S(a,b) material number read: "+str(material_number_sab))

    sab_file.close()

#====================================== Begin writing input
njoy_input = open("NJOY_INPUT.txt","w")

njoy_input.write("-- Processing ENDF to PENDF\n")
njoy_input.write("moder\n")
njoy_input.write("20 -21/\n")

njoy_input.write("reconr\n")
njoy_input.write("-21 -22/\n")
njoy_input.write("'pendf tape'/\n")
njoy_input.write(str(material_number)+" 0/\n")
njoy_input.write("0.001/\n")
njoy_input.write("0/\n")

njoy_input.write("broadr\n")
njoy_input.write("-21 -22 -23\n")
njoy_input.write(str(material_number)+" 1 0 0 0/\n")
njoy_input.write("0.001/\n")
njoy_input.write(str(args.temperature)+"/\n")
njoy_input.write("0/\n")

njoy_input.write("unresr\n")
njoy_input.write("-21 -23 -24/\n")
njoy_input.write(str(material_number)+" 1 1 0/\n")
njoy_input.write(str(args.temperature)+"/\n")
njoy_input.write("0.0/\n")
njoy_input.write("0/\n")

njoy_input.write("heatr\n")
njoy_input.write("-21 -24 -90/\n")
njoy_input.write(str(material_number)+" 0/\n")

#====================================== Thermr free-gas
njoy_input.write("thermr\n")
njoy_input.write("0 -90 -25/\n")
njoy_input.write("0 "+str(material_number)+" 16 1 ")
njoy_input.write("1 ")  #inelastic options 
                        # 0=none
                        # 1=compute as free gas
                        # 2=read s(a,b) and compute matrix
njoy_input.write("0 ")     #elastic options
                           # 0=none
                           # 1=compute using ENDF6 format data (if supplied)
njoy_input.write("0 1 221 1/\n")    
njoy_input.write(str(args.temperature)+"/\n") 
njoy_input.write("0.05 100.0/\n")   

#====================================== Thermr S(a,b)
if (with_sab):
    njoy_input.write("thermr\n")
    njoy_input.write("51 -25 -26/\n")
    njoy_input.write(str(material_number_sab)+" "+str(material_number)+" 16 1 ")
    njoy_input.write("2 ")  #inelastic options 
                            # 0=none
                            # 1=compute as free gas
                            # 2=read s(a,b) and compute matrix
    njoy_input.write("1 ")     #elastic options
                            # 0=none
                            # 1=compute using ENDF6 format data (if supplied)
    njoy_input.write("0 1 249 1/\n")    
    njoy_input.write(str(args.temperature)+"/\n") 
    njoy_input.write("0.05 100.0/\n") 

#====================================== Groupr
njoy_input.write("groupr\n")
if (with_sab):
    njoy_input.write("-21 -26 0 -30\n")
else:
    njoy_input.write("-21 -25 0 -30\n")
njoy_input.write(str(material_number)+" ")
njoy_input.write(str(args.neutron_group_structure)+" ")
njoy_input.write(str(args.gamma_group_structure)+" ")
njoy_input.write(str(args.neutron_weight_function)+" ")
njoy_input.write("7 ")              # Legendre order
njoy_input.write("1 1 1 1/\n")
njoy_input.write(isotope_symbol+str(isotope_mass_number)+isotope_metastable+"/\n")
njoy_input.write(str(args.temperature)+"/\n") 
njoy_input.write("0.0/\n") 

njoy_input.write("3/\n") #All MF3 1D cross-sections
njoy_input.write("3 259 'inverse velocity'/\n") 
njoy_input.write("3 221 'free-gas therm scat'/\n") 

njoy_input.write("6/\n") #All MF6 transfer matrices
njoy_input.write("6 221 'free-gas transfer matrix'/\n")

if (with_sab):
    njoy_input.write("6 249 'inelasting s(a,b) matrix'/\n")
    njoy_input.write("6 250 'elasting s(a,b) matrix'/\n")

njoy_input.write("0/\n") #Terminate reaction list

njoy_input.write("0/\n") #Terminate groupr

njoy_input.write("moder\n")
njoy_input.write("-30 31/\n")

njoy_input.write("stop\n")
njoy_input.close()


#====================================== Run NJOY
os.system("njoy < NJOY_INPUT.txt")
os.system("rm tape*")

print("Copying outputfile to "+args.output_directory+output_filename)

os.system("cp output "+args.output_directory+output_filename)