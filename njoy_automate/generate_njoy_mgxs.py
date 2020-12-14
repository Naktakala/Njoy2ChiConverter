'''Prepares input/output files for NJOY2016 and executes NJOY2016'''
import sys 
import argparse
import os

#====================================== Check version
if sys.version_info[0] < 3:
    raise ValueError(
      "\n Error: This script requires python3 but was executed "
      "with version:\n\n"+sys.version+"\n"
    )

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
argparser.add_argument("--path_to_sab",
                       help="Path to the isotope's thermal scattering data",
                       default="")
argparser.add_argument("--inelastic_thermal_number",
                       help="MT number to use for incoherent inelastic scattering",
                       default="")
argparser.add_argument("--inelastic_thermal_num_atoms",
                       help="MT number to use for incoherent inelastic scattering",
                       default="1")                       
argparser.add_argument("--elastic_thermal_number",
                       help="MT number to use for coherent/incoherent elastic scattering",
                       default="")
                       
argparser.add_argument("--path_to_gamma_endf",
                       help="Path to the isotope's endf file for "
                       "incident gammas",
                       default="")
argparser.add_argument("--temperature",
                       help="NJOY temperature output",type=float,
                       default=296.0)                       
argparser.add_argument("--neutron_group_structure",
                       help="ign in the NJOY2016 manual. 1 for user specified "
                       "(requires --custom_neutron_gs_file to be set)",
                       choices=ign,
                       type=int,
                       default=9)
argparser.add_argument("--neutron_weight_function",
                       help="Neutron weight function. 1 for user specified "
                       "(requires --custom_neutron_wt_file to be set)",
                       choices=niwt,
                       type=int,
                       default=8)
argparser.add_argument("--output_directory",
                       help="Directory where output file needs to be stored."
                       " If not supplied, will default to current directory.",
                       default="")
argparser.add_argument("--output_filename",
                       help="Output file name. If not supplied then a default "
                       "file name will be constructed.",
                       default="")                                              
argparser.add_argument("--gamma_group_structure",
                       help="igg in the NJOY2016 manual. 1 for user specified "
                       "(requires --custom_gamma_gs_file to be set)",
                       choices=igg,
                       type=int,
                       default=0)
argparser.add_argument("--gamma_weight_function",
                       help="Gamma weight function. 1 for user specified "
                       "(requires --custom_gamma_wt_file to be set)",
                       choices=giwt,
                       type=int,
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
with_gamma = False
if (args.path_to_gamma_endf != ""):
  with_gamma = True

if (args.neutron_group_structure == 1 and 
    args.custom_neutron_gs_file == ""):
  msg =  "Error: When specifying \"1\" for neutron group structure "
  msg += "then argument \"custom_neutron_gs_file\" must be supplied."
  raise ValueError(msg)

if (with_gamma and args.gamma_group_structure == 1 and
    args.custom_gamma_gs_file == ""):
  msg  = "Error: When specifying \"1\" for gamma group structure "
  msg += "then argument \"custom_gamma_gs_file\" must be supplied."
  raise ValueError(msg)

if (args.neutron_weight_function == 1 and \
    args.custom_neutron_wt_file == ""):
  msg  = "Error: When specifying \"1\" for neutron weight function "
  msg += "then argument \"custom_neutron_wt_file\" must be supplied."
  raise ValueError(msg)

if (with_gamma and args.gamma_weight_function == 1 and 
    args.custom_gamma_wt_file == ""):
  msg  = "Error: When specifying \"1\" for gamma weight function "
  msg += "then argument \"custom_gamma_wt_file\" must be supplied."
  raise ValueError(msg)

#====================================== Tape reservations
# The following tape numbers are used:
# neutron endf        20
# neutron sab         51
# gamma endf          52
# moder neutron input 20    output 21
# moder gamma   input 52    output 53
# reconr neutro input 21    output 22
# reconr gamma  input 53    output 54
# broadr        input 21,22 output 23
# unresr        input 21,23 output 24
# heatr         input 21,24 output 90
# thermr freeg  input 90    output 25
# thermr sab    input 51,25 output 26
# groupr        input 21,26 output 30
# gaminr        input 52,54 output 55
# moder print   input 30    output 31

#====================================== Extract info from neutron file
os.system("ln -fs "+args.path_to_neutron_endf+" tape20")
has_fission = False
with open('tape20','r') as endf_file:
  lines = endf_file.readlines()
  for l, line in enumerate(lines):
    line = line.split()
    if l == 4 and "zr" in line[0]:
      atomic_number = int(line[0].split('-')[0])
      atomic_symbol = line[0].split('-')[1]
      isotope = "".join([atomic_symbol,str(atomic_number)])
      mass_number = 91.2
      material_number = int(line[-3])
      isotope_metastable = ""
      break

    if l == 5 and "-" in line[0]:
      atomic_number = int(line[0].split('-')[0])
      atomic_symbol = line[0].split('-')[1]
      isotope = "".join([atomic_symbol,str(atomic_number)])
      mass_number = int(line[2]) if line[1]=='-' else line[1].strip('-')
      isotope_metastable = ""
      material_number = int(lines[l+2].split()[2])
      break
  
  for l, line in enumerate(lines[:200]):
    if "fission" in line:
      has_fission = True
      break
  print("neutron material number read: ", material_number)

#====================================== Extract info from gamma file
material_number_gamma = 0
if (args.path_to_gamma_endf != ""):
    os.system("ln -fs "+args.path_to_gamma_endf+" tape52")
    gamma_file = open("tape52","r")
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline() #skip line
    gamma_file.readline(11) #skip field
    gamma_file.readline(11) #skip field

    material_number_gamma = int(gamma_file.readline(22)[9:])

    print("gamma material number read: "+str(material_number_gamma))

    gamma_file.close()

#====================================== Extract info from S(a,b) file
has_sab = False 
material_number_sab = 0
if (args.path_to_sab != ""):
  has_sab = True
  os.system("ln -fs "+args.path_to_sab+" tape51")

  with open("tape51","r") as sab_file:
    lines = sab_file.readlines()
    for l, line in enumerate(lines):
      line = line.split()
      if "MATERIAL" in line:
        material_number_sab = int(line[2])
        break

  print("S(a,b) material number read: "+str(material_number_sab))

#====================================== Define output filename
output_filename = args.output_filename
if (output_filename == ""):
  output_filename = atomic_symbol+str(mass_number)+".njoy"

output_directory = args.output_directory
if (not os.path.isdir(output_directory)):
  os.makedirs(output_directory)

#====================================== Begin writing input
with open("NJOY_INPUT.txt","w") as njoy_input:
  njoy_input.write("-- Processing ENDF to PENDF\n")
  njoy_input.write("moder\n")
  njoy_input.write("20 -21/\n")

  if (args.path_to_gamma_endf != ""):
    njoy_input.write("moder\n")
    njoy_input.write("52 -53/\n")

  njoy_input.write("reconr\n")
  njoy_input.write("-21 -22/\n")
  njoy_input.write("'pendf tape'/\n")
  njoy_input.write(str(material_number)+" 0/\n")
  njoy_input.write("0.001/\n")
  njoy_input.write("0/\n")

  if (args.path_to_gamma_endf != ""):
    njoy_input.write("reconr\n")
    njoy_input.write("-53 -54/\n")
    njoy_input.write("'pendf gamma tape'/\n")
    njoy_input.write(str(material_number_gamma)+" 0/\n")
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
  njoy_input.write("0 ")  #elastic options
                          # 0=none
                          # 1=compute using ENDF6 format data (if supplied)
  njoy_input.write("0 1 221 1/\n")    
  njoy_input.write(str(args.temperature)+"/\n") 
  njoy_input.write("0.005 5.0/\n")   

  #====================================== Thermr S(a,b)
  if (has_sab):
    njoy_input.write("thermr\n")
    njoy_input.write("51 -25 -26/\n")
    njoy_input.write(str(material_number_sab)+" "+str(material_number)+" 16 1 ")
    njoy_input.write("2 ")  #inelastic options 
                            # 0=none
                            # 1=compute as free gas
                            # 2=read s(a,b) and compute matrix
    njoy_input.write("1 ")  #elastic options
                            # 0=none
                            # 1=compute using ENDF6 format data (if supplied)
    njoy_input.write("0 "+args.inelastic_thermal_num_atoms)      #Number principal atoms
    njoy_input.write(" "+args.inelastic_thermal_number+" 0/\n")  #MT number for inelastic  
    njoy_input.write(str(args.temperature)+"/\n") 
    njoy_input.write("0.005 5.0/\n") 

  #====================================== Groupr
  njoy_input.write("groupr\n")
  if (has_sab):
    njoy_input.write("-21 -26 0 -30\n")
  else:
    njoy_input.write("-21 -25 0 -30\n")
  njoy_input.write(str(material_number)+" ")
  njoy_input.write(str(args.neutron_group_structure)+" ")
  njoy_input.write(str(args.gamma_group_structure)+" ")
  njoy_input.write(str(args.neutron_weight_function)+" ")
  njoy_input.write("7 ")              # Legendre order
  njoy_input.write("1 1 1 1/\n")
  njoy_input.write("'"+atomic_symbol+str(mass_number)+isotope_metastable+"'/\n")
  njoy_input.write(str(args.temperature)+"/\n") 
  njoy_input.write("0.0/\n") 

  # Custom neutron gs
  if (args.neutron_group_structure == 1):
    with open(args.custom_neutron_gs_file,"r") as n_gs_file:
      lines = n_gs_file.readlines()
      for line in lines:
        if (line[0] != '#'):
          njoy_input.write(line)
    njoy_input.write("\n")

  # Custom gamma gs
  if (args.gamma_group_structure == 1):
    with open(args.custom_gamma_gs_file,"r") as g_gs_file:
      lines = g_gs_file.readlines()
      for line in lines:
        if (line[0] != '#'):
          njoy_input.write(line)
    njoy_input.write("\n")

  # Custom neutron weight function
  if (args.neutron_weight_function == 1): 
    with open(args.custom_neutron_wt_file,"r") as n_wt_file:
      lines = n_wt_file.readlines()
      for line in lines:
        if (line[0] != '#'):
          njoy_input.write(line)
    njoy_input.write("\n")

  njoy_input.write("3/\n") #All MF3 1D cross-sections
  njoy_input.write("3 259 'inverse velocity'/\n") 
  njoy_input.write("3 221 'free-gas'/\n") 
  if (has_sab):
    njoy_input.write("3 " + args.inelastic_thermal_number + \
                     " 'inelastic s(a,b)'/\n")
    if (args.elastic_thermal_number != ""):
      njoy_input.write("3 " + args.elastic_thermal_number + \
                       " 'elastic s(a,b)'/\n")
  if (has_fission):
    njoy_input.write("3 452 'total nubar'/\n")
    njoy_input.write("3 456 'prompt nubar'/\n")
    njoy_input.write("3 455 'delayed nubar'/\n")
    njoy_input.write("5 18 'prompt chi'/\n")
    njoy_input.write("5 455 'delayed chi'/\n")
  njoy_input.write("6/\n") #All MF6 transfer matrices
  njoy_input.write("6 221 'free-gas matrix'/\n")

  if (has_sab):
    njoy_input.write("6 " + args.inelastic_thermal_number + \
                     " 'inelastic_s(a,b) matrix'/\n")
    if (args.elastic_thermal_number != ""):
      njoy_input.write("6 " + args.elastic_thermal_number + \
                       " 'elastic_s(a,b) matrix'/\n")

  #Gamma matrices
  if (args.path_to_gamma_endf != ""):
      njoy_input.write("16/\n")
      # njoy_input.write("17/\n")
      # njoy_input.write("18/\n")

  njoy_input.write("0/\n") #Terminate reaction list
  njoy_input.write("0/\n") #Terminate groupr

  #====================================== GAMINR
  if (args.path_to_gamma_endf != ""):
    njoy_input.write("gaminr\n")
    njoy_input.write("-53 -54 0 -55/\n")
    #Card 2
    njoy_input.write(str(material_number_gamma)+" ")
    njoy_input.write(str(args.gamma_group_structure)+" ")
    njoy_input.write(str(args.gamma_weight_function)+" ")
    njoy_input.write("7 ")              # Legendre order
    njoy_input.write("1/\n")              # Print option
    njoy_input.write("'Photon interactions'/\n")

    # Custom gamma gs
    if (args.gamma_group_structure == 1):
      with open(args.custom_gamma_gs_file,"r") as g_gs_file:
        lines = g_gs_file.readlines()
        for line in lines:
          if (line[0] != '#'):
            njoy_input.write(line)

    # Custom gamma weight function
    if (args.gamma_weight_function == 1):
      with open(args.custom_gamma_wt_file,"r") as g_wt_file:
        lines = g_wt_file.readlines()
        for line in lines:
          if (line[0] != '#'):
            njoy_input.write(line)

    njoy_input.write("-1/\n") #All data
    njoy_input.write("0/\n")  #Terminate gaminr

  #====================================== Convert groupr output to ASCII
  njoy_input.write("moder\n")
  njoy_input.write("-30 31/\n")

  #====================================== Convert gaminr output to ASCII
  if (args.path_to_gamma_endf != ""):
    njoy_input.write("moder\n")
    njoy_input.write("-55 56/\n")

  njoy_input.write("stop\n")
  njoy_input.close()

#====================================== Run NJOY
os.system("njoy < NJOY_INPUT.txt")
os.system("rm tape* NJOY_INPUT.txt")
print("Copying outputfile to "+output_directory+'/'+output_filename)

# if (args.path_to_gamma_endf != ""):
#     os.system("cat tape56 >>")

os.system("cp output "+output_directory+'/'+output_filename)
os.system("rm output")