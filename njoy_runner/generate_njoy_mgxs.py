'''Prepares input/output files for NJOY2016 and executes NJOY2016/2021'''
import sys 
import argparse
from argparse import RawTextHelpFormatter
import os

# ===================================== Check python version
if sys.version_info[0] < 3:
    print("\n Error: This script requires python3 but was executed "
          "with version:\n\n"+sys.version+"\n")
    sys.exit(1)

# ===================================== Setup arguments
ign  = [*range(1,35)]  # Allowed numeric options for --neutron_group_structure
niwt = [*range(1,13)]  # Allowed numeric options for --neutron_weight_function
igg  = [*range(0,11)]  # Allowed numeric options for --gamma_group_structure
giwt = [*range(1,4)]   # Allowed numeric options for --gamma_weight_function
argparser = argparse.ArgumentParser(
    description="Prepares input/output files "
                "for NJOY2016 and executes NJOY2016",
    formatter_class=RawTextHelpFormatter,
    epilog='''Additional information:
  Neutron group structure options (ign):
   1   arbitrary structure (read in)
   2   csewg 239-group structure
   3   lanl 30-group structure
   4   anl 27-group structure
   5   rrd 50-group structure
   6   gam-i 68-group structure
   7   gam-ii 100-group structure
   8   laser-thermos 35-group structure
   9   epri-cpm 69-group structure
   10  lanl 187-group structure
   11  lanl 70-group structure
   12  sand-ii 620-group structure
   13  lanl 80-group structure
   14  eurlib 100-group structure
   15  sand-iia 640-group structure
   16  vitamin-e 174-group structure
   17  vitamin-j 175-group structure
   18  xmas nea-lanl
   19  ecco 33-group structure
   20  ecco 1968-group structure
   21  tripoli 315-group structure
   22  xmas lwpc 172-group structure
   23  vit-j lwpc 175-group structure
   24  shem cea 281-group structure
   25  shem epm 295-group structure
   26  shem cea/epm 361-group structure
   27  shem epm 315-group structure
   28  rahab aecl 89-group structure
   29  ccfe 660-group structure (30 MeV)
   30  ukaea 1025-group structure (30 MeV)
   31  ukaea 1067-group structure (200 MeV)
   32  ukaea 1102-group structure (1 GeV)
   33  ukaea 142-group structure (200 MeV)
   34  lanl 618-group structure
   
   Neutron group structure weighting options (iwt):
   1   read in smooth weight function
   2   constant
   3   1/e
   4   1/e + fission spectrum + thermal maxwellian
   5   epri-cell lwr
   6   (thermal) -- (1/e) -- (fission + fusion)
   7   same with t-dep thermal part
   8   thermal--1/e--fast reactor--fission + fusion
   9   claw weight function
   10  claw with t-dependent thermal part
   11  vitamin-e weight function (ornl-5505)
   12  vit-e with t-dep thermal part
   
  Gamma group structure options (igg):
   0   none
   1   arbitrary structure (read in)
   2   csewg 94-group structure
   3   lanl 12-group structure
   4   steiner 21-group gamma-ray structure
   5   straker 22-group  structure
   6   lanl 48-group structure
   7   lanl 24-group structure
   8   vitamin-c 36-group structure
   9   vitamin-e 38-group structure
   10  vitamin-j 42-group structure
   
   Gamma group structure weighting options (giwt):
   1   read in smooth weight function
   2   constant
   3   1/e + rolloffs
   
Example Custom group structure file:
------
# Number of groups
31/
# Group boundaries
1.0000e-05 4.6589e+05 9.3178e+05 1.3977e+06 1.8636e+06 2.3295e+06
2.7953e+06 3.2612e+06 3.7271e+06 4.1930e+06 4.6589e+06 5.1248e+06
5.5907e+06 6.0566e+06 6.5225e+06 6.9884e+06 7.4775e+06 7.9434e+06
8.4093e+06 8.8752e+06 9.3411e+06 9.8070e+06 1.0273e+07 1.0739e+07
1.1205e+07 1.1671e+07 1.2136e+07 1.2602e+07 1.3068e+07 1.3534e+07
1.3999e+07 1.4000e+07/
------

Example Custom smooth weight function:
See spectrum_file.txt
   
''')

argparser.add_argument("--njoy_exec_name",
                       help="Name of the NJOY executable",
                       default="njoy")                                              
argparser.add_argument("--path_to_neutron_endf",
                       help="Path to the isotope's endf file for "
                       "incident neutrons",
                       default="")
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
                       help="NJOY temperature output", type=float,
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

# ===================================== Check argument rules
if args.path_to_neutron_endf != "":
  if not os.path.isfile(args.path_to_neutron_endf):
    print("Error: Value supplied to --path_to_neutron_endf does not "
          "point to an existing file.")

if (args.neutron_group_structure == 1 and 
        args.custom_neutron_gs_file == ""):
    print("Error: When specifying \"1\" for neutron group structure "
          "then argument \"custom_neutron_gs_file\" must be supplied.")

if (args.neutron_group_structure == 1 and
        args.custom_neutron_gs_file != ""):
    if not os.path.isfile(args.custom_neutron_gs_file):
        print("Error: Value supplied to --custom_neutron_gs_file does not "
              "point to an existing file.")

if args.path_to_sab != "":
  if not os.path.isfile(args.path_to_sab):
    print("Error: Value supplied to --path_to_sab does not "
          "point to an existing file.")

with_gamma = False
if args.path_to_gamma_endf != "":
    with_gamma = True
    if not os.path.isfile(args.path_to_gamma_endf):
        print("Error: Value supplied to --path_to_gamma_endf does not "
              "point to an existing file.")

if (with_gamma and args.gamma_group_structure == 1 and 
        args.custom_gamma_gs_file == ""):
    print("Error: When specifying \"1\" for gamma group structure "
          "then argument \"custom_gamma_gs_file\" must be supplied.")

if (with_gamma and args.gamma_group_structure == 1 and
        args.custom_gamma_gs_file != ""):
    if not os.path.isfile(args.custom_gamma_gs_file):
        print("Error: Value supplied to --custom_gamma_gs_file does not point to"
              "an existing file.")

if (args.neutron_weight_function == 1 and 
        args.custom_neutron_wt_file == ""):
    print("Error: When specifying \"1\" for neutron weight function "
          "then argument \"custom_neutron_wt_file\" must be supplied.")

if (args.neutron_weight_function == 1 and
        args.custom_neutron_wt_file != ""):
    if not os.path.isfile(args.custom_neutron_wt_file):
        print("Error: Value supplied to --custom_neutron_wt_file does not "
              "point to an existing file.")

if (with_gamma and args.gamma_weight_function == 1 and 
        args.custom_gamma_wt_file == ""):
    print("Error: When specifying \"1\" for gamma weight function "
          "then argument \"custom_gamma_wt_file\" must be supplied.")

if (with_gamma and args.gamma_weight_function == 1 and
        args.custom_gamma_wt_file != ""):
    if not os.path.isfile(args.custom_gamma_wt_file):
        print("Error: Value supplied to --custom_gamma_wt_file does not "
              "point to an existing file.")

if not os.path.isdir(args.output_directory):
    print("Error: Value passed to --output_directory does not point to "
          "an existing folder")

# ===================================== Tape reservations
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

# ================================ Display type of problem ======
if args.path_to_neutron_endf != "":
  if args.path_to_gamma_endf != "":
    print("Neutron + Gamma Problem")
  else:
    print("Neutron only Problem")
elif args.path_to_gamma_endf != "":
  print("Gamma only Problem")
else:
  raise ValueError("Either a neutron or a gamma ENDF must be provided") 

# ===================================== Develop input tapes
if args.path_to_neutron_endf != "":
  os.system("ln -fs "+args.path_to_neutron_endf+" tape20")

if args.path_to_gamma_endf != "":
    os.system("ln -fs "+args.path_to_gamma_endf+" tape52")
    


# ===================================== Extract isotope information from
#                                       endf neutron file
if args.path_to_neutron_endf != "":
  endf_file = open("tape20", "r")
  endf_file.readline()    # skip line
  endf_file.readline()    # skip line
  endf_file.readline()    # skip line
  endf_file.readline()    # skip line
  endf_file.readline()    # skip line
  
  isotope_text = endf_file.readline(11); endf_file.readline()
  isotope_atomic_number = int(isotope_text[0:3])
  isotope_symbol        = isotope_text[4:6].strip()
  isotope_mass_number   = int(isotope_text[7:10])
  isotope_metastable    = isotope_text[10:].strip()
  
  endf_file.readline()    # skip line
  endf_file.readline(11)  # skip field
  endf_file.readline(11)  # skip field
  material_number = int(endf_file.readline(22)[9:])
  
  output_filename = args.output_filename
  if args.output_filename == "":
      output_filename = isotope_symbol + str(isotope_mass_number) + \
                        isotope_metastable+".txt"
  
  endf_file.close()

# ===================================== Extract gamma number information from
#                                       endf gamma file, if required
material_number_gamma = 0
if args.path_to_gamma_endf != "":
  os.system("ln -fs "+args.path_to_gamma_endf+" tape52")
  endf_gamma_file = open("tape52", "r")
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline()    # skip line
  endf_gamma_file.readline(11)  # skip field
  endf_gamma_file.readline(11)  # skip field

  material_number_gamma = int(endf_gamma_file.readline(22)[9:])

  print("gamma material number read: "+str(material_number_gamma))
  
  endf_gamma_file.close()
  
  # ============== Only do for gamma only problem ==========
  if args.path_to_neutron_endf =="":
    os.system("ln -fs "+args.path_to_gamma_endf+" tape52")
    endf_gamma_file = open("tape52", "r")
    endf_gamma_file.readline()    # skip line
    endf_gamma_file.readline()    # skip line
    endf_gamma_file.readline()    # skip line
    endf_gamma_file.readline()    # skip line
    endf_gamma_file.readline()    # skip line
    
    isotope_text = endf_gamma_file.readline(11); endf_gamma_file.readline()
    isotope_atomic_number = int(isotope_text[0:3])
    isotope_symbol        = isotope_text[4:6].strip()
    isotope_mass_number   = int(isotope_text[7:10])
    
    output_filename = args.output_filename
    if args.output_filename == "":
      output_filename = isotope_symbol + str(isotope_mass_number) + \
                          isotope_metastable+".txt"
    
    endf_gamma_file.close()

# ===================================== Extract information from S(a,b) file
#                                       if required
with_sab = False 
material_number_sab = 0
if args.path_to_sab != "":
  with_sab = True

  sab_file = open(args.path_to_sab, "r")
  os.system("ln -fs "+args.path_to_sab+" tape51")

  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline()    # skip line
  sab_file.readline(11)  # skip field
  sab_file.readline(11)  # skip field
  material_number_sab = int(sab_file.readline(22)[9:])

  print("S(a,b) material number read: "+str(material_number_sab))

  sab_file.close()

# ===================================== Begin writing input
njoy_input = open("NJOY_INPUT.txt", "w")
njoy_input.write("-- Processing ENDF to PENDF\n")

# =============================== Neutron-dependent input
if args.path_to_neutron_endf != "":

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

  # ===================================== Thermr free-gas
  njoy_input.write("thermr\n")
  njoy_input.write("0 -90 -25/\n")
  njoy_input.write("0 "+str(material_number)+" 16 1 ")
  njoy_input.write("1 ")     # inelastic options
                             #  0=none
                             #  1=compute as free gas
                             #  2=read s(a,b) and compute matrix
  njoy_input.write("0 ")     # elastic options
                             #  0=none
                             #  1=compute using ENDF6 format data (if supplied)
  njoy_input.write("0 1 221 1/\n")    
  njoy_input.write(str(args.temperature)+"/\n") 
  njoy_input.write("0.005 5.0/\n")   

# ==================================== Gamma-dependent Input
if args.path_to_gamma_endf != "":
  njoy_input.write("moder\n")
  njoy_input.write("52 -53/\n")
  
  njoy_input.write("reconr\n")
  njoy_input.write("-53 -54/\n")
  njoy_input.write("'pendf gamma tape'/\n")
  njoy_input.write(str(material_number_gamma)+" 0/\n")
  njoy_input.write("0.001/\n")
  njoy_input.write("0/\n")

# ===================================== Thermr S(a,b)
if with_sab:
  njoy_input.write("thermr\n")
  njoy_input.write("51 -25 -26/\n")
  njoy_input.write(str(material_number_sab)+" "+str(material_number)+" 16 1 ")
  njoy_input.write("2 ")  # inelastic options
                          #  0=none
                          #  1=compute as free gas
                          #  2=read s(a,b) and compute matrix
  njoy_input.write("0 ")  # elastic options
                          #  0=none
                          #  1=compute using ENDF6 format data (if supplied)
  njoy_input.write("0 "+args.inelastic_thermal_num_atoms)      #Number principal atoms
  njoy_input.write(" "+args.inelastic_thermal_number+" 0/\n")  #MT number for inelastic  
  njoy_input.write(str(args.temperature)+"/\n") 
  njoy_input.write("0.005 5.0/\n") 

# ===================================== Groupr (Only for Neutron included problems)
if args.path_to_neutron_endf != "":
  njoy_input.write("groupr\n")
  if with_sab:
      njoy_input.write("-21 -26 0 -30\n")
  else:
      njoy_input.write("-21 -25 0 -30\n")
  njoy_input.write(str(material_number)+" ")
  njoy_input.write(str(args.neutron_group_structure)+" ")
  njoy_input.write(str(args.gamma_group_structure)+" ")
  njoy_input.write(str(args.neutron_weight_function)+" ")
  njoy_input.write("7 ")              # Legendre order
  njoy_input.write("1 1 1 1/\n")
  njoy_input.write("'"+isotope_symbol+str(isotope_mass_number)+isotope_metastable+"'/\n")
  njoy_input.write(str(args.temperature)+"/\n") 
  njoy_input.write("0.0/\n") 
  
  # Custom neutron gs
  if args.neutron_group_structure == 1:
      n_gs_file = open(args.custom_neutron_gs_file, "r")
      lines = n_gs_file.readlines()
      for line in lines:
          if line[0] != '#':
              njoy_input.write(line)
      n_gs_file.close()
  
  # Custom gamma gs
  if args.gamma_group_structure == 1:
      g_gs_file = open(args.custom_gamma_gs_file, "r")
      lines = g_gs_file.readlines()
      for line in lines:
          if line[0] != '#':
              njoy_input.write(line)
      g_gs_file.close()
  
  # Custom neutron weight function
  if args.neutron_weight_function == 1:
      n_wt_file = open(args.custom_neutron_wt_file, "r")
      lines = n_wt_file.readlines()
      for line in lines:
          if line[0] != '#':
              njoy_input.write(line)
      n_wt_file.close()
  
  njoy_input.write("3/\n")  # All MF3 1D cross-sections
  njoy_input.write("3 259 'inverse velocity'/\n") 
  njoy_input.write("3 221 'free-gas therm scat'/\n") 
  
  njoy_input.write("6/\n")  # All MF6 transfer matrices
  njoy_input.write("6 221 'free-gas matrix'/\n")
  
  if with_sab:
      njoy_input.write("6 "+args.inelastic_thermal_number+" 'inelastic_s(a,b) matrix'/\n")
      if args.elastic_thermal_number != "":
          njoy_input.write("6 "+args.elastic_thermal_number+" 'elastic_s(a,b) matrix'/\n")
  
  # Gamma matrices
  if args.path_to_gamma_endf != "":
      njoy_input.write("16/\n")
      # njoy_input.write("17/\n")
      # njoy_input.write("18/\n")
  
  njoy_input.write("0/\n")  # Terminate reaction list
  
  njoy_input.write("0/\n")  # Terminate groupr
  
  # ===================================== Convert groupr output to ASCII
  njoy_input.write("moder\n")
  njoy_input.write("-30 31/\n")


# ===================================== GAMINR
if args.path_to_gamma_endf != "":
    njoy_input.write("gaminr\n")
    njoy_input.write("-53 -54 0 -55/\n")
    # Card 2
    njoy_input.write(str(material_number_gamma)+" ")
    njoy_input.write(str(args.gamma_group_structure)+" ")
    njoy_input.write(str(args.gamma_weight_function)+" ")
    njoy_input.write("7 ")                # Legendre order
    njoy_input.write("1/\n")              # Print option
    njoy_input.write("'Photon interactions'/\n")

    # Custom gamma gs
    if args.gamma_group_structure == 1:
        g_gs_file = open(args.custom_gamma_gs_file, "r")
        lines = g_gs_file.readlines()
        for line in lines:
            if line[0] != '#':
                njoy_input.write(line)
        g_gs_file.close()

    # Custom gamma weight function
    if args.gamma_weight_function == 1:
        g_wt_file = open(args.custom_gamma_wt_file, "r")
        lines = g_wt_file.readlines()
        for line in lines:
            if line[0] != '#':
                njoy_input.write(line)
        g_wt_file.close()

    njoy_input.write("-1/\n") #All data
    njoy_input.write("0/\n")  #Terminate gaminr

# ===================================== Convert gaminr output to ASCII
    njoy_input.write("moder\n")
    njoy_input.write("-55 56/\n")

njoy_input.write("stop\n")
njoy_input.close()


# ===================================== Run NJOY
# remove output file
os.system('rm -f output')

if args.njoy_exec_name == "njoy21":
    cmd_line = args.njoy_exec_name + " -i NJOY_INPUT.txt -o output"
else:
    cmd_line = args.njoy_exec_name  + " < NJOY_INPUT.txt"
print('command line =',cmd_line)

os.system(cmd_line)
os.system("rm tape*")

if not os.path.isdir(args.output_directory):
  print("Error: Value passed to --output_directory does not point to "
        "an existing folder")

forward_slash = ""
if args.output_directory[-1] != "/":
  forward_slash = "/"

print("Copying outputfile to "+args.output_directory+forward_slash+output_filename)

os.system("cp output "+args.output_directory+forward_slash+output_filename)
# os.system("pwd")
# os.system("cp NJOY_INPUT.txt "+args.output_directory+"NJOY_INPUT.txt")