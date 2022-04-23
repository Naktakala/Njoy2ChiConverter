"""Executes the main steps of conversion"""
import Utils_ChiTechCombiner
import Utils_CombinedXSWriter
#import Utils_Info
#import Utils_ChiTechPlotter

import sys 
import argparse

# ===================================== Check version
if sys.version_info[0] < 3:
    print("\n Error: This script requires python3 but was executed "
          "with version:\n\n"+sys.version+"\n")
    sys.exit(1)

# ===================================== Setup arguments
argparser = argparse.ArgumentParser(
    description="Executes the main steps of conversion")

argparser.add_argument("--output_path",
                       help="Complete path where to store the output xs",
                       default="", required=True)
argparser.add_argument("--chixs_filename",
                       help="Name of Chi XS file",
                       default="", required=True)
argparser.add_argument("--plot",
                       help="If included, will produce xs plots",
                       action='store_true', required=False)

# ====================================== Argument for composite
argparser.add_argument("--chixs_filename_list",
                       help="List of output files produced by ChiTech. Must use same gs",
                       default="", required=True)
argparser.add_argument("--atomic_density",
                       help="List of the atomic densities",
                       default="", required=True)
# ============================= Argument for MCNP

args = argparser.parse_args()    

forward_slash = ""
if args.output_path[-1] != "/":
    forward_slash = "/"
    
# Get atomic density
N_density=[]
N_density_string=args.atomic_density.split(",")
print(N_density_string)
for i in N_density_string:
    N_density.append(float(i))

# ===================================== Combine data from ChiTech files
print("Combining ChiTech data")
chixs_filename_list = args.chixs_filename_list.split(",")
if len(chixs_filename_list) != len(N_density):
    raise ValueError('number of filenames inconsistent with at.density data')

chixs_fullpath_list = []
for chixs_filename in chixs_filename_list:
    print(chixs_filename)
    chixs_fullpath_list.append(args.output_path + forward_slash + chixs_filename)

data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath_list, N_density)

# ===================================== Write cross section

chi_output_complete_path = args.output_path + forward_slash + args.chixs_filename
print("Creating chi-cross-section in file " + chi_output_complete_path)
Utils_CombinedXSWriter.WriteCombinedChiTechFile(data, chi_output_complete_path)

"""
argparser.add_argument("--mcnp",
                       help="If included, will produce plots of mcnp data",
                       action='store_true', required=False)
argparser.add_argument("--mcnp_path",
                       help="Complete path where the mcnp output is stored",
                       default="", required=False)
argparser.add_argument("--mcnp_filename",
                       help="Name of output file produced by MCNP",
                       default="", required=False)
"""
"""
argparser.add_argument("--source_term",
                       help="Description of the source term in string form: particle type, energy value in MeV and if fission is included, separated by a comma",
                       default="", required=True)
"""

"""
# ===================================== Create the source definition
source_def={}
source_string=args.source_term.split(",")
print(source_string)
source_def["Particle type"] = source_string[0]
source_def["Energy value"] = source_string[1]
if len(source_string) > 2:
    source_def["If fission"] = True
else:
    source_def["If fission"] = False

print("Source definition: ", end = "")
for keys in source_def.keys():
    print(str(keys) + ": " + str(source_def[keys]), end = ", ")
print()

# ===================================== Solve the infinite medium equation
big_processed_data = []
for data in big_data:
    processed_data = Utils_Info.InfiniteMediumSpectrum(data, source_def,path=args.output_path, plot=args.plot)
    big_processed_data.append(processed_data)

# ===================================== Generate spectrum
#================================== Get MCNP data if asked
if args.mcnp:
    mcnp_output_complete_path = args.mcnp_path + args.mcnp_filename
else:
    mcnp_output_complete_path = ""

if args.plot: 
    print("Plotting the Spectra")
    Utils_ChiTechPlotter.PlotSpectra(big_processed_data, source_def, path=args.output_path, mcnp_filename=mcnp_output_complete_path)
"""
