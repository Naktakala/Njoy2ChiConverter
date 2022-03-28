"""Executes the main steps of conversion"""
import Utils_ChiTechCombiner
import Utils_Info
import Utils_CombinedXSWriter
import Utils_ChiTechPlotter
import Utils_Info
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
argparser.add_argument("--source_term",
                       help="Description of the source term in string form: particle type, energy value in MeV and if fission is included, separated by a comma",
                       default="", required=True)
argparser.add_argument("--plot",
                       help="If included, will produce xs plots",
                       action='store_true', required=False)

# ====================================== Argument for composite
argparser.add_argument("--composite_chi_output_filename",
                       help="List of output files produced by ChiTech for the composite separated by group structure with || in between each gs",
                       default="", required=True)
argparser.add_argument("--atomic_fraction",
                       help="List of the atomics fraction of all elements in the composite",
                       default="", required=True)
argparser.add_argument("--atomic_density",
                       help="List of the atomics density of all elements in the composite",
                       default="", required=True)
# ============================= Argument for MCNP
argparser.add_argument("--mcnp",
                       help="If included, will produce plots of mcnp data",
                       action='store_true', required=False)
argparser.add_argument("--mcnp_path",
                       help="Complete path where the mcnp output is stored",
                       default="", required=False)
argparser.add_argument("--mcnp_filename",
                       help="Name of output file produced by MCNP",
                       default="", required=False)

args = argparser.parse_args()    


#Get atomic fraction
atom_fraction=[]
atom_fraction_string=args.atomic_fraction.split(",")
print(atom_fraction_string)
for i in atom_fraction_string:
    atom_fraction.append(float(i))

#Get atomic density
N_density=[]
N_density_string=args.atomic_density.split(",")
print(N_density_string)
for i in N_density_string:
    N_density.append(float(i))

# ===================================== Combine disjoint data from ChiTech file
big_data = []
print("Combining disjoint data")
chitech_gs_list = args.composite_chi_output_filename.split("||")
for chi_gs in chitech_gs_list:
    composite_Chifilename = chi_gs.split(",")
    print(composite_Chifilename)
    for i in range (len(composite_Chifilename)):
        composite_Chifilename[i] = args.output_path + "/" + composite_Chifilename[i]
    data = Utils_ChiTechCombiner.BuildCombinedData(composite_Chifilename, atom_fraction, N_density)
    big_data.append(data)

# ===================================== Write cross-section
# FIXME: Possibly add an option for user to pick which group structure to write in chitech format
chi_output_complete_path = args.output_path + "/" + args.chixs_filename
print("Creating chi-cross-section in file " + chi_output_complete_path)
data = big_data[0]
Utils_CombinedXSWriter.WriteCombinedChiTechFile(data, chi_output_complete_path, comment = "# Output")

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

