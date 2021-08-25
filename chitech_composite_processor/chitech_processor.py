"""Executes the main steps of conversion"""
import Utils_ChiTechCombiner
import Utils_Info
#import Utils_ReadNJOYOutput
#import Utils_GsProcessor
import Utils_CombinedXSWriter
import Utils_ChiTechPlotter
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
argparser.add_argument("--molar_mass",
                       help="The molar mass for the element or composite",
                       type=float, required=True)

# ====================================== Argument for composite
argparser.add_argument("--composite_chi_output_filename",
                       help="List of output files produced by ChiTech for the all elements in the composite",
                       default="", required=True)
argparser.add_argument("--composite_atomic_fraction",
                       help="List of the atomics fraction of all elements in the composite",
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

# ============================= Argument for additional Njoy plots
argparser.add_argument("--more_njoy",
                       help="If included, will produce plots of additional njoy data",
                       action='store_true', required=False)
argparser.add_argument("--more_composite_chi_output_filename",
                       help="List of output files produced by ChiTech for the all elements in the composite",
                       default="", required=True)

args = argparser.parse_args()    


#Get atomic fraction
atom_fraction=[]
atom_fraction_string=args.composite_atomic_fraction.split(",")
print(atom_fraction_string)
for i in atom_fraction_string:
    atom_fraction.append(int(i))

# ===================================== Combine disjoint data from ChiTech file
print("Combining disjoint data")
composite_Chifilename = args.composite_chi_output_filename.split(",")
print(composite_Chifilename)
for i in range (len(composite_Chifilename)):
    composite_Chifilename[i] = args.output_path + "/" + composite_Chifilename[i]
data = Utils_ChiTechCombiner.BuildCombinedData(composite_Chifilename, atom_fraction)

# ===================================== Write cross-section
chi_output_complete_path = args.output_path + "/" + args.chixs_filename
print("Creating chi-cross-section in file " + chi_output_complete_path)
Utils_CombinedXSWriter.WriteChiTechFile(data, chi_output_complete_path, comment = "# Output")

# ===================================== Generate spectrum
# FIXME: Check if this works with gamma-only problem later

# Check if neutron problem
neutron_gs = data["neutron_gs"]
if neutron_gs == []:
    print("No spectrum for gamma only problem")
else:
    import Utils_Info
    print(args.mcnp)
    if args.mcnp:
        mcnp_output_complete_path = args.mcnp_path + args.mcnp_filename
    else:
        mcnp_output_complete_path = ""

    # ================================== Process additional njoy data if needed ======
    more_njoy_dict = {}
    if args.more_njoy:
        # ===================================== Combine disjoint data from ChiTech file
        print("Combining disjoint data")
        composite_Chifilename = args.more_composite_chi_output_filename.split(",")
        print(composite_Chifilename)
        for i in range (len(composite_Chifilename)):
            composite_Chifilename[i] = args.output_path + "/" + composite_Chifilename[i]
        extra_data = Utils_ChiTechCombiner.BuildCombinedData(composite_Chifilename, atom_fraction)

    processed_data = Utils_Info.InfiniteMediumSpectrum(data, path=args.output_path, plot=True)
    # FIXME: replaced with args.plot when everything is fixed
    if args.plot: 
        print("Plotting the Spectra")
        Utils_ChiTechPlotter.PlotSpectra(processed_data, path=args.output_path, molar_mass=args.molar_mass, mcnp_filename=mcnp_output_complete_path, more_njoy=args.more_njoy, extra_njoy_data = extra_data)

