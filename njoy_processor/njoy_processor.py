"""Executes the main steps of conversion"""
import Utils_ReadNJOYOutput
import Utils_Combiner
import Utils_XSWriter

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

argparser.add_argument("--path_to_njoy_output",
                       help="Path to output produced by NJOY",
                       default="", required=True)
argparser.add_argument("--output_file_path",
                       help="Complete path and name of output xs",
                       default="", required=True)
argparser.add_argument("--plot",
                       help="If included, will produce xs plots",
                       action='store_true', required=False)

args = argparser.parse_args()                                            


# ===================================== Read NJOY output
print("Reading NJOY output")
raw_njoy_data = Utils_ReadNJOYOutput.ReadNJOYfile(
    args.path_to_njoy_output, verbose=args.plot)

# ===================================== Combine disjoint data
print("Combining disjoint data")
data = Utils_Combiner.BuildCombinedData(raw_njoy_data, plot=args.plot)

# ===================================== Write cross-section
print("Creating chi-cross-section in file "+args.output_file_path)
Utils_XSWriter.WriteChiTechFile(data, args.output_file_path)
