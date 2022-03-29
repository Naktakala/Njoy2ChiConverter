"""Executes the main steps of conversion"""
import Utils_ReadNJOYOutput
import Utils_Combiner
import Utils_XSWriter
# import Utils_NjoySpectrumPlotter
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
argparser.add_argument("--njoy_output_filename",
                       help="Name of output file produced by NJOY",
                       default="", required=True)
argparser.add_argument("--chixs_filename",
                       help="Name of Chi XS file",
                       default="", required=False)
argparser.add_argument("--plot",
                       help="If included, will produce xs plots",
                       default=False,
                       action='store_true', required=False)
args = argparser.parse_args()                                            

forward_slash = ""
if args.output_path[-1] != "/":
    forward_slash = "/"
  
# ===================================== Read NJOY output
njoy_output_complete_path = args.output_path + "/" + args.njoy_output_filename
print("Reading NJOY output located at " + njoy_output_complete_path)
raw_njoy_data = Utils_ReadNJOYOutput.ReadNJOYfile(
    njoy_output_complete_path, verbose=True)

# ===================================== Combine disjoint data
print("Combining disjoint data")
data, problem_description = Utils_Combiner.BuildCombinedData(raw_njoy_data, plot=args.plot)

filename =  args.njoy_output_filename.split("_")
### FIMXME: this is not a good way of recovering the isotope if the user supplies its own output filename
problem_description['isotope'] = filename[0]
print('\n\n Problem description =\n',problem_description,'\n\n')

# ===================================== Write cross-section file in Chi format
if len(args.chixs_filename)==0:
    # generate name automatically
    chixs_filename = args.njoy_output_filename + ".csx"
else:
    chixs_filename = args.chixs_filename
   

chi_output_complete_path = args.output_path + forward_slash + chixs_filename
print("Creating chi cross section in file " + chi_output_complete_path)
Utils_XSWriter.WriteChiTechFile(data, chi_output_complete_path, problem_description)

# copy the plots
if args.plot:
    import os, glob
    filename =  chixs_filename.split(".")  
    # os.system('pwd')    
    # os.system('ls -lasrt')    
    for file in glob.glob("*.png"):
        #print(file)
        os.system("mv "+file+" "+args.output_path+forward_slash+filename[0]+"_"+file)

# ===================================== Done

# argparser.add_argument("--source_term",
                       # help="Description of the source term in string form: particle type and energy value in MeV separated by a comma",
                       # default="", required=True)
# Get the source definition
# source_def={}
# source_string=args.source_term.split(",")
# print(source_string)
# source_def["Particle type"] = source_string[0]
# source_def["Energy value"] = float(source_string[1])
# if len(source_string) > 2:
    # source_def["If fission"] = True
# else:
    # source_def["If fission"] = False

# print("Source definition: ", end = "")
# for keys in source_def.keys():
    # print(str(keys) + ": " + str(source_def[keys]), end = ", ")
# print()

# ===================================== Generate spectrum
# sys.path.insert(1,'../chitech_composite_processor')
# import Utils_Info

# out_p = Utils_Info.InfiniteMediumSpectrum(data, source_def, path=args.output_path, plot=args.plot)

# FIXME: Remove if dont want to plot completely
# if args.plot:
    # Utils_NjoyPlotter.Njoy_plotter(out_p, path=args.output_path)


