#
# Template to run NJOY and then convert to Chi 
# cross-sections.
#
CWD=$PWD

#================================= Set properties here
export NEUTRON_ROOT=~/endf/neutron
export SAB_ROOT=~/endf/neutron_thermal
neutron_file="isotope.endf"
sab_file="isomol.endf"
njoy_directory="${CWD}/njoy_xs/gsname/Tname/"
chi_directory="${CWD}/chi_xs/gsname/Tname/"
output_file_prefix="outfile"

#================================= Run converter
cd njoy_converter

python njoy_converter.py \
--path_to_njoy_output="${njoy_directory}${output_file_prefix}.njoy" \
--output_file_path="${chi_directory}${output_file_prefix}.csx"

cd $CWD