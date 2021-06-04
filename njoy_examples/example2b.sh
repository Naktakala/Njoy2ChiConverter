#
# Simple cross-section for O16 using mostly defaults and 
# the lanl30 neutron group structure and the gamma lanl 12 group structure
#
CWD=$PWD

#================================= Set properties here
export ENDF_ROOT=/Users/janv4/Desktop/Projects/ENDF/ENDF-B-VII.1
neutron_file="n-008_O_016.endf"
gamma_file="photoat-008_O_000.endf"
# sab_file="tsl-graphite.endf"

output_directory="../output/ENDF-B-VII-1/xmas172lanl48/"
output_file_prefix="O16"

#================================= Run NJOY
cd njoy_automate2

python generate_njoy_mgxs.py \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
--path_to_gamma_endf=$ENDF_ROOT/photoat/$gamma_file \
--temperature=293.6 \
--neutron_group_structure=22 \
--custom_neutron_wt_file=$output_directory/spectrum_file.txt \
--gamma_group_structure=6 \
--gamma_weight_function=3 \
--output_directory=$output_directory \
--output_filename=$output_file_prefix.njoy

# --path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
# --path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
# --inelastic_thermal_number=229 \
# --inelastic_thermal_num_atoms=1 \
# --elastic_thermal_number=230 \
# --path_to_gamma_endf= \
# --temperature=296.0 \
# --neutron_group_structure=22 \
# --neutron_weight_function=8 \
# --output_directory=$output_directory \
# --output_filename=$output_file_prefix.njoy \
# --gamma_group_structure=0 \
# --gamma_weight_function=2 \
# --custom_neutron_gs_file="" \
# --custom_gamma_gs_file="" \
# --custom_neutron_wt_file="" \
# --custom_gamma_wt_file="" \

cd $CWD

#================================= Run converter
cd njoy_converter2

python njoy_converter2.py \
--path_to_njoy_output=$output_directory/$output_file_prefix.njoy \
--output_file_path=$output_directory/$output_file_prefix.csx 
# --plot

cd $CWD