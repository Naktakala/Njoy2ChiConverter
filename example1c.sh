#
# Simple cross-section for graphite using mostly defaults and 
# the xmas 172 group structure with a custom weighting spectrum
#
CWD=$PWD

if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

#================================= Set properties here
neutron_file="n-006_C_000.endf"
sab_file="tsl-graphite.endf"

output_directory="../output/ENDF-B-VII-1/xmas172/"
output_file_prefix="Cnat_graphite"

#================================= Run NJOY
cd njoy_runner || exit

python generate_njoy_mgxs.py \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/$neutron_file \
--path_to_sab=$ENDF_ROOT/thermal_scatt/$sab_file \
--inelastic_thermal_number=229 \
--inelastic_thermal_num_atoms=1 \
--elastic_thermal_number=230 \
--temperature=296.0 \
--neutron_group_structure=22 \
--neutron_weight_function=1 \
--custom_neutron_wt_file="spectrum_file.txt" \
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

cd "$CWD" || exit

#================================= Run converter
cd njoy_processor || exit

python njoy_processor.py \
--path_to_njoy_output=$output_directory/$output_file_prefix.njoy \
--output_file_path=$output_directory/$output_file_prefix.csx \
--plot

cd "$CWD" || exit