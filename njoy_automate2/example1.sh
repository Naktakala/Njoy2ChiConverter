export ENDF_ROOT=/Users/janv4/Desktop/Projects/ENDF/ENDF-B-VII.1

python generate_njoy_mgxs.py \
--path_to_neutron_endf=$ENDF_ROOT/neutrons/n-006_C_000.endf \
--path_to_sab=$ENDF_ROOT/thermal_scatt/tsl-graphite.endf \
--temperature=296.0 \
--output_directory="../njoy_xs/epri-cpm69/" \
--output_filename="Cnat_graphite.txt" 
#--path_to_gamma_endf=$ENDF_ROOT/gammas/g-006_C_012.endf \
# --neutron_group_structure=9 \
# --neutron_weight_function=8 \
# --gamma_group_structure=0 
# --gamma_weight_function=2
# --custom_neutron_gs_file=somefile
# --custom_gamma_gs_file=somefile
# --custom_neutron_wt_file=somefile
# --custom_gamma_wt_file=somefile
