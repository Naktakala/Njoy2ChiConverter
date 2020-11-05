#
# Template to run NJOY and then convert to Chi 
# cross-sections.
#
CWD=$PWD

#================================= Set properties here
export NEUTRON_ROOT=neutron_root
export SAB_ROOT=sab_root
neutron_file="isotope.endf"
sab_file="isomol.endf"
njoy_directory="${CWD}/njoy_xs/gsname/Tname"
chi_directory="${CWD}/chi_xs/gsname/Tname"
output_file_prefix="outfile"

if [[ $1 == '0' ]] || [[ $1 == '1' ]]
then
  #================================= Run NJOY
  cd njoy_automate

  python generate_njoy_mgxs.py \
  --path_to_neutron_endf=$NEUTRON_ROOT/$neutron_file \
  --path_to_sab=$SAB_ROOT/$sab_file \
  --neutron_group_structure=gsnum \
  --neutron_weight_function=5 \
  --custom_neutron_gs_file=gsfile \
  --custom_neutron_wt_file="" \
  --temperature=Tval \
  --inelastic_thermal_number=sab_inel \
  --elastic_thermal_number=sab_el \
  --inelastic_thermal_num_atoms=n_atoms \
  --output_directory=$njoy_directory \
  --output_filename=$output_file_prefix.njoy

  cd $CWD
fi

#================================= Run converter
if [[ $1 == '0' ]] || [[ $1 == '2' ]]
then
  cd njoy_converter

  python njoy_converter.py \
  --path_to_njoy_output="${njoy_directory}/${output_file_prefix}.njoy" \
  --output_file_path="${chi_directory}/${output_file_prefix}.csx" \
  --plot

  cd $CWD
fi