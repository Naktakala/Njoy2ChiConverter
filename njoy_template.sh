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
