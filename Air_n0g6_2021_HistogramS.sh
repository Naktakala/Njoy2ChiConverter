#
# Simple cross-section for N14 using manual energy groups
#
CWD=$PWD

if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

#================================= Set properties here
neutron_file="n-007_N_014.endf"
gamma_file="photoat-007_N_000.endf"
# sab_file="tsl-graphite.endf"

output_directory="../output/ENDF-B-VII-1/LANL0_LANL48_Air/"
output_file_prefix="Air_n0g48_njoy2021"
source_description="Gamma,Histogram"
#================================= Get the MCNP output files
mcnp_filename="DTRA_Air_n0_g48_100BinHistoS.out"
mcnp_directory="../MCNP/Air/Gamma_only_PointS/"

#================================ Get the arguments for composite
composite_chi_output_prefix="N14_n0g48.csx,O16_n0g48.csx"
#composite_chi_output_prefix="Chi_n0g48/N14_n0g48.csx,Chi_n0g48/O16_n0g48.csx"
atom_fraction="0.79,0.21"
atom_density="5.11533E-05,5.11533E-05"

#================================= Run converter
cd chitech_composite_processor || exit

python3 chitech_processor.py \
--output_path=$output_directory \
--chixs_filename=$output_file_prefix.csx \
--composite_chi_output_filename=$composite_chi_output_prefix \
--atomic_fraction=$atom_fraction \
--atomic_density=$atom_density \
--source_term=$source_description \
--plot \
--mcnp \
--mcnp_path=$mcnp_directory \
--mcnp_filename=$mcnp_filename \

# --output_path=$output_directory \
# --chixs_filename=$output_file_prefix.csx \
# --composite_chi_output_filename=$composite_chi_output_prefix \
# --composite_atomic_fraction=$atom_fraction \
# --source_term=$source_description \
# --molar_mass=$molar_mass \
# --mcnp \
# --mcnp_path=$mcnp_directory \
# --mcnp_filename=$mcnp_filename \
# --more_chitech \
# --more_composite_chi_output_filename=$more_composite_chi_output_prefix \
# --plot \

cd "$CWD" || exit
