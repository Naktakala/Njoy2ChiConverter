#
# Simple cross section for air using manual energy groups
#
CWD=$PWD

if [[ -z "${ENDF_ROOT}" ]]; then
  echo "ENDF_ROOT not set. Please set or specify custom paths."
  exit
fi

output_directory="../output/testing/LANL30_LANL12/"
output_file_prefix="Air_n30g12"
#================================ Get the arguments for composite
chixs_filename_list="N14_n30g12.csx,O16_n30g12.csx"
atom_density="1,0"

#================================= Run converter
cd chitech_composite_processor || exit

python3 chitech_processor.py \
--output_path=$output_directory \
--chixs_filename=$output_file_prefix.csx \
--chixs_filename_list=$chixs_filename_list \
--atomic_density=$atom_density \

cd "$CWD" || exit
