#! /usr/bin/env bash

isotopes=( "H-1_ZrH" "Zr-nat" "C-nat" "U-235" "U-238" )
temperatures=( "293.6" "400" "600" "800" "1000" )
grp_structs=( "1" "3" "5" "31" "lanl30" "lanl70" "lanl80" "lanl187" "lanl618" "xmas172" )
nummom="7"

echo > generate_xs.log
# Loop over isotopes
for isotope in ${isotopes[*]}
do
	# Loop over temperatures
	for temp in ${temperatures[*]}
	do

		# Get correct room temp for isotope
		if 	[[ $isotope == "C-nat" ]] || \
				[[ $isotope == *"Zr"* ]]  && \
				[[ $temp == "293.6" ]]
		then
			temp="296"
		fi

		# Loop over group structures
		for grp_struct in ${grp_structs[*]}
		do	
			printf "\n*** STARTING %s %sg %sm %sK***\n" $isotope $grp_struct $nummom $temp

			# Get correct template to use
			if	[[ $grp_struct == *"lanl"* ]] || \
				 	[[ $grp_struct == *"xmas"* ]] 
			then
				template="default"
			else
				template="${grp_struct}g"
			fi
			template_path="templates/runNJOY_${template}_template.sh"

			# Isotope parsing
			IFS='_'; read -ra ADDR <<< "$isotope"
			isosub=${ADDR[0]}; molsub=${ADDR[1]}
			isonospacesub=${isosub/-/''} # Take '-' out of isosub	
			nfile=../endf/neutron/"$isonospacesub"_endf.txt # endf file for isotope
			IFS='-'; read -ra ADDR <<< "$isosub"; IFS=' '
			elementsub=${ADDR[0]}
			matnum=`./get_mat_ID.py $nfile` # get the material number

			# Copy template file for modification
			cp $template_path runNJOY.sh

			# Substitutions
			source substitutions/thermal_scatter_sub.sh $isotope runNJOY.sh
			sed -i -e "s/isosub/${isosub}/g" runNJOY.sh
			sed -i -e "s/isonospacesub/${isonospacesub}/g" runNJOY.sh
			sed -i -e "s/elementsub/${elementsub}/g" runNJOY.sh
			sed -i -e "s/matsub/${matnum}/g" runNJOY.sh
			sed -i -e "s/nummoms/${nummom}/g" runNJOY.sh
			sed -i -e "s/temp/${temp}/g" runNJOY.sh
			source substitutions/grp_structure_sub.sh $grp_struct runNJOY.sh
			source substitutions/fission_rxn_sub.sh $nfile runNJOY.sh
			
			# Run NJOY
			echo $isotope >> generate_xs.log
			. ./runNJOY.sh 2>&1 >> generate_xs.log
			rm runNJOY.sh

			# Create output directories
			path="../njoy_xs/${grp_struct}g"
			if [ ! -d ${path} ]; then
				mkdir ${path}
			fi
			if [[ $temp == "293.6" ]] || [[ $temp == "296" ]]
			then
				path="${path}/room"
			else
				path="${path}/${temp}k"
			fi
			if [ ! -d ${path} ]; then
				mkdir ${path}
			fi

			# Create unique name for the isotope/molecule
			if [[ $molsub != "" ]]; then
				filename="${isonospacesub}_${molsub}"
			else
				filename="${isonospacesub}"
			fi
			# Create directories for output
			mv output "${path}/${filename}".txt
			
			printf "*** FINISHED %s %sg %sm %sK ***\n" $isotope $grp_struct $nummom $temp
		done
	done
done
