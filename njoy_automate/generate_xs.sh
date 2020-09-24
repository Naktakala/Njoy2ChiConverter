#! /usr/bin/env bash

isotopes=( "H-1" "O-16" "C-nat" "U-235" "U-238" "Pu-239" )
temperatures=( "293.6" "400" "500" "600" "800" )
grp_structs=( "1" "3" "5" "31" "lanl30" "lanl70" "lanl80" "lanl187" "lanl618" "xmas172")
nummoms=( "7" )

echo > generate_xs.log
for isotope in ${isotopes[*]}
do
	for temp in ${temperatures[*]}
	do
		for grp_struct in ${grp_structs[*]}
		do	
			for nummom in ${nummoms[*]}
			do
				printf "\n*** STARTING %s %sg %sm %sK***\n" $isotope $grp_struct $nummom $temp

				# Determine the correct template to use
				if [[ $grp_struct == *"lanl"* ]] || [[ $grp_struct == *"xmas"* ]]
				then
					template="default"
				else
					template="${grp_struct}g"
				fi
				template_path="templates/runNJOY_${template}_template.sh"

				# Correct temp for S(\alpha,\beta) in graphite
				if [[ $isotope == "C-nat" ]] || [[ $isotope == "Zr-nat" ]]
				then
					if [[ $temp == "293.6" ]]
					then
						temp="296"
					fi
				fi

				# Isotope parsing
				isosub=$isotope
				isonospacesub=${isosub/-/''} # Take '-' out of isosub			
				nfile=../endf/neutron/"$isonospacesub"_endf.txt # endf file for isotope
				IFS='-'
				read -ra ADDR <<< "$isotope"  # Split isotope name on '-'
				IFS=' '
				elementsub=${ADDR[0]}
				matnum=`./get_mat_ID.py $nfile` # get the material number

				# Copy template file for modification
				cp $template_path runNJOY.sh 

				# Substitutions
				source substitutions/thermal_scatter_sub.sh $isosub runNJOY.sh
				sed -i -e "s/isosub/${isosub}/g" runNJOY.sh
				sed -i -e "s/isonospacesub/${isonospacesub}/g" runNJOY.sh
				sed -i -e "s/elementsub/${elementsub}/g" runNJOY.sh
				sed -i -e "s/matsub/${matnum}/g" runNJOY.sh
				sed -i -e "s/nummoms/${nummom}/g" runNJOY.sh
				sed -i -e "s/temp/${temp}/g" runNJOY.sh
				source substitutions/grp_structure_sub.sh $grp_struct runNJOY.sh
				source substitutions/fission_rxn_sub.sh $nfile runNJOY.sh
				if [[ $isosub == "H-1" ]]
				then
					if [[ $temp == "293.6" ]]
					then
						sed -i -e "s/alt/296/g" runNJOY.sh
					else
						sed -i -e "s/alt/${temp}/g" runNJOY.sh
					fi
				fi
				
				# Run NJOY
				echo $isotope >> generate_xs.log
				. ./runNJOY.sh 2>&1 >> generate_xs.log
				rm runNJOY.sh

				# Write to unique filename and move to correct dir
				if [[ $temp == "293.6" ]] || [[ $temp == "296" ]]
				then
					description="${grp_struct}g${nummom}m_room"
				else
					description="${grp_struct}g${nummom}m${temp}k"
				fi
				mv output ../njoy_xs/"$isonospacesub"_"$description".txt
				
				printf "*** FINISHED %s %sg %sm %sK ***\n" $isotope $grp_struct $nummom $temp
			done
		done
	done
done
