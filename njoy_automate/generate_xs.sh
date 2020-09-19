#! /usr/bin/env bash

function default_gs_mapping(){
	if [[ $1 == *"lanl"* ]] || [[ $1 == *"xmas"* ]]
	then
		if [[ $1 == "lanl30" ]]
		then
			sed -i -e 's/gs/3/g' $2
		elif [[ $1 == "lanl70" ]]
		then
			sed -i -e 's/gs/11/g' $2
		elif [[ $1 == "lanl80" ]]
		then
			sed -i -e 's/gs/13/g' $2
		elif [[ $1 == "lanl187" ]]
		then
			sed -i -e 's/gs/10/g' $2
		elif [[ $1 == "lanl618" ]]
		then
			sed -i -e 's/gs/34/g' $2
		elif [[ $1 == "xmas172" ]]
		then
			sed -i -e 's/gs/22/g' $2
		elif [[ $1 == "xmasnea" ]]
		then
			sed -i -e 's/gs/18/g' $2
		else
			echo "Unrecognized default group structure."
		fi
	fi
}

function fission_reaction_sub(){
	if grep -Eq "fission" $1 
	then
		fission_rxns=$( cat substitutions/fission_rxns.txt )
		sed -i -e "s/fission_rxns/${fission_rxns}/g" $2
	else
		sed -i -e "/fission_rxns/d" $2
	fi
	
}

function s_alpha_beta_sub(){
	if [[ $1 == "H-1" ]]
	then
		sab=$( cat substitutions/s_alpha_beta.txt )
		sed -i -e "s/s_alpha_beta_input/${sab}/g" $2
		sed -i -e "s/s_alpha_beta_rxns/6 222 \'h2o_therm matrix\'/g" $2
		sed -i -e "s/matid/1/g" $2
		sed -i -e "s/icoh/1/g" $2
		sed -i -e "s/natom/2/g" $2
		sed -i -e "s/mtref/222/g" $2
		sed -i -e "s/-22 -26/-22 -27/g" $2
	elif [[ $1 == "C-nat" ]]
	then
		sab=$( cat substitutions/s_alpha_beta.txt )
		rxns=$( cat substitutions/graphite_rxns.txt )
		sed -i -e "s/s_alpha_beta_input/${sab}/g" $2
		sed -i -e "s/s_alpha_beta_rxns/${rxns}/g" $2
		sed -i -e "s/matid/31/g" $2
		sed -i -e "s/icoh/1/g" $2
		sed -i -e "s/natom/1/g" $2
		sed -i -e "s/mtref/229/g" $2
		sed -i -e "s/-22 -26/-22 -27/g" $2
	else
		sed -i -e "/neutron_thermal/d" $2
		sed -i -e "/s_alpha_beta_input/d" $2
		sed -i -e "/s_alpha_beta_rxns/d" $2
	fi
}

# temperatures=( "296" "400" "500" "600" "800" )
temperatures=( "296" )

isotopes=( "C-nat" )
# isotopes=( "H-1" "O-16" "C-nat" "U-235" "U-238" "Pu-239" )

# grp_structs=( "1" "3" "5" "31" "lanl30" "lanl70" "lanl80" "lanl187" "lanl618" "xmas172")
# grp_structs=( "lanl187" "lanl618" "xmas172" )
grp_structs=( "xmas172" )

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
				s_alpha_beta_sub $isosub runNJOY.sh
				sed -i -e "s/isosub/${isosub}/g" runNJOY.sh
				sed -i -e "s/isonospacesub/${isonospacesub}/g" runNJOY.sh
				sed -i -e "s/elementsub/${elementsub}/g" runNJOY.sh
				sed -i -e "s/matsub/${matnum}/g" runNJOY.sh
				sed -i -e "s/nummoms/${nummom}/g" runNJOY.sh
				sed -i -e "s/temp/${temp}/g" runNJOY.sh
				default_gs_mapping ${grp_struct} runNJOY.sh
				fission_reaction_sub $nfile runNJOY.sh
				
				# Run NJOY
				echo $isotope >> generate_xs.log
				. ./runNJOY.sh 2>&1 >> generate_xs.log
				rm runNJOY.sh

				# Write to unique filename and move to correct dir
				if [[ $temp == "293.6" ]]
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
