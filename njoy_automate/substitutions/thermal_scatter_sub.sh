#! /usr/bin/env bash

# thermr called for H i	n H2O
if [[ $1 == "H-1_H2O" ]] 
then
	# Load up the thermal data
	cat > tmp.txt <<-END
		ln -fs ../endf/neutron_thermal/H_in_H2O_endf.txt tape41
	END
	DATA=$( cat tmp.txt )
	sed -i -e "s|thermal_data|${DATA}|g" $2

	# Construct the THERMR input
	cat > tmp.txt <<-END
		thermr\\
		41 -26 -27/\\
		1 matsub 16 1 2 1 0 2 222 1/\\
		temp/\\
		0.001 10.0
	END
	THERMR=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_thermr|${THERMR}|g" $2

	# Construct the S(alpha,beta) reactions
	cat > tmp.txt <<- END
		6 222 \'H_in_H2O_inelastic matrix\'
	END
	RXN=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_rxns|${RXN}|g" $2
	sed -i -e "s/-22 -26/-22 -27/g" $2
	rm tmp.txt

# thermr called for H in ZrH
elif [[ $1 == "H-1_ZrH" ]]
then
	# Load up the thermal data
	cat > tmp.txt <<-END
		ln -fs ../endf/neutron_thermal/H_in_ZrH_endf.txt tape41
	END
	DATA=$( cat tmp.txt )
	sed -i -e "s|thermal_data|${DATA}|g" $2

	# Construct the THERMR input
	cat > tmp.txt <<-END
		thermr\\
		41 -26 -27/\\
		7 matsub 16 1 2 1 0 1 225 1/\\
		temp/\\
		0.001 10.0
	END
	THERMR=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_thermr|${THERMR}|g" $2

	# Construct the S(alpha,beta) reactions
	cat > tmp.txt <<- END
		6 225 \'H_in_ZrH_inelastic matrix\'/\\
		6 226 \'H_in_ZrH_elastic matrix\'
	END
	RXN=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_rxns|${RXN}|g" $2
	sed -i -e "s/-22 -26/-22 -27/g" $2
	rm tmp.txt

# thermr called for graphite
elif [[ $1 == "C-nat" ]] 
then
	# Load up the thermal data
	cat > tmp.txt <<-END
		ln -fs ../endf/neutron_thermal/graphite_endf.txt tape41
	END
	DATA=$( cat tmp.txt )
	sed -i -e "s|thermal_data|${DATA}|g" $2

	# Construct the THERMR input
	cat > tmp.txt <<-END
		thermr\\
		41 -26 -27/\\
		31 matsub 16 1 2 1 0 1 229 1/\\
		temp/\\
		0.001 10.0
	END
	THERMR=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_thermr|${THERMR}|g" $2
	
	# Construct the S(alpha,beta) reactions
	cat > tmp.txt <<-END
		6 229 \'graphite_inelastic matrix\'/\\
		6 230 \'graphite_elastic matrix\'
	END
	RXN=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_rxns|${RXN}|g" $2
	sed -i -e "s/-22 -26/-22 -27/g" $2
	rm tmp.txt
	
# thermr called for ZrH
elif [[ $1 == "Zr-nat" ]] 
then
	# Load up the thermal data
	cat > tmp.txt <<-END
		ln -fs ../endf/neutron_thermal/Zr_in_ZrH_endf.txt tape41
	END
	DATA=$( cat tmp.txt )
	sed -i -e "s|thermal_data|${DATA}|g" $2

	# Construct the THERMR input
	cat > tmp.txt <<-END
		thermr\\
		41 -26 -27/\\
		58 matsub 16 1 2 1 0 1 235 1/\\
		temp/\\
		0.001 10.0
	END
	THERMR=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_thermr|${THERMR}|g" $2

	# Construct the S(alpha,beta) reactions
	cat > tmp.txt <<-END
		6 235 \'Zr_in_ZrH_inelastic matrix\'/\\
		6 236 \'Zr_in_ZrH_elastic matrix\'
	END
	RXN=$( cat tmp.txt )
	sed -i -e "s|s_alpha_beta_rxns|${RXN}|g" $2
	sed -i -e "s/-22 -26/-22 -27/g" $2
	rm tmp.txt

# delete substitution lines otherwise
else
	sed -i -e "/thermal_data/d" $2
	sed -i -e "/s_alpha_beta_thermr/d" $2
	sed -i -e "/s_alpha_beta_rxns/d" $2
fi
