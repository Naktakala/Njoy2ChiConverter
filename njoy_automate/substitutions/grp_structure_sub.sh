#! /usr/bin/env bash

if [[ $1 == *"lanl"* ]] || [[ $1 == *"xmas" ]] 
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
