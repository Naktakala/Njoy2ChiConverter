# !/usr/bin/env bash

if grep -Eq "fission" $1 
then
	cat > tmp.txt <<-END
		3 452 \'total nubar\'/\\
		3 456 \'prompt nubar\'/\\
		3 455 \'delayed nubar\'/\\
		5 18 \'prompt chi\'/\\
		5 455 \'delayed chi\'/\\	
	END
	RXN=$( cat tmp.txt )
	sed -i -e "s|fission_rxns|${RXN}|g" $2
	rm tmp.txt
else
	sed -i -e "/fission_rxns/d" $2
fi

