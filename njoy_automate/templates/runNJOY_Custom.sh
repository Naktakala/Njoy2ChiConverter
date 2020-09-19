#! /usr/bin/env bash
echo 'Getting ENDF input tapes and NJOY executable'
ln -fs ../endf/neutron/isonospacesub_endf.txt tape20
ln -fs ../endf/neutron_thermal/isonospacesub_endf8.txt tape26                            
echo 'Running NJOY'
cat>input_isosub <<EOF
moder
20 -21/
reconr
-21 -22/
'pendf tape for isosub from endf/b-vii1'/
matsub 3/
0.005/
'isosub from endf/b-vii1 tape'/
'processed by the njoy nuclear data processing system'/
'see original endf/b-vii1 tape for details of evaluation'/
0/
broadr
-21 -22 -23/
matsub 1 0 0 0/     
0.005/
temp/
0/
heatr
-21 -23 -22/
matsub 0/                                                        
thermr
0 -22 -24/
0 matsub 8 1 1 0 0 1 221 2/
temp/
0.01 100.0/
thermr
26 -24 -23
matid matsub 8 1 2 1 0 1 229 2
temp/
0.05 100.0/
groupr
-21 -23 0 -24/
matsub gs 0 2 nummoms 1 1 1 1/                   GROUP STRUCTURE IS 2ND ENTRY, LEGENDRE MOMENTS IS 5TH
'isosub'/
temp/
0.0
3/
3 259 'inverse velocity'/
3 221 'therm scat'/
6/
s_alpha_beta_rxns/
0/
0/
moder
-23 25/
stop
EOF
/Users/janv4/Desktop/Projects/NJOY21/NJOY21/bin/njoy21<input_isosub -o output
echo 'Cleaning up'
rm input_isosub
rm -f tape*
