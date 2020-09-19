#! /usr/bin/env bash
echo 'Getting ENDF input tapes and NJOY executable'
ln -fs ../endf/neutron/isonospacesub_endf.txt tape21
ln -fs ../endf/neutron_thermal/isonospacesub_endf.txt tape41                            
echo 'Running NJOY'
cat>input_isosub <<EOF
moder
21 -22/
reconr
-22 -23/
'pendf tape for isosub from endf/b-vii1'/
matsub 3/
0.001/
'isosub from endf/b-vii1 tape'/
'processed by the njoy nuclear data processing system'/
'see original endf/b-vii1 tape for details of evaluation'/
0/
broadr
-22 -23 -24/
matsub 1 0 0 0/     
0.001/
temp/
0/
unresr
-22 -24 -25/
matsub 1 1 0/       
temp/
0.0/
0/
heatr
-22 -25 -90/
matsub 0/                                                                
thermr
0 -90 -26/
0 matsub 16 1 1 0 0 1 221 1/
temp/
0.05 10/
s_alpha_beta_input/
groupr
-22 -26 0 -30/
matsub gs 0 5 nummoms 1 1 1 1/                   GROUP STRUCTURE IS 2ND ENTRY, LEGENDRE MOMENTS IS 5TH
'isosub'/
temp/
0.0
3/
3 259 'inverse velocity'/
3 221 'therm scat'/
fission_rxns/
6/
6 221 'transfer matrix'/
s_alpha_beta_rxns/
0/
0/
moder
-30 31/
stop
EOF
njoy<input_isosub
echo 'Cleaning up'
rm input_isosub
rm -f tape*
