#! /usr/bin/env bash
echo 'Getting ENDF input tapes and NJOY executable'
ln -fs ../endf/neutron/isonospacesub_endf.txt tape21
thermal_data            
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
0.001 10.0/
s_alpha_beta_thermr/
groupr
-22 -26 0 -30/
matsub 1 0 5 nummoms 1 1 1 1/                   GROUP STRUCTURE IS 2ND ENTRY, LEGENDRE MOMENTS IS 5TH
'isosub'/
temp/
0.0/
31/                                       NUMBER OF GROUPS (GROUP STRUCTURE = 1 = READ-IN, ONLY)
1.0000e-05 4.6589e+05 9.3178e+05 1.3977e+06 1.8636e+06 2.3295e+06
2.7953e+06 3.2612e+06 3.7271e+06 4.1930e+06 4.6589e+06 5.1248e+06
5.5907e+06 6.0566e+06 6.5225e+06 6.9884e+06 7.4775e+06 7.9434e+06
8.4093e+06 8.8752e+06 9.3411e+06 9.8070e+06 1.0273e+07 1.0739e+07
1.1205e+07 1.1671e+07 1.2136e+07 1.2602e+07 1.3068e+07 1.3534e+07
1.3999e+07 1.4000e+07/                    GROUP BOUNDARIES IN EV (GROUP STRUCTURE = 1 = READ-IN, ONLY)
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
