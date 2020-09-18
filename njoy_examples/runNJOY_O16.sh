#! /usr/bin/env bash
echo 'Getting ENDF input tapes and NJOY executable'
ln -fs endf/neutron/O16_endf.txt tape21
echo 'Running NJOY'
cat>input_O-16 <<EOF
moder
21 -22
reconr
-22 -23
'pendf tape for O-16 from endf/b-vii1'/
825 3/ MATERIAL NUMBER
0.001/
'O-16 from endf/b-vii1 tape'/
'processed by the njoy nuclear data processing system'/
'see original endf/b-vii1 tape for details of evaluation'/
0/
broadr
-22 -23 -24
825 1 0 0 0/ MATERIAL NUMBER
0.001/
293.6/ TEMPERATURE
0/
unresr
-22 -24 -25
825 1 1 0/ MATERIAL NUMBER        
293.6/ TEMPERATURE
0.0
0/
heatr
-22 -25 -90/
825 0/ MATERIAL NUMBER                                                
thermr
0 -90 -26
0 825 16 1 1 0 0 1 221 1/ MATERIAL NUMBER
293.6/ TEMPERATURE
0.001 100
groupr
-22 -26 0 0
825 1 3 2 1 1 1 1 1/ MATERIAL NUMBER, GROUP STRUCTURE IS 2ND ENTRY, LEGENDRE MOMENTS IS 5TH
'O-16'/
293.6/ TEMPERATURE
0.0
31/ 31 groups, followed by 32 group boundaries (eV)
1.0000e-05 4.6589e+05 9.3178e+05 1.3977e+06 1.8636e+06 2.3295e+06
2.7953e+06 3.2612e+06 3.7271e+06 4.1930e+06 4.6589e+06 5.1248e+06
5.5907e+06 6.0566e+06 6.5225e+06 6.9884e+06 7.4775e+06 7.9434e+06
8.4093e+06 8.8752e+06 9.3411e+06 9.8070e+06 1.0273e+07 1.0739e+07
1.1205e+07 1.1671e+07 1.2136e+07 1.2602e+07 1.3068e+07 1.3534e+07
1.3977e+07 1.4000e+07 / energy group bounds
3/
6/
16/
0/
0/
stop
EOF
njoy<input_O-16
echo 'Cleaning up'
rm input_O-16
rm -f tape*
mv output output_O-16.txt
