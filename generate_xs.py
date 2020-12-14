import os
import sys
import functions as fcts

opt = 0
neutron_root = '/Users/zachhardy/Projects/endf/neutron'
sab_root = '/Users/zachhardy/Projects/endf/neutron_thermal'

#===== Parse any command line arguments
if len(sys.argv) > 1:
  opt = int(sys.argv[1])
if len(sys.argv) > 2:
  neutron_root = sys.argv[2]
if len(sys.argv) > 3:
  sab_root = sys.argv[3]

#===== Check inputs
assert opt in [0,1,2], "Invalid opt. Must be 0, 1, or 2."
assert os.path.isdir(neutron_root), "Invalid neutron root."
assert os.path.isdir(sab_root), "Invalid S(alpha,beta) root."

isotopes = [('H1','ZrH'), ('H1','H2O'), ('H1',''), ('O16',''), 
            ('Zrnat','ZrH'), ('U235',''), ('U238','')]
grp_structs = ['1g','3g','5g','6g','31g','lanl30g','lanl70g',
               'lanl80g','lanl187g','lanl618g','xmas172g']
temperatures = [293.6,400,500,600,800]
plot = False

#===== Loop over isotopes
for isotope,molecule in isotopes:
  outp = fcts.ThermScatInfo(isotope, molecule)
  isomol, sab_inel, sab_el, n_atoms = outp

  #===== Loop over group structures
  for gs in grp_structs:
    gsnum, gsfile = fcts.ParseGroupStructure(gs)

    #===== Loop over temperatures
    for temp in temperatures:
      Tname = fcts.ParseTemperature(temp)
      if molecule in ['graphite','ZrH']:
        if temp == 293.6: temp = 296.0

      msg  = "*"*40 + "\n"
      msg += "STARTING: {:^10}{:^10}{:^10}\n"
      msg += "*"*40 + "\n"
      print(msg.format(isotope+" "+molecule,gs,Tname))

      with open('template.sh','rt') as fin:
        with open('run.sh', 'wt') as fout:
          lines = fin.readlines()
          for l, line in enumerate(lines):
            if line[0] != "#":
              #===== Replace neutron root
              if 'neutron_root' in line:
                line = line.replace('neutron_root',neutron_root)
              #===== Replace S(alpha,beta) root
              if 'sab_root' in line:
                line = line.replace('sab_root',sab_root)
              #===== Replace isotope
              if 'isotope' in line:
                line = line.replace('isotope',isotope)
              #===== Replace output filename
              if 'outfile' in line:
                if isomol is None: line = line.replace('outfile',isotope)
                else: line = line.replace('outfile',isomol)
              if 'path_to_sab' in line:
                if isomol is None: line = ''
              #===== Replace group structure name
              if 'gsname' in line:
                line = line.replace('gsname', str(gs))
              #===== Replace temperature name
              if 'Tname' in line:
                line = line.replace('Tname',Tname)
              #===== Replace group structure number
              if 'gsnum' in line:
                line = line.replace('gsnum',str(gsnum))
              #===== Replace group structure file
              if 'gsfile' in line:
                if gsnum != 1: line = ''
                else: line = line.replace('gsfile',gsfile)
              #===== Replace temperature
              if 'Tval' in line:
                line = line.replace('Tval',str(temp))
              #===== Replace S(a,b) file
              if 'isomol' in line:
                if isomol is None: line = ''
                else: line = line.replace('isomol',isomol)
              #===== Replace S(a,b) inelastic MT number
              if 'sab_inel' in line:
                if sab_inel is None: line = ''
                else: line = line.replace('sab_inel',str(sab_inel))
              #===== Replace S(a,b) elastic MT number
              if 'sab_el' in line:
                if sab_el is None: line = ''
                else: line = line.replace('sab_el',str(sab_el))
              #===== Replace n_atoms for S(a,b)
              if 'n_atoms' in line:
                if n_atoms is None: line = ''
                else: line = line.replace('n_atoms',str(n_atoms))
              #===== Replace plot
              if 'plot' in line:
                if not plot: line = ''
            fout.write(line)
      os.system(". ./run.sh "+str(opt))
      os.system("rm run.sh")

      msg  = "\n" + "*"*40 + "\n"
      msg += "FINISHED: {:^10}{:^10}{:^10}\n"
      msg += "*"*40 + "\n"
      print(msg.format(isotope+" "+molecule,gs,Tname))
      