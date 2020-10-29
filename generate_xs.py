import os
import sys

#=================================================
def ParseGroupStructure(gs):
  #========== Check group structure
  valid_gs = [
    '1g', '3g', '5g', '6g', '31g',
    'lanl30g',  'lanl70g', 'lanl80g',
    'lanl187g', 'lanl618g', 'xmas172g'
  ]
  if gs not in valid_gs:
    raise ValueError("Invalid group structure.")

  #========== Return grp stucture info
  if gs in ['1g', '3g', '5g', '6g', '31g']:
    return 1, '../grp_structs/custom'+gs+'.txt'
  elif 'lanl30g' in gs:
    return 3, ''
  elif 'lanl70g' in gs:
    return 11, ''
  elif 'lanl80g' in gs:
    return 13, ''
  elif 'lanl187g' in gs:
    return 10, ''
  elif 'lanl618g' in gs:
    return 34, ''
  elif 'xmas172g' in gs:
    return 22, ''

#=================================================
def ParseTemperature(temp):
  if temp in [296.0,293.6]: return 'room'
  else: return str(temp)+'k'

#=================================================
def ThermScatInfo(isotope, molecule):
  if molecule == '':
    return None, None, None, None
  elif isotope == 'H1' and molecule == 'H2O':
    return 'H1_H2O', 222, None, 2
  elif isotope == 'H1' and molecule == 'ZrH':
    return 'H1_ZrH', 225, 226, 1
  elif isotope == 'Zrnat' and molecule == 'ZrH':
    return 'Zrnat_ZrH', 235, 236, 1
  elif isotope == 'Cnat' and molecule == 'graphite':
    return 'Cnat_graphite', 229, 230, 1

#=================================================
isotopes = [('H1','ZrH'),('Zrnat','ZrH'),('U235',''),('U238','')]
grp_structs = ['1g','3g','5g','6g','31g','lanl70g','lanl80g','lanl618g']
temperatures = [293.6,400,500,600,800,1000]

if len(sys.argv) > 1:
  if int(sys.argv[1]) == 0:
    template = "njoy2chi_template.sh"
  elif int(sys.argv[1]) == 1:
    template = "njoy_template.sh"
  elif int(sys.argv[1]) == 2:
    template = "chi_template.sh"
else:
  template = "njoy2chi_template.sh"



#===== Loop over isotopes
for isotope,molecule in isotopes:
  outp = ThermScatInfo(isotope, molecule)
  isomol, sab_inel, sab_el, n_atoms = outp

  #===== Loop over group structures
  for gs in grp_structs:
    gsnum, gsfile = ParseGroupStructure(gs)

    #===== Loop over temperatures
    for temp in temperatures:
      Tname = ParseTemperature(temp)
      if molecule in ['graphite','ZrH']:
        if temp == 293.6: temp = 296.0

      msg  = "*"*40 + "\n"
      msg += "STARTING: {:^10}{:^10}{:^10}\n"
      msg += "*"*40 + "\n"
      print(msg.format(isotope+" "+molecule,gs,Tname))

      with open(template,'rt') as fin:
        with open('run.sh', 'wt') as fout:
          lines = fin.readlines()
          for l, line in enumerate(lines):
            if line[0] != "#":
              #===== Replace isotope
              if 'isotope' in line:
                line = line.replace('isotope',isotope)
              #===== Replace output filename
              if 'outfile' in line:
                if isomol is None: line = line.replace('outfile',isotope)
                else: line = line.replace('outfile',isomol)
              if 'path_to_sab' in line:
                if isomol is None:
                  line = ''
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
            fout.write(line)
      os.system(". ./run.sh")
      os.system("rm run.sh")

      msg  = "\n" + "*"*40 + "\n"
      msg += "FINISHED: {:^10}{:^10}{:^10}\n"
      msg += "*"*40 + "\n"
      print(msg.format(isotope+" "+molecule,gs,Tname))
      