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