'''
Cross section parser for analysis of data.

Author: Zachary Hardy
Date: 12/2020
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

import njoy_converter.Utils_ReadNJOYOutput as Utils_ReadNJOYOutput
import njoy_converter.Utils_Combiner as Utils_Combiner
import njoy_converter.Utils_ChiWriter as Utils_ChiWriter
import njoy_converter.Utils_Info as Utils_Info

############################################################
def ParseGroupStructure(val, xs):
  G         = len(val)
  E_bndrys  = np.zeros(G+1)
  E_midpts  = np.zeros(G)
  dE        = np.zeros(G)
  for i in range(G):
    E_bndrys[i]  = val[i][1]
    E_midpts[i]  = 0.5*(val[i][2] + val[i][1])
    dE[i]        = val[i][2] - val[i][1]
  E_bndrys[G-1] = val[G-1][2]

  if 'G' not in xs:
    xs['G']        = G
    xs['E_bndrys'] = E_bndrys
    xs['E_midpts'] = E_midpts
    xs['dE']       = dE
  else:
    try:
      if G != xs['G']:
        msg = "Number of groups do not match."
        raise ValueError(msg)
      for g in range(G+1):
        if E_bndrys[g] != xs['E_bndrys'][g]:
          msg = "Group boundaries do not match."
          raise ValueError(msg)
    except ValueError as err:
      print(err.args[0])
      sys.exit(0)

############################################################
if __name__ == "__main__":

  rxns = ['sigma_t', 'sigma_a', 'sigma_s', 'sigma_f', 
          'sigma_s_el', 'sigma_s_inel', 'sigma_s_nxn', 
          'sigma_s_freegas', 'sigma_s_sab', 'sigma_s_sab_el', 
          'sigma_s_sab_inel']

  isotopes = [('U235','',0.001078), ('U238','',0.004259), 
              ('Zrnat','ZrH',0.031854), ('H1','',0.050966)]

  root = 'njoy_xs'
  grp_struct = 'xmas172g'
  temperature = 'room'
  
  
  xs_dir = os.path.join(root, grp_struct, temperature)
  assert os.path.isdir(xs_dir), "Invalid directory."

  xs = {}
  for iso,mol,density in isotopes:
    isomol = iso+'_'+mol if mol != '' else iso
    filepath = os.path.join(xs_dir, isomol+'.njoy')
    assert os.path.isfile(filepath), "Invalid filepath."

    raw_njoy_data = Utils_ReadNJOYOutput.ReadNJOYfile(filepath)
    data = Utils_Combiner.BuildCombinedData(raw_njoy_data)

    ParseGroupStructure(data['neutron_gs'], xs)

    #===== Iterate over items in dictionary
    for key,val in data.items():

      #===== Parse reactions, scale by densities
      if key in rxns:
        if key not in xs:
          xs[key] = np.zeros(xs['G'])
        xs[key] += density * np.array(val)
  

  E_midpts = xs['E_midpts'][::-1]
  sig_t = xs['sigma_t']
  sig_a = xs['sigma_a']
  sig_s = xs['sigma_s']

  fig = plt.figure(figsize=(6,6))
  plt.title("UZrH$_{1.6}$ without $S(\\alpha, \\beta)$")
  plt.xlabel("Neutron Energy (eV)")
  plt.ylabel("Cross Section (cm$^{-1}$)")
  plt.semilogx(E_midpts,sig_t,label=r"$\sigma_t$")
  plt.semilogx(E_midpts,sig_a,label=r"$\sigma_a$")
  plt.semilogx(E_midpts,sig_s,label=r"$\sigma_s$")
  plt.legend()
  plt.grid(True)
  plt.show()