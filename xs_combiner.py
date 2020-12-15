'''
Cross section parser for analysis of data.

Author: Zachary Hardy
Date: 12/2020
'''

import os
from os.path import isfile
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.arraysetops import isin

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
  E_bndrys[G] = val[G-1][2]

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
def CombineXS(isotopes, grp_struct, temperature):
  assert isinstance(isotopes, list), "'isotopes' must be a list."
  assert all([isinstance(iso, tuple) for iso in isotopes]), \
    "All entries of isotopes must be a tuple."
  assert all([len(iso)==3 for iso in isotopes]), \
    "All entries of isotopes must be a tuple of length 3."
  assert isinstance(grp_struct, str), "'grp_struct' must be a string."
  assert isinstance(temperature, str), "'temperature' must be a string."

  root = 'njoy_xs'
  xs_dir = os.path.join(root, grp_struct, temperature)
  assert os.path.isdir(xs_dir), "Invalid directory path."

  #===== Get the cross sections
  data = []
  for iso,mol,density in isotopes:
    isomol = iso+'_'+mol if mol != '' else iso
    filepath = os.path.join(xs_dir, isomol+'.njoy')
    assert os.path.isfile(filepath), "Invalid filepath."

    #===== Parse NJOY files
    raw_njoy_data = Utils_ReadNJOYOutput.ReadNJOYfile(filepath)
    data_ = Utils_Combiner.BuildCombinedData(raw_njoy_data)
    data += [data_]

  #===== Preliminary calculations
  Nf = 0.0
  for i in range(len(data)):
    if np.sum(data[i]['sigma_f']) > 0.0:
      Nf += isotopes[i][-1]

  #===== Combine data
  xs = {}
  for i in range(len(data)):
    density = isotopes[i][-1]

    #===== Parse group structure
    ParseGroupStructure(data[i]['neutron_gs'], xs)
    xs['neutron_gs'] = data[i]['neutron_gs']

    #===== Parse total cross section
    if 'sigma_t' not in xs:
      xs['sigma_t'] = np.zeros(xs['G'])
    xs['sigma_t'] += density * np.array(data[i]['sigma_t'])

    #===== Parse absorption cross section
    if 'sigma_a' not in xs:
      xs['sigma_a'] = np.zeros(xs['G'])
    xs['sigma_a'] += density * np.array(data[i]['sigma_a'])

    #===== Parse scattering cross section
    if 'sigma_s' not in xs:
      xs['sigma_s'] = np.zeros(xs['G'])
    xs['sigma_s'] += density * np.array(data[i]['sigma_s'])

    #===== Parse prompt fission quantities
    if 'sigma_f' not in xs:
      xs['sigma_f'] = np.zeros(xs['G'])
      xs['nu_sigma_f'] = np.zeros(xs['G'])
      xs['chi_prompt'] = np.zeros(xs['G'])
    nu = np.array(data[i]['nu_total'])
    xs['sigma_f'] += density * np.array(data[i]['sigma_f'])
    xs['nu_sigma_f'] += density * nu * np.array(data[i]['sigma_f'])
    xs['chi_prompt'] += density/Nf * np.array(data[i]['chi_prompt'])

    #===== Parse transfer matrices
    tr_mat = np.array(data[i]['transfer_matrices'])
    if 'transfer_matrices' not in xs:
      xs['transfer_matrices'] = np.zeros(tr_mat.shape)
    xs['transfer_matrices'] += density * tr_mat
  return xs


############################################################
if __name__ == "__main__":

  isotopes = [('U235','',0.01), ('U238','',0.01)]
  grp_struct = 'xmas172g'
  temperature = 'room'
  plot_xs = False
  
  xs = CombineXS(isotopes, grp_struct, temperature)
  Utils_Info.ComputeKinf(xs)

  if plot_xs:
    E_midpts = xs['E_midpts'][::-1]
    sig_t = xs['sigma_t']
    sig_a = xs['sigma_a']
    sig_s = xs['sigma_s']
    sig_f = xs['sigma_f']
    nu_sig_f = xs['nu_sigma_f']

    fig = plt.figure(figsize=(6,6))
    plt.title("UZrH$_{1.6}$ with $S(\\alpha, \\beta)$")
    plt.xlabel("Neutron Energy (eV)")
    plt.ylabel("Cross Section (cm$^{-1}$)")
    plt.semilogx(E_midpts,sig_t,label=r"$\sigma_t$")
    plt.semilogx(E_midpts,sig_a,label=r"$\sigma_a$")
    plt.semilogx(E_midpts,sig_s,label=r"$\sigma_s$")
    plt.semilogx(E_midpts,sig_f,label=r"$\sigma_f$")
    plt.semilogx(E_midpts,nu_sig_f,label=r"$\nu \sigma_f$")
    plt.legend()
    plt.grid(True)
    plt.savefig('/Users/zachhardy/Desktop/withSaB.png')
    plt.show()

