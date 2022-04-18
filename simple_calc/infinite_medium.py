# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 00:55:42 2022

@author: jean.ragusa
"""
import numpy as np

import sys, os
sys.path.insert(0, os.path.realpath('../chitech_composite_processor'))
import Utils_ChiTechCombiner
import Utils_Info
import Utils_NjoySpectrumPlotter
import matplotlib.pyplot as plt
plt.close('all')

chixs_fullpath = []
ID = 'Ar'
if ID =='Ar':
    chixs_fullpath.append('../output/testing/XMAS_172/Ar36_n172.csx')
if ID =='Al':
    chixs_fullpath.append('../output/testing/XMAS_172/Al27_n172.csx')
if ID =='N':
    chixs_fullpath.append('../output/testing/XMAS_172/N14_n172.csx')
if ID =='H':
    chixs_fullpath.append('../output/testing/XMAS_172/H1_n172.csx')
N_density = []
N_density.append(1.)
data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
print(data.keys())


source_def = {}
source_def["particle_type"] = 'neutron'
source_def["energy"] = 18.

outp = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)

Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp, './')

if ID =='Ar':
    A = np.loadtxt('../output/testing/XMAS_172/mcnp_Ar36_flx_tally.txt')
if ID =='Al':
    A = np.loadtxt('../output/testing/XMAS_172/mcnp_Al27_flx_tally.txt')
if ID =='N':
    A = np.loadtxt('../output/testing/XMAS_172/mcnp_N14_flx_tally.txt')
if ID =='H':
    A = np.loadtxt('../output/testing/XMAS_172/mcnp_H1_flx_tally.txt')
mcnpE = A[:,0]
mcnpF = A[:,1] * 4e9

# dE = np.diff(np.insert(mcnpE,0,0.))
dE = np.diff(mcnpE)
spectrum = mcnpF[1:]/dE

E = []
F = []
for g in range(A.shape[0]-1):
    E += [mcnpE[g], mcnpE[g+1]]
    F += [spectrum[g], spectrum[g]]

E = np.asarray(E)
F = np.asarray(F)

fig_n = plt.gcf().number
print(fig_n)
fig = plt.figure(fig_n)
ax_list = fig.axes
print(ax_list)
ax_list[0].semilogy(E, F, '+--', label='mcnp')
ax_list[1].loglog(E, F,'+--', label='mcnp')
plt.legend()
plt.show()