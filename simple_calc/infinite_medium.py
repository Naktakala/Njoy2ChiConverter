# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 00:55:42 2022

@author: jean.ragusa
"""
import sys, os
sys.path.insert(0, os.path.realpath('../chitech_composite_processor'))
import Utils_ChiTechCombiner
import Utils_Info
import Utils_NjoySpectrumPlotter
import matplotlib.pyplot as plt
plt.close('all')

chixs_fullpath = []
chixs_fullpath.append('../output/testing/LANL187_LANL48/N14_n187g48.csx')
N_density = []
N_density.append(1.)
data = Utils_ChiTechCombiner.BuildCombinedChiTechData(chixs_fullpath, N_density)
print(data.keys())


source_def = {}
source_def["particle_type"] = 'neutron'
source_def["energy"] = 'fission'

outp = Utils_Info.InfiniteMediumSpectrum(data, source_def, plot=True)

Utils_NjoySpectrumPlotter.Njoy_spectrum_plotter(outp, './')