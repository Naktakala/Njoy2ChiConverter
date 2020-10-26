''' Converts neutron gamma njoy output to chitech-format '''
import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import Utils_ReadNJOYOutput as Reader
import Utils_Combiner as Combiner
import Utils_ChiWriter as Writer
import Utils_Info as Info

#####################################################################
# Stand-alone usage
if __name__ == "__main__":
  from shutil import copyfile
  
  njoy_dir = "njoy_xs"
  chi_dir = "chi_xs"

  assert os.path.isdir(njoy_dir)
  assert os.path.isdir(chi_dir)

  # mt1 = np.loadtxt('graphite_xs/graphite_mt1.txt')
  # mt2 = np.loadtxt('graphite_xs/graphite_mt2.txt')
  # mt230 = np.loadtxt('graphite_xs/graphite_mt230.txt')
  # mt229 = np.loadtxt('graphite_xs/graphite_mt229.txt')

  gs_dirs = []
  temp_dirs = []
  isotopes = []

  plt.close('all')
  for gs in os.listdir(njoy_dir):
    if gs in gs_dirs or len(gs_dirs) == 0:
      njoy_gs = "/".join([njoy_dir,gs])
      chi_gs = "/".join([chi_dir,gs])

      for temp in os.listdir(njoy_gs):
        if temp in temp_dirs or len(temp_dirs) == 0:
          njoy_temp = "/".join([njoy_gs,temp])
          chi_temp = "/".join([chi_gs,temp])
          if not os.path.isdir(chi_temp):
            os.makedirs(chi_temp)
          
          for isotope in os.listdir(njoy_temp):
            if isotope.split(".")[0] in isotopes or len(isotopes) == 0:
              njoy_file = "/".join([njoy_temp,isotope])
              chi_file = "/".join([chi_temp,isotope]).replace(".txt",".csx")
              
              print("\nPARSING FILE: " + njoy_file)
              raw_njoy_data = Reader.ReadNJOYfile(njoy_file)
              data = Combiner.BuildCombinedData(raw_njoy_data,plot=False)
              Writer.WriteChiTechFile(data, chi_file)
              # Info.InfiniteMediumSpectrum(data)
              # Info.ComputeKinf(data)

  # plt.show()
