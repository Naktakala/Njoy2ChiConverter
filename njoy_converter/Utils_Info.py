import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#====================================================================
def InfiniteMediumSpectrum(data):
  neutron_gs = data["neutron_gs"]
  gamma_gs = data["gamma_gs"]
  sig_t = data["sigma_t"]
  transfer_mats = data["transfer_matrices"]
  transfer_mats_nonzeros = data["transfer_matrices_sparsity"]

  #======================================= Convert data to numpy data
  # for infinite medium
  G         = np.size(sig_t)
  G_neutron = len(neutron_gs)
  G_gamma   = len(gamma_gs)

  v_sig_t = np.array(sig_t)
  M_sig_gp_to_g = np.zeros([G,G])

  for gprime in range(0,G):
    for g in transfer_mats_nonzeros[0][gprime]:
      M_sig_gp_to_g[gprime,g] = transfer_mats[0][gprime,g]

  S = M_sig_gp_to_g.transpose()
  T = np.diag(v_sig_t)

  A = T - S
  A_inv = np.linalg.inv(A)

  v_src = np.zeros(G)
  v_src[0] = 1.0

  v_psi = np.matmul(A_inv,v_src)

  #======================================= Build data/energy
  neutron_group_bndries = []
  neutron_bndry_edge_values = []
  gamma_group_bndries = []
  gamma_bndry_edge_values = []
  for g in range(0,len(neutron_gs)):
    gprime = G_neutron - g -1

    lo_bound = neutron_gs[g][1]/1e6
    hi_bound = neutron_gs[g][2]/1e6

    bin_width = hi_bound-lo_bound

    neutron_group_bndries.append(lo_bound)
    neutron_group_bndries.append(hi_bound)

    neutron_bndry_edge_values.append(v_psi[gprime]/bin_width)
    neutron_bndry_edge_values.append(v_psi[gprime]/bin_width)

  for g in range(0,len(gamma_gs)):
    gprime = (G_gamma - g -1) + G_neutron

    lo_bound = gamma_gs[g][1]/1e6
    hi_bound = gamma_gs[g][2]/1e6

    bin_width = hi_bound-lo_bound

    gamma_group_bndries.append(lo_bound)
    gamma_group_bndries.append(hi_bound)

    gamma_bndry_edge_values.append(v_psi[gprime]/bin_width)
    gamma_bndry_edge_values.append(v_psi[gprime]/bin_width)

  last_nval = neutron_bndry_edge_values[len(neutron_bndry_edge_values)-1]
  if (gamma_bndry_edge_values != []):
    last_gval = gamma_bndry_edge_values[len(gamma_bndry_edge_values)-1]

  for i in range(0,len(neutron_bndry_edge_values)):
    neutron_bndry_edge_values[i] /= last_nval

  if (gamma_bndry_edge_values != []):
    for i in range(0,len(gamma_bndry_edge_values)):
      gamma_bndry_edge_values[i] /= last_gval

  plt.figure(figsize=(6,6))
  plt.plot(neutron_group_bndries, neutron_bndry_edge_values)
  plt.yscale("log")
  plt.xlabel("Energy (MeV)")
  plt.ylabel("$\phi(E)$")

  if (gamma_bndry_edge_values != []):
    plt.figure(figsize=(6,6))
    plt.plot(gamma_group_bndries, gamma_bndry_edge_values)
    plt.yscale("log")


#====================================================================
def ComputeKinf(data):
  #======================================= Get relevant data
  neutron_gs = data["neutron_gs"]
  sig_t = np.array(data["sigma_t"])
  sig_a = np.array(data["sigma_a"])
  chi_p = np.array(data["chi_prompt"])
  nu_p = np.array(data["nu_prompt"])
  sig_f = np.array(data["sigma_f"])
  transfer_mats = np.array(data["transfer_matrices"])
  transfer_mats_nonzeros = data["transfer_matrices_sparsity"]
  G = np.size(sig_t)

  #======================================= Format data for k-problem
  M_sig_gp_to_g = np.zeros([G,G])
  for gprime in range(0,G):
    for g in transfer_mats_nonzeros[0][gprime]:
      M_sig_gp_to_g[gprime,g] = transfer_mats[0][gprime,g]

  T = np.diag(sig_t)
  S = M_sig_gp_to_g.transpose()
  psi = np.linalg.solve(T-S, chi_p)

  #======================================= Compute k
  k = np.sum(nu_p*sig_f*psi) / np.sum(sig_a*psi)
  print("k_inf:\t{:.6f}".format(k))

  #======================================= Build data/energy
  neutron_group_bndries = []
  neutron_bndry_edge_values = []
  for g in range(0,len(neutron_gs)):
    gprime = G - g -1

    lo_bound = neutron_gs[g][1]/1e6
    hi_bound = neutron_gs[g][2]/1e6

    bin_width = hi_bound-lo_bound

    neutron_group_bndries.append(lo_bound)
    neutron_group_bndries.append(hi_bound)

    neutron_bndry_edge_values.append(psi[gprime]/bin_width)
    neutron_bndry_edge_values.append(psi[gprime]/bin_width)

  plt.figure(figsize=(6,6))
  plt.plot(neutron_group_bndries, neutron_bndry_edge_values)
  plt.yscale("log")
  plt.xlabel("Energy (MeV)")
  plt.ylabel("$\phi(E)$")
