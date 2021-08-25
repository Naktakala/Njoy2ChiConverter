import numpy as np 
import matplotlib.pyplot as plt
#====================================================================
def GenerateSpectrumData(neutron_gs, psi, sigma_heat, gamma_gs=[]):
  G_neutron = len(neutron_gs)
  G_gamma   = len(gamma_gs)
  G         = G_neutron + G_gamma
  assert len(psi) == G, \
    "Total neutron+gamma groups not compatible with psi."
  assert len(sigma_heat) == G, \
    "Total neutron+gamma groups not compatible with heating XS"

  #==================== Generate neutron spectrum data
  n_bndrys, n_vals = [], []
  n_heating = []
  for g in range(G_neutron):
    gprime    = G_neutron - g - 1
    lo_bound  = neutron_gs[g][1]*1.0e-6
    hi_bound  = neutron_gs[g][2]*1.0e-6
    bin_width = hi_bound-lo_bound
    spectrum  = psi[G_neutron-g-1] / bin_width
    n_bndrys += [lo_bound, hi_bound]
    n_vals   +=  [spectrum, spectrum]

    heating_spectrum  = sigma_heat[G_neutron-g-1]
    n_heating += [heating_spectrum, heating_spectrum]

  n_bndrys = np.array(n_bndrys)
  n_vals   = np.array(n_vals)
  n_heating = np.array(n_heating)

  #==================== Generate gamma spectrum data
  g_bndrys, g_vals = [], []
  g_heating = []
  for g in range(G_gamma):
    gprime = (G_gamma - g - 1) + G_neutron
    lo_bound  = gamma_gs[g][1]*1.0e-6
    hi_bound  = gamma_gs[g][2]*1.0e-6
    bin_width = hi_bound-lo_bound
    spectrum  = psi[gprime] / bin_width
    g_bndrys += [lo_bound, hi_bound]
    g_vals   += [spectrum, spectrum]

    heating_spectrum  = sigma_heat[gprime]
    g_heating += [heating_spectrum, heating_spectrum]

  g_bndrys = np.array(g_bndrys)
  g_vals   = np.array(g_vals)
  g_heating = np.array(g_heating)

  return [(n_bndrys, n_vals), (g_bndrys, g_vals), (n_heating, g_heating)]

#====================================================================
def InfiniteMediumSpectrum(data, path, plot=False):
  neutron_gs = data["neutron_gs"]
  gamma_gs = data["gamma_gs"]
  sig_t = data["sigma_t"]
  transfer_mats = data["transfer_matrices"]
  transfer_mats_nonzeros = data["transfer_matrices_sparsity"]
  sig_heat = data["sigma_heat"]

  #======================================= Get data from dictionary
  G_neutron = len(neutron_gs)
  G_gamma   = len(gamma_gs)
  coupled_txt = ''
  if G_gamma>0:
    coupled_txt = '_coupled_ng'

  G = np.size(sig_t)
  
  v_sig_t = np.array(sig_t)
  M_sig_gp_to_g = np.zeros([G,G])

  for gprime in range(G):
    for g in transfer_mats_nonzeros[0][gprime]:
      M_sig_gp_to_g[gprime,g] = transfer_mats[0][gprime,g]

  #======================================= Solve for psi
  S = M_sig_gp_to_g.transpose()
  T = np.diag(v_sig_t)

  if plot:
    fig = plt.figure()
    plt.matshow(np.log(S))
    if G_gamma>0:
      g_txt = '_g'+str(G_gamma)
    else:
      g_txt = ''
    filename = path+'transfert_matrix_n'+str(G_neutron)+g_txt+'.png'
    plt.savefig(filename)

  A = T - S
  A_inv = np.linalg.inv(A)

  v_src = np.zeros(G)
  v_src[0] = 1.0

  v_psi = np.matmul(A_inv,v_src)
  print("Norm spectrum: ")
  print(np.linalg.norm(v_psi))

  #======================================= Build data/energy
  outp = GenerateSpectrumData(neutron_gs, v_psi, sig_heat, gamma_gs)
  neutron_group_bndries = outp[0][0]
  neutron_spectrum = outp[0][1]
  gamma_group_bndries = outp[1][0]
  gamma_spectrum = outp[1][1]
  neutron_heating_spectrum = outp[2][0]
  gamma_heating_spectrum = outp[2][1]

  #========================== Compute the heating rate for njoy
  for i in range (0, len(neutron_heating_spectrum)):
    neutron_heating_spectrum[i] *= neutron_spectrum[i]
  for i in range (0, len(gamma_heating_spectrum)):
    gamma_heating_spectrum[i] *= gamma_spectrum[i]

  
  #========================== Check for type of problems ==================#
  #First NJOY normalization
  #======================================= Plot the spectra
  if plot:
    #================================= Plot energy spectrum
    #================================= Flux
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax = fig.add_subplot(1, 2, 1)
    ax.semilogy(neutron_group_bndries, neutron_spectrum)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("$\phi(E)$")
    ax.grid('on')
    ax = fig.add_subplot(1, 2, 2)
    ax.loglog(neutron_group_bndries, neutron_spectrum)
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("$\phi(E)$")
    ax.grid('on')
    plt.suptitle('neutron spectrum')
    plt.savefig(path+'Neutron_spectrum_'+str(G_neutron)+coupled_txt+'.png')

    #================================= Heating XS
    fig = plt.figure(figsize=(12,6))
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    ax = fig.add_subplot(1, 2, 1)
    ax.semilogy(neutron_group_bndries, neutron_heating_spectrum,color='r')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("H(E) (eV/s)")
    ax.grid('on')
    ax = fig.add_subplot(1, 2, 2)
    ax.loglog(neutron_group_bndries, neutron_heating_spectrum, color='r')
    ax.set_xlabel("Energy (MeV)")
    ax.set_ylabel("H(E) (eV/s)")
    ax.grid('on')
    plt.suptitle('neutron heating')
    plt.savefig(path+'Neutron_Heating_spectrum_'+str(G_neutron)+coupled_txt+'.png')

    #Testing
    if (gamma_spectrum != []):
    #================================= Flux
      fig = plt.figure(figsize=(12,6))
      fig.subplots_adjust(hspace=0.4, wspace=0.4)
      ax = fig.add_subplot(1, 2, 1)
      ax.semilogy(gamma_group_bndries, gamma_spectrum)
      ax.set_xlabel("Energy (MeV)")
      ax.set_ylabel("$\phi(E)$")
      ax.grid('on')
      ax = fig.add_subplot(1, 2, 2)
      ax.loglog(gamma_group_bndries, gamma_spectrum)
      ax.set_xlabel("Energy (MeV)")
      ax.set_ylabel("$\phi(E)$")
      ax.grid('on')
      plt.suptitle('gamma spectrum')
      plt.savefig(path+'Gamma_spectrum_'+str(G_gamma)+'.png')

    #================================= Heating
      fig = plt.figure(figsize=(12,6))
      fig.subplots_adjust(hspace=0.4, wspace=0.4)
      ax = fig.add_subplot(1, 2, 1)
      ax.semilogy(gamma_group_bndries, gamma_heating_spectrum, color ='r')
      ax.set_xlabel("Energy (MeV)")
      ax.set_ylabel("H(E) (eV/s)")
      ax.grid('on')
      ax = fig.add_subplot(1, 2, 2)
      ax.loglog(gamma_group_bndries, gamma_heating_spectrum, color='r')
      ax.set_xlabel("Energy (MeV)")
      ax.set_ylabel("H(E) (eV/s)")
      ax.grid('on')
      plt.suptitle('gamma heating')
      plt.savefig(path+'Gamma_heating_spectrum_'+str(G_gamma)+'.png')


#====================================================================
def ComputeKinf(data):
  #======================================= Get relevant data
  neutron_gs = data["neutron_gs"]
  sig_t = np.array(data["sigma_t"])
  ## sig_a = np.array(data["sigma_a"])
  chi_p = np.array(data["chi_prompt"])
  try:
    nu_p = np.array(data["nu_prompt"])
    sig_f = np.array(data["sigma_f"])
    nu_sig_f = nu_p * sig_f
  except:
    nu_sig_f = np.array(data['nu_sigma_f'])
  transfer_mats = np.array(data["transfer_matrices"])
  G = np.size(sig_t)
  
  #======================================= Solve the eigenproblem
  M_sig_gp_to_g = np.zeros([G,G])
  for gprime in range(G):
    for g in range(G):
      M_sig_gp_to_g[gprime,g] = transfer_mats[0][gprime,g]

  T = np.diag(sig_t)
  S = M_sig_gp_to_g.transpose()
  psi = np.linalg.solve(T-S, chi_p)

  #======================================= Compute k
  k = np.sum(nu_sig_f*psi)
  print("k_inf:\t{:.6g}".format(k))
  print("rho:\t{:.8g}".format(1.0e5*(k-1)/k))

  #======================================= Build data/energy
  outp = GenerateSpectrumData(neutron_gs, psi)
  neutron_group_bndries = outp[0][0]
  neutron_spectrum = outp[0][1]

  fig = plt.figure(figsize=(6,6))
  plt.loglog(neutron_group_bndries, neutron_spectrum)
  # plt.xlabel("Energy (MeV)")
  plt.ylabel("$\phi(E)$")
  plt.grid(True)
  plt.show()
