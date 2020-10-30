import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#====================================================================
def BuildCombinedData(raw_njoy_data, plot=False):
  ''' Combines enjoy raw data into a dictionary of vectors and matrices '''
  # ================================= Get dictionaries
  group_structures = raw_njoy_data["group_structures"]
  cross_sections = raw_njoy_data["cross_sections"]
  transfer_matrices = raw_njoy_data["transfer_matrices"]

  # ================================= Determine # of groups
  neutn_gs = group_structures["neutron"]
  gamma_gs = group_structures["gamma"] if "gamma" in group_structures else []

  G_n = len(neutn_gs)
  G_g = len(gamma_gs)
  G   = G_n + G_g 
  
  # ================================= Combine sig_t
  sig_t = np.zeros(G)
  
  if ("(n,total)" in cross_sections):
    sig_t_ndata = cross_sections["(n,total)"]
    for entry in sig_t_ndata:
      g = entry[0]
      v = entry[1]
      sig_t[G_n-g-1] += v

  if ("(g,total)" in cross_sections):
    sig_t_gdata = cross_sections["(g,total)"]
    for entry in sig_t_gdata:
      g = entry[0]
      v = entry[1]
      sig_t[G_n + G_g-g-1] += v

  # ================================= Inverse velocity term
  inv_v = np.zeros(G)

  if ("inv_velocity" in cross_sections):
    inv_v_ndata = cross_sections["inv_velocity"]
    for entry in inv_v_ndata:
      g = entry[0]
      v = entry[1]
      inv_v[G_n-g-1] += v

  # ================================= Fission data
  sig_f = np.zeros(G)
  if ("(n,fission)" in cross_sections):
    sig_f_data = cross_sections["(n,fission)"]
    for entry in sig_f_data:
      g = entry[0]
      v = entry[1]
      sig_f[G_n-g-1] += v

  nu_total = np.zeros(G)
  if ("total_nubar" in cross_sections):
    total_nubar_data = cross_sections["total_nubar"]
    for entry in total_nubar_data:
      g = entry[0]
      v = entry[1]
      nu_total[G_n-g-1] += v
  
  nu_prompt = np.zeros(G)
  if ("prompt_nubar" in cross_sections):
    prompt_nubar_data = cross_sections["prompt_nubar"]
    for entry in prompt_nubar_data:
      g = entry[0]
      v = entry[1]
      nu_prompt[G_n-g-1] += v

  chi_prompt = np.zeros(G)
  if ("prompt_chi" in cross_sections):
    prompt_chi_data = cross_sections["prompt_chi"]
    for entry in prompt_chi_data:
      g = entry[0]
      v = entry[1]
      chi_prompt[G_n-g-1] += v

  # ================================= Delayed neutron data
  decay_const=[]
  if ("decay_constants" in cross_sections):
    for entry in cross_sections["decay_constants"]:
      decay_const.append(entry[1])
    decay_const = np.asarray(decay_const)
  J = len(decay_const) # number of precursor groups


  nu_delayed = np.zeros(G)
  if ("delayed_nubar" in cross_sections):
    delayed_nubar_data = cross_sections["delayed_nubar"]
    for entry in delayed_nubar_data:
      g = entry[0]
      v = entry[1]
      nu_delayed[G_n-g-1] += v

  chi_delayed = np.zeros((G,J))
  if ("delayed_chi" in cross_sections):
    delayed_chi_data = cross_sections["delayed_chi"]
    for entry in delayed_chi_data:
      g = entry[0]
      v = entry[1:]
      chi_delayed[G_n-g-1] += v

  gamma = np.zeros(J)
  if (np.sum(nu_delayed)>0 and np.sum(chi_delayed)>0):
    gamma = np.sum(chi_delayed,axis=0)
  chi_delayed /= gamma

  # ================================= Combine transfer matrices
  # Keys available to all isotopes
  n_to_n_elastic_keys = ["(n,elastic)"]
  n_to_n_freegas_keys = ["mt221"]
  n_to_g_transfer_keys = [ "(n,g)", "(n,inel)", "(n,np)", 
                           "(n,nd)", "(n,p)", "(n,d)", 
                           "(n,t)", "(n,a)"]
  g_to_g_transfer_keys = ["(g,coherent)", "(g,incoherent)", \
                          "(g,pair_production)"]
  # graphite, H in ZrH, Zr in ZrH
  n_to_n_sab_elastic_keys = ["mt230", "mt226", "mt236"] 
  # H in H2O, graphite, H in ZrH, Zr in ZrH
  n_to_n_sab_inelastic_keys = ["mt222", "mt229", "mt225", "mt235"] 
  sab_treatment = False # has S(alpha, beta)

  # ===== Get the transfer matrices
  # Adding all the elastic scattering data
  nranges_to_nranges_elastic = []
  for rxn in n_to_n_elastic_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_elastic.append(transfer_matrices[rxn])

  # Adding all the (n,nxx) data, inelastic data
  nranges_to_nranges_inelastic = []
  for nn in range(1,40+1):
    rx_name = "(n,n{:02d})".format(nn)
    if rx_name in transfer_matrices:
      mat = transfer_matrices[rx_name]
      nranges_to_nranges_inelastic.append(mat)
  for nxn in range(2,4+1):
    rx_name = "(n,{:01d}n)".format(nxn)
    if rx_name in transfer_matrices:
      mat = transfer_matrices[rx_name]
      nranges_to_nranges_inelastic.append(mat)
  if "(n,nc)" in transfer_matrices:
    mat = transfer_matrices["(n,nc)"]
    nranges_to_nranges_inelastic.append(mat)

  # Adding all the free gas elastic scattering data
  nranges_to_nranges_freegas = []
  for rxn in n_to_n_freegas_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_freegas.append(transfer_matrices[rxn])

  # Adding all the elastic scattering S(\alpha,\beta) data
  nranges_to_nranges_sab_el = []
  for rxn in n_to_n_sab_elastic_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_sab_el.append(transfer_matrices[rxn])

  # Adding all the inelastic scattering S(\alpha,\beta) data
  nranges_to_nranges_sab_inel = []
  for rxn in n_to_n_sab_inelastic_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_sab_inel.append(transfer_matrices[rxn])
      sab_treatment = True

  # Adding all the neutron to gamma data
  nranges_to_granges = []
  for rxn in n_to_g_transfer_keys:
    if rxn in transfer_matrices:
      nranges_to_granges.append(transfer_matrices[rxn])

  # Adding all the gamma to gamma data
  granges_to_granges = []
  for rxn in g_to_g_transfer_keys:
    if rxn in transfer_matrices:
      granges_to_granges.append(transfer_matrices[rxn])

  # ===== Computing the max number of moments
  max_num_moms = 0
  for range_data in nranges_to_nranges_elastic:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_inelastic:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_freegas:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_sab_el:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_sab_inel:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_granges:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in granges_to_granges:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)

  # ===== Initializing the transfer matrices
  transfer_mats = []
  for m in range(0,max_num_moms):
    transfer_mats.append(np.zeros((G,G)))
  transfer_mats_nonzeros = []

  #=======================================
  # Lambda-ish to add neutron data
  def AddTransferNeutron(data_vals,offset=0):
    for entry in data_vals:
      num_moms = len(entry)-2
      gprime = G_n - entry[0] - 1
      g      = G_n - entry[1] - 1 + offset
    
      for m in range(0,num_moms):
        v = entry[m+2]
        transfer_mats[m][gprime,g] += v

  #=======================================
  # Lambda-ish to add gamma data
  def AddTransferGamma(data_vals):
    # (g,coherent)
    for entry in data_vals:
      gprime = G_n + G_g - entry[0] - 1
      g      = G_n + G_g - entry[1] - 1

      # Determine number of moments
      num_moms = len(entry)-2
      for m in range(0,num_moms):
        v = entry[m+2]
        transfer_mats[m][gprime,g] += v

  # ===== Format the transfer matrices, store scattering
  transfer_el = np.copy(transfer_mats)
  transfer_inel = np.copy(transfer_mats)
  transfer_sab_el = np.copy(transfer_mats)
  transfer_sab_inel = np.copy(transfer_mats)
  transfer_freegas = np.copy(transfer_mats)
  transfer_std = np.copy(transfer_mats)
  transfer_sab = np.copy(transfer_mats)

  # Regular elastic scatter (MT-2)
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_elastic:
    AddTransferNeutron(range_data)
  transfer_el = np.copy(transfer_mats)
  sig_s_el = np.sum(transfer_el[0],axis=1)

  # Regular inelastic scattering
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_inelastic:
    AddTransferNeutron(range_data)
  transfer_inel = np.copy(transfer_mats)
  sig_s_inel = np.sum(transfer_inel[0],axis=1)

  # Freegas thermal scattering
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_freegas:
    AddTransferNeutron(range_data)
  transfer_freegas = np.copy(transfer_mats)
  sig_s_freegas = np.sum(transfer_freegas[0],axis=1)

  # Elastic S(alpha, beta)
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_sab_el:
    AddTransferNeutron(range_data)
  transfer_sab_el = np.copy(transfer_mats)
  sig_s_sab_el = np.sum(transfer_sab_el[0],axis=1)
  
  # Inelastic S(alpha, beta)
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_sab_inel:
    AddTransferNeutron(range_data)
  transfer_sab_inel = np.copy(transfer_mats)
  sig_s_sab_inel = np.sum(transfer_sab_inel[0],axis=1)

  # Form elastic+inelastic matrices
  for m in range(0,max_num_moms):
    transfer_std[m] = transfer_el[m] + transfer_inel[m]
    transfer_sab[m] = transfer_sab_el[m] + transfer_sab_inel[m]
  
  # Replace MT-2 with thermal treatment
  transfer_mats = np.copy(transfer_std)
  for m in range(0,max_num_moms):
    if (sab_treatment): mat_m = transfer_sab[m]
    else: mat_m = transfer_freegas[m]
    for gprime in range(0,G):
      for g in range(0,G):
        if (np.abs(mat_m[gprime,g]) > 1.0e-18):
          transfer_mats[m][gprime,g] = mat_m[gprime,g]
        
  # Regular n,\gamma transfer <- add to above
  for range_data in nranges_to_granges:
    AddTransferNeutron(range_data,offset=G_g)

  # Regular \gamma,\gamma transfer <- add to above
  for range_data in granges_to_granges:
    AddTransferGamma(range_data)

  # ===== Compute total scattering and absorption
  sig_s_sab = sig_s_sab_el + sig_s_sab_inel
  sig_s = np.sum(transfer_mats[0],axis=1)
  sig_a = sig_t - sig_s
      
  # ===== Determine sparsity of the transfer matrices
  for m in range(0,max_num_moms):
    mat_non_zeros = []
    for gprime in range(0,G):
      non_zeros = []
      for g in range(0,G):
        if (abs(transfer_mats[m][gprime,g]) > 1.0e-18):
          non_zeros.append(g)
          mat_non_zeros.append(non_zeros)
    transfer_mats_nonzeros.append(mat_non_zeros)

  # ===== Plot the matrix
  if plot:
    Atest = transfer_mats[0]
    nz = np.nonzero(Atest)
    Atest[nz] = np.log10(Atest[nz]) + 10.0
    
    plt.figure(figsize=(6,6))
    im = plt.imshow(Atest, cmap=cm.Greys)
    plt.xticks(np.arange(0,G,10), [str(g) for g in range(0,G,10)])
    plt.yticks(np.arange(0,G,10), [str(g) for g in range(0,G,10)])
    plt.xlabel('Destination energy group')
    plt.ylabel('Source energy group')
    plt.gca().xaxis.set_ticks_position('top')
    plt.gca().xaxis.set_label_position('top')
    # plt.savefig("SERPENTTransferMatrix.png")

    #================================== Build group structures
    np_neutn_gs = np.matrix(neutn_gs)
    nbin_lo = np_neutn_gs[:,1]
    nbin_hi = np_neutn_gs[:,2]
    nbin_center = (0.5*(nbin_lo+nbin_hi))[::-1]
  
    #================================== Plot Cross sections
    plt.figure(figsize=(6,6))
    plt.semilogx(nbin_center,sig_t,label=r"$\sigma_t$")
    plt.semilogx(nbin_center,sig_a,label=r"$\sigma_a$")
    plt.semilogx(nbin_center,sig_s,label=r"$\sigma_s$")
    plt.legend()
    plt.show()
  
    #================================== Plot scattering
    plt.semilogx(nbin_center,sig_s,label=r"$\sigma_s$")
    plt.semilogx(nbin_center,sig_s_el,label=r"$\sigma_s$ elastic")
    plt.semilogx(nbin_center,sig_s_inel,label=r"$\sigma_s$ inelastic")
    if np.sum(sig_s_freegas) > 0.0:
      plt.semilogx(nbin_center,sig_s_freegas,label=r"$\sigma_s$ freegas")
    if np.sum(sig_s_sab) > 0.0:
      plt.semilogx(nbin_center,sig_s_sab,label=r"$\sigma_s$ total SAB")
      plt.semilogx(nbin_center,sig_s_sab_el,label=r"$\sigms_s$ elastic SAB")
      plt.semilogx(nbin_center,sig_s_sab_inel,label=r"$\sigms_s$ inelastic SAB")
    plt.legend()
    plt.show()

  #================================== Build return data
  return_data = {}
  return_data["neutron_gs"] = neutn_gs
  return_data["gamma_gs"] = gamma_gs
  return_data["sigma_t"] = sig_t
  return_data["sigma_a"] = sig_a
  return_data["sigma_s"] = sig_s
  return_data["sigma_s_freegas"] = sig_s_freegas
  return_data["sigma_s_sab"] = sig_s_sab
  return_data["sigma_s_sab_el"] = sig_s_sab_el
  return_data["sigma_s_sab_inel"] = sig_s_sab_inel
  return_data["sigma_s_el"] = sig_s_el
  return_data["sigma_s_inel"] = sig_s_inel
  return_data["sigma_f"] = sig_f
  return_data["nu_total"] = nu_total
  return_data["nu_prompt"] = nu_prompt
  return_data["nu_delayed"] = nu_delayed
  return_data["chi_prompt"] = chi_prompt
  return_data["chi_delayed"] = chi_delayed
  return_data["decay_constants"] = decay_const
  return_data["gamma"] = gamma
  return_data["inv_velocity"] = inv_v
  return_data["transfer_matrices"] = transfer_mats
  return_data["transfer_matrices_sparsity"] = transfer_mats_nonzeros
  return return_data 