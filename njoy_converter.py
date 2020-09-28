''' Converts neutron gamma njoy output to chitech-format '''
import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#====================================================================
def PrintLine(line_num,line):
  ''' Prints a file line without the line feed. '''
  print("Line: ",line_num,line,end='') 

#====================================================================
def ProcessGroupStructure(nL_i,lines):
  ''' Processes a group structure.
    params:
      nL_i  File line number before data starts.
      lines An array of files lines.
      
    return:
      returns the group_structure '''
  group_struct = []
  nL = nL_i + 1
  words = lines[nL].split()
  while (True):
    if (len(words) == 0):
      nL += 1; words = lines[nL].split()
    if (words[0] == "1"): break
    nL += 1; words = lines[nL].split()

  while (True):
    words = lines[nL].split()
    if not words: 
      break 

    group_struct.append( [int(words[0])-1,
                          float(words[1]),
                          float(words[3])])
    nL += 1
  return group_struct

#====================================================================
def ProcessCrossSection(nL_i, lines,header_size=5,line_incr=1):
  ''' Reads 1D neutron cross-sections from the lines. 
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.
    
    return:
      A table containing the group wise xs. '''
  xs=[]
  # Skip header
  nL = nL_i + header_size
  words = lines[nL].split()
  num_words = len(words)

  while (num_words>=2):
    number_word = words[1]
    number_word = number_word.replace("+","E+")
    number_word = number_word.replace("-","E-")
    xs.append([int(words[0])-1,float(number_word)])
    nL += line_incr

    words = lines[nL].split()
    num_words = len(words)

  return xs

#====================================================================
def ProcessPromptChi(nL_i, lines, header_size=4, line_incr=1):
  ''' Reads 1D prompt fission spectrum from the lines. 
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.
    
    return:
      A table containing the group wise xs. '''
  xs=[]
  # Skip header
  nL = nL_i + header_size
  words = lines[nL].split()
  num_words = len(words)

  g = 0
  while (num_words>=2):
    if not words[0].isdigit():
      nL += 2
      words = lines[nL].split()
      num_words = len(words)

    for number_word in words[1:]:
      number_word = number_word.replace("+","E+")
      number_word = number_word.replace("-","E-")
      xs.append([g, float(number_word)])
      g += 1
      
    nL += line_incr
    words = lines[nL].split()
    num_words = len(words)
  return xs

#====================================================================
def ProcessDecayConstants(nL_i, lines):
  ''' Reads delayed neutron decay constants from delayed chi.
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.

    return:
      A table containing delayed neutron group decay constants.
  '''
  decay_const = []
  nL = nL_i + 2
  nL += 1 if lines[nL].split()[0]=="normalized" else 3
  words = lines[nL].split()
  for j, number_word in enumerate(words[1:]):
    number_word = number_word.replace("+","E+")
    number_word = number_word.replace("-","E-")
    decay_const.append([j, float(number_word)])
  return decay_const

#====================================================================
def ProcessDelayedChi(nL_i, lines):
  ''' Reads the delayed neutron spectra from the lines.
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.
    
    return:
      A table containing the group wise delayed 
      neutron precursor spectrum coefficients. '''
  matrix = []
  nL = nL_i + 2
  nL += 3 if lines[nL].split()[0]=="normalized" else 5
  words = lines[nL].split()
  num_words = len(words)
  while (num_words>=2):
    g = int(words[0])-1
    entry = [g]
    for j, number_word in enumerate(words[1:]):
      number_word = number_word.replace("+","E+")
      number_word = number_word.replace("-","E-")
      entry.append(float(number_word))
    matrix.append(entry)

    nL += 1
    words = lines[nL].split()
    num_words = len(words)

  return matrix

#====================================================================
def ProcessTransferMatrix(nL_i, lines):    
  ''' Reads transfer matrices from the lines. 
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.
    
    return:
      A table containing the group wise transfer coefficients. '''    
  matrix = []
  # Skip 4 lines
  add_skip = 0
  if lines[nL_i+1].find("particle emission")>=0 or \
       lines[nL_i+2].find("spectrum constant below")>=0:
    add_skip = 1
  nL = nL_i + 5 + add_skip
  words = lines[nL].split()
  num_words = len(words)

  while (num_words>=2):
    # Transform words 0.000-00 to 0.000E+00
    for i in range(0,num_words):
      word = words[i]
      loc_of_sign = word.find("-",2)
      if loc_of_sign == -1:
        loc_of_sign = word.find("+",2)
      first_part = word[0:loc_of_sign]
      secnd_part = word[loc_of_sign:]

      secnd_part = secnd_part.replace("+","E+")
      secnd_part = secnd_part.replace("-","E-")

      words[i] = first_part+secnd_part

    # Develop table entry
    if words[0] != "normalization":
      entry = []
      entry.append(int(words[0])-1) # gprime
      entry.append(int(words[1])-1) # g
      for i in range(2,num_words):
         entry.append(float(words[i]))
      matrix.append(entry)

    nL += 1
    words = lines[nL].split()
    num_words = len(words)
  return matrix  

#====================================================================
def ProcessTransferMatrixB(nL_i, lines, process_overflow=False):    
  ''' Reads transfer matrices from the lines. 
    params:
      nL_i    File line number before data starts.
      lines   An array of file lines.
      line_incr Line increment to use
    
    return:
      A table containing the group wise transfer coefficients. '''    
  matrix = []
  # Skip 3 lines
  nL = nL_i + 4
  words = lines[nL].split()
  num_words = len(words)

  while (num_words>=2):
    # Transform words 0.000-00 to 0.000E+00
    for i in range(0,num_words):
      word = words[i]
      loc_of_sign = word.find("-",2)
      if loc_of_sign == -1:
        loc_of_sign = word.find("+",2)
      first_part = word[0:loc_of_sign]
      secnd_part = word[loc_of_sign:]

      secnd_part = secnd_part.replace("+","E+")
      secnd_part = secnd_part.replace("-","E-")

      words[i] = first_part+secnd_part

    # Develop table entry
    entry = []
    entry.append(int(words[0])-1) # gprime
    entry.append(int(words[1])-1) # g
    for i in range(2,num_words):
      entry.append(float(words[i]))

    # Process overflow
    if process_overflow:
      words = lines[nL+1].split()

      word = words[0]
      loc_of_sign = word.find("-",2)
      if loc_of_sign == -1:
        loc_of_sign = word.find("+",2)
      first_part = word[0:loc_of_sign]
      secnd_part = word[loc_of_sign:]

      secnd_part = secnd_part.replace("+","E+")
      secnd_part = secnd_part.replace("-","E-")

      words[0] = first_part+secnd_part

      entry.append(float(words[0]))
      nL += 1

    # Complete matrix entry
    matrix.append(entry)

    # Move to next line
    nL += 1
    words = lines[nL].split()
    num_words = len(words)

    # Skip lines with xsec and heat
    while (num_words>=2 and not words[1].isnumeric()):
      nL += 1
      words = lines[nL].split()
      num_words = len(words)
      
  return matrix  

#====================================================================
def ReadNJOYfile(njoy_filename="output"):
  ''' Reads an NJOY output file
    params:
      njoy_filename Name of the NJOY file to process. 
      
    return:
      Returns a complex dictionary of raw data. '''

  njoy_file = open(njoy_filename, 'r')
  file_lines = njoy_file.readlines()

  njoy_raw_data = {}

  group_structures = {}
  cross_sections = {}
  transfer_matrices = {}

  flag_run_processed             = False
  flag_gamma_structure_processed = False

  nL=-1
  while (nL<(len(file_lines)-1)):
    nL += 1
    line = file_lines[nL]
    words = file_lines[nL].split()
    num_words = len(words)

    # if (line.find("neutron group structure") != -1 and words[3] != "option"):
    #   group_structures["neutron"] = ProcessGroupStructure(nL,file_lines)

    if (line.find("sigma zeroes") != -1):
      group_structures["neutron"] = ProcessGroupStructure(nL,file_lines)

    if (line.find("gamma group structure") != -1):
      if not flag_gamma_structure_processed:
        group_structures["gamma"] = ProcessGroupStructure(nL,file_lines)
        flag_gamma_structure_processed = True 

    if (line.find("for mf") != -1 and
        line.find("mt") != -1 ):

      if (words[2] == "3" and words[5] == "1" ):
        cross_sections["(n,total)"] = \
          ProcessCrossSection(nL,file_lines,line_incr=2)  

      if (words[2] == "3" and words[5] == "2" ):
        cross_sections["(n,elastic)"] = \
          ProcessCrossSection(nL,file_lines)        

      if (words[2] == "3" and words[5] == "4" ):
        cross_sections["(n,inelastic)"] = \
          ProcessCrossSection(nL,file_lines) 

      if (words[2] == "3" and words[5] == "16" ):
        cross_sections["(n,2n)"] = \
          ProcessCrossSection(nL,file_lines) 

      if (words[2] == "3" and words[4] == "mt259"):
        cross_sections["inv_velocity"] = \
          ProcessCrossSection(nL,file_lines)
      
      if (words[2] == "3" and words[5] == "18"):
        cross_sections["(n,fission)"] = \
          ProcessCrossSection(nL,file_lines)

      if (words[2] == "3" and words[4] == "mt452"):
        cross_sections["total_nubar"] = \
          ProcessCrossSection(nL,file_lines)

      if (words[2] == "3" and words[4] == "mt456"):
        cross_sections["prompt_nubar"] = \
          ProcessCrossSection(nL,file_lines)

      if (words[2] == "3" and words[4] == "mt455"): 
        cross_sections["delayed_nubar"] = \
          ProcessCrossSection(nL,file_lines)

      if (words[2] == "5" and words[5] == "18"):
        cross_sections["prompt_chi"] = \
          ProcessPromptChi(nL,file_lines,4)

      if (words[2] == "5" and words[4] == "mt455"):
        cross_sections["decay_constants"] = \
          ProcessDecayConstants(nL,file_lines)
        cross_sections["delayed_chi"] = \
          ProcessDelayedChi(nL,file_lines)

      if (words[num_words-1] == "matrix"):
        transfer_matrices[words[num_words-3]] = \
          ProcessTransferMatrix(nL,file_lines)

      if (words[1] == "mf23" and words[3] == "mt501" ):
        cross_sections["(g,total)"] = \
          ProcessCrossSection(nL,file_lines,header_size=4) 

      if (words[1] == "mf23" and words[3] == "mt502" ):
        cross_sections["(g,coherent)"] = \
          ProcessCrossSection(nL,file_lines,header_size=4)   

      if (words[1] == "mf23" and words[3] == "mt504" ):
        cross_sections["(g,incoherent)"] = \
          ProcessCrossSection(nL,file_lines,header_size=4)  

      if (words[1] == "mf23" and words[3] == "mt516" ):
        cross_sections["(g,pair_production)"] = \
          ProcessCrossSection(nL,file_lines,header_size=4)      

      if (words[1] == "mf26" and words[3] == "mt502"):
        transfer_matrices["(g,coherent)"] = \
          ProcessTransferMatrixB(nL,file_lines,process_overflow=True)

      if (words[1] == "mf26" and words[3] == "mt504"):
        transfer_matrices["(g,incoherent)"] = \
          ProcessTransferMatrixB(nL,file_lines,process_overflow=True)

      if (words[1] == "mf26" and words[3] == "mt516"):
        transfer_matrices["pair_production"] = \
          ProcessTransferMatrixB(nL,file_lines,process_overflow=False)         

  njoy_file.close()

  njoy_raw_data["group_structures"] = group_structures
  njoy_raw_data["cross_sections"] = cross_sections
  njoy_raw_data["transfer_matrices"] = transfer_matrices

  return njoy_raw_data

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
  precursor_fraction = np.zeros(J)
  if (np.sum(nu_delayed)>0 and np.sum(chi_delayed)>0):
    nu_bar_delayed = np.mean(nu_delayed)
    precursor_fraction = np.sum(chi_delayed,axis=0)
    gamma = nu_bar_delayed*precursor_fraction
  chi_delayed /= np.sum(chi_delayed,axis=0)

  # ================================= Combine transfer matrices
  # Keys available to all isotopes
  n_to_n_elastic_keys = ["(n,elastic)"]
  n_to_g_transfer_keys = [ "(n,g)", "(n,inel)", "(n,np)", 
                           "(n,nd)", "(n,p)", "(n,d)", 
                           "(n,t)", "(n,a)"]
  g_to_g_transfer_keys = ["(g,coherent)", "(g,incoherent)", \
                          "(g,pair_production)"]
  n_to_n_freegas_keys = ["mt221"]

  # Keys for S(alpha, beta) thermal corrections. Note that this 
  # will NOT be correct for inelastic scattering if data for 
  # H2O and ZrH exist in H1 files. As it stands, all reaction 
  # keys are additive. For isotopes with only one S(alpha, beta) 
  # material, all other reactions will be skipped over because 
  # those reactions do not live in the file. A flag could be used to 
  # choose the correct reaction numbers. 

  # graphite, H in ZrH, Zr in ZrH
  n_to_n_sab_elastic_keys = ["mt230", "mt226", "mt236"] 
  # H in H2O, graphite, H in ZrH, Zr in ZrH
  n_to_n_sab_inelastic_keys = ["mt222", "mt229", "mt226", "mt235"] 

  # ===== Get the transfer matrices
  # Adding all the elastic scattering data
  nranges_to_nranges_elastic = []
  for rxn in n_to_n_elastic_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_elastic.append(transfer_matrices[rxn])

  # Adding all the (n,nxx) data, inelastic data
  nranges_to_nranges_inelastic = []
  for nn in range(1,24+1):
    rx_name = "(n,n{:02d})".format(nn)
    if rx_name in transfer_matrices:
      mat = transfer_matrices[rx_name]
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
  transfer_mats_standard = np.copy(transfer_mats)
  transfer_mats_freegas  = np.copy(transfer_mats)
  transfer_mats_sab_el   = np.copy(transfer_mats)
  transfer_mats_sab_inel = np.copy(transfer_mats)

  # Regular elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_elastic:
    AddTransferNeutron(range_data)
  sig_s_el = np.sum(transfer_mats[0],axis=1)
  transfer_mats_standard = np.copy(transfer_mats)

  # Regular inelastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_inelastic:
    AddTransferNeutron(range_data)
  sig_s_inel = np.sum(transfer_mats[0],axis=1)
  for m in range(0,max_num_moms):
    transfer_mats_standard[m] += transfer_mats[m]

  # Regular n,\gamma transfer
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_granges:
    AddTransferNeutron(range_data,offset=G_g)
  for m in range(0,max_num_moms):
    transfer_mats_standard[m] += transfer_mats[m]

  # Regular \gamma,\gamma transfer
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in granges_to_granges:
    AddTransferGamma(range_data)
  for m in range(0,max_num_moms):
    transfer_mats_standard[m] += transfer_mats[m]
  
  # Freegas elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_freegas:
    AddTransferNeutron(range_data)
  sig_s_freegas = np.sum(transfer_mats[0],axis=1)
  transfer_mats_freegas = np.copy(transfer_mats)

  # S(alpha,beta) elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_sab_el:
    AddTransferNeutron(range_data)
  sig_s_el_sab = np.sum(transfer_mats[0],axis=1)
  transfer_mats_sab_el = np.copy(transfer_mats)
  
  # S(alpha,beta) inelastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_inelastic:
    AddTransferNeutron(range_data)
  sig_s_inel_sab = np.sum(transfer_mats[0],axis=1)
  transfer_mats_sab_inel = np.copy(transfer_mats)

  # ===== Store stuff
  # Uncorrected quantities
  sig_t_uncorr = sig_t
  sig_s_uncorr = np.sum(transfer_mats_standard[0],axis=1)
  sig_a = sig_t_uncorr - sig_s_uncorr
  # Correction terms
  sig_s_sab = sig_s_el_sab + sig_s_inel_sab
  
  # ===== Make the termal scattering corrections
  # This is a bit complex. First, the freegas transfer
  # matrices are overlaid on the standard transfer matrices.
  # Then, the same process is carried out with the 
  # inelastic S(alpha, beta) terms. Finally, the elastic
  # S(alpha, beta) terms are added to the result.
  transfer_mats = np.copy(transfer_mats_standard)
  for m in range(0,max_num_moms):
    m_transfer_mat_freegas = transfer_mats_freegas[m]
    m_transfer_mat_sab_inel = transfer_mats_sab_inel[m]
    for gprime in range(0,G):
      for g in range(0,G):
        if (np.abs(m_transfer_mat_freegas[gprime,g]) > 1.0e-18):
          transfer_mats[m][gprime,g] = m_transfer_mat_freegas[gprime,g]
        if (np.abs(m_transfer_mat_sab_inel[gprime,g]) > 1.0e-18):
          transfer_mats[m][gprime,g] = m_transfer_mat_sab_inel[gprime,g]
    transfer_mats[m] += transfer_mats_sab_el[m]

  # ===== Update cross section vectors
  sig_s = np.sum(transfer_mats[0],axis=1)
  sig_t = sig_a + sig_s
      
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

  #================================== Build return data
  return_data = {}
  return_data["neutron_gs"] = neutn_gs
  return_data["gamma_gs"] = gamma_gs
  return_data["sigma_t"] = sig_t
  return_data["sigma_t_uncorr"] = sig_t_uncorr
  return_data["sigma_a"] = sig_a
  return_data["sigma_s"] = sig_s
  return_data["sigma_s_uncorr"] = sig_s_uncorr
  return_data["sigma_s_freegas"] = sig_s_freegas
  return_data["sigma_s_sab"] = sig_s_sab
  return_data["sigma_s_el"] = sig_s_el
  return_data["sigma_s_inel"] = sig_s_inel
  return_data["sigma_s_sab_el"] = sig_s_el_sab
  return_data["sigma_s_sab_inel"] = sig_s_inel_sab
  return_data["sigma_f"] = sig_f
  return_data["nu_total"] = nu_total
  return_data["nu_prompt"] = nu_prompt
  return_data["nu_delayed"] = nu_delayed
  return_data["chi_prompt"] = chi_prompt
  return_data["chi_delayed"] = chi_delayed
  return_data["decay_constants"] = decay_const
  return_data["gamma"] = gamma
  return_data["precursor_fraction"] = precursor_fraction
  return_data["inv_velocity"] = inv_v
  return_data["transfer_matrices"] = transfer_mats
  return_data["transfer_matrices_sparsity"] = transfer_mats_nonzeros
  return return_data 

#====================================================================
def WriteChiTechFile(data,chi_filename="output.cxs",comment="# Output"):
  cf = open(chi_filename,'w')

  sig_t = data["sigma_t"]
  sig_s = data["sigma_s"]
  sig_f = data["sigma_f"]
  nu_total = data["nu_total"]
  nu_prompt = data["nu_prompt"]
  nu_delayed = data["nu_delayed"]
  chi_prompt = data["chi_prompt"]
  chi_delayed = data["chi_delayed"]
  decay_const = data["decay_constants"]
  precursor_fraction = data["precursor_fraction"]
  gamma = data["gamma"]
  ddt_coeff = data["inv_velocity"] 
  transfer_mats = data["transfer_matrices"]
  transfer_mats_nonzeros = data["transfer_matrices_sparsity"]

  G = np.size(sig_t)
  M = len(transfer_mats)
  J = len(decay_const)

  cf.write(comment+"\n")
  cf.write("NUM_GROUPS "+str(G)+"\n")
  cf.write("NUM_MOMENTS "+str(M)+"\n")
  cf.write("NUM_PRECURSORS "+str(J)+"\n")

  cf.write("SIGMA_T_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(sig_t[g]))
    cf.write("\n")
  cf.write("SIGMA_T_END"+"\n")

  cf.write("SIGMA_S_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(sig_s[g]))
    cf.write("\n")
  cf.write("SIGMA_S_END"+"\n")

  cf.write("SIGMA_F_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(sig_f[g]))
    cf.write("\n")
  cf.write("SIGMA_F_END"+"\n")

  cf.write("NU_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(nu_total[g]))
    cf.write("\n")
  cf.write("NU_END"+"\n")

  cf.write("NU_PROMPT_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(nu_prompt[g]))
    cf.write("\n")
  cf.write("NU_PROMPT_END"+"\n")

  cf.write("CHI_PROMPT_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(chi_prompt[g]))
    cf.write("\n")
  cf.write("CHI_PROMPT_END"+"\n")

  cf.write("DDT_COEFF_BEGIN"+"\n")
  for g in range(0,G):
    cf.write("{:<4d}".format(g)+ " ")
    cf.write("{:<g}".format(ddt_coeff[g]/100))
    cf.write("\n")
  cf.write("DDT_COEFF_END"+"\n")

  cf.write("TRANSFER_MOMENTS_BEGIN"+"\n")
  for m in range(0,M):
    cf.write("# l = " + str(m) + "\n")
    for gprime in range(0,G):
      for g in transfer_mats_nonzeros[m][gprime]:
        cf.write("M_GPRIME_G_VAL"+ " ")
        cf.write("{:<4d}".format(m)+ " ")
        cf.write("{:<4d}".format(gprime)+ " ")
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<g}".format(transfer_mats[m][gprime,g]))
        cf.write("\n")
  cf.write("TRANSFER_MOMENTS_END"+"\n")

  if J > 0:
    cf.write("NU_DELAYED_BEGIN"+"\n")
    for g in range(0,G):
      cf.write("{:<4d}".format(g)+ " ")
      cf.write("{:<g}".format(nu_delayed[g]))
      cf.write("\n")
    cf.write("NU_DELAYED_END"+"\n")

    cf.write("PRECURSOR_LAMBDA_BEGIN"+"\n")
    for j in range(0,J):
      cf.write("{:<4d}".format(j)+ " ")
      cf.write("{:<g}".format(decay_const[j]))
      cf.write("\n")
    cf.write("PRECURSOR_LAMBDA_END"+"\n")

    cf.write("PRECURSOR_FRACTION_BEGIN"+"\n")
    for j in range(0,J):
      cf.write("{:<4d}".format(j)+ " ")
      cf.write("{:<g}".format(precursor_fraction[j]))
      cf.write("\n")
    cf.write("PRECURSOR_FRACTION_END"+"\n")

    cf.write("PRECURSOR_GAMMA_BEGIN"+"\n")
    for j in range(0,J):
      cf.write("{:<4d}".format(j)+ " ")
      cf.write("{:<g}".format(gamma[j]))
      cf.write("\n")
    cf.write("PRECURSOR_GAMMA_END"+"\n")

    cf.write("CHI_DELAYED_BEGIN"+"\n")
    for g in range(0,G):
      for j in range(0,J):
        cf.write("G_PRECURSORJ_VAL"+ " ")
        cf.write("{:<4d}".format(g)+ " ")
        cf.write("{:<4d}".format(j)+ " ")
        cf.write("{:<g}".format(chi_delayed[g][j]))
        cf.write("\n")
    cf.write("CHI_DELAYED_END"+"\n")

  cf.close()

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

  # plt.figure(figsize=(6,6))
  # plt.plot(neutron_group_bndries, neutron_bndry_edge_values)
  # plt.yscale("log")
  # plt.xlabel("Energy (MeV)")
  # plt.ylabel("$\phi(E)$")
  # plt.show()

  if (gamma_bndry_edge_values != []):
    plt.figure(figsize=(6,6))
    plt.plot(gamma_group_bndries, gamma_bndry_edge_values)
    plt.yscale("log")
    plt.show()

#####################################################################
# Stand-alone usage
if __name__ == "__main__":
  from shutil import copyfile
  
  njoy_dir = "njoy_xs"
  chi_dir = "chi_xs"

  assert os.path.isdir(njoy_dir)
  assert os.path.isdir(chi_dir)

  mt1 = np.loadtxt('graphite_xs/graphite_mt1.txt')
  mt2 = np.loadtxt('graphite_xs/graphite_mt2.txt')
  mt230 = np.loadtxt('graphite_xs/graphite_mt230.txt')
  mt229 = np.loadtxt('graphite_xs/graphite_mt229.txt')

  plt.close('all')
  for njoy_file in os.listdir(njoy_dir):
    if njoy_file.endswith(".txt"):
      njoy_path = os.path.join(njoy_dir, njoy_file)

      # Parse NJOY cross section files
      print("\nPARSING FILE: " + njoy_path)
      with open(njoy_path, 'r') as xs_file:
        raw_njoy_data = ReadNJOYfile(njoy_path)
    

      # Reformat raw NJOY data
      data = BuildCombinedData(raw_njoy_data,plot=False)

      # Write to ChiTech cross section file
      chi_file = njoy_file.replace(".txt",".csx")
      chi_path = os.path.join(chi_dir, chi_file)
      WriteChiTechFile(data, chi_path)

  plt.show()
  # print(raw_njoy_data["transfer_matrices"]["(n,elastic)"][0])
  # print(raw_njoy_data["transfer_matrices"]["coherent"][11])
  # print(raw_njoy_data["transfer_matrices"]["incoherent"][11])
  # print(raw_njoy_data["transfer_matrices"]["pair_production"][8])