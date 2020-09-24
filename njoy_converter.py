''' Converts neutron gamma njoy output to chitech-format '''
import sys
import os
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#====================================================================
def Format1DCrossSection(xs,gs):
  '''Formats a 1D group average cross-section as 
     a blocky bargraph type plottable array '''
  
  gsr = gs[::-1]

  xs_bar = []
  e_bar  = []
  
  g=0
  for gsv in gsr:
    xs_bar.append(xs[g])
    e_bar.append(gsr[g][2]/1e6)

    xs_bar.append(xs[g])
    e_bar.append(gsr[g][1]/1e6)
    g+=1
  
  return np.array(xs_bar), e_bar


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
    for entry in total_nubar_data:
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
    nu_bar_delayed = np.mean(nu_delayed)
    delayed_frac = np.sum(chi_delayed,axis=0)
    gamma = nu_bar_delayed*delayed_frac

  # Normalize chi_delayed spectrum to sum to 1
  chi_delayed /= np.sum(chi_delayed,axis=0)

  # ================================= Combine transfer matrices

  n_to_n_elastic_keys = ["(n,elastic)"]
  n_to_n_elastic_sab_keys = ["mt230", "mt236"]
  n_to_n_elastic_freegas_keys = ["mt221"]
  n_to_n_inelastic_sab_keys = ["mt222", "mt229", "mt235"]
  n_to_n_inelastic_keys = ["(n,2n)"]
  n_to_g_transfer_keys = [ "(n,g)", "(n,inel)", "(n,np)", 
                           "(n,nd)", "(n,p)", "(n,d)", 
                           "(n,t)", "(n,a)"]
  g_to_g_transfer_keys = ["(g,coherent)", "(g,incoherent)", \
                          "(g,pair_production)"]

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
  nranges_to_nranges_elastic_freegas = []
  for rxn in n_to_n_elastic_freegas_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_elastic_freegas.append(transfer_matrices[rxn])

  # Adding all the elastic scattering S(\alpha,\beta) data
  nranges_to_nranges_elastic_sab = []
  for rxn in n_to_n_elastic_sab_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_elastic_sab.append(transfer_matrices[rxn])

  # Adding all the inelastic scattering S(\alpha,\beta) data
  nranges_to_nranges_inelastic_sab = []
  for rxn in n_to_n_inelastic_sab_keys:
    if rxn in transfer_matrices:
      nranges_to_nranges_inelastic_sab.append(transfer_matrices[rxn])

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
  for range_data in nranges_to_nranges_elastic_freegas:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_elastic_sab:
    if range_data:
      max_num_moms = max(max_num_moms,len(range_data[0])-2)
  for range_data in nranges_to_nranges_inelastic_sab:
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
  transfer_mats_sab      = np.copy(transfer_mats)

  # Regular elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_elastic:
    AddTransferNeutron(range_data)
  sig_s_el = np.sum(transfer_mats[0],axis=1)
  transfer_mats_standard += transfer_mats

  # Regular inelastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_inelastic:
    AddTransferNeutron(range_data)
  sig_s_inel = np.sum(transfer_mats[0],axis=1)
  transfer_mats_standard += transfer_mats

  # Regular n,\gamma transfer
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_granges:
    AddTransferNeutron(range_data,offset=G_g)
  transfer_mats_standard += transfer_mats

  # Regular \gamma,\gamma transfer
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in granges_to_granges:
    AddTransferGamma(range_data)
  transfer_mats_standard += transfer_mats
  
  # Freegas elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_elastic_freegas:
    AddTransferNeutron(range_data)
  sig_s_el_freegas = np.sum(transfer_mats[0],axis=1)
  transfer_mats_freegas += transfer_mats

  # S(alpha,beta) elastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_elastic_sab:
    AddTransferNeutron(range_data)
  sig_s_el_sab = np.sum(transfer_mats[0],axis=1)
  transfer_mats_sab_el = np.copy(transfer_mats)
  
  # S(alpha,beta) inelastic scatter
  transfer_mats = [0.0*mat for mat in transfer_mats]
  for range_data in nranges_to_nranges_inelastic_sab:
    AddTransferNeutron(range_data)
  sig_s_inel_sab = np.sum(transfer_mats[0],axis=1)
  transfer_mats_sab_inel = np.copy(transfer_mats)

  for m in range(max_num_moms):
    transfer_mats_sab[m] = transfer_mats_sab_inel[m] + transfer_mats_sab_el[m]

  # ===== Store stuff
  # Uncorrected quantities
  sig_t_uncorr = sig_t
  sig_s_uncorr = np.sum(transfer_mats_standard[0],axis=1)
  sig_a = sig_t_uncorr - sig_s_uncorr
  # Correction terms
  sig_s_freegas = sig_s_el_freegas
  sig_s_sab = sig_s_el_sab + sig_s_inel_sab
  
  # ===== Make the termal scattering corrections
  # This is a bit complex. The S(\alpha,\beta) thermal scattering
  # corrections should replace the standard transfer matrix where 
  # there is overlap, and in all other places, simply be placed
  # into the transfer matrix.
  transfer_mats = np.copy(transfer_mats_standard)
  for m in range(0,max_num_moms):
    m_transfer_mat_sab = transfer_mats_sab[m]
    # for nz in np.nonzero(m_transfer_mat_sab):
    #   transfer_mats[m][nz] = m_transfer_mat_sab[nz]
    for gprime in range(G):
      for g in range(G):
        if (m_transfer_mat_sab[gprime,g] > 1.0e-18):
          transfer_mats[m][gprime,g] = m_transfer_mat_sab[gprime,g]

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
    Atest = transfer_mats_sab_inel[0]
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
    plt.show()

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
  return_data["sigma_s_el_sab"] = sig_s_el_sab
  return_data["sigma_s_inel_sab"] = sig_s_inel_sab
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
    cf.write("PRECURSOR_LAMBDA_BEGIN"+"\n")
    for j in range(0,J):
      cf.write("{:<4d}".format(j)+ " ")
      cf.write("{:<g}".format(decay_const[j]))
      cf.write("\n")
    cf.write("PRECURSOR_LAMBDA_END"+"\n")

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

elasto=np.matrix([[1.0E-5,0.0],
[1.0625E-5,0.0],
[1.125E-5,0.0],
[1.1875E-5,0.0],
[1.25E-5,0.0],
[1.375E-5,0.0],
[1.5E-5,0.0],
[1.625E-5,0.0],
[1.75E-5,0.0],
[1.875E-5,0.0],
[2.0E-5,0.0],
[2.1875E-5,0.0],
[2.375E-5,0.0],
[2.5625E-5,0.0],
[2.75E-5,0.0],
[2.9375E-5,0.0],
[3.125E-5,0.0],
[3.3125E-5,0.0],
[3.5E-5,0.0],
[3.875E-5,0.0],
[4.25E-5,0.0],
[4.625E-5,0.0],
[5.0E-5,0.0],
[5.3125E-5,0.0],
[5.625E-5,0.0],
[5.9375E-5,0.0],
[6.25E-5,0.0],
[6.875E-5,0.0],
[7.5E-5,0.0],
[8.125E-5,0.0],
[8.75E-5,0.0],
[9.375E-5,0.0],
[1.0E-4,0.0],
[1.0625E-4,0.0],
[1.125E-4,0.0],
[1.1875E-4,0.0],
[1.25E-4,0.0],
[1.375E-4,0.0],
[1.5E-4,0.0],
[1.625E-4,0.0],
[1.75E-4,0.0],
[1.875E-4,0.0],
[2.0E-4,0.0],
[2.1875E-4,0.0],
[2.375E-4,0.0],
[2.5625E-4,0.0],
[2.75E-4,0.0],
[2.9375E-4,0.0],
[3.125E-4,0.0],
[3.3125E-4,0.0],
[3.5E-4,0.0],
[3.875E-4,0.0],
[4.25E-4,0.0],
[4.555488E-4,0.0],
[4.55549E-4,0.0],
[4.625E-4,0.0],
[5.0E-4,0.0],
[5.3125E-4,0.0],
[5.625E-4,0.0],
[5.9375E-4,0.0],
[6.25E-4,0.0],
[6.875E-4,0.0],
[7.5E-4,0.0],
[8.125E-4,0.0],
[8.75E-4,0.0],
[9.375E-4,0.0],
[0.001,0.0],
[0.0010625,0.0],
[0.001125,0.0],
[0.0011875,0.0],
[0.00125,0.0],
[0.001375,0.0],
[0.0015,0.0],
[0.001625,0.0],
[0.00175,0.0],
[0.001822195,0.0],
[0.001822197,7.394727],
[0.001875,7.18648],
[0.002,6.737325],
[0.0021875,6.15984],
[0.002375,5.673537],
[0.0025625,5.2584],
[0.00275,4.899873],
[0.0029375,4.587115],
[0.003125,4.311888],
[0.0033125,4.067819],
[0.0035,3.8499],
[0.003875,3.477329],
[0.004099939,3.286549],
[0.004099941,3.286547],
[0.00425,3.170506],
[0.004515511,2.984081],
[0.004515513,3.331187],
[0.004625,3.252329],
[0.00497106,3.025918],
[0.004971062,4.821734],
[0.005,4.793828],
[0.0053125,4.511838],
[0.005625,4.26118],
[0.0059375,4.036908],
[0.00625,3.835062],
[0.006337706,3.78199],
[0.006337708,4.192849],
[0.006875,3.865171],
[0.007288781,3.645747],
[0.007288783,4.526682],
[0.0075,4.3992],
[0.008125,4.0608],
[0.008615451,3.829631],
[0.008615453,4.591872],
[0.00875,4.521264],
[0.009375,4.219846],
[0.01,3.956106],
[0.010625,3.723394],
[0.01125,3.516539],
[0.01138871,3.473709],
[0.01138873,3.473703],
[0.01180428,3.351417],
[0.0118043,3.505454],
[0.011875,3.484584],
[0.0125,3.310354],
[0.01354653,3.054615],
[0.01354655,4.041785],
[0.01375,3.981981],
[0.01400207,3.910296],
[0.01400209,3.910291],
[0.015,3.650149],
[0.01536872,3.562576],
[0.01536874,5.170426],
[0.01590422,4.996342],
[0.01590424,5.281365],
[0.01625,5.168991],
[0.01639975,5.121791],
[0.01639977,5.362697],
[0.0175,5.025542],
[0.01764647,4.983829],
[0.01764649,4.983823],
[0.01806204,4.869162],
[0.01806206,4.907669],
[0.01851759,4.786941],
[0.01851761,5.008652],
[0.01875,4.946574],
[0.01988423,4.664413],
[0.01988425,4.730033],
[0.02,4.702658],
[0.02083531,4.514123],
[0.02083533,5.484886],
[0.02091526,5.463924],
[0.02091528,5.5242],
[0.021325,5.418063],
[0.02216198,5.213442],
[0.022162,5.377434],
[0.02232189,5.338916],
[0.02232191,5.338911],
[0.02265,5.261576],
[0.023975,4.97079],
[0.02493525,4.779367],
[0.02493527,4.779363],
[0.0253,4.710462],
[0.02535082,4.70102],
[0.02535084,4.744462],
[0.0268374,4.48166],
[0.02683742,4.599764],
[0.02684375,4.598679],
[0.0283875,4.348597],
[0.02915512,4.234104],
[0.02915514,4.324946],
[0.02945076,4.281533],
[0.02945078,4.381935],
[0.02993125,4.311594],
[0.02994629,4.309429],
[0.02994631,4.829403],
[0.031475,4.594847],
[0.03160857,4.57543],
[0.03160859,4.604963],
[0.03206412,4.539541],
[0.03206414,4.712286],
[0.03343077,4.519651],
[0.03343079,4.57309],
[0.03367063,4.540515],
[0.03367065,4.566891],
[0.0344618,4.462048],
[0.03446182,4.487346],
[0.0345625,4.474274],
[0.03570851,4.330679],
[0.03570853,4.473026],
[0.03586842,4.453087],
[0.03586844,4.453085],
[0.03689945,4.328661],
[0.03689947,4.328658],
[0.03765,4.242369],
[0.03889736,4.106325],
[0.03889738,4.146904],
[0.04038393,3.994255],
[0.04038395,4.051048],
[0.0406396,4.025564],
[0.04063962,4.175248],
[0.0407375,4.165216],
[0.04109515,4.128967],
[0.04109517,4.128965],
[0.04141496,4.097082],
[0.04141498,4.151273],
[0.04246179,4.048932],
[0.04246181,4.324787],
[0.04270166,4.300496],
[0.04270168,4.57345],
[0.0429973,4.542006],
[0.04299732,4.643045],
[0.043825,4.555357],
[0.04473954,4.462239],
[0.04473956,4.462237],
[0.04555488,4.382374],
[0.0455549,4.422635],
[0.0469125,4.294649],
[0.04721717,4.266937],
[0.04721719,4.281036],
[0.04792838,4.217512],
[0.0479284,4.436741],
[0.04800833,4.429354],
[0.04800835,4.456668],
[0.05,4.279146],
[0.05007039,4.27313],
[0.05007041,4.285721],
[0.053125,4.0393],
[0.05393047,3.978971],
[0.05393049,4.131209],
[0.0549615,4.053712],
[0.05496152,4.085176],
[0.05512141,4.073326],
[0.05512143,4.402868],
[0.05625,4.314532],
[0.05828607,4.163815],
[0.05828609,4.163813],
[0.05870164,4.134338],
[0.05870166,4.14353],
[0.05910142,4.115503],
[0.05910144,4.315005],
[0.059375,4.295124],
[0.05963692,4.27626],
[0.05963694,4.478642],
[0.0625,4.27348],
[0.06280158,4.252958],
[0.0628016,4.301045],
[0.06296149,4.290122],
[0.06296151,4.290121],
[0.06361693,4.245922],
[0.06361695,4.253725],
[0.06559903,4.125198],
[0.06559905,4.202343],
[0.06875,4.009741],
[0.06979473,3.949721],
[0.06979475,4.091154],
[0.07011454,4.072495],
[0.07011456,4.221066],
[0.07407038,3.995635],
[0.0740704,4.029299],
[0.075,3.979357],
[0.07698775,3.876614],
[0.07698777,3.9744],
[0.08102354,3.776436],
[0.08102356,3.804389],
[0.08125,3.793786],
[0.08150327,3.781997],
[0.08150329,3.886708],
[0.08579472,3.692296],
[0.08579474,3.696396],
[0.08619449,3.679253],
[0.08619451,3.812859],
[0.0875,3.755971],
[0.08928758,3.680775],
[0.0892876,3.72049],
[0.09375,3.543398],
[0.09380309,3.541393],
[0.09380311,3.769374],
[0.09892568,3.574188],
[0.0989257,3.621197],
[0.1,3.582294],
[0.1014032,3.532723],
[0.1014034,3.629731],
[0.1024984,3.590955],
[0.1024986,3.692831],
[0.10625,3.562447],
[0.1081165,3.500946],
[0.1081167,3.515228],
[0.1085962,3.499707],
[0.1085964,3.513835],
[0.1091476,3.49609],
[0.1091478,3.598479],
[0.1125,3.491254],
[0.1147099,3.423995],
[0.1147101,3.430154],
[0.1166204,3.373966],
[0.1166206,3.44724],
[0.11875,3.385425],
[0.1226941,3.276598],
[0.1226943,3.314597],
[0.1239808,3.280202],
[0.123981,3.334171],
[0.1242764,3.326245],
[0.1242766,3.346237],
[0.125,3.326872],
[0.1273695,3.264981],
[0.1273697,3.299362],
[0.1299271,3.23442],
[0.1299273,3.29404],
[0.1316535,3.250849],
[0.1316537,3.290581],
[0.1375,3.150671],
[0.1383185,3.132026],
[0.1383187,3.152209],
[0.1399808,3.114781],
[0.139981,3.115997],
[0.1403805,3.10713],
[0.1403807,3.162842],
[0.1420428,3.125832],
[0.142043,3.129341],
[0.143138,3.105402],
[0.1431382,3.132891],
[0.1475977,3.038234],
[0.1475979,3.087047],
[0.15,3.037611],
[0.1555894,2.928487],
[0.1555896,2.930301],
[0.1563805,2.915481],
[0.1563807,2.917267],
[0.15726,2.900956],
[0.1572602,3.005347],
[0.1625,2.908439],
[0.1644531,2.873898],
[0.1644533,2.930301],
[0.173947,2.770371],
[0.1739472,2.771674],
[0.1747466,2.758994],
[0.1747468,2.817809],
[0.175,2.813732],
[0.1822195,2.702253],
[0.1822197,2.774525],
[0.1875,2.69639],
[0.1917135,2.637128],
[0.1917137,2.68278],
[0.2,2.571628],
[0.2002815,2.568014],
[0.2002817,2.568431],
[0.200897,2.560564],
[0.2008972,2.601305],
[0.2112063,2.474333],
[0.2112065,2.479965],
[0.2114462,2.477154],
[0.2114464,2.489077],
[0.213828,2.461354],
[0.2138282,2.462027],
[0.2144435,2.454962],
[0.2144437,2.45496],
[0.2153862,2.444217],
[0.2153864,2.459637],
[0.21875,2.421816],
[0.218959,2.419505],
[0.2189592,2.420748],
[0.2204856,2.40399],
[0.2204858,2.445674],
[0.2325056,2.31924],
[0.2325058,2.320755],
[0.2326487,2.31933],
[0.2326489,2.33692],
[0.2367012,2.296912],
[0.2367014,2.302112],
[0.2375,2.294371],
[0.2385476,2.284295],
[0.2385478,2.286826],
[0.2409853,2.263695],
[0.2409855,2.297212],
[0.2530926,2.187321],
[0.2530928,2.187319],
[0.2544677,2.175501],
[0.2544679,2.183098],
[0.25625,2.167916],
[0.2590473,2.144506],
[0.2590475,2.157406],
[0.2623961,2.129874],
[0.2623963,2.157817],
[0.275,2.058921],
[0.2759017,2.052192],
[0.2759019,2.052992],
[0.2759426,2.052689],
[0.2759428,2.054822],
[0.2763814,2.051561],
[0.2763816,2.057255],
[0.2773409,2.050139],
[0.2773411,2.053918],
[0.2804581,2.031091],
[0.2804583,2.035944],
[0.284718,2.005484],
[0.2847182,2.028988],
[0.29375,1.966604],
[0.2996869,1.927644],
[0.2996871,1.928217],
[0.3003814,1.92376],
[0.3003816,1.927356],
[0.30278,1.912089],
[0.3027802,1.921966],
[0.307951,1.889695],
[0.3079512,1.906609],
[0.3125,1.878856],
[0.3237027,1.813832],
[0.3237029,1.816038],
[0.3237585,1.815727],
[0.3237587,1.816896],
[0.3253576,1.807968],
[0.3253578,1.808573],
[0.326013,1.804938],
[0.3260132,1.809459],
[0.33125,1.780853],
[0.3320951,1.776321],
[0.3320953,1.789373],
[0.3492934,1.70127],
[0.3492936,1.702938],
[0.35,1.699501],
[0.3501571,1.698738],
[0.3501573,1.701541],
[0.3571502,1.668225],
[0.3571504,1.682844],
[0.3751872,1.601943],
[0.3751874,1.602012],
[0.3752123,1.601906],
[0.3752125,1.606021],
[0.3831165,1.572887],
[0.3831167,1.581757],
[0.3875,1.563865],
[0.4022802,1.506407],
[0.4022804,1.506456],
[0.4023769,1.506095],
[0.4023771,1.506144],
[0.4026558,1.505101],
[0.402656,1.505895],
[0.4034795,1.502821],
[0.4034797,1.502894],
[0.4035437,1.502656],
[0.4035439,1.503048],
[0.4043432,1.500076],
[0.4043434,1.502353],
[0.4099939,1.481648],
[0.4099941,1.488994],
[0.425,1.43642],
[0.4313312,1.415336],
[0.4313314,1.415336],
[0.4313554,1.415257],
[0.4313556,1.417441],
[0.4365905,1.400446],
[0.4365907,1.400722],
[0.4377824,1.396909],
[0.4377826,1.402212],
[0.4597993,1.335069],
[0.4597995,1.335093],
[0.4603264,1.333565],
[0.4603266,1.333638],
[0.4608229,1.332201],
[0.4608231,1.332201],
[0.4612859,1.330864],
[0.4612861,1.330936],
[0.4622454,1.328174],
[0.4622456,1.329115],
[0.466482,1.317044],
[0.4664822,1.321694],
[0.4898897,1.258542],
[0.4898899,1.258781],
[0.4946535,1.246659],
[0.4946537,1.24708],
[0.4960043,1.243684],
[0.4960045,1.243707],
[0.4960927,1.243486],
[0.4960929,1.247523],
[0.5,1.237774],
[0.5211785,1.187476],
[0.5211787,1.187535],
[0.5228817,1.183667],
[0.5228819,1.183917],
[0.5266144,1.175526],
[0.5266146,1.178303],
[0.5532742,1.121526],
[0.5532744,1.121591],
[0.5533551,1.121428],
[0.5533553,1.121697],
[0.5580473,1.112266],
[0.5580475,1.114312],
[0.5860989,1.060979],
[0.5860991,1.061069],
[0.586242,1.06081],
[0.5862422,1.061034],
[0.5903913,1.053578],
[0.5903915,1.055102],
[0.6202241,1.004352],
[0.6202243,1.004356],
[0.6202883,1.004252],
[0.6202885,1.004252],
[0.6210878,1.002959],
[0.621088,1.003073],
[0.6236464,0.9989576],
[0.6236466,1.000051],
[0.625,0.997885],
[0.6552549,0.9518099],
[0.6552551,0.9518432],
[0.6575224,0.948561],
[0.6575226,0.9485809],
[0.6578125,0.9481629],
[0.6578127,0.9490805],
[0.691577,0.9027443],
[0.6915772,0.9027493],
[0.6925365,0.9014988],
[0.6925367,0.9015011],
[0.6928898,0.9010417],
[0.69289,0.901692],
[0.7278936,0.8583307],
[0.7278938,0.8583491],
[0.7285964,0.8575214],
[0.7285966,0.8575217],
[0.7288781,0.8571905],
[0.7288783,0.8576157],
[0.75,0.8334633],
[0.7657776,0.8162912],
[0.7657778,0.8166345],
[0.8051454,0.7767052],
[0.8051456,0.7767091],
[0.8057776,0.7760999],
[0.8057778,0.7761009],
[0.8064172,0.7754856],
[0.8064174,0.7754898],
[0.8081037,0.7738716],
[0.8081039,0.7740532],
[0.8423098,0.7426192],
[0.84231,0.7427902],
[0.875,0.7150395],
[0.8850402,0.7069279],
[0.8850404,0.7069371],
[0.8864581,0.7058065],
[0.8864583,0.7059043],
[0.9224864,0.6783349],
[0.9224866,0.6784087],
[0.9687442,0.6460146],
[0.9687444,0.6460147],
[0.9689589,0.6458717],
[0.9689591,0.6458809],
[0.9743278,0.642322],
[0.974328,0.6423255],
[0.9781656,0.6398055],
[0.9781658,0.6398082],
[0.9818825,0.6373863],
[0.9818827,0.6373864],
[0.9820034,0.6373081],
[0.9820036,0.6373312],
[1.0,0.6258615],
[1.006307,0.6219389],
[1.006309,0.6219775],
[1.056853,0.5922315],
[1.056855,0.5922305],
[1.058047,0.5915633],
[1.058049,0.5915629],
[1.05851,0.5913052],
[1.058512,0.5913042],
[1.058645,0.5912299],
[1.058647,0.5912326],
[1.063569,0.5884965],
[1.063571,0.5884963],
[1.067646,0.5862501],
[1.067648,0.5862608],
[1.093772,0.5722584],
[1.093774,0.5722762],
[1.148602,0.5449588],
[1.148604,0.5449579],
[1.14882,0.5448554],
[1.148822,0.5448548],
[1.150314,0.5441481],
[1.150316,0.5441473],
[1.152418,0.5431547],
[1.15242,0.5431542],
[1.15597,0.5414862],
[1.155972,0.5414854],
[1.156933,0.5410357],
[1.156935,0.5410403],
[1.184882,0.5282792],
[1.184884,0.5282857],
[1.244866,0.5028311],
[1.244868,0.5028303],
[1.245258,0.5026728],
[1.24526,0.5026721],
[1.247698,0.5016899],
[1.2477,0.5016894],
[1.249865,0.5008204],
[1.249867,0.500822],
[1.25,0.5007687],
[1.279636,0.4891711],
[1.279638,0.4891739],
[1.343921,0.4657755],
[1.343923,0.4657748],
[1.344691,0.4655088],
[1.344693,0.4655082],
[1.346442,0.4649035],
[1.346444,0.4649036],
[1.378034,0.4542461],
[1.378036,0.4542471],
[1.447426,0.4324703],
[1.447428,0.4324697],
[1.448865,0.4320408],
[1.448867,0.4320402],
[1.449984,0.4317074],
[1.449986,0.4317068],
[1.450283,0.4316184],
[1.450285,0.431618],
[1.46924,0.4260496],
[1.469242,0.4260491],
[1.480077,0.4229302],
[1.480079,0.4229303],
[1.5,0.4173135],
[1.55411,0.4027837],
[1.554112,0.4027832],
[1.554139,0.4027762],
[1.554141,0.4027757],
[1.554719,0.4026259],
[1.554721,0.4026254],
[1.555216,0.4024973],
[1.555218,0.4024968],
[1.561436,0.400894],
[1.561438,0.4008935],
[1.564074,0.4002178],
[1.564076,0.4002173],
[1.565523,0.3998474],
[1.565525,0.3998469],
[1.565964,0.3997348],
[1.565966,0.3997343],
[1.566709,0.3995447],
[1.566711,0.3995442],
[1.568629,0.3990557],
[1.568631,0.3990552],
[1.573105,0.3979202],
[1.573107,0.3979198],
[1.585765,0.3947435],
[1.585767,0.3947432],
[1.6652,0.3759133],
[1.665202,0.3759128],
[1.66544,0.3758591],
[1.665442,0.3758587],
[1.666998,0.3755078],
[1.667,0.3755074],
[1.667027,0.3755013],
[1.667029,0.3755008],
[1.680615,0.3724653   ],
[1.680617,0.3724649],
[1.695096,0.3692834],
[1.695098,0.3692831],
[1.75,0.3576977],
[1.779931,0.3516827],
[1.779933,0.3516823],
[1.791769,0.3493592],
[1.791771,0.3493588],
[2.0,0.3129855],
[2.1875,0.2861582],
[2.375,0.2635667],
[2.5625,0.2442814],
[2.75,0.2276258],
[3.125,0.2003107],
[3.5,0.1788489],
[3.875,0.1615409],
[4.25,0.1472873],
[4.625,0.1353451],
[4.999999,0.1251942],
[5.0,0.1251942],
[5.000001,0.0],
[2.0E7,0.0]])

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
      if "Cnat" in njoy_file and "172g" in njoy_file:
        print("\nPARSING FILE: " + njoy_path)
        with open(njoy_path, 'r') as xs_file:
          raw_njoy_data = ReadNJOYfile(njoy_path)
      

        # Reformat raw NJOY data
        data = BuildCombinedData(raw_njoy_data,plot=True)

        # Write to ChiTech cross section file
        chi_file = njoy_file.replace(".txt",".csx")
        chi_path = os.path.join(chi_dir, chi_file)
        WriteChiTechFile(data, chi_path)

        plt.figure(figsize=(12,6))
        plt.semilogx(elasto[:,0]*1e-6,elasto[:,1],"-r",label=r"MT230 Janis")

        xs = data["sigma_t"]
        gs = data["neutron_gs"]
        xs, gs = Format1DCrossSection(xs, gs)
        plt.semilogx(gs,xs,"-k",linewidth=0.5,label=r"MT1")

        xs = data["sigma_s_el"]
        gs = data["neutron_gs"]
        xs, gs = Format1DCrossSection(xs, gs)
        plt.semilogx(gs,xs,"-b",linewidth=0.5,label=r"MT2")

        xs = data["sigma_s_inel_sab"]
        gs = data["neutron_gs"]
        xs229, gs = Format1DCrossSection(xs, gs)
        # plt.semilogx(gs,xs229,"-y",linewidth=0.5,label=r"MT229")

        xs = data["sigma_s_el_sab"]
        gs = data["neutron_gs"]
        xs230, gs = Format1DCrossSection(xs, gs)
        # plt.semilogx(gs,xs230,"-b",linewidth=0.5,label=r"MT230")

        plt.semilogx(gs,xs229+xs230,"-g",linewidth=0.5,label=r"MT229+MT230")

        plt.xlabel("Energy (Mev)")
        plt.ylabel("Cross section (b)")
        plt.legend()
        plt.grid(True)
        plt.ylim([0.0,10.0])  
        plt.show()

  # plt.show()
  # print(raw_njoy_data["transfer_matrices"]["(n,elastic)"][0])
  # print(raw_njoy_data["transfer_matrices"]["coherent"][11])
  # print(raw_njoy_data["transfer_matrices"]["incoherent"][11])
  # print(raw_njoy_data["transfer_matrices"]["pair_production"][8])