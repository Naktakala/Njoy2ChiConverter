"""Contains all the necessary functionality to
   read NJOY cross-section files."""


# ===================================================================
def PrintLine(line_num, line):
    """ Prints a file line without the line feed. """
    print("Line: ", line_num, line, end = '')


# ===================================================================
def ProcessGroupStructure(nL_i, lines):
    """ Processes a group structure.
      params:
        nL_i  File line number before data starts.
        lines An array of files lines.
        
      return:
        returns the group_structure """
    group_struct = []
    nL = nL_i + 1
    words = lines[nL].split()
    while True:
        if len(words) == 0:
            nL += 1
            words = lines[nL].split()
        if words[0] == "1":
            break
        nL += 1
        words = lines[nL].split()

    while True:
        words = lines[nL].split()
        if not words:
            break

        group_struct.append([int(words[0]) - 1,
                             float(words[1]),
                             float(words[3])])
        nL += 1
    return group_struct


# ===================================================================
def ProcessCrossSection(nL_i, lines, header_size = 5, line_incr = 1, weighting_spectrum=False):
    """ Reads 1D neutron cross-sections from the lines.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
      
      return:
        A table containing the group wise xs. """
    xs = []
    # Skip header
    add_skip = 0
    if lines[nL_i + 1].find("particle emission") >= 0:
        add_skip = 1

    nL = nL_i + header_size + add_skip
    words = lines[nL].split()
    num_words = len(words)

    while num_words >= 2:
        number_word = words[1]
        number_word = number_word.replace("+", "E+")
        number_word = number_word.replace("-", "E-")
        if not weighting_spectrum:
            xs.append([int(words[0]) - 1, float(number_word)])
        else:
            if words[0]=='flx':
                xs.append([float(number_word)])
        nL += line_incr

        words = lines[nL].split()
        num_words = len(words)

    return xs


# ===================================================================
def ProcessPromptChi(nL_i, lines, header_size = 4, line_incr = 1):
    """ Reads 1D prompt fission spectrum from the lines.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
      
      return:
        A table containing the group wise xs. """
    xs = []
    # Skip header
    nL = nL_i + header_size
    words = lines[nL].split()
    num_words = len(words)

    g = 0
    while num_words >= 2:
        if not words[0].isdigit():
            nL += 2
            words = lines[nL].split()
            num_words = len(words)

        for number_word in words[1:]:
            number_word = number_word.replace("+", "E+")
            number_word = number_word.replace("-", "E-")
            xs.append([g, float(number_word)])
            g += 1

        nL += line_incr
        words = lines[nL].split()
        num_words = len(words)
    return xs


# ===================================================================
def ProcessDecayConstants(nL_i, lines):
    """ Reads delayed neutron decay constants from delayed chi.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
    
      return:
        A table containing delayed neutron group decay constants.
    """
    decay_const = []
    nL = nL_i + 2
    nL += 1 if lines[nL].split()[0] == "normalized" else 3
    words = lines[nL].split()
    for j, number_word in enumerate(words[1:]):
        number_word = number_word.replace("+", "E+")
        number_word = number_word.replace("-", "E-")
        decay_const.append([j, float(number_word)])
    return decay_const


# ===================================================================
def ProcessDelayedChi(nL_i, lines):
    """ Reads the delayed neutron spectra from the lines.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
      
      return:
        A table containing the group wise delayed 
        neutron precursor spectrum coefficients. """
    matrix = []
    nL = nL_i + 2
    nL += 3 if lines[nL].split()[0] == "normalized" else 5
    words = lines[nL].split()
    num_words = len(words)
    while num_words >= 2:
        g = int(words[0]) - 1
        entry = [g]
        for j, number_word in enumerate(words[1:]):
            number_word = number_word.replace("+", "E+")
            number_word = number_word.replace("-", "E-")
            entry.append(float(number_word))
        matrix.append(entry)

        nL += 1
        words = lines[nL].split()
        num_words = len(words)

    return matrix


# ===================================================================
def ProcessTransferMatrix(nL_i, lines):
    """ Reads transfer matrices from the lines.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
      
      return:
        A table containing the group wise transfer coefficients. """
    matrix = []
    # Skip 4 lines
    add_skip = 0
    if lines[nL_i + 1].find("particle emission") >= 0 or \
       lines[nL_i + 2].find("spectrum constant below") >= 0:
        add_skip = 1
    nL = nL_i + 5 + add_skip
    words = lines[nL].split()
    num_words = len(words)

    while num_words >= 2:
        # Transform words 0.000+00 to 0.000E+00
        for i in range(0, num_words):
            word = words[i]
            loc_of_sign = word.find("-", 2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+", 2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+", "E+")
            secnd_part = secnd_part.replace("-", "E-")

            words[i] = first_part + secnd_part

        # Develop table entry
        if words[0] != "normalization":
            entry = [int(words[0]) - 1, int(words[1]) - 1]
            for i in range(2, num_words):
                entry.append(float(words[i]))
            matrix.append(entry)

        nL += 1
        words = lines[nL].split()
        num_words = len(words)
    return matrix


# ===================================================================
def ProcessTransferMatrixB(nL_i, lines, process_overflow = False):
    """ Reads transfer matrices from the lines.
      params:
        nL_i    File line number before data starts.
        lines   An array of file lines.
        line_incr Line increment to use
      
      return:
        A table containing the group wise transfer coefficients. """
    matrix = []
    # Skip 3 lines
    nL = nL_i + 4
    words = lines[nL].split()
    num_words = len(words)

    while num_words >= 2:
        # Transform words 0.000-00 to 0.000E+00
        for i in range(0, num_words):
            word = words[i]
            loc_of_sign = word.find("-", 2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+", 2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+", "E+")
            secnd_part = secnd_part.replace("-", "E-")

            words[i] = first_part + secnd_part

        # Develop table entry
        entry = [int(words[0]) - 1, int(words[1]) - 1]
        for i in range(2, num_words):
            entry.append(float(words[i]))

        # Process overflow
        if process_overflow:
            words = lines[nL + 1].split()

            word = words[0]
            loc_of_sign = word.find("-", 2)
            if loc_of_sign == -1:
                loc_of_sign = word.find("+", 2)
            first_part = word[0:loc_of_sign]
            secnd_part = word[loc_of_sign:]

            secnd_part = secnd_part.replace("+", "E+")
            secnd_part = secnd_part.replace("-", "E-")

            words[0] = first_part + secnd_part

            entry.append(float(words[0]))
            nL += 1

        # Complete matrix entry
        matrix.append(entry)

        # Move to next line
        nL += 1
        words = lines[nL].split()
        num_words = len(words)

        # Skip lines with xsec and heat
        while num_words >= 2 and not words[1].isnumeric():
            nL += 1
            words = lines[nL].split()
            num_words = len(words)

    return matrix


# ===================================================================
def ReadNJOYfile(njoy_filename = "output", verbose = False):
    """ Reads an NJOY output file
      params:
        njoy_filename Name of the NJOY file to process. 
        
      return:
        Returns a complex dictionary of raw data. """

    njoy_file = open(njoy_filename, 'r')
    file_lines = njoy_file.readlines()

    njoy_raw_data = {}

    group_structures = {}
    cross_sections = {}
    transfer_matrices = {}
    transfer_matrices['neutron'] ={}
    transfer_matrices['gamma'] ={}
    weighting_spectrum = {}

    # flag_run_processed             = False
    flag_gamma_structure_processed = False

    nL = -1
    while nL < (len(file_lines) - 1):
        nL += 1
        line = file_lines[nL]
        words = file_lines[nL].split()
        num_words = len(words)

        # if (line.find("neutron group structure") != -1 and words[3] != "option"):
        #   group_structures["neutron"] = ProcessGroupStructure(nL,file_lines)

        if line.find("sigma zeroes") != -1:
            group_structures["neutron"] = ProcessGroupStructure(nL, file_lines)

        if line.find("gamma group structure") != -1:
            if not flag_gamma_structure_processed:
                group_structures["gamma"] = ProcessGroupStructure(nL, file_lines)
                flag_gamma_structure_processed = True

###---------lines containing both mf and mt -----------------------------------
        if line.find("for mf") != -1 and line.find("mt") != -1:

###---------cross section -----------------------------------------------------
            if words[num_words - 2] == "cross" and words[num_words - 1] == "section":
                # incr=2 for n,total
                if words[2] == "3" and words[5] == "1":
                    cross_sections[words[num_words - 3]] = \
                        ProcessCrossSection(nL, file_lines, line_incr = 2)
                    # this is to get the weighting spectrum
                    weighting_spectrum[words[num_words - 3]] = \
                        ProcessCrossSection(nL, file_lines, weighting_spectrum=True)
                else:
                    cross_sections[words[num_words - 3]] = \
                        ProcessCrossSection(nL, file_lines)

###---------other MF 3 --------------------------------------------------------
            if words[2] == "3" and words[4] == "mt259":
                cross_sections["inv_velocity"] = \
                    ProcessCrossSection(nL, file_lines)

            if words[2] == "3" and words[4] == "mt221":
                cross_sections["free_gas"] = \
                    ProcessCrossSection(nL, file_lines)

            if words[2] == "3" and words[4] == "mt452":
                cross_sections["total_nubar"] = \
                    ProcessCrossSection(nL, file_lines)

            if words[2] == "3" and words[4] == "mt456":
                cross_sections["prompt_nubar"] = \
                    ProcessCrossSection(nL, file_lines)

            if words[2] == "3" and words[4] == "mt455":
                cross_sections["delayed_nubar"] = \
                    ProcessCrossSection(nL, file_lines)

###---------MF 5 --------------------------------------------------------------
            if words[2] == "5" and words[5] == "18":
                cross_sections["prompt_chi"] = \
                    ProcessPromptChi(nL, file_lines, 4)

            if words[2] == "5" and words[4] == "mt455":
                cross_sections["decay_constants"] = \
                    ProcessDecayConstants(nL, file_lines)
                cross_sections["delayed_chi"] = \
                    ProcessDelayedChi(nL, file_lines)

###---------matrix ------------------------------------------------------------
            # caveat: the pp values from transfer(mf26) = 2x the pp from the xsec(mf23)
            # caveat: the n,2n values from transfer(mf8/mt16) = 2x the n,2n from the xsec(mf3/mt16)
            if words[num_words - 1] == "matrix":
                particle_type = words[num_words - 2]
                reaction_type  = words[num_words - 3]
                if particle_type=="free-gas":
                    particle_type = "neutron"
                if particle_type=="inelastic_s(a,b)":
                    particle_type = "neutron"
                if particle_type=="elastic_s(a,b)":
                    particle_type = "neutron"
                transfer_matrices[particle_type][reaction_type] = \
                    ProcessTransferMatrix(nL, file_lines)

###---------MF 23 -------------------------------------------------------------
            # Total photon interaction
            if words[1] == "mf23" and words[3] == "mt501":
                cross_sections["(g,total)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

            # Photon coherent scattering
            if words[1] == "mf23" and words[3] == "mt502":
                cross_sections["(g,coherent)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

            # Photon incoherent scattering
            if words[1] == "mf23" and words[3] == "mt504":
                cross_sections["(g,incoherent)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

            # 515: Pair production, electron field
            # 517: Pair production, nuclear field
            # 516: Pair production; sum of MT=515, 517.
            if words[1] == "mf23" and words[3] == "mt516":
                cross_sections["(g,pair_production)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

            # Photoelectric absorption
            if words[1] == "mf23" and words[3] == "mt522":
                cross_sections["(g,abst)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

            if words[1] == "mf23" and words[3] == "mt525":
                cross_sections["(g,heat)"] = \
                    ProcessCrossSection(nL, file_lines, header_size = 4)

###---------MF 26 -------------------------------------------------------------
            if words[1] == "mf26" and words[3] == "mt502":
                transfer_matrices['gamma']["(coherent)"] = \
                    ProcessTransferMatrixB(nL, file_lines,
                                           process_overflow = True)

            # we do not see to save the xsec(g), it is just the sum_k xs(g->k)
            # the incoh heat value is not saved either,
            # the total heat xs (mf23/mt525) is the sum of the heat from
            # inch(mf23/mt504) + abst(mf23/mt522) + pp(mf23/mt516)
            if words[1] == "mf26" and words[3] == "mt504":
                transfer_matrices['gamma']["(incoherent)"] = \
                    ProcessTransferMatrixB(nL, file_lines,
                                           process_overflow = True)

            # caveat: the pp values from transfer(mf26) = 2x the pp from the xsec(mf23)
            if words[1] == "mf26" and words[3] == "mt516":
                transfer_matrices['gamma']["(pair_production)"] = \
                    ProcessTransferMatrixB(nL, file_lines,
                                           process_overflow = False)

    njoy_file.close()

    njoy_raw_data["group_structures"] = group_structures
    njoy_raw_data["cross_sections"] = cross_sections
    njoy_raw_data["transfer_matrices"] = transfer_matrices

    if verbose:
        print("Cross-sections extracted:")
        xss = njoy_raw_data["cross_sections"]
        for k in xss:
            print("     " + k)

        print("Transfer matrices extracted")
        mats = njoy_raw_data["transfer_matrices"]
        for ptype in mats:
            print("     " + ptype)
            for k in mats[ptype]:
                print("        " + k)


    return njoy_raw_data
