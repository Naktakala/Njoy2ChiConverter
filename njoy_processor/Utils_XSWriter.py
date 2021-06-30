"""Writes cross-sections to file"""
import numpy as np


# ===================================================================
def WriteChiTechFile(data, chi_filename = "output.cxs", comment = "# Output"):
    cf = open(chi_filename, 'w')

    sig_t = data["sigma_t"]
    sig_a = data["sigma_a"]
    sig_s = data["sigma_s"]
    sig_f = data["sigma_f"]
    sig_heat = data["sigma_heat"]
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

    cf.write(comment + "\n")
    cf.write("NUM_GROUPS " + str(G) + "\n")
    cf.write("NUM_MOMENTS " + str(M) + "\n")
    if J > 0:
        cf.write("NUM_PRECURSORS " + str(J) + "\n")

    cf.write("\n")

    cf.write("SIGMA_T_BEGIN" + "\n")
    for g in range(0, G):
        cf.write("{:<4d}".format(g) + " ")
        cf.write("{:<g}".format(sig_t[g]))
        cf.write("\n")
    cf.write("SIGMA_T_END" + "\n\n")
    
    cf.write("SIGMA_A_BEGIN" + "\n")
    for g in range(0, G):
        cf.write("{:<4d}".format(g) + " ")
        cf.write("{:<g}".format(sig_a[g]))
        cf.write("\n")
    cf.write("SIGMA_A_END" + "\n\n")

    cf.write("SIGMA_S_BEGIN" + "\n")
    for g in range(0, G):
        cf.write("{:<4d}".format(g) + " ")
        cf.write("{:<g}".format(sig_s[g]))
        cf.write("\n")
    cf.write("SIGMA_S_END" + "\n\n")
    
    cf.write("SIGMA_HEAT_BEGIN" + "\n")
    for g in range(0, G):
        cf.write("{:<4d}".format(g) + " ")
        cf.write("{:<g}".format(sig_heat[g]))
        cf.write("\n")
    cf.write("SIGMA_HEAT_END" + "\n\n")
    
    if np.linalg.norm(sig_f) > 1.0e-20:
        cf.write("SIGMA_F_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(sig_f[g]))
            cf.write("\n")
        cf.write("SIGMA_F_END" + "\n\n")

    if np.linalg.norm(nu_total) > 1.0e-20:
        cf.write("NU_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(nu_total[g]))
            cf.write("\n")
        cf.write("NU_END" + "\n\n")

    if np.linalg.norm(nu_prompt) > 1.0e-20:
        cf.write("NU_PROMPT_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(nu_prompt[g]))
            cf.write("\n")
        cf.write("NU_PROMPT_END" + "\n\n")

    if np.linalg.norm(chi_prompt) > 1.0e-20:
        cf.write("CHI_PROMPT_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(chi_prompt[g]))
            cf.write("\n")
        cf.write("CHI_PROMPT_END" + "\n\n")

    if np.linalg.norm(ddt_coeff) > 1.0e-20:
        cf.write("DDT_COEFF_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(ddt_coeff[g] / 100))
            cf.write("\n")
        cf.write("DDT_COEFF_END" + "\n\n")

    cf.write("TRANSFER_MOMENTS_BEGIN" + "\n")
    for m in range(0, M):
        cf.write("# l = " + str(m) + "\n")
        for gprime in range(0, G):
            for g in transfer_mats_nonzeros[m][gprime]:
                cf.write("M_GPRIME_G_VAL" + " ")
                cf.write("{:<4d}".format(m) + " ")
                cf.write("{:<4d}".format(gprime) + " ")
                cf.write("{:<4d}".format(g) + " ")
                cf.write("{:<g}".format(transfer_mats[m][gprime, g]))
                cf.write("\n")
    cf.write("TRANSFER_MOMENTS_END" + "\n\n")

    if J > 0:
        cf.write("PRECURSOR_LAMBDA_BEGIN" + "\n")
        for j in range(0, J):
            cf.write("{:<4d}".format(j) + " ")
            cf.write("{:<g}".format(decay_const[j]))
            cf.write("\n")
        cf.write("PRECURSOR_LAMBDA_END" + "\n\n")

        cf.write("PRECURSOR_GAMMA_BEGIN" + "\n")
        for j in range(0, J):
            cf.write("{:<4d}".format(j) + " ")
            cf.write("{:<g}".format(gamma[j]))
            cf.write("\n")
        cf.write("PRECURSOR_GAMMA_END" + "\n\n")

        cf.write("NU_DELAYED_BEGIN" + "\n")
        for g in range(0, G):
            cf.write("{:<4d}".format(g) + " ")
            cf.write("{:<g}".format(nu_delayed[g]))
            cf.write("\n")
        cf.write("NU_DELAYED_END" + "\n\n")

        cf.write("CHI_DELAYED_BEGIN" + "\n")
        for g in range(0, G):
            for j in range(0, J):
                cf.write("G_PRECURSORJ_VAL" + " ")
                cf.write("{:<4d}".format(g) + " ")
                cf.write("{:<4d}".format(j) + " ")
                cf.write("{:<g}".format(chi_delayed[g][j]))
                cf.write("\n")
        cf.write("CHI_DELAYED_END" + "\n\n")

    cf.close()
