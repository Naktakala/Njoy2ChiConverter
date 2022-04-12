import numpy as np 
import matplotlib.pyplot as plt
import Histogram_Source as Histo
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
    
        #heating_spectrum  = sigma_heat[G_neutron-g-1]
        heating_spectrum  = sigma_heat[gprime]*spectrum
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
    
        heating_spectrum  = sigma_heat[gprime]*spectrum
        g_heating += [heating_spectrum, heating_spectrum]
  
    g_bndrys = np.array(g_bndrys)
    g_vals   = np.array(g_vals)
    g_heating = np.array(g_heating)
  
    return [(n_bndrys, n_vals), (g_bndrys, g_vals), (n_heating, g_heating)]

#====================================================================
def InfiniteMediumSpectrum(data, source_def, plot=False):
    neutron_gs = data["neutron_gs"]
    gamma_gs = data["gamma_gs"]
    sig_t = data["sigma_t"]
    transfer_mats = data["transfer_matrices"]
    sig_heat = data["sigma_heat"]
  
    #======================================= Get data from dictionary
    G_neutron = len(neutron_gs)
    G_gamma   = len(gamma_gs) 
    G = np.size(sig_t)
    
    v_sig_t = np.array(sig_t)
    M_sig_gp_to_g = np.zeros([G,G])
    M_sig_gp_to_g[:,:] = transfer_mats[0][:,:]
  
    #======================================= Solve for psi
    S = M_sig_gp_to_g.transpose()
    T = np.diag(v_sig_t)
  
    if plot:
        plt.matshow(np.log(S))
        if G_gamma>0:
            g_txt = '_g'+str(G_gamma)
        else:
            g_txt = ''
        filename = 'transfert_matrix_n'+str(G_neutron)+g_txt+'.png'
        plt.savefig(filename)
  
    A = T - S
    A_inv = np.linalg.inv(A)
  
    v_src = np.zeros(G)
    # Get the source description
    particle = source_def["particle_type"].lower()
    if source_def["energy"]=='fission':
        # Get the source term:
        if particle == "neutron":
            src_term = create_source_spectrum(neutron_gs, 1., fission = True)
            for i in range (len(src_term)):
                v_src[i] = src_term[i]
        else:
            raise Exception('for fission, particle must be neutron')
    else:
        if particle == "neutron":
            src_term = create_source_spectrum(neutron_gs, source_def["energy"])
            for i in range (len(src_term)):
                v_src[i] = src_term[i]
        elif particle == "gamma":
            src_term = create_source_spectrum(gamma_gs  , source_def["energy"])
            for i in range (len(src_term)):
                index = i + len(neutron_gs)
                v_src[index] = src_term[i]
        else:
            raise Exception('particle must be neutron or gamma')
    print(v_src)        

    v_psi = np.matmul(A_inv,v_src)
    # print("Norm spectrum: ")
    # print(np.linalg.norm(v_psi))
  
    #======================================= Build data/energy
    outp = GenerateSpectrumData(neutron_gs, v_psi, sig_heat, gamma_gs)
  
    return outp

# ===================================================
def create_source_spectrum(group_structure,myE,fission=False):
    #Convert to eV. myE must be entered in MeV
    myE *= pow(10,6)

    Gn = len(group_structure)
    v_src = np.zeros(Gn) # this is a multigroup src spectrum (sum_g = 1)
    ###
    # find the idx in the neutron-gs when myE is located
    index = []
    for i in range (len(group_structure)):
        #If between bounds
        if ((myE < group_structure[i][2]) and (myE > group_structure[i][1])):
            index = [i]
        #If equals to upper bound
        elif ( (myE == group_structure[i][2]) and (i != (len(group_structure) - 1)) ):
            index = [i, i+1]
        #If equals to lower bound
        elif ( (myE == group_structure[i][1]) and (i != 0 ) ):
            index = [i, i-1]

    # if by chance, myE = existing bound, spread 0.5 to both bins
    ###
    if not fission:
        for idx in index:
            if len(index) == 1:
                v_src[idx] = 1.0 # from 13.94-14.2 MeV
            elif len(index) > 1:
                v_src[idx] = 0.5
            else:
                print("Error: The source energy is outside of the bounds")
            # v_src[0] = 1.0 # highest energy group
        v_src = np.flip(v_src)
    else:
        import scipy.integrate as integrate
        def chi(E):
            a = 0.965 # MeV
            b = 2.29  # 1/MeV
            return np.exp(-E/a) * np.sinh(np.sqrt(b*E))
        for i in range(Gn):
            Elow = group_structure[i][1]/1e6
            Eupp = group_structure[i][2]/1e6
            result = integrate.quad(lambda x: chi(x), Elow, Eupp)
            v_src[i] = result[0]
        v_src = np.flip(v_src)

    v_src /= np.sum(v_src)
    return v_src

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

