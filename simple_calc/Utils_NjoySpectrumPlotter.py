import numpy as np 
import matplotlib.pyplot as plt

#======================================= Plot the spectra
def Njoy_spectrum_plotter(outp, path):
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
    #================================= Plot energy spectrum
    #================================= Flux
    fig = plt.figure(figsize=(12,6), )
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
    plt.savefig(path+'Neutron_spectrum.png')

    #================================= Heating XS
    # fig = plt.figure(figsize=(12,6))
    # fig.subplots_adjust(hspace=0.4, wspace=0.4)
    # ax = fig.add_subplot(1, 2, 1)
    # ax.semilogy(neutron_group_bndries, neutron_heating_spectrum,color='r')
    # ax.set_xlabel("Energy (MeV)")
    # ax.set_ylabel("H(E) (eV/s)")
    # ax.grid('on')
    # ax = fig.add_subplot(1, 2, 2)
    # ax.loglog(neutron_group_bndries, neutron_heating_spectrum, color='r')
    # ax.set_xlabel("Energy (MeV)")
    # ax.set_ylabel("H(E) (eV/s)")
    # ax.grid('on')
    # plt.suptitle('neutron heating')
    # plt.savefig(path+'Neutron_Heating_spectrum.png')

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
      plt.savefig(path+'Gamma_spectrum.png')

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
      plt.savefig(path+'Gamma_heating_spectrum.png')

