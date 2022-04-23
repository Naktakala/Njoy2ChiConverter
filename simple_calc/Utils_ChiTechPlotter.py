"""Combines NJOY raw data into comprehensible transport data"""
import numpy as np 
import matplotlib.pyplot as plt
import dtra_n_g_mcnp_post_process as mcnp_reader
from matplotlib.ticker import MultipleLocator
import matplotlib.dates as mdates

import Utils_Info as processor


def PlotSpectra(big_processed_data, source, path, mcnp_filename = ""):
  #========================== Check for type of plots ==================#
  #from dtra_n_g_mcnp_post_process import n_spectrum
  if (mcnp_filename != ""):
    tally_lists = ['24','26','214','216']
    tally, tally_dict, atom_density, gram_density = mcnp_reader.ReadMcnpFile(mcnp_filename, tally_lists)
    print("Gram density = %.4f g/cc"%(gram_density))
    for particle_type in tally_dict.keys():
      for reaction_type in tally_dict[particle_type]:
        if len(tally_dict[particle_type][reaction_type]) > 0:
          fig, ax = plt.subplots(dpi = 750)
          ax.set_xlabel('E (MeV)')
          ax.xaxis.set_label_coords(0.5, -0.07)
          if reaction_type == "Flux":
            ax.set_ylabel('$\phi$(E) ($n$ $cm^{-2}$ $MeV^{-1}$)')
          elif reaction_type == "Heating":
            ax.set_ylabel('H(E) (eV $s^{-1}$)')
          ax.yaxis.set_label_coords(-0.09, 0.5)
          plt.suptitle("MCNP vs NJOY: " + particle_type + ' ' + reaction_type + " Spectrum", y=1)
          ax.spines['right'].set_visible(False)
          ax.spines['top'].set_visible(False)
          plt.setp(ax, xticks=[], yticks=[])
          fig.tight_layout(pad=4.0)
          row, col = 1,2
          lw=1 # Linewidth
          for r in range (len(tally_dict[particle_type][reaction_type])):
            spectrum = []
            Sn_E = []
            Sn_Flx = []
            #Get MCNP data
            tally_name = tally_dict[particle_type][reaction_type][r]
            E, Flx = [], []
            for tally_value in tally:
              if tally_value['tallyID'] == tally_name:
                E = tally_value['values'][:,0]
                Flx = tally_value['values'][:,1]
                spectrum.append(Flx)
                spectrum.append(E)
                E, Flx = [], []
                E, Flx = mcnp_reader.convert_spectrum(spectrum)
                
                # Flip MCNP spectrum
                E = np.flip(E)
                Flx = np.flip(Flx)

                if reaction_type == "Flux":
                  #If flux
                  #Flx = [i*pow(10,12) for i in Flx]
                  Flx = [i*pow(10,12) for i in Flx]
                elif reaction_type == "Heating":
                    #If heating:
                  #Flx = [i*pow(10,18)*gram_density for i in Flx]
                  Flx = [i*pow(10,18)*gram_density for i in Flx]
                # # Normalization
                # max = np.max(Flx)
                # Flx /= max
            #G = len(list(set(E)))
            #Start plotting
            ax = fig.add_subplot(row,col,r+1)
            ax.plot(E, Flx, lw=lw, linewidth = 2, label = "$MCNP_{48}$")
            #Get Njoy data
            for outp in big_processed_data:
              # If deals with neutrons 
              neutron_group_bndries = outp[0][0]
              neutron_spectrum = outp[0][1]
              gamma_group_bndries = outp[1][0]
              gamma_spectrum = outp[1][1]
              neutron_heating_spectrum = outp[2][0]
              gamma_heating_spectrum = outp[2][1]

              #Get the group structure
              G_n = len(list(set(neutron_group_bndries)))
              G_g = len(list(set(gamma_group_bndries)))
              if G_n > 0: 
                G_n -= 1
              if G_g > 0:
                G_g -= 1

              if particle_type == "Neutron":
                # If deals with flux:
                if reaction_type == "Flux":
                  Sn_E = neutron_group_bndries
                  Sn_Flx = neutron_spectrum
                # If deals with heating:
                elif reaction_type == "Heating":
                  Sn_E = neutron_group_bndries
                  Sn_Flx = neutron_heating_spectrum
              elif particle_type == "Gamma":
                if (gamma_spectrum != []):
                  # If deals with flux:
                  if reaction_type == "Flux":
                    Sn_E = gamma_group_bndries
                    Sn_Flx = gamma_spectrum
                  # If deals with heating:
                  elif reaction_type == "Heating":
                    Sn_E = gamma_group_bndries
                    Sn_Flx = gamma_heating_spectrum
                
                # max = np.max(Sn_Flx)
                # Sn_Flx /= max
                    
              ax.plot(Sn_E, Sn_Flx, lw=lw,linewidth = 1.5, label = "$NJOY_{n%.dg%.d}$"%(G_n, G_g), color ='r')
            #ax.plot(E, Flx, lw=lw, label = "$MCNP_{%.d}$"%(G_n), color = 'y')

            #ax.set_xlim(right = 1.2e1)
            # if reaction_type == 'Heating':
            #   ax.set_ylim(bottom = 10e3)
            # else:
            #   ax.set_ylim(bottom = 10e1)

            if float(tally_name) > 210:
              ax.set_title("Fine MCNP GS", fontsize=11)
            else:
              ax.set_title("Coarse MCNP GS", fontsize=11)
            ax.set_yscale('log')
            ax.set_xscale('log')
            plt.xticks(fontsize=8)
            plt.yticks(fontsize=8)
            plt.grid()
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            ax.xaxis.set_major_locator(MultipleLocator(3))
        
            #ax.xaxis.set_minor_locator(MultipleLocator(1))
            #ax.yaxis.set_minor_locator(MultipleLocator(0.2))
            # if r == (len(tally_dict[particle_type][reaction_type]) - 1):
              #plt.legend(bbox_to_anchor=(0.5, 0.0, 2.0, 0.75), loc=8, borderaxespad=0., frameon=False, fontsize=10)
          handle_list, label_list = [], []
          for ax in fig.axes:
            handles, labels = ax.get_legend_handles_labels()
            for handle, label in zip(handles, labels):
                if label not in label_list:
                    handle_list.append(handle)
                    label_list.append(label)
          plt.legend(handle_list, label_list, loc = 0)
          fig.subplots_adjust(wspace=0.585, hspace=0.485, right = 0.994, left = 0.12, top = 0.90, bottom = 0.1)
          plt.savefig(path + '/' + particle_type + '_' + reaction_type + '_SpectraSource_Histogram'+ '.png')

  else:
    for outp in big_processed_data:
      # If deals with neutrons 
      neutron_group_bndries = outp[0][0]
      neutron_spectrum = outp[0][1]
      gamma_group_bndries = outp[1][0]
      gamma_spectrum = outp[1][1]
      neutron_heating_spectrum = outp[2][0]
      gamma_heating_spectrum = outp[2][1]

      #Get the group structure
      G_n = len(list(set(neutron_group_bndries)))
      G_g = len(list(set(gamma_group_bndries)))
      if G_n > 0: 
        G_n -= 1
      if G_g > 0:
        G_g -= 1

      if (neutron_spectrum != []):
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
        fig.legend()
        plt.suptitle('Neutron spectrum')
        plt.savefig(path+'Njoy Neutron_spectrum_n%.dg%.d.png'%(G_n, G_g))

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
        fig.legend()
        plt.suptitle('Neutron heating')
        plt.savefig(path+'Njoy Neutron_Heating_spectrum_n%.dg%.d.png'%(G_n, G_g))

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
        fig.legend()
        plt.suptitle('Gamma spectrum')
        plt.savefig(path+'Njoy Gamma_spectrum_n%.dg%.d.png'%(G_n, G_g))

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
        fig.legend()
        plt.suptitle('Gamma heating')
        plt.savefig(path+'Njoy Gamma_heating_spectrum_n%.dg%.d.png'%(G_n, G_g))


