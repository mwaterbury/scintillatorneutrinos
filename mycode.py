import sys
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.image as img
import random
import time

from scipy.stats import poisson
from scipy.stats import norm
from scipy.integrate import simps

from os import listdir
import re


class Analyzer():
    
    ########################################################
    # Initialization and Prepare Weights
    ########################################################
    
    def __init__(self, prepare_muon=True, prepare_numu=True, prepare_nuel=True):
        
        # Logging
        print ('Prepare Data Manually')
        
        if prepare_numu: self.prepare_numu()
        else: self.numu_data = []
        if prepare_nuel: self.prepare_nuel()
        else: self.nuel_data = []
        if prepare_muon: self.prepare_muon()
        else: self.muon_data = []
    
    def weight_muon_neutrinos(self,data):
    
        # Import data for the number of expected muon neutrino interactions during run 3
        expectedData1 = np.loadtxt('FluenceFiles/FASER_-14.txt')
        expectedData2 = np.loadtxt('FluenceFiles/FASER_14.txt')
        
        expEnergy = expectedData1[:,0]
        expTotal  = (expectedData1[:,1] + expectedData2[:,1]) / 150 # To normalize to 1 fb^-1
        
        # Setup variables for weighting neutrino events.
        volFlu  = 25 * 25 * 100  # Volume of tungsten target used for FASER_-14.txt file
        rhoTung = 19.3           # Density of tungsten
        volTung = 25 * 30 * 106  # Volume of tungsten target used in FLUKA simulation
        rhoLead = 11.35          # Density of lead
        volLead = 30 * 40 * 10   # Volume of lead shielding used in FLUKA simulation
        
        # Weight by the masses of the targets
        tungWeight = (volTung/volFlu)
        leadWeight = volLead/volFlu * rhoLead/rhoTung
        
        # Weight per primary
        numSimTung = 13 * 1000
        numSimLead = 5 * 1000
        
        # Simulation energies are close to those in FASER_-14.txt
        # Use a small buffer to ensure there is a match
        eps  = 0.01 # Within 10 MeV
        for event in data:
            energy = event['primaryEnergy']
            zpos   = event['zposition']
            
            mask = (expEnergy < (energy + eps)) & (expEnergy > (energy - eps))
            numExp = expTotal[mask]
            if zpos < 120: # If interacted in tungsten
                weight = tungWeight * numExp / numSimTung
            else:
                weight = leadWeight * numExp / numSimLead
            
            event['weight'] = weight[0] # Method with the mask returns a numpy array with one-element

    def weight_elec_neutrinos(self, data):
        
        # Import data for the number of expected muon neutrino interactions during run 3
        expectedData1 = np.loadtxt('FluenceFiles/FASER_-12.txt')
        expectedData2 = np.loadtxt('FluenceFiles/FASER_12.txt')
        
        expEnergy = expectedData1[:,0]
        expTotal  = (expectedData1[:,1] + expectedData2[:,1]) / 150 # To normalize to 1 fb^-1
        
        # Setup variables for weighting neutrino events.
        volFlu  = 25 * 25 * 100  # Volume of tungsten target used for FASER_-14.txt file
        rhoTung = 19.3           # Density of tungsten
        volTung = 25 * 30 * 106  # Volume of tungsten target used in FLUKA simulation
        rhoLead = 11.35          # Density of lead
        volLead = 30 * 40 * 10   # Volume of lead shielding used in FLUKA simulation
        
        # Weight by the masses of the targets
        tungWeight = (volTung/volFlu)
        leadWeight = volLead/volFlu * rhoLead/rhoTung
        
        # Weight per primary
        numSimTung = 13 * 300
        numSimLead = 5 * 300
        
        # Simulation energies are close to those in FASER_-14.txt
        # Use a small buffer to ensure there is a match
        eps  = 0.01 # Within 10 MeV
        for event in data:
            energy = event['primaryEnergy']
            zpos   = event['zposition']
            
            mask = (expEnergy < (energy + eps)) & (expEnergy > (energy - eps))
            numExp = expTotal[mask]
            if zpos < 120: # If interacted in tungsten
                weight = tungWeight * numExp / numSimTung
            else:
                weight = leadWeight * numExp / numSimLead
            
            event['weight'] = weight[0] # Method with the mask returns a numpy array with one-element

    def weight_muons(self, files, veto):
        
        data = np.array([])
        
        # Load fluence files for muons and anti-muons and combine them
        expectedData1 = np.loadtxt('FluenceFiles/negative_muon_flux.csv')
        expectedData2 = np.loadtxt('FluenceFiles/positive_muon_flux.csv')
        
        expEnergy = expectedData1[:,0]
        expFlu    = expectedData1[:,1] + expectedData2[:,1]
        
        # expFlu has units of GeV^-1 cm^-2 s^-1
        # Multiply by energy bin width (200 GeV/bin) to change to cm^-2 s^-1 /bin
        expFlu = expFlu * 200
        
        LHCflu = 2 * (10 ** -5) # fb^-1 s^-1 from LHC
        expFlu = expFlu / LHCflu # fb (at the LHC) * cm^-2 (per unit area at FASER) per bin
        
        A_lr = 6.5 * 35 # left/right simulation area in cm^2
        A_tb = 43  * 4 # top/bot simulation area in cm^2
        A_c = 35 * 30 # center cimulation area in cm^2

        expTot_lr = expFlu * A_lr
        expTot_tb = expFlu * A_tb
        expTot_c = expFlu * A_c * veto
        
        # Number simulated in each region for each input file
        simTot = self.countSimTotals(files)
        
        # Controls accepted range for matching primary energy - Anything below ~200 GeV should be identical
        eps  = 0.01 # Within 10 MeV
        
        for file in files:
            tmp = np.load('NumpyArrays/' + file, allow_pickle=True)
            side = re.match(r'.*_([a-z]*)[0-9][0-9][0-9].npy', file)[1]
            if (side == 'left') or (side == 'right'):
                expTot = expTot_lr
            elif (side == 'top') or (side == 'bot'):
                expTot = expTot_tb
            elif (side == 'center'):
                expTot = expTot_c
            
            # Add the appropriate weight to each event
            for event in tmp:
                energy = event['primaryEnergy']
                simTotCur = simTot[side][energy]
                
                # Find which energy this primary falls into
                mask = (expEnergy < (energy + eps)) & (expEnergy > (energy - eps))
                # Number of expected muons at this energy
                numExp = expTot[mask]
                # weight = number expected / number simulated
                event['weight'] = numExp[0] / simTotCur # Weight for 1 fb^-1 run at the LHC
            
            data = np.append(data,tmp)
        
        return data

    def countSimTotals(self,files):
        
        simTot = {}
        for side in ['left', 'right', 'top', 'bot']:
            sideDict = {}
            sideFiles = [file for file in files if re.search(f'.*{side}.*', file)]
            for file in sideFiles:
                energy = np.load('NumpyArrays/' + file, allow_pickle=True)[0]['primaryEnergy']
                if not(energy in sideDict.keys()):
                    sideDict[energy] = 21527
                else:
                    sideDict[energy] += 21527
            simTot[side] = sideDict
        # Determine simTot for center
        centerDict = {}
        centerFile = [file for file in files if re.search(f'.*center.*', file)]
        for file in centerFile:
            energy = np.load('NumpyArrays/' + file, allow_pickle=True)[0]['primaryEnergy']
            if not(energy in centerDict.keys()):
                centerDict[energy] = 20000
            else:
                centerDict[energy] += 20000
        simTot['center'] = centerDict
        
        return simTot

    def prepare_numu(self,):
        
        # Prepare muon neutrinos
        start_time = time.time()
        ls = listdir('NumpyArrays/')

        files = [file for file in ls if re.search(r'muon_neutrinos',file)]
        self.numu_data = np.array([])
        for file in files:
            tmp = np.load('NumpyArrays/' + file, allow_pickle=True)
            self.weight_muon_neutrinos(tmp)
            np.save('NumpyArrays/' + file, tmp)
            
            self.numu_data = np.append(self.numu_data, tmp)

        print (' ... found', len(self.numu_data), 'numu events in', round(time.time() - start_time,2), 'seconds')

    def prepare_nuel(self,):
        # Prepare electron neutrinos
        start_time = time.time()
        ls = listdir('NumpyArrays/')

        files = [file for file in ls if re.search(r'elec_neutrinos',file)]
        self.nuel_data = np.array([])
        for file in files:
            tmp = np.load('NumpyArrays/' + file, allow_pickle=True)
            self.weight_elec_neutrinos(tmp)
            np.save('NumpyArrays/' + file, tmp)
            
            self.nuel_data = np.append(self.nuel_data, tmp)

        print (' ... found', len(self.nuel_data), 'nue events in', round(time.time() - start_time,2), 'seconds')

    def prepare_muon(self, veto=10**-6):
        # Prepare muon
        start_time = time.time()
        ls = listdir('NumpyArrays/')

        files = [file for file in ls if re.search(r'muons',file)]
        self.muon_data = self.weight_muons(files, veto)
        
        print (' ... found', len(self.muon_data), 'muon events in', round(time.time() - start_time,2), 'seconds')

    def get_muon_statistics(self,):

        # get weights
        objects = []
        for event in self.muon_data:
            tmp = {'energy': event['primaryEnergy'],
                'weight': event['weight']}
            if not(tmp in objects):
                objects.append(tmp)
        sorted(objects, key=lambda d: d['energy'])
        
        # peroperly weight lr and tb
        A_lr = 6.5 * 35 # left/right simulation area in cm^2
        A_tb = 43  * 4 # top/bot simulation area in cm^2

        weights = {}
        for obj in objects:
            energy = obj['energy']
            weight = obj['weight']
            if not(energy in weights.keys()):
                weights[energy] = weight
            else:
                if weight > weights[energy]:
                    weight_lr = weight
                    weight_tb = weights[energy]
                else:
                    weight_lr = weights[energy]
                    weight_tb = weight
                weights[energy] = (weight_lr * A_lr + weight_tb * A_tb)/(A_lr + A_tb)

        # print
        weights = dict(sorted(weights.items(), key=lambda item: item[0]))
        for item in weights:
            print (item, ":", round(weights[item],3))

    ########################################################
    # Event Displays
    ########################################################

    def display_event(self, event, filename=None):
        
        # initialize
        multiplier=2
        matplotlib.rcParams.update({'font.size': 15*multiplier})
        matplotlib.rcParams['axes.linewidth'] = multiplier
        fig = plt.figure(figsize=(16*multiplier,7*multiplier))
        
        # Setup some nice axes
        x0, x1 = 0.01, 0.94
        width = (x1-x0)/4
        h_tr, h_pic, h_sc = 0.55, 0.35, 0.08
        ax0 = fig.add_axes([x0          , 0.02+h_tr      , (x1-x0)   , h_pic])
        ax1 = fig.add_axes([x0 + 0*width, 0.02           , width     , h_tr ])
        ax2 = fig.add_axes([x0 + 1*width, 0.02           , width     , h_tr ])
        ax3 = fig.add_axes([x0 + 2*width, 0.02           , width     , h_tr ])
        ax4 = fig.add_axes([x0 + 3*width, 0.02           , 1.25*width, h_tr ])
        ax5 = fig.add_axes([x0          , 0.02+h_tr+h_pic, 0.8       , 0.05  ])
        ax6 = fig.add_axes([x0+0.8      , 0.02+h_tr+h_pic, x1-x0-0.8 , 0.05  ])
        
        #plot image of detector
        image = img.imread('Fig_Layout.jpg')
        ax0.imshow(image)
        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.spines.clear()
        
        # Plot the tracker images
        bins = (np.linspace(-12.5, 12.5, 26),np.linspace(-12.5, 12.5, 26))
        for ax, layer in zip([ax1,ax2,ax3,ax4],['hits1','hits2','hits3','hits4']):
            image = np.array(event[layer])
            if len(image)==0: image = np.array([[0,0,0]])
            hist = ax.hist2d(image.T[0], image.T[1], weights=image.T[2], bins=bins,
                             norm=matplotlib.colors.LogNorm(vmin=1,vmax=100),cmap='cool')
            ax.set_xticks([])
            ax.set_yticks([])
        fig.colorbar(hist[3], ax=ax4)
        
        # Plot scintillators
        ax5.hist2d(np.linspace(0.5,7.5,8), np.zeros(8) ,
                   weights=event['scintillator'][1:9],
                   bins=(np.linspace(0,8,9),[-1,1]),
                   norm=matplotlib.colors.LogNorm(vmin=0.1,vmax=1),
                   cmap='rainbow')
        for x in np.linspace(0,7,8):
            ax5.plot([x,x],[-1,1],lw=multiplier,c="k")
            ax5.text(x+0.06,-0.75,str(int(x+1)))
        ax5.set_xticks([])
        ax5.set_yticks([])
        
        # Plot calorimeter
        ax6.text(0.5,0.5,"E="+str(event['calorimeter'])+" GeV",ha="center", va="center")
        ax6.fill([0,1,1,0], [0,0,1,1], c="lightgreen")
        ax6.set_xlim(0,1)
        ax6.set_ylim(0,1)
        ax6.set_xticks([])
        ax6.set_yticks([])
        
        # Return figure so user can manipulate or save the plots
        if filename is not None:
            plt.savefig(filename)
        plt.show()


    def display_random_event(self, particle="muon", requirement="True"):
        
        #check particle type:
        if particle not in ['muon', 'numu', 'nuel']:
            print ("particle must be either 'muon', 'numu' or 'nuel'")
            return 1
        
        found=False
        while found is not True:
            if particle == 'muon':
                length = len(self.muon_data)
                ievent = random.randrange(length)
                filename = "EventDisplays/Event_muon_"+str(ievent)+".pdf"
                event = self.muon_data[ievent]
                if eval(requirement): found=True
            elif particle == 'numu':
                length = len(self.numu_data)
                ievent = random.randrange(length)
                filename = "EventDisplays/Event_numu_"+str(ievent)+".pdf"
                event = self.numu_data[ievent]
                if eval(requirement): found=True
            elif particle == 'nuel':
                length = len(self.nuel_data)
                ievent = random.randrange(length)
                filename = "EventDisplays/Event_nuel_"+str(ievent)+".pdf"
                event = self.nuel_data[ievent]
                if eval(requirement): found=True
        self.display_event(event,filename=filename)

    
    ########################################################
    # Observables
    ########################################################

    def define_observable_from_function(self, name, function):
        for event in self.muon_data: event[name] = function(event)
        for event in self.numu_data: event[name] = function(event)
        for event in self.nuel_data: event[name] = function(event)

    def define_observable(self, name, definition):
        self.define_observable_from_function(name,
                eval('lambda event : ' + definition))
#         for event in self.muon_data: event[name] = eval(definition)
#         for event in self.numu_data: event[name] = eval(definition)
#         for event in self.nuel_data: event[name] = eval(definition)
    
    def get_histodata(self, particle, observable, bins=None, requirement="True"):
    
        #check particle type, assign data:
        particles = {'muon': self.muon_data, 'numu': self.numu_data, 'nuel': self.nuel_data}
        if particle in particles.keys(): data = particles[particle]
        else: print ("Error: particle must be either 'muon', 'numu' or 'nuel'")
        
        # loop through events
        
        # Writing requirement as a lambda function runs much faster (~5x)
        # - I presume this has to do with interpretting the string as code
        #   before running the loop rather than looping and then intrepetting
        #   the string as code each iteration
        func = eval('lambda event : ' + requirement)

#         selected_data = []
#         for event in data:
#             if func(event) == False: continue
#             value = event[observable]
#             weight = event['weight']
#             selected_data.append([value, weight, weight*weight])
#         selected_data = np.array(selected_data)

        # Python's List comprehension is faster than looping when appending
        selected_data = [[event[observable],
                          event['weight'],
                          event['weight']**2]
                             for event in data if func(event)]
        selected_data = np.array(selected_data)

        # bin data get sum of weights, and squared sum of weights for each bin
        yvals , _ = np.histogram(selected_data.T[0], weights=selected_data.T[1], bins=bins)
        yvals2, _ = np.histogram(selected_data.T[0], weights=selected_data.T[2], bins=bins)

        # reformat binned data
        # for uncertainties, see https://www.zeuthen.desy.de/~wischnew/amanda/discussion/wgterror/working.html
        xvals = (bins[:-1] + bins[1:]) / 2.
        binned_data = np.array([xvals, yvals, np.sqrt(yvals2)])

        return binned_data
    

    def plot_histogram(self, dataset, logx=False, logy=False, xlim=None, ylim=None, xlabel="Observable", filename=None):
        
        # prepare plot
        matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=(8,6))
        
        # plot binned data
        ax = plt.subplot(1,1,1)
        for data, color, label in dataset:
            for x,y,yerr in data.T: ax.errorbar(x, y, yerr=yerr, fmt='.', color=color)
            ax.scatter(data.T[0][0], data.T[0][1], color=color, label=label, marker='.')
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Expected Event Rate [fb/bin]")
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
        if xlim is not None: ax.set_xlim(xlim[0], xlim[1])
        if ylim is not None: ax.set_ylim(ylim[0], ylim[1])

        # finalize
        ax.legend(frameon=False, labelspacing=0, fontsize=13)
        ax.grid()
        if filename is not None: fig.savefig(filename)
        plt.show()


    def plot_2dhistogram(self,
        observable1, observable2,
        bins1, bins2,
        requirement=None,
        xlabel="Observable1", ylabel="Observable2",
        logx=False, logy=False,
        filename=None
    ):

        # check if data exists:
        if len(self.numu_data)==0: print ("Error: no numu signal data")
        if len(self.muon_data)==0: print ("Error: no muon background data")

        # prepare plot
        matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=(8*3,6))
        
        # loop through events
        selected_data_s, selected_data_b = [], []
        for event in self.numu_data:
            if eval(requirement) == False: continue
            selected_data_s.append([event[observable1], event[observable2], event['weight']])
        for event in self.muon_data:
            if eval(requirement) == False: continue
            selected_data_b.append([event[observable1], event[observable2], event['weight']])
        selected_data_s = np.array(selected_data_s)
        selected_data_b = np.array(selected_data_b)
        
        # plot signal data
        ax = plt.subplot(1,3,1)
        hs = ax.hist2d(
            selected_data_s.T[0], selected_data_s.T[1], weights=selected_data_s.T[2],
            bins=[bins1,bins2], range=[[bins1[0], bins1[-1]],[bins2[0], bins2[-1]]],
            norm=matplotlib.colors.LogNorm(),cmap='jet',
        )
        plt.colorbar(hs[3])
        ax.set_title("Muon Neutrino [fb/bin")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
                   
        ax = plt.subplot(1,3,2)
        hb = ax.hist2d(
            selected_data_b.T[0], selected_data_b.T[1], weights=selected_data_b.T[2],
            bins=[bins1,bins2], range=[[bins1[0], bins1[-1]],[bins2[0], bins2[-1]]],
            norm=matplotlib.colors.LogNorm(),cmap='jet',
        )
        plt.colorbar(hb[3])
        ax.set_title("Muons [fb/bin]")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
                   
        selected_data_r = []
        vmin, vmax = 10e10,0
        for ix in range(len(bins1)-1):
            x = (bins1[ix+1]+bins1[ix])/2.
            for iy in range(len(bins2)-1):
                y = (bins2[iy+1]+bins2[iy])/2.
                if hb[0][ix][iy]==0: r = 0
                else: r=hs[0][ix][iy]/hb[0][ix][iy]
                selected_data_r.append([x,y,r])
                if r>0 and r<vmin: vmin=r
                if r>vmax: vmax=r
        selected_data_r = np.array(selected_data_r)

        ax = plt.subplot(1,3,3)
        hr = ax.hist2d(
            selected_data_r.T[0], selected_data_r.T[1], weights=selected_data_r.T[2],
            bins=[bins1,bins2], range=[[bins1[0], bins1[-1]],[bins2[0], bins2[-1]]],
            norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax),cmap='jet',
        )
        plt.colorbar(hr[3])
        ax.set_title("Ratio: Muon Neutrino / Muon")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if logx: ax.set_xscale("log")
        if logy: ax.set_yscale("log")
    
        if filename is not None: fig.savefig(filename)
        plt.show()

    def plot_2dhistogram_primary(self, observable, bins, requirement=None,
        label="Observable", log=False, filename=None ):
    
        # check if data exists:
        if len(analyser.numu_data)==0: print ("Error: no numu signal data")
        if len(analyser.muon_data)==0: print ("Error: no muon background data")
        
        # prepare plot
        matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=(8*2,6))
        bins_b = np.linspace(0,3600,18+1)
        bins_s = np.logspace(1,4,30+1)
        
        # loop through events
        selected_data_s, selected_data_b = [], []
        for event in self.numu_data:
            if eval(requirement) == False: continue
            selected_data_s.append([event['primaryEnergy'], event[observable], event['weight']])
        for event in self.muon_data:
            if eval(requirement) == False: continue
            selected_data_b.append([event['primaryEnergy'], event[observable], event['weight']])
        selected_data_s = np.array(selected_data_s)
        selected_data_b = np.array(selected_data_b)

        # plot signal data
        ax = plt.subplot(1,2,1)
        hs = ax.hist2d(
            selected_data_s.T[0], selected_data_s.T[1], weights=selected_data_s.T[2],
            bins=[bins_s,bins], range=[[bins_s[0], bins_s[-1]],[bins[0], bins[-1]]],
            norm=matplotlib.colors.LogNorm(),cmap='jet',
        )
        plt.colorbar(hs[3])
        ax.set_title("Muon Neutrino [fb/bin")
        ax.set_xlabel('Neutrino Energy [GeV]')
        ax.set_ylabel(label)
        ax.set_xscale("log")
        if log: ax.set_yscale("log")
                   
        ax = plt.subplot(1,2,2)
        hb = ax.hist2d(
            selected_data_b.T[0], selected_data_b.T[1], weights=selected_data_b.T[2],
            bins=[bins_b,bins], range=[[bins_b[0], bins_b[-1]],[bins[0], bins[-1]]],
            norm=matplotlib.colors.LogNorm(),cmap='jet',
        )
        plt.colorbar(hb[3])
        ax.set_title("Muons [fb/bin]")
        ax.set_xlabel('Muon Energy [GeV]')
        ax.set_ylabel(label)
        if log: ax.set_yscale("log")
                   
        if filename is not None: fig.savefig(filename)
        plt.show()
    
    ########################################################
    # Count events after cuts
    ########################################################
    
    def count_events(self, particle, requirement="True"):
        
        # check particle type:
        if particle not in ['muon', 'numu', 'nuel']:
            print ("particle must be either 'muon', 'numu' or 'nuel'")
            return 0
        
        # asign data
        if particle == 'muon':   data = self.muon_data
        elif particle == 'numu': data = self.numu_data
        elif particle == 'nuel': data = self.nuel_data
        
        # evaulate weights and check requirements
        func = eval('lambda event: ' + requirement)
        weights = [event['weight'] for event in data if func(event)]
        weights = np.array(weights)
        
        # obtain count and error
        count = np.sum(weights)
        error = np.sqrt(np.sum(weights ** 2))
        
        return count, error
    

    def calculate_stats(self, requirement="True", luminosity = 1):
        
        def integrate_with_gaussian(f, mu, sigma):
            if sigma <= 0: return f
            xvals = np.linspace(mu - 10 * sigma, mu + 10 * sigma)
            yvals = f * norm.pdf((xvals-mu)/sigma) / sigma
            return simps(yvals, xvals)

        # get signal and background rate
        muon_mean, muon_err = self.count_events(particle='muon', requirement=requirement)
        numu_mean, numu_err = self.count_events(particle='numu', requirement=requirement)
        nuel_mean, nuel_err = self.count_events(particle='nuel', requirement=requirement)

        bkg_mean, bkg_err = muon_mean, muon_err
        sig_mean, sig_err = numu_mean + nuel_mean, np.sqrt(numu_err**2+nuel_err**2)
        
        # multiply luminosity
        sig_mean, sig_err = sig_mean*luminosity, sig_err*luminosity
        bkg_mean, bkg_err = bkg_mean*luminosity, bkg_err*luminosity

        x_max = int(np.round(sig_mean + sig_err + bkg_mean + bkg_err))

        p_func  = lambda x,y : integrate_with_gaussian(poisson.pmf(x,y), bkg_mean, bkg_err)

        xvals   = list(range(5*x_max))
        p_xvals = [p_func(n,bkg_mean) for n in xvals]
        pvals   = [sum(p_xvals[n:]) for n in xvals]

        BF_func = lambda x, y : p_func(x, bkg_mean + y)
        BF = [pvals[n] / integrate_with_gaussian(BF_func(n, sig_mean), sig_mean, sig_err) for n in xvals]

        # summarize results
        analysis = {}
        analysis['nsignal']     = (sig_mean*luminosity, sig_err*luminosity)
        analysis['nbackground'] = (bkg_mean*luminosity, bkg_err*luminosity)
        
        analysis['x-values'] = xvals
        analysis['p-values'] = pvals
        analysis['Bayes Factors'] = BF

        analysis['p-expected'] = sum([pvals[n] * poisson.pmf(n, sig_mean + bkg_mean) for n in xvals])
        analysis['BF-expected'] = sum([BF[n] * poisson.pmf(n, sig_mean + bkg_mean) for n in xvals])

        return analysis

class Poisson_Analysis():
    
    def marg_like(self, No, Ns, Nb, db):
        bs = np.linspace(max(Nb-10*db,0), Nb + 10*db, 200)
        integrand = [poisson.pmf(No, Ns + b) * norm.pdf(b, Nb, db)
                     for b in bs]
        return simps(integrand, bs, axis=0)
    
    def __init__(self, Ns, Nb, db, alpha = 10**-8):
        lo = int(np.floor(poisson.ppf(alpha*10**-6, Nb)))
        hi = int(np.ceil(poisson.ppf(1-alpha*10**-6, Ns + Nb)))
        No = np.arange(lo-1,hi+1)
        
        mus = np.linspace(0,3,301)
        null_mask = mus == 0
        alt_mask = mus == 1
        
        likes = [poisson.pmf(No, Nb + mu * Ns) for mu in mus]
        pvals = [poisson.sf(No, Nb + mu * Ns) for mu in mus]
        sig = [np.abs(norm.ppf(pval)) for pval in pvals]

        marg_likes = [self.marg_like(No, mu*Ns, Nb, db) for mu in mus]
        marg_pvals = [np.array([1-np.sum(like[:n]) for n in range(len(No))[1:]]) for like in marg_likes]
        marg_sig = [np.abs(norm.ppf(pval)) for pval in marg_pvals]
        
        # Make everything a numpy array for good measure
        likes = np.array(likes)
        pvals = np.array(pvals)
        sig = np.array(sig)
        marg_likes = np.array(marg_likes)
        marg_pvals = np.array(marg_pvals)
        marg_sig = np.array(marg_sig)
        
        self.data = {'Ns' : Ns, 'Nb' : Nb, 'db' : db}
        self.accuracy = alpha
        self.mus = mus
        self.No = No
        self.likelihood = likes
        self.pvals = pvals
        self.significance = sig
        self.marginal_likelihood = marg_likes
        self.marginal_pvals = marg_pvals
        self.marginal_significance = marg_sig
        self.bayes_factor = likes[null_mask] / likes[alt_mask]
        self.marginal_bayes_factor = marg_likes[null_mask] / marg_likes[alt_mask]
        
        pvals = pvals[null_mask]
        sig   = sig[null_mask]
        like1 = likes[alt_mask]
        exp_pval = np.sum([p*like for (p,like) in zip(pvals, like1)])
        exp_sig  = np.sum([s*like for (s,like) in zip(sig, like1)])
        
        marg_pvals = marg_pvals[null_mask]
        marg_sig   = marg_sig[null_mask]
        marg_like1 = marg_likes[alt_mask]
        exp_marg_pval = np.sum([p*like[1:]
                           for (p,like) in zip(marg_pvals, marg_like1)])
        exp_marg_sig  = np.sum([s*like[1:]
                           for (s,like) in zip(marg_sig, marg_like1)])
        self.expected = {
            'pval' : exp_pval,
            'significance' : exp_sig,
            'marginal_pval' : exp_marg_pval,
            'marginal_significance' : exp_marg_sig
        }

    def plot_marginal_comparison(self):
       
        # Pull data from self
        Ns    = self.data['Ns']
        Nb    = self.data['Nb']
        db    = self.data['db']
        No = self.No
        lo = No[0]+1
        hi = No[len(No)-1]-1
        like = self.likelihood
        pvals_with_mus = self.pvals
        marg_like = self.marginal_likelihood
        marg_pval_with_mus = self.marginal_pvals
        
        # Setup variables for plotting 
        alpha = self.accuracy
        yticks = 10**np.arange(np.log10(alpha),0.1,2)
        Zs = np.arange(1,8.1)
        pval_sigmas = np.array([1-norm.cdf(Z) for Z in Zs])
        Zs = Zs[pval_sigmas >= alpha]
        pval_sigmas = pval_sigmas[pval_sigmas >= alpha]
        zeros = np.zeros(np.shape(No))
        
        # Masks for extracting null and alternative hypotheses
        mus = self.mus
        null_mask = mus == 0
        alt_mask  = mus == 1

        #Setup variables for cleaner code
        like0 = like[null_mask][0]
        like1 = like[alt_mask][0]
        pvals  = pvals_with_mus[null_mask][0]
        
        marg_like0 = marg_like[null_mask][0]
        marg_like1 = marg_like[alt_mask][0]
        marg_pvals = marg_pval_with_mus[null_mask][0]

        matplotlib.rcParams.update({'font.size': 14})
        fig = plt.figure(figsize=(8*2,6))

        ax = plt.subplot(1,2,1)
        
        # Mark mean value of No
        ax.plot([Ns + Nb, Ns + Nb], [alpha*10**-0.5,10**0.5], 'r--')
        ax.text( Ns + Nb + 0.01*(hi-lo), 10**0.125, '$N_s + N_b$', color='red')
        
        # Plot null hypothesis likelihoods
        ax.scatter(No, like0,
                   color='blue', marker='<', s=16, label='$f(N_o, N_b)$')
        ax.scatter(No, marg_like0,
                   color='orange', marker='<', s=16, label='$f_\mathrm{marg}(N_o, N_b, \sigma_b)$')

        # Plot alterative hypothesis likelihoods
        ax.scatter(No, like1, color='blue', marker='>', s=16, label='$f(N_o, N_b + N_s)$')
        ax.scatter(No, marg_like1, color='orange', marker='>', s=16, label='$f_\mathrm{marg}(N_o, N_s, N_b, \sigma_b)$')
        for (Z,pval_sigma) in zip(Zs,pval_sigmas):
            if pval_sigma >= alpha:
                mask = pvals <= pval_sigma
                ax.fill_between(No[mask], zeros[mask], like0[mask], color = 'blue', alpha=0.05)
                ax.text(No[mask][0],alpha*10**-0.375,f'{int(Z)}$\sigma$', fontsize='x-small')

                mask = np.append(False, marg_pvals <= pval_sigma)
                ax.fill_between(No[mask], zeros[mask], marg_like0[mask], color = 'orange', alpha=0.05)

        ax.legend(loc="upper right")

        ax.set_xlabel('$N_o$')
        ax.set_ylabel('Likelihood')
        ax.set_xlim([lo,hi])
        ax.set_ylim([alpha*10**-0.5,10**0.5])
        ax.set_yscale('log')
        ax.set_yticks(yticks);

        ax = plt.subplot(1,2,2)
        
        # Mark mean value of No
        ax.plot([Ns + Nb, Ns + Nb], [alpha*10**-0.5,10**0.5], 'r--')
        ax.text(Ns+Nb+0.01*(hi-lo),10**0.125,'$N_s + N_b$', color='red')
        
        # Plot p-values
        ax.scatter(No, pvals, s=8, label='$\hat{p}(N_o)$')
        ax.scatter(No[1:], marg_pvals, s=8, label='$\hat{p}_\mathrm{marg}(N_o)$')
        
        # Mark p-values corresponding to significance levels
        for (Z,pval_sigma) in zip(Zs,pval_sigmas):
            ax.plot([lo-1, hi+1], [pval_sigma, pval_sigma], 'k--')
            ax.text(hi-0.06*(hi-lo), pval_sigma*10**0.125, f'{int(Z)}$\sigma$')
        ax.legend(loc="lower left")

        ax.set_xlabel('$N_o$')
        ax.set_ylabel('$p$-value')
        ax.set_xlim([lo,hi])
        ax.set_ylim([alpha*10**-0.5,10**0.5])
        ax.set_yscale('log')
        ax.set_yticks(yticks)