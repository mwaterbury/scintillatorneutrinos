import sys
import os
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.image as img
import random
import time

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

    def weight_muons(self, files):
        
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
        
        expTot_lr = expFlu * A_lr
        expTot_tb = expFlu * A_tb
        
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

    def prepare_muon(self,):
        # Prepare muon
        start_time = time.time()
        ls = listdir('NumpyArrays/')

        files = [file for file in ls if re.search(r'muons',file)]
        self.muon_data = self.weight_muons(files)
        
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
    
    def define_observable(self, name, definition):
        for event in self.muon_data: event[name] = eval(definition)
        for event in self.numu_data: event[name] = eval(definition)
        for event in self.nuel_data: event[name] = eval(definition)

    def define_observable_from_function(self, name, function):
        for event in self.muon_data: event[name] = function(event)
        for event in self.numu_data: event[name] = function(event)
        for event in self.nuel_data: event[name] = function(event)
    
    def get_histodata(self, particle, observable, bins=None, requirement=None):
    
        #check particle type, assign data:
        particles = {'muon': self.muon_data, 'numu': self.numu_data, 'nuel': self.nuel_data}
        if particle in particles.keys(): data = particles[particle]
        else: print ("Error: particle must be either 'muon', 'numu' or 'nuel'")
        
        # loop through events
        selected_data = []
        for event in data:
            if eval(requirement) == False: continue
            value = event[observable]
            weight = event['weight']
            selected_data.append([value, weight, weight*weight])
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

    ########################################################
    # Binning
    ########################################################
    
    
    

    def binCalEnergy(self, particle, bins = None, requirement=None, weighted=True):

        #check particle type:
        if particle not in ['muon', 'nu_mu', 'nu_el']:
            print ("particle must be either 'muon', 'nu_mu' or 'nu_el'")
            return 1
        
        if particle == 'muon':
            data = self.muon_data
        elif particle == 'nu_mu':
            data = self.numu_data
        elif particle == 'nu_el':
            data = self.nuel_data
        
        func = eval('lambda event: ' + requirement)
        matches = [event for event in data if func(event)]
        if bins is None:
            deposits = [event['calorimeter'] for event in matches]
            minDep = min(deposits) - 1
            maxDep = max(deposits) + 1
            low = np.log10(max(10**-2, minDep)) # Sets 0.01 as the smallest bin
            high = max(3, np.log10(maxDep))     # Sets 1000 as the largest bin
            bins = np.logspace(low, high, 100)

        primEnergies = sorted(set([event['primaryEnergy'] for event in matches]))

        binned = {'bins': bins}

        for energy in primEnergies:
            subdata = [(event['calorimeter'], event['weight']) for event in matches if event['primaryEnergy'] == energy]
            cals    = [x[0] for x in subdata]
            weights = [x[1] for x in subdata]

            counts = np.zeros(np.size(bins))
            uncertainties = np.zeros(np.size(bins))
            for weight in set(weights):
                mask = weights == weight

                cal = np.array(cals)[mask]
                ids = np.searchsorted(bins,cal)

                count = np.zeros(np.size(bins))
                raw_count = np.zeros(np.size(bins))
                for id_ in ids:
                    raw_count[id_] += 1
                    count[id_]     += weight
                    
                if weighted:
                    counts += count
                    mask = raw_count > 0
                    uncertainties[mask] += count[mask] / np.sqrt(raw_count[mask])
                else:
                    counts += raw_count
                    uncertainties += np.sqrt(raw_count)

            binned[energy] = (counts, uncertainties)

        return binned
    
    def sumCalEnergyBins(self, particle, bins = None, requirement = None, weighted = True, primaryEnergies = None):
        
        data = self.binCalEnergy(particle, bins, requirement, weighted)
        
        if primaryEnergies is None:
            primaryEnergies = [key for key in data.keys() if isinstance(key,float)]
        
        bins = data['bins']
        counts = np.zeros(np.shape(bins))
        uncertainties = np.zeros(np.shape(bins))
        for energy in primaryEnergies:
            counts += data[energy][0]
            uncertainties += data[energy][1]
        
        return bins, counts, uncertainties
    
    ########################################################
    # Count events after cuts
    ########################################################
    
    def countEvents(self, particle, requirement="True"):
        
        #check particle type:
        if particle not in ['muon', 'nu_mu', 'nu_el']:
            print ("particle must be either 'muon', 'nu_mu' or 'nu_el'")
            return 0
        
        if particle == 'muon':
            data = self.muon_data
        elif particle == 'nu_mu':
            data = self.numu_data
        elif particle == 'nu_el':
            data = self.nuel_data
        
        func = eval('lambda event: ' + requirement)
        return len([True for event in data if func(event)])
    
    def integrateCounts(self, particle, calEnergies, requirement="True", weighted = True):
    
        #check particle type:
        if particle not in ['muon', 'nu_mu', 'nu_el']:
            print ("particle must be either 'muon', 'nu_mu' or 'nu_el'")
            return 0
        
        if particle == 'muon':
            data = self.muon_data
        elif particle == 'nu_mu':
            data = self.numu_data
        elif particle == 'nu_el':
            data = self.nuel_data
        
        func = eval('lambda event : ' + requirement)
        # Apply an early cut with the scintillator cut so we don't search the same events repeatedly
        subdata = [(event['calorimeter'], event['weight']) for event in data
                   if func(event) and
                      event['calorimeter'] > min(calEnergies)]

        raw_counts = [len(subdata)]
        weighted_counts = [np.sum([w for (e, w) in subdata])]

        # Sort and pop off minimum energy since cut was already made
        calEnergies = sorted(calEnergies)
        calEnergies = calEnergies[1:]

        for energy in calEnergies:
            subdata = [(e,w) for (e,w) in subdata if e > energy]
            raw_counts.append(len(subdata))
            weighted_counts.append(np.sum([w for (e, w) in subdata]))

        raw_counts = np.array(raw_counts)
        weighted_counts = np.array(weighted_counts)
        if weighted:
            mask = raw_counts > 0
            counts = weighted_counts
            uncertainties = np.zeros(np.shape(counts))
            uncertainties[mask] = counts[mask] / np.sqrt(raw_counts[mask])
            return weighted_counts, uncertainties
        else:
            return raw_counts, np.sqrt(raw_counts)
