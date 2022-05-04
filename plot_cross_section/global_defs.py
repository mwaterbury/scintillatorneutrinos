# -*- coding: utf-8 -*-

#################################################################################
#										
# global_defs.py						
#										
# Contains global definitions, constants, and tools.
#										
# Created: 21/03/2016 13:47				
# Last modified: 22/03/2016 01:33
#										
#################################################################################


# User-defined colors

myOrange = '#FFAA00'
myBlue = '#5e81b5'
myLightBlue = '#DEF0FF'
myPink = '#FF8C66'
myRed = '#EB6235'
myLightRed = '#FFAAAA'
myDarkRed = '9c402b'
myGreen = '#8FB032'
myDarkGreen = '#4A6500'
myLightPurple = '#c681ad'
myPurple = '#ac4283'
myBrightPurple = '#713059'
myYellow = '#F8DB00'

# Earth's radius

REarth = 6371.0 # [km]

# Hubble horizon

lengthHubble = 3.89 # [Gpc]

# Average electron fraction inside the Earth

avgElectronFraction = 0.494

# Masses

massElectron = 0.510998928e-3 # [GeV]
massMuon = 105.6583715e-3 # [GeV]
massProton = 938.272046e-3 # [GeV]
massNeutron = 939.565379e-3 # [GeV]
massIsoNucleon = ( massProton+massNeutron ) / 2.0 # [GeV]
massW = 80.385 # [GeV]

# Decay widths

widthW = 2.085 # [GeV]

# Branching ratios

branchWToHad = 0.6741
branchWToMuon = 0.1063

# Fermi constant

GF = 1.1663787e-5 # [GeV^-2]

# Light speed 

speedLightKmPerS = 299792458.e-3 # [km s^-1]

# Average inelasticity <y> for neutrino-nucleon DIS, above Enu ~ 100 TeV

avgInelasticity = 0.25

# Density of ice

densityIce = 0.92 # [g cm^-3]

# Molar density of ice

molar_density_ice = 1./18.01528 # [mol g^-1]

# Effective volume of IceCube

Veff = 1.e15 # [cm^3]

# Conversion factors

convGeVToGram =  1.782661907e-24 # GeV --> g
convErgToGeV = 624.151 # erg --> GeV
convInvGeVToCm = 1.98e-14 # GeV^-1 --> cm
convGeVToInvSec = 2.41799e23 # GeV --> s^-1
convGpcToInvGeV = 2.48876e40 # Gpc --> GeV^-1

# Avogadro number

NAv = 6.022140857e23 # [mol^-1]

# Random seed

ranSeed = 1234

# Cosmology

OmegaM = 0.27 # Adimensional matter energy fraction
OmegaL = 0.73 # Adimensional vacii, energy fraction
OmegaK = 0.00 # Adimensional curvature energy fraction
h = 0.678 
H0 = 100.0 * h # Hubble constant [km s^-1 Mpc^-1]
# lengthHubble = speedLightKmPerS / H0 / 1.e3 # [Gpc] # need to delete the definition of lengthHubble at the top


