import re
import numpy as np

numTotPartPerRun = 100000
numRuns = 1
beamWidth = 50.0
RNG = 0.0

# Generates a string with the appropriate amount of whitespace
def str2WHAT(input):
    WHAT = input
    while len(WHAT) < 10:
        WHAT = ' ' + WHAT
    return WHAT

# Generates a string with the correct formatting for the beam information
def genBeamInfoLine(energy,dx,dy,particle):
    BEAM = 'BEAM      '
    # - here sets energy instead of momentum of the particle
    WHAT2 = str2WHAT('-' + str(energy))
    WHAT34 = '                    '
    WHAT5 = str2WHAT(str(dx))
    WHAT6 = str2WHAT(str(dy))
    # 1. here sets a rectangle
    WHAT7 = str2WHAT('1.')
    WHAT8 = particle + '\n'
    newLine = BEAM + WHAT2 + WHAT34 + WHAT5 + WHAT6 + WHAT7 + WHAT8
    return newLine

def genPartNumLine(numPart):
    return 'START     ' + str2WHAT(str(int(numPart))) + '\n'

for filename in ['fasernu_bot', 'fasernu_left', 'fasernu_right', 'fasernu_top']:
# for filename in ['fasernu_bot']:
    inputFile = open(filename + '.inp')
    nmuData = np.loadtxt('../MuonFluenceFiles/negative_muon_flux.csv')
    # pmuData = np.loadtxt('../MuonFluenceFiles/positive_muon_flux.csv')

    nmuEn = nmuData[:,0]
    nmuFl = nmuData[:,1]
    # nmuDist = nmuFl/np.sum(nmuFl)

    # Create a few file for each muon charge and energy
    nmuOut = []
    for en in nmuEn:
        nmuOut.append(open(str(int(np.floor(en))) + '_' + filename + '.inp','w'))

    # pmuEn = pmuData[:,0]
    # pmuFl = pmuData[:,1]
    # pmuDist = pmuFl/np.sum(pmuFl)
    #
    # pmuOut = []
    # for en in pmuEn:
    #     pmuOut.append(open('pos_mu_' + str(int(np.floor(en))) + '.inp','w'))

    # Write each new file line by line
    for line in inputFile:
        for n in range(np.size(nmuEn)):
            outputFile = nmuOut[n]
            if re.search(r'BEAM ',line):
                m = re.split(r'  *',line)
                energy = nmuEn[n]
                dx = m[2]
                dy = m[3]
                newLine = genBeamInfoLine(energy, dx, dy, 'MUON-')
            # elif re.search('START',line):
                # numPart = numTotPartPerRun * nmuDist[n]
                # newLine = genPartNumLine(numPart)
            elif re.search(r'RANDOMIZ ',line):
                newLine = 'RANDOMIZ  ' + str2WHAT(str(1.0)) + str2WHAT(str(RNG)) + '\n'
                RNG = RNG + 1.0
            else:
                newLine = line
            outputFile.write(newLine)

        # for n in range(np.size(pmuEn)):
        #     outputFile = pmuOut[n]
        #     if re.search(r'BEAM ',line):
        #         energy = pmuEn[n]
        #         newLine = genBeamInfoLine(energy, beamWidth, 'MUON+')
        #     elif re.search('START',line):
        #         numPart = numTotPartPerRun * pmuDist[n]
        #         newLine = genPartNumLine(numPart)
        #     else:
        #         newLine = line
        #     outputFile.write(newLine)

    inputFile.close()
