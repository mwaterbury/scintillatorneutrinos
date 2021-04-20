import re
import numpy as np
from os import listdir



# Find appropriate files in directory
dir = './'
ls = listdir(dir)
files = [dir + fname for fname in ls if re.match(r'.*001_fort.24',fname)]
files.sort()

numSc = 10
evRecord_initiated = False

for fname in files:
    print(fname)
    m = re.match(r'./fasernu-([0-9][0-9]*)001.*',fname)
    filenum = int(m[1])
    # print(filenum)
    if filenum < 10:
        zpos = 13*filenum - 2
    elif (filenum > 9) & (filenum < 15):
        zpos = 136 + 2*(filenum-10)
    elif filenum > 14:
        zpos = 573 + 2*(filenum-15)
    # print(zpos)


    try:
        f = open(fname)

        line = f.readline()
        while line:
            # if re.search(r'Energy=',line):
            #     # Extract primary energy
            #     m = re.match(r'.*\[Energy=(.*)\].*',line)
            #     primaryEn = float(m[1])

            if re.search(r'Event',line):
                # Extract scintillator and event from line
                m = re.match(r'.*s([0-9]*)EnBin.*Event #: *([0-9]*).*',line)
                try:
                    sNum = int(m[1])
                    eNum = int(m[2])
                except TypeError:
                    print(line)
                    quit()

                # Read whether bin was hit
                line = f.readline()
                m = re.split(r'  *',line)
                hit = int(re.split(r'\n',m[len(m)-1])[0])

                if hit:
                    # Skip lines to energy line
                    line = f.readline()
                    # Extract energy from line
                    m = re.split(r'  *',line)
                    # Separate amplitude from exponent from scientific notation
                    if re.search(r'E',m[len(m)-1]):
                        [amp, exp] = re.split(r'E',m[len(m)-1])
                        # Remove newline character from exponent
                        exp = re.split(r'\n',exp)[0]
                        E = float(amp) * 10 ** float(exp)
                    else:
                        E = float(m[len(m)-2])
                else:
                    # Set energy deposited to 0 if there is no hit
                    E = float(0)

                # Record detector into event history
                if sNum == 1:
                    evHist = np.array((zpos, hit))
                    # enHist = np.array((primaryEn,E))
                    enHist = np.array(E)
                else:
                    evHist = np.append(evHist,hit)
                    enHist = np.append(enHist,E)
                # If last scintillator, add event history to record
                if sNum == numSc:
                    eventData = np.append(evHist,enHist)
                    if not(evRecord_initiated):
                        evRecord = eventData
                        evRecord_initiated = True
                    else:
                        evRecord = np.concatenate((evRecord, eventData))
            line = f.readline()
    finally:
        f.close()

np.save("fasernu_neutrinoHistories.npy",evRecord)
