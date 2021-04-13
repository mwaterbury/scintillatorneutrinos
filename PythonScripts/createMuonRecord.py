import re
import numpy as np
from os import listdir



# Find appropriate files in directory
dir = './'
ls = listdir(dir)
files = [dir + fname for fname in ls if re.match(r'.*001_fort.25',fname)]
files.sort()

numSc = 10
evRecord = False

for fname in files:
    try:
        f = open(fname)

        line = f.readline()
        while line:
            if re.search(r'Event',line):
                # Extract scintillator and event from line
                m = re.match(r'.*s([0-9]*)MuBin.*Event #: *([0-9]*).*',line)
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

                # Record detector into event history
                if sNum == 1:
                    evHist = np.array(hit)
                else:
                    evHist = np.append(evHist,hit)
                # If last scintillator, add event history to record
                if sNum == numSc:
                    if not(evRecord):
                        evRecord = [evHist]
                    else:
                        evRecord = np.concatenate((evRecord, [evHist]))
            line = f.readline()
    finally:
        f.close()

np.save("fasernu_muonRecord.npy",evRecord)
