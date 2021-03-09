import re
import numpy as np
from os import listdir



# Find appropriate files in directory
ls = listdir()
files = [fname for fname in ls if re.match(r'.*00[2-3]_fort.24',fname)]

for fname in files:

    numSc = 9
    numEv = 100000
    evRecord_tmp = np.zeros([numEv,numSc])
    evEnergy_tmp = np.zeros([numEv,numSc])

    try:
        f = open(fname)

        line = f.readline()
        while line:
            if re.search(r'Event',line):
                # Extract scintillator and event from line
                m = re.match(r'.*s([0-9])evBin.*Event #: *([0-9][0-9]*).*',line)
                sNum = int(m[1])
                eNum = int(m[2])

                # Read whether bin was hit
                line = f.readline()
                m = re.split(r'  *',line)
                hit = int(re.split(r'\n',m[len(m)-1])[0])
                # Mark as binned in event record
                evRecord_tmp[eNum-1,sNum-1] = hit

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
                    evEnergy_tmp[eNum-1,sNum-1] = E
            line = f.readline()
    finally:
        f.close()

    try:
        evRecord = np.concatenate((evRecord, evRecord_tmp))
        evEnergy = np.concatenate((evEnergy, evEnergy_tmp))
    except NameError:
        evRecord = evRecord_tmp
        evEnergy = evEnergy_tmp

np.save("fasernu_eventRecord.npy",evRecord)
np.save("fasernu_eventEnergy.npy",evEnergy)
