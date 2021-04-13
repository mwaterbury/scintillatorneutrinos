import re
import numpy as np
from os import listdir



# Find appropriate files in directory
dir = './NeutrinoRun-0330/'
ls = listdir(dir)
files = [dir + fname for fname in ls if re.match(r'.*00[0-9]_fort.25',fname)]
files.sort()

numSc = 10
numEv = 1000

for fname in files:

    evRecord_tmp = np.zeros([numEv,numSc])

    try:
        f = open(fname)

        line = f.readline()
        while line:
            if re.search(r'Event',line) and re.search(r'mu',line):
                # Extract scintillator and event from line
                m = re.match(r'.*s([0-9]*)MuBin.*Event #: *([0-9]*).*',line)
                try:
                    sNum = int(m[1])
                    eNum = int(m[2])
                except TypeError:
                    print(line)
                    quit()

                if sNum == 10:
                    sNum = 7

                # Read whether bin was hit
                line = f.readline()
                m = re.split(r'  *',line)
                hit = int(re.split(r'\n',m[len(m)-1])[0])
                # Mark as binned in event record
                evRecord_tmp[eNum-1,sNum-1] = hit
            line = f.readline()
    finally:
        f.close()

    try:
        evRecord = np.concatenate((evRecord, evRecord_tmp))
    except NameError:
        evRecord = evRecord_tmp

np.save("fasernu_muonRecord.npy",evRecord)
