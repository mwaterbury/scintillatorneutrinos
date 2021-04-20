#!/bin/zsh

# Set FLUKA directory
export FLUPRO=~/.local/fluka/

# Get input files for runs
# FILES=$(ls . | grep -G '[0-9][0-9]*_fasernu_.*inp')

# Input files for specific energies
ENERGIES=(327 3321 3521 524 721 924)
ENERGIES=(121 327 524 721 924 1124 1324 1724 2124 2324 2521 2921 3321 3521)
for ENERGY in $ENERGIES; do
  FILES=($FILES $(ls . | grep -G '[0-9][0-9]*_fasernu_.*inp' | grep $ENERGY))
done

MAXCPU=5
echo $FILES
for FILE in $FILES; do
    ((i=i%$MAXCPU)); ((i++==0)) && wait
    rfluka -N3 -M5 $FILE &
done
