#!/bin/zsh

# Set FLUKA directory
export FLUPRO=~/.local/fluka/

# Get input files for runs
# FILES=$(ls . | grep -G 'fasernu-[0-9][0-9]*.inp')

# Input files for specific energies
#ENERGIES=(327 3321 3521 524 721 924)
#for ENERGY in $ENERGIES; do
#  FILES=($FILES $(ls . | grep -G '[0-9][0-9]*_fasernu_.*inp' | grep $ENERGY))
#done

# Input files
INTS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21)
for INT in $INTS; do
   FILES=($FILES $(ls . | grep -G '.*fasernu.*inp' | grep $INT))
done

MAXCPU=5
echo $FILES
for FILE in $FILES; do
    ((i=i%$MAXCPU)); ((i++==0)) && wait
    rfluka -N0 -M1 $FILE &
done
