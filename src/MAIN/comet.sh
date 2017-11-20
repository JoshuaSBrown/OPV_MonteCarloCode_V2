#!/bin/bash
#SBATCH --output=shout.txt
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:20:00
#SBATCH --mail-user=joshbro42867@yahoo.com

present=$(pwd)
folder=$(pwd | awk -F "/" '{print $NF}')
WORK=/oasis/scratch/comet/joshbro/temp_project/
cd $WORK
echo $present

rm -rf $folder
mkdir $folder
rm -rf PARAMETERS
mkdir PARAMETERS
cp $present/../PARAMETERS/parameters.txt PARAMETERS/
cd $folder
cp $present/* .
./run > out.txt
gprof ./run > out2.txt
cp *txt $present
cp *out $present
cd $present
