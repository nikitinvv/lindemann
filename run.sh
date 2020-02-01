#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -J rec
#SBATCH -e err
#SBATCH --mem 160G
#SBATCH --exclusive
#SBATCH -N 1
#SBATCH -p v100
#SBATCH --exclude gn1
#load the modules required for you program - customise for your program

source ~/.bashrc-2019-04-27

for k in {0..5025..1005};
do
 python rectv.py /data/staff/tomograms/viknik/Lindemann/ 1150 $k $(($k+1005))
 python rectv.py /data/staff/tomograms/viknik/Lindemann/ 1149 $k $(($k+1005))
 python rectv.py /data/staff/tomograms/viknik/Lindemann/ 1136 $k $(($k+1005))
 python rectv.py /data/staff/tomograms/viknik/Lindemann/ 1137 $k $(($k+1005))
done