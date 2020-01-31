#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -J rec
#SBATCH -e err
#SBATCH --mem 160G
#SBATCH --exclusive
#SBATCH -N 1
#load the modules required for you program - customise for your program

source ~/.bashrc-2019-04-27

# for k in {0..5025..1005};
# do
     #python preprocess.py /data/staff/tomograms/viknik/Lindemann/Continuous_2D_8400eV_500ms_2_1137.h5 $k
python preprocess.py /data/staff/tomograms/viknik/Lindemann/Continuous_1D_8400eV_100ms_2_1157.h5 0
# done