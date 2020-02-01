#!/bin/bash
#SBATCH -t 2-00:00:00
#SBATCH -J rec
#SBATCH -e err2
#SBATCH --mem 160G
#SBATCH --exclusive
#SBATCH -N 1
#SBATCH -p v100
#SBATCH --exclude gn1
#load the modules required for you program - customise for your program

source ~/.bashrc-2019-04-27

python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1127 180 1185 1169
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1127 1185 2190 1173
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1127 2190 3195 1169
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1127 3195 4200 1169
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1127 4200 5205 1169

python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1128 180 1185 1164
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1128 1185 2190 1172
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1128 2190 3195 1164
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1128 3195 4200 1172
python rectv_1127_1128.py /data/staff/tomograms/viknik/Lindemann/ 1128 4200 5205 1164