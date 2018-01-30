#!/bin/bash -l

#SBATCH --job-name=lidarMC
#SBATCH --output=lidarMC.txt
#SBATCH --ntasks=1

./lidar
