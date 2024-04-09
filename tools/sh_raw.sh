#!/bin/bash
#SBATCH_job -J sub_1
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH_core -c 3

source /public1/soft/modules/module.sh
module load anaconda
cd ~/johnson/RigidMultiblobsWall-master/multi_bodies/examples/static_massive_height_test
python