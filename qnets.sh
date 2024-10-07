#!/bin/bash

#SBATCH --nodes=1                    # 1 node
#SBATCH --ntasks-per-node=1          # 1 task per node
#SBATCH --cpus-per-task=32
#SBATCH --mem=494000
#SBATCH --time=24:00:00              # time limits: 24 hours
#SBATCH --account=INF24_npqcd_1      # account name
#SBATCH --partition=boost_usr_prod   # partition name
#SBATCH --qos=normal                 # quality of service
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=andrea.degirolamo@poliba.it

export OMP_NUM_THREADS=32
export JULIA_NUM_THREADS=32
julia run.jl $1 $2 $3 -p 32 --threads 32
