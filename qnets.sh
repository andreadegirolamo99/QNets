#!/bin/bash

#SBATCH --nodes=1                    # 1 node
#SBATCH --ntasks-per-node=1          # 1 task per node
#SBATCH --cpus-per-task=24
#SBATCH --mem=375300
#SBATCH --time=24:00:00              # time limits: 24 hours
#SBATCH --account=INF24_npqcd_0      # account name
#SBATCH --partition=g100_usr_prod   # partition name
#SBATCH --qos=normal                 # quality of service
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=andrea.degirolamo@poliba.it

export OMP_NUM_THREADS=24
export JULIA_NUM_THREADS=24
julia run.jl $1 $2 $3 $4 $5 $6 -p 24 --threads 24
