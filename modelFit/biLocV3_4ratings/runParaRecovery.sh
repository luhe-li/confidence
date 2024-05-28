#!/bin/bash
#
#SBATCH --job-name=ParaRecovery
#SBATCH -a 1-10 # modify based on looping variable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=32GB
#SBATCH --time=40:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=ll3981@nyu.edu
#SBATCH --output=slurm%j.out

module purge
module load matlab/2022a

matlab <<EOF

paraRecovery

EOF
