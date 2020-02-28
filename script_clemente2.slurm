#!/bin/bash
#SBATCH --job-name=total
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/multipole_density_profile

export OMP_NUM_THREADS=20
export MKL_NUM_THREADS=20

source activate py2env

srun python -u run_quadrupole_analyse.py -file 'table_mcmc_total.out' -profile 'profile_total.cat'
