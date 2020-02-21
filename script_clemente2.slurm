#!/bin/bash
#SBATCH --job-name=total_t
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/multipole_density_profile

export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

source activate py2env

srun python -u run_quadrupole_analyse.py -file 'table_mcmc_centred.out' \
-profile 'profile_total.cat' -misscentred 'False'

##srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_total.cat' -ncores 10 -ang 't' -misscentred 'True' -nit 100
