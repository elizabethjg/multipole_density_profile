#!/bin/bash
#SBATCH --job-name=fit_tcos
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

srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_total.cat' -ncores 10 -ang 'twl'
