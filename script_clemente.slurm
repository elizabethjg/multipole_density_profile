#!/bin/bash
#SBATCH --job-name=qfit_b1tp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/multipole_density_profile

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

source activate py2env

srun python -u make_profile_redMapper.py

#srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper' \
-file 'profile_median_bin1.cat' -ncores 10 -ang 'tp' -misscentred 'True' -nit 100 -RIN 1
