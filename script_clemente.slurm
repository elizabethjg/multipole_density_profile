#!/bin/bash
#SBATCH --job-name=0_10 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-0:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /mnt/clemente/lensing/multipole_density_profile

export OMP_NUM_THREADS=10
export MKL_NUM_THREADS=10

source activate py2env

srun python -u fit_profile_monopole_misscentred.py -folder '/mnt/clemente/lensing/redMaPPer/' \\
 -file 'profile_pru_highmass.cat' -ncores 10
 
srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' -file  'profile_pru_highmass.cat' \\
-ncores 10 -ang 'tp' -misscentred 'False' -nit 200 -ROUT 700

##srun python -u run_fit.py -folder '/mnt/clemente/lensing/redMaPPer/test/' -ini 0 -fini 10

##srun python -u make_profile_redMapper.py

##srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper' -file 'profile_median_bin1.cat' -ncores 10 -ang 'tp' -misscentred 'True' -nit 100 -RIN 1
