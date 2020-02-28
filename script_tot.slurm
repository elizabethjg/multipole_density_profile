#!/bin/bash
#SBATCH --job-name=t_int_both
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /home/elizabeth/multipole_density_profile

export OMP_NUM_THREADS=10

module load softwares/anaconda2/5.0-intel-2018.0
source activate py2env

srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 't' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'twl' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'twd' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'tp' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'tpwl' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'tpwd' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
srun python -u fit_profile_quadrupole.py -folder '/home/elizabeth/profiles_redMapper/' \
-file 'profile_total.cat' -ncores 10 -ang 'control' -misscentred 'False' -nit 200 -ROUT 5000 -component 'both'
