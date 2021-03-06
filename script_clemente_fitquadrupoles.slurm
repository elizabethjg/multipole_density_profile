#!/bin/bash
#SBATCH --job-name=fit_total
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

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 't' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'twl' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'twd' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tp' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwl' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwd' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwdl' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'control' -misscentred 'False' -component 'tcos'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 't' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'twl' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'twd' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tp' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwl' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwd' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'tpwdl' -misscentred 'False' -component 'xsin'

srun python -u fit_profile_quadrupole.py -folder '/mnt/clemente/lensing/redMaPPer/' \
-file 'profile_terciles_total.cat' -ncores 10 -ang 'control' -misscentred 'False' -component 'xsin'
