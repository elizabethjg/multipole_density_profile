#!/bin/bash
#SBATCH --job-name=misscentred
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
## An email is sent when the job is canceled or finished
## Change "root" for your e-mail address
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=ejgonzalez@unc.edu.ar


cd /home/elizabeth/multipole_density_profile

export OMP_NUM_THREADS=20

module load softwares/anaconda2/5.0-intel-2018.0
source activate py2env

srun python -u run_parallel.py
