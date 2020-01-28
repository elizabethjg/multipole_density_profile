#!/bin/bash
#SBATCH --job-name=monopole_miss
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

srun python -u fit_profile_monopole_misscentred_parallelr.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_original_bin4.cat' -ncores 20

##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_twl.cat' -ncores 56
##
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_twl.cat' -ncores 56
##
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_tcos.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_twl.cat' -ncores 56
##
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin1_twl.cat' -ncores 56
##
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin2_twl.cat' -ncores 56
##
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_t.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tp.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwdl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_tpwl.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_twd.cat' -ncores 56
##srun python -u fit_profile_quadrupole_xsin.py -folder '/mnt/clemente/lensing/redMaPPer/' -file 'profile_bin3_twl.cat' -ncores 56
