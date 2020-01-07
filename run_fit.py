import sys
sys.path.append('/home/elizabeth/multipole_density_profile')
from fit_profile import *

fit_profile_monopole_onlymass('profile_bin4.cat',folder='/home/elizabeth/',ncores=30)
fit_profile_monopole_misscentred('profile_bin4.cat',folder='/home/elizabeth/',ncores=30)
