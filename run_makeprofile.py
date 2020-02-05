import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
from make_profile_redMapper import *

'''
samples  = ['original_bin1','original_bin2','original_bin3','original_bin4']
lmin     = np.array([20.,23.42,28.3,39.7])
lmax     = np.array([23.42,28.3,39.7,145.])
zmin     = np.ones(len(lmin))*0.1
zmax     = np.ones(len(lmin))*0.33
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*(100./0.7)
ROUT     = np.ones(len(lmin))*(10000./0.7)
ndots    = np.ones(len(lmin))*20

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,RIN,ROUT,ndots]).T

for j in range(3):
    makeprofile_unpack(entrada[j+1])


'''

samples  = ['terciles_total','terciles_bin1','terciles_bin2','terciles_bin3']
lmin     = np.array([20.,20.,24.87,33.24])
lmax     = np.array([150.,24.87,33.24,150.])
percentil = [False,[0,1./3.],[1./3.,2./3.],[2./3.,1.]]
zmin     = np.ones(len(lmin))*0.1
zmax     = np.ones(len(lmin))*0.4
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*100.
ROUT     = np.ones(len(lmin))*5000.
ndots    = np.ones(len(lmin))*10

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,percentil,
                    RIN,ROUT,ndots]).T

for j in range(len(entrada)):
    makeprofile_unpack(entrada[j])

'''

name_cat = ['CFHT']*4 + ['CS82']*4 + ['KiDS']*4
samples  = ['original_bin1','original_bin2','original_bin3','original_bin4']*3
lmin     = np.tile(np.array([20.,23.42,28.3,39.7]),3)
lmax     = np.tile(np.array([23.42,28.3,39.7,145.]),3)
zmin     = np.ones(len(lmin))*0.1
zmax     = np.ones(len(lmin))*0.33
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*(100./0.7)
ROUT     = np.ones(len(lmin))*(10000./0.7)
ndots    = np.ones(len(lmin))*20
zlim     = np.append(np.ones(8)*1.3,np.ones(4)*0.9)


name_cat = name_cat + ['CFHT']*4 + ['CS82']*4 + ['KiDS']*4
samples  = samples +  ['terciles_total','terciles_bin1','terciles_bin2','terciles_bin3']*3
lmin     = np.append(lmin,np.tile(np.array([20.,20.,28.,40.]),3))
lmax     = np.append(lmax,np.tile(np.array([150.,28.,40.,150.]),3))
zmin     = np.append(zmin,np.ones(12)*0.1)
zmax     = np.append(zmax,np.ones(12)*0.33)
z_back   = np.append(z_back,np.ones(12)*0.1)
odds_min = np.append(odds_min,np.ones(12)*0.5)
RIN      = np.append(RIN,np.ones(12)*200.)
ROUT     = np.append(ROUT,np.ones(12)*5000.)
ndots    = np.append(ndots,np.ones(12)*10)
zlim     = np.append(zlim,np.append(np.ones(8)*1.3,np.ones(4)*0.9))

for j in range(3):
    ini = j*8
    fin = (j+1)*8
    makeindprofile_parallel(name_cat[ini:fin],samples[ini:fin],lmin[ini:fin],
                            lmax[ini:fin],zmin[ini:fin],zmax[ini:fin],
                            z_back[ini:fin],odds_min[ini:fin],RIN[ini:fin],
                            ROUT[ini:fin],ndots[ini:fin],zlim[ini:fin])
#'''
