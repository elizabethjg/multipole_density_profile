import sys
sys.path.append('/home/eli/Documentos/PostDoc/halo-elongation/multipole_density_profile')
sys.path.append('/mnt/clemente/lensing/multipole_density_profile')
sys.path.append('/home/elizabeth/multipole_density_profile')
sys.path.append('/home/eli/Documentos/Astronomia/posdoc/halo-elongation/multipole_density_profile')
from make_profile_redMapper import *

# profile_redMapper('total_z04_v3',20.,150.,zmin = 0.1, zmax = 0.4,
                      # z_back = 0.1, odds_min = 0.5, RIN = 100., ROUT = 5000., ndots = 10.)

'''
samples  = ['total_z033','total_z04','total_z045']
lmin     = np.array([20.,20.,20.])
lmax     = np.array([150.,150.,150.])
percentil = [False,False,False]
zmin     = np.ones(len(lmin))*0.1
zmax     = np.array([0.33,0.4,0.45])
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*(100.)
ROUT     = np.ones(len(lmin))*(5000.)
ndots    = np.ones(len(lmin))*10

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,percentil,
                    RIN,ROUT,ndots]).T

for j in range(3):
    makeprofile_unpack(entrada[j])
'''


'''
samples  = ['original_withR_bin1','original_withR_bin2','original_withR_bin3','original_withR_bin4']
lmin     = np.array([20.,23.42,28.3,39.7])
lmax     = np.array([23.42,28.3,39.7,145.])
percentil = [False,False,False,False]
zmin     = np.ones(len(lmin))*0.1
zmax     = np.ones(len(lmin))*0.33
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*(100./0.7)
ROUT     = np.ones(len(lmin))*(10000./0.7)
ndots    = np.ones(len(lmin))*20

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,percentil,
                    RIN,ROUT,ndots]).T

for j in range(4):
    makeprofile_unpack(entrada[j])


#'''

# '''

samples  = ['total','median_bin1','median_bin2',
            'terciles_bin1','terciles_bin2','terciles_bin3',
            'medianz1','medianz2',
            'terciles_z1','terciles_z2','terciles_z3']
lmin     = np.array([20.,20.,27.982,20.,24.763,33.214]+[20.]*5)
lmax     = np.array([150.,27.982,150.,24.763,33.214,150.]+[150.]*5)
zmin     = np.array([0.1]*6+[0.1,0.313,0.1,0.268,0.346])
zmax     = np.array([0.4]*6+[0.313,0.4,0.268,0.346,0.4])
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*100.
ROUT     = np.ones(len(lmin))*5000.
ndots    = np.ones(len(lmin))*10

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,RIN,ROUT,ndots]).T

for j in range(len(entrada)):
    makeprofile_unpack(entrada[j])



# '''
'''
samples  = ['terciles_wR_z33_total','terciles_wR_z33_bin1','terciles_wR_z33_bin2','terciles_wR_z33_bin3']
lmin     = np.array([20.,20.,24.87,33.24])
lmax     = np.array([150.,24.87,33.24,150.])
percentil = [False,[0,1./3.],[1./3.,2./3.],[2./3.,1.]]
zmin     = np.ones(len(lmin))*0.1
zmax     = np.ones(len(lmin))*0.33
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
RIN      = np.ones(len(lmin))*100.
ROUT     = np.ones(len(lmin))*5000.
ndots    = np.ones(len(lmin))*12

entrada = np.array([samples,lmin,lmax,zmin,zmax,
                    z_back,odds_min,percentil,
                    RIN,ROUT,ndots]).T

for j in range(4):
    makeprofile_unpack(entrada[j])

'''
'''
name_cat = ['CFHT']*3 + ['CS82']*3 + ['RCSL']*3 + ['KiDS']*3
samples  = ['total_z033','total_z04','total_z045']*4
lmin     = np.array([20.]*12)
lmax     = np.array([150.]*12)
zmin     = np.ones(len(lmin))*0.1
zmax     = np.array([0.33,0.4,0.45]*4)
z_back   = np.ones(len(lmin))*0.1
odds_min = np.ones(len(lmin))*0.5
percentil = [False,False,False]*4
RIN      = np.ones(len(lmin))*(100.)
ROUT     = np.ones(len(lmin))*(5000.)
ndots    = np.ones(len(lmin))*10
zlim     = np.append(np.ones(9)*1.3,np.ones(3)*0.9)


#name_cat = name_cat + ['CFHT']*4 + ['CS82']*4 + ['KiDS']*4
#samples  = samples +  ['terciles_total','terciles_bin1','terciles_bin2','terciles_bin3']*3
#lmin     = np.append(lmin,np.tile(np.array([20.,20.,28.,40.]),3))
#lmax     = np.append(lmax,np.tile(np.array([150.,28.,40.,150.]),3))
#zmin     = np.append(zmin,np.ones(12)*0.1)
#zmax     = np.append(zmax,np.ones(12)*0.33)
#z_back   = np.append(z_back,np.ones(12)*0.1)
#odds_min = np.append(odds_min,np.ones(12)*0.5)
#RIN      = np.append(RIN,np.ones(12)*200.)
#ROUT     = np.append(ROUT,np.ones(12)*5000.)
#ndots    = np.append(ndots,np.ones(12)*10)
#zlim     = np.append(zlim,np.append(np.ones(8)*1.3,np.ones(4)*0.9))

for j in range(2):
    ini = j*6
    fin = (j+1)*6
    makeindprofile_parallel(name_cat[ini:fin],samples[ini:fin],lmin[ini:fin],
                            lmax[ini:fin],zmin[ini:fin],zmax[ini:fin],
                            z_back[ini:fin],odds_min[ini:fin],percentil[ini:fin],
                            RIN[ini:fin],ROUT[ini:fin],ndots[ini:fin],zlim[ini:fin])
#'''
