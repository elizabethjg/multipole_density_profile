import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-folder', action='store', dest='folder',default='./')
parser.add_argument('-ini', action='store', dest='ini',default=0)
parser.add_argument('-fini', action='store', dest='fini', default=10)
args = parser.parse_args()

ini  = args.ini
ini  = int(ini)
fini = args.fini
fini = int(fini)

folder    = args.folder

for j in range(ini,fini):
    
    pname  = 'profile_medianas_'+str(j)+'.cat'
    
    print '---------------------------'
    print pname
    print '###########################'
    
    # correr = 'python -u fit_profile_monopole_misscentred.py -folder '+folder+' -file '+pname+' -ncores 10'
    # os.system(correr)

    f = open(folder+pname,'r')
    lines = f.readlines()
    pcc = float((lines[-1][1:-2]))
    M200 = float((lines[-2][1:-2]))
    f.close()
    
    f=open(folder+pname,'a')
    f.write('#'+str('%.2f' % (10**(M200)/1.e14))+'   \n')
    f.write('#'+str('%.2f' % (pcc))+'   \n')
    f.close()
    
    # os.system('mv /mnt/clemente/lensing/redMaPPer/test/profile_medianas_tp_'+str(j)+'.cat /mnt/clemente/lensing/redMaPPer/test/profile_medianas_'+str(j)+'_tp.cat')
    correr2 = 'python -u fit_profile_quadrupole.py -folder '+folder+' -file  '+pname+' -ncores 10 -ang \'tp\' -misscentred \'False\' -nit 200 -ROUT 700 -component both'
    os.system(correr2)

    print '###########################'
