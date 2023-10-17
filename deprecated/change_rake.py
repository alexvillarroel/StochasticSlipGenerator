import os
import numpy as np
os.chdir('../Slab/')
var=np.genfromtxt('sam_rake.xyz',delimiter=' ')
file=open('sam_rake2.xyz','w')
# this program change the order of columns and change the lon in -180,180
for i in range(np.shape(var)[0]):
    file.write(str(np.round(var[i,1]-360,decimals=2))+','+str(var[i,0])+','+str(var[i,2])+'\n')