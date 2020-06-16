import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
exps = [pd.read_csv('experiment_01.csv'),
        pd.read_csv('experiment_02.csv'),
        pd.read_csv('experiment_03.csv'),
        pd.read_csv('experiment_04.csv'),
        pd.read_csv('experiment_05.csv'),
        pd.read_csv('experiment_06.csv'),
        pd.read_csv('experiment_07.csv'),
        pd.read_csv('experiment_08.csv'),
        pd.read_csv('experiment_09.csv'),
        pd.read_csv('experiment_10.csv'),
        pd.read_csv('experiment_11.csv'),
        pd.read_csv('experiment_12.csv'),
        pd.read_csv('experiment_13.csv'),
        pd.read_csv('experiment_14.csv'),
        pd.read_csv('experiment_15.csv'),
        pd.read_csv('experiment_16.csv'),
        pd.read_csv('experiment_17.csv'),
        pd.read_csv('experiment_18.csv')]
unworn = [0,1,2,3,4,10,11,16]
worn = [5,6,7,8,9,12,13,14,15,17]
mean_x, mean_y, mean_z ,mean_t= [], [], [], []
mean_x_unworn, mean_y_unworn, mean_z_unworn = [], [], []
mean_x_worn, mean_y_worn, mean_z_worn, = [], [], []
mean_x1=0
mean_y1=0
mean_z1=0

dfw=np.zeros(18)
a=0
b=0
c=0
len1=0
len2=0
len3=0
sst=0
ssw=0

for i in range(0,18):
    err_pos_x = np.asarray(exps[i]['X1_ActualPosition']-exps[i]['X1_CommandPosition'])
    err_pos_y = np.asarray(exps[i]['Y1_ActualPosition']-exps[i]['Y1_CommandPosition'])
    err_pos_z = np.asarray(exps[i]['Z1_ActualPosition']-exps[i]['Z1_CommandPosition'])
    mean_x.append(abs(err_pos_x))
    mean_y.append(abs(err_pos_y))
    mean_z.append(abs(err_pos_z))

for i in range(len(mean_x)):
    for j in range(mean_x[i].size):
        a=a+mean_x[i][j]
        b=b+mean_y[i][j]
        c=c+mean_z[i][j]
    len1=len1+mean_x[i].size
    len2=len2+mean_y[i].size
    len3=len3+mean_z[i].size
#xgm= (18exps . (x)k obs)/(total numebr of obs of 18 exps)
xgm=(a+b+c)/(len1+len2+len3)
#x1,x2,x3 la mean cua delta pos x,y,z (18exps . (x)k obs/ itself)
x1=a/len1
x2=b/len2
x3=c/len3
for i in range(18):
    for j in range(mean_x[i].size):
        sst=sst+((mean_x[i][j]-xgm)**2)+((mean_y[i][j]-xgm)**2)+((mean_z[i][j]-xgm)**2)
        ssw=ssw+((mean_x[i][j]-x1)**2)+((mean_y[i][j]-x2)**2)+((mean_z[i][j]-x3)**2)
ssb=len1*((x1-xgm)**2 + (x2-xgm)**2 + (x3-xgm)**2)
dfb=2
dfw=len1*3-3
dft=len1*3-1
msb=ssb/dfb
msw=ssw/dfw

print('xgm=', xgm)
print('SSB=',ssb)
print('SST=',sst)
print('dfb=',dfb)
print('dfw=',dfw)
print('dft=',dft)
print('MSB=',msb)
print('MSW=',msw)
print('F=',msb/msw)



