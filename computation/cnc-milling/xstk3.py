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
err_xs, err_ys, err_zs ,err_t= [], [], [], []
ubs = []
u_x_unworn, mean_y_unworn, mean_z_unworn = [], [], []
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
xcs, ycs, zcs = np.array([]), np.array([]), np.array([])
for exp in exps:
    err_pos_x = np.asarray(exp['X1_ActualPosition']-exp['X1_CommandPosition'])
    err_pos_y = np.asarray(exp['Y1_ActualPosition']-exp['Y1_CommandPosition'])
    err_pos_z = np.asarray(exp['Z1_ActualPosition']-exp['Z1_CommandPosition'])
    err_xs.append(err_pos_x)
    err_ys.append(err_pos_y)
    err_zs.append(err_pos_z)
    xcs = np.concatenate((xcs, abs(err_pos_x)))
    ycs = np.concatenate((ycs, abs(err_pos_y)))
    zcs = np.concatenate((zcs, abs(err_pos_z)))
uws = np.asarray([np.mean(xcs), np.mean(ycs), np.mean(zcs)])
ut = np.mean([xcs, ycs, zcs])
for err_x, err_y, err_z in zip(err_xs,err_ys,err_zs):
    ubs.append([np.mean(err_x), np.mean(err_y), np.mean(err_z)])
ubs = np.asarray(ubs)
# for
# #xgm= (18exps . (x)k obs)/(total numebr of obs of 18 exps)
# xgm=(a+b+c)/(len1+len2+len3)
# #x1,x2,x3 la mean cua delta pos x,y,z (18exps . (x)k obs/ itself)
# x1=a/len1
# x2=b/len2
# x3=c/len3
# for i in range(18):
#     for j in range(mean_x[i].size):
#         sst+=((mean_x[i][j]-xgm)**2)+((mean_y[i][j]-xgm)**2)+((mean_z[i][j]-xgm)**2)
#         ssw+=((mean_x[i][j]-x1)**2)+((mean_y[i][j]-x2)**2)+((mean_z[i][j]-x3)**2)
# ssb=(x1-xgm)**2+(x2-xgm)**2+(x3-xgm)**2
# dfb=2
# dfw=54-3
# dft=54-1
# msb=ssb/dfb
# msw=ssw/dfw

# print('xgm=', xgm)
# print('SSB=',ssb)
# print('SST=',sst)
# print('dfb=',dfb)
# print('dfw=',dfw)
# print('dft=',dft)
# print('MSB=',msb)
# print('MSW=',msw)
# print('F=',msb/msw)



