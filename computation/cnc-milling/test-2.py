import pandas as pd
import numpy as np

df = pd.read_csv('experiment_01.csv')
a=0
b=0
variance = 0
size = df.shape[0]
c= np.random.rand(size,3)
d= np.random.rand(size,3)
for i in range (size):
    a=df.at[i,'X1_ActualPosition'] +a
mean = a/df.shape[0]   
for i in range (size):
    c[i,0] = mean - df.at[i,'X1_ActualPosition']
for i in range (size-1):
    d[i,0]=(c[i+1,0]-c[i,0])**2
    variance = d[i,0] +variance
variance = variance/(size-1)    
#print('X1_ActualPosition =',mean)    
print('variance of X1_ActualPosition=',variance)   

a=0
b=0
variance = 0
#print('a =',a)
for i in range (size):
    a=df.at[i,'Y1_ActualPosition'] +a
mean = a/df.shape[0]   
for i in range (size):
    c[i,1] = mean - df.at[i,'Y1_ActualPosition']
for i in range (size-1):
    d[i,1]=(c[i+1,1]-c[i,1])**2
    variance = d[i,1] +variance
variance = variance/(size-1)    
#print('Y1_ActualPosition =',mean)    
print('variance of Y1_ActualPosition=',variance)   

a=0
b=0
variance = 0
for i in range (size):
    a=df.at[i,'Z1_ActualPosition'] +a
mean = a/df.shape[0]   
for i in range (size):
    c[i,2] = mean - df.at[i,'Z1_ActualPosition']
for i in range (size-1):
    d[i,2]=(c[i+1,2]-c[i,2])**2
    variance = d[i,2] +variance
variance = variance/(size-1)    
#print('Z1_ActualPosition =',mean)    
print('variance of Z1_ActualPosition=',variance)   
