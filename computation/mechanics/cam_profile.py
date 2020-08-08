from lib_cam_profile import *

b1, b2, b3, b4 = 50, 120, 90, 100
num, n = 500, 600
omg = n*pi/30
h = 5

n1 = int(num*b1/360)
n2 = int(num*b2/360)
n3 = int(num*b3/360)
n4 = int(num*b4/360)

t1 = np.linspace(0,b1,n1,endpoint=0)*pi/180
t2 = np.linspace(0,b2,n2,endpoint=0)*pi/180
t3 = np.linspace(0,b3,n3,endpoint=0)*pi/180
t4 = np.linspace(0,b4,n4,endpoint=1)*pi/180
b1r, b2r, b3r, b4r = np.array([b1,b2,b3,b4])*np.pi/180
xs = np.hstack([t1, b1r+t2, b1r+b2r+t3, b1r+b2r+b3r+t4])

#%%
s1, ds1, dds1, ddds1 = FullModifiedTrapezoidal(t1, h, b1)
s2 = np.ones(n2) * h
s3, ds3, dds3, ddds3 = h-np.array(FullCycloidal(t3, h, b3))
s4 = np.zeros(n4)
ds2 = np.zeros(n2)

ys = np.hstack([s1, s2, s3, s4])
dys = np.hstack([ds1, ds2, ds3, s4])
ddys = np.hstack([dds1, ds2, dds3, s4])
dddys = np.hstack([ddds1, ds2, ddds3, s4])

#%%
rho_min = 1.5
Rb_min = rho_min - min(ys + ddys)
print(Rb_min)
Rb = 30
print(max(dys)-min(dys))

s, x, theta = ys, np.hstack((0,np.diff(ys))), xs

RA = np.sqrt((Rb+s)**2 + x**2)
sigA = theta + np.arctan(x/(Rb + s))

xn = RA * np.cos(sigA)
yn = -np.sign(omg) * RA * np.sin(sigA)

fig, ax = plt.subplots(figsize=(10,8))
ax.plot(xn,yn,'k',label='cam profile')
ax.grid()
ax.set_aspect('equal')
ax.legend()
plt.tight_layout()

# np.savetxt('cam_profile.txt',np.array([xn,yn,np.zeros(len(xn))]).T,delimiter='\t')