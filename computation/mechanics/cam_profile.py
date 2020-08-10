from lib_cam_profile import *

"""
USABLE ONLY FOR DOUBLE DWELL CAMS
"""

u = UnitRegistry()
slant = 6 * u.deg
b1, b2, b3 = np.array([55,5,55])*u.deg
b4 = 360*u.deg - b1 - b2 - b3
n = 500*u.rpm

h = 3*np.cos(slant.to(u.rad)) * u.mm

num = 1500
n1 = int(num*b1/u.deg/360)
n2 = int(num*b2/u.deg/360)
n3 = int(num*b3/u.deg/360)
n4 = int(num*b4/u.deg/360)

t1 = np.linspace(0,b1,n1,endpoint=0)
t2 = np.linspace(0,b2,n2,endpoint=0)
t3 = np.linspace(0,b3,n3,endpoint=0)
t4 = np.linspace(0,b4,n4,endpoint=1)
xs = np.hstack([t1, b1+t2, b1+b2+t3, b1+b2+b3+t4])

s1, ds1, dds1, ddds1 = RDFD(t1, h, b1, 'modified trapezoidal')
s2 = np.ones(n2) * h
s3, ds3, dds3, ddds3 = RDFD(t3, h, b3, 'modified trapezoidal', rise=False)
s4 = np.zeros(n4)
ds2 = np.zeros(n2)

ys = np.hstack([s1, s2, s3, s4])
dys = np.hstack([ds1, ds2, ds3, s4])/u.rad
ddys = np.hstack([dds1, ds2, dds3, s4])/u.rad**2
dddys = np.hstack([ddds1, ds2, ddds3, s4])/u.rad**3

#%%
PlotSVAJ(xs, ys, dys, ddys, dddys, n)
# RadialCamTranslatingFlatFaced(xs, ys, dys, ddys, n, 1, 15, 21, savefig=True)
# RadialCamTranslatingRoller(xs, ys, n, 16.5, 0, show=True)
