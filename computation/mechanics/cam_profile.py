from lib_cam_profile import *

"""
USABLE ONLY FOR DOUBLE DWELL CAMS
"""

u = UnitRegistry()
slant = 6 * u.deg
b1, b2, b3 = np.array([65,20,65])*u.deg
b4 = 360*u.deg - b1 - b2 - b3
n = 500*u.rpm
h = 2*np.cos(slant.to(u.rad)) * u.mm

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
s4 = np.zeros(n4) * h.units

ds2 = (s2-s2)/u.rad
dds2, ddds2 = ds2/u.rad**2, ds2/u.rad**3
ds4, dds4, ddds4 = s4/u.rad, s4/u.rad**2, s4/u.rad**3

ys = np.hstack([s1, s2, s3, s4])
dys = np.hstack([ds1, ds2, ds3, ds4])
ddys = np.hstack([dds1, dds2, dds3, dds4])
dddys = np.hstack([ddds1, ddds2, ddds3, ddds4])

#%%
PlotSVAJ(xs, ys, dys, ddys, dddys, n, savefig=True)
# RadialCamTranslatingFlatFaced(xs, ys, dys, ddys, n, 1, 12, 10, show=True, savefig=True)
# RadialCamTranslatingRoller(xs, ys, n, 16.5, 0, show=True)

#%%

# including some errors for signal analysis (butterworth, etc.)
# m = 1e3*1.5
# ys = np.hstack([s1, s2, s3, s4])
# e1 = np.random.default_rng().standard_normal(len(ys))/m * u.mm
# ys = ys + e1
# e2 = np.random.default_rng().standard_normal(len(ys)-1)/m**2 * u.mm
# dys = (np.diff(ys) + e2*0)/np.diff(xs.to(u.rad))
# e3 = np.random.default_rng().standard_normal(len(dys)-1)/m * u.mm / u.rad
# ddys = (np.diff(dys) + e3*0)/np.diff(xs.to(u.rad))[1:]
# e4 = np.random.default_rng().standard_normal(len(ddys)-1)/m * u.mm / u.rad**2
# dddys = (np.diff(ddys) + e4*0)/np.diff(xs.to(u.rad))[2:]
# xs = xs[3:]
# ys = ys[3:]
# dys = dys[2:]
# ddys = ddys[1:]
# np.savetxt('cam_profile.txt',np.array([xs.magnitude, ys.magnitude, dys.magnitude,
#                         ddys.magnitude, dddys.magnitude]).T,delimiter='\t')