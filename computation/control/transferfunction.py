from control import *
import numpy as np
import matplotlib.pyplot as plt
def Damp(zeta, lower_bound=-100, upper_bound=100):
    x = np.linspace(lower_bound,upper_bound,100)
    y = np.tan(np.pi-np.arccos(zeta))*x
    return x,y
#%%
s = tf('s')
sys = lambda k, g: step_info(feedback(k*g,1))
b0,b1,b2,b3 = 0,-1,-6,-9
g = tf(1,np.poly([b0,b1,b2,b3]))
sisotool(g)
# g2 = tf([1,3.064],np.poly([0,-4,-6]))
# f = feedback(44.59*g,1)
# # rlist, klist = rlocus(g)

# # plt.plot(*Damp(0.4923))

# p = -1.9748+5.404j
# zc = -5
# ps, zs = np.array([0,-4,-6]), np.array([zc])
# zeros = np.sum(np.angle(p-zs))
# poles = np.sum(np.angle(p-ps))
# pc = np.real(p) - np.imag(p)/np.tan(np.pi*np.arange(-5,5)-(poles-zeros))