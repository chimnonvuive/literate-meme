from control import *
from numpy import *
from matplotlib.pyplot import *

f = tf(poly([-8]),poly([-3,-6,-10]))
z = 0.4559
root_locus(f)
x = linspace(-10,10,1000)
y = tan(pi-arccos(z))*x
plot(x,y)

def FindP(s,ps,zs,k):
    from sympy import solve, arg, symbols
    pc = symbols('p', real=True)
    eq = 0
    for z in zs: eq += angle(z-s)
    for p in ps: eq -= angle(p-s)
    return list(solve(eq - arg(pc-s) - pi*(2*k+1),pc))

# pc = FindP(-2.0138+5.2541j,[0,-4,-6],[-5], 0)

# f = tf(poly([-5]),poly([-42.9598,0,-4,-6]))
# rlocus(f)
# x = linspace(-10,10,1000)
# y = tan(pi-arccos(z))*x
# plot(x,y)