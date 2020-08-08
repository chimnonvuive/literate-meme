from sympy.logic.boolalg import And
from numpy import pi
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import sympy.physics.mechanics as me

np.set_printoptions(precision=2, suppress=True)
me.init_vprinting()

t = sp.symbols('theta') # radian

def Full345(tn, h, beta):
    beta = beta * pi / 180
    s = h * (10*(t/beta)**3 - 15*(t/beta)**4 + 6*(t/beta)**5)
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def Full4567(tn, h, beta):
    beta = beta * pi / 180
    s = h * (35*(t/beta)**4 - 84*(t/beta)**5 + 70*(t/beta)**6 - 20*(t/beta)**7)
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def Full8thOrder(tn, h, beta):
    beta = beta * pi / 180
    s = h * (6.09755*(t/beta)**3 - 20.78040*(t/beta)**5 + 26.73155*(t/beta)**6 -\
          13.60965*(t/beta)**7 + 2.56095*(t/beta)**8)
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def FullCycloidal(tn, h, beta):
    beta = beta * pi / 180
    s = h * (t/beta - 1/2/pi*sp.sin(2*pi*t/beta))
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def FullDoubleHarmonic(tn, h, beta):
    beta = beta * pi / 180
    s = h/2* ((1-sp.cos(pi*t/beta)) - (1-sp.cos(2*pi*t/beta))/4)
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def FullModifiedTrapezoidal(tn, h, beta):
    beta = beta * pi / 180
    s = sp.Piecewise(
        (h * (.38898448*t/beta - .0309544*sp.sin(4*pi*t/beta)), 
         And(0<=t, t<beta/8)),
        (h * (2.44406184*(t/beta)**2 - .22203097*t/beta + .00723407), 
         And(beta/8<=t, t<3*beta/8)),
        (h * (1.6110154*t/beta - .0309544*sp.sin(4*pi*t/beta-pi) - .3055077),
         And(3*beta/8<=t, t<5*beta/8)), 
        (h * (-2.44406184*(t/beta)**2 + 4.6660917*t/beta - 1.2292648),
         And(5*beta/8<=t, t<7*beta/8)), 
        (h * (.6110154 + .38898448*t/beta + .0309544*sp.sin(4*pi*t/beta - 3*pi)),
         And(7*beta/8<=t, t<=beta))
    )
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)
    
def FullSimpleHarmonic(tn, h, beta):
    beta = beta * pi / 180
    s = h/2 * (1 - sp.cos(pi*t/beta))
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def FullModifiedSinusoidal(tn, h, beta):
    beta = beta * pi / 180
    s = sp.Piecewise(
        (h * (.43990085*t/beta - .0350062*sp.sin(4*pi*t/beta)), 
        And(0<=t, t<beta/8)),
        (h * (.28004957 + .43990085*t/beta - .31505577*sp.cos(4*pi/3*t/beta-pi/6)), 
        And(beta/8<=t, t<7*beta/8)),
        (h * (.56009915 + .43990085*t/beta - .0350062*sp.sin(4*pi*t/beta-2*pi)), 
        And(7*beta/8<=t, t<=beta))
    )
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def FullSCCA(tn, h, beta, a, b, c):
    if a+b+c != 0:
        raise ValueError('a + b + c does not equal 0!')
    beta = beta * pi / 180
    x = t / beta
    Ca = 4*pi**2 / ((pi**2-8)*(a**2-c**2) - 2*pi*(pi-2)*a + pi**2)
    Cv = Ca * ((a+c)/pi + b/2)
    Cj = Ca * pi / a
    s = h*sp.Piecewise(
        (Ca * (a/pi*x - (a/pi)**2*sp.sin(pi/a*x)),
         And(0<=x, x<=a/2, a!=0)),
        (Ca * (x**2/2 + a*(1/pi-1/2)*x + a**2*(1/8-1/pi**2)),
         And(a/2<=x, x<=(1-c)/2)),
        (Ca * ((a/pi+b/2)*x + (c/pi)**2 + a**2*(1/8-1/pi**2) - (1-c)**2/8 - (c/pi)**2*sp.cos(pi/c*(x-(1-c)/2))),
         And((1-c)/2<=x, x<=(1+c)/2, c!=0)),
        (Ca * (-x**2/2 + (a/pi+1-a/2)*x + (2*c**2-a**2)*(1/pi**2-1/8) - 1/4),
         And((1+c)/2<=x, x<=1-a/2)),
        (Ca * (a/pi*x + 2*(c**2-a**2)/pi**2 + ((1-a)**2-c**2)/4 - (a/pi)**2*sp.sin(pi/a*(x-1))),
         And(1-a/2<=x, x<=1, a!=0)),
        (1,
         And(x>1))
    )
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)
    
def HalfSimpleHarmonic1(tn, h, beta):
    beta = beta * pi / 180
    s = h * (1 - sp.cos(pi*t/2/beta))
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)

def HalfSimpleHarmonic2(tn, h, beta):
    beta = beta * pi / 180
    s = h * sp.sin(pi*t/2/beta)
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)
    
def HalfCycloidal1(tn, h, beta):
    beta = beta * pi / 180
    s = h * (t/beta + 1/pi*sp.sin(pi*t/beta))
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    s = sp.lambdify(t, s)
    return s(tn), v(tn), a(tn), j(tn)
    
def HalfCycloidal2(h, beta):
    beta = beta * pi / 180
    s = h * (t/beta - 1/pi*sp.sin(pi*t/beta))
    v = sp.lambdify(t, s.diff())
    a = sp.lambdify(t, s.diff(t,2))
    j = sp.lambdify(t, s.diff(t,3))    
    return sp.lambdify(t, s), v, a, j

def PlotSVAJ(xs, ys, dys, ddys, dddys, omg, ticks=4, savefig=False):
    
    fig, ax = plt.subplots(2,2,figsize=(12,8))
    
    ax[0][0].plot(xs*180/pi,ys,label='S',c='k')
    ax[0][0].set_xticks(np.linspace(0, 360, ticks))
    ax[0][0].set_xticklabels(np.linspace(0, 360, ticks))
    ax[0][0].set_xlabel(r'$\theta$, deg')
    ax[0][0].set_ylabel(r'displacement S, mm')
    ax[0][0].grid()
    
    ax[0][1].plot(xs*180/pi,dys*omg,label='V',c='k')
    ax[0][1].set_xticks(np.linspace(0, 360, ticks))
    ax[0][1].set_xticklabels(np.linspace(0, 360, ticks))
    ax[0][1].set_xlabel(r'$\theta$, deg')
    ax[0][1].set_ylabel(r'velocity V, mm/s')
    ax[0][1].grid()
    
    ax[1][0].plot(xs*180/pi,ddys*omg**2/1e3,label='A',c='k')
    ax[1][0].set_xticks(np.linspace(0, 360, ticks))
    ax[1][0].set_xticklabels(np.linspace(0, 360, ticks))
    ax[1][0].set_xlabel(r'$\theta$, deg')
    ax[1][0].set_ylabel(r'acceleration A, m/s$^2$')
    ax[1][0].grid()
    
    ax[1][1].plot(xs*180/pi,dddys*omg**3/1e3,label='J',c='k')
    ax[1][1].set_xticks(np.linspace(0, 360, ticks))
    ax[1][1].set_xticklabels(np.linspace(0, 360, ticks))
    ax[1][1].set_xlabel(r'$\theta$, deg')
    ax[1][1].set_ylabel(r'jerk J, m/s$^3$')
    ax[1][1].grid()
    
    plt.tight_layout()
    print(r'Peak displacement: {0:.3f} (mm)'.format(np.max(np.abs(ys))))
    print(r'Peak velocity: {0:.3f} (mm/s)'.format(np.max(abs(dys*omg))))
    print(r'Peak acceleration: {0:.3f} (m/s^2)'.format(np.max(abs(ddys*omg**2)/1e3)))
    print(r'Peak jerk: {0:.3f} (m/s^3)'.format(np.max(abs(dddys*omg**3/1e3))))
    
    if savefig is True:
        plt.savefig('SVAJ.png',dpi=600)