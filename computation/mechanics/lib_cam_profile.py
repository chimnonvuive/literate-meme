from sympy.logic.boolalg import And
from sympy import symbols, sin, cos, Piecewise
from numpy import pi
from pint import *
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import sympy.physics.mechanics as me

def compact(arr, ud):
    return np.asarray(arr.to(ud)/ud, dtype=float)

np.set_printoptions(precision=2, suppress=True)
me.init_vprinting()

t, h, beta = sp.symbols('theta h beta')

# For SCCA
a, b, c = sp.symbols('a b c')
x = t / beta
Ca = 4*pi**2 / ((pi**2-8)*(a**2-c**2) - 2*pi*(pi-2)*a + pi**2)

DOUBLE_DWELL = {
        '345'  : h * (10*x**3 - 15*x**4 + 6*x**5),
        '4567' : h * (35*x**4 - 84*x**5 + 70*x**6 - 20*x**7),
        'cycloidal' : h * (x - 1/2/pi*sin(2*pi*x)),
        'modified trapezoidal' : h * Piecewise(
        (.38898448*x - .0309544*sin(4*pi*x), 
         And(0<=x, x<1/8)),
        (2.44406184*x**2 - .22203097*x + .00723407, 
         And(1/8<=x, x<3/8)),
        (1.6110154*x - .0309544*sin(4*pi*x-pi) - .3055077,
         And(3/8<=x, x<5/8)), 
        (-2.44406184*x**2 + 4.6660917*x - 1.2292648,
         And(5/8<=x, x<7/8)), 
        (.6110154 + .38898448*x + .0309544*sin(4*pi*x - 3*pi),
         And(7/8<=x, x<=1))),
        'modified sine' : h * Piecewise(
        (.43990085*x - .0350062*sin(4*pi*x), 
        And(0<=x, x<beta/8)),
        (.28004957 + .43990085*x - .31505577*cos(4*pi/3*x-pi/6), 
        And(beta/8<=x, x<7/8)),
        (.56009915 + .43990085*x - .0350062*sin(4*pi*x-2*pi), 
        And(7/8<=x, x<=1))),
        'SCCA' : Ca * h * Piecewise(
        (a/pi*x - (a/pi)**2*sin(pi/a*x),
         And(0<=x, x<=a/2, a!=0)),
        (x**2/2 + a*(1/pi-1/2)*x + a**2*(1/8-1/pi**2),
         And(a/2<=x, x<=(1-c)/2)),
        ((a/pi+b/2)*x + (c/pi)**2 + a**2*(1/8-1/pi**2) - (1-c)**2/8 \
             - (c/pi)**2*cos(pi/c*(x-(1-c)/2)),
         And((1-c)/2<=x, x<=(1+c)/2, c!=0)),
        (-x**2/2 + (a/pi+1-a/2)*x + (2*c**2-a**2)*(1/pi**2-1/8) - 1/4,
         And((1+c)/2<=x, x<=1-a/2)),
        (a/pi*x + 2*(c**2-a**2)/pi**2 + ((1-a)**2-c**2)/4 \
             - (a/pi)**2*sin(pi/a*(x-1)),
         And(1-a/2<=x, x<=1, a!=0)),
        (1,
         And(x>1))),
        '8th order'  : h * (6.09755*x**3 - 20.78040*x**5 + 26.73155*x**6 \
                            - 13.60965*x**7 + 2.56095*x**8),
        'harmonic'  : h/2 * (1 - cos(pi*x)),
        'double harmonic' : h/2* ((1-cos(pi*x)) - (1-cos(2*pi*x))/4)
        }

LIST_TYPE = list(DOUBLE_DWELL.keys())
LIST_S = list(DOUBLE_DWELL.values())

def RDFD(tn, hn, betan, curve, *args, rise=True):
    
    u = tn._REGISTRY
    betan = float(betan.to(u.rad) / u.rad)
    tn = np.asarray(tn.to(u.rad) / u.rad, dtype=float)
    h_unit = hn.units
    hn = float(hn / h_unit)
    try:
        S = DOUBLE_DWELL.get(curve).subs({h: hn, beta: betan})
        V = sp.lambdify(t, S.diff())
        A = sp.lambdify(t, S.diff(t,2))
        J = sp.lambdify(t, S.diff(t,3))    
        S = sp.lambdify(t, S)
    except ValueError:
        an, bn, cn = args
        if an+bn+cn != 1:
            raise ValueError('a + b + c does not equal 0!')
        Can = Ca.subs({a:an, c:cn})
        S = DOUBLE_DWELL.get(curve).subs({h: hn, beta: betan,
                                          Ca:Can, a:an, b:bn, c:cn})
        V = sp.lambdify(t, S.diff())
        A = sp.lambdify(t, S.diff(t,2))
        J = sp.lambdify(t, S.diff(t,3))    
        S = sp.lambdify(t, S)
    
    if rise is True:
        return np.array([S(tn), V(tn), A(tn), J(tn)]) * u.dimensionless * h_unit
    else:
        return -np.array([S(tn)-hn, V(tn), A(tn), J(tn)]) * u.dimensionless * h_unit

def RadialCamTranslatingFlatFaced(xs, ys, dys, ddys, n, rho_min, Rb, b, show=False, savefig=False):
    
    u = xs._REGISTRY
    s, x, theta = ys, dys, xs.to(u.rad)
    rho_min = rho_min * ys.units
    Rb = Rb * ys.units
    b = b * ys.units
    Rb_min = rho_min - min(ys + ddys)
    b_min = max(dys)-min(dys)
    
    if Rb <= Rb_min:
        print('Rb < {0:.4f}. Choose base circle radius = {0:.4f} \
                                                  instead'.format(Rb_min))
        Rb = Rb_min
    if b.magnitude <= b_min.magnitude:
        print('b < {0:.4f}. Choose face width = {0:.4f} instead'.format(
                                                b_min.magnitude * ys.units))
        b = b_min.magnitude * ys.units
    
    RA = np.sqrt((Rb+s)**2 + x**2)
    sigA = theta + np.arctan(x/(Rb + s))
    
    xn = RA * np.cos(sigA)
    yn = -np.sign(n.magnitude) * RA * np.sin(sigA)
    
    if show is True:
        fig, ax = plt.subplots()
        ax.plot(xn, yn, 'k', label='cam profile')
        ax.plot(0, 0, 'ro')
        ax.set_xlabel('abscissa (mm)'), ax.set_ylabel('ordinal (mm)')
        ax.grid(), ax.set_aspect('equal'), ax.legend()
        plt.tight_layout()
        plt.savefig('cam profile.png', dpi=600)
        
    if savefig is True:
        np.savetxt('cam_profile.txt',np.array([xn.magnitude, yn.magnitude,
                                   np.zeros(len(xn))]).T,delimiter='\t')

def RadialCamTranslatingRoller(xs, ys, n, Rp, e, show=False, savefig=False):
    
    u = xs._REGISTRY
    e, s, theta, Rp = e * ys.units, ys, xs.to(u.rad), Rp * ys.units
    d = np.sqrt(Rp**2 - e**2)
    R = np.sqrt((d+s)**2 + e**2)
    Lambda = np.arctan(e/(d+s)) + theta
    
    xn = R * np.cos(Lambda)
    yn = -np.sign(n.magnitude) * R * np.sin(Lambda)
    
    if show is True:
        fig, ax = plt.subplots()
        ax.plot(xn, yn, 'k', label='cam profile')
        ax.set_xlabel('abscissa (mm)')
        ax.set_ylabel('ordinal (mm)')
        ax.grid()
        ax.set_aspect('equal')
        ax.legend()
        plt.tight_layout()
        plt.savefig('cam profile.png', dpi=600)
        
    if savefig is True:
        np.savetxt('cam_profile.txt',np.array([xn.magnitude, yn.magnitude,
                                   np.zeros(len(xn))]).T,delimiter='\t')


def PlotSVAJ(xs, ys, dys, ddys, dddys, omg, ticks=5, savefig=False):
    
    u = xs._REGISTRY
    omg = omg.to(u.rad/u.s)
    dys = dys*omg
    ddys = ddys*omg**2
    dddys = dddys*omg**3
    
    xn = compact(xs, u.deg)
    yn = compact(ys, u.mm)
    dyn = compact(dys, u.mm/u.s)
    ddyn = compact(ddys, u.mm/u.s**2)
    dddyn = compact(dddys, u.mm/u.s**3)
    
    # fig, ax = plt.subplots(2,2,figsize=(12,8))
    
    # ax[0][0].plot(xn,yn,label='S',c='k')
    # ax[0][0].set_xticks(np.linspace(0, 360, ticks))
    # ax[0][0].set_xticklabels(np.linspace(0, 360, ticks))
    # ax[0][0].set_xlabel(r'$\theta$, deg')
    # ax[0][0].set_ylabel(r'displacement S, mm')
    # ax[0][0].grid()
    
    # ax[0][1].plot(xn,dyn,label='V',c='k')
    # ax[0][1].set_xticks(np.linspace(0, 360, ticks))
    # ax[0][1].set_xticklabels(np.linspace(0, 360, ticks))
    # ax[0][1].set_xlabel(r'$\theta$, deg')
    # ax[0][1].set_ylabel(r'velocity V, mm/s')
    # ax[0][1].grid()
    
    # ax[1][0].plot(xn,ddyn,label='A',c='k')
    # ax[1][0].set_xticks(np.linspace(0, 360, ticks))
    # ax[1][0].set_xticklabels(np.linspace(0, 360, ticks))
    # ax[1][0].set_xlabel(r'$\theta$, deg')
    # ax[1][0].set_ylabel(r'acceleration A, mm/s$^2$')
    # ax[1][0].grid()
    
    # ax[1][1].plot(xn,dddyn,label='J',c='k')
    # ax[1][1].set_xticks(np.linspace(0, 360, ticks))
    # ax[1][1].set_xticklabels(np.linspace(0, 360, ticks))
    # ax[1][1].set_xlabel(r'$\theta$, deg')
    # ax[1][1].set_ylabel(r'jerk J, mm/s$^3$')
    # ax[1][1].grid()

    # plt.tight_layout()
    
    # Tc = 493.66e-3*ddyn*dyn/omg.magnitude
    print(500*max(yn))
    Tc = 500*yn*dyn/omg.magnitude
    fig, ax = plt.subplots()
    ax.plot(xn,Tc*1e-3,label='Torque',c='k')
    ax.set_title('Torque on camshaft, $n=250rpm$')
    ax.set_xticks(np.linspace(0, 360, ticks))
    ax.set_xticklabels(np.linspace(0, 360, ticks))
    ax.set_xlabel(r'$\theta$, deg')
    ax.set_ylabel(r'$T_c$, $N\cdot m$')
    ax.legend()
    ax.grid()
    # print(max(abs(Tc)))
    # plt.plot(-50*yn/dyn)
    # plt.plot(-ddyn/yn)
    
    # fig, ax = plt.subplots()
    # a = ys.magnitude+ddys.magnitude
    # ax.plot(xn,a,label='A',c='k')
    # ax.plot(xn[np.where(a==min(a))],min(a),'ro')
    # ax.set_xticks(np.linspace(0, 360, ticks))
    # ax.set_xticklabels(np.linspace(0, 360, ticks))
    # ax.set_xlabel(r'$\theta$, deg')
    # ax.set_ylabel(r'$s+\frac{ds^2}{d\phi}$, mm')
    # ax.grid()
    
    plt.tight_layout()
    
    print(r'Peak displacement: {0:.4f}'.format(max(abs(ys)).to_compact()))
    print(r'Peak velocity: {0:.4f}'.format(max(abs(dys)).to_compact()))
    print(r'Peak acceleration: {0:.4f}'.format(max(abs(ddys)).to_compact()))
    print(r'Peak jerk: {0:.4f}'.format(max(abs(dddys)).to_compact()))
    
    if savefig is True:
        plt.savefig('SVAJ.png',dpi=600)