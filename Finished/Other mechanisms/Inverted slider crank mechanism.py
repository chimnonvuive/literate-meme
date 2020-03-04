import numpy as np
import matplotlib.pyplot as plt
import kinematicAnalysis as kA
from sympy import *
from sympy.abc import x, y, z
from symbolicAnalysis import *

def pos(phi, AC=.15, BC=.2, AD=0.35):
    
    rC = AC
    BA1, C_angle1 = kA.CPA2(phi, BC, rC)[0]
    BA2, C_angle2 = kA.CPA2(phi, BC, rC)[1]
    rB1, rB2 = kA.toExp(BA1, phi), kA.toExp(BA2, phi)
    rBC1, rBC2 = -kA.toExp(BC, C_angle1), -kA.toExp(BC, C_angle2)
    
    # BA must be positive and the angle ABC has to be smaller than 90 deg
    if np.abs(np.angle(rBC1)-np.angle(rB1))<np.pi/2 and BA1>0:
        rB, rBC = rB1, rBC1
    else: rB, rBC = rB2, rBC2
    rD = kA.toExp(AD, np.angle(rB))
    return rB, rD, rBC
    
    # rC, rB, rD = vec([AC, 0, 0]), x*vec([cos(phi), sin(phi), 0]), AD*vec([cos(phi), sin(phi), 0])
    # sol = solve((rB-rC)[0]**2 + (rB-rC)[1]**2 - BC**2, x)
    # if sol[0]>0: rB=rB.subs(x, sol[0])
    # else: rB=rB.subs(x, sol[1])
    # return rB, rD
print(pos(np.pi))
def vel(phi, omg1, AC=0.15, BC=0.2, AD=0.35): # AC < BC
    
    "find rB, rC, rD"
    rA, rC = 0, AC
    rB, rD, rBC = pos(phi, AC, BC, AD)
    
    "find velocities"
    # find vD
    vD, vB1 = kA.toVel(rD, omg1), kA.toVel(rB, omg1)
    
    # find vB2=vB3, omg2=omg3
    # vB2 =     vB1 +       vB2B1
    # _|_BC     _|_AB       //AB
    # omg2?     known       |VB2B1|?
    vB2_coeff = kA.toExp(rB-rC, np.pi/2)
    # using vB2 = vB1 + vB2B1 and the fact that vB2B1 only has 2 directions:
    # one moving along rBC and the other is along rCB, we get the following
    # code:
    omg2, vB2B1 = kA.solve_w_angle(vB1, vB2_coeff, np.angle(rB))
    vB2B1 = kA.toExp(vB2B1, np.angle(rB))
    vB2=vB3 = omg2*vB2_coeff
    return vD, vB1, vB2, vB2B1, omg2
    
print(vel(pi/3, 1))
def acc(phi1, omg1, alp1, AC=0.15, BC=0.2, AD=0.35):
    "find accelerations"
    # find aD
    rA, rC = 0, AC
    rB, rD, rBC = pos(phi1, AC, BC, AD)
    vB2B1, omg2 = vel(phi1, omg1, AC, BC, AD)[-2:]
    aD = -omg1**2*rD
    
    # find aB1, aB2=aB3, alp2=alp3
    # aB2n +    aB2t =      aB1 +   aB2B1c +    aB2B1r
    # //BC      _|_BC       known   _|_vB2B1    //AB
    # known     alp2?       known   known       |aB2B1r|?
    aB1 = kA.totAcc(rB, omg1, alp1)
    aB2n = -omg2**2*(rB-rC)
    aB2B1c = 2*omg2*kA.toExp(vB2B1, np.pi/2)
    # b = [np.real(aB1 + aB2B1c - aB2n), np.imag(aB1 + aB2B1c - aB2n)]
    # print(aB2B1c)
    aB2t_coeff = kA.toExp(rB-rC, np.pi/2)
    # aB2B1r_angle = phi1 + np.pi/2
    # # print(aB2B1r_coeff)
    # for i in range(2):
    #     A = [[np.real(aB2t_coeff), -np.cos(aB2B1r_angle)[i]],
    #          [np.imag(aB2t_coeff), -np.sin(aB2B1r_angle)[i]]]
    #     alp2, norm_aB2B1r = np.linalg.solve(A, b)
    #     # aB2_check = aB1 + aB2B1c + np.exp(1j*aB2B1r_angle[i])*norm_aB2B1r
    #     # print(aB2, aB2_check)
    #     # if aB2 == aB2_check: break
    #     if np.abs(np.linalg.solve(A, b)[0])>10e-4:
    #         alp2, norm_aB2B1r = np.linalg.solve(A, b)
    #     #     print(norm_aB2B1r)
    #         break
    alp2, aB2B1r = kA.solve_w_angle(aB1+aB2B1c-aB2n, aB2t_coeff, np.angle(rB))
    aB2B1r = kA.toExp(aB2B1r, np.angle(rB))
    aB2 = aB2n + alp2*aB2t_coeff
    # print(aB2)
    return aB1, aD, aB2, alp2


def animate(phi1s=np.linspace(0, 360, 300)*np.pi/180, AC=.15, BC=.2, AD=0.35):
    
    from matplotlib.animation import ArtistAnimation
    fig, ax = plt.subplots()
    container = []
    rA, rC = [0,0,0], [AC,0,0]
    
    for phi in phi1s:
        rB, rD = pos(phi, AC, BC, AD)
        link1, = ax.plot([rA[0], rD[0]], [rA[1], rD[1]], 'ro-')
        link3, = ax.plot([rB[0], rC[0]], [rB[1], rC[1]], 'bo-')
        # A = ax.text(np.real(rA), np.imag(rA)+.01, 'A')
        # B = ax.text(np.real(rB)-.01, np.imag(rB)+.01, 'B')
        # C = ax.text(np.real(rC), np.imag(rC)+.01, 'C')
        # D = ax.text(np.real(rD), np.imag(rD)+.01, 'D')
        container.append([link1, link3])#, A, B, C, D])
        
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    # ani.save('Inverted_R_RTR.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel(omg1, phi1s=np.linspace(0, 360, 300)*np.pi/180,
              AC=.15, BC=.2, AD=.35):
    
    """
    Plot velocity graph of slider 2
    """
    
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    v2x = []
    for phi in phi1s:
        v2x.append(np.real(vel(phi, omg1, AC=.15, BC=.2, AD=.35)[2]))
        
    fig, ax = plt.subplots()
    ax.plot(phi1s, v2x)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{AB}$ ($radian$)')
    ax.set_ylabel(r'Velocity of slider 2 ($m/s$)')
    ax.grid()
    
    fig.savefig('Velocity_diagram_Inverted_R_RTR', dpi=600, papertype='a4')
    
def plot_acc(omg1, alp1, phi1s=np.linspace(0, 360, 300)*np.pi/180,
              AC=.15, BC=.2, AD=0.35):
    
    """
    Plot acceleration graph of slider 2
    """
    
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    aCxs = []
    for phi in phi1s:
        aCxs.append(np.real(acc(phi, omg1, alp1, AC=.15, BC=.2, AD=.35)[2]))
    
    fig, ax = plt.subplots()
    ax.plot(phi1s, aCxs)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{AB}$ ($radian$)')
    ax.set_ylabel(r'Acceleration of slider 2 ($m/s^2$)')
    ax.grid()
    
    fig.savefig('Acceleration_diagram_Inverted_R_RTR', dpi=600, papertype='a4')    