import numpy as np
import matplotlib.pyplot as plt
import kinematicAnalysis as kA
import dynamicAnalysis as dA

from matplotlib import rc


def pos(phi1, AB=1, BC=1):
        
    """
    Calculate displacements
    """
    
    rB = kA.toExp(AB, phi1)
    B_angle1 = kA.CPA2(0, BC, rB)[0][1]
    B_angle2 = kA.CPA2(0, BC, rB)[1][1]
    rBC1 = kA.toExp(BC, B_angle1)
    rBC2 = kA.toExp(BC, B_angle2)
    rC1 = rB - rBC1
    rC2 = rB - rBC2
    
    if AB==BC:
        if np.real(rB) > 0:
            return rB, rC2, rBC2
        if np.real(rB) < 0:
            return rB, rC1, rBC1
    
    return (rB, rC1, rBC1), (rB, rC2, rBC2)

def vel(phi1, omg1, variation=0, AB=1, BC=1):    
    
    """
    Calculate velocities and angular velocities
    """
    
    rB, rC, rBC = pos(phi1, AB, BC)[variation]
    vB = kA.vel(rB, omg1, 0)
    vC, omg2 = kA.vCPA2(vB, -rC, rBC, 0, 0)
    return vB, vC, omg2

def acc(phi1, omg1, alp1, variation=0, AB=1, BC=1):
    
    """
    Calculate accelerations and angular acceleration
    """
    
    rB, rC, rBC = pos(phi1, AB, BC)[variation]
    vC, omg2 = vel(phi1, omg1, variation, AB, BC)[-2:]
    aB = kA.acc(rB, omg1, 0, 0, alp1)
    aC, alp2 = kA.aCPA2(aB, -rC, rBC, 0, omg2, vC, 0, 0, 0)
    return aB, aC, alp2

def reactions(phi1, omg1, alp1, ms, Is, F_ext,
            g=-10j, variation=0, AB=1, BC=1):
    rB, rC, rBC = pos(phi1, AB, BC)[variation]
    aB, aC, alp2 = acc(phi1, omg1, alp1, variation, AB, BC)
    a_Cs = np.array([aB/2, (aB+aC)/2, aC])
    vC = vel(phi1, omg1, variation, AB, BC)[1]
    coeffs = -ms*g + ms*a_Cs
    
    coeffs[2] -= np.sign(-vC)*F_ext
    A = [[1, 0,0,0,0,0,0,0],
          [-1,0,1,0,0,0,0,0],
          [0,-1,0,1,0,0,0,0],
          [0,0,-1,0,1,0,0,0],
          [0,0,0,-1,0,1,0,0],
          [0,0,0,0,-1,0,1,0],
          [0,0,0,0,0,-1,0,1],
          [0,0,-np.imag(rBC)/2, np.real(rBC)/2, -np.imag(rBC)/2, np.real(rBC)/2,0,0]]
    b = [0, 
          np.real(-np.sign(-vC)*F_ext+coeffs[2]),
          np.imag(-np.sign(-vC)*F_ext+coeffs[2]),
          np.real(coeffs[1]),
          np.imag(coeffs[1]),
          np.real(coeffs[0]),
          np.imag(coeffs[0]),
          Is[1]*alp2]
    res = np.linalg.solve(A, b)
        
    F30 = res[0] + res[1]*1j
    F23 = res[2] + res[3]*1j
    F12 = res[4] + res[5]*1j
    F01 = res[6] + res[7]*1j
    M_equilibrium = Is[0]*alp1 - dA.toCross(rB/2, -F12) - dA.toCross(-rB/2, F01)
    return F30, F23, F12, F01, M_equilibrium

def dyad(phi1, omg1, alp1, ms, Is, F_ext,
         g=-10j, variation=0, AB=1, BC=1):
    rB, rC, rBC = pos(phi1, AB, BC)[variation]
    aB, aC, alp2 = acc(phi1, omg1, alp1, variation, AB, BC)
    a_Cs = np.array([aB/2, (aB+aC)/2, aC])
    Gs = ms*g

    vC = vel(phi1, omg1, variation, AB, BC)[1]
    F_ext *= -np.sign(vC)
    
    F12x = np.real((ms*a_Cs)[2] + (ms*a_Cs)[1]) - F_ext
    F12y = (Is[1]*alp2 + dA.toCross(rBC/2, (ms*a_Cs)[1]) + np.imag(rBC)*F12x -
            dA.toCross(rBC/2, Gs[1]))/np.real(rBC)
    F12 = F12x + F12y*1j
    F03 = (ms*a_Cs)[2] + (ms*a_Cs)[1] - (Gs[1] + Gs[2] + F12 + F_ext)
    F23 = (ms*a_Cs)[1] - (Gs[1] + F12)
    F01 = (ms*a_Cs)[0] - (Gs[0] + -F12)
    M_equilibrium = Is[0]*alp1 - dA.toCross(rB/2, -F12) - dA.toCross(-rB/2, F01)
    return F03, F23, F12, F01, M_equilibrium

def virtual(phi1, omg1, alp1, ms, F_ext, g=-10j, variation=0, AB=1, BC=1):
    rB, rC, rBC = pos(phi1, AB, BC)[variation]

def animate(phi1s, variation=0, fps=60, AB=1, BC=1):
    
    """
    Animate movement of the mechanism
    """
    
    rA = 0
    fig, ax = plt.subplots()
    container = []
    
    for phi in phi1s:
        rB, rC, rBC = pos(phi, AB, BC)[variation]
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rB, rC]), np.imag([rB, rC]), 'go-')
        A = ax.text(np.real(rA), np.imag(rA)+.07, 'A')
        B = ax.text(np.real(rB), np.imag(rB)+.07, 'B')
        C = ax.text(np.real(rC), np.imag(rC)+.07, 'C')
        container.append([link1, link2, A, B, C])
        
        ax.plot(np.real(rC), np.imag(rC), 'b,-')

    from matplotlib.animation import ArtistAnimation
    ani = ArtistAnimation(fig, container, interval=fps, blit=True)
    ax.grid()
    ax.set_aspect('equal')
    ani.save('R_RRT.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel(phi1s, omg1, variation=0, AB=1, BC=1):
    
    """
    Plot velocity graph of slider C
    """
    
    vCx = []
    for phi in phi1s:
        vCx.append(np.real(vel(phi, omg1, variation, AB, BC)[1]))
        
    fig, ax = plt.subplots()
    ax.plot(phi1s, vCx)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax.set_ylabel(r'Velocity of slider 3 ($m/s$)')
    ax.grid()
    
    fig.savefig('Velocity_diagram_R_RRT', dpi=600, papertype='a4')
    
def plot_acc(phi1s, omg1, alp1, variation=0, AB=1, BC=1):
    
    """
    Plot acceleration graph of slider C
    """
    
    aCxs = []
    for phi in phi1s:
        aCxs.append(np.real(acc(phi, omg1, alp1, variation, AB, BC)[1]))
    
    fig, ax = plt.subplots()
    ax.plot(phi1s, aCxs)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax.set_ylabel(r'Acceleration of slider 3 ($m/s^2$)')
    ax.grid()
    
    fig.savefig('Acceleration_diagram_R_RRT', dpi=600)
    
def plot_vel_acc_M(phi1s, omg1, alp1, ms, Is, F_ext,
            g=-10j, variation=0, AB=1, BC=1):
    
    # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    vCx, aCx, Meq = [], [], []
    
    for phi in phi1s:
        vCx.append(np.real(vel(phi, omg1, variation, AB, BC))[1])
        aCx.append(np.real(acc(phi, omg1, alp1, variation, AB, BC))[1])
        Meq.append(reactions(phi, omg1, alp1, ms, Is, F_ext, g, variation, AB, BC)[-1])
    fig, ax = plt.subplots(nrows=2, figsize=(6,8))
    ax[0].plot(phi1s, vCx, label='Velocity of slider 3')
    ax[0].plot(phi1s, aCx, label='Acceleration of slider 3')
    ax[1].plot(phi1s, Meq, label='Drive moment on crank 1')
    
    xlabels = np.linspace(0, 2*np.pi, 4, endpoint=False)
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax[0].set_xticks(xlabels)
    ax[0].set_xticklabels(xticks)
    ax[0].set_ylabel(r'v - a diagram of slider 3 ($m/s^2$)')
    ax[0].legend()
    ax[0].grid()
    
    ax[1].set_xticks(xlabels)
    ax[1].set_xticklabels(xticks)
    ax[1].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[1].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[1].legend()
    ax[1].grid()
    plt.tight_layout()
    
    fig.savefig('Diagram_R_RRT', dpi=600)

AB, BC = 1, 1.5
hs = np.array([1e-2, 1e-2, 1e-2])
ws = np.array([AB, BC, 1e-2])
ds = np.array([1e-3, 1e-3, 1e-3])
rhos = np.array([1e5, 1e5, 1e7])
ms = rhos*hs*ws*ds
Is = ms*(ws**2 + hs**2)/12
F_ext = 100
animate(np.linspace(0,2*np.pi,200),variation=1, BC=BC)
plot_vel_acc_M(np.linspace(0,2*np.pi,200), 1, 0, ms, Is, F_ext, variation=1, BC=BC)