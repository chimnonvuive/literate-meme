import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

def pos(theta1, AB, AD, BC, CD):
    
    rB = kA.toExp(AB, theta1)
    rDB = AD - rB
    
    # find angle of rCB, rDC
    res1, res2 = kA.CPA4(BC, CD, rDB)
    CB_angle1 = res1[0]
    CB_angle2 = res2[0]
    rC1 = kA.toExp(BC, CB_angle1) + rB
    rC2 = kA.toExp(BC, CB_angle2) + rB
    return (rB, rC1), (rB, rC2)

def vel(theta1, omg1, AB, AD, BC, CD, variation=0):
    rB, rC = pos(theta1, AB, AD, BC, CD)[variation]
    vB = kA.vel(rB, omg1)
    omg2, omg3 = kA.vCPA4(vB, rC-rB, rC-AD)
    vC = kA.vel(rC-AD, omg3)
    return vB, vC, omg2, omg3

def acc(theta1, omg1, alp1, AB, AD, BC, CD, variation=0):
    rB, rC = pos(theta1, AB, AD, BC, CD)[variation]
    omg2, omg3 = vel(theta1, omg1, AB, AD, BC, CD, variation)[2:]
    aB = kA.acc(rB, omg1, alp1)
    alp2, alp3 = kA.aCPA4(aB, rC-rB, rC-AD, omg2, omg3)
    aC = kA.acc(rC-AD, omg3, alp3)
    return aB, aC, alp2, alp3

def reactions(theta1, omg1, alp1, ms, Is, M_ext,
              AB, AD, BC, CD, g=-9.807j, variation=0):
    rB, rC = pos(theta1, AB, AD, BC, CD)[variation]
    aB, aC, alp2, alp3 = acc(theta1, omg1, alp1, AB, AD, BC, CD, variation)
    aGs = np.array([aB/2, (aB + aC)/2, aC/2])
    
    M_ext *= -np.sign(vel(theta1, omg1, AB, AD, BC, CD, variation)[3])
    Gs, FGs = ms*g, ms*aGs
    
    # find F03, F12 using dyad method (link 2-3)
    # equilibrium equations, 
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                dA.moment(AD-rC), [0, 0],
                [0, 0], dA.moment(rB-rC)], (4, 4))
    b = dA.mat([kA.toList(FGs[1]+FGs[2]-Gs[1]-Gs[2]),
               [Is[2]*alp3 + dA.toCross((AD-rC)/2, FGs[2]-Gs[2]) - M_ext,
                Is[1]*alp2 + dA.toCross((rB-rC)/2, FGs[1]-Gs[1])]], 4)
    F03, F12 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F23, F01
    F23 = FGs[2] - Gs[2] - F03
    F01 = FGs[0] - Gs[0] - (-F12)
    
    # find drive moment
    M_eq = Is[0]*alp1 - (dA.toCross(rB/2, -F12) + dA.toCross(-rB/2, F01))
    return F03, F23, F12, F01, M_eq

def animate(theta1s, AB, AD, BC, CD, variation=0):
    
    rA, rD = 0, AD
    fig, ax = plt.subplots()
    container = []
    
    for theta1 in theta1s:
        rB, rC = pos(theta1, AB, AD, BC, CD)[variation]
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rB, rC]), np.imag([rB, rC]), 'go-')
        link3, = ax.plot(np.real([rD, rC]), np.imag([rD, rC]), 'bo-')
        A = ax.text(np.real(rA), np.imag(rA)+.01, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.01, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.01, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.01, r'$D$')
        container.append([link1, link2, link3, A, B, C, D])


    from matplotlib.animation import ArtistAnimation
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_ylim(ymax=CD+.05)
    ax.grid()
    ax.set_aspect('equal')
    ani.save('mechanism_1_{0:d}.mp4'.format(variation), 
              fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(omg1, alp1, ms, Is, M_ext,
              AB, AD, BC, CD, g=-9.807j, variation=0, frames=1000):
    
    # plot velocity and acceleration of point G2
    theta1s = np.linspace(0, 2*np.pi, frames)
    vCx, aCx, Meq = [], [], []
    vCy, aCy = [], []
    for theta1 in theta1s:
        vCx.append(np.real(vel(theta1, omg1, AB, AD, BC, CD, variation)[1]))
        vCy.append(np.imag(vel(theta1, omg1, AB, AD, BC, CD, variation)[1]))
        # aCx.append(np.real(pos(theta1, AB, AD, BC, CD)[variation][1]))
        # aCy.append(np.imag(pos(theta1, AB, AD, BC, CD)[variation][1]))
        aCx.append(np.real(acc(theta1, omg1, alp1, AB, AD, BC, CD, variation)[1]))
        aCy.append(np.imag(acc(theta1, omg1, alp1, AB, AD, BC, CD, variation)[1]))
        Meq.append(reactions(theta1, omg1, alp1, ms, Is, M_ext, AB, AD, BC, CD, g, variation)[-1])
    
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(theta1s, vCx, label='x-axis')
    ax[0].plot(theta1s, vCy, label='y-axis')
    ax[1].plot(theta1s, aCx, label='x-axis')
    ax[1].plot(theta1s, aCy, label='y-axis')
    ax[2].plot(theta1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of C($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of C ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('mechanism_1_{0:d}'.format(variation), dpi=600)

AB, AD, BC, CD, theta1, n = .08, .19, .21, .12,  60, 500
variation=0

ds = np.array([1e-3, 1e-3, 1e-3])
ws = np.array([1e-2, 1e-2, 1e-2])
hs = np.array([AB, BC, CD])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

g, M_ext = -9.807j, 3000
theta1, omg1 = theta1*np.pi/180, n*np.pi/30

# animate(np.linspace(0,2*np.pi,200), AB, AD, BC, CD, variation)
plot_vel_acc_M(omg1, 30, ms, Is, M_ext, AB, AD, BC, CD, variation)