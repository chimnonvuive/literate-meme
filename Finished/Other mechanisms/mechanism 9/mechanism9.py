import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(theta1, case):
    rB = kA.toExp(case[0], theta1)
    rC = case[1]
    rD2 = rC + kA.toExp(case[0]+case[1], np.angle(rB-rC))
    rD = rC + kA.toExp(case[2], np.angle(rB-rC)-np.pi/2)
    A, B_angle = kA.CPA2(0, case[3], rD)[1]
    rE = rD - kA.toExp(case[3], B_angle)
    return rB, rC, rD, rE, rD2

def vel(theta1, case):
    rB, rC, rD, rE, rD2 = pos(theta1, case)
    
    vB = kA.vel(rB, case[-1])
    v2_T, omg3 = kA.vCPA1(0, rB-rC, rB, case[-1])
    vD, vD2 = kA.vel(rD-rC, omg3), kA.vel(rD2-rC, omg3)
    vE, omg4 = kA.vCPA2(-vD, rE-rC, rE-rD)
    return vB, vD, vE, vD2, omg3, omg4, v2_T

def acc(theta1, case):
    rB, rC, rD, rE, rD2 = pos(theta1, case)
    vB, vD, vE, vD2, omg3, omg4, v2_T = vel(theta1, case)
    aB = kA.acc(rB, case[-1])
    a2_T, alp3 = kA.aCPA1(0, rB-rC, rB, omg3, case[-1], v1_T=v2_T)
    aD, aD2 = kA.acc(rD-rC, omg3, alp3), kA.acc(rD2-rC, omg3, alp3)
    aE, alp4 = kA.aCPA2(-aD, rE-rC, rE-rD, 0, omg4, v1_T=vE)
    return aB, aD, aE, aD2, alp3, alp4

def virtual(theta1, case, ms, F_ext, g=-9.807j):
    vB, vD, vE, vD2 = vel(theta1, case)[:4]
    F = F_ext*-np.sign(vE)
    vGs = np.array([vB/2, vB, (vD+vD2)/3, (vD+vE)/2, vE])
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F, vE))/case[-1]
    return M_eq

def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    theta1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for theta1 in theta1s:
        rB, rC, rD, rE, rD2 = pos(theta1, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link31, = ax.plot(np.real([rC, rD]), np.imag([rC, rD]), 'go-')
        link32, = ax.plot(np.real([rC, rD2]), np.imag([rC, rD2]), 'g-')
        
        link4, = ax.plot(np.real([rD, rE]), np.imag([rD, rE]), 'bo-')
        ax.plot(np.real(rE), np.imag(rE), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.02, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.02, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.02, r'$C$')
        D1 = ax.text(np.real(rD), np.imag(rD)+.02, r'$D_1$')
        D2 = ax.text(np.real(rD2), np.imag(rD2)+.02, r'$D_2$')
        E = ax.text(np.real(rE), np.imag(rE)+.02, r'$E$')
        
        container.append([link1, link31, link32, link4, A, B, D2, C, D1, E])
    link3, = ax.plot(np.real(rC), np.imag(rC), 'ro')
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    plt.show()
    ani.save('mechanism_9.mp4', fps=60, bitrate=5000, dpi=400)
    
def plot_vel_acc_M(case, ms, Is, F_ext=2000, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    theta1s = np.linspace(0, 2*np.pi, frames)
    v0E, a0E, Meq = [], [], []
    
    for theta1 in theta1s:
        v0E.append(vel(theta1, case)[2])
        a0E.append(acc(theta1, case)[2])
        Meq.append(virtual(theta1, case, ms, F_ext, g))
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))

    ax[0].plot(theta1s, v0E, label='x-axis')
    ax[1].plot(theta1s, a0E, label='x-axis')
    ax[2].plot(theta1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of joint E ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of joint E ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('mechanism_9', dpi=600)

#         AB  AC   CD   DE   theta1  n
case1 = [.18, .08, .2, .5,   30,     600*np.pi/30]

ds = 1e-3
ws = np.array([1e-2, 5e-2, 1e-2, 1e-2, 5e-2])
hs = np.array([.18, 2e-2, .18+.08+.2, .5, 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

F_ext = 2000

# animate(case1)
plot_vel_acc_M(case1, ms, Is, F_ext)