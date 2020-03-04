import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(theta1, case):
    rB = kA.toExp(case[0], theta1)
    rE = -case[1]
    rC = rB + kA.toExp(case[2], np.angle(rB-rE))
    rC2 = rB + kA.toExp(case[0]+case[1], np.angle(rE-rB))
    CD_angle = kA.CPA2(-np.pi/2, case[3], rC-case[4])[1][1]
    rD = rC - kA.toExp(case[3], CD_angle)
    return rB, rC, rD, rE, rC2

def vel(theta1, case):
    rB, rC, rD, rE, rC2 = pos(theta1, case)
    vB = kA.vel(rB, case[-1])
    v2_T, omg2 = kA.vCPA1(0, rB-rE, rB, case[-1])
    vC, vC2 = kA.vel(rC-rE, omg2, v2_T), kA.vel(rC2-rE, omg2, v2_T)
    v5_T, omg4 = kA.vCPA2(-vC, 1j, rD-rC)
    vD = kA.vel(1j, v_T=v5_T)
    return vB, vC, vD, vC2, omg2, omg4, v2_T, v5_T

def acc(theta1, case):
    rB, rC, rD, rE, rC2 = pos(theta1, case)
    omg2, omg4, v2_T, v5_T = vel(theta1, case)[-4:]
    aB = kA.acc(rB, case[-1])
    a2_T, alp2 = kA.aCPA1(0, rB-rE, rB, omg2, case[-1], v1_T=v2_T)
    aC = kA.acc(rC-rE, omg2, alp2, v2_T, a2_T)
    aC2 = kA.acc(rC2-rE, omg2, alp2, v2_T, a2_T)
    a5_T, alp4 = kA.aCPA2(-aC, 1j, rD-rC, omg2=omg4, v1_T=v5_T)
    aD = kA.acc(1j, v_T=v5_T, a_T=a5_T)
    return aB, aC, aD, aC2, alp2, alp4

def virtual(theta1, case, ms, F_ext, g=-9.807j):
    F_ext*=-np.sign(vel(theta1, case)[1])
    vB, vC, vD, vC2 = vel(theta1, case)[:4]
    vGs = [vB/2, (vC+vC2)/2, 0, (vC+vD)/2, vD]
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F_ext, vD))/case[-1]
    return M_eq

def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    theta1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for theta1 in theta1s:
        rB, rC, rD, rE, rC2 = pos(theta1, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rC, rC2]), np.imag([rC, rC2]), 'go-')
        link4, = ax.plot(np.real([rC, rD]), np.imag([rC, rD]), 'bo-')
        ax.plot(np.real(rD), np.imag(rD), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.02, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.02, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.02, r'$C_1$')
        C2 = ax.text(np.real(rC2), np.imag(rC2)+.02, r'$C_2$')
        D = ax.text(np.real(rD), np.imag(rD)+.02, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.02, r'$E$')
        
        container.append([link1, link2, link4, A, B, C2, C, D, E])
    ax.plot(np.real(rE), np.imag(rE), 'ro-')
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    plt.show()
    ani.save('mechanism_10.mp4', fps=60, bitrate=5000, dpi=400)
    
def plot_vel_acc_M(case, ms, Is, F_ext=2000j, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    theta1s = np.linspace(0, 2*np.pi, frames)
    vD, aD, Meq = [], [], []
    
    for theta1 in theta1s:
        
        vD.append(np.imag(vel(theta1, case)[2]))
        aD.append(np.imag(acc(theta1, case)[2]))
        Meq.append(virtual(theta1, case, ms, F_ext, g))
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(theta1s, vD, label='y-axis')
    ax[1].plot(theta1s, aD, label='y-axis')
    ax[2].plot(theta1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of joint D ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of joint D ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('mechanism_10', dpi=600)

#         AB   AE   BC   CD   a   theta1  n
case1 = [.10, .25, .17, .3, .09,  120,    600*np.pi/30]

ds = 1e-3
ws = np.array([1e-2, 1e-2, 5e-2, 1e-2, 5e-2])
hs = np.array([.10, (.10+.25+.17), 2e-2, .3, 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

animate(case1)
plot_vel_acc_M(case1, ms, Is)