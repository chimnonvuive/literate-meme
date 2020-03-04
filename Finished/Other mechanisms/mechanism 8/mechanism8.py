import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(theta1, case):
    rB = kA.toExp(case[0], theta1)
    rD = case[1]
    rC = rB + kA.toExp(case[2], np.angle(rD-rB))
    rB_max = rB + kA.toExp(case[2]*2, np.angle(rD-rB))
    rE = rC + kA.toExp(case[3], np.angle(rC-rB)-np.pi/2)
    rEFx = np.real(rE)-case[-3]
    A, B_angle = kA.CPA2(-np.pi/2, case[4], rEFx)[1]
    rF = rE - kA.toExp(case[4], B_angle)
    return rB, rC, rE, rD, rE, rF, rB_max

def vel(theta1, case):
    rB, rC, rE, rD, rE, rF, rB_max = pos(theta1, case)
    vB = kA.vel(rB, case[-1])
    v2_T, omg2 = kA.vCPA1(0, rB-rD, rB, case[-1])
    vB2 = kA.vel(rB_max-rD, omg2, v2_T)
    vC = kA.vel(rC-rD, omg2, v2_T)
    vE = vC + kA.vel(rE-rC, omg2)
    vF, omg4 = kA.vCPA2(-vE, rF-case[-3], rF-rE)
    return vB, vC, vE, vF, vB2, omg2, omg4, v2_T

def acc(theta1, case):
    rB, rC, rE, rD, rE, rF, rB_max = pos(theta1, case)
    omg2, omg4, v2_T = vel(theta1, case)[-3:]
    aB = kA.acc(rB, case[-1])
    a2_T, alp2 = kA.aCPA1(0, rB-rD, rB, omg2, case[-1], v1_T=v2_T)
    aB2 = kA.acc(rB_max-rD, omg2, alp2, v2_T, a2_T)
    aC = kA.acc(rC-rD, omg2, alp2, v2_T, a2_T)
    aE = aC + kA.acc(rE-rC, omg2, alp2)
    aF, alp4 = kA.aCPA2(-aE, rF-case[-3], rF-rE, 0, omg4, 0, vel(theta1, case)[3])
    return aB, aC, aE, aF, aB2, alp2, alp4

def virtual(theta1, case, ms, F_ext, g=-9.807j):
    F = F_ext*-np.sign(vel(theta1, case)[1])
    vB, vC, vE, vF, vB2, omg2, omg4, v2_T = vel(theta1, case)
    vGs = np.array([vB/2, (vB+vB2+vE)/3, 0, (vE+vF*1j)/2, vF*1j])
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F, vF))/case[-1]
    return M_eq

def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    theta1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for theta1 in theta1s:
        rB, rC, rE, rD, rE, rF, rB_max = pos(theta1, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link21, = ax.plot(np.real([rB, rB_max]), np.imag([rB, rB_max]), 'go-')
        link22, = ax.plot(np.real([rC, rE]), np.imag([rC, rE]), 'g-')
        
        link4, = ax.plot(np.real([rE, rF]), np.imag([rE, rF]), 'bo-')
        ax.plot(np.real(rF), np.imag(rF), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.02, r'$A$')
        B1 = ax.text(np.real(rB), np.imag(rB)+.02, r'$B_1$')
        B2 = ax.text(np.real(rB_max), np.imag(rB_max)+.02, r'$B_2$')
        C = ax.text(np.real(rC), np.imag(rC)+.02, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.02, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.02, r'$E$')
        F = ax.text(np.real(rF), np.imag(rF)+.02, r'$F$')
        
        container.append([link1, link21, link22, link4, A, B1, B2, C, D, E, F])
    link3, = ax.plot(np.real(rD), np.imag(rD), 'ro')
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    plt.show()
    ani.save('mechanism_8.mp4', fps=60, bitrate=5000, dpi=400)
    
def plot_vel_acc_M(case, ms, Is, F_ext=2000j, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    theta1s = np.linspace(0, 2*np.pi, frames)
    v0F, a0F, Meq = [], [], []
    
    for theta1 in theta1s:
        
        v0F.append(vel(theta1, case)[3])
        a0F.append(acc(theta1, case)[3])
        Meq.append(virtual(theta1, case, ms, F_ext, g))
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(theta1s, v0F, label='y-axis')
    ax[1].plot(theta1s, a0F, label='y-axis')
    ax[2].plot(theta1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of joint F ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of joint F ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('mechanism_8', dpi=600)

#         AB   AD   BC   CE   EF   a   theta1  n
case1 = [.22, 1.0, .65, .18, .70, .60, 300,    600*np.pi/30]

ds = 1e-3
ws = np.array([1e-2, 1e-2, 5e-2, 1e-2, 5e-2])
hs = np.array([.22, (.65*2+.18), 2e-2, .70, 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

F_ext = 2000j

# animate(case1)
plot_vel_acc_M(case1, ms, Is, F_ext)