import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(phi1, case):
    
    AB, BC, CE, CD, EF, a, b, c = case[:8]
    rB, rD = kA.toExp(AB, phi1), a + b*1j

    # rC-rD + rB-rC = rBD
    CD_angle = kA.CPA4(CD, BC, rB-rD)[0][0]    
    rC = rD + kA.toExp(CD, CD_angle)
    rE = rC + kA.toExp(CE, np.angle(rC-rB))
    
    # find rFA', rE-rF
    res1, res2 = kA.CPA2(0, EF, rE-c*1j)
    
    # rF1 = rD + rFD1
    rF = case[7]*1j + kA.toExp(res2[0], 0)
    
    return rB, rC, rD, rE, rF

def vel(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg1 = case[-1]*np.pi/30
    vB = kA.vel(rB, omg1)
    omg2, omg3 = kA.vCPA4(vB, rC-rB, rC-rD)    
    
    # find vE = vB + vEB
    vC, vE = kA.vel(rC-rD, omg3), vB + kA.vel(rE-rB, omg2)
    v5_T, omg4 = kA.vCPA2(-vE, 1, rF-rE)
    vF = kA.vel(1, v_T=v5_T)
    return vB, vC, vE, vF, omg2, omg3, omg4, v5_T

def acc(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg2, omg3, omg4, v5_T = vel(phi1, case)[-4:]
    omg1 = case[-1]*np.pi/30
    aB = kA.acc(rB, omg1)
    alp2, alp3 = kA.aCPA4(aB, rC-rB, rC-rD, omg2, omg3)
    aC = kA.acc(rC-rD, omg3, alp3)
    aE = aB + kA.acc(rE-rB, omg2, alp2)
    a5_T, alp4 = kA.aCPA2(-aE, 1, rF-rE, omg2=omg4, v1_T=v5_T)
    aF = kA.acc(1, v_T=v5_T, a_T=a5_T)
    return aB, aC, aE, aF, alp2, alp3, alp4

def reactions(phi1, case, ms, Is, F_ext, g=-9.807j):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    aB, aC, aE, aF, alp2, alp3, alp4 = acc(phi1, case)
    F_ext *= -np.sign(np.real(vel(phi1, case)[3]))
    
    aGs = np.array([aB/2, (aB+aE)/2, aC/2, (aE+aF)/2, aF])
    Gs, FGs = ms*g, ms*aGs
    
    # find F05, F34 using dyad method
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0], dA.moment(rE-rF),
                [1, 0, 0, 0]], (4, 4))
    b = dA.mat([kA.toList(FGs[3]+FGs[4]-Gs[3]-Gs[4]-F_ext),
                [Is[3]*alp4 + dA.toCross((rE-rF)/2, FGs[3]-Gs[3]), 0]], 4)
    F05, F34 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F45
    F45 = FGs[4] - Gs[4] - F05 - F_ext
    
    # find F03, F12 using dyad method, moment at C
    rG3, rG2 = rC/2, (rB+rE)/2
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                dA.moment(rD-rC), [0, 0],
                [0, 0], dA.moment(rB-rC)], (4, 4))
    b = dA.mat([kA.toList(FGs[2]+FGs[1]-Gs[2]-Gs[1]-(-F34)),
                [Is[2]*alp3 + dA.toCross(rG3-rC, FGs[2]-Gs[2]),
                Is[1]*alp2 + dA.toCross(rG2-rC, FGs[1]-Gs[1]) \
                    - dA.toCross(rE-rC, -F34)]], 4)
    F03, F12 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F23
    F23 = FGs[2] - Gs[2] - F03
    
    # find F01
    F01 = FGs[0] - Gs[0] - (-F12)
    
    # find drive moment of crank 1
    M_eq = dA.toCross(rB/2, FGs[0]-Gs[0]) - dA.toCross(rB, -F12)
    
    return F05, F45, F23, F03, F12, F01, M_eq

def virtual_work(theta1, case, ms, F_ext, g=-9.807j):
    omg1 = case[-1]*np.pi/30
    rB, rC, rD, rE, rF = pos(theta1, case)
    vB, vC, vE, vF = vel(theta1, case)[:4]
    F_ext*=-np.sign(np.real(vF))
    
    vGs = [vB/2, (vB+vE)/2, vC/2, (vE+vF)/2, vF]
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F_ext, vF))/omg1
    return M_eq

def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    phi1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for phi in phi1s:
        rB, rC, rD, rE, rF = pos(phi, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rB, rE]), np.imag([rB, rE]), 'go-')
        link3, = ax.plot(np.real([rD, rC]), np.imag([rD, rC]), 'bo-')
        link4, = ax.plot(np.real([rE, rF]), np.imag([rE, rF]), 'ko-')
        ax.plot(np.real(rD), np.imag(rD), 'bo-')
        ax.plot(np.real(rF), np.imag(rF), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.01, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.01, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.01, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.01, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.01, r'$E$')
        F = ax.text(np.real(rF), np.imag(rF)+.01, r'$F$')
        
        container.append([link1, link2, link3, link4, A, B, C, D, E, F])
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    # ax.set_ylim(ymax=case[1]+.02)
    ax.grid()
    plt.show()
    
    ani.save('mechanism_4.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(case, ms, Is, F_ext, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    phi1s = np.linspace(0, 2*np.pi, frames)
    vF, aF, Meq = [], [], []
    
    for phi in phi1s:
        
        vF.append(np.real(vel(phi, case)[3]))
        aF.append(np.real(acc(phi, case)[3]))
        Meq.append(reactions(phi, case, ms, Is, F_ext, g)[-1])
        # Meq.append(virtual_work(phi, case, ms, F_ext, g))
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(phi1s, vF, label='x-axis')
    ax[1].plot(phi1s, aF, label='x-axis')
    ax[2].plot(phi1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram - link 5 ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram - link 5 ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($rad$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('Mechanism4', dpi=600)

#         AB   BC   CE   CD   EF   a    b    c   phi  n
case1 = [.09, .40, .25, .12, .21, .22, .35, .40, 60,  400]
case2 = [.15, .67, .45, .22, .32, .33, .60, .65, 135, 600]
case3 = [.22, 1.0, .65, .35, .60, .55, .90, 1.2, 240, 800]
case4 = [.16, .70, .50, .25, .48, .40, .60, .70, 300, 1e3]

ds = 1e-3
ws = np.array([1e-2, 1e-2, 1e-2, 1e-2, 5e-2])
hs = np.array([case1[0], case1[1]+case1[2], case1[3], case1[4], 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

F_ext = 2000

animate(case1)
plot_vel_acc_M(case1, ms, Is, F_ext)