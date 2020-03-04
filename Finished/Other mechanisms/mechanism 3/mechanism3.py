import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(phi1, case):
    
    rB = kA.toExp(case[0], phi1)
    rD = case[5] + case[6]*1j

    # rCD + rBC = rBD
    CD_angle = kA.CPA4(case[2], case[1], rB-rD)[0][0]
    rC = rD + kA.toExp(case[2], CD_angle)
    rE = rD + kA.toExp(case[3], np.pi+CD_angle)
    
    res1, res2 = kA.CPA2(0, case[4], rE-rD)
    rFD = kA.toExp(res2[0], 0)
    rF = rD + rFD
    
    return rB, rC, rD, rE, rF

def vel(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg1 = case[-1]*np.pi/30
    vB = kA.vel(rB, omg1)
    omg2, omg3 = kA.vCPA4(vB, rC-rB, rC-rD)
    vC, vE = kA.vel(rC-rD, omg3), kA.vel(rE-rD, omg3)
    v5_T, omg4 = kA.vCPA2(-vE, rF-rD, rF-rE)
    vF = kA.vel(rF-rD, v_T=v5_T)
    return vB, vC, vE, vF, omg2, omg3, omg4, v5_T

def acc(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg2, omg3, omg4, v5_T = vel(phi1, case)[-4:]
    omg1 = case[-1]*np.pi/30
    aB = kA.acc(rB, omg1)
    alp2, alp3 = kA.aCPA4(aB, rC-rB, rC-rD, omg2, omg3)
    aC, aE = kA.acc(rC-rD, omg3, alp3), kA.acc(rE-rD, omg3, alp3)
    a5_T, alp4 = kA.aCPA2(-aE, rF-rD, rF-rE, omg2=omg4, v1_T=v5_T)
    aF = kA.acc(rF-rD, v_T=v5_T, a_T=a5_T)
    return aB, aC, aE, aF, alp2, alp3, alp4

def reactions(phi1, case, ms, Is, F_ext, g=-9.807j):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    aB, aC, aE, aF, alp2, alp3, alp4 = acc(phi1, case)
    F = F_ext*-np.sign(np.real(vel(phi1, case)[1]))
    Gs, FGs = ms*g, ms*np.array([aB/2, (aB+aC)/2, (aC+aE)/2, (aE+aF)/2, aF])
    
    # find F05, F34 using dyad method
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0], dA.moment(rE-rF),
                [1, 0, 0, 0]], (4, 4))
    b = dA.mat([kA.toList(FGs[3]+FGs[4]-Gs[3]-Gs[4]-F),
                [Is[3]*alp4 + dA.toCross((rE-rF)/2, FGs[3]-Gs[3]), 0]], 4)
    F05, F34 = kA.toComplexes(np.linalg.solve(A, b))
    # find F45
    F45 = FGs[4] - Gs[4] - F05 - F
    
    # find F03, F12 using dyad method, moment at C
    rG3, rG2 = (rE+rC)/2, (rB+rC)/2
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                dA.moment(rD-rC), [0, 0],
                [0, 0], dA.moment(rB-rC)], (4, 4))
    b = dA.mat([kA.toList(FGs[2]+FGs[1]-Gs[2]-Gs[1]-(-F34)),
                [Is[2]*alp3 + dA.toCross(rG3-rC, FGs[2]-Gs[2]) - dA.toCross(rE-rC, -F34),
                 Is[1]*alp2 + dA.toCross(rG2-rC, FGs[1]-Gs[1])]], 4)
    F03, F12 = kA.toComplexes(np.linalg.solve(A, b))
    # find F23
    F23 = FGs[2] - Gs[2] - (-F34) - F03
    F01 = FGs[0] - Gs[0] - (-F12)  
    
    # find drive moment of crank 1
    M_eq = dA.toCross(rB/2, FGs[0]-Gs[0]) - dA.toCross(rB, -F12)
    return F05, F45, F23, F03, F12, F01, M_eq

def virtual_work(theta1, case, ms, F_ext, g=-9.807j):
    omg1 = case[-1]*np.pi/30
    rB, rC, rD, rE, rF = pos(theta1, case)
    vB, vC, vE, vF = vel(theta1, case)[:4]
    F_ext*=-np.sign(np.real(vF))
    
    vGs = [vB/2, (vB+vC)/2, (vC+vE)/2, (vE+vF)/2, vF]
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
        link2, = ax.plot(np.real([rB, rC]), np.imag([rB, rC]), 'go-')
        link3, = ax.plot(np.real([rC, rE]), np.imag([rC, rE]), 'bo-')
        link4, = ax.plot(np.real([rE, rF]), np.imag([rE, rF]), 'ko-')
        ax.plot(np.real(rD), np.imag(rD), 'bo-')
        ax.plot(np.real(rF), np.imag(rF), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.002, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.002, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.002, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.002, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.002, r'$E$')
        F = ax.text(np.real(rF), np.imag(rF)+.002, r'$F$')
        
        container.append([link1, link2, link3, link4, A, B, C, D, E, F])
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.set_ylim(ymax=case[1]+.02)
    ax.grid()
    plt.show()
    
    ani.save('mechanism_3.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(case, ms, Is, F_ext=2000, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    phi1s = np.linspace(0, 2*np.pi, frames)
    vF, aF, Meq = [], [], []
    
    for phi in phi1s:
        vF.append(np.real(vel(phi, case)[3]))
        aF.append(np.real(acc(phi, case)[3]))
        # Meq.append(reactions(phi, case, ms, Is, F_ext, g)[-1])
        Meq.append(virtual_work(phi, case, ms, F_ext))
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
    
    fig.savefig('Mechanism3', dpi=600)

#         AB   BC   CD   DE   EF   a    b   phi  n
case1 = [.02, .03, .03, .03, .06, .03, .01, 30,  600]
case2 = [.05, .20, .20, .31, .30, .10, .06, 150, 700]
case3 = [.09, .11, .07, .11, .12, .035,.025,240, 900]
case4 = [.22, .27, .17, .25, .28, .09, .055,330, 1100]

ds = 1e-3
ws = np.array([1e-2, 1e-2, 1e-2, 1e-2, 5e-2])
hs = np.array([case1[0], case1[1], case1[2]+case1[3], case1[4], 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

animate(case1)
plot_vel_acc_M(case1, ms, Is)
