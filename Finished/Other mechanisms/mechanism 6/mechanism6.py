import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(phi1, case):
    
    AB, AC, BD, DE, EF, a, b = case[:-2]
    rB, rC, rF = kA.toExp(AB, phi1), -AC, (a-AC) + b*1j
    
    # find rBD
    rBD = kA.toExp(BD, np.angle(rB-rC))
    rD = rB - rBD
    
    # find rDE and rE
    angle_DE = kA.CPA4(EF, DE, rD-rF)[1][1]
    rDE = kA.toExp(DE, angle_DE)
    rE = rD-rDE
    
    return rB, rC, rD, rE, rF

def vel(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg1 = case[-1]*np.pi/30
    
    # find vB
    vB = kA.vel(rB, omg1)
    
    # find v2_T, omg2
    v2_T, omg2 = kA.vCPA1(0, rB-rC, rB, omg1)
    
    # find vE = vB + vEB
    vD = vB+kA.vel(rD-rB, omg2)#kA.vel(rD-rC, omg2, v2_T)
    omg4, omg5 = kA.vCPA4(vD, rE-rD, rE-rF)
    vE = kA.vel(rE-rF, omg5)
    return vB, vD, vE, omg2, omg4, omg5, v2_T

def acc(phi1, case):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    omg2, omg4, omg5, v2_T = vel(phi1, case)[3:]
    omg1 = case[-1]*np.pi/30
    
    # find aB
    aB = kA.acc(rB, omg1)
    
    # find angular acceleration of link2, link3
    a2_T, alp2 = kA.aCPA1(0, rB-rC, rB, omg2, omg1, v1_T=v2_T)
    
    # find aC, aE
    aD = aB+kA.acc(rD-rB, omg2, alp2)#kA.acc(rD-rC, omg2, alp2, v2_T, a2_T)
    alp4, alp5 = kA.aCPA4(aD, rE-rD, rE-rF, omg4, omg5)
    aE = kA.acc(rE-rF, omg5, alp5)
    return aB, aD, aE, alp2, alp4, alp5

def reactions(phi1, case, ms, Is, M_ext, g=-9.807j):
    
    rB, rC, rD, rE, rF = pos(phi1, case)
    aB, aD, aE, alp2, alp4, alp5 = acc(phi1, case)
    M_ext *= -np.sign(vel(phi1, case)[-2])
    
    aGs = np.array([aB/2, (aB+aD)/2, 0, (aE+aD)/2, aE/2])
    Gs, FGs = ms*g, ms*aGs

    # find F05, F34 using dyad method, moment at E
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                dA.moment(rF-rE), [0, 0],
                [0, 0], dA.moment(rD-rE)], (4, 4))
    b = dA.mat([kA.toList(FGs[3]+FGs[4]-Gs[3]-Gs[4]),
                [Is[4]*alp5 + dA.toCross((rF-rE)/2, FGs[4]-Gs[4]) - M_ext,
                 Is[3]*alp4 + dA.toCross((rD-rE)/2, FGs[3]-Gs[3])]], 4)
    F05, F24 = kA.toComplexes(np.linalg.solve(A, b))

    # find F45
    F45 = FGs[4] - Gs[4] - F05
    
    # find F03, F12 using dyad method, moment at C
    rG2 = (rB+rD)/2
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0], dA.moment(rB-rC),
                kA.toList(rB-rD), [0, 0]], (4, 4))
    b = dA.mat([kA.toList(FGs[2]+FGs[1]-Gs[2]-Gs[1]-(-F24)),
                [Is[1]*alp2 + dA.toCross([rG2-rC, -(rD-rC)], 
                                         [FGs[1]-Gs[1], -F24]), 0]], 4)
    F03, F12 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F23
    F23 = FGs[2] - Gs[2] - F03
    
    # find F01
    F01 = FGs[0] - Gs[0] - (-F12)
    
    # find drive moment of crank 1
    M_eq = dA.toCross([rB/2, -rB], [FGs[0]-Gs[0], -F12])# - dA.toCross(rB, -F12)

    return F05, F45, F23, F03, F12, F01, M_eq

def virtual_work(phi1, case, ms, M_ext, g=-9.807j):
    
    omg1, omg5 = case[-1]*np.pi/30, vel(phi1, case)[-2]
    vB, vD, vE = vel(phi1, case)[:3]
    M_ext *= -np.sign(omg5)
    
    vGs = [vB/2, (vB+vD)/2, 0, (vD+vE)/2, vE/2]
    M_eq = -(dA.toDot(ms*g, vGs) + M_ext*omg5)/omg1
    return M_eq
    
def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    phi1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for phi in phi1s:
        rB, rC, rD, rE, rF = pos(phi, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rB, rD]), np.imag([rB, rD]), 'go-')
        link4, = ax.plot(np.real([rD, rE]), np.imag([rD, rE]), 'bo-')
        link5, = ax.plot(np.real([rE, rF]), np.imag([rE, rF]), 'ko-')
        ax.plot(np.real(rD), np.imag(rD), 'b,-')
        ax.plot(np.real(rE), np.imag(rE), 'b,-')
        A = ax.text(np.real(rA), np.imag(rA)+.01, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.01, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.01, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.01, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.01, r'$E$')
        # F = ax.text(np.real(rF), np.imag(rF)+.01, r'$F$')
        container.append([link1, link2, link4, link5, A, B, C, D, E])
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.set_ylim(ymax=case[1])
    ax.grid()
    plt.show()
    ani.save('mechanism_6.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(case, ms, Is, M_ext, frames=500):
    
    phi1s = np.linspace(0, 2*np.pi, frames)
    Meq, vEx, aEx, vEy, aEy = [],[],[],[],[]
    
    for phi in phi1s:
        
        vEx.append(np.real(vel(phi, case)[2]))
        aEx.append(np.real(acc(phi, case)[2]))
        vEy.append(np.imag(vel(phi, case)[2]))
        aEy.append(np.imag(acc(phi, case)[2]))
        Meq.append(reactions(phi, case, ms, Is, M_ext)[-1])
        # Meq.append(virtual_work(phi, case, ms, M_ext))
    
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(phi1s, vEx, label='x-axis')
    ax[0].plot(phi1s, vEy, label='y-axis')
    ax[1].plot(phi1s, aEx, label='x-axis')
    ax[1].plot(phi1s, aEy, label='y-axis')
    ax[2].plot(phi1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of point E ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of point E ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($rad$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('Mechanism6', dpi=600)

#          AB   AC   BD   DE   EF   a    b   phi  n
cases = [[.06, .18, .27, .12, .08, .02, .01, 225, 600],
         [.08, .20, .30, .16, .11, .025, .13, 45, 750],
         [.05, .12, .18, .10, .065, .017, .085, 135, 950],
         [.04, .10, .16, .09, .063, .013, .065, 315, 1150]]
n = 0

ds = 1e-3
ws = np.array([1e-2, 1e-2, 5e-2, 1e-2, 1e-2])
hs = np.array([cases[n][0], cases[n][2], 2e-2, 
               cases[n][3], cases[n][4]])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

M_ext = 3000

animate(cases[n])
plot_vel_acc_M(cases[n], ms, Is, M_ext)