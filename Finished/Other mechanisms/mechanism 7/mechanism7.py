import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(phi1, case):
    
    AB, AD, BC, a = case[:-2]
    CD_max = np.sqrt((AB/AD*(a+AD))**2 + (a+AD)**2)
    rB, rD = kA.toExp(AB, phi1), AD
    rBD = rB - rD
    # rCD + rBC = rBD
    rCD_max, rCD = kA.toExp([CD_max, np.abs(rBD)+BC], np.angle(rBD))
    rD_max = rD - (rCD_max-rCD)
    
    rC = rD + rCD
    rE = -a + np.imag(rC)*1j
    
    rC_max = rE + a
    return rB, rC, rD, rD_max, rE, rC_max

def vel(phi1, case):
    
    rB, rC, rD, rD_max, rE, rC_max = pos(phi1, case)
    omg1 = case[-1]*np.pi/30
    
    # find vB
    vB = kA.vel(rB, omg1)
    
    # angular velocity of link 2
    v2_T, omg2 = kA.vCPA1(0, rB-rD, rB, omg1)
    vC = kA.vel(rC-rD, omg2, v2_T)
    vE = vC_max = np.imag(vC)*1j
    vD_max = kA.vel(rD_max-rD, omg2, v2_T)
    return vB, vC, vE, vD_max, vC_max, omg2, v2_T

def acc(phi1, case):
    
    rB, rC, rD, rD_max, rE, rC_max = pos(phi1, case)
    vB, vC, vE, vD_max, vC_max, omg2, v2_T = vel(phi1, case)
    omg1 = case[-1]*np.pi/30
    
    # find aB
    aB = kA.acc(rB, omg1)
    
    # find angular acceleration of link2, link3
    a2_T, alp2 = kA.aCPA1(0, rB-rD, rB, omg2, omg1, v1_T=v2_T)
    
    # find aC, aE
    aC = kA.acc(rC-rD, omg2, alp2, v2_T, a2_T)
    aE = aC_max = np.imag(aC)*1j
    aD_max = kA.acc(rD_max-rD, omg2, alp2, v2_T, a2_T)
    return aB, aC, aE, aD_max, aC_max, alp2, a2_T

def reactions(phi1, case, ms, Is, F_ext, g=-9.807j):
    
    rB, rC, rD, rD_max, rE, rC_max = pos(phi1, case)
    aB, aC, aE, aD_max, aC_max, alp2, a2_T = acc(phi1, case)
    F = F_ext*-np.sign(vel(phi1, case)[2])
    aGs = np.array([aB/2, (aD_max+aC)/2, 0, aC, (aE+aC_max)/2])
    Gs, FGs = ms*g, ms*aGs
    
    # find F05, F24 using dyad method, moment at C
    rG5 = (rE+rC_max)/2
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                dA.moment(rE-rC), [0, 0],
                [0, 0, 1, 0]], (4, 4))
    b = dA.mat([kA.toList(FGs[3]+FGs[4]-Gs[3]-Gs[4]-F),
                [dA.toCross(rG5-rC, FGs[4]-Gs[4]), 0]], 4)
    F05, F24 = kA.toComplexes(np.linalg.solve(A, b))
        
    # find F45
    F45 = FGs[4] - Gs[4] - F05 - F
    
    # find F32, F12 using dyad method, moment at D
    rG2 = (rC+rD_max)/2
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                kA.toList(rB-rD), [0, 0],
                [0, 0], dA.moment(rB-rD)], (4, 4))
    b = dA.mat([kA.toList(FGs[1]-Gs[1]-(-F24)),
                [0, dA.toCross(rG2-rD, FGs[1]-Gs[1]) - dA.toCross(rC-rD, -F24)]], 4)
    F32, F12 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F03
    F03 = FGs[2] - Gs[2] - (-F32)
    
    # find F01
    F01 = FGs[0] - Gs[0] - (-F12)
    
    # find drive moment of crank 1
    M_eq = dA.toCross(rB/2, FGs[0]-Gs[0]) - dA.toCross(rB, -F12)
    
    return F05, F45, F32, F03, F12, F01, M_eq

def virtual_work(phi1, case, ms, F_ext, g=-9.807j):
    
    omg1 = case[-1]*np.pi/30
    vB, vC, vE, vD_max, vC_max, omg2, v2_T = vel(phi1, case)
    F = F_ext*-np.sign(vel(phi1, case)[2])
        
    vGs = [vB/2, (vD_max+vC)/2, 0, vC, (vE+vC_max)/2]
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F, vE))/omg1
    return M_eq
    
def animate(case):
    
    from matplotlib.animation import ArtistAnimation
    
    fig, ax = plt.subplots()
    rA, container = 0, []
    phi1s = np.linspace(0, 2*np.pi, 500)
        
    for phi in phi1s:
        rB, rC, rD, rD_max, rE, rC_max = pos(phi, case)
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rD_max, rC]), np.imag([rD_max, rC]), 'go-')
        link5, = ax.plot(np.real([rE, rC_max]), np.imag([rE, rC_max]), 'bo-')
        ax.plot(np.real(rE), np.imag(rE), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.01, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.01, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.01, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.01, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.01, r'$E$')
        
        container.append([link1, link2, link5, A, B, C, D, E])
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    plt.show()
    ani.save('mechanism_7.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(case, ms, Is, F_ext, 
                   g=-9.807j, frames=1000):
    
    # plot velocity and acceleration of point G2
    phi1s = np.linspace(0, 2*np.pi, frames)
    vE, aE, Meq = [], [], []
    
    for phi in phi1s:
        
        vE.append(np.imag(vel(phi, case)[2]))
        aE.append(np.imag(acc(phi, case)[2]))
        # Meq.append(reactions(phi, case, ms, Is, F_ext, g)[-1])
        Meq.append(virtual_work(phi, case, ms, F_ext, g))
    
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(phi1s, vE, label='y-axis')
    ax[1].plot(phi1s, aE, label='y-axis')
    ax[2].plot(phi1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram of E ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram of E ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($rad$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('Mechanism7', dpi=600)

#          AB   AD   BC   a   phi  n
cases = [[.08, .20, .10, .2, 300, 1000],
         [.06, .15, .08, .50, 210, 1200],
         [.10, .25, .12, .50, 150, 1300],
         [.12, .30, .15, .50,  30, 1500]]

n = 0

CD_max = np.sqrt((cases[n][0]/cases[n][1]*(cases[n][3]+cases[n][1]))**2 +
                 (cases[n][3]+cases[n][1])**2)

ds = 1e-3
ws = np.array([1e-2, 1e-2, 5e-2, 5e-2, 1e-2])
hs = np.array([cases[n][0], CD_max, 2e-2, 2e-2, cases[n][3]])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

F_ext = 2000j

animate(cases[n])
plot_vel_acc_M(cases[n], ms, Is, F_ext)