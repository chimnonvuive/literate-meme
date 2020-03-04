import numpy as np
import kinematicAnalysis as kA
import dynamicAnalysis as dA
import matplotlib.pyplot as plt
from matplotlib import rc


rc('text', usetex=True)

def pos(theta1, case):
    
    rB = kA.toExp(case[0], theta1)
    rD = case[5] + case[6]*1j
    rBD = rB - rD

    # rCD + rBC = rBD
    CD_angle = kA.CPA4(case[2], case[1], rBD)[0][1]
    rCD = kA.toExp(case[2], CD_angle)
    rED = kA.toExp(case[2]+case[3], CD_angle)
    rC, rE = rD + rCD, rD + rED
    
    res1, res2 = kA.CPA2(np.pi/2, case[4], case[7]+rD+rED)
    r0F_1, rFE_1 = kA.toExp(res1[0], np.pi/2), kA.toExp(case[4], res1[1])
    r0F_2, rFE_2 = kA.toExp(res2[0], np.pi/2), kA.toExp(case[4], res2[1])
    
    return (rB, rC, rD, rE, rCD, rED, r0F_1, rFE_1), \
           (rB, rC, rD, rE, rCD, rED, r0F_2, rFE_2)

def vel(theta1, case, variation=0):
    
    rB, rC, rD, rE, rCD, rED, r0F, rFE = pos(theta1, case)[variation]
    omg1 = case[-1]*np.pi/30
    vB = kA.vel(rB, omg1)
    omg2, omg3 = kA.vCPA4(vB, rC-rB, rCD)
    vE = kA.vel(rED, omg3, 0)
    v05, omg4 = kA.vCPA2(-vE, r0F, rFE)
    v5 = kA.vel(r0F, v_T=v05)
    return vB, v05, v5, omg2, omg3, omg4

def acc(theta1, case, variation=0):
    
    rB, rC, rD, rE, rCD, rED, r0F, rFE = pos(theta1, case)[variation]
    v05, v5, omg2, omg3, omg4 = vel(theta1, case, variation)[1:]
    omg1 = case[-1]*np.pi/30
    
    aB = kA.acc(rB, omg1)
    alp2, alp3 = kA.aCPA4(aB, rC-rB, rCD, omg2, omg3)
    aC = kA.acc(rCD, omg3, alp3)
    aE = kA.acc(rED, omg3, alp3)
    a05, alp4 = kA.aCPA2(-aE, r0F, rFE, omg2=omg4, v1_T=v05)
    a5 = kA.acc(r0F, v_T=v05, a_T=a05)
    return aB, aC, aE, a5, a05, alp2, alp3, alp4

def reactions(theta1, case, ms, Is, F_ext, g=-9.807j, variation=0):
    
    rB, rC, rD, rE, rCD, rED, r0F, rFE = pos(theta1, case)[variation]
    aB, aC, aE, a5, a05, alp2, alp3, alp4 = acc(theta1, case, variation)
    F = F_ext*-np.sign(vel(theta1, case, variation)[2])
    
    aGs = np.array([aB/2, (aB+aC)/2, aE/2, (aE+a5)/2, a5])
    Gs, FGs = ms*g, ms*aGs
    
    # find F05, F34 using dyad method
    A = dA.mat([[1, 0, 1, 0],
                [0, 1, 0, 1],
                [0, 0], dA.moment(-rFE),
                [0, 1, 0, 0]], (4, 4))
    b = dA.mat([kA.toList(FGs[3]+FGs[4]-Gs[3]-Gs[4]-F),
                [Is[3]*alp4 + dA.toCross(-rFE/2, FGs[3]-Gs[3]), 0]], 4)
    F05, F34 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F45
    F45 = FGs[4] - Gs[4] - F05 - F_ext
    
    # find F03, F12 using dyad method at C
    rG3, rG2 = (rE+rD)/2, (rB+rC)/2
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
    
    # find F01
    F01 = FGs[0] - Gs[0] - (-F12)
    
    # find drive moment of crank 1
    M_eq = dA.toCross(rB/2, FGs[0]-Gs[0]) - dA.toCross(rB, -F12)
     
    return F05, F45, F23, F03, F12, F01, M_eq
def virtual(theta1, case, ms, F_ext, g=-9.807j, variation=0):
    omg1 = case[-1]*np.pi/30
    rB, rC, rD, rE, rCD, rED, r0F, rFE = pos(theta1, case)[variation]
    vB, v05, v5, omg2, omg3, omg4 = vel(theta1, case, variation)
    vC, vE = kA.vel(rC-rD, omg3), kA.vel(rE-rD, omg3)
    F = F_ext*-np.sign(np.imag(v5))
    vGs = [vB/2, (vB+vC)/2, vE/2, (vE+v5)/2, v5]
    M_eq = -np.sum(dA.toDot(ms*g, vGs) + dA.toDot(F, v5))/omg1
    return M_eq
    
def animate(case, variation=0):
    
    from matplotlib.animation import ArtistAnimation
    
    theta1s = np.linspace(0, 2*np.pi, 300)
    fig, ax = plt.subplots()
    rA, container = 0, []
    
    for theta1 in theta1s:
        rB, rC, rD, rE, rCD, rED, r0F, rFE = pos(theta1, case)[variation]
        rF = -case[7] + r0F
        
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link2, = ax.plot(np.real([rB, rC]), np.imag([rB, rC]), 'go-')
        link3, = ax.plot(np.real([rD, rE]), np.imag([rD, rE]), 'bo-')
        link4, = ax.plot(np.real([rE, rF]), np.imag([rE, rF]), 'ko-')
        ax.plot(np.real(rF), np.imag(rF), 'b,-')
        
        A = ax.text(np.real(rA), np.imag(rA)+.01, r'$A$')
        B = ax.text(np.real(rB), np.imag(rB)+.01, r'$B$')
        C = ax.text(np.real(rC), np.imag(rC)+.01, r'$C$')
        D = ax.text(np.real(rD), np.imag(rD)+.01, r'$D$')
        E = ax.text(np.real(rE), np.imag(rE)+.01, r'$E$')
        F = ax.text(np.real(rF), np.imag(rF)+.01, r'$F$')
        # G2 = ax.text(np.real(rG2), np.imag(rG2)+.01, r'$G_2$')
        
        container.append([link1, link2, link3, link4, A, B, C, D, E, F])
    
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
    ax.set_aspect('equal')
    ax.grid()
    plt.show()
    
    ani.save('mechanism_2_{0:d}.mp4'.format(variation), 
              fps=60, bitrate=5000, dpi=400)

def plot_vel_acc_M(case, ms, Is, F_ext=2000j, 
                   g=-9.807j, variation=0, frames=1000):
    
    # plot velocity and acceleration of point G2
    theta1s = np.linspace(0, 2*np.pi, frames)
    v0F, a0F, Meq = [], [], []
    
    for theta1 in theta1s:
        
        v0F.append(np.imag(vel(theta1, case, variation)[2]))
        a0F.append(np.imag(acc(theta1, case, variation)[3]))
        # Meq.append(reactions(theta1, case, ms, Is, F_ext, g, variation)[-1])
        Meq.append(virtual(theta1, case, ms, F_ext, g, variation))
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
    
    fig.savefig('mechanism_2_{0:d}'.format(variation), dpi=600)

#         AB   BC   CD   CE   EF   a    b    c   theta1  n
case1 = [.15, .40, .37, .23, .23, .30, .45, .37, 330, 500]
# case2 = [.10, .27, .25, .15, .15, .20, .30, .25, 120, 800]
# case3 = [.12, .32, .30, .18, .18, .24, .36, .30, 240, 900]
# case4 = [.20, .55, .50, .30, .30, .40, .60, .50, 300, 1e3]
variation=1
ds = 1e-3
ws = np.array([1e-2, 1e-2, 1e-2, 1e-2, 5e-2])
hs = np.array([case1[0], case1[1], case1[2]+case1[3], case1[4], 2e-2])
rho = 8000
ms = ds*ws*hs*rho
Is = ms*(ws**2 + hs**2)/12

# animate(case1, variation)
plot_vel_acc_M(case1, ms, Is, variation=variation)