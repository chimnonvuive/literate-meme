import numpy as np
import matplotlib.pyplot as plt
import kinematicAnalysis as kA
import dynamicAnalysis as dA
from matplotlib import rc


def pos(phi1, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    
    rA = 0
    rB = kA.toExp(AB, phi1)
    rC = AC*1j
    rD = rC + kA.toExp(CD, np.angle(rC-rB))
    rF = rC + kA.toExp(DF-CD, np.angle(rB-rC))
    rG = kA.toExp(AG, np.angle(rD))
    
    return rA, rB, rC, rD, rF, rG

def vel(phi1, omg1, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    
    rA, rB, rC, rD, rF, rG = pos(phi1, AB, AC, CD, DF, AG)
    v2 = kA.vel(rB, omg1, 0)
    v3_T, omg3 = kA.vCPA1(0, rB-rC, rB, omg1, 0)
    v5_T, omg5 = kA.vCPA1(0, rD, rD-rC, omg3, 0)
    v4 = kA.vel(rD-rC, omg3, 0)
    vF = kA.vel(rF-rC, omg3, 0)
    vG = kA.vel(rG, omg5, 0)
    return v2, v3_T, v5_T, v4, vF, vG, omg3, omg5

def acc(phi1, omg1, alp1, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):

    rA, rB, rC, rD, rF, rG = pos(phi1, AB, AC, CD, DF, AG)
    v3_T, v5_T = vel(phi1, omg1, AB, AC, CD, DF, AG)[1:3]
    omg3, omg5 = vel(phi1, omg1, AB, AC, CD, DF, AG)[-2:]
    a2 = kA.acc(rB, omg1, 0, 0, alp1)
    a3_T, alp3 = kA.aCPA1(0, rB-rC, rB, omg3, omg1, v3_T, 0, 0, alp1)
    a5_T, alp5 = kA.aCPA1(0, rD, rD-rC, omg5, omg3, v5_T, 0, 0, alp3)
    a4 = kA.acc(rD-rC, omg3, 0, 0, alp3)
    aF = kA.acc(rF-rC, omg3, 0, 0, alp3)
    aG = kA.acc(rG, omg5, 0, 0, alp5)
    return a2, a3_T, a5_T, a4, aF, aG, alp3, alp5

def reactions(phi1, omg1, alp1, ms, Is, M5_ext,
            g=-9.807j, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    
    rA, rB, rC, rD, rF, rG = pos(phi1, AB, AC, CD, DF, AG)
    M5_ext = -np.sign(vel(phi1, omg1, AB, AC, CD, DF, AG)[-1])*M5_ext
    a2, a3_T, a5_T, a4, aF, aG, alp3, alp5 = acc(phi1, omg1, alp1,
                                                 AB, AC, CD, DF, AG)
    aGs = np.array([a2/2, a2, (a4+aF)/2, a4, aG/2])
    Gs, FGs = ms*g, ms*aGs

    # find F45, F05 with matrix solve AX=B
    # moment equation, force equation, dot product of F45 and rD = 0
    A = dA.mat([dA.moment(rD-rG/2), dA.moment(rA-rG/2),
          [1, 0, 1, 0], 
          [0, 1, 0, 1], 
          [np.real(rD), np.imag(rD), 0, 0]], (4,4))
    b = dA.mat([[Is[-1]*alp5-M5_ext], kA.toList(FGs[-1] - Gs[-1]), [0]], 4)
    F45, F05 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F34
    F34 = FGs[-2] - (Gs[-2] + (-F45))
        
    # find F03, F23 with matrix solve AX=B
    # moment equation, force equation, dot product of F23 and rF-rC = 0
    rG3 = (rD+rF)/2
    
    A = dA.mat([dA.moment(rC-rG3), dA.moment(rB-rG3), 
                [1, 0, 1, 0], 
                [0, 1, 0, 1], 
                [0, 0, np.real(rB-rG3), np.imag(rB-rG3)]], (4,4))
    b = dA.mat([[Is[2]*alp3 - dA.toCross(rD-rG3, -F34)], 
                kA.toList(FGs[2] - Gs[2] - (-F34)), [0]], 4)
    F03, F23 = kA.toComplexes(np.linalg.solve(A, b))
    
    # find F12, F01
    F12 = FGs[1] - Gs[1] - (-F23)
    F01 = FGs[0] - Gs[0] - (-F12)
    M_eq = Is[0]*alp1 - dA.toCross(rB/2, -F12) - dA.toCross(-rB/2, F01)
    return F05, F45, F34, F03, F23, F12, F01, M_eq

# def virtual(phi1, omg1, alp1, ms, M5_ext, g=-9.807j, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    

def animate(phi1s=np.linspace(0, 360, 300)*np.pi/180,
            AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    
    from matplotlib.animation import ArtistAnimation
    
    fig, ax = plt.subplots()
    
    container = []
    for phi in phi1s:
        
        rA, rB, rC, rD, rF, rG = pos(phi, AB, AC, CD, DF, AG)
        
        pinA, = ax.plot(np.real(rC), np.imag(rC), 'bo')
        pinC, = ax.plot(np.real(rC), np.imag(rC), 'bo')
        link1, = ax.plot(np.real([rA, rB]), np.imag([rA, rB]), 'ro-')
        link3, = ax.plot(np.real([rD, rF]), np.imag([rD, rF]), 'go-')
        link5, = ax.plot(np.real([rA, rG]), np.imag([rA, rG]), 'bo-')
        A = ax.text(np.real(rA)*1.05, np.imag(rA)*1.05, 'A', fontsize=12)
        B = ax.text(np.real(rB)*1.05, np.imag(rB)*1.05, 'B', fontsize=12)
        C = ax.text(np.real(rC)*1.05, np.imag(rC)*1.05, 'C', fontsize=12)
        D = ax.text(np.real(rD)*1.05, np.imag(rD)*1.05, 'D', fontsize=12)
        F = ax.text(np.real(rF)*1.05, np.imag(rF)*1.05, 'F', fontsize=12)
        G = ax.text(np.real(rG)*1.05, np.imag(rG)*1.05, 'G', fontsize=12)
        
        container.append([pinA, pinC, link1, link3, link5, A, B, C, D, F, G])

    ax.grid()
    ax.set_aspect('equal')
    ax.set_xlim(xmin=-np.max([AG, DF]), xmax=np.max([AG, DF]))
    ax.set_ylim(ymin=-np.max([AG, DF]), ymax=np.max([AG, DF]))
    ani = ArtistAnimation(fig, container, interval=15, blit=True)
    ani.save('R_RTR_RTR.mp4', fps=60, bitrate=5000, dpi=400)

def plot_vel(omg1, phi1s=np.linspace(0, 360, 1000)*np.pi/180, AB=1, BC=1):
    
    """
    Plot velocity graph of point G
    """
    
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    vGxs = []
    for phi in phi1s:
        vGxs.append(np.real(vel(phi, omg1, AB, BC)[-3]))
        
    fig, ax = plt.subplots()
    ax.plot(phi1s, vGxs)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{BA}$ ($radian$)')
    ax.set_ylabel(r'Velocity of point G ($m/s$)')
    ax.grid()
    
    fig.savefig('Velocity_diagram_R_RTR_RTR', dpi=600, papertype='a4')
   
def plot_acc(omg1, alp1, phi1s=np.linspace(0, 360, 1000)*np.pi/180, AB=1, BC=1):
    
    """
    Plot acceleration graph of point G
    """
    
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    aGxs = []
    for phi in phi1s:
        aGxs.append(np.real(acc(phi, omg1, alp1, AB, BC)[-1]))
    
    fig, ax = plt.subplots()
    ax.plot(phi1s, aGxs)
    
    xlabels = np.linspace(0, 2, 5)*np.pi
    xticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    ax.set_xticks(xlabels)
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'Angle of $r_{BA}$ ($radian$)')
    ax.set_ylabel(r'Acceleration of point G ($m/s^2$)')
    ax.grid()
    
    fig.savefig('Acceleration_diagram_R_RTR_RTR', dpi=600, papertype='a4')

def plot_vel_acc_M(phi1s, omg1, alp1,  ms, Is, M5_ext,
            g=-9.807j, AB=0.15, AC=0.1, CD=0.15, DF=0.4, AG=0.3):
    
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    vGx, aGx, Meq = [], [], []
    vGy, aGy = [], []
    
    for phi in phi1s:
        vGx.append(np.real(vel(phi, omg1, AB, AC, CD, DF, AG)[-3]))
        aGx.append(np.real(acc(phi, omg1, alp1, AB, AC, CD, DF, AG)[-3]))
        vGy.append(np.imag(vel(phi, omg1, AB, AC, CD, DF, AG)[-3]))
        aGy.append(np.imag(acc(phi, omg1, alp1, AB, AC, CD, DF, AG)[-3]))
        Meq.append(reactions(phi, omg1, alp1, ms, Is, M5_ext, 
                             g, AB, AC, CD, DF, AG)[-1])
    fig, ax = plt.subplots(nrows=3, figsize=(8, 11))
    
    ax[0].plot(phi1s, vGx, label='x-axis')
    ax[0].plot(phi1s, vGy, label='y-axis')
    ax[1].plot(phi1s, aGx, label='x-axis')
    ax[1].plot(phi1s, aGy, label='y-axis')
    ax[2].plot(phi1s, Meq, label='Drive moment')
    
    labels = np.linspace(0, 2*np.pi, 5)
    ticks  = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
    
    ax[0].set_xticks(labels)
    ax[0].set_xticklabels(ticks)
    ax[0].set_ylabel(r'Velocity diagram ($m/s^2$)')
    ax[0].grid()
    ax[0].legend()
        
    ax[1].set_xticks(labels)
    ax[1].set_xticklabels(ticks)
    ax[1].set_ylabel(r'Acceleration diagram ($m/s^2$)')
    ax[1].grid()
    ax[1].legend()

    ax[2].set_xticks(labels)
    ax[2].set_xticklabels(ticks)
    ax[2].set_xlabel(r'Angle of $r_{B}$ ($radian$)')
    ax[2].set_ylabel(r'Drive moment ($N\cdot m$)')
    ax[2].grid()
    ax[2].legend()

    plt.tight_layout()
    
    fig.savefig('Diagram_R_RTR_RTR', dpi=600)

hs = np.array([1, 2, 1, 2, 1])*1e-2
ws = np.array([0, 5, 0, 5, 0])*1e-2
ds = np.ones(5)*1e-3
ls = np.array([1.5, 0, 4, 0, 3])*1e-1
n = 50
omg1 = n*np.pi/30
rho = 8000
ms = hs*(ws+ls)*ds*rho
Is = ms*(ws**2 + ls**2 + hs**2)/12
M5_ext = 100

animate(np.linspace(0,2*np.pi,200))
plot_vel_acc_M(np.linspace(0, 2*np.pi, 400), omg1, 0, ms, Is, M5_ext)