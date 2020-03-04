from sympy import *
from symbolicAnalysis import *
from sympy.abc import x, y, z, t
import matplotlib.pyplot as plt
import numpy as np
init_printing(pretty_print=True, )

theta = Function('theta')(t)
AB, AC, DC, DF, AG = .15, .1, .15, .4, .3
rA, rC, omg1, alp1 = vec([0,0,0]), vec([0,AC,0]), vec([0,0,5/3*pi]), vec([0,0,0])

def pos(phi1):
    
    rB = AB*vec([cos(phi1), sin(phi1), 0])
    pm = sign(rB[0])
    
    # find rD
    BC_angle = atan((rB-rC)[1]/(rB-rC)[0]) + pi
    rD = rC + pm*DC*vec([cos(BC_angle), sin(BC_angle), 0])
    
    # find rF
    BC_angle = atan((rC-rD)[1]/(rC-rD)[0])
    rF = rC + pm*(DF-DC)*vec([cos(BC_angle), sin(BC_angle), 0])
    
    # find rG
    D_angle = atan(rD[1]/rD[0])
    rG = -pm*AG*vec([cos(D_angle), sin(D_angle), 0])

    return rB, rD, rF, rG

def vel(phi1, omg1):
    
    rB, rD, rF, rG = pos(phi1)
    
    # find omg3, vB3B2
    omg3, vB3B2 = vec([0,0,x]), vec([y, z, 0])
    rCB_perp = GramSchmidt([rC-rB, vec([0,1,0]), vec([0,0,1])])[1]
    eq1 = totVel(rB-rC, omg3) - totVel(rB, omg1, v_rel=vB3B2)
    eq2 = rCB_perp.dot(vB3B2)
    sol = linsolve([eq1[0], eq1[1], eq2], [x, y, z])
    omg3 = omg3.subs(x, float(sol.args[0][0]))
    vB3B2 = vB3B2.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    
    # find omg5, vD5D4
    omg5, vD5D4 = vec([0,0,x]), vec([y, z, 0])
    rG_perp = GramSchmidt([rG, vec([0,1,0]), vec([0,0,1])])[1]
    eq1 = totVel(rD, omg5) - totVel(rD-rC, omg3, v_rel=vD5D4)
    eq2 = rG_perp.dot(vD5D4)
    sol = linsolve([eq1[0], eq1[1], eq2], [x, y, z])
    omg5 = omg5.subs(x, float(sol.args[0][0]))
    vD5D4 = vD5D4.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    
    # find vD, vF, vG
    vB = totVel(rB, omg1)
    vD = totVel(rD-rC, omg3)
    vF = totVel(rF-rC, omg3)
    vG = totVel(rG, omg5)
    
    return vB, vD, vF, vG, omg3, omg5, vB3B2, vD5D4
    
def acc(phi1, omg1, alp1):
    
    rB, rD, rF, rG = pos(phi1)
    omg3, omg5, vB3B2, vD5D4 = vel(phi1, omg1)[-4:]
    aB = totAcc(rB, omg1, alp1)

    
    # find alp3
    alp3, aB3B2r = vec([0,0,x]), vec([y,z,0])
    rCB_perp = GramSchmidt([rC-rB, vec([0,1,0]), vec([0,0,1])])[1]
    eq1 = totAcc(rB-rC, omg3, alp3) - totAcc(vec([0,0,0]), omg3, alp3, aB, 
                                             v_rel=vB3B2, a_rel=aB3B2r)
    eq2 = rCB_perp.dot(aB3B2r)

    sol = linsolve([eq1[0], eq1[1], eq2], [x,y,z])
    alp3 = alp3.subs(x, float(sol.args[0][0]))
    aB3B2r = aB3B2r.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    
    # find aD    
    aD = totAcc(rD-rC, omg3, alp3)
    
    # find alp5
    alp5, aD5D4r = vec([0,0,x]), vec([y, z, 0])
    rG_perp = GramSchmidt([rG, vec([0,1,0]), vec([0,0,1])])[1]
    eq1 = totAcc(rD, omg5, alp5) - totAcc(vec([0,0,0]), omg5, alp5, aD,
                                          v_rel=vD5D4, a_rel=aD5D4r)
    eq2 = rG_perp.dot(aD5D4r)
    sol = linsolve([eq1[0], eq1[1], eq2], [x, y, z])
    alp5 = alp5.subs(x, float(sol.args[0][0]))
    aD5D4r = aD5D4r.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    
    # find aD, aF, aG
    aF = totAcc(rF-rC, omg3, alp3)
    aG = totAcc(rG, omg5, alp5)
    
    return aB, aD, aF, aG, alp3, alp5, aB3B2r, aD5D4r


def d_pos(phi1):
    
    rB = AB*vec([cos(theta), sin(theta), 0])
    pm = sign(rB[0].subs(theta, phi1))
    
    # find rD
    CB_angle = atan((rC-rB)[1]/(rC-rB)[0])
    rD = rC + pm*DC*vec([cos(CB_angle), sin(CB_angle), 0])
    
    # find rF
    BC_angle = atan((rB-rC)[1]/(rB-rC)[0])
    rF = rC + pm*(DF-DC)*vec([cos(BC_angle), sin(BC_angle), 0])
    
    # find rG
    D_angle = atan(rD[1]/rD[0])
    rG = -pm*AG*vec([cos(D_angle), sin(D_angle), 0])
    
    return rB, rD, rF, rG
        # rB.subs(theta, phi1), rD.subs(theta, phi1),\
        # rF.subs(theta, phi1), rG.subs(theta, phi1)
        
def d_vel(phi1, omg1):
    
    rB, rD, rF, rG = d_pos(phi1)
    vB, vD, vF, vG = rB.diff(), rD.diff(), rF.diff(), rG.diff()
    
    # find omg3
    omg3 = vec([0,0,x])
    eq = vD - omg3.cross(rD-rC)
    sol = solveset(eq[0], x)
    omg3 = omg3.subs(x, sol.args[0])
    
    # find omg5
    omg5 = vec([0,0,x])
    eq = vG - omg5.cross(rG)
    sol = solveset(eq[0], x)
    omg5 = omg5.subs(x, sol.args[0])
    
    return vB, vD, vF, vG, omg3, omg5,\
        vB.subs([(theta.diff(), omg1[2]), (theta, phi1)]),\
            vD.subs([(theta.diff(), omg1[2]), (theta, phi1)]),\
                vF.subs([(theta.diff(), omg1[2]), (theta, phi1)]),\
                    vG.subs([(theta.diff(), omg1[2]), (theta, phi1)]),\
                        omg3.subs([(theta.diff(), omg1[2]), (theta, phi1)]),\

def d_acc(phi1, omg1, alp1):
    
    rB, rD, rF, rG = d_pos(phi1)
    aB, aD, aF, aG = rB.diff(t, 2), rD.diff(t, 2), rF.diff(t, 2), rG.diff(t, 2)
    omg3, omg5 = d_vel(phi1, omg1)[-2:]
    # alp3, alp5 = omg3.diff(), omg5.diff()
    
    # find alp3
    alp3 = vec([0,0,x])
    eq = aD - totAcc(rD-rC, omg3, alp3)
    sol = solveset(eq[0], x)
    alp3 = alp3.subs(x, sol.args[0])
    
    # find alp5
    alp5 = vec([0,0,x])
    eq = aG - totAcc(rG, omg5, alp5)
    sol = solveset(eq[0], x)
    alp5 = alp5.subs(x, sol.args[0])
    
    return aB, aD, aF, aG, alp3, alp5,\
        aB.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\
            aD.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\
                aF.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\
                    aG.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\
                        alp3.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\
                            alp5.subs([(theta.diff(t, 2), alp1[2]), (theta.diff(), omg1[2]), (theta, phi1)]),\


def c_vel(phi1, omg1):
    
    rB, rD, rF, rG = pos(phi1)
    rCB_perp = GramSchmidt([rC-rB, vec([0,1,0]), vec([0,0,1])])[1]
    rG_perp = GramSchmidt([rG, vec([0,1,0]), vec([0,0,1])])[1]
    
    # contour 1
    # find omg3, vB3B2
    omg12, v23, omg30 = vec([0,0,x]), vec([y,z,0]), vec([0,0,t])
    eq1 = omg1 + omg12 + omg30
    eq2 = rB.cross(omg12) + rC.cross(omg30) + v23
    eq3 = rCB_perp.dot(v23)
    sol = linsolve([eq1[2], eq2[0], eq2[1], eq3], [x, y, z, t])
    omg12 = omg12.subs(x, float(sol.args[0][0]))
    v23 = v23.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    omg30 = omg30.subs(t, float(sol.args[0][3]))
    
    # contour 2
    # find omg5, v45
    omg50, omg34, v45 = vec([0,0,x]), vec([0,0,y]), vec([z,t,0])
    eq1 = -omg30 + omg34 + omg50
    eq2 = rC.cross(-omg30) + rD.cross(omg34) + v45
    eq3 = rG_perp.dot(v45)
    sol = linsolve([eq1[2], eq2[0], eq2[1], eq3], [x,y,z,t])
    omg50 = omg50.subs(x, float(sol.args[0][0]))
    omg34 = omg34.subs(y, float(sol.args[0][1]))
    v45 = v45.subs([(z, float(sol.args[0][2])), (t, float(sol.args[0][3]))])
    
    # return omg3, omg5
    return -omg30, -omg50, v23, v45

def c_acc(phi1, omg1, alp1):
    
    rB, rD, rF, rG = pos(phi1)
    omg3, omg5, v23, v45 = c_vel(phi1, omg1)
    rCB_perp = GramSchmidt([rC-rB, vec([0,1,0]), vec([0,0,1])])[1]
    rG_perp = GramSchmidt([rG, vec([0,1,0]), vec([0,0,1])])[1]
    
    # contour 1
    # find alp3, aB3B2r
    alp12, a23, alp30 = vec([0,0,x]), vec([y,z,0]), vec([0,0,t])
    eq1 = alp1 + alp12 + alp30
    eq2 = rB.cross(alp12) + rC.cross(alp30) + a23 + 2*omg3.cross(v23) +\
        omg1.cross(omg1.cross(rB)) + omg3.cross(omg3.cross(rB-rC))
    eq3 = rCB_perp.dot(a23)
    sol = linsolve([eq1[2], eq2[0], eq2[1], eq3], [x,y,z,t])
    alp12 = alp12.subs(x, float(sol.args[0][0]))
    a23 = a23.subs([(y, float(sol.args[0][1])), (z, float(sol.args[0][2]))])
    alp30 = alp30.subs(t, float(sol.args[0][3]))
    
    # contour 2
    # find alp5, a45
    alp50, alp34, a45 = vec([0,0,x]), vec([0,0,y]), vec([z,t,0])
    eq1 = -alp30 + alp34 + alp50
    eq2 = rC.cross(-alp30) + rD.cross(alp34) + a45 + 2*omg5.cross(v45) +\
        omg3.cross(omg3.cross(rD-rC)) + omg5.cross(omg5.cross(rD))
    eq3 = rG_perp.dot(a45)
    sol = linsolve([eq1[2], eq2[0], eq2[1], eq3], [x,y,z,t])
    alp50 = alp50.subs(x, float(sol.args[0][0]))
    alp34 = alp34.subs(y, float(sol.args[0][1]))
    a45 = a45.subs([(z, float(sol.args[0][2])), (t, float(sol.args[0][3]))])
    
    # return alp3, alp5
    return -alp30, -alp50, a23, a45
    
def animate(phi1s):
    
    from matplotlib.animation import ArtistAnimation
    
    fig, ax = plt.subplots()
    
    container = []
    for phi in phi1s:
        
        rB, rD, rF, rG = pos(phi)
        
        pinA, = ax.plot(rA[0], rA[1], 'bo')
        pinC, = ax.plot(rC[0], rC[1], 'bo')
        link1, = ax.plot([rA[0], rB[0]], [rA[1], rB[1]], 'ro-')
        link3, = ax.plot([rF[0], rD[0]], [rF[1], rD[1]], 'go-')
        link5, = ax.plot([rA[0], rG[0]], [rA[1], rG[1]], 'bo-')
        A = ax.text(rA[0]+.05, rA[1]+.05, 'A', fontsize=12)
        B = ax.text(rB[0]+.05, rB[1]+.05, 'B', fontsize=12)
        C = ax.text(rC[0]+.05, rC[1]+.05, 'C', fontsize=12)
        D = ax.text(rD[0]+.05, rD[1]+.05, 'D', fontsize=12)
        F = ax.text(rF[0]+.05, rF[1]+.05, 'F', fontsize=12)
        G = ax.text(rG[0]+.05, rG[1]+.05, 'G', fontsize=12)
        
        container.append([pinA, pinC, link1, link3, link5, A, B, C, D, F, G])

    ax.grid()
    ax.set_aspect('equal')
    ani = ArtistAnimation(fig, container, interval=20, blit=True)
def plot_vel_acc_M(phi1s, omg1, alp1):
    from matplotlib import rc
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    
    vGx, aGx = [], []
    vGy, aGy = [], []
    
    for phi in phi1s:
        vGx.append(np.real(vel(phi, omg1)[-1]))
        aGx.append(np.real(acc(phi, omg1, alp1)[-3]))
        vGy.append(np.imag(vel(phi, omg1)[-1]))
        aGy.append(np.imag(acc(phi, omg1, alp1, AB, AC, CD, DF, AG)[-3]))

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
    
animate(np.linspace(0,2*np.pi,200))
plot_vel_acc_M(np.linspace(0,2*np.pi,200), 1,-1)