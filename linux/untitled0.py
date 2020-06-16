import numpy as np

from scipy.integrate import solve_ivp

from sympy import lambdify, symbols, sin, cos
from sympy.solvers.solveset import linsolve
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import dynamicsymbols, init_vprinting

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.animation import ArtistAnimation



init_vprinting()
rc('text', usetex=True)

O = ReferenceFrame('O')
l1, l2, m1, m2, g, t = symbols('L_1 L_2 m_1 m_2 g t')
q1, q2 = dynamicsymbols('q_1 q_2')
interval = [0, 5]
variables = [q1, q1.diff(), q2, q2.diff()]
num = [(l1, 1), (l2, 0.5), (m1, 1), (m2, 1), (g, 10)]

omg1, omg2 = q1.diff()*O.z, q2.diff()*O.z
alp1, alp2 = omg1.diff(t, O), omg2.diff(t, O)
IC1, IC2 = m1*l1**2/12, m2*l2**2/12

"PVA of rC1, rC2"
# position
rC1 = l1/2 * (cos(q1)*O.x + sin(q1)*O.y)
rC2 = 2*rC1 + l2/2 * (cos(q2)*O.x + sin(q2)*O.y)

# velocity
vC1 = rC1.diff(t, O)
vC2 = rC2.diff(t, O)

# acceleration
aC1 = vC1.diff(t, O)
aC2 = vC2.diff(t, O)

"Force Analysis"
G1 = -m1*g * O.y
G2 = -m2*g * O.y

F12 = m2*aC2 - G2
F01 = m1*aC1 - G1 - (-F12)

eq1 = (-IC1*alp1 + (rC1 ^ -F12) + (-rC1 ^ F01)).to_matrix(O)[2].subs(num)
eq2 = (-IC2*alp2 + ((2*rC1-rC2) ^ F12)).to_matrix(O)[2].subs(num)

sol = linsolve((eq1, eq2), (q1.diff(t, 2), q2.diff(t, 2)))

# "plotting graphs"
   
# dq1 = lambdify(variables, q1.diff(), 'numpy')
# dq2 = lambdify(variables, q2.diff(), 'numpy')
# ddq1 = lambdify(variables, list(sol)[0][0], 'numpy')
# ddq2 = lambdify(variables, list(sol)[0][1], 'numpy')

# fun = lambda t1, y: [dq1(*y), ddq1(*y), dq2(*y), ddq2(*y)]

# res = solve_ivp(fun, interval, [-np.pi/4, 0, -np.pi/3, 0], 
#                 t_eval=np.linspace(*interval, 1000))

# t = res.t
# q1_arr, q2_arr = res.y[0], res.y[2]

# fig, ax = plt.subplots()
# ax.plot(t, q1_arr*180/np.pi, label=r'$q_1$ (deg)')
# ax.plot(t, q2_arr*180/np.pi, label=r'$q_2$ (deg)')
# ax.set_xlabel('t (s)')
# ax.set_ylabel('angle of 2 links (deg)')
# ax.legend()
# ax.grid()
# plt.tight_layout()
# plt.show
# plt.savefig('Figure1.png')

# fig2, ax2 = plt.subplots()
# container = []

# for i in range(len(q1_arr)):
#     qB = 2*rC1.to_matrix(O).subs([(q1, q1_arr[i]), *num])
#     xC = qB[0] + (l2 * cos(q2)).subs([(q2, q2_arr[i]), *num])
#     yC = qB[1] + (l2 * sin(q2)).subs([(q2, q2_arr[i]), *num])
#     link1, = ax2.plot([0, qB[0]], [0, qB[1]], 'ro-')
#     link2, = ax2.plot([qB[0], xC], [qB[1], yC], 'bo-')
#     ax2.plot(xC, yC, 'b,')
#     container.append([link1, link2])
    
# ani = ArtistAnimation(fig2, container, interval=20, blit=True)
# ax2.set_aspect('equal')
# plt.show()
# ani.save('Figure2.mp4')