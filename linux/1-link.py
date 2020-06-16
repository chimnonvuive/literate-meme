import numpy as np
from sympy import symbols, solve, sin, cos
from sympy.physics.vector import ReferenceFrame
from sympy.physics.mechanics import dynamicsymbols, init_vprinting

init_vprinting

O = ReferenceFrame('O')
t, m, L, g = symbols('t m L g')
theta = dynamicsymbols('theta')


rC = L/2 * cos(theta) * O.x + L/2 * sin(theta) * O.y

IC = m*L**2/12
IO = IC + m*L**2/4
omg = theta.diff() * O.z
alp = omg.diff(t, O)

vC = rC.diff(t, O)
aC = vC.diff(t, O)

G = - m * g * O.y
F01 = m * aC - G
eq = IO * alp - (rC ^ G)

# output is -16.1⋅cos(θ(t))
rhs = solve(eq.to_matrix(O)[2], theta.diff(t,2))[0].subs([(g, 32.2), (L, 3)])

# solve the equations

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
from matplotlib import rc

rc('text', usetex=True)

fun = lambda t, y: [y[1], -16.1*np.cos(y[0])]
res = solve_ivp(fun, [0, 10], [np.pi/4, 0], t_eval=np.linspace(0, 5, 1000))

t = res.t
y = res.y

# plot the graphs

fig, ax = plt.subplots(nrows=3, figsize=(4,5))
plots = [[t, y[0]], [t, y[1]], [y[0], y[1]]]
ylabels = [r'$\theta$', r'$\theta$', r'$\omega$']
xlabels = [r'$t$', r'$t$', r'$\theta$']

for i in range(3):
    ax[i].plot(*plots[i], 'b,-')
    ax[i].grid()
    ax[i].set_xlabel(xlabels[i])
    ax[i].set_ylabel(ylabels[i])
plt.tight_layout()
plt.show

# animate the body

fig2, ax2 = plt.subplots()
container = []

for angle in y[0]:
    qC = rC.to_matrix(O).subs([(theta, angle), (L, 3)])
    link, = ax2.plot([0, qC[0]], [0, qC[1]], 'ro-')
    A = ax2.text(0, 0, r'$A$')
    C = ax2.text(qC[0], qC[1], r'$C$')
    container.append([link, A, C])

ani = ArtistAnimation(fig2, container, interval=20, blit=True)
ax2.grid()
ax2.set_aspect('equal')
plt.show
