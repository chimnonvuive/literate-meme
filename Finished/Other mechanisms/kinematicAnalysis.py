import numpy as np


def toExp(rXs, thetas):
    
    """
    Every complex number can be written in polar form:
        z = a + jb = sqrt(a**2 + b**2)*exp(arctan(b/a))
    """
    
    return np.abs(rXs)*np.exp(1j*(np.angle(rXs) + thetas))

def toList(rX):
    
    """
    Turn a complex number into column matrix:
        z = a + jb = [a, b]
    """
    
    return np.array([np.real(rX), np.imag(rX)])

def toComplexes(results):
    
    """
    For convenience, turn a list of numbers into complex numbers
    [real(r1), imag(r1), real(r2), imag(r2), ... , real(rn), imag(rn)]
    = r1, r2, ... , rn
    """
    
    return [results[i]+results[i+1]*1j for i in range(0, len(results), 2)]

"""
Displacement vector problem
Consider a triangle ABC with a + b = c
Every planar problem derives from 4 cases below. In these cases, only 2
unknowns are considered since a complex vector has 2 components: real and 
imaginary. We will look further into finding norm and argument of the
vectors based on these 4 cases.
"""

def CPA1(a, b):
    
    """
    Complex Polar Algebra, method 1, Theory of Machines and Mechanisms
    The unknowns are C, C_angle
    """
    
    A, A_angle = np.abs(a), np.angle(a)
    B, B_angle = np.abs(b), np.angle(b)
    C = np.sqrt(A**2 + B**2 + 2*A*B*np.cos(B_angle-A_angle))
    C_angle = np.arctan((A*np.sin(A_angle) + B*np.sin(B_angle))/\
                        (A*np.cos(A_angle) + B*np.cos(B_angle)))
    return C, C_angle

def CPA2(A_angle, B, c):
    
    """
    Complex Polar Algebra, method 2, Theory of Machines and Mechanisms
    The unknowns are A, B_angle
    Using this case leads to 2 different results. Be sure to check out for
    appropriate conditions to come up with the right solution.    
    """
    
    C, C_angle = np.abs(c), np.angle(c)
    arcsin_res = np.array([np.arcsin(C*np.sin(C_angle-A_angle)/B),
                           np.pi-np.arcsin(C*np.sin(C_angle-A_angle)/B)])
    results = []
    for arcsin in arcsin_res:
        B_angle = A_angle + arcsin
        A = C*np.cos(C_angle-A_angle) - B*np.cos(B_angle-A_angle)
        results.append(np.array([A, B_angle]))
    return results[0], results[1]

def CPA3(A_angle, B_angle, c):
    
    """
    Complex Polar Algebra, method 3, Theory of Machines and Mechanisms
    The unknowns are A, B
    """
    
    C, C_angle = np.abs(c), np.angle(c)
    B = C*np.sin(C_angle-A_angle)/np.sin(B_angle-A_angle)
    A = C*np.sin(C_angle-B_angle)/np.sin(A_angle-B_angle)
    return A, B

def CPA4(A, B, c):
    
    """
    Complex Polar Algebra, method 4, Theory of Machines and Mechanisms
    The unknowns are A_angle, B_angle
    Using this case leads to 2 different results. Be sure to check out for
    appropriate conditions to come up with the right solution.
    """
    
    C, C_angle = np.abs(c), np.angle(c)
    pm = np.array([1, -1])
    B_angle = C_angle + pm*np.arccos((C**2+B**2-A**2)/(2*B*C))
    A_angle = C_angle + pm*np.arccos((C**2+A**2-B**2)/(2*A*C))
    results = []
    if A*np.sin(A_angle[0] - C_angle) == -B*np.sin(B_angle[0] - C_angle):
        results.append([A_angle[0], B_angle[0]])
        results.append([A_angle[1], B_angle[1]])
        return np.array(results[0]), np.array(results[1])
    else:
        results.append([A_angle[0], B_angle[1]])
        results.append([A_angle[1], B_angle[0]])
        return np.array(results[0]), np.array(results[1])


"Velocity Analysis"

def vel(r, omg=0, v_T=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms,
    returns total velocity (wrt to fixed frame of reference)
    
    Take first derivative of displacement vector:
        d
        -- (r) = v = R_T*exp(j*arg(r)) + j*omg*r
        dt
    
    It can also be written in matrix form as:
    [cos(angle(r)), -imag(r)] [ v_T ]
    [sin(angle(r)),  real(r)] [ omg ]
    
    The first row is real. The second row is imaginary.
    
    v_T: translational velocity
    omg: angular velocity
    r: displacement vector
    """

    return v_T*np.exp(1j*np.angle(r)) + 1j*omg*r

def vCPA1(v_known, r1, r2, omg2=0, v2_T=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms
    The unknowns are v1_T, omg1
    v_known + v_r1 = v_r2
    """
    
    A = [[-np.cos(np.angle(r1)),  np.imag(r1)],
         [-np.sin(np.angle(r1)), -np.real(r1)]]
    b = [np.real(v_known) - (np.cos(np.angle(r2))*v2_T - np.imag(r2)*omg2), 
         np.imag(v_known) - (np.sin(np.angle(r2))*v2_T + np.real(r2)*omg2)]
    return np.linalg.solve(A, b)

def vCPA2(v_known, r1, r2, omg1=0, v2_T=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms
    The unknowns are v1_T, omg2
    v_known + v_r1 = v_r2
    """
    
    A = [[-np.cos(np.angle(r1)), -np.imag(r2)],
         [-np.sin(np.angle(r1)),  np.real(r2)]]
    b = [np.real(v_known) - ( np.imag(r1)*omg1 + np.cos(np.angle(r2))*v2_T),
         np.imag(v_known) - (-np.real(r1)*omg1 + np.sin(np.angle(r2))*v2_T)]
    return np.linalg.solve(A, b)

def vCPA3(v_known, r1, r2, omg1=0, omg2=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms
    The unknowns are V1_T, V2_T
    v_known + v_r1 = v_r2
    """
    
    A = [[-np.cos(np.angle(r1)), np.cos(np.angle(r2))], 
         [-np.sin(np.angle(r1)), np.sin(np.angle(r2))]]
    b = [np.real(v_known) - ( np.imag(r1)*omg1 - np.imag(r2)*omg2),
         np.imag(v_known) - (-np.real(r1)*omg1 + np.real(r2)*omg2)]
    return np.linalg.solve(A, b)

def vCPA4(v_known, r1, r2, v1_T=0, v2_T=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms,
    The unknowns are omg1, omg2
    v_known + v_r1 = v_r2
    """
    
    A = [[ np.imag(r1), -np.imag(r2)],
         [-np.real(r1),  np.real(r2)]]
    b = [np.real(v_known) - (-np.cos(np.angle(r1))*v1_T + np.cos(np.angle(r2))*v2_T),
         np.imag(v_known) - ( np.sin(np.angle(r1))*v1_T + np.sin(np.angle(r2))*v2_T)]
    return np.linalg.solve(A, b)


"Acceleration Analysis"

"v_T: translational velocity"
"a_T: translational acceleration"

def coriolis(theta, omg=0, v_T=0): return 2j*v_T*omg*np.exp(1j*theta)

def centripetal(rX, omg=0): return -rX*omg**2

def acc(rX, omg=0, alp=0, v_T=0, a_T=0):
    
    """
    Complex-algebraic velocity analysis, Theory of Machines and Mechanisms,
    returns total acceleration (wrt to fixed frame of reference)
    
    Convetionally, there are 4 elements in an acceleration vector:
        - Centripetal acceleration: -omg**2*r
        - Tangential acceleration: alp*r
        - Coriolis acceleration: norm = |2*v_T*omg|
        - Relative acceleration: depends on many things
    
    Take first derivative of velocity vector:
        d         d
        -- (v) =  -- (R_T*exp(j*arg(r)) + j*omg*r)
        dt        dt
    
    Every acceleration vector on a plane can be written as
    [cos(angle(r)), -imag(r)] [ a_T ] _ [real(r)*omg**2 + 2*omg*v_T*sin(angle(r))]
    [sin(angle(r)),  real(r)] [ alp ]   [imag(r)*omg**2 - 2*omg*v_T*cos(angle(r))]
    
    
    The first row is real. The second row is imaginary.
    
    r: displacement vector
    omg: angular velocity
    v_T: translational velocity
    alp: angular acceleration
    a_T: translational acceleration
    """
    
    return a_T*np.exp(1j*np.angle(rX)) + 1j*alp*rX +\
        coriolis(np.angle(rX), omg, v_T) + centripetal(rX, omg)

def aCPA1(a_known, r1, r2, omg1=0, omg2=0, alp2=0, v1_T=0, v2_T=0, a2_T=0):
    
    """
    a_known + a_1 = a_2
    return a1_T, alp1
    """
    
    a1_known = centripetal(r1, omg1) + coriolis(np.angle(r1), omg1, v1_T)
    a2_known = centripetal(r2, omg2) + coriolis(np.angle(r2), omg2, v2_T)
    free_coeff = a_known + a1_known - a2_known
    
    A = [[-np.cos(np.angle(r1)),  np.imag(r1)],
         [-np.sin(np.angle(r1)), -np.real(r1)]]
    b = [np.real(free_coeff) - (np.cos(np.angle(r2))*a2_T - np.imag(r2)*alp2),
         np.imag(free_coeff) - (np.sin(np.angle(r2))*a2_T + np.real(r2)*alp2)]
    return np.linalg.solve(A, b)
    
def aCPA2(a_known, r1, r2, omg1=0, omg2=0, alp1=0, v1_T=0, v2_T=0, a2_T=0):
    
    """
    a_known + a_1 = a_2
    return a1_T, alp2
    """
    
    a1_known = centripetal(r1, omg1) + coriolis(np.angle(r1), omg1, v1_T)
    a2_known = centripetal(r2, omg2) + coriolis(np.angle(r2), omg2, v2_T)
    free_coeff = a_known + a1_known - a2_known
    
    A = [[-np.cos(np.angle(r1)), -np.imag(r2)],
         [-np.sin(np.angle(r1)),  np.real(r2)]]
    b = [np.real(free_coeff) - ( np.imag(r1)*alp1 + np.cos(np.angle(r2))*a2_T),
         np.imag(free_coeff) - (-np.real(r1)*alp1 + np.sin(np.angle(r2))*a2_T)]
    return np.linalg.solve(A, b)

def aCPA3(a_known, r1, r2, omg1=0, omg2=0, alp1=0, alp2=0, v1_T=0, v2_T=0):
    
    """
    a_known + a_1 = a_2
    return a1_T, a2_T
    """
    
    a1_known = centripetal(r1, omg1) + coriolis(np.angle(r1), omg1, v1_T)
    a2_known = centripetal(r2, omg2) + coriolis(np.angle(r2), omg2, v2_T)
    free_coeff = a_known + a1_known - a2_known
    
    A = [[-np.cos(np.angle(r1)), np.cos(np.angle(r2))],
         [-np.sin(np.angle(r1)), np.sin(np.angle(r2))]]
    b = [np.real(free_coeff) - ( np.imag(r1)*alp1 - np.imag(r2)*alp2),
         np.imag(free_coeff) - (-np.real(r1)*alp1 + np.real(r2)*alp2)]
    return np.linalg.solve(A, b)

def aCPA4(a_known, r1, r2, omg1=0, omg2=0, v1_T=0, v2_T=0, a1_T=0, a2_T=0):
    
    """
    a_known + a_1 = a_2
    return alp1, alp2
    """
    
    a1_known = centripetal(r1, omg1) + coriolis(np.angle(r1), omg1, v1_T)
    a2_known = centripetal(r2, omg2) + coriolis(np.angle(r2), omg2, v2_T)
    free_coeff = a_known + a1_known - a2_known
    
    A = [np.imag([r1, -r2]),
         np.real([-r1, r2])]
    b = [np.real(free_coeff) - (-np.cos(np.angle(r1))*a1_T + np.cos(np.angle(r2))*a2_T),
         np.imag(free_coeff) - (-np.sin(np.angle(r1))*a1_T + np.sin(np.angle(r2))*a2_T)]
    return np.linalg.solve(A, b)