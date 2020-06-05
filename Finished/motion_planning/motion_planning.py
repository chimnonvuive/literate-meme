#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import modern_robotics as mr
import numpy as np
from scipy.spatial.distance import euclidean

# import csv file as an array
ARRAY = np.genfromtxt('obstacles.csv', delimiter=',')

# generate coordinates from array in the form of [x,y,radius]
CIRCLES = np.array([[x,y,d/2] for element in ARRAY for x,y,d in [element]])

# set arrays of xs, ys, radii individually
XS, YS, RS = CIRCLES.T


# In[534]:


def IsCollisionFree(A, B, h):
    
    # find MN perp AB, where AB is a finite line, M is center, N in BA
    BA = B - A
    matA = [[BA[0], BA[1]], [BA[1], -BA[0]]]
    matB = [XS[h]*BA[0] + YS[h]*BA[1],
            B[0]*BA[1] - B[1]*BA[0]]
    
    # N in BA, MN perp AB
    N = np.linalg.solve(matA, matB)
    
    # center of circle #h
    M = (XS[h],YS[h])
    
    # if MB < r or MA < r or (AB=AN+BN & MN < r) return False (proved by geometry)
    if euclidean(M, B) < RS[h] or       euclidean(M, A) < RS[h] or       (euclidean(A,B) == (euclidean(A,N) + euclidean(B,N)) and
           euclidean(N, M) < RS[h]):
        
        return False
    
    else:
        
        if h == 0:
            return True
        else:
            return IsCollisionFree(A, B, h-1)

        
def x_newPosition(x_nearest, x_samp, d):
    
    # assign variables to the coordinates
    xa, ya = x_nearest
    xc, yc = x_samp
    
    """
    Solving the equations symbolically yields 2 results of x_new
    
    eq1 = (ya-yb)/(xa-xb) - (ya-yc)/(xa-xc)
    eq2 = (xa-xb)**2 + (ya-yb)**2 - d**2
    
    In this case I use sympy and copy the results. The code is as follows:
    
    from sympy import symbols, solve, simplify
    xa, xb, xc, ya, yb, yc, d = symbols('xa xb xc ya yb yc d')
    eq1 = (ya-yb)/(xa-xb) - (ya-yc)/(xa-xc)
    eq2 = (xa-xb)**2 + (ya-yb)**2 - d**2
    xb1, yb1 = solve([eq1, eq2], [xb, yb])[0]
    xb2, yb2 = solve([eq1, eq2], [xb, yb])[1]
    """
    
    xb1 = (-d*xa + d*xc + xa*np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2))/np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2)
    yb1 = -d*(ya - yc)/np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2) + ya
    xb2 = ( d*xa - d*xc + xa*np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2))/np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2)
    yb2 =  d*(ya - yc)/np.sqrt(xa**2 - 2*xa*xc + xc**2 + ya**2 - 2*ya*yc + yc**2) + ya
    x_new_1 = np.array([xb1, yb1])
    x_new_2 = np.array([xb2, yb2])
    
    # return the correct x_new if ab + bc = ac algebraically
    if (euclidean(x_nearest, x_new_1) + euclidean(x_new_1, x_samp)) == euclidean(x_nearest, x_samp):
        return x_new_1
    else:
        return x_new_2


def Path(edges):
    start, goal = edges[0][0], edges[-1][1]
    path = [goal]
    tmp_edges = edges[::-1]
    i, n = 0, len(edges)
    while i<n-2:
        for j in range(i+1, n):
            if np.array_equal(tmp_edges[i][0], tmp_edges[j][1]):
                path.append(tmp_edges[j][1])
                i = j
                break
    path.append(start)
    return path[::-1]
    
    
def RRT(x_start, goal, r_goal, nodes_max, d):
    
    "Plot RRT on a graph"
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_aspect('equal')
    
    # plot the circles (obstacles)
    circles_plot = []
    for x, y, r in zip(XS, YS, RS):
        circles_plot.append(plt.Circle((x,y), r, color='blue'))
    for circle in circles_plot:
        ax.add_artist(circle)
    
    "begin the algorithm"
    tree_size = 0
    
    # initialize T with x_start
    nodes = [x_start]
    edges = []
    
    # generate random point (x, y) in the interval [-0.5, 0.5] with 10% bias around the goal
    # you can change the bias to see the differences
    bias = 0.1
    x_uniform = np.random.rand(int(nodes_max*(1-bias)), 2) - 0.5
    x_biased = np.random.rand(int(nodes_max*bias), 2)*2*r_goal + (0.5 - 2*r_goal)
    
    # connect 2 arrays together as a new array and then shuffle it
    x_samps = np.concatenate((x_uniform, x_biased), axis=0)
    np.random.shuffle(x_samps)

    while tree_size < nodes_max:

        # compute distances of all points relative to rand_pt
        all_distances = np.linalg.norm(np.asarray(nodes) - x_samps[tree_size], axis=1)

        # choose the point with smallest distance to rand_pt based on index of T
        smallest_distance = all_distances.argmin()
        x_nearest = nodes[smallest_distance]
        x_new = x_newPosition(x_nearest, x_samps[tree_size], d)
        
        if IsCollisionFree(x_new, x_nearest, len(XS)-1):
    
            nodes.append(x_new)
            edges.append((x_nearest, x_new))
            ax.plot([x_nearest[0], x_new[0]], [x_nearest[1], x_new[1]], color='red', linewidth=0.5)
            
            if euclidean(x_new, goal) < r_goal:
                
                # save figure as png file and display it along with the path
                path = np.asarray(Path(edges))
                ax.plot(path[:,0], path[:,1], color='green')
                plt.savefig('RRT_algorithm.png', dpi=600)
                plt.show()
                return "SUCCESS", np.asarray(nodes), np.asarray(edges), path
        
        tree_size += 1
    
    return "FAILURE", np.asarray(nodes), np.asarray(edges)

# declare inputs
nodes_max = 1500
x_start = np.array([-0.5, -0.5])
d = 0.025
goal = np.array([0.5, 0.5])
r_goal = 0.05
isSuccess, nodes, edges, path = RRT(x_start, goal, r_goal, nodes_max, d)


# In[554]:


np.savetxt("nodes.csv", nodes, fmt='%.5f', delimiter=',')
np.savetxt("edges.csv", np.c_[edges[:,0,:],edges[:,1,:]], header='x1, y1, x2, y2', fmt=['%.5f']*4, delimiter=',')
np.savetxt("path.csv", path, fmt='%.5f', delimiter=',')

