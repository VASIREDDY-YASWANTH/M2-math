import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')        #path to my scripts
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import *

#if using termux
import subprocess
import shlex
#end if



#Triangle vertices
A = np.array([0,0]).reshape(-1,1)
B = np.array([3,3]).reshape(-1,1) 

#Incircle parameters
h  = A
V = np.eye(2)
r=2
u = -B
f = LA.norm(B)**2-r**2

#contact points
[Q,P] = contact(V,u,f,h)

#printing CONTACT POINTS
print("P=",P)
print("Q=",Q)

#Generating all lines
x_AP = line_gen(A,P)
x_AQ = line_gen(A,Q)
x_BP = line_gen(B,P)
x_BQ = line_gen(B,Q)
x_AB = line_gen(A,B)

#generating incircle
x_circ = circ_gen(B,2)

#Plotting all lines
plt.plot(x_AP[0,:],x_AP[1,:],label='$AP$')
plt.plot(x_AQ[0,:],x_AQ[1,:],label='$AQ$')
plt.plot(x_BP[0,:],x_BP[1,:],label='$BP$')
plt.plot(x_BQ[0,:],x_BQ[1,:],label='$BQ$')
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_circ[0,:],x_circ[1,:],label='$incircle$')

#Labeling the coordinates
A = A.reshape(-1,1)
B = B.reshape(-1,1)
P = P.reshape(-1,1)
Q = Q.reshape(-1,1)

tri_coords = np.block([[A,B,P,Q]])
plt.scatter(tri_coords[0, :], tri_coords[1, :])
vert_labels = ['A','B','P','Q']
for i, txt in enumerate(vert_labels):
    offset = 10 if txt == 'F' else -10
    plt.annotate(txt,
                 (tri_coords[0, i], tri_coords[1, i]),
                 textcoords="offset points",
                 xytext=(0, offset),
                 ha='center')
plt.xlabel('$X-axis$')
plt.ylabel('$Y-axis$')
plt.legend(loc='upper left')
plt.grid() # minor
plt.axis('equal')
plt.savefig('eigen.png')
plt.show()
