import numpy as np
from sympy import init_printing
init_printing(use_unicode=True)
from sympy.matrices import Matrix

import sys                            
sys.path.insert(0, '/home/gadepall/github/geometry/codes/CoordGeo')        #path to my scripts
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
import numpy.linalg as LA

def line_intersect(n1,A1,n2,A2):
  N=np.vstack((n1,n2))
  p = np.zeros(2)
  p[0] = n1@A1
  p[1] = n2@A2
  #Intersection
  P=np.linalg.inv(N)@p
  return P

A=np.array([-3,-5])
B=np.array([3,-5])
C=np.array([-4,-3])

# Triangle midpoints
D = (B + C) / 2
E = (C + A) / 2
F = (A + B) / 2

O = line_intersect(AB,F,AC,E)
G=(B+C)/2
#O=np.array([0,-2.25])
X = A - O
radius = np.linalg.norm(X)
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_OA = line_gen(O,A)
x_OG=line_gen(O,G)
X_OE=line_gen(O,E)
X_OF=line_gen(O,F)


#Generating the circumcirclecircle
[O,r] = ccircle(A,B,C)
x_ccirc= circ_gen(O,radius)
#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
#Plotting line OG
#plt.plot(x_OG[0,:],x_OG[1,:],label='$OG$')
#Plotting lines OE and OF
#plt.plot(x_OE[0,:],x_OE[1,:],label='$OE$')
#plt.plot(x_OF[0,:],x_OF[1,:],label='$OF$')

#Plotting the circumcircle
plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')


A = A.reshape(-1,1)
B = B.reshape(-1,1)
C = C.reshape(-1,1)
O = O.reshape(-1,1)
G = G.reshape(-1,1)
E = E.reshape(-1,1)
F = F.reshape(-1,1)

#Labeling the coordinates
tri_coords = np.block([[A,B,C,O,G,E,F]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','O','G','E','F']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() 
plt.axis('equal')
plt.savefig('figs/triangle/perp_bisect5.png')
