import sys                                          #for path to external scripts
sys.path.insert(0, '/home/gadepall/github/geometry/codes/CoordGeo')        #path to my scripts
import numpy as np
import mpmath as mp
import numpy.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen
#if using termux
import subprocess
import shlex
#end if

A = np.array([-3, -5])  # Use floating-point numbers
B = np.array([3, -5])  # Use floating-point numbers
C = np.array([-4, -3])  # Use floating-point numbers

# Midpoint of each line
def midpoint(P, Q):
    return (P + Q) / 2  
#normal vector 
def norm_vec(A,B):
  omat = np.array([[0,1],[-1,0]]) 
  return omat.T@(A-B)
#to find the coefficients and constant of the equation of perpendicular bisector of BC
def perpendicular_bisector(B, C):
    midBC=midpoint(B,C)
    dir=B-C
    constant = -dir.T @ midBC
    return dir,constant
equation_coeff1,const1 = perpendicular_bisector(A, B)
equation_coeff2,const2 = perpendicular_bisector(B, C)
equation_coeff3,const3 = perpendicular_bisector(C, A)
print(f'Equation for perpendicular bisector of AB: ({equation_coeff1[0]:.2f})x + ({equation_coeff1[1]:.2f})y + ({const1:.2f}) = 0')
print(f'Equation for perpendicular bisector of  BC: ({equation_coeff2[0]:.2f})x + ({equation_coeff2[1]:.2f})y + ({const2:.2f}) = 0')
print(f'Equation for perpendicular bisector of  CA: ({equation_coeff3[0]:.2f})x + ({equation_coeff3[1]:.2f})y + ({const3:.2f}) = 0')
#circumcentre of triangle ABC
def ccircle(A,B,C):
  p = np.zeros(2)
  n1 = equation_coeff1[:2]
  p[0] = 0.5*(np.linalg.norm(A)**2-np.linalg.norm(B)**2)
  n2 = equation_coeff2[:2]
  p[1] = 0.5*(np.linalg.norm(B)**2-np.linalg.norm(C)**2)
  #Intersection
  N=np.block([[n1],[n2]])
  O=np.linalg.solve(N,p)
  return O
O=ccircle(A,B,C)
# Generate points along a line
def line_gen(A, B):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(0, 1, len)
    for i in range(len):
        temp1 = A + lam_1[i] * (B - A)
        x_AB[:, i] = temp1.T
    return x_AB
# Perpendicular bisector
def line_dir_pt(m, A, k1=0, k2=1):
    len = 10
    dim = A.shape[0]
    x_AB = np.zeros((dim, len))
    lam_1 = np.linspace(k1, k2, len)
    for i in range(len):
        temp1 = A + lam_1[i] * m
        x_AB[:, i] = temp1.T
    return x_AB
# Calculate the perpendicular vector and plot arrows
def perpendicular(B, C, label):
    perpendicular=norm_vec(B,C)
    mid = midpoint(B, C)
    x_D = line_dir_pt(perpendicular, mid, 0, 1)
    plt.arrow(mid[0], mid[1], perpendicular[0], perpendicular[1], color='blue', head_width=0.4, head_length=0.4, label=label)
    plt.arrow(mid[0], mid[1], -perpendicular[0], -perpendicular[1], color='blue', head_width=0.4, head_length=0.4)
    return x_D

#Perpendicular Bisector
#x_D = perpendicular(A, B, 'OD')
#x_E = perpendicular(B, C, 'OE')
#x_F = perpendicular(C, A, 'OF')


from triangle.funcs import *

  
def dir_vec(A,B):
  return B-A
def midpoint(P, Q):
    return (P + Q) / 2  

#Triangle vertices
A = A.reshape(-1,1)
B = B.reshape(-1,1) 
C = C.reshape(-1,1) 

D,E,F = tri_mid_pt(A,B,C)
G=(C+B)/2
# direction vector along line joining A & B
AB = dir_vec(A,B)
# direction vector along line joining A & C
AC = dir_vec(A,C)
#O = line_intersect(AB,F,AC,E)
print(O)

#Generating the circumcircle
[O,R] = ccircle(A,B,C)
x_circ= circ_gen(O,R)

#angle BOC = 2 angle BAC
m1 = dir_vec(A,B)
m2 = dir_vec(A,C)
m3 = dir_vec(B,C)
v1 = dir_vec(O,B)
v2 = dir_vec(O,C)

print(ang_vec(m1,m2),2*mp.pi-ang_vec(v1,v2))

#OD perpendicular to BC
m3 = dir_vec(B,C)
print(m3.T@(O-D))


#Plotting the circumcircle
plt.plot(x_circ[0,:],x_circ[1,:],label='$circumcircle$')
#Generating all lines
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A)
x_OD = line_gen(O,D)
x_OE = line_gen(E,O)
x_OF = line_gen(O,F)
x_OG = line_gen(O,G)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:])
plt.plot(x_BC[0,:],x_BC[1,:])
plt.plot(x_CA[0,:],x_CA[1,:])
#plt.plot(x_OD[0,:],x_OD[1,:],'--',label='$OD$')
#plt.plot(x_OE[0,:],x_OE[1,:],'--',label='$OE$')
#plt.plot(x_OF[0,:],x_OF[1,:],'--',label='$OF$')
#plt.plot(x_OG[0,:],x_OG[1,:],label='$OG$')

#Labeling the coordinates
tri_coords = np.block([[A,B,C,D,E,F,O]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D','E','F','O']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')
plt.savefig('figs/triangle/perp-bisect.png')
plt.show()
