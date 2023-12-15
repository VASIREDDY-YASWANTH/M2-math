
import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')        #path to my scripts
import numpy as np
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd
import math as m

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen

#if using termux
import subprocess
import shlex
#end if

A = np.array([0,0]).reshape(-1,1)
B = np.array([5,0]).reshape(-1,1) 
r=4                 #length of AP and AQ
theta=m.pi/6           # angle QAB and angle PAB  30 degrees
P = np.array([r*m.cos(theta),-r*m.sin(theta)]).reshape(-1,1)
Q = np.array([r*m.cos(theta),r*m.sin(theta)]).reshape(-1,1)

print(P,"\n",Q,"\n")

a = LA.norm(A-B)    #length of AB
x1 = LA.norm(B-P)    #length of BP
x2 = LA.norm(B-Q)    #length of BQ
#print(x1,"\n",x2,"\n")
print("1.By the definition of AAS congruency rule if any two pair of angles and one pair of corresponding sides are equal then the two triangles are said to be congruent \n\n")


if(x1==x2):
    print("2.BP equals BQ")
else:
    print("2.BP not equals to BQ")


#Generating all lines
x_AB = line_gen(A,B)
x_AQ = line_gen(A,Q)
x_AP = line_gen(A,P)
x_BQ = line_gen(B,Q)
x_BP = line_gen(B,P)


#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_AQ[0,:],x_AQ[1,:],label='$AQ$')
plt.plot(x_AP[0,:],x_AP[1,:],label='$AP$')
plt.plot(x_BQ[0,:],x_BQ[1,:],label='$BQ$')
plt.plot(x_BP[0,:],x_BP[1,:],label='$BP$')

#Labeling the coordinates
tri_coords = np.block([[A,B,P,Q]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','P','Q']
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

#if using termux
plt.savefig('fig_mat_comp.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()

