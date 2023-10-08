#Code by GVV Sharma
#September 7, 2023
#Revised October 1, 2023
#Revised October 2, 2023
#released under GNU GPL
#Matrix Algebra


import sys                                          #for path to external scripts
sys.path.insert(0, './CoordGeo')        #path to my scripts
import numpy as np
import numpy.linalg as LA
import scipy.linalg as SA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen


#if using termux
import subprocess
import shlex
#end if



#-----------------Vectors-------------------------------

#Input parameters from excel file
df= pd.read_excel('vertices.xlsx')
#print(df ,"\n")

#Triangle Vertices
G_v= df.to_numpy()[:,:]

#Direction vector circulant matrix
C_m= SA.circulant([1,0,-1]).T
#print(G_v, C_m ,"\n")


#Direction vector Matrix
G_dir = G_v@C_m
#print(G_dir,"\n")

#Normal vector matrix
G_n = R_o@G_v@C_m
#print(G_n,"\n")

#Find the line constants
cmat = np.diag(G_n.T@G_v).reshape(-1,1)

#line matrix
linmat = np.block([G_n.T,cmat])
#print(linmat,"\n")

#sides vector
dis = np.linalg.norm(G_dir, axis=0).reshape(-1,1)
#print(dis,"\n")

'''
#Finding the angles of the triangle
dmat = np.diag(1/d)
G_dnorm = G_dir@dmat
G_dgram = G_dnorm.T@G_dnorm
#print(np.degrees(np.arccos(G_dgram)))
'''
#print("vector ends")
#-----------------Vectors Ends-------------------------------

#-----------------Medians-------------------------------
#Median circulant matrix
C_mid = SA.circulant([0,1,1]).T
#print(C_mid,"\n")

#Mid point matrix
G_mid = 0.5*G_v@C_mid
#print(G_mid,"\n")

#Median direction circulant matrix
C_mid_dir = SA.circulant([1,-0.5,-0.5])

#Median direction matrix
G_med_dir = G_v @ C_mid_dir
#print(G_med_dir,"\n")

#Normal vector matrix
G_n_med = R_o@G_med_dir

#Find the line constants
cmat_med = np.diag(G_n_med.T@G_v).reshape(-1,1)

#median  matrix
linmat_med = np.block([G_n_med.T,cmat_med])
#print(linmat_med ,"\n")

#Find the centroid
G_G=LA.lstsq(G_n_med.T,cmat_med)
#print("G_G=",G_G)
#print(LA.lstsq(G_n_med.T,cmat_med),"\n")
#print("median ends")
#-----------------Median Ends-------------------------------

#-----------------Altitude-------------------------------

#Circulant matrix
C_alt= SA.circulant([0,-1,1]).T
#print(C_alt ,"\n")

#Normal Matrix
G_dir_alt = G_v@C_alt
#print(G_dir_alt,"\n")

#Find the line constants
cmat_alt = np.diag(G_dir_alt.T@G_v).reshape(-1,1)
#print( np.block([G_dir_alt.T,cmat_alt]),"\n")

#altitude matrix
linmat_alt= np.block([G_dir_alt.T,cmat_alt])
#print(linmat_alt,"\n")

#Find the orthocentre
G_H=LA.lstsq(G_dir_alt.T,cmat_alt)
print(LA.lstsq(G_dir_alt.T,cmat_alt))
#print("altitude ends")


#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
#Find the line constants
cmat_perp_bis= np.diag(G_dir_alt.T@G_mid).reshape(-1,1)
#print( np.block([G_dir_alt.T,cmat_perp_bis]))


#Find the Circumcentre
#G_O = A.lstsq(G_dir_alt.T,cmat_perp_bis)
G_O = LA.lstsq(G_dir_alt.T, cmat_perp_bis, rcond=None)
#print("G_O",G_O)
#print("\n",A.lstsq(G_dir_alt.T,cmat_perp_bis))
#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
#Incircle circulant matrix
C_in = SA.circulant([1,1,0]).T
#print("c_in=",C_in,"\n")
#m,n,p
secvec = LA.inv(C_in)@dis
cont_mat =np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])])
G_incir = G_v @ cont_mat
#print("incircle=" ,G_incir)
#print(cont_mat,"\n")
#np.block(np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]]))
#print(np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])]))
cont_mat = np.array([secvec[1]/dis[1],secvec[0]/dis[2],0],dtype=object)
#print("\n",cont_mat)
#print("\n",secvec[1]/dis[0])
#print("\n",C_in,"\n","\n",C_in.T)
tvec = np.array([1,2,3]).reshape(-1,1)
tC = SA.circulant([0,1,0])
#print("\n",tC,"\n",tvec)
#print(tC@tvec)





A = G_v[:, 0]
B = G_v[:, 1]
C = G_v[:, 2]
D = G_mid[:, 0]
E = G_mid[:, 1]
F = G_mid[:, 2]
G = G_G[0]
H = G_H[0]
O = G_O[0]
D3=G_incir[:,0]
E3=G_incir[:,1]
F3=G_incir[:,2]
print("A=",A,"B=",B,"C=",C,"D=",D,"E=",E,"F=",F,"G=",G,"H=",H,"O=",O,"D3=",D3,"E3=",E3,"F3=",F3)


A=A.reshape(-1,1)
B=B.reshape(-1,1)
C=C.reshape(-1,1)
D=D.reshape(-1,1)
E=E.reshape(-1,1)
F=F.reshape(-1,1)
G=G.reshape(-1,1)
H=H.reshape(-1,1)
D3=D3.reshape(-1,1)
E3=E3.reshape(-1,1)
F3=F3.reshape(-1,1)
print("A=",A,"B=",B,"C=",C,"D=",D,"E=",E,"F=",F,"G=",G,"H=",H,"O=",O,"D3=",D3,"E3=",E3,"F3=",F3)

radius = np.linalg.norm(A-O)
#print("r=",r)

#Generating the circumcirclecircle
[O,r] = ccircle(A,B,C)
x_ccirc= circ_gen(O,radius)


#Generating all lines 
x_AB = line_gen(A,B)
x_BC = line_gen(B,C)
x_CA = line_gen(C,A) 
x_AD = line_gen(A,D)
x_BE = line_gen(B,E)
x_CF = line_gen(C,F)
x_OA = line_gen(O,A)




#Plotting all lines
plt.plot(x_AB[0,:],x_AB[1,:],label='$AB$')
plt.plot(x_BC[0,:],x_BC[1,:],label='$BC$')
plt.plot(x_CA[0,:],x_CA[1,:],label='$CA$')
plt.plot(x_AD[0,:],x_AD[1,:],label='$AD$')
#plt.plot(x_BE[0,:],x_BE[1,:],label='$BE$')
#plt.plot(x_CF[0,:],x_CF[1,:],label='$CF$')
#plt.plot(x_OA[0,:],x_OA[1,:],label='$OA$')
#Plotting the circumcircle
#plt.plot(x_ccirc[0,:],x_ccirc[1,:],label='$circumcircle$')



#Labeling the coordinates
tri_coords = np.block([[A,B,C,D]])
plt.scatter(tri_coords[0,:], tri_coords[1,:])
vert_labels = ['A','B','C','D',]
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
plt.savefig('mat_alt1.png')
#subprocess.run(shlex.split("termux-open ./figs/tri_sss.pdf"))
#else
plt.show()

