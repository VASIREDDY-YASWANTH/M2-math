#Code by GVV Sharma
#September 7, 2023
#Revised October 1, 2023
#Revised October 2, 2023
#released under GNU GPL
#Matrix Algebra


import sys                                          #for path to external scripts
sys.path.insert(0, '/home/VASIREDDY_YASWANTH/github/geometry/codes/CoordGeo')        #path to my scripts
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
df= pd.read_excel('tables/vertices.xlsx')
print(df ,"\n")

#Triangle Vertices
G_v= df.to_numpy()[:,:]
print(G_v, C_m ,"\n")

#Direction vector circulant matrix
C_m= SA.circulant([1,0,-1]).T


#Direction vector Matrix
G_dir = G_v@C_m
print(G_dir,"\n")

#Normal vector matrix
G_n = R_o@G_v@C_m
print(G_n,"\n")

#Find the line constants
cmat = np.diag(G_n.T@G_v).reshape(-1,1)

#line matrix
linmat = np.block([G_n.T,cmat])
print(linmat,"\n")

#sides vector
dis = np.linalg.norm(G_dir, axis=0).reshape(-1,1)
print(dis,"\n")

'''
#Finding the angles of the triangle
dmat = np.diag(1/d)
G_dnorm = G_dir@dmat
G_dgram = G_dnorm.T@G_dnorm
#print(np.degrees(np.arccos(G_dgram)))
'''
print("vector ends")
#-----------------Vectors Ends-------------------------------

#-----------------Medians-------------------------------
#Median circulant matrix
C_mid = SA.circulant([0,1,1]).T
print(C_mid,"\n")

#Mid point matrix
G_mid = 0.5*G_v@C_mid
print(G_mid,"\n")

#Median direction circulant matrix
C_mid_dir = SA.circulant([1,-0.5,-0.5])

#Median direction matrix
G_med_dir = G_v @ C_mid_dir
print(G_med_dir,"\n")

#Normal vector matrix
G_n_med = R_o@G_med_dir

#Find the line constants
cmat_med = np.diag(G_n_med.T@G_v).reshape(-1,1)

#median  matrix
linmat_med = np.block([G_n_med.T,cmat_med])
print(linmat_med ,"\n")

#Find the centroid
#print(LA.lstsq(G_n_med.T,cmat_med),"\n")
print("median ends")
#-----------------Median Ends-------------------------------

#-----------------Altitude-------------------------------

#Circulant matrix
C_alt= SA.circulant([0,-1,1]).T
print(C_alt ,"\n")

#Normal Matrix
G_dir_alt = G_v@C_alt
print(G_dir_alt,"\n")

#Find the line constants
cmat_alt = np.diag(G_dir_alt.T@G_v).reshape(-1,1)
#print( np.block([G_dir_alt.T,cmat_alt]),"\n")

#altitude matrix
linmat_alt= np.block([G_dir_alt.T,cmat_alt])
print(linmat_alt,"\n")

#Find the orthocentre
print(LA.lstsq(G_dir_alt.T,cmat_alt))
print("altitude ends")
#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
#Find the line constants
cmat_perp_bis= np.diag(G_dir_alt.T@G_mid).reshape(-1,1)
print( np.block([G_dir_alt.T,cmat_perp_bis]))


#Find the Circumcentre
#print("\n",A.lstsq(G_dir_alt.T,cmat_perp_bis))
#-----------------Altitude Ends-------------------------------

#-----------------Perpendicular Bisector-------------------------------
#Incircle circulant matrix
C_in = SA.circulant([1,1,0]).T
print(C_in,"\n")
#m,n,p
secvec = LA.inv(C_in)@dis
cont_mat =np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])])
print(cont_mat,"\n")
#np.block(np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]]))
print(np.array([np.block([secvec[1]/dis[1],secvec[0]/dis[2],0]), np.block([0, secvec[2]/dis[2],secvec[1]/dis[0]]),np.block([secvec[2]/dis[1],0,secvec[0]/dis[0]])]))
cont_mat = np.array([secvec[1]/dis[1],secvec[0]/dis[2],0],dtype=object)
print("\n",cont_mat)
print("\n",secvec[1]/dis[0])
print("\n",C_in,"\n","\n",C_in.T)
tvec = np.array([1,2,3]).reshape(-1,1)
tC = SA.circulant([0,1,0])
print("\n",tC,"\n",tvec)
#print(tC@tvec)

