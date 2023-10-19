#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"coeffs-mat.h"
int main() {
    int m = 2;  // You can set the dimensions of your matrix
    int n = 3;

// Create a matrix
	double **G_v = createMat(m, n);
	double **C_m= createMat(3,3);
	double **C_mid= createMat(3,3);
	double **R_o=createMat(2,2);
	double **C_mid_dir= createMat(3,3);

// Reading matrix data from .dat file
G_v=loadtxt("vert.dat",m,n);
C_m=loadtxt("C.dat",3,3);
C_mid=loadtxt("C_mid.dat",3,3);
C_mid_dir=loadtxt("C_mid_dir.dat",3,3);
R_o=loadtxt("R.dat",2,2);



	// **********************VECTORS*****************************
	double **G_dir=matmul(G_v,C_m,2,3,3);
	double **G_n=matmul(R_o,G_dir,2,2,3);
	double **G_con=diag(matmul(transpose(G_n,2,3),G_v,3,2,3),3,3);
	double **G_dis=sqrt_diag(matmul(transpose(G_dir,2,3),G_dir,3,2,3),3,3);
	double **G_line=h_concat(transpose(G_n,2,3),transpose(G_con,1,3),3,2,3,1);
	//***********************MEDIANS*******************************
	double **G_mid=nmatmul(0.5,matmul(G_v,C_mid,2,3,3),2,3);
	double **G_med_dir=matmul(G_v,C_mid_dir,2,3,3);
	double **G_n_med = matmul(R_o,G_med_dir,2,2,3);
	double **cmat_med=diag(matmul(transpose(G_n_med,2,3),G_v,3,2,3),3,3);
	double **linemat_med=h_concat(transpose(G_n_med,2,3),transpose(cmat_med,1,3),3,2,3,1);	
	double **G_G=line_intersect(linemat_med,3,3);

//printf("\n Vectors \n");
//printf("direction matrix= \n");
//print(G_dir,2,3);
//printf("normal matrix= \n");
//print(G_n,2,3);
//printf("constant matrix= \n");
//print(G_con,1,3);
//printf("distance matrix= \n");
//print(G_dis,1,3);
//printf("line  matrix= \n");
//print(G_line,3,3);

printf("\n Medians \n");
printf("midpoint matrix= \n");
print(G_mid,2,3);
printf("median direction  matrix= \n");
print(G_med_dir,2,3);
printf("median normal matrix= \n");
print(G_n_med,2,3);
printf("median constant matrix= \n");
print(cmat_med,1,3);
printf("median line matrix= \n");
print(linemat_med,3,3);
printf("Centroid = \n");
print(G_G,1,2);


}



