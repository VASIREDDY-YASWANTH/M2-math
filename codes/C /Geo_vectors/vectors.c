#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"coeffs-mat.h"
int main() {
    int m = 2;  // You can set the dimensions of your matrix
    int n = 3;

    // Create a matrix
    double **A,**R_o,**C_m,**G_dir,**G_n,**G_con,**G_dis,**G_line;
	A = createMat(m, n);
	C_m= createMat(3,3);
	R_o=createMat(2,2);

A=loadtxt("vert.dat",m,n);
C_m=loadtxt("C.dat",3,3);
R_o=loadtxt("R.dat",2,2);

//print(A,2,3);

	G_dir=matmul(A,C_m,2,3,3);
	G_n=matmul(R_o,G_dir,2,2,3);
	G_con=diag(matmul(transpose(G_n,2,3),A,3,2,3),3,3);
	G_dis=sqrt_diag(matmul(transpose(G_dir,2,3),G_dir,3,2,3),3,3);
	G_line=h_concat(transpose(G_n,2,3),transpose(G_con,1,3),3,2,3,1);
	//G_con=diag(G_c,3,3);


printf("direction matrix= \n");
print(G_dir,2,3);
printf("normal matrix= \n");
print(G_n,2,3);
printf("constant matrix= \n");
print(G_con,1,3);
printf("distance matrix= \n");
print(G_dis,1,3);
printf("line  matrix= \n");
print(G_line,3,3);
}




