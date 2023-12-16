#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"matfun.h"
int main() {
	double **h= createMat(2,1);
	double **u= createMat(2,1);
h=loadMat("A.dat",2,1); //vertex A
u=loadMat("B.dat",2,1); //centre of circle B

	double r=2;  // radius of circle 
        double **V=Mateye(2);  //2x2 identity matrix
        u[0][0]=-u[0][0];    u[1][0]=-u[1][0];
        double **f=createMat(1,1); f[0][0]=pow(Matnorm(u,2),2)-r*r;


double **contactPts=contactPoints(V,h,u,f);

// Printing matrices . 	
printf("contactpts = \n");
printMat(contactPts,2,2);
}

