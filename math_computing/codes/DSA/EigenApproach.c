#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "listgen.h"
#include "listfun.h"
int  main(){
//******************Eigen value approach to find contact points********************************
//load matrix from file
avyuh *h= loadList("A.dat", 2, 1);		    // vertex A
avyuh *u=loadList("B.dat",2,1);			    // Incenter I
double r= 2; 			    // incircle radius
u->vector->data=-u->vector->data;	u->next->vector->data=-u->next->vector->data;  // incentre u=-I  
avyuh *V=Listeye(2);			  			 	// 2x2 identity matrix for circle
avyuh *f=createList(1,1); f->vector->data=pow(Listnorm(u),2)-r*r;	//     from circle equatin x^2+y^2=f
//
avyuh *contactpts=ListContactPts(V,h,u,f);     // function to find contact points 
//
// Printing lists . 	
printf("contact points = \n");	printList(contactpts);
return 0;
}
