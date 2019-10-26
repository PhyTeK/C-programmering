static char RCSindx[]= "$Id: bio.c,v 1.2 2000/01/11 17:27:05 root Exp root $";

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr.h"
#include "nrutil.h"
#include "lasar.h"

/* 
   Lasar method applied to a bioremediation case
   P. Martinet Jan 2000

   Compile:
   -------

   gcc -o bio bio.c -lm -lrecipes_c

   Revisions:
   ---------

   $Log: bio.c,v $
   Revision 1.2  2000/01/11 17:27:05  root
   Gaussian quadrature

   Revision 1.1  2000/01/10 18:00:22  root
   Initial revision


*/

static float U = 2;    /* mean velocity (m/day) */
static float s2 = .01; /* s2=sigma^2 logconductivity variance */
static float Iyh = 3;  /* horizontal log condictivity integral scale */
static float e = 0.1;  /* anisotropy ratio */
static float dtau=0.1;
static float A=100;    /* Upper limit of integration domain */
float TX;              /* global time variable */

/* Guassian Quadrature weights */
static float x[]={0.0,0.1488743389,0.4333953941,
		  0.6794095682,0.8650633666,0.9739065285};
static float w[]={0.0,0.2955242247,0.2692667193,
		  0.2190863625,0.1494513491,0.0666713443};

float bfun() 
{
  float t1,t2,t3;
  t1=e*e;
  t2=sqrt(1-t1);
  t3=(t1-1)*(t1-1);

  return 1+(19*t1-10*t1*t1)/16/t3-e*(13-4*t1)*asin(t2)/16/t2/t3;

}

float X11(float tp)
{
  return 2/bfun()*Iyh*Iyh*s2*(tp*bfun()+exp(-tp*bfun())-1);
}

float D0(float t)
{
  return Iyh*Iyh*s2*(1-exp(-t*bfun()));  /* dX11/dt */
}

float gfun(float tau, float x)
{
  float t0,t1,t2,t3,X11t;

  t0=tau*U/Iyh;  /* t=tau */
  t1=x-U*tau;
  t2=t1*t1;
  X11t = X11(t0);
  t3=(U*X11t + t1*D0(t0))/pow(X11t,1.5);

  return exp(-t2/2/X11t)*t3/sqrt(2*M_PI);
}

float ak(int k,int l,float t,float s)
{
  /* Kernel */

  float C,dtau2;
  dtau2 = dtau*dtau;      
  C = 1.0/sqrt(2*M_PI*dtau2);
  
  return  exp(-(t-s)*(t-s)/2/dtau2)*C;
}



void voltra1(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(float,float), float (*ak)(int, int, float, float))
{
	void lubksb(float **a, int n, int *indx, float b[]);
	void ludcmp(float **a, int n, int *indx, float *d);
	int i,j,k,l,*indx;
	float d,sum,**a,*b;

	indx=ivector(1,m);
	a=matrix(1,m,1,m);
	b=vector(1,m);
	t[1]=t0;
 	for (k=1;k<=m;k++) f[k][1]= -(*g)(TAU,t[1]);
	for (i=2;i<=n;i++) {
		t[i]=t[i-1]+h;
		for (k=1;k<=m;k++) {
			sum= -(*g)(TAU,t[i]);
			for (l=1;l<=m;l++) {
				sum += 0.5*h*(*ak)(k,l,t[i],t[1])*f[l][1];
				for (j=2;j<=n;j++)  /* replace j<i with j<=n */
					sum += h*(*ak)(k,l,t[i],t[j])*f[l][j];
				a[k][l]= -0.5*h*(*ak)(k,l,t[i],t[i]);
			}
			b[k]=sum;
		}
		ludcmp(a,m,indx,&d);
		lubksb(a,m,indx,b);
		for (k=1;k<=m;k++) f[k][i]=b[k];
	}
	free_vector(b,1,m);
	free_matrix(a,1,m,1,m);
	free_ivector(indx,1,m);
}

float intcg(float t)
{
  /*
    Return the integral of Cn(t,tau) gn(tau;x1) dtau 
    between 0 and A=100 using Gaussian quadrature at the points:
    tau_i (i=1..10) defined in qgaus.
    The result is the expected value <Cn(t;x1)> that should be equal to <C(t;x1)>.

    1) The Cn(t,tau_i) are obtained from ODE. The results is stored in the matrix
    CN[1:10][1:NSTEP].

    2) We form the product Cn(t,tau_i) gn(tau_i;x1) in the matrix CNG[1:10][1:NSTEP]

    3) We applied Gaussian quadrature for each distribution CNG[k][1:NSTEP]

    4) Return the result CNG[k][l] that corresponds with the time t

  */
  int i,j,k,var;
  float xr,xm,dx,s,*dtau_i,*tau_i,**CN,**GN;
  float a=0.0,b=A;

  tau_i=vector(1,10);

  /* The tau_i */
  
  xm=0.5*(b+a);
  xr=0.5*(b-a);
  s=0;
  for(j=1,k=1;j<=5;j++,k += 2) {
    dx=xr*x[j];
    tau_i[k]=xm-dx;
    tau_i[k+1]=xm+dx;
  }

  
  /* 
     We need 10 run on Cn for each tau_i and the 
     corresponding dtau_i such that tau_i is included in the results and 
     is the last one
  */

  dtau_i=vector(1,10);
  
  for (i=1;i<=10;i++){
    dtau_i[i] = tau_i[i]/MAXTAU;
  
  }

  /* CN for the variable var
   var=1 for C
   var=2 for O
   var=3 for N
   var=4 for X
  */

  var=1;

  CN=matrix(1,10,1,NSTEP);
  
  for (i=1;i<=10;i++)
/*     biode(tau_i[i],dtau_i[i],var,CN[i]); */
    

  /* Obtain by inversion the corresponding gn(tau_i;x1) */

  GN=matrix(1,10,1,NSTEP);

  /* To be continue ... */


  free_matrix(GN,1,10,1,NSTEP);
  free_vector(dtau_i,1,10);
  free_matrix(CN,1,10,1,NSTEP);
  free_vector(tau_i,1,10);
  
}

int main() {

  int i,nn;
  float t0=0.0,*t,**f,x,tau;
  
  t=vector(1,N);
  f=matrix(1,M,1,N);

  voltra1(N,M,t0,H,t,f,gfun,ak);
  
/*   for (tau=0;tau<=2;tau += 0.005){ */
/*     printf("%f %f %f %f %f\n",tau,gfun(tau,.5),gfun(tau,1),gfun(tau,2),gfun(tau,3)); */
/*   } */

  for (nn=1;nn<=N;nn++)
    printf("%f %f %f\n",t[nn],f[1][nn]/2,gfun(TAU,t[nn])); /* Warning f/2 */

  

  free_vector(t,1,N);
  free_matrix(f,1,M,1,N);
  
  
  return 0;

}
  
  
