/* Exact solution for the kinetic model (Thomas 48 & Amundsom 50)
   p. 308 Cvetkovic & Dagan 1996

   compile with : gcc -O2 -o kinex kinex.c -lm -lrecipes_c

   P. Martinet 12 Mars 2000
*/

#include <stdio.h>
#include <math.h>
#include "nr.h"

float fphi(float ,float );
float ffphi(float ,float );
float Gamma(float,float,float,float,float,float);
float kern(float,float);
float qgausp(float (*func)(float,float), float , float ,float );
float Nsol(float,float,float);

float fphi(float u,float v){
  int m,n;
  float phi=0.0;

  for (m=0;m<30;m++)
    for (n=0;n<m;n++)
      phi += pow(u,m)*pow(v,n)/factrl(m)/factrl(n);
  
return phi;
}


float kern (float x,float v){
  return bessi0(2*sqrt(v*x))*exp(-x); 
}

float ffphi(float u, float v){
  return exp(u)*qgausp(kern,0,u,v);
}

float Gamma(float t,float tau, float k1,float k2,float Nm,float C0){
  float r,s,xi,x,bsx,fpr,num;
  if (t<=tau) return 0.0;
  r=C0*k1+k2;
  s=k1*k2*Nm/r;
  xi=t-tau;
  x=2*sqrt(k1*k2*Nm*tau*xi);
  bsx=bessi0(x);
  fpr=fphi(r*xi,s*tau);
  num=bsx+fpr;
  return num/(num+fphi(k1*Nm*tau,k2*xi));
}

float Nsol(float C,float B,float Nm){
  return Nm*B*C/(1+B*C);
}

int main(){
  float t,tau;
  float C0=1.0,Nm=0.2,k1=100.0,k2=1.0;
  float r,s,C,B;
  r=C0*k1+k2;
  s=k1*k2*Nm/r;
  B=k1/k2;
  printf("# Cin=%f C0=%f Nm=%f k1=%f k2=%f\n",0,C0,Nm,k1,k2);
  printf("# %5s %15s %15s %15s\n","t","tau","C(t,tau)","N(t,tau)");
  for (t = 0;t <= 2;t += .01)
    for (tau = 0;tau <= 1.1;tau += .1){
      if(tau>=1.0) {
	C=C0*Gamma(t,tau,k1,k2,Nm,C0);
	printf("%e %e %e %e\n",t,tau,C,Nsol(C,B,Nm));
      }
    }
return 0;
}



float qgausp(float (*func)(float,float), float a, float b, float v)
{
	int j;
	float xr,xm,dx,s;
	static float x[]={0.0,0.1488743389,0.4333953941,
		0.6794095682,0.8650633666,0.9739065285};
	static float w[]={0.0,0.2955242247,0.2692667193,
		0.2190863625,0.1494513491,0.0666713443};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*((*func)(xm+dx,v)+(*func)(xm-dx,v));
	}
	return s *= xr;
}

