/* Numerical solution of the kinetic model

   compile with : gcc -O2 -o kinum kinum.c -lm -lrecipes_c

   P. Martinet 12 Mars 2000
*/

#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#define TMAX 2.0
#define TAUMAX 2.0

typedef struct parm{
  float k1,k2,Nm,C0;
  float dt,dtau;
} PARM;

float newC(float,float,float,PARM *);
float newN(float,float,PARM *);

float newC(float C,float N,float oldC,PARM *P){
  float Cn;
  
  Cn = P->k1*C*(N-P->Nm);
  Cn += P->k2*N-(C-oldC)/P->dtau;
  Cn *= P->dt;
  Cn += C;

  return Cn;
}

float newN(float C,float N,PARM *P){
  float Nn;
  
  Nn = -P->k1*C*(N-P->Nm);
  Nn -= P->k2*N;
  Nn *= P->dtau;
  Nn += N;

  return Nn;
}

int main(){
  
  PARM P,*p=&P;

  float t,tau;
  float **C,**N,B;

  int Nt,Ntau,i,j;

  p->k1=50;
  p->k2=10;
  p->Nm=0.2;
  p->dt=0.01;
  p->dtau=0.01;
  p->C0=1.0;

  Nt = (int)ceil(TMAX/p->dt)+1;
  Ntau = (int)ceil(TAUMAX/p->dtau)+1;


  C = matrix(1,Nt,1,Ntau);
  N = matrix(1,Nt,1,Ntau);

  /* Initial & Boundary values */
  
  B=p->k1/p->k2;

  for (i=1;i <= Ntau;i++){
    C[1][i] = p->C0;
    N[1][i] = p->Nm*B*p->C0/(1.0 + B*p->C0);
  }
  for (i=1;i <= Nt;i++)
    C[i][1] = N[i][1] = 0.0;

  for (tau=p->dtau,i=2;tau < TAUMAX + p->dtau;tau += p->dtau,i++)
    for (t=p->dt,j=2;t < TMAX + p->dt;t += p->dt,j++){
      C[j][i] = newC(C[j-1][i],N[j-1][i],C[j-1][i-1],p);
      N[j][i] = newN(C[j-1][i],N[j-1][i],p);
  
    }

  /* Write out the results */
  for (i=Ntau;i<Ntau;i++)
    for (j=1;j<=Nt;j++)
      printf("%e %e %e %e\n",(j-1)*p->dt,(i-1)*p->dtau,C[j][i],N[j][i]);

  free_matrix(C,1,Nt,1,Ntau);
  free_matrix(N,1,Nt,1,Ntau);

  return 0;
}
      







