static char rcsid[]= "$Id: bac.c,v 1.10 2000/01/09 22:07:23 root Exp $";

/*    BIOREMEDIATION
      --------------

      We solve the general transport equation:


      RX dX/dt + dX/dtau = Psi_X

      where, Psi_X are the sink/source term for the quantity X

      Psi_C = Psi_0/F0 + Psi_N/FN
      Psi_O = -mu0/Y0 C/(KCO + C) O/(KO + O) X FO S(O)
      Psi_N = -muN/YN C/(KCN + C) N/(KN + N) X FN (1-S(O))
      Psi_X = -YO/FO Psi_O - YN/FN Psi_N - lambda X
      
      S(0) = 0 if O < OT
      S(O) = 1 if O > OT    S is called the switching function
      
      A1 = mu_O.X/Y0 and A2 = mu_O.X.F/YO, where mu_O is the maximum
      microbial growth rate for aerobic metabolism, X is the
      concentration of bacterial, F is the stoichiometric ratio between
      oxygen and contaminant. Finally, Y= is the yield coefficient for
      aerobic metabolim. 
      
      Kco = the half-saturation constant of the contaminant for aerobic
      metabolism.  Ko = the half-saturation constant of oxygen for
      aerobic metabolism.
      
      with `Def' option we defined the following parameters:

      Microcolony Kinetic Parameters (Molz et al. 1986)
      ------------------------------
      Kco = 0.120 mg/cm3
      Ko = 0.00077 mg/cm3
      mu_0 = 4.34 1/days
      Y0 = 0.278
      X = 90 mg/cm3
      
      Then we have A1=1405 mg/cm3/day and A2=A1.F(O)
      
      with `KCB' option the parameters are defined in the paper from
      Kaluarachchi, Cvetkovic and Berglund (1999).
      
      
      modified: 990710
      990917: loop in tau values

      Revisions
      ---------
      
      $Log: bac.c,v $
      Revision 1.10  2000/01/09 22:07:23  root
      Working version: Good results for dtau=0.01
      For dtau=1 somes differences with KCB results

      Revision 1.9  2000/01/07 17:45:39  root
      2D alpha version
      good results for dtau=1
      quite strange for dtau=0.01

      Revision 1.8  2000/01/05 03:42:38  root
      2D results (t,tau)

      Revision 1.7  2000/01/03 15:16:56  root
      Bugs in derivs

      Revision 1.6  2000/01/03 04:19:02  root
      KCB with tau-coordinate

      Revision 1.5  2000/01/02 18:27:27  root
      KCB model

      Revision 1.4  1999/12/27 16:34:53  root
      Working version

      Revision 1.3  1999/12/23 04:28:26  root
      Used the recipes_c library

      Revision 1.2  1999/12/22 16:21:15  root
      Description of the problem

      Revision 1.1  1999/12/22 13:17:06  root
      Initial revision
      
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "nr.h"
#include "lasar.h"

typedef struct PARM {
  float Ao,An,mu0,muN,Y0,YN,Fo,Fn;
  float Kco,Ko,Kcn,Kn,C,O,N,X,Xcrit;
  float lambda;
} parm;

parm P,*pm=&P;
float dxsav,*xp,**yp,**yold,h;
int kmax,kount,nvar,idx;

/* Switching function S(O) */

int Swt(float x){
  if (x <= pm->Xcrit) return 0;
  if (x >  pm->Xcrit) return 1;
}

/* Derivate dC/dt, dO/dt, dN/dt, dX/dt */

void derivs(float x,float *y,float *dydx){
  float C,O,N,X,Fo,Fn;

  float Ao,An,Kco,Ko,Kcn,Kn,CO,CN;
  float R=1,Rc=R,Ro=1,Rn=1,Rx=R;   /* 1<R<3 */
  float Gc,Go,Gn,Gx;
  float Cold,Oold,Nold,Xold,er;
  int i;

  Ao = -pm->mu0/pm->Y0; 
  An = -pm->muN/pm->YN;

  Fo = pm->Fo;
  Fn = pm->Fn;
  Kco = pm->Kco;
  Ko  = pm->Ko;
  Kcn = pm->Kcn;
  Kn  = pm->Kn;

  C = (*(y+1));
  O = (*(y+2));
  N = (*(y+3));
  X = (*(y+4));

  /* Mono-kinetic source/sink terms */

  CO = Ao*C*O/(Kco+C)/(Ko+O)*X*Swt(O);
  CN = An*C*N/(Kcn+N)/(Kn+N)*X*(1-Swt(O));

  Cold = yold[2][idx];
  Oold = yold[3][idx];
  Nold = yold[4][idx];
  Xold = yold[5][idx];

  /* Gradient dC/dtau */

  Gc = (C-Cold)/dtau;
  Go = (O-Oold)/dtau;
  Gn = (N-Nold)/dtau;
  Gx = (X-Xold)/dtau;

  *(dydx+1) = (CO + CN - Gc)/Rc;
  *(dydx+2) = (CO*Fo - Go)/Ro;
  *(dydx+3) = (CN*Fn - Gn)/Rn;
  *(dydx+4) = (- pm->Y0*CO - pm->YN*CN - pm->lambda*X - Gx)/Rx;
}


int biode(float taui,float dtau,int var,float *cn)
{
  /* Solve the ODE for bioremediation of 4 independent variables
     C,O,N,X (var=1,2,3,4)
  */
  float eps=1e-6,x1=T1,x2=T2,*y;
  float h1=0.01,hmin=0.0;
  float h,*dydxp,*yout,x,tau;
  int nok,nbad;
  int i=0,j=0,m,k;
  char *s;

  div_t result1,result2;

  /* Init the function parameters */
  /* Kaluarachchi, Cvetkovic, Berglund example */
  printf("# KCB case\n");

  pm->Kco= 0.008;
  pm->Kcn= 0.008;
  pm->Kn = 0.0002;
  pm->Ko = 0.0002; 
  pm->C = 0.010;   /* for 10 days */
  pm->O = 0.010;
  pm->N = 0.015;
  pm->X = 0.0005;
  pm->Xcrit = 0.0002;
  pm->mu0= 1.98;
  pm->muN= 1.98;
  pm->Y0 = 0.1;
  pm->YN = 0.1;
  pm->Fo = 3.13;
  pm->Fn = 4.72;
  pm->lambda=0.01;

  xp = vector(1,MAXSTP);
  yp = matrix(1,NVAR,1,MAXSTP);
  dydxp = vector(1,NVAR);
  y = vector(1,NVAR);
  yold = matrix(1,NVAR+1,1,MAXSTP);  /* xCONX at i-1 */
  yout = vector(1,NVAR);

  kmax=MAXSTP;

  dxsav=(x2-x1)/kmax;
  h=(x2-x1)/MAXSTP;

  /* Initial values for yold and yp at tau=0*/

  for(i=1;i<=MAXSTP;i++){
    yold[1][i] = h*(i-1);
    yold[2][i] = yp[1][i] = pm->C;
    if (i >= 10.0/h)
      yold[2][i] = yp[1][i] = 0.0;
    yold[3][i] = yp[2][i] = pm->O;
    yold[4][i] = yp[3][i] = pm->N;
    yold[5][i] = yp[4][i] = pm->X;
  }


  /* Apply Runge Kutta's method */

  printf("# Runge Kutta\n");

  j=1;

  for(tau=0;tau <= TAUMAX + dtau;tau += dtau) {
    j++;

    /* Injected concentrations */

    y[1] = pm->C;
    if(tau >= 5.0) y[1] = 0.0; /* After 10 days tau=5*/

    y[2] = pm->O;
    y[3] = pm->N;
    y[4] = pm->X;

    xp[1]=x1;
    x=x1;

    for (i=1;i<MAXSTP;i++) {
      idx=i;
     
      derivs(x,y,dydxp);

      rk4(y,dydxp,NVAR,x,h,y,derivs);

      if (x+h == x)
	nrerror("Step size too small in routine bac");
      x += h;
      xp[i+1] = x;

      /* insert the result for the variable var in the vector cn for the time i*/

      if(STP(taui)) {
	cn[i] = y[var]; 
      }
      

      yp[1][i] = y[1];
      yp[2][i] = y[2];
      yp[3][i] = y[3];
      yp[4][i] = y[4];
      
    }


    for(i=1;i<MAXSTP;i++){
      yold[1][i] = xp[i];
      yold[2][i] = yp[1][i];
      yold[3][i] = yp[2][i];
      yold[4][i] = yp[3][i];
      yold[5][i] = yp[4][i];
    }
    
  }


  free_vector(xp,1,MAXSTP);
  free_matrix(yp,1,NVAR,1,MAXSTP);
  free_vector(y,1,NVAR);
  free_matrix(yold,1,NVAR+1,1,MAXSTP);
  free_vector(dydxp,1,NVAR);
  free_vector(yout,1,NVAR);

  return 0;

}
















































  
