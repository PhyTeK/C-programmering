/* 011016: Maple BJ jacobians 
   011022: News BJ normalised
   011023: New reaction part
   011119: 2 bacteria
 */

/* File doc/PDMtemplate; TRANSIENT V2.51; 14 Aug. 2001 */
/* (C) 2001 P. Martinet */

/* Template for the Problem Definition Module source file(s)	*/
#define trueSteady 0
#define PDM_H "demo/skb/skbbio.h"

/* PDMpath is the path of the PDMheader.h relative to the
		TRD (TRANSIENT Root Directory);
	   PDMheader.h is the name of the PDM header file that goes with
		this source file	*/
#include "../../transien.h"
/*  TRDpath is the path of the TRD relative to the directory
		where this source file will be located	*/

/* The above two lines are mandatory, as well as all the function
   declarations bellow with the exception of reaction_part().
   However, what is in the body of individual functions is more or less
   just a recommendation (except for the return statement of get_neqns_etc(),
   and the 'nnodes =  ...' line in generate_grid().
   And something to hold the addresses of c, c_o, F, A, ...., Y must
   always be declared.
*/
#ifndef Y2S
#define Y2S 31556926. /* Mean Astronomical Year */
#endif
#define NEQNS	5  /* Number of PDE */


typedef double ROW[NEQNS];


/* Declaration of other global variables not listed in PDMheader.h
	(not being input parameters) will come here */

static double New,eNew,aNew,InveNew,InvtNew,Ma,Maw;
static double MaNew,MawNew,MaweNew,alphaeNew,aalphaeNew,Idh,Idhh,NewIdh;
static double alphae,kwbeta,aalphae,kkwbeta,TauIdh,Newh2,dIdhh;
static double LvmS0,LvmwS0,LsV;

static ROW *c, *c_o; /* new and old "concentration" as treated
			as 2-dimensional matrices here */

static ROW *AJ, *BJ, *DJ, *XJ, *YJ; /* Jacobian submatrices  - 2-DIM */

static double *FD; /* Discretized F_k */

static double *x,h; /* grid */


inline void BJacob (double *);


int get_neqns_etc(void)
/*===================*/
{
    title = "   Deep microbial live";
    time_units = "";
    tEnd = Y2S*1000; /* 1000 years */
    tOutInit = Y2S*10; /* 10 years */
    tStepInit = 0.5*Y2S;
    return NEQNS;
}


void generate_grid(void)
/*====================*/
{
    nnodes = interv+1;   /* total number of nodes */
    h = L/(double)interv; 
}


void Initialize(double *c_vec, double *c_o_vec)
/*===========================================*/
{
  double fw,ffw,X, XX;
  int j;

  c = (ROW *)c_vec;
  c_o = (ROW *)c_o_vec;
  
  fw = (L/v)*mw*S0/(aw + S0);
  ffw = (L/v)*mmw*S0/(aaw + S0);
  X = (kw + beta)/fw;
  XX = (kkw + bbeta)/ffw;
  w0 = ((1 + aG)*X - 1.)/(X - 1.);
  ww0 = ((1 + aaG)*XX - 1.)/(XX - 1.);

    for(j=0; j<nnodes; j++) {
	c_o[j][0] = c[j][0] = Si; 
	c_o[j][1] = c[j][1] = ui;
  	c_o[j][2] = c[j][2] = wi; 
	c_o[j][3] = c[j][3] = uui;
  	c_o[j][4] = c[j][4] = wwi;
    }

}


void initNL3BANDV(double***LBp, double***UBp, double*F,
		  double*A, double*B, double*D)
/*===================================================*/
{
/* Declare arrays of pointers to BDA arrays: */
    static double *LB[NEQNS],*UP[NEQNS];
/* Or if there are e.g. no upper bounds for all unknowns,
   there is no need to declare UB
*/
    int j;
/* Declare some BDA, for EXAMPLE: */
    static double Lower[2],Upper[2];

    *LBp = LB;
    *UBp = UP;

    /* or if there are no upper bounds (no UB): *UBp = NULL */

    Lower[0] = 1.;  /* Bound type: constant */
    Lower[1] = 0.;  /* Bound value */
    Upper[0] = 1.;
    Upper[1] = 1.;

    for(j=0; j<NEQNS; j++) {
      LB[j] = Lower; 
      UP[j] = Upper;
    }

/* Save addresses of F, Jacobian submatrices: */

    FD= F;

    AJ = (ROW *)A;
    BJ = (ROW *)B;
    DJ = (ROW *)D;
    XJ= AJ;
    YJ= DJ;

}

#define MAXLETT 80
/* xxx should be larger than the length of the longest line + 1 */

#define GETLINE  p= fgets(text,MAXLETT,infile);\
if(p == NULL) break;\
if(text[strlen(text)-1] != 10) {\
PE"\nInput file has lines longer then %d characters!\n", MAXLETT-1);\
PE"Please, increase MAXLETT in %s!\n\n", __FILE__); exitTR(13);}


FILE *infile;

void read_profiles(char *filename)
/*==============================*/
/* It is assumed that infile was created by pr_profiles */
{
char *p, text[MAXLETT], comnt;
int j=0, k;

if((infile = fopen(filename,"r")) == NULL) {
PE"Failed to open  %s\n", filename);
PE"Initializations made in Initialize() unchanged\n");
return;
}

for(;;) {
GETLINE
comnt= *p; /* CmntC may differ between runs */
break;
}
do {
GETLINE
if(*p != '\n' && *p != comnt) break;
/* Process file header */
} while(1);

do {
sscanf(text, " %*e %le %le %le", &c_o[j][0], &c_o[j][1], &c_o[j][2]);
j++;
GETLINE

} while(1);
(void)fclose(infile);

if(!(steady & 1)) 
for(j= 0; j<nnodes; j++)
for(k=0; k<NEQNS; k++) c[j][k]= c_o[j][k];
}

void auxiliary_parameters_for_evalf(void)
/*=====================================*/
{
Idh = 1./(2*h);
Idhh = 1./(h*h);
TauIdh = Tau*Idh;
New = 1. - Tau;
dIdhh = 2.0*Idhh*New;
NewIdh = New*Idh;
Newh2 = New/(h*h);
eNew = epsilon*New;
aNew = a*New;
InveNew = New/epsilon;
InvtNew = New*tStepInv;
LsV = L/v;
Ma = m*a;
Maw = mw*aw;
MaNew = Ma*New;
MawNew = Maw*New;
MaweNew = epsilon*MawNew;
alphae = alpha/epsilon;
aalphae = aalpha/eepsilon;
aalphaeNew = aalphae*New;
alphaeNew = alphae*New;
kwbeta = kw + beta;
kkwbeta = kkw + bbeta;
LvmS0 = LsV*m*S0;
LvmwS0 = LsV*mw*S0;

}

inline void BJacob(double *C){
double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17;
double t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32;
double t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45, t46, t47;
double t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60, t61, t62;
double t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75, t76, t77;
double t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90, t91, t92;
double t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105;
double t106, t107, t108, t109, t110, t111, t112, t113, t114, t115, t116, t117;
double t118, t119, t120, t121, t122, t123, t124, t125, t126, t127, t128, t129;
double t130, t131, t132, t133, t134, t135, t136, t137, t138, t139, t140, t141;
double t142, t143, t144, t145, t146, t147, t148, t149, t150, t151, t152, t153;
double t154, t155, t156, t157, t158, t159, t160, t161, t162, t163, t164, t165;
double t166, t167, t168, t169, t170, t171, t172, t173, t174, t175, t176, t177;
double t178, t179, t180, t181, t182, t183, t184, t185, t186, t187, t188, t189;
double t190, t191, t192, t193, t194, t195, t196, t197, t198, t199;
double t200, t201, t202, t203, t204, t205, t206, t207, t208, t209, t210;
double t211, t212, t213, t214, t215, t216, t217, t218, t219, t220;
double t221, t222, t223, t224, t225, t226, t227, t228, t229, t230;
double t231, t232, t233, t234, t235, t236, t237, t238, t239, t240;
double t241, t242, t243, t244, t245, t246, t247, t248, t249, t250;


      t1 = aw+C[0];
      t2 = 1/t1;
      t3 = t1*t1;
      t4 = 1/t3;
      t5 = C[0]*t4;
      t10 = m*C[1];
      t11 = a+C[0];
      t12 = 1/t11;
      t14 = t11*t11;
      t15 = 1/t14;
      t18 = C[3]*mm;
      t19 = aa+C[0];
      t20 = 1/t19;
      t22 = t19*t19;
      t23 = 1/t22;
      t26 = eepsilon*C[4];
      t27 = aaw+C[0];
      t28 = 1/t27;
      t29 = mmw*t28;
      t31 = mmw*C[0];
      t32 = t27*t27;
      t33 = 1/t32;
      t37 = 1/v;
      t39 = L*New;
      t41 = t39*t37;
      t42 = m*C[0];
      t43 = t42*t12;
      t47 = t37*mw;
      t48 = C[0]*t2;
      t51 = mm*C[0];
      t52 = t51*t20;
      t57 = C[0]*t28;
      t62 = 1/winf;
      t63 = (C[2]+C[4])*t62;
      t64 = 1.0-t63;
      t65 = aG+1.0-t63;
      t66 = 1/t65;
      t67 = t64*t66;
      t68 = 1.0-t67;
      t70 = C[2]*C[0];
      t83 = L*t37;
      t85 = alpha*t62;
      t87 = C[4]*t62;
      t88 = 1.0-t87;
      t93 = t62*t66;
      t94 = t65*t65;
      t97 = t64/t94*t62;
      t99 = t2*(t93-t97);
      t103 = t47*L;
      t105 = t85*C[1];
      t124 = alpha/epsilon;
      t125 = t62*C[2];
      t133 = C[2]*(-t48*t93+t48*t97);
      t139 = t124*t62*C[1];
      t149 = aaG+1.0-t63;
      t150 = 1/t149;
      t151 = t64*t150;
      t152 = 1.0-t151;
      t155 = t26*mmw;
      t164 = t149*t149;
      t167 = t64/t164*t62;
      t168 = t62*t150-t167;
      t173 = aalpha*C[3]*t62;
      t206 = C[4]*(-t31*t28*t62*t150+t31*t28*t167);
      t209 = aalpha/eepsilon;
      t211 = t209*C[3]*t62;
      
      BJ[0][0] = ((-t2+t5)*epsilon*C[2]*mw-t10*t12+t10*C[0]*t15-t18*t20+\
      t18*C[0]*t23-t26*t29+t26*t31*t33)*t37*t39;
      BJ[0][1] = -t41*t43;
      BJ[0][2] = -epsilon*New*L*t47*t48;
      BJ[0][3] = -t41*t52;
      BJ[0][4] = -eepsilon*New*L*t37*mmw*t57;
      BJ[1][0] = ((C[2]*t2*t68-t70*t4*t68)*epsilon*mw+C[1]*(m*t12-t42*t15))*t37*t39;
      BJ[1][1] = (t83*t43+t85*C[2]-kf-alpha*t88)*New;
      BJ[1][2] = ((t48*t68+t70*t99)*epsilon*t103+t105+epsilon*beta)*New;
      BJ[1][4] = (epsilon*C[2]*t83*mw*C[0]*t99+t105)*New;
      BJ[2][0] = C[2]*(t2*t64*t66-t5*t67)*mw*t41;
      BJ[2][1] = (-t124*t125+t124*t88)*New;
      BJ[2][2] = ((t133+t48*t67)*mw*t83-t139-kw-beta)*New;
      BJ[2][4] = (t133*t103-t139)*New;
      BJ[3][0] = (C[3]*(mm*t20-t51*t23)+t26*t29*t152-t155*C[0]*t33*t152)*t37*t39;
      BJ[3][2] = (t26*t83*t31*t28*t168+t173)*New;
      BJ[3][3] = (t83*t52+aalpha*C[2]*t62-kkf-aalpha*t88)*New;
      BJ[3][4] = ((eepsilon*mmw*t57*t152+t155*t57*t168)*t37*L+t173+eepsilon*bbeta)*New;
      BJ[4][0] = C[4]*(t29*t151-t31*t33*t64*t150)*t41;
      BJ[4][2] = (t206*t83-t211)*New;
      BJ[4][3] = (-t209*t125+t209*(1.0-t87))*New;
      BJ[4][4] = ((t206+t31*t28*t64*t150)*t37*L-t211-kkw-bbeta)*New;
}


void reaction_part(double *C)
/*=========================*/
/* Here C is the mixed old-nex solution for a single node */
{
  double f0,ff0,fw,ffw,G2,GG2;
  double W;

  W = (C[2]+C[4])/winf;

    f0 = LsV*m*C[0]/(a + C[0]);
    ff0 = LsV*mm*C[0]/(aa + C[0]);
    fw = LsV*mw*C[0]/(aw + C[0]);
    ffw = LsV*mmw*C[0]/(aaw + C[0]);
    G2 = (1. - W)/(aG + 1.- W);
    GG2 = (1. - W)/(aaG + 1.- W);

    FD[0] = -C[1]*f0-epsilon*C[2]*fw -C[3]*ff0-eepsilon*C[4]*ffw;
    FD[1] = C[1]*(f0-kf) + epsilon*C[2]*fw*(1.-G2);
    FD[1] += -alpha*C[1]*(1.- W) + epsilon*beta*C[2];
    FD[2] = C[2]*(fw*G2 - kwbeta) + alphae*C[1]*(1.0 - W);
    FD[3] = C[3]*(ff0-kkf) + eepsilon*C[4]*ffw*(1.-GG2);
    FD[3] += -aalpha*C[3]*(1.- W) + eepsilon*bbeta*C[4];
    FD[4] = C[4]*(ffw*GG2 - kkwbeta) + aalphae*C[3]*(1.0 - W);

    BJacob(C);

}

void evalf(int j)
/*===============
  Remember, that in the calling function, XJ=AJ and YJ=DJ. All elements
  of AJ, BJ, DJ (and thus XJ and YJ) are set to zero in NL3BANDV before
  evalf() is called. Thus only their nonzero elements need to be set here.
*/
{
    static ROW VM, V0, VP;
    static double *cBarM, *cBar0, *cBarP;
    int k;

    if(j == 0){
	/*
	  NL3BANDV first calls evalf() with j=0 and
	  on each subsequent call j is increased by 1.
	*/
	
	for(k=0; k<NEQNS; k++) {
	    V0[k] = Tau * c_o[0][k] + New * c[0][k];
	    VP[k] = Tau * c_o[1][k] + New * c[1][k];
	}
	cBarM= (double *)VM;
	cBar0= (double *)V0;
	cBarP= (double *)VP;

	/* Now do the x=0 boundary conditions: 
	   All concentrations constants */
	if (BC0 == 1){

	    FD[0] = c[0][0] - S0;
	    FD[1] = c[0][1] - u0;
	    FD[2] = c[0][2] - w0;
	    FD[3] = c[0][3] - uu0;
	    FD[4] = c[0][4] - ww0;

	    BJ[0][0] = 1.0;
	    BJ[1][1] = 1.0;
	    BJ[2][2] = 1.0;
	    BJ[3][3] = 1.0;
	    BJ[4][4] = 1.0;

	} else if (BC0 == 2 || BC0 == 3){
	  /* First derivatives = 0 at x=0 */
	  for(k=0;k<NEQNS;k++){
	    FD[k] = -3.* c[0][k] + 6.*c[1][k] -1.*c[2][k];
	    BJ[k][k] = -3.;
	    DJ[k][k] = 6.;
	    XJ[k][k] = -1.;
	  }

	} if (BC0 == 3) {
	    /* Constant in flow at x=0 for the diffusive species (H2) */
	    FD[0] = D0*Idh*(-3.*c[0][0] + 6.*c[1][0] - c[2][0]) + Jin;
	    BJ[0][0] = -3.*D0*Idh;
	    DJ[0][0] =  6.*D0*Idh;
	    XJ[0][0] =  1.*D0*Idh;
	}
    } else if (j<interv) {

	{ double *p;
	register k;
	p = cBarM;
	cBarM = cBar0;
	cBar0 = cBarP;
	cBarP = p;
	for(k=0; k<NEQNS; k++)
	    cBarP[k] = Tau * c_o[j+1][k] + New * c[j+1][k];
	}

	reaction_part(cBar0);

	FD[0] += Idh*(cBarM[0]-cBarP[0]);
	FD[1] += Idh*(cBarM[1]-cBarP[1]);
	FD[3] += Idh*(cBarM[3]-cBarP[3]);
	
	FD[0] += D0*Idhh*(cBarP[0]-2.*cBar0[0]+cBarM[0]);
	FD[1] += D1*Idhh*(cBarP[1]-2.*cBar0[1]+cBarM[1]);
	FD[2] += D2*Idhh*(cBarP[2]-2.*cBar0[2]+cBarM[2]);
	FD[3] += D3*Idhh*(cBarP[3]-2.*cBar0[3]+cBarM[3]);
	FD[4] += D4*Idhh*(cBarP[4]-2.*cBar0[4]+cBarM[4]);

	BJ[0][0] -= dIdhh*D0;
	BJ[1][1] -= dIdhh*D1;
	BJ[2][2] -= dIdhh*D2;
	BJ[3][3] -= dIdhh*D3;
	BJ[4][4] -= dIdhh*D4;

	/* time derivative: */
	for (k=0; k<NEQNS; k++){
	    FD[k] += (c_o[j][k] - c[j][k]) * tStepInv;
	    BJ[k][k] -= tStepInv;
	}    
	
	AJ[0][0] = NewIdh + D0*Newh2;
	AJ[1][1] = NewIdh + D1*Newh2;
	AJ[2][2] = D2*Newh2;
	AJ[3][3] = NewIdh + D3*Newh2;
	AJ[4][4] = D4*Newh2;
	DJ[0][0] = -NewIdh + D0*Newh2;
	DJ[1][1] = -NewIdh + D1*Newh2;
	DJ[2][2] = D2*Newh2;
	DJ[3][3] = -NewIdh + D3*Newh2;
	DJ[4][4] = D4*Newh2;
 
    } else if (j == interv) { /* x=L boundary conditions:*/

	if(BCL == 1){
	    FD[0] = c[j+1][0] - Si;
	    FD[1] = c[j+1][1] - ui;
	    FD[2] = c[j+1][2] - wi;
	    FD[3] = c[j+1][3] - uui;
	    FD[4] = c[j+1][4] - wwi;

	    BJ[0][0] = 1.0;
	    BJ[1][1] = 1.0;
	    BJ[2][2] = 1.0;
	    BJ[3][3] = 1.0;
	    BJ[4][4] = 1.0;

	}else if (BCL == 2 || BCL == 3){
	    /* All differential are zero at x=L (j=N-1) */
	    for(k=0;k<NEQNS;k++){
		FD[k] = 3.* c[j][k] -4.* c[j-1][k] + c[j-2][k];
		BJ[k][k] = 3.;
		AJ[k][k] = -4.;
		YJ[k][k] = 1.;
	    }
	    
	} if (BCL == 3) {
	    /* Constant out flux of H2 at x=L */
	    FD[0] = D0*Idh*(3.*c[j][0] -4.*c[j-1][0] + c[j-2][0]) + Jout;
	    AJ[0][0] = -4.*D0*Idh;
	    BJ[0][0] = 3.*D0*Idh;
	    YJ[0][0] = 1.*D0*Idh;
	}
    }
}


void pr_profiles(int MaxIter)
/*=========================*/
{
    int j, i;

    if(trueSteady) {
	PR"%c\t\tSteady State results\n%c\n", CmntC, CmntC);
    PR"%c No. of iterations done = %d\n", CmntC, MaxIter);
} else {
PR"%c Transient results for t = %.10e sec = %e years\n", CmntC, t, t/Y2S);
PR"%c Last tStep = %e sec = %.10e years\n", CmntC, tStepLast, tStepLast/Y2S);
PR"%c Max. # of iterations done since last profile = %d\n",CmntC, MaxIter);
PR"%c Stability Condition: fw(S0)G(0)-kw-beta ( %2.4f ) > 0\n", CmntC, \
mw*S0/(aw+S0)/(a+1.)-kw-beta);
}
report_max_resid(out, CmntC);

PR"%c\t\t All concentrations in mol/dm^3\n",CmntC);
PR"%c   x(m)      t(s)        [S]          [u]          [w]         [u2]        [w2]\n", CmntC);
for(j= 0;j<nnodes;j++) {
PR" %.5f ", j*h);
PR"%e %e %e %e %e %e",t,c_o[j][0],c_o[j][1],c_o[j][2],c_o[j][3],c_o[j][4]);
PR"\n");
}
(void)fflush(out);
}


void init_Other(void)
/*=================*/
{ }

int update_Other(void)
/*==================*/
{ }

void pr_Other(int only_Other)
/*===========================*/
{ }










