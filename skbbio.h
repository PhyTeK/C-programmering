
IN_PAR(L,double,1.0,%g,Model length)
CHECK_PAR(L>0)

IN_PAR(v,double,1.0,%g,Advective velocity)
CHECK_PAR(v>=0)

IN_PAR(interv,int,512,%d,Number of intervals)
CHECK_PAR(interv>0)

IN_PAR(D0,double,0.0,%g,Diffusion coefficient for S (H2))
CHECK_PAR(D0>=0)

IN_PAR(D1,double,0.0,%g,Motility for u)
CHECK_PAR(D1>=0)

IN_PAR(D2,double,0.0,%g,Motility for w)
CHECK_PAR(D2>=0)

IN_PAR(D3,double,0.0,%g,Motility for u2)
CHECK_PAR(D3>=0)

IN_PAR(D4,double,0.0,%g,Motility for w2)
CHECK_PAR(D4>=0)

IN_PAR(aw,double,0.1,%g,Monod function coef for wall bacteria)
CHECK_PAR(aw>=0)

IN_PAR(a,double,0.1,%g,Monod function coef for free bacteria)
CHECK_PAR(a>=0)

IN_PAR(mw,double,1.0,%g,Monod function coef for wall bacteria)
CHECK_PAR(mw>=0)

IN_PAR(aaw,double,0.1,%g,Monod function coef for wall bacteria 2)
CHECK_PAR(aaw>=0)

IN_PAR(aa,double,0.1,%g,Monod function coef for free bacteria 2)
CHECK_PAR(aa>=0)

IN_PAR(aG,double,0.1,%g,Wall-bound function coef for free bacteria 1)
CHECK_PAR(aG>=0)

IN_PAR(aaG,double,0.1,%g,Wall-bound function coef for free bacteria 2)
CHECK_PAR(aaG>=0)

IN_PAR(mmw,double,1.0,%g,Monod function coef for wall bacteria 2)
CHECK_PAR(mmw>=0)

IN_PAR(m,double,1.0,%g,Monod function coef for free bacteria)
CHECK_PAR(m>=0)

IN_PAR(mm,double,1.0,%g,Monod function coef for free bacteria 2)
CHECK_PAR(mm>=0)

IN_PAR(S0,double,1.0,%g,Concentration of Hydrogen at x=0)
CHECK_PAR(S0>=0)

IN_PAR(u0,double,0.0,%g,Concentration of free bacteria at x=0)
CHECK_PAR(u0>=0)

IN_PAR(w0,double,1.0,%g,Concentration of wall bacteria at x=0)
CHECK_PAR(w0>=0)

IN_PAR(uu0,double,0.0,%g,Concentration of free bacteria 2 at x=0)
CHECK_PAR(uu0>=0)

IN_PAR(ww0,double,0.0,%g,Concentration of wall bacteria 2 at x=0)
CHECK_PAR(ww0>=0)

IN_PAR(Si,double,0.0,%g,Initial concentration of Hydrogen)
CHECK_PAR(Si>=0)

IN_PAR(ui,double,0.0,%g,Initial concentration of free bacteria)
CHECK_PAR(ui>=0)

IN_PAR(wi,double,0.0,%g,Initial concentration of wall bacteria)
CHECK_PAR(wi>=0)

IN_PAR(uui,double,0.0,%g,Initial concentration of free bacteria 2)
CHECK_PAR(uui>=0)

IN_PAR(wwi,double,0.0,%g,Initial concentration of wall bacteria 2)
CHECK_PAR(wwi>=0)

IN_PAR(winf,double,1.0,%g,Maximum occupation upper bound)
CHECK_PAR(winf>=0)

IN_PAR(kf,double,0.2,%g,Dead ratio for free bacteria)
CHECK_PAR(kf>=0.0)

IN_PAR(kw,double,0.2,%g,Dead ratio for wall bacteria)
CHECK_PAR(kw>=0.0)

IN_PAR(kkf,double,0.2,%g,Dead ratio for free bacteria 2)
CHECK_PAR(kkf>=0.0)

IN_PAR(kkw,double,0.2,%g,Dead ratio for wall bacteria 2)
CHECK_PAR(kkw>=0.0)

IN_PAR(alpha,double,0.2,%g,Attachement constant)
CHECK_PAR(alpha>=0)

IN_PAR(beta,double,0.2,%g,Detachement constant)
CHECK_PAR(beta>=0)

IN_PAR(aalpha,double,0.2,%g,Attachement constant for bacteria 2)
CHECK_PAR(aalpha>=0)

IN_PAR(bbeta,double,0.2,%g,Detachement constant for bacteria 2)
CHECK_PAR(bbeta>=0)

IN_PAR(epsilon,double,10.0,%g,Non dimensional parameter)
CHECK_PAR(epsilon>=0)

IN_PAR(eepsilon,double,10.0,%g,Non dimensional parameter 2)
CHECK_PAR(eepsilon>=0)

IN_PAR(BCL,int,1,%d,Boundary Condition at x=L: \
 BCL=1  Diricklet BCL=2 Neumann BCL=3 +Fluxes)
CHECK_PAR(BCL>0)
CHECK_PAR(BCL<4)

IN_PAR(BC0,int,1,%d,Boundary Condition at x=0: \
 BC0=1  Diricklet BC0=2 Neumann BC0=3 +Fluxes)
CHECK_PAR(BC0>0)
CHECK_PAR(BC0<4)
  
IN_PAR(Jin,double,0.1,%g,In flow of H2 at x=0)
CHECK_PAR(Jin >= 0)

IN_PAR(Jout,double,0.1,%g,Out flow of H2 at x=L)
CHECK_PAR(Jout >= 0)













