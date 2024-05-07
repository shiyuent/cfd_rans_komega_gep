#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "tinyexpr.h"


#define NO_SLIP
#define EXPLICIT
#define TINY 1.0e-8

void EquallySpacedGrid(double *y,int N);
void StretchedGrid(double *y,int N);
void ReadGridAndInitialGuessData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N);
void WriteGridAndInitialGuessData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N);
double Calc_Utau(double *U,double *y,double nu);
void LinearInterpolate(double *func,double *y,double yp,double *value,int N);
void ReadRestartData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N);

int main(void)
{
  double betStar,a1,*alp, alp1,alp2,*bet,bet1,bet2,*sigo,sigo1,sigo2,*sigk,sigk1,sigk2,*arg1;
  double *A,*B,*C,*RHS,temp; /* Pointer to arrays used for banded matrices */
  int N,i,iter,maxiter;
  double *y,*k,*om,*U,*d,dx,*T,*S,*F1,*F2,*CDkom,nu,utau, *wallDist;/* utau is the friction velocity*/
  double xi,ubulk;
  double temp1,temp2,temp3;
  double *U_old,*om_old,*k_old;
  double *Temperature,*Temperature_old;
  double L2Norm_U,L2Norm_k,L2Norm_om,iacc,L2Norm_Temperature;
  double y_cut,er_v2,om_cut,er_om;
  double dy,dyplus1;
  double uplus,ThetaPlus;
  double gravity,RayleighNumber,beta,GrasshofNumber,PrandtlNumber,kappa;
  double prodkTilde,prod,Gb,fudge,fudge_om,Prt;

  double *f1Gep, *f2Gep, *I1,*J1, *tempGrad, *velGrad, *uvGrad, *vcGrad, *ajex, *dd;
  double *uc, *vc, *vv, *uv, *uu, *prodk, *L;
  double *alphaGep;

  char InFilename[100],OutFilename[100];
  char DefaultGuess,bc_iacc;

  int stop=0,n;
  FILE *fp;

// reference  
// NASA Langley Research Center. The Menter Shear Stress Transport Turbulence Model, 2015.
// F.R. Menter, M. Kuntz, and R. Langtry. Ten years of industrial experience with the SST turbulence model. In Proceedings of the fourth international symposium on turbulence, heat and mass transfer, pages 625–632, Antalya, Turkey, 2003. Begell House.
// Menter, F. L. O. R. I. A. N. R. (1993, July). Zonal two equation kw turbulence models for aerodynamic flows. In 23rd fluid dynamics, plasmadynamics, and lasers conference (p. 2906).
  /*************************
    Simulation parameters
    ************************/
  //N=769;
   N=201;
 //N=251;
  maxiter=200000;
  dx=1.0; // for the so-called time steps in time discretisation
  
  /*********************
    Model parameters
    ********************/
//  RayleighNumber=1.0e+5;
// RayleighNumber=5.4e+5;
//  RayleighNumber=2.e+6;
//  RayleighNumber=5.e+6;
  RayleighNumber=2.0e+7;
 // RayleighNumber=1.0e+8;
 // RayleighNumber=1.0e+9;

  PrandtlNumber=0.709;
  GrasshofNumber=RayleighNumber/PrandtlNumber;
  
 gravity=9.81;
 kappa=1./RayleighNumber;
 beta=1./(GrasshofNumber*gravity);
 nu=1./GrasshofNumber;

//gravity=1.0;
//beta=1.0;
//nu=sqrt(1/GrasshofNumber);
//kappa=nu/PrandtlNumber;

/*gravity=1.0;
  beta=RayleighNumber;
  nu=1./PrandtlNumber;
  kappa=1.0;*/


  betStar=0.09; //https://en.wikipedia.org/wiki/SST_(Menter%E2%80%99s_Shear_Stress_Transport)#cite_note-1
  a1=0.31; //https://www.openfoam.com/documentation/cpp-guide/html/guide-turbulence-ras-k-omega-sst.html

  alp1=5.0/9.0;
  alp2=0.44;
  bet1=3.0/40.0;
  bet2=0.0828;
  sigo1=0.5;
  sigo2=0.856;
  sigk1=0.85;
  sigk2=1.0;
// simple gradient diffusion hypothesis
//fudge=-1;
fudge=0;
fudge_om=0.0;
Prt=0.9;


  y=dvector(0,N-1);
  wallDist=dvector(0,N-1);
  om=dvector(0,N-1);
  U=dvector(0,N-1);
  Temperature=dvector(0,N-1);
  k=dvector(0,N-1);
  d=dvector(0,N-1);/*d is the eddy viscosity*/

  T=dvector(0,N-1);/*t is the turbulence time scale*/
  velGrad=dvector(0,N-1);/*velGrad is the velocity gradient*/
  S=dvector(0,N-1);/*S is sqrt(2*SijSij) the deformation rate*/
  F1=dvector(0,N-1);/*F1 is the blending function for k-om or k-ep*/
  F2=dvector(0,N-1);/*F2 is for eddy viscosity*/
  alp=dvector(0,N-1);/*the blending alpha*/
  bet=dvector(0,N-1);/*the blending beta*/
  arg1=dvector(0,N-1);/*arg1 in CDkomega and in omega equation*/
  sigk=dvector(0,N-1);/*blending sigk in k equation*/
  sigo=dvector(0,N-1);/*bledning sigo in om equation*/
  CDkom=dvector(0,N-1);/*bledning coefficients*/

  // extra term for GEP 
T=dvector(0,N-1);/*t is the turbulence time scale*/
L=dvector(0,N-1);/*l is the turbulence legth scale*/

f1Gep=dvector(0,N-1);
f2Gep=dvector(0,N-1);

uc=dvector(0,N-1);
vc=dvector(0,N-1);
uu=dvector(0,N-1);
uv=dvector(0,N-1);
vv=dvector(0,N-1);

//  Gb=dvector(0,N-1);
prodk=dvector(0,N-1);

velGrad=dvector(0,N-1);
tempGrad=dvector(0,N-1);
I1=dvector(0,N-1);
J1=dvector(0,N-1);

uvGrad=dvector(0,N-1); // exGrad is the gradient uv-2 nut Sij
vcGrad=dvector(0,N-1); // 
dd=dvector(0,N-1);
ajex=dvector(0,N-1);

alphaGep=dvector(0,N-1);


  A=dvector(1,N);
  B=dvector(1,N);
  C=dvector(1,N);
  RHS=dvector(1,N);

  U_old=dvector(0,N-1);
  Temperature_old=dvector(0,N-1);
  om_old=dvector(0,N-1);
  k_old=dvector(0,N-1);

  /*EquallySpacedGrid(y,N);*/
  /*StretchedGrid(y,N);*/

  sprintf(InFilename,"input.data");
  ReadRestartData(InFilename,y,U,Temperature,k,om,N);

  for(i=0;i<=N-1;i++){
   if (i<=(int)(N/2)) {
     wallDist[i]=y[i];
   } else{
     wallDist[i]=y[N-1]-y[i];
   }
  }

// make a interface for run_rans_in_loop concept
// input models from a string
// const char *expression = argv[1];
// const char *expression = "1.116+ 0.205*p - 12.*q";
// char *expression = "1.116+ 0.205*p - 12.*q";
// char expression[] = "1.116+ 0.205*p - 12.*q";
  int model_length=1000;

  char expression0[model_length];
  FILE *fpp0=fopen("heat.gep", "r");
  if (fpp0 == NULL) {
    printf("there was an error opening gep.model\n");
    return 1;
  }
// printf("The import models is: ");
  if( fgets (expression0, model_length, fpp0)!=NULL ) {
  /* writing content to stdout */
  //     puts(expression);
    printf("The heat-flux model is %s\n", expression0);
   }
  fclose(fpp0);

  char expression1[model_length];
  FILE *fpp1=fopen("velof1.gep", "r");
  if (fpp1 == NULL) {
    printf("there was an error opening gep.model\n");
    return 1;
  }
  if( fgets (expression1, model_length, fpp1)!=NULL ) {
    printf("The stress-flux f1 model is %s\n", expression1);
   }
  fclose(fpp1);

  char expression2[model_length];
  FILE *fpp2=fopen("velof2.gep", "r");
  if (fpp2 == NULL) {
    printf("there was an error opening gep.model\n");
    return 1;
  }
  if( fgets (expression2, model_length, fpp2)!=NULL ) {
    printf("The stress-flux f2 model is %s\n", expression2);
   }
  fclose(fpp2);

  double p, q;
  te_variable vars[] = {{"p", &p}, {"q", &q}};
  /* This will compile the expression and check for errors. */
  int err;
  te_expr *heat_gep = te_compile(expression0, vars, 2, &err);
  te_expr *velof1_gep = te_compile(expression1, vars, 2, &err);
  te_expr *velof2_gep = te_compile(expression2, vars, 2, &err);



  stop=0;
  iter=0;
  
  while((iter<maxiter) && (stop==0)){
    iter++;
    U[0]=0.0;

    memcpy(U_old,U,N*sizeof(double));
    memcpy(Temperature_old,Temperature,N*sizeof(double));
    memcpy(k_old,k,N*sizeof(double));
    memcpy(om_old,om,N*sizeof(double));
    
    /*
     *Calculating eddy viscosity
     */
    velGrad[0]=(U[1]-U[0])/(y[1]-y[0]);
    CDkom[0]=DMAX(2.0*sigo2/om[0]*((k[1]-k[0])/(y[1]-y[0]))*((om[1]-om[0])/(y[1]-y[0])),1e-10);
    for(i=1;i<=N-2;i++){
      velGrad[i]=((U[i+1]-U[i-1])/(y[i+1]-y[i-1]));
      CDkom[i]=DMAX(2.0*sigo2/om[i]*((k[i+1]-k[i-1])/(y[i+1]-y[i-1]))*((om[i+1]-om[i-1])/(y[i+1]-y[i-1])),1e-10);
    }
    velGrad[N-1]=(U[N-1]-U[N-2])/(y[N-1]-y[N-2]);
    CDkom[N-1]=DMAX(2.0*sigo2/om[N-1]*((k[N-1]-k[N-2])/(y[N-1]-y[N-2]))*((om[N-1]-om[N-2])/(y[N-1]-y[N-2])),1e-10);

    for(i=0;i<=N-1;i++){
      S[i]=sqrt(velGrad[i]*velGrad[i]);
      F2[i]=tanh(pow(DMAX(2.0*sqrt(k[i])/(betStar*om[i]*wallDist[i]),500*nu/(pow(wallDist[i],2.)*om[i])),2.)); //blending function 2
      d[i]=fabs(a1*k[i]/DMAX(a1*om[i],S[i]*F2[i]));  //eddy visocity 
      F1[i]=tanh(pow(DMIN(DMAX(sqrt(k[i])/(betStar*om[i]*wallDist[i]),500*nu/(pow(wallDist[i],2.)*om[i])),4.0*sigo2*k[i]/(CDkom[i]*wallDist[i]*wallDist[i])),4.));
    }

// extra term for GEP implementation
  for(i=0;i<=N-1;i++){
    T[i]=fabs(1/om[i]);
    L[i]=fabs(T[i]*sqrt(k[i])); 
// d[i]=fabs(0.09*k[i]*k[i]/ep[i]);
//    d[i]=1;
  }

  for(i=0;i<=N-1;i++){
	 I1[i]=0.5*velGrad[i]*velGrad[i]*T[i]*T[i];
  }

    tempGrad[0]=(Temperature[1]-Temperature[0])/(y[1]-y[0]);
    for(i=1;i<=N-2;i++){
      tempGrad[i]=((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));
    }
    tempGrad[N-1]=(Temperature[N-1]-Temperature[N-2])/(y[N-1]-y[N-2]);

    for(i=0;i<=N-1;i++){
       J1[i]=tempGrad[i]*tempGrad[i]*L[i]*L[i];
    }
    
 //       for(i=0;i<=N-1;i++){
 //      alphaGep[i]=0.868436 + I1[i]*(2.23575 - 4.45*J1[i]) + J1[i]*(-14.4161 +J1[i]);    // 1e5
//        alphaGep[i]=0.911 + I1[i]* (2.205 - J1[i]) + (-4.444 + J1[i])*J1[i];    
//        alphaGep[i]=0.992 - 6.15*J1[i] - 3.* I1[i]* (-0.616667 + I1[i] + J1[i]);    
//       alphaGep[i]=1. - I1[i]*I1[i] - 2.* (1.089 + I1[i]* (-2. + J1[i]))*(I1[i] - 4.097*J1[i])*(-1. + I1[i] + 0.5*J1[i]);    
 //         alphaGep[i]=1.116 + 0.205*I1[i] - 12.*J1[i]; //2e7
//       alphaGep[i]=1.097 + 3.*(-1.66667 + J1[i])*J1[i];    
//        alphaGep[i]=1.11492 - 6.12727*J1[i] + I1[i]* (-0.36365 + (3.84654 + 10.2245*I1[i] - 2.0449*J1[i])* J1[i]);    
//        }

for(i=0;i<=N-1;i++){
   // alphaGep[i]=1.116 + 0.205*I1[i] - 12.*J1[i]; 
   // f1Gep[i]= 0.205 - 2.43*I1[i] + 4.86 * (-1.57 + (-1. + I1[i]) * I1[i]) *J1[i] ;   
   // f2Gep[i]= -0.57 + J1[i] * (-2.23 + (-6.615 - 3. * J1[i])* J1[i])  ;   
	  p=I1[i];
	  q=J1[i];
    alphaGep[i]=te_eval(heat_gep); 
    f1Gep[i]=te_eval(velof1_gep); 
    f2Gep[i]=te_eval(velof2_gep); 
  }
 
// define uc, vc, vv, uv respectively please.
 //   for(i=0;i<=N-1;i++){
 //     //uv[i]=f1Gep[i]*T[i]*vv[i]*velGrad[i]+f2Gep[i]*gravity*beta*T[i]*vc[i];
 //     //uu[i]=f2Gep[i]*gravity*beta*T[i]*uc[i];
 //     //dd[i]=fabs(f1Gep[i]*T[i]*vv[i]); // modify eddy visocity in a gep velocity models
 //     dd[i]=d[i]; // baseline-cases
 //     exGrad[i] = 0; // baseline case
 //   }
// baseline results  
  for(i=0;i<=N-1;i++){
    // alphaGep[i]=1.111;
    uc[i]=0;
    vc[i]=-1*d[i]*alphaGep[i]*tempGrad[i];

    uu[i]=2/3*k[i];
    // uv[i]=-1*d[i]*f1Gep[i]*velGrad[i];
    uv[i]=-1*d[i]*f1Gep[i]*velGrad[i]+f2Gep[i]*gravity*beta*T[i]*vc[i];
    vv[i]=2/3*k[i];
    dd[i]=d[i]*f1Gep[i]; // baseline-cases
    // ajex[i]=-1*f2Gep[i]*gravity*beta*T[i]*vc[i];
    ajex[i] = uv[i]+dd[i]*velGrad[i]; // baseline case
  }

// uvGrad is the extra term except the (2 nut Sij)
  uvGrad[0]=(ajex[1]-ajex[0])/(y[1]-y[0]);
  for(i=1;i<=N-2;i++){
    uvGrad[i]=((ajex[i+1]-ajex[i-1])/(y[i+1]-y[i-1]));
  }
  uvGrad[N-1]=(ajex[N-1]-ajex[N-2])/(y[N-1]-y[N-2]);

  vcGrad[0]=(vc[1]-vc[0])/(y[1]-y[0]);
  for(i=1;i<=N-2;i++){
    vcGrad[i]=((vc[i+1]-vc[i-1])/(y[i+1]-y[i-1]));
  }
  vcGrad[N-1]=(vc[N-1]-vc[N-2])/(y[N-1]-y[N-2]);


    /********************************
     * k-equation
     ********************************/
    /*
     * Set up boundary problem close to the wall (k=0)
     */
    A[1]=0.0;
    B[1]=1.0;
    C[1]=0.0;
    RHS[1]=0.0;
    for(i=1;i<=N-2;i++){
      /*
       *Setting up coefficients for k-equation
       */
      sigk[i]=sigk1*F1[i]+sigk2*(1.0-F1[i]);

      A[i+1]=(1./((y[i]-y[i-1])*(y[i+1]-y[i-1])))*(2.*nu+(d[i]*sigk[i])+(d[i-1]*sigk[i-1]));
      C[i+1]=(1./((y[i+1]-y[i])*(y[i+1]-y[i-1])))*(2.*nu+(d[i]*sigk[i])+(d[i+1]*sigk[i+1]));
      B[i+1]=-A[i+1]-C[i+1]-(1./dx)-betStar*om[i];
   //   Gb=beta*gravity*d[i]*(k[i]/om[i])*((U[i+1]-U[i-1])/(y[i+1]-y[i-1]))*((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));
    // Gb=(gravity*beta)*(d[i]/Prt)*((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));
    Gb=(gravity*beta)*(d[i]*alphaGep[i])*((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));
      prodkTilde=DMIN(d[i]*pow((U[i+1]-U[i-1])/(y[i+1]-y[i-1]),2.),10*betStar*k[i]*om[i]); // how to deal with om[i] and k[i]
      prod=prodkTilde+fudge*Gb;
      RHS[i+1]=-prod-(k[i]/dx);
    }
    /*
     *Set up boundary problem close to the other wall
     */
    A[N]=0.0;
    B[N]=1.0;
    C[N]=0.0;
    RHS[N]=0.0;

    tridag(A,B,C,RHS,k-1,N); // the reason for -1 is we get k_{i+1}, make the pointer-1 to get k_{i}


    /***************************
     *omega equation
     **************************/
    A[1]=0.0;
    B[1]=1.0;
    C[1]=0.0;
    //RHS[1]=k[1]*2.*nu/pow(y[1],2.);
  RHS[1]=60*nu/(pow(y[1],2.)*bet1); //https://en.wikipedia.org/wiki/SST_(Menter%E2%80%99s_Shear_Stress_Transport)#cite_note-1
//     RHS[1]=2.6221e-17;
//    RHS[1]=10*betStar*k[1]*60*nu/(pow(y[1],2.)*bet1); //https://en.wikipedia.org/wiki/SST_(Menter%E2%80%99s_Shear_Stress_Transport)#cite_note-1
    //RHS[1]=0.0; 

    for(i=1;i<=N-2;i++){
      /*
       *Setting up coefficients for omega-equation
       */
      sigo[i]=sigo1*F1[i]+sigo2*(1.0-F1[i]);
      arg1[i]=2.0*(1.0-F1[i])*sigo2/om[i]*((k[i+1]-k[i-1])/(y[i+1]-y[i-1]))/(y[i+1]-y[i-1]);
      alp[i]=alp1*F1[i]+alp2*(1.0-F1[i]);
      bet[i]=bet1*F1[i]+bet2*(1.0-F1[i]);

      A[i+1]=(1./((y[i]-y[i-1])*(y[i+1]-y[i-1])))*(2.*nu+(d[i]*sigo[i])+(d[i-1]*sigo[i-1]))-arg1[i];
      C[i+1]=(1./((y[i+1]-y[i])*(y[i+1]-y[i-1])))*(2.*nu+(d[i]*sigo[i])+(d[i+1]*sigo[i+1]))+arg1[i];
      B[i+1]=-A[i+1]-C[i+1]-(1./dx)-(bet[i]*om[i]);

      prodkTilde=DMIN(d[i]*pow((U[i+1]-U[i-1])/(y[i+1]-y[i-1]),2.),10*betStar*1.0*k[i]*om[i]); // how to deal with om[i] and k[i]
      //Gb=beta*gravity*d[i]*(k[i]/om[i])*((U[i+1]-U[i-1])/(y[i+1]-y[i-1]))*((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));
      Gb=(gravity*beta)*(d[i]*alphaGep[i])*((Temperature[i+1]-Temperature[i-1])/(y[i+1]-y[i-1]));

      prod=alp[i]*prodkTilde/d[i]+fudge_om*Gb;
      RHS[i+1]=-prod-om[i]/dx;
    }
    A[N]=0.0;
    B[N]=1.0;
    C[N]=0.0;
   // RHS[N]=0.0; 
   //  RHS[N]=2.6221e-17;
   RHS[N]=60*nu/(pow((y[N-1]-y[N-2]),2.)*bet1); 
   // RHS[N]=10*betStar*k[N-2]*60*nu/(pow((y[N-1]-y[N-2]),2.)*bet1); 
   // RHS[N]=k[N-2]*2.*nu/pow(y[N-1]-y[N-2],2.);
    tridag(A,B,C,RHS,om-1,N);

    /**************************************
      Solving for U
      **************************************/
    A[1]=0.0;
    B[1]=1.0;
    C[1]=0.0;
    RHS[1]=0.0;
    for(i=1;i<=N-2;i++){
      A[i+1]=(1./((y[i]-y[i-1])*(y[i+1]-y[i-1])))*(2.*nu+dd[i]+dd[i-1]);
      C[i+1]=(1./((y[i+1]-y[i])*(y[i+1]-y[i-1])))*(2.*nu+dd[i]+dd[i+1]);
      B[i+1]=-A[i+1]-C[i+1]-(1./dx);
      RHS[i+1]=-U[i]/dx-gravity*beta*(Temperature[i]-0.5)+uvGrad[i];
    }
    A[N]=0.0;
    B[N]=1.0;
    C[N]=0.0;
    RHS[N]=0.0;
    tridag(A,B,C,RHS,U-1,N);

    /**************************************
      Solving for Temperature
      **************************************/
    A[1]=0.0;
    B[1]=1.0;
    C[1]=0.0;
    RHS[1]=1.0;
    for(i=1;i<=N-2;i++){
      // A[i+1]=(1./((y[i]-y[i-1])*(y[i+1]-y[i-1])))*(2.*kappa+(d[i]*alphaGep[i]+d[i-1]*alphaGep[i-1]));
      // C[i+1]=(1./((y[i+1]-y[i])*(y[i+1]-y[i-1])))*(2.*kappa+(d[i]*alphaGep[i]+d[i+1]*alphaGep[i+1]));
      A[i+1]=(1./((y[i]-y[i-1])*(y[i+1]-y[i-1])))*(2.*kappa); 
      C[i+1]=(1./((y[i+1]-y[i])*(y[i+1]-y[i-1])))*(2.*kappa);
      B[i+1]=-A[i+1]-C[i+1]-(1./dx);
      RHS[i+1]=-Temperature[i]/dx+vcGrad[i];
    }
    A[N]=0.0;
    B[N]=1.0;
    C[N]=0.0;
    RHS[N]=0.0;
    tridag(A,B,C,RHS,Temperature-1,N);

    L2Norm_U=0.0;
    L2Norm_Temperature=0.0;
    L2Norm_k=0.0;
    L2Norm_om=0.0;
    for(i=0;i<=N-1;i++){
      L2Norm_U+=pow(U[i]-U_old[i],2.);
      L2Norm_Temperature+=pow(Temperature[i]-Temperature_old[i],2.);
      L2Norm_k+=pow(k[i]-k_old[i],2.);
      L2Norm_om+=pow(om[i]-om_old[i],2.);
    }
    L2Norm_U=pow(L2Norm_U/N,0.5);
    L2Norm_Temperature=pow(L2Norm_Temperature/N,0.5);
    L2Norm_k=pow(L2Norm_k/N,0.5);
    L2Norm_om=pow(L2Norm_om/N,0.5);
    
    if(iter%5000==0)
      fprintf(stdout,"%d %e %e %e %e\n",iter,L2Norm_U,L2Norm_Temperature,L2Norm_k,L2Norm_om);

    if((L2Norm_Temperature < TINY) && (L2Norm_U < TINY) && (L2Norm_k < TINY) && (L2Norm_om < TINY))
      stop=1;
  }

  
  utau=pow(nu*fabs((U[1]-U[0])/(y[1]-y[0])),0.5);
  uplus=pow(nu*(U[1]-U[0])/(y[1]-y[0]),0.5);
  ThetaPlus=fabs((kappa/uplus)*(Temperature[1]-Temperature[0])/(y[1]-y[0]));
  fprintf(stderr,"nu=%e  utau=%e thetaplus=%e\n",nu,utau,ThetaPlus);

  fp=fopen("output.dat","w");
  for(i=1;i<=N;i++)
    fprintf(fp,"%e %e %e %e %e %e %e %e %e %e %e %e\n",y[i-1],U[i-1]/sqrt(gravity*beta),Temperature[i-1]-0.5,k[i-1]/(gravity*beta),om[i-1]*nu/(gravity*beta),d[i-1]/sqrt(gravity*beta),\
    y[i-1],U[i-1],Temperature[i-1],k[i-1],om[i-1],d[i-1]);
  fclose(fp);

  sprintf(OutFilename,"LatestGridAndGuess.data");
  WriteGridAndInitialGuessData(OutFilename,y,U,Temperature,k,om,N);

  return 0;
}

// sub-functions 
void EquallySpacedGrid(double *y,int N)
{
  int i;
  for(i=0;i<=N-1;i++)
    y[i]=(2./(N-1))*i;
}

void StretchedGrid(double *y,int N)
{
  int i;
  double x,dx;
  dx=pow(log(2),0.5)/(N/2+1);
  x=0.0;
  y[0]=0.0;

  for(i=1;i<=(N/2)-1;i++){
    x+=dx;
    y[i]=exp(x*x)-1;
  }

  for(i=1;i<=N/2;i++)
    y[N/2+(i-1)]=2.-y[N/2-i];
}

void ReadGridAndInitialGuessData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N)
{
  FILE *fp;
  int i;
  fp=fopen(Filename,"r");
  for(i=0;i<=N-1;i++){
    fscanf(fp,"%lf %lf %lf %lf",&y[i],&U[i],&k[i],&om[i]);
    Temperature[i]=1.0;
    U[i]=sin(2*3.14159265359*y[i]);
  }
  fclose(fp);
}

void ReadRestartData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N)
{
  FILE *fp;
  int i;
  fp=fopen(Filename,"r");
  for(i=0;i<=N-1;i++)
    fscanf(fp,"%lf %lf %lf %lf %lf",&y[i],&U[i],&Temperature[i],&k[i],&om[i]);
  fclose(fp);
}

void WriteGridAndInitialGuessData(char *Filename,double *y,double *U,double *Temperature,double *k,double *om,int N)
{
  FILE *fp;
  int i;
  fp=fopen(Filename,"w");
  for(i=0;i<=N-1;i++)
    fprintf(fp,"%15.12e %15.12e %15.12e %15.12e %15.12e\n",y[i],U[i],Temperature[i],k[i],om[i]);
  
  fclose(fp);
}

double Calc_Utau(double *U,double *y,double nu)
{
  return pow(nu*(U[2]-U[0])/(y[2]-y[0]),0.5);
}

void LinearInterpolate(double *func,double *y,double yp,double *value,int N)
{
  int i;
  i=0;
  while(y[i]<yp)
    i++;

  *value=func[i-1]+((func[i]-func[i-1])/(y[i]-y[i-1]))*(yp-y[i-1]);
}
