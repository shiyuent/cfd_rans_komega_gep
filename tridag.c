#define NRANSI
#include "nrutil.h"

void tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n)
{
  unsigned long j;
  double bet,*gam;

  gam=dvector(1,n);
  if (b[1] == 0.0) nrerror("Error 1 in tridag");
  u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
	  gam[j]=c[j-1]/bet;
	  bet=b[j]-a[j]*gam[j];
	  if (bet == 0.0)	nrerror("Error 2 in tridag");
	  u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
	  u[j] -= gam[j+1]*u[j+1];
	free_dvector(gam,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software =$j7. */
