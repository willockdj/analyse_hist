/******************************************************/
/**** Routine to give values of the error function ****/
/**** Dave Willock May 7th 1996                    ****/
/******************************************************/

#include <stdio.h>
#include <math.h>

double erf_mac(double x, double accuracy);
 
double erf_asymptotic(double x, double accuracy);
 
double own_erf(double x, double accuracy)
{
#include "own_maths.h"
double erf_ans;

erf_ans = erf(x);

/*
if (x < 4.0)
  {
erf_ans = erf_mac(x, accuracy);
  }
else 
  {
erf_ans = erf_asymptotic(x, accuracy);
  }

printf ("Them: %10.6f Us : %10.6f\n", erf_ans2, erf_ans);
*/
return erf_ans;
}
