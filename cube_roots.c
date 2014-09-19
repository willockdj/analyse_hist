#include <stdio.h>
#include <math.h>
#include "own_maths.h"

/*********************************************************/
/*** Get roots of a cubic equation if they are real ******/
/*** Formulae from Abramowitz and Stegun handbook   ******/
/*** Formula of form x3 + a[2]x2 + a[1]x + a[0] = 0 ******/
/*** Note that cubic term has unity coeff. and that ******/
/*** coefficients are in order in p_coeffs          ******/
/*** Dave Willock July 1999                         ******/
/*********************************************************/

void cube_roots(double *p_coeffs, double *p_roots)
  {
    double q, r, r2, sr, si, s1ps2, is1ms2;
    double ur, ui, theta;
    double a0, a1, a2;
    double discrim, mag_z;

    a0= *p_coeffs;
    p_coeffs++;
    a1= *p_coeffs;
    p_coeffs++;
    a2= *p_coeffs;

    q= a1/3.0 - (a2 * a2) / 9.0;
    
    r= (a1*a2 - 3.0*a0)/6.0 - a2 * a2 * a2 / 27.0;
    r2 = r*r;

/* printf("In cube roots a2= %10.6f,  a1= %10.6f,  a0= %10.6f, q= %10.6f, r= %10.6f\n",
                           a2, a1, a0, q, r);                                                 */

    discrim= q * q * q + r2;
/* printf("discrim= %10.6f\n", discrim);   */
    if (discrim <= 1.0E-9 )
      {
        discrim= sqrt(-discrim); 
        mag_z= r2 + discrim*discrim;  
        mag_z= sqrt(mag_z);

/**** ur and ui are the real and imaginary parts of r * exp(i(theta)) representation ***/
/**** of s1 and s2                                                                   ***/

        ur = r/mag_z;
        ui = discrim/mag_z;

/*        printf("ur= %10.6f, ui= %10.6f\n", ur, ui);      */

        if (ur <= 1.0 && ur >= -1.0 )
          {
             theta = acos(ur); 
          }
        else
          {
             printf("ur out of range for acos in cube_roots\n");
          }

/*** get cube roots for s1 and s2 *****************************************************/
        ur= cos (theta/3.0);
        ui= sin (theta/3.0);

        mag_z= pow(mag_z, one_third);         

        sr= mag_z*ur;
        si= mag_z*ui;

        s1ps2= 2.0* sr;

        *p_roots= s1ps2 - a2/3.0;
        p_roots++;

        *p_roots= -sr - a2/3.0 - si * sqrt(3.0);
        p_roots++;
 
        *p_roots= -sr - a2/3.0 + si * sqrt(3.0);
         
      }
    else
      {
        printf("Error imaginary roots in cube_roots routine\n");
        /* exit(0); */
      }
  
    return;
  }


