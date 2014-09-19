#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

/*********************************************************************************/
/********  Routine to work out the minimum image vector for   ********************/
/********  periodic systems                                   ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *frac_a, double *frac_b, double *frac_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void min_image( double *x, double *y, double *z)
  {
#include "header.h"
double fract_a, fract_b, fract_c;

/******* find the fractional co-ords of this vector *************/
 
    cart_to_fract( *x, *y, *z, &fract_a, &fract_b, &fract_c, &recip_latt_vec[0]);

/******* Altered to attempt speed up Dave March 98 **************/
/******* now all images in fractional range 0 to 1 **************/

    if (fract_a > 0)
      {
        fract_a= fract_a- floor(fract_a);
      }
    else
      {
        fract_a= fract_a- ceil(fract_a);
      }
    if (fract_a > 0)
      {
        fract_b= fract_b- floor(fract_b);
      }
    else
      {
        fract_b= fract_b- ceil(fract_b);
      }
    if (fract_c > 0)
      {
        fract_c= fract_c- floor(fract_c);
      }
    else
      {
        fract_c= fract_c- ceil(fract_c);
      }

/****** convert back to cartessian ******************************/

    fract_to_cart( x, y, z, fract_a, fract_b, fract_c, &latt_vec[0] );

    return;   

  }
