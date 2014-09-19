/*******************************************************************/
/* routine to generate a rotation matrix from an axis and an angle */
/* started 1th May 1997, Dave Willock                              */
/*******************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"


void gen_rot_matrix(double *p_quarternion, double *p_axis, double theta)
{

  int iloop,icoloum;

  double theta_on_2, e[4], product[10];

/* Build up the quarternion matrix */

  theta_on_2 = theta/2.0;

  e[0] = cos(theta_on_2); 
  e[1] = (*p_axis)* sin(theta_on_2);
  p_axis++;
  e[2] = (*p_axis)* sin(theta_on_2);
  p_axis++;
  e[3] = (*p_axis)* sin(theta_on_2);

/* arrange quarternion elements :
 *
 *     0  1  2
 *     3  4  5
 *     6  7  8
 *
 * products set as
 *
 *    0 = e0e0,  1 = e1e1,  2 = e2e2,  3 = e3e3,
 *    4 = e0e1,  5 = e0e2,  6 = e0e3,
 *    7 = e1e2,  8 = e1e3,
 *    9 = e2e3
 */

     for (iloop=0; iloop <= 3; iloop++) product[iloop]= e[iloop]*e[iloop];
     for (iloop=4; iloop <= 6; iloop++) product[iloop]= e[0]*e[iloop-3];
     for (iloop=7; iloop <= 8; iloop++) product[iloop]= e[1]*e[iloop-5];
     product[9]= e[2]*e[3];

/* 0 */
     *p_quarternion= product[0] + product[1] - product[2] - product[3];
     p_quarternion++;

/* 1 */
     *p_quarternion= 2.0*(product[7] - product[6]);
     p_quarternion++;

/* 2 */
     *p_quarternion= 2.0*(product[8] + product[5]);
     p_quarternion++;

/* 3 */
     *p_quarternion= 2.0*(product[7] + product[6]);
     p_quarternion++;

/* 4 */
     *p_quarternion= product[0] - product[1] + product[2] - product[3];
     p_quarternion++;

/* 5 */
     *p_quarternion= 2.0*(product[9] - product[4]);
     p_quarternion++;

/* 6 */
     *p_quarternion= 2.0*(product[8] - product[5]);
     p_quarternion++;

/* 7 */
     *p_quarternion= 2.0*(product[9] + product[4]);
     p_quarternion++;

/* 8 */
     *p_quarternion= product[0] - product[1] - product[2] + product[3];
     p_quarternion++;

 return;
}

