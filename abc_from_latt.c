/*********************************************************************************/
/********  Routine to work out                                ********************/
/********  a,b,c,alpha,beta,gamma from cartessian lattice     ********************/
/********  vectors as in p_latt_vec                           ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

#include <math.h>
#include <stdio.h>
#include "constants.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void abc_from_latt( double *p_abc, double *p_latt_vec )
  {
    double alpha, beta, gamma;
    double dot;
    double *p_this_comp, *p_this_abc;
    double *p_a, *p_b, *p_c;

    int iloop;

/********* work out the sizes for a b c ********************************/
    
    p_a= p_latt_vec;
    p_b= p_latt_vec+3;
    p_c= p_latt_vec+6;

    *p_abc     = sqrt(vec_dot(p_a, p_a));
    *(p_abc+1) = sqrt(vec_dot(p_b, p_b));
    *(p_abc+2) = sqrt(vec_dot(p_c, p_c));

/********* work out the angles for alpha, beta and gamma ***************/
 
    alpha = vec_dot(p_b, p_c);
    beta  = vec_dot(p_a, p_c);
    gamma = vec_dot(p_a, p_b);

    *(p_abc+3) = RAD_TO_DEG* acos(alpha/(*(p_abc+1) * *(p_abc+2)));
    *(p_abc+4) = RAD_TO_DEG* acos( beta/(*p_abc     * *(p_abc+2)));
    *(p_abc+5) = RAD_TO_DEG* acos(gamma/(*p_abc     * *(p_abc+1)));

    return;
  }
