#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

void unit_vector(double *p_vector);

void vec_cross(double *p_A, double *p_B, double *p_cross);

/* routine to carry out a unitary transform a symmetric matrix using (UT) H U */

void mat_transform( double *p_matrix, double *p_u )
{
  int iloop,icoloum, i, j, row_index, column_index;

  double *p_origin_comp;
  double theta, theta_on_2, e[4], quarternion[9], product[10], ut_u[9];
  double mat_u[9], ut_mat_u[9];
  double eee[3], axis[3], vec[3], dot;

/*******************************************************************/
/*** Written to diagonalise 3x3 matricies with upper diagonal ******/
/*** form on entry  Dave Willock Fri Jul 23 1999              ******/
/*** Works by generating rotation matricies which align the   ******/
/*** coluomns of the matrix with the axes                     ******/
/*******************************************************************/
/***DEBUG TEST *****/
  *p_matrix    = 1;
  *(p_matrix+1)= 1;
  *(p_matrix+2)= 0;
  *(p_matrix+3)= -1;
  *(p_matrix+4)= 0;
  *(p_matrix+5)= 1;
/***DEBUG TEST *****/

  eee[0]= 1.0;
  eee[1]= 0.0;
  eee[2]= 0.0;

  vec[0]=*p_matrix;
  vec[1]=*(p_matrix+1);
  vec[2]=*(p_matrix+2);

  unit_vector(&vec[0]);

  vec_cross(&vec[0], &eee[0], &axis[0]);

  unit_vector(&axis[0]);

/*  printf("Axis: %10.6f %10.6f %10.6f\n", axis[0], axis[1], axis[2]);  */

  theta= acos(axis[0]);

/*  printf("Axis: %10.6f %10.6f %10.6f theta: %10.6f\n", axis[0], axis[1], axis[2], theta*RAD_TO_DEG); */
  
/* Build up the quarternion matrix */

  theta_on_2 = theta/2.0;

  e[0] = cos(theta_on_2); 
  e[1] = axis[0]* sin(theta_on_2);
  e[2] = axis[1]* sin(theta_on_2);
  e[3] = axis[2]* sin(theta_on_2);

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

/* diagonals */

     quarternion[0]= product[0] + product[1] - product[2] - product[3];
     quarternion[4]= product[0] - product[1] + product[2] - product[3];
     quarternion[8]= product[0] - product[1] - product[2] + product[3];

/* upper right off diagonals */

     quarternion[1]= 2.0*(product[7] - product[6]);
     quarternion[2]= 2.0*(product[8] + product[5]);
     quarternion[5]= 2.0*(product[9] - product[4]);

/* lower left off diagonals */

     quarternion[3]= 2.0*(product[7] + product[6]);
     quarternion[6]= 2.0*(product[8] - product[5]);
     quarternion[7]= 2.0*(product[9] + product[4]);

/* Test with multiplication */

    ut_u[0] =   quarternion[0] * quarternion[0]
              + quarternion[1] * quarternion[1]
              + quarternion[2] * quarternion[2];

    ut_u[1] =   quarternion[0] * quarternion[3]
              + quarternion[1] * quarternion[4]
              + quarternion[2] * quarternion[5];

    ut_u[2] =   quarternion[0] * quarternion[6]
              + quarternion[1] * quarternion[7]
              + quarternion[2] * quarternion[8];

    ut_u[4] =   quarternion[3] * quarternion[3]
              + quarternion[4] * quarternion[4]
              + quarternion[5] * quarternion[5];

    ut_u[5] =   quarternion[3] * quarternion[6]
              + quarternion[4] * quarternion[7]
              + quarternion[5] * quarternion[8];

    ut_u[8] =   quarternion[6] * quarternion[6]
              + quarternion[7] * quarternion[7]
              + quarternion[8] * quarternion[8];

    ut_u[3] = ut_u[1];
    ut_u[6] = ut_u[2];
    ut_u[7] = ut_u[5];

/** Print out matrix ***/

/*    printf("quarternion (using as u!)\n");                                           */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n", 
                               quarternion[i], quarternion[i+1], quarternion[i+2]);    */
/*    printf("\n");                                                                    */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n", 
                               ut_u[i], ut_u[i+1], ut_u[i+2]);                         */

/** Apply transformation to original matrix **/

    mat_u[0] =   *p_matrix     * quarternion[0]
               + *(p_matrix+1) * quarternion[3]
               + *(p_matrix+2) * quarternion[6];

    mat_u[1] =   *p_matrix     * quarternion[1]
               + *(p_matrix+1) * quarternion[4]
               + *(p_matrix+2) * quarternion[7];

    mat_u[2] =   *p_matrix     * quarternion[2]
               + *(p_matrix+1) * quarternion[5]
               + *(p_matrix+2) * quarternion[8];

    mat_u[3] =   *(p_matrix+1) * quarternion[0]
               + *(p_matrix+3) * quarternion[3]
               + *(p_matrix+4) * quarternion[6];

    mat_u[4] =   *(p_matrix+1) * quarternion[1]
               + *(p_matrix+3) * quarternion[4]
               + *(p_matrix+4) * quarternion[7];

    mat_u[5] =   *(p_matrix+1) * quarternion[2]
               + *(p_matrix+3) * quarternion[5]
               + *(p_matrix+4) * quarternion[8];

    mat_u[6] =   *(p_matrix+2) * quarternion[0]
               + *(p_matrix+4) * quarternion[3]
               + *(p_matrix+5) * quarternion[6];

    mat_u[7] =   *(p_matrix+2) * quarternion[1]
               + *(p_matrix+4) * quarternion[4]
               + *(p_matrix+5) * quarternion[7];

    mat_u[8] =   *(p_matrix+2) * quarternion[2]
               + *(p_matrix+4) * quarternion[5]
               + *(p_matrix+5) * quarternion[8];

/** Do uT part                              **/

    ut_mat_u[0] =   quarternion[0] * mat_u[0]
                  + quarternion[3] * mat_u[3]
                  + quarternion[6] * mat_u[6];

    ut_mat_u[1] =   quarternion[0] * mat_u[1]
                  + quarternion[3] * mat_u[4]
                  + quarternion[6] * mat_u[7];

    ut_mat_u[2] =   quarternion[0] * mat_u[2]
                  + quarternion[3] * mat_u[5]
                  + quarternion[6] * mat_u[8];

    ut_mat_u[3] =   quarternion[1] * mat_u[0]
                  + quarternion[4] * mat_u[3]
                  + quarternion[7] * mat_u[6];

    ut_mat_u[4] =   quarternion[1] * mat_u[1]
                  + quarternion[4] * mat_u[4]
                  + quarternion[7] * mat_u[7];

    ut_mat_u[5] =   quarternion[1] * mat_u[2]
                  + quarternion[4] * mat_u[5]
                  + quarternion[7] * mat_u[8];

    ut_mat_u[6] =   quarternion[2] * mat_u[0]
                  + quarternion[5] * mat_u[3]
                  + quarternion[8] * mat_u[6];

    ut_mat_u[7] =   quarternion[2] * mat_u[1]
                  + quarternion[5] * mat_u[4]
                  + quarternion[8] * mat_u[7];

    ut_mat_u[8] =   quarternion[2] * mat_u[2]
                  + quarternion[5] * mat_u[5]
                  + quarternion[8] * mat_u[8];

/** Print out matrix ***/
/*    printf("mat_u\n");                                                   */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n",
                               mat_u[i], mat_u[i+1], mat_u[i+2]);          */

/*    printf("ut_mat_u\n");                                                */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n",
                               ut_mat_u[i], ut_mat_u[i+1], ut_mat_u[i+2]); */

exit(0);

 return;
}

