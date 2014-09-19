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
  double theta, theta_on_2, e[4], product[10], ut_u[9];
  double mat_u[9], ut_mat_u[9];
  double eee[3], axis[3], vec[3], dot;

/** Apply transformation to original matrix **/

    mat_u[0] =   *p_matrix     * *p_u
               + *(p_matrix+1) * *(p_u+3)
               + *(p_matrix+2) * *(p_u+6);

    mat_u[1] =   *p_matrix     * *(p_u+1)
               + *(p_matrix+1) * *(p_u+4)
               + *(p_matrix+2) * *(p_u+7);

    mat_u[2] =   *p_matrix     * *(p_u+2)
               + *(p_matrix+1) * *(p_u+5)
               + *(p_matrix+2) * *(p_u+8);

    mat_u[3] =   *(p_matrix+1) * *p_u
               + *(p_matrix+3) * *(p_u+3)
               + *(p_matrix+4) * *(p_u+6);

    mat_u[4] =   *(p_matrix+1) * *(p_u+1)
               + *(p_matrix+3) * *(p_u+4)
               + *(p_matrix+4) * *(p_u+7);

    mat_u[5] =   *(p_matrix+1) * *(p_u+2)
               + *(p_matrix+3) * *(p_u+5)
               + *(p_matrix+4) * *(p_u+8);

    mat_u[6] =   *(p_matrix+2) * *p_u
               + *(p_matrix+4) * *(p_u+3)
               + *(p_matrix+5) * *(p_u+6);

    mat_u[7] =   *(p_matrix+2) * *(p_u+1)
               + *(p_matrix+4) * *(p_u+4)
               + *(p_matrix+5) * *(p_u+7);

    mat_u[8] =   *(p_matrix+2) * *(p_u+2)
               + *(p_matrix+4) * *(p_u+5)
               + *(p_matrix+5) * *(p_u+8);

/** Do uT part                              **/

    ut_mat_u[0] =   *p_u    * mat_u[0]
                  + *(p_u+3) * mat_u[3]
                  + *(p_u+6) * mat_u[6];

    ut_mat_u[1] =   *p_u    * mat_u[1]
                  + *(p_u+3) * mat_u[4]
                  + *(p_u+6) * mat_u[7];

    ut_mat_u[2] =   *p_u    * mat_u[2]
                  + *(p_u+3) * mat_u[5]
                  + *(p_u+6) * mat_u[8];

    ut_mat_u[3] =   *(p_u+1) * mat_u[0]
                  + *(p_u+4) * mat_u[3]
                  + *(p_u+7) * mat_u[6];

    ut_mat_u[4] =   *(p_u+1) * mat_u[1]
                  + *(p_u+4) * mat_u[4]
                  + *(p_u+7) * mat_u[7];

    ut_mat_u[5] =   *(p_u+1) * mat_u[2]
                  + *(p_u+4) * mat_u[5]
                  + *(p_u+7) * mat_u[8];

    ut_mat_u[6] =   *(p_u+2) * mat_u[0]
                  + *(p_u+5) * mat_u[3]
                  + *(p_u+8) * mat_u[6];

    ut_mat_u[7] =   *(p_u+2) * mat_u[1]
                  + *(p_u+5) * mat_u[4]
                  + *(p_u+8) * mat_u[7];

    ut_mat_u[8] =   *(p_u+2) * mat_u[2]
                  + *(p_u+5) * mat_u[5]
                  + *(p_u+8) * mat_u[8];

/** Print out matrix ***/
/*    printf("mat_u\n");                                                     */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n",
                               mat_u[i], mat_u[i+1], mat_u[i+2]);            */

/*    printf("ut_mat_u\n");                                                  */
/*    for ( i=0; i<9; i+=3) printf("%10.6f  %10.6f  %10.6f\n",
                               ut_mat_u[i], ut_mat_u[i+1], ut_mat_u[i+2]);   */

 return;
}

