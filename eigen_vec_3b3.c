#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

void unit_vector(double *p_vector);

void vec_cross(double *p_A, double *p_B, double *p_cross);

void eigen_vec_3b3( double *p_matrix, double eigen_val, double *p_eigen_vec, int eigen_no )
{
  int iloop,icoloum, i, j, row_index, column_index;

  double *p_origin_comp;
  double theta, theta_on_2, e[4], quarternion[9], product[10], ut_u[9];
  double mat_u[9], ut_mat_u[9];
  double eee[3], axis[3], vec[3], dot;
  double divider1, divider2, divider3;

/***************************************/
/*** Debug variables *******************/

  double bit1, bit2, bit3, det;

/*** Debug variables *******************/
/***************************************/

/************************************************/
/*** Take eigenvalues off diagonal elements *****/
/************************************************/

  eee[0]= *p_matrix    -eigen_val;
  eee[1]= *(p_matrix+3)-eigen_val;
  eee[2]= *(p_matrix+5)-eigen_val;

/***************************************/
/*** Debug work out determinant ********/
/*  printf("Considering Eigenvalue %10.6f\n", eigen_val);  */

  bit1= eee[1]*eee[2]                - *(p_matrix+4)* *(p_matrix+4);
  bit2= *(p_matrix+1)*eee[2]         - *(p_matrix+2)* *(p_matrix+4);
  bit3= *(p_matrix+1)* *(p_matrix+4) - *(p_matrix+2)* eee[1];

/*  printf("body diag: %10.6f  %10.6f  %10.6f \n",eee[0],  eee[1], eee[2]);        */
/*  printf("Cofactors: %10.6f  %10.6f  %10.6f \n",bit1,bit2,bit3);                 */

  det = eee[0]*bit1 -  *(p_matrix+1) *bit2 + *(p_matrix+2) *bit3;

/*  printf("Determinant = %10.6f\n",det);                                          */

/*** Debug variables *******************/
/***************************************/

/************************************************/
/*** Check for axes already aligned *************/
/************************************************/

  if (fabs(eee[0]) < 1E-6 && fabs(*(p_matrix+1)) < 1E-6 && fabs(*(p_matrix+2)) < 1E-6 )
    {
      *p_eigen_vec     = 1.0;
      *(p_eigen_vec+1) = 0.0;
      *(p_eigen_vec+2) = 0.0;
      return;
    }
  else if (fabs(eee[1]) < 1E-6 && fabs(*(p_matrix+1)) < 1E-6 && fabs(*(p_matrix+4)) < 1E-6 )
    {
      *p_eigen_vec     = 0.0;
      *(p_eigen_vec+1) = 1.0;
      *(p_eigen_vec+2) = 0.0;
      return;
    }
  else if (fabs(eee[2]) < 1E-6 && fabs(*(p_matrix+2)) < 1E-6 && fabs(*(p_matrix+4)) < 1E-6 )
    {
      *p_eigen_vec     = 0.0;
      *(p_eigen_vec+1) = 0.0;
      *(p_eigen_vec+2) = 1.0;
      return;
    }

/**********************************************************************************************/
/*** Check for x-alignment already i.e. avoid using x=1 as start point for axes in zy plane ***/
/**********************************************************************************************/
  if ( fabs(*(p_matrix+1)) > 1E-6 && fabs(*(p_matrix+2))  > 1E-6) 
    {
/*********************************************************/
/*** Solutions based on setting x component to 1 first ***/
/*********************************************************/

      divider1= *(p_matrix+2)*eee[1]- *(p_matrix+1)* *(p_matrix+4);
      divider2= *(p_matrix+2)* *(p_matrix+4) - *(p_matrix+1)*eee[2];
      divider3= *(p_matrix+4)* *(p_matrix+4) - eee[1]*eee[2];

      if ( fabs(eee[1]) > 1E-6 && fabs(divider1) > 1E-6 )
        {
/*          printf("Using divider 1\n");       */
          *p_eigen_vec = 1.0;
          *(p_eigen_vec+2) = ( *(p_matrix+1) * *(p_matrix+1) - eee[0] * eee[1] ) / divider1;
          *(p_eigen_vec+1) = ( -*(p_matrix+4) * *(p_eigen_vec+2) - *(p_matrix+1) ) / eee[1];
        }
      else if ( fabs(*(p_matrix+4)) > 1E-6 && fabs(divider2) > 1E-6 )
        {
/*          printf("Using divider 2\n");        */
          *p_eigen_vec = 1.0;
          *(p_eigen_vec+2) = ( *(p_matrix+1) * *(p_matrix+2) - *(p_matrix+4) * eee[0])/ divider2;
          *(p_eigen_vec+1)=  ( -*(p_matrix+2) - eee[2]* *(p_eigen_vec+2) ) / *(p_matrix+4);
        }
      else if (fabs(*(p_matrix+4)) > 1E-6 && fabs(divider3) > 1E-6 )
        {
/*          printf("Using divider 3\n");        */
          *p_eigen_vec = 1.0;
          *(p_eigen_vec+2) = (*(p_matrix+2) *eee[1] -  *(p_matrix+1) * *(p_matrix+4))/divider3; 
          *(p_eigen_vec+1)=  ( -*(p_matrix+2) - eee[2]* *(p_eigen_vec+2) ) / *(p_matrix+4);
        }
    }
/*********************************************************/
/*** Deal with yz plane cases, set y=1 first! ************/
/*** In the knowledge that x is aligned       ************/
/*********************************************************/

 else if ( fabs(*(p_matrix+4)) > 1E-6 )
    {
      *p_eigen_vec     = 0.0;
      *(p_eigen_vec+2) = 1.0;

      *(p_eigen_vec+1) = -eee[2]/ *(p_matrix+4);
    }
/*********************************************************/
/*** Deal with xy plane cases, set x=1 first! ************/
/*** In the knowledge that z is aligned       ************/
/*********************************************************/
  else if ( fabs(*(p_matrix+1)) > 1E-6 )
    {
      *p_eigen_vec     = 1.0;
      *(p_eigen_vec+2) = 0.0;

      *(p_eigen_vec+1) = -eee[0]/ *(p_matrix+1);
    }

 else
    {
       printf("eigen_vec_3b3 cannot find a solution for the matrix with eigenvalue %10.6f:\n",
                                                                                        eigen_val);
       printf("%10.6f %10.6f %10.6f \n", *p_matrix,     *(p_matrix+1), *(p_matrix+2) );
       printf("%10.6f %10.6f %10.6f \n", *(p_matrix+1), *(p_matrix+3), *(p_matrix+4) );
       printf("%10.6f %10.6f %10.6f \n", *(p_matrix+2), *(p_matrix+4), *(p_matrix+5) );
       printf("Exiting\n");
       exit(0);
    } 

/*  printf("Sending %10.6f %10.6f %10.6f to unit vector\n", 
                         *p_eigen_vec, *(p_eigen_vec+1), *(p_eigen_vec+2));   */

  unit_vector(p_eigen_vec);

/*   printf("Eigenvector for %11.6f : %11.6f %11.6f %11.6f \n", eigen_val, 
                                       *p_eigen_vec, *(p_eigen_vec+1), *(p_eigen_vec+2) ); */

 return;
}

