#include <stdio.h>
#include <time.h>

void matrix_x_ket( double *p_matrix, double *p_a, double *p_b,
                   int dimension );

void ket_bra_to_matrix( double *p_matrix, double *p_ai, double *p_bi,
                        double factor, int dimension );

double bra_x_ket( double *p_a, double *p_b, int dimension );

double bra_mat_ket( double *p_bra, double *p_matrix, double *p_ket, int dimension );

void bfgs_update(double *p_G, double *p_gamma, double *p_delta, int dimension,
                 int triang_length)
{

/***********************************************************************/
/**** Test triangular matrix programs with arrays upto 10x10 ***********/
/***********************************************************************/

double temp_bra[1000], temp_mat[1000];

/***********************************************************************/
/**** Test times for triagular matrix and square matrix ****************/
/***********************************************************************/

double denominator, denominator_2, factor;

int iloop; 

/*************************************************************************/
/***** BFGS update formula:                                  *************/
/*****                                                       *************/
/***** G => G+ (<d|g> + <g|G|g>) |d><d|                      *************/
/*****             |<d|g>| 2                                 *************/
/*****                                                       *************/
/*****       - ( |d><g|G + G|g><d| )                         *************/
/*****              <d|g>                                    *************/
/*************************************************************************/

 denominator   = bra_x_ket(p_delta, p_gamma, dimension);
 denominator_2 =  denominator * denominator;

 matrix_x_ket(p_G, p_gamma, &temp_bra[0], dimension);

 factor = denominator 
        + bra_mat_ket(p_gamma, p_G, p_gamma, dimension);

 factor = 0.5*factor/ denominator_2;

 for (iloop=0; iloop < triang_length; iloop++) temp_mat[iloop] = 0.0; 

 ket_bra_to_matrix(&temp_mat[0], p_delta, p_delta, factor, dimension);

 factor = -1/ denominator;

 ket_bra_to_matrix(p_G, p_delta, &temp_bra[0], factor, dimension); 

/***** Add in first bit **********/

printf("DEBUG>> %d\n",triang_length);
 for (iloop=0; iloop < triang_length; iloop++)
   {
      *p_G += temp_mat[iloop];
      p_G++;
   }

return;
}
