#include <stdio.h>
#include <time.h>

void print_product(double *p_matrix, double *p_vector, int dimension);

void matrix_x_ket( double *p_matrix, double *p_a, double *p_b,
                   int dimension );

void ket_bra_to_matrix( double *p_matrix, double *p_ai, double *p_bi,
                        double factor, int dimension );

double bra_x_ket( double *p_a, double *p_b, int dimension );

double bra_mat_ket( double *p_bra, double *p_matrix, double *p_ket, int dimension );

main()
{

/***********************************************************************/
/**** Test triangular matrix programs with arrays upto 10x10 ***********/
/***********************************************************************/

double tri_matrix[55], g_vector[10], delta[10], temp_bra[10];
double temp_mat[55];

/***********************************************************************/
/**** Test times for triagular matrix and square matrix ****************/
/***********************************************************************/

double mat_contents, denominator, denominator_2, factor;

int iloop, jloop, index, dimension;

dimension=3;
mat_contents=1;
index = 0;
for (iloop=0; iloop < dimension; iloop++)
  {
     g_vector[iloop]= iloop+1.0;
     delta[iloop]= 0.0;

     for (jloop=iloop; jloop < dimension; jloop++)
        {
           tri_matrix[index]= mat_contents;
           temp_mat[index]= 0;
           index++;
           mat_contents++;
        }
  }

/*************************************************************************/
/**** Specific example ***************************************************/
/*************************************************************************/

g_vector[0]= -2.0;
g_vector[1]=  1.0;
g_vector[2]= -1.0;

delta[0]= 1.0;
delta[1]=-1.0;
delta[2]= 2.0;

tri_matrix[0]= 1.0;
tri_matrix[1]= 1.0;
tri_matrix[2]=-1.0;
tri_matrix[3]=-2.0;
tri_matrix[4]= 2.0;
tri_matrix[5]= 2.0;

printf("Starting matrix and g_vector\n");
print_product(&tri_matrix[0], &g_vector[0], dimension);

matrix_x_ket(&tri_matrix[0], &g_vector[0], &temp_bra[0], dimension);

printf("delta vector:\n");
for (iloop=0; iloop < dimension; iloop++)
  {
     printf("%d => %10.6f\n", iloop, delta[iloop]);
  }

/**************************************************************************/
/******* n-dimensional dot product ****************************************/
/**************************************************************************/

printf("\n<delta|g_vector> = %10.6f\n", bra_x_ket(&delta[0], &g_vector[0], dimension));

printf("\n<g_vector|M|g_vector> = %10.6f\n",
 bra_mat_ket(&g_vector[0], &tri_matrix[0],  &g_vector[0], dimension));

/*************************************************************************/
/******* Now have all the bits to try dummy update step!!!!! *************/
/*************************************************************************/

/*************************************************************************/
/***** BFGS update formula:                                  *************/
/*****                                                       *************/
/***** G => G+ (<d|g> + <g|G|g>) |d><d|                      *************/
/*****             |<d|g>| 2                                 *************/
/*****                                                       *************/
/*****       - ( |d><g|G + G|g><d| )                         *************/
/*****              <d|g>                                    *************/
/*************************************************************************/

printf("\nInitial matrix:\n\n");
print_product(&tri_matrix[0], &g_vector[0], dimension);

denominator   = bra_x_ket(&delta[0], &g_vector[0], dimension);
denominator_2 =  denominator * denominator;

printf("\ndenominator= %10.6f %10.6f\n",denominator, denominator_2);
matrix_x_ket(&tri_matrix[0], &g_vector[0], &temp_bra[0], dimension);

printf("\nTemp bra :\n");
for (iloop=0; iloop < dimension; iloop++)
  {
     printf("%10.6f\n",temp_bra[iloop]);
  }



factor = denominator 
       + bra_mat_ket(&g_vector[0], &tri_matrix[0],  &g_vector[0], dimension);

factor = 0.5*factor/ denominator_2;

ket_bra_to_matrix(&temp_mat[0], &delta[0], &delta[0], factor, dimension);

printf("\nFirst time round factor= %10.6f\n\n",factor);
print_product(&temp_mat[0], &g_vector[0], dimension);


factor = -1/ denominator;
printf("Second time round factor= %10.6f\n\n",factor);


 ket_bra_to_matrix(&tri_matrix[0], &delta[0], &temp_bra[0], factor, dimension); 


/***** Add in first bit **********/
index = 0;
for (iloop=0; iloop < dimension; iloop++)
  {
     for (jloop=iloop; jloop < dimension; jloop++)
        {
           tri_matrix[index] += temp_mat[index];
           index++;
        }
  }

printf(" Result of BFGS update:- \n");
print_product(&tri_matrix[0], &g_vector[0], dimension);


}
