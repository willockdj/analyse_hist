#include <stdio.h>
/***************************************************************/
/****** To calculate the value of the determinant |A|     ******/
/****** where A symmetric matrix sent as its top diagonal ******/
/***************************************************************/
/****** Dave Willock June 1999 *********************************/
/***************************************************************/

double calc_det( double *p_matrix, int dimension )
{
int i,j,k;
int index, index_step;

double term, det;

/***************************************************************/
/*** Calculation uses the top left to bottom right products  ***/
/*** minus bottom left to top right products method bottom   ***/
/*** to allow determinants of arbitrary size to be evaluated ***/
/***************************************************************/

/***************************************************************/
/*** top left to bottom right set of products ******************/
/***************************************************************/

for (i=0; i < dimension; i++)
  {
    index=i;
    for (j=0; j < dimension-i; j++)
      {
        if (j!=0) index+= dimension-j+1;
        printf(" %d ", index);
      }
    index= dimension-i;
    j=0;
    for (k=dimension-i; k < dimension; k++)
      {
        printf(" %d ", index);
        j++;
        index+= dimension-j+1; 
      }
    printf("\n");
  }

/***************************************************************/
/*** bottom left to top right set of products   ****************/
/*** Get to central index and then work outward ****************/
/***************************************************************/

for (i=0; i < dimension; i++)
  {
    index= i;
    index_step= dimension-2;
    for (j=0; j <= i; j++)
      {
         printf(" %d ", index);
         index+= index_step;
         index_step-=1;
      }

    printf("\n");
  }



return det;
}
