/***************************************************************************/
/*** Apply a true square matrix operator to a set of vectors with  *********/
/*** a total number of elements given by num_elem (i.e. num_elem/3 *********/
/*** vectors)                                                      *********/
/*** started Dave Willock May 10th 1997                            *********/
/***************************************************************************/
#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void apply_matrix_to_vecs(double *p_vecs, double *p_matrix, int num_elems)
{

  int iloop,icoloum, num_vecs;

  double hold_coord[3], *p_elem;

  num_vecs= (num_elems+1)/3;

printf("DEBUG>> Have %d vectors in apply_matrix_to_vecs\n",num_vecs);
printf("DEBUG>> They have %d elements in total\n",num_elems+1);

  for (iloop= 0; iloop <= num_vecs; iloop++)
  {
     p_elem= p_matrix;
     for (icoloum=0; icoloum< 3; icoloum++)
     {
        hold_coord[icoloum]= 0.0;
        hold_coord[icoloum]+= *p_elem*(*p_vecs);
        p_elem++;
        hold_coord[icoloum]+= *p_elem*(*(p_vecs+1));
        p_elem++;
        hold_coord[icoloum]+= *p_elem*(*(p_vecs+2));
        p_elem++;
     } 
     *p_vecs = hold_coord[0];
     p_vecs++;
     *p_vecs = hold_coord[1];
     p_vecs++;
     *p_vecs = hold_coord[2];
     p_vecs++;
  }

 return;
}

