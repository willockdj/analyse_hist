/***********************************************************************/
/* print_dashes stolen from zebedde April 1997                         */
/*                                                                     */
/* Code by Dewi taken by Dave                                          */
/***********************************************************************/
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
/**************************************************************************/
/* print_dashes                                                           */
/*  prints n_dashes dashes followed by a new line                         */
/**************************************************************************/

void print_dashes(int ndashes,FILE *fp)
	{
#include "header.h"
	int i;
	for (i=0;i<ndashes;i++)
		 fprintf(fp,"-");
    fprintf(fp,"\n");
	return;
}
