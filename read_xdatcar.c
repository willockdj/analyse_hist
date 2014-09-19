/***************************************************************************/
/*** Read in a VASP XDATCAR file for MD trajectories ***********************/
/*** Dave Willock May 06 ***************************************************/
/*** Modified on Feb 2011 to read in the Cerium atom label from potcar *****/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

int read_line(FILE *fp, int *p_ichar);

int read_atom_data_vasp( FILE *fp, atom *p_atom, int *p_mol_number,
                          coord_flags *p_fix_flags );

double get_double(FILE *input_fp, int skip_lines, int *p_error);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

void put_string (FILE *fp, int *p_ichar, int length);

/*---------------------------------------------------------------------------*/

/* read in a xdatcar vasp formatted file */

int read_xdatcar(FILE *fp, int *p_num_frames, atom *p_molecule, int num_atoms,
                 coord_flags *p_fix_flags, int *p_iconf, int *p_is_fract, int *p_is_cart )
{
  int ichar[LINESIZ];
  int iatom,idummy;
  int skip, lower_case=FALSE;
  int error, good_line;

  char *p_key;
  char *tok;

/* assume xdatcar file is just a read after atoms have been sorted out */
/* On first call the number of frames is assumed to be negative        */

   printf("Trying to read XDATCAR file format\n");

   if (*p_num_frames < 0)
     {
        
/*** Skip title line ****/
        skip = TRUE;
        tok= tok_get(fp, skip, lower_case);
        printf("Read: %s\n", tok);

/*** Skip the integer and the lattice vectors ****/
        idummy = get_integer( fp, skip, &error ); 
        printf("The integer = %d\n",idummy);

/*** Jump lines to set up for first frame read next time in here ***/
        skip=TRUE;
        tok= tok_get(fp, skip, lower_case);
        printf("TOK: %s\n",tok);
        tok= tok_get(fp, skip, lower_case);
        printf("TOK: %s\n",tok);
        tok= tok_get(fp, skip, lower_case);
        printf("TOK: %s\n",tok);
        tok= tok_get(fp, skip, lower_case);
        printf("TOK: %s\n",tok);   
        tok= tok_get(fp, skip, lower_case);
        printf("TOK: %s\n",tok);   

        return 0;
    }
  else
    {
      printf("Back in read_xdatcar to get another frame\n");

/*** Should be left with a line to skip between frames ***/
/*** If this moves to end of file we can stop....      ***/

     skip=TRUE;
     tok= tok_get(fp, skip, lower_case);
     printf("TOK: %s",tok);   

     if (strcmp(tok, END_OF_INPUT) == 0) return FALSE;

     skip=FALSE; 
     tok= tok_get(fp, skip, lower_case);
     printf("  TOK: %s",tok);   
     if (strcmp(tok, "Direct"))
     {
	printf("This is a direct lattice fractional coordinates.\n");
	*p_is_fract=TRUE;
	*p_is_cart=FALSE;
     }
     else
     {
	printf("This is a direct lattice fractional coordinates.\n");
	*p_is_cart=TRUE;
        *p_is_fract=FALSE;
     }
     tok= tok_get(fp, skip, lower_case);
     printf("  TOK: %s\n",tok);   

     *p_iconf= atoi(tok);
     printf("XDATCAR says configuration %d\n", *p_iconf);
     

 /**** Must already have number of frames and am now just reading one in ****/
      skip = TRUE;    
 /**** First skip blank line ***/
/*      while (!(tok= tok_get(fp, skip, lower_case))); */

/*      fscanf(fp,"\n"); */
/*      printf("TOK: %s\n",tok);  */

      for (iatom=0; iatom < num_atoms; iatom++)
        {
          good_line = read_atom_data_vasp( fp, p_molecule, &idummy, p_fix_flags );

          printf("Read: %10.6f  %10.6f  %10.6f good= %d\n", p_molecule->x,
                                                            p_molecule->y,
                                                            p_molecule->z, good_line);
          p_molecule++;
        }
     printf("Frame read so off back to main.....\n");
    }

  return good_line;
}
