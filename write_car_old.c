
#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"

/* ------Prototype-list---------------------------------------- */

void write_atom_data(atom *p_atom);

void put_string (int *p_ichar, int length);

/* ------------------------------------------------------------ */

void write_car( int *p_header_line, int *p_title_line, int *p_date_line,
                atom *p_molecule, int pbc, double *p_abc, int num_atoms,
                int do_header)

{
   int iloop, this_atom;

   if (do_header)
     {
       put_string(  p_header_line,100);

/* check if periodic boundaries were set */

       if (pbc) 
        {
/*          printf("PBC=ON\n");                   */
          put_string(  p_title_line,100);
          put_string(  p_date_line,100);
/*          printf("PBC");                   */
          for (iloop=0; iloop < 6; iloop++)
           {
/*              printf("%10.4f",*(p_abc+iloop));                   */
           }
/*          printf(" (P1)\n");                   */
        }
       else
        {
/*          printf("PBC=OFF\n");                   */
          put_string(  p_title_line,100);
          put_string(  p_date_line,100);
        }
    }

  for (this_atom=0; this_atom < num_atoms; this_atom++)
   {
      write_atom_data(p_molecule);
      p_molecule++;
   }
/* printf("end\n");                   */

return;
}

