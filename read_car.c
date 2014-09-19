#include <stdio.h>
#include <limits.h>
#include "structures.h"
#include "global_values.h"

/* protype list for this routine */

int read_line(FILE *fp, int *p_ichar);

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int read_atom_data(int *p_ichar, int *p_num_atoms, int at_end, int num_of_chars,
                   int *p_num_of_mols, atom *p_atom );

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

void put_string (int *p_ichar, int length);

/*---------------------------------------------------------------------------*/


/* read in a .car insight formatted file */

int read_car( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
              int *p_num_of_mols, int *p_num_mol_members, 
              double *p_abc, int *p_been_before)

{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end;
  int start,start2,start3,start4;

  char *p_key;

  atom *p_atom;

/* deal with top of file on first pass only */

 if (!*p_been_before)
   {

/* read in first line and check for BIOSYM header */

     *p_been_before= 1;

     num_of_chars = read_line( fp, p_header_line);
     p_key= "BIOSYM";
     is_biosym= locate_string( p_key, p_header_line, num_of_chars);

     if (!is_biosym)
       {
          printf("The BIOSYM header is missing from this file check it really is a .car file \n");
          return 1;
       }

/* check for periodic boundary conditions */

     num_of_chars = read_line( fp, &ichar[0]);
     p_key = "PBC=ON";
     *p_pbc= locate_string( p_key, &ichar[0], num_of_chars);

   }

/* Read in lines common to all frame entries */

/* get title and date lines */

     num_of_chars= read_line ( fp, p_title_line);

/* return Rnegative values (read errors) for parent program to deal with */

     if (num_of_chars < 0) return num_of_chars;

     num_of_chars= read_line ( fp, p_date_line); 


/* if pbc set get abc alpha beta gamma */

   if (*p_pbc == 1) 
     {
        num_of_chars= read_line(fp, &ichar[0]);

        place= -1;
        for (iloop=0; iloop <= 5; iloop++) 
          {
             *(p_abc+iloop)= get_doub(&ichar[0], num_of_chars, &place, &itsanum); 
             if (!itsanum)
	       {
                 printf("failure whilst trying to read in lattice vectors from line: \n");
                 put_string( &ichar[0],100);
                 return -15;
               }
          }
     }

/* Read in the atomic data */

   at_end=0;
   p_atom= p_molecule;
   *p_num_mol_members=0;

   while (at_end != 2) 
     {
       num_of_chars= read_line( fp, &ichar[0]);

       at_end= read_atom_data(&ichar[0],p_num_atoms, at_end, num_of_chars,
                          p_num_of_mols, p_atom);

       if (at_end == 0)
         {
           ++*p_num_mol_members;
           p_atom++;
         }
     
       if (at_end == 1)
         {
           p_num_mol_members++;
           *p_num_mol_members=0;
         }
    }

  *p_num_atoms--;
  *p_num_of_mols--;

  return 0;
}
