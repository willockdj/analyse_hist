#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"
#include "mopac_driver.h"

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int find_field(int *p_ichar, int field_desired, int start, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

void read_mopac_out(FILE *file_fp, atom *p_molecule, int num_atoms,
                    double *p_heat, double *p_grad, double *p_reverse_mopac,
                    int num_grads, int have_host, int need_grads )
{

int iloop, good_read, good_second_read, place;
int num_of_chars, num_sec_chars;
int done_heat, done_grad, is_heat, is_grad, is_mullik, done_mullik;
int ichar[LINESIZ];
int skip, index, icomp;

div_t fract;

double *p_this_grad;

char *p_key;

atom *p_atom;

/*********************************************************************/
/**** Read in the bits we want from the MOPAC run ********************/
/*********************************************************************/

good_read= TRUE;
p_atom= p_molecule;

while (good_read)
  {
    num_of_chars= read_line(file_fp, &ichar[0]);

    good_read= num_of_chars != -10;   

    if (good_read)
      {
         p_key= "FINAL HEAT OF FORMATION";
         is_heat= locate_string( p_key, &ichar[0], num_of_chars);
  
         p_key= "FINAL  POINT  AND  DERIVATIVES";
         is_grad= locate_string( p_key, &ichar[0], num_of_chars);

         p_key= "MULLIKEN POPULATIONS AND CHARGES";
         is_mullik= locate_string( p_key, &ichar[0], num_of_chars);

         place=0;
         if (is_heat)
           {
             *p_heat= get_doub( &ichar[0], num_of_chars, &place, &done_heat );
           }
         else if (is_grad)
           {
             num_sec_chars = read_line(file_fp, &ichar[0]);
             num_sec_chars = read_line(file_fp, &ichar[0]);

/**********************************************************************************/
/**** If we have a host then the 6-elements of gradient ignored by MOPAC **********/
/**** will be non-zero due to host...guest interactions but here MOPAC   **********/
/**** will not have written values into them, so we need to put the      **********/
/**** elements in place correctly!! MOPAC gives                          **********/
/****                                atom 1   :  ignr  ignr  ignr        **********/
/****                                atom 2   :  x     ignr  ignr        **********/
/****                                atom 3   :  x     y     ignr        **********/
/****                                atom rest:  x     y     z   etc     **********/
/****                                                                    **********/
/**** We wish to use the molecule as an adsorbate so what zero c of m    **********/
/**** gradient. To do this we add up all current forces, this is the     **********/
/**** current centre of mass force, then remove this from atom 1!!       **********/
/****                                                                    **********/
/**********************************************************************************/

             if (have_host)
               {

                  skip=3; 
                  index=0;
                  p_this_grad= p_grad;

                  for (iloop=0; iloop < num_grads-6; iloop++)                     
                    {
                       num_sec_chars = read_line(file_fp, &ichar[0]); 
                       good_second_read= num_sec_chars != -10;

                       if (good_second_read)
                          {
                             place = find_field(&ichar[0], 7, 0, num_sec_chars);
                             if (skip >0 && iloop != 2) 
                               {
                                  p_this_grad+= skip;
                                  index+= skip;
                                  skip--;
                               }
                             *p_this_grad = get_doub( &ichar[0], num_sec_chars, &place, &done_grad);

                             p_this_grad++;
                             index++;
                          } 
                    }
               }
             else
               {
                  for (iloop=0; iloop < num_grads; iloop++)
                    {
                      num_sec_chars = read_line(file_fp, &ichar[0]);
                      good_second_read= num_sec_chars != -10;
 
                      if (good_second_read)
                        {
                           place = find_field(&ichar[0], 7, 0, num_sec_chars);
                           *p_grad = get_doub( &ichar[0], num_sec_chars, &place, &done_grad);

                           p_grad++;
                        }
                    }
               }
           }
         else if (is_mullik)
           {
              for (iloop=0; iloop < num_atoms; iloop++)
                {
                  num_sec_chars = read_line(file_fp, &ichar[0]);
                  good_second_read= num_sec_chars != -10;

                  if (good_second_read)
                    {
                       place = find_field(&ichar[0], 3, 0, num_sec_chars);
                       p_atom->part_chge= get_doub( &ichar[0], num_sec_chars, &place, &done_mullik);

                       p_atom++;
                    }
                }
           }
      }
  }

fclose(file_fp);

if (!done_heat) 
  {
    printf("Error mopac.out does not contain heat of formation\n");
    exit(0);
  }
if (!done_grad && need_grads) 
  {
    printf("Error mopac.out does not contain Gradiants\n");
    exit(0);
  }
if (!done_mullik) 
  {
    printf("Error mopac.out does not contain mulliken charges\n");
    exit(0);
  }
return;

}

