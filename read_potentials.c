/****************************************************/
/********* Routines for reading .frc files      *****/
/********* begun 27/11/95 Dave Willock          *****/
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "data.h"

void open_file(FILE **p_file, char *p_filename, char *p_status);

int read_line(FILE *fp, int *p_ichar);

int locate_string_in_string( char *p_key, char *p_char, int num_of_chars );
 
int string_from_int(int *p_int, char *p_string);

void get_nonbond_info(FILE *fp, 
                      int *p_line, char *p_line_string,
                      int *p_num_of_chars);

void get_nonbond_params(FILE *fp, int *p_line, char *p_line_string,
                        int *p_num_of_chars);

void get_stretch_params(FILE *fp, int *p_line, char *p_line_string,
                                int *p_num_of_chars, int which);

void get_equivalences(FILE *fp, int *p_line, char *p_line_string,
                                int *p_num_of_chars);

void read_potentials()
{

#include "header.h"

int num_of_chars,iloop;
int ilength1, ilength2;
int line[BUFFER],num_lines;
int grep_int[81], ipoint, iend, itsanum;
int grep_second_int[81];
int idummy, version_found, version_in_use;
int sign, num_digi, ref;
int which_stretch, found_stretch;

double version, a2;

char line_string[BUFFER];
char label[3];
char file_name[81];
char grep_second_string[81];
char non_bond_vdw[BUFFER], bond_stretch[BUFFER], equivalence[BUFFER];

FILE *file_fp;

open_file(&file_fp, "cff91_czeo.frc", "r");

/**************************************************/
/*** Work out what type of forcefield it is *******/
/**************************************************/

num_of_chars=1;
num_lines=0;
version_found=FALSE;

while (num_of_chars != END_OF_FILE && !version_found )
  {
    num_lines++;
    num_of_chars=read_line(file_fp, &line[0]);
    idummy= string_from_int(&line[0], &line_string[0]);

    if ( locate_string_in_string(VERSION, &line_string[0], num_of_chars) )
       {
          if ( locate_string_in_string(PCFF_STRING, &line_string[0], num_of_chars) )  
             {
                fprintf(output_fp,"The potential file contains terms for the pcff potential type.\n");
                version_in_use= PCFF;
                version_found=TRUE;
                strcpy(&non_bond_vdw[0], NON_BOND_VDW_PCFF);
                strcpy(&equivalence[0], EQUIVALENCE_PCFF);
             }
          else if ( locate_string_in_string(CVFF_STRING, &line_string[0], num_of_chars) ) 
             {
                fprintf(output_fp,"The potential file contains terms for the cvff potential type.\n");
                version_in_use= CVFF;
                version_found=TRUE;
                strcpy(&non_bond_vdw[0], NON_BOND_VDW_CVFF);
                strcpy(&equivalence[0], EQUIVALENCE_CVFF);
             } 
          else if ( locate_string_in_string(CFF91_STRING, &line_string[0], num_of_chars) )
             {
                fprintf(output_fp,"The potential file contains terms for the cff91 potential type.\n");
                version_in_use= CFF91;
                version_found=TRUE;
                strcpy(&non_bond_vdw[0], NON_BOND_VDW_CFF91);
                strcpy(&equivalence[0], EQUIVALENCE_CFF91);
             }
          else
             {
                fprintf(output_fp,"ERROR: Cannot interpret potential type from frc file.\n");
                fprintf(output_fp,"Offending line : %s\n", line_string);
                fflush(output_fp);
                fflush(stdout);
             }
       }
  }

/*************************************************************************/
/***** Get the potential information *************************************/
/*************************************************************************/
num_stretches=-1;

while (num_of_chars != END_OF_FILE )
  {
    num_lines++;
    num_of_chars=read_line(file_fp, &line[0]);
    idummy= string_from_int(&line[0], &line_string[0]);

/**************************************************/
/******* Now try and read info from these *********/
/******* lines !!!                        *********/
/**************************************************/
/******* First see if we have one of the  *********/
/******* Myriad of possible stretch pots  *********/
/**************************************************/

    found_stretch= FALSE;
     if ( locate_string_in_string(MORSE_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= MORSE_STRETCH;
        }
     else if (locate_string_in_string(QUARTIC_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= QUARTIC_STRETCH;
        }
     else if (locate_string_in_string(QUADRATIC_STRETCH_STRING, &line_string[0], num_of_chars))
        {
           found_stretch= TRUE;
           which_stretch= QUADRATIC_STRETCH;
        }

    if ( locate_string_in_string(&non_bond_vdw[0], &line_string[0], num_of_chars)  )
      {

/***************************************************/
/***** So we are at the grepped for title  *********/
/***** line !!!                            *********/
/***************************************************/
       
/***************************************************/
/***** Read info lines starting with @ symbols *****/
/***************************************************/

        while (!locate_string_in_string(INFO_LINE, &line_string[0], num_of_chars))
           {
              num_of_chars=read_line(file_fp, &line[0]);
              idummy= string_from_int(&line[0], &line_string[0]);
           }

        get_nonbond_info(file_fp, &line[0], 
                                 &line_string[0], &num_of_chars);  

/*****************************************************/
/***** Now get the values !!! ************************/
/*****************************************************/

        get_nonbond_params(file_fp, &line[0], &line_string[0], &num_of_chars); 

/*****************************************************/
/**** Check the non-bond type is sensible ************/
/*****************************************************/

        if (strcmp(pot_info.type, R_EPS) == 0)
          {
            fprintf(output_fp,"type checks out: potentials in r, eps form\n");
          }
        else if (strcmp(pot_info.type, A_B) == 0)
          {
            fprintf(output_fp,"type checks out: potentials in A, B form\n");
          }
        else
          {
             fprintf(output_fp,"ERROR: Do not understand the non-bond potential type >>%s<< found in the frc file\n",
                                                                                                      pot_info.type);
             fprintf(output_fp,"Currently available types are: r-eps and A-B\n");
             fflush(output_fp);
             fflush(stdout);
             exit(0);
          }

/*****************************************************/
/**** Prepare values required by combining rules *****/
/**** Set bits that aren't required to zero      *****/
/*****************************************************/

       if (strcmp(pot_info.combination, SIXTH_POWER) == 0)
           {
             for  (iloop=0; iloop <= num_potential_types; iloop++)
               {
                  a2 =  potent[iloop].a* potent[iloop].a;
                  potent[iloop].a3 = potent[iloop].a * a2;
                  potent[iloop].a6 = potent[iloop].a3 * potent[iloop].a3;
                  potent[iloop].sqrt_b = sqrt(potent[iloop].b);
                  potent[iloop].sqrt_a = 0;
               }  
           }
       else if (strcmp(pot_info.combination, GEOMETRIC) == 0)
           {
             for  (iloop=0; iloop <= num_potential_types; iloop++)
               {
                 potent[iloop].sqrt_a = sqrt(potent[iloop].a);
                 potent[iloop].sqrt_b = sqrt(potent[iloop].b); 
                 potent[iloop].a3     = 0;
                 potent[iloop].a6     = 0;
               }
           }
       else
           {
             fprintf(output_fp,
                     "ERROR: Do not understand the non-bond potential combining rules >>%s<< found in the frc file\n",
                                                                                                pot_info.combination);
             fprintf(output_fp,"Currently available types are: sixth-power and geometric\n");
             fflush(output_fp);
             fflush(stdout);
           }
      }
/*****************************************************/
/******* Deal with bond stretch param. reading *******/
/*****************************************************/
    else if ( found_stretch )
      {

/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

        while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars)) 
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_stretch_params( file_fp, &line[0], 
                            &line_string[0], &num_of_chars, which_stretch );

      }

/*****************************************************/
/******* Read in the equivalence table ***************/
/*****************************************************/

    else if ( locate_string_in_string(&equivalence[0], &line_string[0], num_of_chars) )
      {
/*****************************************************/
/******* move to the title lines *********************/
/*****************************************************/

         while (!locate_string_in_string(TITLE_LINE, &line_string[0], num_of_chars))
          {
            num_of_chars=read_line(file_fp, &line[0]);
            idummy= string_from_int(&line[0], &line_string[0]);
          }

        get_equivalences( file_fp, &line[0],
                                             &line_string[0], &num_of_chars );

      }
/*****************************************************/
/******* End of titles recognition if blocks *********/
/*****************************************************/

   }
/*****************************************************/
/******* End of parameter reading while loop *********/
/******* Close the file behind you!!!        *********/
/*****************************************************/

fclose(file_fp);
return;
}

