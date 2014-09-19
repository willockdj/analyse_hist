#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int read_line(FILE *fp, int *p_ichar);

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum );

void read_output(FILE *file_fp, double *p_eng_tot, 
                 double *p_temp_tot,  double *p_eng_cfg, 
                 double *p_eng_vdw,   double *p_eng_cou, 
                 double *p_eng_bnd,   double *p_eng_ang, 
                 double *p_eng_dih,   double *p_eng_tet, 
                 double *p_time,      double *p_eng_pv,   
                 double *p_temp_rot,  double *p_vir_cfg,  
                 double *p_vir_vdw,   double *p_vir_cou,   
                 double *p_vir_bnd,   double *p_vir_ang,   
                 double *p_vir_con,   double *p_vir_tet,   
                 double *p_volume,    double *p_temp_shl,
                 double *p_eng_shl,   double *p_vir_shl,
                 double *p_alpha,     double *p_beta,
                 double *p_gamma,     double *p_vir_pmf,
                 double *p_press,     int *p_num_in_list, 
                 int want_rdf,        rdf *p_rdf,
                 int *p_num_rdfs,     int skip_out,
                 analysis *p_anal_flags )
{

int good_read, num_of_chars, itsanum;
int ichar[LINESIZ], ndigi, sign;

int is_start, i, done_time, done_all, first;
int last_start, place, num, occurances, iloop;
int first_time, done, found, jump;

char *p_key;
char dummy[100];

/*********************************************************************/
/*** Prepare the ground **********************************************/
/*********************************************************************/

last_start= 0;
jump=0;
*p_num_in_list=-1;
first_time= FALSE;
done_time= FALSE;

/*********************************************************************/
/**** Read in the bits we want from the FIELD file *******************/
/*********************************************************************/

done = FALSE;
while ( !done && !done_time)
  {
    while ( strncmp(dummy,"-----",5) != 0 && !done )
                             done = fscanf( file_fp, "%90s",dummy) == EOF;

/*** skip heading characters ***/
    fscanf( file_fp, "%90s", dummy);

    if (strncmp(dummy,"step",4) == 0)
      {
          while ( strncmp(dummy,"-----",5) != 0 && !done )
                             done = fscanf( file_fp, "%90s",dummy) == EOF;

          done = fscanf( file_fp, "%90s", dummy) == EOF;
      }
    else if (strncmp(dummy,"----",4) == 0)
      {
          done = fscanf( file_fp, "%90s", dummy) == EOF;

          while ( strncmp(dummy,"-----",5) != 0 && !done )
                             done = fscanf( file_fp, "%90s",dummy) == EOF;

          done = fscanf( file_fp, "%90s", dummy) == EOF;
      }
    else if (strncmp(dummy,"switching",9) == 0)
      {
          done = fscanf( file_fp, "%90s", dummy) == EOF;

          while ( strncmp(dummy,"-----",5) != 0 && !done )
                             done = fscanf( file_fp, "%90s",dummy) == EOF;

          done = fscanf( file_fp, "%90s", dummy) == EOF;
      }
    else if (strncmp(dummy,"run",3) == 0)
      {
          printf("Finished time data reading\n");
          done_time=TRUE;
      }
    

    if (!done && !done_time) 
      {
        jump++;
        if ( jump == skip_out)
          {
            jump=0;
            (*p_num_in_list)++;

        if (*p_num_in_list > MAX_OUT_LIST)
          {
             printf("ERROR: Data exceeds maximum in read_output\n");
             printf("       Hit OUTPUT line %d but MAX_OUT_LIST = %d\n",
                              *p_num_in_list, MAX_OUT_LIST);
             exit(0);
          }

        if ( p_anal_flags->eng_tot )
          {
            fscanf( file_fp, "%le", p_eng_tot);
            p_eng_tot++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->temp_tot )
          {
            fscanf( file_fp, "%le", p_temp_tot);
            p_temp_tot++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_cfg )
          {
            fscanf( file_fp, "%le", p_eng_cfg);
            p_eng_cfg++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_vdw )
          {
            fscanf( file_fp, "%le", p_eng_vdw);
            p_eng_vdw++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_cou )
          {
            fscanf( file_fp, "%le", p_eng_cou);
            p_eng_cou++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_bnd )
          {
            fscanf( file_fp, "%le", p_eng_bnd);
            p_eng_bnd++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_ang )
          {
            fscanf( file_fp, "%le", p_eng_ang);
            p_eng_ang++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_dih )
          {
            fscanf( file_fp, "%le", p_eng_dih);
            p_eng_dih++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->eng_tet )
          {
            fscanf( file_fp, "%le", p_eng_tet);
            p_eng_tet++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

/*******************************************/
/*** Always Read in the time ***************/
/*******************************************/

        fscanf( file_fp, "%le", p_time);
        p_time++;

        if ( p_anal_flags->eng_pv )
          {
            fscanf( file_fp, "%le", p_eng_pv);
            p_eng_pv++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->temp_rot )
          {
            fscanf( file_fp, "%le", p_temp_rot);
            p_temp_rot++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_cfg )
          {
            fscanf( file_fp, "%le", p_vir_cfg );
            p_vir_cfg++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_vdw )
          {
            fscanf( file_fp, "%le", p_vir_vdw );
            p_vir_vdw++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_cou )
          {
            fscanf( file_fp, "%le", p_vir_cou );
            p_vir_cou++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_bnd )
          {
            fscanf( file_fp, "%le", p_vir_bnd );
            p_vir_bnd++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_ang )
          {
            fscanf( file_fp, "%le", p_vir_ang );
            p_vir_ang++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_con )
          {
            fscanf( file_fp, "%le", p_vir_con );
            p_vir_con++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

        if ( p_anal_flags->vir_tet )
          {
            fscanf( file_fp, "%le", p_vir_tet );
            p_vir_tet++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          }

/**** Don't bother with cpu time *****/

        fscanf( file_fp, "%*s");

        if ( p_anal_flags->volume )
          {
            fscanf( file_fp, "%le", p_volume );
            p_volume++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->temp_shl )
          {
            fscanf( file_fp, "%le", p_temp_shl );
            p_temp_shl++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->eng_shl )
          {
            fscanf( file_fp, "%le", p_eng_shl );
            p_eng_shl++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->vir_shl )
          {
            fscanf( file_fp, "%le", p_vir_shl );
            p_vir_shl++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->alpha )
          {
            fscanf( file_fp, "%le", p_alpha );
            p_alpha++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->beta )
          {
            fscanf( file_fp, "%le", p_beta );
            p_beta++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->gamma )
          {
            fscanf( file_fp, "%le", p_gamma );
            p_gamma++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->vir_pmf )
          {
            fscanf( file_fp, "%le", p_vir_pmf );
            p_vir_pmf++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 

        if ( p_anal_flags->press )
          {
            fscanf( file_fp, "%le", p_press );
            p_press++;
          }   
        else
          {
            fscanf( file_fp, "%*s");
          } 
      }
    }
  }

/*************************************************/
/**** Look for radial distribution functions *****/
/*************************************************/
    if ( want_rdf )
      {
        done= FALSE;
        while ( strncmp(dummy,"RADIAL",6) != 0 && !done )
                             done = fscanf( file_fp, "%90s",dummy) == EOF;

        if ( done )
          {
             printf("ERROR: End of file found before RDF data encountered\n");
             exit(0);
          }
        
/************************************************/
/*** Skip to titles *****************************/
/************************************************/
        
         for (i=1; i<=6; i++) fscanf( file_fp, "%*s");

         done_all = FALSE;
         first= TRUE;
         *p_num_rdfs = -1;
         p_rdf--;

         while (!done_all)
           {
             done = FALSE;

             while ( !done )
              {
                fscanf( file_fp, "%90s",dummy);
                done_all = strcmp(dummy, "time") == 0; 
                done     = done_all || strcmp(dummy, "g(r)") == 0;  

                if (done && !done_all)
                  {
                    p_rdf++;
                    (*p_num_rdfs)++;
                    p_rdf->len = -1;

/**** Read atom types for this rdf remove colon from first **********/

                    fscanf( file_fp, "%90s",dummy);
                    strncpy( p_rdf->atom1, &dummy[1], 4);

                    fscanf( file_fp, "%90s",dummy);
                    strncpy( p_rdf->atom2, dummy, 4);

/*** Skip column headings ***********************/
                    for (i=1; i<=4; i++) fscanf( file_fp, "%90s",dummy);
                  }

                if ( !done_all)
                  {
                    (p_rdf->len)++;
                    p_rdf->r[p_rdf->len] = atof(dummy);

                    fscanf( file_fp, "%le%le", &(p_rdf->g[p_rdf->len]),
                                               &(p_rdf->n[p_rdf->len]));
     
                  }
              }
          }
      }

fclose(file_fp);
return;

}

