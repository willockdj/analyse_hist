/*********************************************************************/
/*** write_car will now just write out the co-ords it is sent   ******/
/*** without alterations. Any messing with the structure should ******/
/*** be done elsewhere! Dave Willock November 2005              ******/
/*********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "maxima.h"
#include "ewald.h"
#include "constants.h"
#include "structures.h"
#include "global_values.h"

#define DEBUG FALSE

/* ------Prototype-list---------------------------------------- */

void write_atom_data(FILE *fp, atom *p_atom, double scale_factor, coord_flags *p_fix_flags);

void put_string (FILE *fp, int *p_ichar, int length);

void min_image( double *x, double *y, double *z, double *p_recip_latt_vec, double *p_latt_vec);

void cart_to_fract( double cart_x,  double cart_y,  double cart_z, 
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z, 
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void vec_cross(double *p_A, double *p_B, double *p_cross);

double size_vector(double *p_vector);

double unit_vector(double *p_vector);
/* ------------------------------------------------------------ */

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags,
                int need_draw, win_details *p_win_geom, int num_windows ) 

{
   int iii, iwin, iloop, mol_current, this_atom;
   int iavec, ibvec, icvec;
   double ta[3],tb[3],t[3], x,y,z;
   double fract_a, fract_b, fract_c;
   double rrr, yyy[3], theta, dtheta;
   char *p_this_char;

   atom current_atom;
   atom draw_atom;
   atom *p_atom;

   coord_flags *p_this_fix;

/*** Drawing variables ***/
   int idraw;

/*** Use the first element of header_line equal -1 as an indicator that */
/*** the header line has not been read                                  */

   if (start_frame)
     {
        if (*p_header_line != -1)
          {
             put_string( fp, p_header_line,100);
          }
        else
          {
             fprintf(fp, "!BIOSYM archive 3\n");
          }
     }

/* check if periodic boundaries were set */

   if (pbc) 
    {
      if (start_frame) fprintf(fp, "PBC=ON\n");

/*
      if (*p_title_line != -1)
        {
          put_string(fp, p_title_line,100);
          fprintf(fp, "\n");
        }
      else
        {
          fprintf(fp, "%s", p_c_title_line);
        }
*/

     fprintf(fp, "analyse_hist generated file\n");

/*** Use the first element of date_line equal -1 as an indicator that */
/*** the header line has not been read                                  */

   if (*p_date_line != -1)
     {
      put_string(fp, p_date_line,100);
     }
   else
     {
      fprintf(fp, "!DATE Mon Oct 12 12:17:26 1998\n");
     }

      fprintf(fp,"PBC");

      for (iloop=0; iloop < 6; iloop++)
       {
         if ( iloop <= 2 )
           {
             fprintf(fp,"%10.4f",*(p_abc+iloop) * scale_factor * *(p_super+iloop));
           }
         else
           {
             fprintf(fp,"%10.4f",*(p_abc+iloop)*RAD_TO_DEG);
           }
       }
      fprintf(fp," (P1)\n");
    }
   else
    {
      if (start_frame) fprintf(fp,"PBC=OFF\n");
      put_string(fp, p_title_line,100);
      put_string(fp, p_date_line,100);
    }

    if (need_draw)
      {
        idraw=0;
      }

    if (pbc)
      {
        mol_current=0;

        p_atom=p_molecule;
        for (iavec=0; iavec<= (*p_super)-1; iavec++)
          {
            ta[0]=iavec * *p_latt_vec;
            ta[1]=iavec * *(p_latt_vec+1);
            ta[2]=iavec * *(p_latt_vec+2);

          for (ibvec=0; ibvec<=*(p_super+1)-1; ibvec++)
            {
              tb[0]=ibvec * *(p_latt_vec+3);
              tb[1]=ibvec * *(p_latt_vec+4);
              tb[2]=ibvec * *(p_latt_vec+5);

              for (icvec=0; icvec<=*(p_super+2)-1; icvec++)
                {
                  t[0]= ta[0] + tb[0] + icvec * *(p_latt_vec+6);
                  t[1]= ta[1] + tb[1] + icvec * *(p_latt_vec+7);
                  t[2]= ta[2] + tb[2] + icvec * *(p_latt_vec+8);

                  p_atom=p_molecule;
                  p_this_fix= p_fix_flags;
                  for (this_atom=0; this_atom < num_atoms; this_atom++)
                    {

                       if (*(p_mol_number+this_atom) != mol_current)
                         {
                            mol_current++;
                            fprintf(fp, "end\n");
                         }

                       current_atom = *p_atom;
                       current_atom.x += t[0];
                       current_atom.y += t[1];
                       current_atom.z += t[2];

                       if (DEBUG) printf("writing data for atom at %10.6f %10.6f %10.6f\n",
                                          current_atom.x,current_atom.y,current_atom.z);

                       write_atom_data(fp, &current_atom, scale_factor, p_this_fix );
                       p_this_fix++;
                       p_atom++;
                   }
/*** Add any features requested by drawing options ****/
                if (need_draw)
                  {   
/*** DEBUG, Check centres *********/
                     sprintf(current_atom.group,"DRAW");
                     sprintf(current_atom.group_no,"1"); 
                     sprintf(current_atom.pot,"dr");
                     sprintf(current_atom.elem,"S");
                     current_atom.part_chge=0.0;   
                     p_this_fix= p_fix_flags;
                     for (iwin=0; iwin<=num_windows; iwin++)
                       {
                         sprintf(current_atom.label,"%s%d","S",iwin);
                         current_atom.x= p_win_geom->centre[iwin][0];
                         current_atom.y= p_win_geom->centre[iwin][1];
                         current_atom.z= p_win_geom->centre[iwin][2];
                         write_atom_data(fp, &current_atom, scale_factor, p_this_fix ); 
                       }
/*** DEBUG, Check centres end *****/
                     printf("Drawing circles\n"); 
                     sprintf(current_atom.group,"DRAW");
                     sprintf(current_atom.group_no,"1"); 
                     sprintf(current_atom.pot,"dr");
                     sprintf(current_atom.elem,"H");
                     current_atom.part_chge=0.0;   
                     p_this_fix= p_fix_flags;

                     for (iwin=0; iwin<=num_windows; iwin++)
                       {
/*** Form a y axis in the plane ***/
                         vec_cross(&(p_win_geom->r_vec[iwin][0]),  
                                   &(p_win_geom->norm[iwin][0]),   
                                   &yyy[0]); 

                         rrr= size_vector(&(p_win_geom->r_vec[iwin][0]));
                         unit_vector(&yyy[0]);

                         printf("yyy: %10.6f %10.6f %10.6f rrr: %10.6f\n",
                                  yyy[0],yyy[1],yyy[2],rrr);  

                         theta= 0.0;
                         dtheta= two_pi/20;
/*** Don't bother with circles just yet ***/
                         for (iii=0; iii<20; iii++) 
                           {
                             idraw++;
                             sprintf(current_atom.label,"%s%d","H",idraw);

                             theta += dtheta;
                             current_atom.x= p_win_geom->centre[iwin][0]
                                          + p_win_geom->r_vec[iwin][0]* sin(theta)
                                          + rrr*yyy[0]*cos(theta); 

                             current_atom.y=  p_win_geom->centre[iwin][1]
                                          + p_win_geom->r_vec[iwin][1]* sin(theta)
                                          + rrr*yyy[1]*cos(theta);  

                             current_atom.z=   p_win_geom->centre[iwin][2]
                                          + p_win_geom->r_vec[iwin][2]* sin(theta)
                                          + rrr*yyy[2]*cos(theta);   

                             write_atom_data(fp, &current_atom, scale_factor, p_this_fix ); 
                           }
                         printf("All done for window %d\n",iwin);
                       }
                  }
              }
          }
      }
  }
else
   {
       mol_current=0;
       p_atom=p_molecule;
       p_this_fix= p_fix_flags;
       for (this_atom=0; this_atom < num_atoms; this_atom++)
         {

            if (*(p_mol_number+this_atom) != mol_current)
              {
                 mol_current++;
                 fprintf(fp, "end\n");
              }

            write_atom_data(fp, p_atom, scale_factor, p_this_fix );
            p_this_fix++;
            p_atom++;
        }

   }

fprintf(fp,    "end\nend\n");

return;
}

