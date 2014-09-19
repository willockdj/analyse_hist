/*********************************************************************************/
/********  Routine to work out cartessian real and reciprocal ********************/
/********  lattice vectors from the a,b,c,alpha,beta,gamma    ********************/
/********  line of a biosym .car vector, assuming that the XYZ********************/
/********  convention is followed in the .mdf file.           ********************/
/********                                                     ********************/
/********  p_latt_vec and p_recip_latt_vec point to 9 comp.   ********************/
/********  arrays to be filled                                ********************/
/********  ax ay az bx by bz cx cy cz                         ********************/
/********                                                     ********************/
/********  Adapted Jan 2000 for DLPOLY like files. Here we    ********************/
/********  get latt_vec matrix from the input file so the     ********************/
/********  orientation of abc is not needed logical abc_form  ********************/
/********  used as a switch for this                          ********************/
/********                                                     ********************/
/********  Dave Willock May 1995, Jan 2000                    ********************/
/*********************************************************************************/

#include <math.h>
#include <stdio.h>
#include "constants.h"
#include "own_maths.h"
#include "ewald.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void cart_latt_vecs( double *p_abc, double *p_latt_vec, double *p_real_latt_sizes, 
                     double *p_recip_latt_vec, double *p_recip_latt_sizes,
                     double *p_cell_volume, int num_atoms, int abc_form)
  {
    double bx_cx, mag_a, mag_b, mag_c;
    double a_cross_b[3], b_cross_c[3], c_cross_a[3]; 
    double alpha, beta, gamma;
    double dot_set[9];
    double cell_length;

    int iloop;

    if (abc_form)
      {

/********* make alpha beta gamma in radians ****************************/

         mag_a = *p_abc;
         mag_b = *(p_abc+1);
         mag_c = *(p_abc+2);
         alpha = *(p_abc+3)/RAD_TO_DEG;
         beta  = *(p_abc+4)/RAD_TO_DEG;
         gamma = *(p_abc+5)/RAD_TO_DEG;

/********* X indicates that a is to be aligned with the x axis *********/

         *p_latt_vec    = mag_a;
         *(p_latt_vec+1)= 0.0;
         *(p_latt_vec+2)= 0.0;

/********* XY indicates that the b axis is to lie in the plane formed ******/
/********* by the y-axis and the x-axis i.e. perp. to Z               ******/

         *(p_latt_vec+3)= mag_b * cos( gamma);
         *(p_latt_vec+4)= mag_b * sin( gamma);
         *(p_latt_vec+5)= 0.0;

/******** XYZ indicates that the c vector position is now defined **********/

         *(p_latt_vec+6)= mag_c * cos( beta); 

         bx_cx = ( *(p_latt_vec+3)) * ( *(p_latt_vec+6));

         *(p_latt_vec+7)= ( mag_b*mag_c*cos( alpha)  - bx_cx )/ ( *(p_latt_vec+4));

         *(p_latt_vec+8)= sqrt ( mag_c*mag_c - *(p_latt_vec+6) * *(p_latt_vec+6) 
                                             - *(p_latt_vec+7) * *(p_latt_vec+7)); 

      }

/****** work out the magnitude of each lattice vector **********************/

    *p_real_latt_sizes = *p_latt_vec * *p_latt_vec +
                         *(p_latt_vec+1) * *(p_latt_vec+1) +
                         *(p_latt_vec+2) * *(p_latt_vec+2);

    *p_real_latt_sizes = sqrt(*p_real_latt_sizes);

    *(p_real_latt_sizes+1) = *(p_latt_vec+3) * *(p_latt_vec+3) +
                             *(p_latt_vec+4) * *(p_latt_vec+4) +
                             *(p_latt_vec+5) * *(p_latt_vec+5);

    *(p_real_latt_sizes+1) = sqrt(*(p_real_latt_sizes+1));
   
    *(p_real_latt_sizes+2) = *(p_latt_vec+6) * *(p_latt_vec+6) +
                             *(p_latt_vec+7) * *(p_latt_vec+7) +
                             *(p_latt_vec+8) * *(p_latt_vec+8);

    *(p_real_latt_sizes+2) = sqrt(*(p_real_latt_sizes+2));

/****** If we have matrix version of lattice vectors work out abc **********/

    if (!abc_form)
     {
       *p_abc    = *p_real_latt_sizes;
       *(p_abc+1)= *(p_real_latt_sizes+1);
       *(p_abc+2)= *(p_real_latt_sizes+2);

/* alpha */

       *(p_abc+3) = *(p_latt_vec+3) * *(p_latt_vec+6) +
                    *(p_latt_vec+4) * *(p_latt_vec+7) +
                    *(p_latt_vec+5) * *(p_latt_vec+8);

       *(p_abc+3) = acos(*(p_abc+3)/(*(p_abc+1)* *(p_abc+2)));

/* beta */

       *(p_abc+4) = *p_latt_vec     * *(p_latt_vec+6) +
                    *(p_latt_vec+1) * *(p_latt_vec+7) +
                    *(p_latt_vec+2) * *(p_latt_vec+8);

       *(p_abc+4) = acos(*(p_abc+4)/(*p_abc * *(p_abc+2)));

/* gamma */

       *(p_abc+5) = *p_latt_vec     * *(p_latt_vec+3) +
                    *(p_latt_vec+1) * *(p_latt_vec+4) +
                    *(p_latt_vec+2) * *(p_latt_vec+5);

       *(p_abc+5) = acos(*(p_abc+5)/(*p_abc * *(p_abc+1)));
     } 
/****** generate reciprocal lattice vectors ********************************/

    vec_cross( p_latt_vec,   p_latt_vec+3, &a_cross_b[0]); 
    vec_cross( p_latt_vec+3, p_latt_vec+6, &b_cross_c[0]); 
    vec_cross( p_latt_vec+6, p_latt_vec  , &c_cross_a[0]); 

    *p_cell_volume= vec_dot(p_latt_vec, &b_cross_c[0]); 

    *p_recip_latt_vec    = b_cross_c[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+1)= b_cross_c[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+2)= b_cross_c[2]/ (*p_cell_volume);

    *(p_recip_latt_vec+3)= c_cross_a[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+4)= c_cross_a[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+5)= c_cross_a[2]/ (*p_cell_volume);

    *(p_recip_latt_vec+6)= a_cross_b[0]/ (*p_cell_volume);
    *(p_recip_latt_vec+7)= a_cross_b[1]/ (*p_cell_volume);
    *(p_recip_latt_vec+8)= a_cross_b[2]/ (*p_cell_volume);

/**************** Work out recip vector magnitudes ********************/

    *p_recip_latt_sizes =  *p_recip_latt_vec     * *p_recip_latt_vec
                         + *(p_recip_latt_vec+1) * *(p_recip_latt_vec+1)
                         + *(p_recip_latt_vec+2) * *(p_recip_latt_vec+2);

    *p_recip_latt_sizes = sqrt(*p_recip_latt_sizes);

    *(p_recip_latt_sizes+1) =  *(p_recip_latt_vec+3) * *(p_recip_latt_vec+3)
                             + *(p_recip_latt_vec+4) * *(p_recip_latt_vec+4)
                             + *(p_recip_latt_vec+5) * *(p_recip_latt_vec+5);
 
    *(p_recip_latt_sizes+1) = sqrt(*(p_recip_latt_sizes+1));

    *(p_recip_latt_sizes+2) =  *(p_recip_latt_vec+6) * *(p_recip_latt_vec+6)
                             + *(p_recip_latt_vec+7) * *(p_recip_latt_vec+7)
                             + *(p_recip_latt_vec+8) * *(p_recip_latt_vec+8);

    *(p_recip_latt_sizes+2) = sqrt(*(p_recip_latt_sizes+2));

/**************** Set up Ewald sum parameters *************************/
    
    cell_length= pow(*p_cell_volume, one_third);

    kappa = pi_tothehalf * pow((num_atoms+1), one_sixth)/ cell_length;

/***DEBUG***/
    kappa = kappa*1.75;
/***END DEBUG***/

    four_pi_over_vol= four_pi/(*p_cell_volume);

    pi_sqrd_over_kappa_sqrd= pi*pi/(kappa*kappa);

/***********************************************************************/
/***** factor of pi takes into account that my reciprocal space vecs ***/
/***** are simply a X b / cell_vol etc and Norgett uses 2pi times this */
/***** According to Catlow and Norgett best sum limits are           ***/
/*****           Gm = 2*kappa*f_param and Rm= f_param/kappa          ***/
/***********************************************************************/
   
    recip_sum_max = f_param * kappa / pi;
    real_sum_max  = f_param / kappa;

    recip_sum_max2 = recip_sum_max*recip_sum_max;
    real_sum_max2  = real_sum_max*real_sum_max;
    
/**************** Check by forming all dot products this should *******/
/**************** give the identity matrix ****************************/

    dot_set[0] = vec_dot(p_recip_latt_vec, p_latt_vec);    
    dot_set[1] = vec_dot(p_recip_latt_vec+3, p_latt_vec);    
    dot_set[2] = vec_dot(p_recip_latt_vec+6, p_latt_vec);    

    dot_set[3] = vec_dot(p_recip_latt_vec, p_latt_vec+3);    
    dot_set[4] = vec_dot(p_recip_latt_vec+3, p_latt_vec+3);    
    dot_set[5] = vec_dot(p_recip_latt_vec+6, p_latt_vec+3);    

    dot_set[6] = vec_dot(p_recip_latt_vec, p_latt_vec+6);    
    dot_set[7] = vec_dot(p_recip_latt_vec+3, p_latt_vec+6);    
    dot_set[8] = vec_dot(p_recip_latt_vec+6, p_latt_vec+6);    

    return;
  }
