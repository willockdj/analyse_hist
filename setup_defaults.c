/*******************************************************************/
/********* Routine to initialise variables that are the ************/
/********* same throughout a run, and to set up DEFAULT ************/
/********* values for variables                         ************/
/********* Dave and Dewi began October 1995             ************/
/*******************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "own_maths.h"
#include "ewald.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void setup_defaults(void)
{ 
#include "header.h"
double rdummy, e_charge, epsilon, avos_number;
double arg;

/*** Turn off debugging output by default ***********************/

  DEBUG= FALSE;

  steric = FALSE;
  non_bonded = TRUE;
  charges = TRUE;

/********Mathematical and Physical constants *******************/
 
  one_sixth= 1.0/6.0;
  one_third= 1.0/3.0;
  
  pi = 2.0*acos(0.0);


  two_pi= 2.0*pi;
  four_pi= 4.0*pi;
  four_pi_sqrd= two_pi*two_pi;
  pi_tothehalf = sqrt(pi);

printf ("Imake pi = %10.6f\n",pi);
printf("two_pi= %10.6f\n four_pi=%10.6f\n four_pi_sqrd=%10.6f\n pi_tothehalf=%10.6f\n",
                      two_pi,four_pi,four_pi_sqrd,pi_tothehalf);

  return;
}
