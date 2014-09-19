#include <stdio.h>
#include <math.h>

#include "maxima.h"
#include "ewald.h"
#include "structures.h"

void build_ellipse(atom *p_ellipse_atoms, int num_ellipse_atoms, 
                   double aaa, double bbb, double ccc,
                   double *p_axis1, double *p_axis2, double *p_centre)
{
atom *p_atom;
int iii, iatom;
double xxx[3], yyy[3], centre[3];
double theta, dtheta;

   printf("\nbuilding ellipse using %d atoms\n", num_ellipse_atoms); 
   printf("aaa= %10.6f bbb= %10.6f ccc=%10.6f\n",aaa,bbb,ccc);
   printf("axis 1: %10.6f  %10.6f  %10.6f \n", *p_axis1, *(p_axis1+1), *(p_axis1+2));
   printf("axis 2: %10.6f  %10.6f  %10.6f \n", *p_axis2, *(p_axis2+1), *(p_axis2+2));

   p_atom = p_ellipse_atoms;

   for (iatom=0; iatom < num_ellipse_atoms; iatom++)
     {
        sprintf(p_atom->group,"DRAW");
        sprintf(p_atom->group_no,"1"); 
        sprintf(p_atom->pot,"dr");
        sprintf(p_atom->elem,"H");
        p_atom->part_chge=0.0;   

        p_atom++;
     }

   xxx[0]= *p_axis1;
   xxx[1]= *(p_axis1+1);
   xxx[2]= *(p_axis1+2);

   yyy[0]= *p_axis2;
   yyy[1]= *(p_axis2+1);
   yyy[2]= *(p_axis2+2);

   centre[0] = *p_centre;
   centre[1] = *(p_centre+1);
   centre[2] = *(p_centre+2);

   printf("Have xxx %10.6f %10.6f %10.6f\n", xxx[0], xxx[1], xxx[2]);
   printf("Have yyy %10.6f %10.6f %10.6f\n", yyy[0], yyy[1], yyy[2]);
   printf("Have centre %10.6f %10.6f %10.6f\n", centre[0], centre[1], centre[2]);
   printf("Have aaa %10.6f bbb %10.6f ccc %10.6f\n", aaa, bbb, ccc);

   p_atom = p_ellipse_atoms;
   theta= 0.0;
   dtheta= two_pi/NUM_ELLIPSE_DOTS;
   for (iii=0; iii< NUM_ELLIPSE_DOTS; iii++)
     {
       sprintf(p_atom->label,"%s%d","H",iii+1);

       theta += dtheta;
       p_atom->x=   xxx[0] * aaa * ccc * cos(theta)
                  + yyy[0] * bbb * ccc * sin(theta); 

       p_atom->y=   xxx[1] * aaa * ccc * cos(theta)
                  + yyy[1] * bbb * ccc * sin(theta); 

       p_atom->z=   xxx[2] * aaa * ccc * cos(theta)
                  + yyy[2] * bbb * ccc * sin(theta); 

       p_atom->x += centre[0];
       p_atom->y += centre[1]; 
       p_atom->z += centre[2]; 

       p_atom++;
     }

return;
}
