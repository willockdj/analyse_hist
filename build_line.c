#include <stdio.h>
#include <math.h>

#include "maxima.h"
#include "ewald.h"
#include "structures.h"

void build_line(atom *p_line_atoms, int num_line_atoms, 
                double aaa, double *p_line, double *p_centre)
{
atom *p_atom;
int iii, iatom;
double xxx[3], centre[3];
double dx, scale;

   printf("\nbuilding line using %d atoms\n", num_line_atoms); 
   printf("aaa= %10.6f \n",aaa);
   printf("direction 1: %10.6f  %10.6f  %10.6f \n", *p_line);

   p_atom = p_line_atoms;

   for (iatom=0; iatom < num_line_atoms; iatom++)
     {
        sprintf(p_atom->group,"DRAW");
        sprintf(p_atom->group_no,"1"); 
        sprintf(p_atom->pot,"dr");
        sprintf(p_atom->elem,"H");
        p_atom->part_chge=0.0;   

        p_atom++;
     }

   xxx[0]= *p_line;
   xxx[1]= *(p_line+1);
   xxx[2]= *(p_line+2);

   centre[0] = *p_centre;
   centre[1] = *(p_centre+1);
   centre[2] = *(p_centre+2);

   p_atom = p_line_atoms;
   dx = aaa/NUM_LINE_DOTS;
   scale=0;
   for (iii=0; iii< NUM_LINE_DOTS; iii++)
     {
       sprintf(p_atom->label,"%s%d","H",iii+1);

       p_atom->x=   xxx[0] * scale;
       p_atom->y=   xxx[1] * scale;
       p_atom->z=   xxx[2] * scale;

       p_atom->x += centre[0];
       p_atom->y += centre[1]; 
       p_atom->z += centre[2]; 

       p_atom++;
       scale += dx;
     }

return;
}
