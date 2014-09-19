/****************************************************************/
/* circumcircle.c: Define the circumcircle for a set of atoms  **/
/******     started April 2011 Dave Willock                    **/
/****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void unit_vector(double *p_vector);

double size_vector(double *p_vector);

void min_image( double *x, double *y, double *z);

double circumcircle(atom *p_molecule, winlist *p_winsets, int iwin, double *p_centre,
                    double *p_norm, double *p_rvec)
{
#include "header.h"
atom *p_atom1, *p_atom2, *p_atom3;

double vec1[3], vec2[3], dot;
double s12[3], s23[3], s31[3];
double b12[3], b23[3], b31[3];
double p1[3], p2[3], p3[3];
double m12[3], m23[3];
double b12b23, b23v, b12v, radius, lambda;

printf("Calculating circumcircle..%d..\n", iwin);

/**** Atoms 1, 2, 3 are set to define the triangle whose circumcircle we need ***/

p_atom1=p_molecule+ p_winsets->iatom[iwin][0];
p_atom2=p_molecule+ p_winsets->iatom[iwin][1];
p_atom3=p_molecule+ p_winsets->iatom[iwin][2];

/**** vectors defining sides *****/

s12[0]= p_atom2->x - p_atom1->x;
s12[1]= p_atom2->y - p_atom1->y;
s12[2]= p_atom2->z - p_atom1->z;

s23[0]= p_atom3->x - p_atom2->x;
s23[1]= p_atom3->y - p_atom2->y;
s23[2]= p_atom3->z - p_atom2->z;

s31[0]= p_atom1->x - p_atom3->x;
s31[1]= p_atom1->y - p_atom3->y;
s31[2]= p_atom1->z - p_atom3->z;

min_image( &s12[0], &s12[1], &s12[2]);
min_image( &s23[0], &s23[1], &s23[2]);
min_image( &s31[0], &s31[1], &s31[2]);

/*** normal to plane of triangle from cross product ***/
vec_cross(&s12[0], &s23[0], p_norm);
unit_vector(p_norm);

/*** Find vectors perpendicular to sides s12 and s23 ***/
/*** Pointing into triangle                          ***/
vec_cross(p_norm, &s12[0], &b12[0]);
vec_cross(p_norm, &s23[0], &b23[0]);
vec_cross(p_norm, &s31[0], &b31[0]);

unit_vector(&b12[0]);
unit_vector(&b23[0]);
unit_vector(&b31[0]);

/**** m vectors get to the centre point of each side ***/
/**** based on atom 1 being at the origin.           ***/
m12[0] = 0.5*s12[0];
m12[1] = 0.5*s12[1];
m12[2] = 0.5*s12[2];

m23[0] = s12[0]+0.5*s23[0];
m23[1] = s12[1]+0.5*s23[1];
m23[2] = s12[2]+0.5*s23[2];

/*** For similtaneous equations need difference vector ***/
vec1[0]=m23[0]-m12[0];
vec1[1]=m23[1]-m12[1];
vec1[2]=m23[2]-m12[2];

b12b23 = vec_dot( &b12[0], &b23[0]);

vec2[0]= b12[0]- b12b23 * b23[0];
vec2[1]= b12[1]- b12b23 * b23[1];
vec2[2]= b12[2]- b12b23 * b23[2];

if ( 1.0 - fabs(b12b23) < 1.0E-6)
  {
    printf("ERROR: Unlikely co-ordinate set in circumcirle, have to stop.\n");
    exit(0);
  }

lambda = vec_dot( &vec1[0], &vec2[0] ) /  ( 1 - b12b23 *b12b23 );

vec1[0] = m12[0] + lambda * b12[0];
vec1[1] = m12[1] + lambda * b12[1];
vec1[2] = m12[2] + lambda * b12[2];

*p_rvec= vec1[0];
p_rvec++;
*p_rvec= vec1[1];
p_rvec++;
*p_rvec= vec1[2];

*p_centre = p_atom1->x + vec1[0];
 p_centre++;
*p_centre = p_atom1->y + vec1[1];
 p_centre++;
*p_centre = p_atom1->z + vec1[2];

radius= size_vector(&vec1[0]);
/**** Check with distances from three corners ****/
printf("distance from 1 : %10.6f ", radius);

vec2[0]= s12[0]-vec1[0];
vec2[1]= s12[1]-vec1[1];
vec2[2]= s12[2]-vec1[2];

printf("from 2 : %10.6f ", size_vector(&vec2[0]));

vec2[0]= -s31[0]-vec1[0];
vec2[1]= -s31[1]-vec1[1];
vec2[2]= -s31[2]-vec1[2];

printf("and from 3 : %10.6f\n ", size_vector(&vec2[0]));

return radius;
}

		
