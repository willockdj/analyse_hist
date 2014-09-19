/********************************************/
/*** Calculate the forster kappa factor   ***/
/*** between two molecules :              ***/
/*** kappa = uD.uA - 3 (uD.rDA)(uD.rDA)   ***/
/*** where uD and uA are unit dipoles     ***/
/*** for donor and acceptor and           ***/
/*** rDA is the donor-acceptor separation ***/
/***  started February 2000 Dave Willock  ***/
/********************************************/
#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

int mol_dipole( atom *p_molecule, int num_atoms, vector *p_mol_dipole );

void unit_vector(double *p_vector);

void forster_kappa(atom *p_donor, int num_donor_atoms, 
                   atom *p_accept, int num_accept_atoms,
                   double *p_kappa)
{
vector u_donor, u_accept;
double usize, uD_dot_uA, uD_dot_rDA, uA_dot_rDA;

double donor_cofm[3], accept_cofm[3], r_DA[3];
double donor_totm, accept_totm;

centre_of_mass( &donor_cofm[0], &donor_totm, p_donor, num_donor_atoms, -1);
centre_of_mass( &accept_cofm[0], &accept_totm, p_accept, num_accept_atoms, -1);

r_DA[0] = accept_cofm[0] - donor_cofm[0];
r_DA[1] = accept_cofm[1] - donor_cofm[1];
r_DA[2] = accept_cofm[2] - donor_cofm[2];

unit_vector(&r_DA[0]);

mol_dipole( p_donor, num_donor_atoms, &u_donor);

usize=  u_donor.x * u_donor.x
      + u_donor.y * u_donor.y
      + u_donor.z * u_donor.z;

usize= sqrt(usize);

u_donor.x = u_donor.x / usize;
u_donor.y = u_donor.y / usize;
u_donor.z = u_donor.z / usize;

mol_dipole( p_accept, num_accept_atoms, &u_accept);

usize=  u_accept.x * u_accept.x
      + u_accept.y * u_accept.y
      + u_accept.z * u_accept.z;

usize= sqrt(usize);

u_accept.x = u_accept.x / usize;
u_accept.y = u_accept.y / usize;
u_accept.z = u_accept.z / usize;

uD_dot_uA   = u_donor.x * u_accept.x
            + u_donor.y * u_accept.y
            + u_donor.z * u_accept.z;

uD_dot_rDA  = u_donor.x * r_DA[0]      
            + u_donor.y * r_DA[1]      
            + u_donor.z * r_DA[2];

uA_dot_rDA  = u_accept.x * r_DA[0]      
            + u_accept.y * r_DA[1]      
            + u_accept.z * r_DA[2];

*p_kappa = uD_dot_uA - 3.0 * uD_dot_rDA * uA_dot_rDA;

return;
}
