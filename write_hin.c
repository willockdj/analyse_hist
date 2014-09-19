#include <stdio.h>
#include "maxima.h"
#include "structures.h"

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms, int which_mol );

void write_atom_data_hin(FILE *p_file, atom *p_atom, int index, int start);


void write_hin(FILE *hin_fp, atom *p_molecule, double *p_abc, int num_of_mols,
               int *p_num_mol_members)
{

atom mol_copy[MAX_ATOMS];

double total_mass, c_of_m[3];

int start_mol, iloop, iatom;
               

/***** TEST WRITING FOR HIN FILE     ******/

   fprintf(hin_fp,"forcefield mm+\nsys 0 0 1\n");
   fprintf(hin_fp,"view 40 0.12234 56.35 16.35 0.9932689");
   fprintf(hin_fp," -0.07986101 0.08389916 0.06445849");
   fprintf(hin_fp," 0.9829017 0.1724799 -0.09623905");
   fprintf(hin_fp," -0.1659109 0.9814335 0.20127 0.42266 -56.419\n");

   fprintf(hin_fp,"box %6.3f  %6.3f  %6.3f\n", *p_abc, *(p_abc+1), *(p_abc+2)); 
   fprintf(hin_fp,"seed -1111\n");

/*********************************************/
/*** Find centre of mass of first molecule ***/
/*********************************************/
   centre_of_mass(&c_of_m[0], &total_mass, p_molecule,
                  (*p_num_mol_members)-1, -1 );

/*****************************************************/
/*** Make copy with centre of mass at screen centre **/
/*****************************************************/

   start_mol=0;
   for (iloop = 0; iloop < num_of_mols; iloop++)
    {
      fprintf(hin_fp,"mol %d\n",iloop+1);
      for (iatom= start_mol; iatom < start_mol+*p_num_mol_members; iatom++)
       {
         mol_copy[iatom]= *p_molecule;
         mol_copy[iatom].x= p_molecule->x - c_of_m[0];
         mol_copy[iatom].y= p_molecule->y - c_of_m[1];
         mol_copy[iatom].z= p_molecule->z - c_of_m[2];
         p_molecule++;
       }

      for (iatom= start_mol; iatom < start_mol+*p_num_mol_members; iatom++)
       {
         write_atom_data_hin(hin_fp, &mol_copy[iatom], iatom, start_mol);
       }

     fprintf(hin_fp,"endmol %d\n",iloop+1);
     start_mol += *p_num_mol_members;
     p_num_mol_members++;
    }

return;
}

