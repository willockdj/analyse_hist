/*************************************************/
/* C Program to Print Output from a HISTORY file */
/*************************************************/

# include <stdio.h>
# include <stdlib.h>
# include "structures.h"

#define MAIN
# include "own_maths.h"
#undef MAIN

void centre_of_mass(double *p_c_of_m, double *p_total_mass, atom *p_molecule,
                    int num_atoms );

void moments_of_inertia(atom *p_molecule, int num_atoms, double *p_c_of_m,
                           double *p_m_of_inertia, double *p_eigenvals );

void eigen_vec_3b3( double *p_matrix, double eigen_val, double *p_eigen_vec );

void write_car( int *p_header_line, int *p_title_line, int *p_date_line,
                atom *p_molecule, int pbc, double *p_abc, int num_atoms,
                int do_header);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void mat_transform( double *p_matrix, double *p_u );

short main ()
{
FILE *fptr;
FILE *fptr_in;
char fname_out[20];
double m_float = 0;
int iloop = 0;
int scan_loop = 0;
char c;
int linecount = 0;
char words[50];
char key[20];
int ikey[20];
int timestep[20];
short num_atoms[20];
float int_timestep[20];
float vectors1[20];
double vectors[20];
char atm_name[20];
short index[20];
float atm_mass[20];
float charge[20];
double velocity[20];
double force[20];
char title_store[50];
int icount;
short i = 0;
atom molecule[1000];

/****************************/
/* Variables for Mom_inert **/
/****************************/

double c_of_m[3], total_mass;
double m_of_inertia[6], moi_eigenvals[3], moi_eigenvecs[9];
double big_eigen, scale, vec[3], temp, abc[6];

int index_new, jloop;
int header_line[100], title_line[100], date_line[100];

atom show_axes[10];

/****************************/
/* Some maths constants *****/
/****************************/

one_third= 1.0/3.0;

/****************************/
/* Gets the Output Filename */
/****************************/

printf ("Enter Output Filename: ");
iloop = 0;
fgets (fname_out, 20, stdin);
while (fname_out[iloop] != '\n') iloop++;
    fname_out[iloop] = '\0';
printf ("Output Filename Entered Was: '%s'\n", fname_out);
printf ("\n");


/**********************************/
/* Opens HISTORY File for Reading */
/**********************************/

if ( (fptr_in = fopen ("HISTORY","r") ) == NULL)
    {
    printf ("ERROR Opening HISTORY FILE\n");
    perror ("open");
    return 1;
    }

/**********************************/
/* Opens OUTPUT File for Writting */
/**********************************/

printf ("Opening HISTORY FILE\n");
if ( (fptr = fopen (fname_out,"w") ) == NULL)
    {
    printf (" ERROR Opening Output File\n"); 
    perror ("open");
    return 2;
    }

/*******************************************************/
/* Scans Information From File, And Then Prints Output */
/*******************************************************/

icount = 0;
while ( (c = fgetc (fptr_in) ) != '\n') 
    {
    title_store[icount] = c;
    icount++; 
    }
title_store[icount] = '\0';
fprintf ( fptr, "%s\n",title_store);
icount = 0;
while ( (c = fgetc (fptr_in) ) != '\n') 
    {
    key[icount] = c;
    icount++; 
    }
key[icount] = '\0';
fprintf ( fptr, "%s\n",key);

    fscanf ( fptr_in, "%s %d %d %d %d %f", &words[0], &timestep[0], &num_atoms[0], &ikey[2], &ikey[3], &int_timestep[0]); 
    fprintf ( fptr, "%s     %d      %d      %d     %d      %f\n", words, timestep[0], num_atoms[0], ikey[2], ikey[3], int_timestep[0]); 
    fscanf ( fptr_in, "%f %le %le", &vectors1[0], &vectors[0], &vectors[1]); 
    fprintf ( fptr, " %.2f      %6.4le  %6.4le\n", vectors1[0], vectors[0], vectors[1]); 
    fscanf ( fptr_in, "%f %f %le", &vectors1[1], &vectors1[2], &vectors[2]); 
    fprintf ( fptr, "%.2f       %.2f      %6.4le\n", vectors1[1], vectors1[2], vectors[2]); 
    fscanf ( fptr_in, "%le %le %f", &vectors[3], &vectors[4], &vectors1[3]); 
    fprintf ( fptr, "%6.4le  %6.4le   %.2f\n", vectors[3], vectors[4], vectors1[3]); 
    i = 0;
    for (i = 0; i < num_atoms[0]; i++)
        {
        fscanf ( fptr_in, "%s %d %lf %lf", &(molecule[i].label[0]), &index[0], &(molecule[i].mass), 
                                         &(molecule[i].part_chge));
        fprintf ( fptr, "%-6s         %3d   %9.6f   %8.7f\n", molecule[i].label, 
                                                              index[0], molecule[i].mass, 
                                                              molecule[i].part_chge);

        fscanf ( fptr_in, "%le %le %le", &(molecule[i].x), &(molecule[i].y), &(molecule[i].z) ); 
        fprintf ( fptr, "%  6.4le %6.4le %6.4le\n", velocity[0], velocity[1], velocity[2]); 
        fscanf ( fptr_in, "%le %le %le", &force[0], &force[1], &force[2]); 
        fprintf ( fptr, "%  6.4le %6.4le %6.4le\n", force[0], force[1], force[2]); 
        }

/***************************/
/* Work out Centre of Mass */
/***************************/

centre_of_mass(&c_of_m[0], &total_mass, &molecule[0], num_atoms[0]-1 );

printf("Calculated:\n");
printf("Centre of mass vector: %10.6f %10.6f %10.6f\n", c_of_m[0], c_of_m[1], c_of_m[2]);
printf("Total Mass           : %10.6f\n", total_mass);
printf("Calculated:\n");

/******************************/
/* Work out Moment of Inertia */
/******************************/

moments_of_inertia(&molecule[0], num_atoms[0]-1, &c_of_m[0],
                           &m_of_inertia[0], &moi_eigenvals[0] );

big_eigen=0.0;
 for (iloop=0; iloop < 3; iloop++)
   {
     eigen_vec_3b3( &m_of_inertia[0], moi_eigenvals[iloop], &moi_eigenvecs[3*iloop] );
     if (moi_eigenvals[iloop] > big_eigen) big_eigen=moi_eigenvals[iloop];
   }

/********************************************************************************/
/*** Put out a show axes molecule as a car file *********************************/
/********************************************************************************/

  scale= 5.0/big_eigen;
  strcpy(show_axes[0].label,"OR");
  strcpy(show_axes[1].label,"X");
  strcpy(show_axes[2].label,"Y");
  strcpy(show_axes[3].label,"Z");
  strcpy(show_axes[0].elem,"H");
  strcpy(show_axes[0].group,"axes");
  strcpy(show_axes[0].group_no,"1");
  strcpy(show_axes[0].pot,"h");
  show_axes[0].part_chge = 0.0;
  show_axes[0].x = 0.0;
  show_axes[0].y = 0.0;
  show_axes[0].z = 0.0;
  jloop=0;
  printf("\n\nMaking up axes:\n");
  for (iloop=1; iloop < 4; iloop++)
   {
     show_axes[iloop].x = scale*moi_eigenvals[iloop-1]*moi_eigenvecs[jloop];
     printf("%10.6f ",moi_eigenvecs[jloop]);
     jloop++;
     show_axes[iloop].y = scale*moi_eigenvals[iloop-1]*moi_eigenvecs[jloop];
     printf("%10.6f ",moi_eigenvecs[jloop]);
     jloop++;
     show_axes[iloop].z = scale*moi_eigenvals[iloop-1]*moi_eigenvecs[jloop];
     printf("%10.6f\n",moi_eigenvecs[jloop]);
     jloop++;
     strcpy(show_axes[iloop].elem,"H");
     strcpy(show_axes[iloop].group,"axes");
     strcpy(show_axes[iloop].group_no,"1");
     strcpy(show_axes[iloop].pot,"h");
     show_axes[iloop].part_chge = 0.0;
   }
  printf("\n\n");

  write_car( &header_line[0], &title_line[0], &date_line[0],
             &show_axes[0], FALSE, &abc[0], 4, TRUE);

/********************************************************/
/*** Write out molecule shifted to c of m ***************/
/********************************************************/
  
  vec[0]= -c_of_m[0];
  vec[1]= -c_of_m[1];
  vec[2]= -c_of_m[2];
  move_molecule(&molecule[0], num_atoms[0]-1, &vec[0]);
 
  write_car( &header_line[0], &title_line[0], &date_line[0],
             &molecule[0], FALSE, &abc[0], num_atoms[0], FALSE);

/* Transpose */
  temp= moi_eigenvecs[3];
  moi_eigenvecs[3]= moi_eigenvecs[1];
  moi_eigenvecs[1]= temp;
  
  temp= moi_eigenvecs[6];
  moi_eigenvecs[6]= moi_eigenvecs[2];
  moi_eigenvecs[2]= temp;

  temp= moi_eigenvecs[7];
  moi_eigenvecs[7]= moi_eigenvecs[5];
  moi_eigenvecs[5]= temp;

 index_new=0;
 printf("Moments of inertia matrix:\n");
 for (iloop=0; iloop < 3; iloop++)
  {
     for (jloop=0    ; jloop < iloop    ; jloop++)
       {
         printf("             ");
       }
     for (jloop=iloop; jloop < 3; jloop++)
       {
         printf("%12.4f ",m_of_inertia[index_new]);
         index_new++;
       }
     printf("\n");
  }

    printf("\nEigenvectors coloumwise\n");
    printf(" %12.4f %12.4f %12.4f \n", moi_eigenvecs[0], moi_eigenvecs[1], moi_eigenvecs[2]);
    printf(" %12.4f %12.4f %12.4f \n", moi_eigenvecs[3], moi_eigenvecs[4], moi_eigenvecs[5]);
    printf(" %12.4f %12.4f %12.4f \n", moi_eigenvecs[6], moi_eigenvecs[7], moi_eigenvecs[8]);

/********************************************************************************/
/*** Trial matrix mult for transforms *******************************************/
/********************************************************************************/

  mat_transform( &m_of_inertia[0], &moi_eigenvecs[0]); 

/***********************/
/* Closes HISTORY File */
/***********************/

if ( (fclose (fptr_in) ) == EOF)
    {
    printf ("%s ERROR Closing Input File\n");
    perror ("close");
    return 3;
    }

/********************************/
/* Closes Specified Output File */
/********************************/

if ( (fclose (fptr) ) == EOF)
    {
    printf ("%s ERROR Closing Output File\n");
    perror ("close");
    return 4;
    }
return 0;
}
