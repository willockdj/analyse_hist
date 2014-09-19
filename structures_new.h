/* structres.h */
/* Mine first then those from orig (tim) after for reference */

typedef struct
{
 int steric;
 double non_bonded;			        /* non_bonding energy */
 double vdw_rep;                                /* van der Waals repulsive */
 double vdw_disp;                               /* van der Waals dispersive */
 double charges;				/* coulombic energy */
 double guest_guest;                            /* part of total energy due to guest-guest interactions */
 double guest_host;                             /* part of total energy due to guest-host  interactions */
 int acceptance;
 double minimizer_init_total;			/* minimizer non_bonding energy */
 double minimizer_end_total;			/* minimizer total energy */
 double minimizer_init_nonbond;			/* minimizer non_bonding energy */
 double minimizer_end_nonbond;			/* minimizer total energy */
} energy;

/* rudimentary statistics structure */
 typedef struct
 {
  int tries;
  int tries_after_build;
  int accepted;
  int accepted_after_build;
 }stats;

typedef struct
{
double stretch;
} internal_energy;

 
/* structure for the atoms */

typedef struct
{
  double x;
  double y;
  double z;
}atom;

typedef struct
{
  char atom1[3];
  char atom2[3];
} bond;

typedef struct
{
  char atom_type[3];
  int  num;
} atom_number;

typedef struct
{
  int start;
  int end;
} links;

typedef struct
{
double matrix[9];
double translation[3];
} symm_ops;

typedef struct
{
int start;
int end;
int num;
} list_partition;

typedef struct
{
char name[3];
} types;

typedef struct
{
   int mole1;
   char label1[7];
   int index1;
   int mole2;
   char label2[7];
   int index2;
} pair_list;


