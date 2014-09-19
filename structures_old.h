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
  char  label[7];
  double x;
  double y;
  double z;
  char pot[6];
  char group[5];
  char group_no[9];
  char elem[3];
  double  part_chge; 
  int nb_list;                  /* index of non-bonding parameters */
  double vdw;			/* van der waals radius of the atom */
  int num_neigh; 		/* number of neighbours */
  int neighb[5]; 		/*initialise -1=not bonded */
  int neighb_stretch_list[5];   /* index for intra stretch potential */
  double vdw_energy;            /* This atoms vdw energy contribution to the total */
  double electrostatic_pot;     /* The electrostatic potential at the atom position */
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


