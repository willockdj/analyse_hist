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
  double vx;
  double vy;
  double vz;
  double ax;
  double ay;
  double az;
  char pot[4];
  char group[5];
  char group_no[9];
  char elem[3];
  double  part_chge; 
  int nb_list;                  /* index of non-bonding parameters */
  double vdw;			/* van der waals radius of the atom */
  int num_neigh; 		/* number of neighbours */
  int neighb[15]; 		/*initialise -1=not bonded */
  int num_images;               /* Number of symmetry related images added Nov 98 DJW */
  int image[10];                /* image indexes initially allow only 10 images */
  int mol;                      /* molecule index for this atom */
  int neighb_stretch_list[10];   /* index for intra stretch potential */
  double theta;                 /* angle for rotation interpolation */
  double electrostatic_pot;     /* The electrostatic potential at the atom position */
  double mass;                  /* For the atomic mass of the atom */
  double bscat;                 /* For the neutron scattering length for this atom */
}atom;


typedef struct
{
  char atom1[3];
  char atom2[3];
} bond;

typedef struct
{
  char atom1[3];
  char atom2[3];
  char atom3[3];
} angle;

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
char name[7];
char name2[7];
} types;

typedef struct
{
   int mole1;
   char label1[7];
   int index1;
   int mole2;
   char label2[7];
   int index2;
   double sep;
} pair_list;

typedef struct
{
int mo_inertia;
int big_med_small[3];
int contact;
int coords;
int angles;
int dipoles;
int forster;
int add_atom;
int bisector;
int sum_comps;
int correlation;
int eng_tot;
int temp_tot;
int eng_cfg;
int eng_vdw;
int eng_cou;
int eng_bnd;
int eng_ang;
int eng_dih;
int eng_tet;
int time;
int eng_pv;
int temp_rot;
int vir_cfg;
int vir_vdw;
int vir_cou;
int vir_bnd;
int vir_ang;
int vir_con;
int vir_tet;
int volume;
int temp_shl;
int eng_shl;
int vir_shl;
int alpha;
int beta;
int gamma;
int vir_pmf;
int press;
int hin;
int pdb;
int msd;
int bonds;
int rgyr;
int windows;
int drawline;
int drawellipse;
int elli_calc;
int elli_vol;
} analysis;

typedef struct
{
double x;
double y;
double z;
} vector;

typedef struct
{
double totx;
double toty;
double totz;
int num;
} tot_vecs;

typedef struct
{
char atom1[10];
char atom2[10];
double r[MAXRDF];
double g[MAXRDF];
double n[MAXRDF];
int len;
} rdf;

/*** Structures for defining and recording windows                  ***/
/*** num_atoms is the number of atoms expected to define the window ***/
/*** num_windows is the number of windows to expect per molecule    ***/
/***                                                                ***/
/*** In winlist indices of the atoms defining the window for this   ***/
/*** molecule are stored.                                           ***/
/*** Structures for defining and recording windows                  ***/
/*** added March 2011, Dave Willock                                 ***/
/***                                                                ***/
/***                                                                ***/
/***                                                                ***/
typedef struct
{
int num_atoms;
int num_windows;
char atoms[MAX_WINDOW_ATOMS][10];
char neigh[MAX_WINDOW_ATOMS][10];
} windef;

typedef struct
{
int num_windows;
int num_atoms[MAX_WIN_TOT];
int iatom[MAX_WIN_TOT][MAX_WINDOW_ATOMS];
}winlist;

typedef struct
{
double radius[MAX_WIN_TOT];
}win_geometry;

typedef struct
{
double r_vec[MAX_WIN_TOT][3];
double norm[MAX_WIN_TOT][3];
double centre[MAX_WIN_TOT][3];
}win_details;

typedef struct
{
int  fx;
int  fy;
int  fz;
} coord_flags;

typedef struct
{
 int bins;
 char atom1[10];
 char atom2[10];
 char atom3[10];
 char atom4[10];
 double min_limit;
 double max_limit;
 double max_limit1;
 double max_limit2;
 double max_limit3;
} monitors;

/*** monitor information for 2D data, e.g. bond length and angle ****/
typedef struct
{
  int bins1;
  int bins2;
  char atom1[10];
  char atom2[10];
  double min_limit1;
  double max_limit1;
  double min_limit2;
  double max_limit2;
  double delta1;
  double delta2;
} twoD_monitors;

