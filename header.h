/************************************/
/* THIS IS MY STUFF                 */
/************************************/

#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

/********************** read variables********************/

EXTERNAL int read_new_line;
EXTERNAL int line_no;


/********************** files pointers********************/

EXTERNAL FILE *input_fp;
EXTERNAL FILE *output_fp;
EXTERNAL FILE *pore_fp; /* pore coordinate file */
EXTERNAL FILE *adsorbed_fp; /* adsorbate coordinate file */
EXTERNAL FILE *fragment_fp; /* fragment coordinate file */
EXTERNAL FILE *template_fp; /* good template file pointer */
EXTERNAL FILE *strategy_fp; /* strategy file for discover */
EXTERNAL FILE *forcefield_fp; /* forcefield_file for discover */
EXTERNAL FILE *discover_fp; /* discover files (car/mdf) for minimizer */
EXTERNAL FILE *gooduns_fp;  /* accepted templates arc file */

/******* Animation file pointers **********************************/

EXTERNAL FILE *anim_file_fp;
EXTERNAL FILE *anim_read_fp;
EXTERNAL FILE *anim_show_fp;

EXTERNAL FILE *gooduns_read_fp;
EXTERNAL FILE *gooduns_show_fp;

/************ files and temporary file read variables ***************/

EXTERNAL char inputfile[FILELEN_MAX];
EXTERNAL char adsorbed_file[FILELEN_MAX];
EXTERNAL char fragment_file[FILELEN_MAX];
EXTERNAL char template_file[FILELEN_MAX];
EXTERNAL char template_strategy_file[FILELEN_MAX];
EXTERNAL char  inpore_strategy_file[FILELEN_MAX];
EXTERNAL char gooduns[FILELEN_MAX];
EXTERNAL char defaults_file[FILELEN_MAX];

EXTERNAL char discover_path[256];            /* path for discover */
EXTERNAL char template_min_car[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char template_min_mdf[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char inpore_min_car[FILELEN_MAX]; /* for intermediate discover use */
EXTERNAL char inpore_min_mdf[FILELEN_MAX]; /* for intermediate discover use */

EXTERNAL char discover_forcefield_name[FILELEN_MAX];
EXTERNAL char forcefield_library[FILELEN_MAX]; /*user potentials file */

/****************************************************************************/
/********** Variables used to keep track of potentials **********************/
/****************************************************************************/

EXTERNAL int num_potential_types;
EXTERNAL int num_stretches;
EXTERNAL int num_equivalences;

EXTERNAL char minimizer_name[20];      	/* name of minimization tech to use */
EXTERNAL char title[80];
EXTERNAL char buffer[256];
EXTERNAL char dummy_head[256];   /* temp storage for biosym headers */
EXTERNAL char pore_title[256];
EXTERNAL char seed_title[256];
EXTERNAL char fragment_title[256];


/* command line for mopac from strategy file*/

EXTERNAL char mopac_cmdline_molecule[256];  
EXTERNAL char mopac_cmdline_inpore[256]; 
EXTERNAL char mopac_root[256];
EXTERNAL char mopac_path[256];  /* path for mopac*/

#ifdef MAIN
EXTERNAL int max_mopac_atoms =99999;  /*silly until we work out what the max is*/
#else
EXTERNAL int max_mopac_atoms;  
#endif

/********************* pore/template/fragment structures*********************/

EXTERNAL atom_number atom_limit[NUM_ELEMENTS];/* struct with max elem. concentration*/
EXTERNAL int have_forbidden_bonds;		/* logical for forbidden bonds */
EXTERNAL int num_forbidden_bonds;			/* number of forbidden bonds */
EXTERNAL int have_conc_limits; 			/* logical for concentration limit */
EXTERNAL int num_conc_limits; 			/* number of   concentration limit */
EXTERNAL links link_atoms[MAX_ATOMS];      /* struct with frag-frag bond list */
EXTERNAL bond forbidden_bond[MAX_ATOMS];   /* struct for forbidden bonds */
EXTERNAL int action_weights[10];  /* weights for action selection */
EXTERNAL int num_links;                 /* number of frag-frag bonds */
EXTERNAL int num_frag_weights;
EXTERNAL int num_action_weights;
EXTERNAL int sum_action_weights;
EXTERNAL int sum_frag_weights;
EXTERNAL int frag_weights_given;
EXTERNAL int action_weights_given;


/*********************periodic structure info ***********************/

/********************************************************************/
/***** pbc    : flags that periodic boundary conditions are *********/
/*****          to be used                                  *********/
/***** abc    : holds a b c alpha beta gamma Angstroms      *********/
/*****                                     and degrees      *********/
/***** latt_vec : holds cartessian vectors for a b c        *********/
/***** recip_latt_vec : holds recip. space a* b* c*         *********/
/********************************************************************/

EXTERNAL int pbc;
EXTERNAL double abc[6];
EXTERNAL double latt_vec[9];
EXTERNAL double real_latt_sizes[3];
EXTERNAL double old_latt_vec[9];
EXTERNAL double old_recip_latt_vec[9];
EXTERNAL double recip_latt_vec[9];
EXTERNAL double recip_latt_sizes[3];
EXTERNAL double cell_volume;


/*********************stuff to do with symmetry *********************/

EXTERNAL symm_ops symm[10]; /* the symmetry operations themselves */
EXTERNAL int symm_set;      /* flag to say symmetry set */
EXTERNAL int num_symm_ops;  /* number of symmetry operations defined (0=1!) */

/*********************Box limits for non periodic systems ***********/
EXTERNAL double box_limits[6];
EXTERNAL int user_box;		/* flag for user supplied or not */
EXTERNAL double box_fraction; /* fraction of extents to use */

/********************* Flags for energy and output *********************/
EXTERNAL int initial_minimize_template;
            /*flag: do a minimize template before we start or not */
EXTERNAL int initial_minimize_inpore;
            /*flag: do a minimize template inpore before we start or not */


EXTERNAL int verbose; /* verbose output flag */
EXTERNAL int steric;
EXTERNAL int non_bonded;
EXTERNAL int charges;
EXTERNAL double vdw_scale;
EXTERNAL double stop_ctf; /*stopping bit */
EXTERNAL double ch_ctf;
EXTERNAL double nb_ctf, nb_ctf_2; /* non-bond and charge cutoff distances */
EXTERNAL double ring_ctf, ring_ctf_2; /* ring formation cutoff distances */

EXTERNAL double max_shake_step;
EXTERNAL double  max_rock_step;

EXTERNAL int  num_modify_attempts;
EXTERNAL int  num_shake_attempts;
EXTERNAL int  num_rock_attempts;

EXTERNAL int animate_flag; /* flag for writing animation */
EXTERNAL int num_anime_frames; /* count of frames in animation */
EXTERNAL int num_goodun_frames; /* count of frames in animation */

EXTERNAL energy interaction_energy;
EXTERNAL internal_energy intra_energy;

/******************** DEBUG flag ******************************/

EXTERNAL int DEBUG;

/******************** Control data ****************************/

EXTERNAL double prob_test;

#undef EXTERNAL
