#define BOOLEAN int
#define TRUE 1
#define FALSE 0

#define NOT_SET "Not_Set"

#define MIN(A,B) (A < B ? A : B)
#define MAX(A,B) (A > B ? A : B)

#define PROB_TEST_DEFAULT 0.01   /* (real) default test probability */
#define MODIFY_TRY_DEFAULT 500   /* (int) default max no of modify attempts*/
#define STOP_CUTOFF_DEFAULT 30.0  /* (real) default stop cutoff distance */
#define RING_CUTOFF_DEFAULT 3.50  /* (real) default ring cutoff distance */
#define CH_CUTOFF_DEFAULT 30.0  /* (real) default charge cutoff distance */
#define NB_CUTOFF_DEFAULT 12.0  /* (real) default non-bond cutoff distance */
#define ATTEMPTS_DEFAULT 10     /* (int) default for shake/rock tries */
#define SHAKE_STEP_DEFAULT 0.1  /* (real) default for shake step */
#define ROCK_STEP_DEFAULT 1.0   /* (real) default for rock step */
#define VDW_SCALE_DEFAULT 1.0   /* (real) default vdw scaling parameter */

#define BOX_FRACTION_DEFAULT 0.8   /* (real) default box fraction parameter */

/****************************************************************/
/***** seed types for seed selection  ***************************/
#define MOLE 1
#define ARCH  2
#define FRAG 3
/***note abbreviated to avoid clash with reader***/


/*****default filenames for the intermediate discover files******/
#define DISCOVER_FORCEFIELD_DEFAULT "$BIOSYM_LIBRARY/cff91_czeo.bin"
#define TEMPLATE_STRATEGY_DEFAULT "temp_template_min.inp"
#define TEMPLATE_CAR_DEFAULT "temp_template_min.car"
#define TEMPLATE_MDF_DEFAULT "temp_template_min.mdf"
#define INPORE_STRATEGY_DEFAULT "temp_inpore_min.inp"
#define INPORE_CAR_DEFAULT "temp_inpore_min.car"
#define INPORE_MDF_DEFAULT "temp_inpore_min.mdf"
#define GOODUNS_DEFAULT    "good_uns"
#define BOX_OUTPUT_DEFAULT "box"

/*****filename for the peek file******/
#define PEEK_FILENAME "zebedde_peek.pcar"

/*****default filenames for analysis run i.e. none *******/
#define NO_ANALYSE "NO_ANALYSE"

/*****default command line for MOPAC*****/
#define DEFAULT_MOPAC_COMMANDLINE "PM3 NOINTER XYZ GEO-OK MMOK Precise "
#define DEFAULT_MOPAC_OUTPUT "mopac_minimize"

/******* defaults for animations ********************************/

#define DEFAULT_ANIMATION_FILE "animation"
#define DEFAULT_ANIMATION_TYPE BIOSYM_ANIMATION

/****************************************************************/
/***** key words for potential reader ***************************/
/****************************************************************/

#define INFO_LINE            "@"
#define TITLE_LINE           "!"
#define ILLUSTRATION_LINE    ">"
#define JUST_RETURN          "\n"

/****************************************************************/
/***** Non-bond potential keywords ******************************/
/****************************************************************/

#define POT_COMBINATION "combination"
#define POT_TYPE        "type"

#define R_EPS           "r-eps"
#define A_B             "A-B"
#define SIXTH_POWER     "sixth-power"
#define GEOMETRIC       "geometric"

/**********************************************************/
/****** Potential file headings for potential types *******/
/****** Can now do pcff or cvff Dave Willock Mar.97 *******/
/**********************************************************/

#define VERSION "#version" 
#define PCFF_STRING "pcff"
#define CVFF_STRING "cvff"
#define CFF91_STRING "cff91"

#define PCFF 0
#define CVFF 1
#define CFF91 2

#define NON_BOND_VDW_PCFF "#nonbond(9-6)"
#define EQUIVALENCE_PCFF "#equivalence"

#define NON_BOND_VDW_CVFF "#nonbond(12-6)"
#define EQUIVALENCE_CVFF "#equivalence"

#define NON_BOND_VDW_CFF91 "#nonbond(9-6)"
#define EQUIVALENCE_CFF91 "#equivalence"

#define QUARTIC_STRETCH    0 
#define MORSE_STRETCH      1
#define QUADRATIC_STRETCH  2

#define QUARTIC_STRETCH_STRING    "#quartic_bond"
#define MORSE_STRETCH_STRING      "#morse_bond"
#define QUADRATIC_STRETCH_STRING  "#quadratic_bond"

/**********************************************************/
/****** own error codes ***********************************/
/**********************************************************/

#define END_OF_INPUT "Thatsit" /* As defined in read_line */
#define END_OF_FILE -10 /* As defined in read_line */
