#define MAX_ATOMS  5000
#define MAXTYPES 10
#define MAX_LINE_LEN 256
#define BUFFER  256
#define PORE_GROUP pore
#define END_CAR "end"
#define NUM_ELEMENTS 104  /* update when more found we found deuterium! */

#define FILELEN_MAX 128

#define LINESIZ 256
#define MAXMOL 1000
#define MAXMOL3  3000
#define MAX_MSD_LIST 800

/****************************************************************/
/***** Maximum number of frames to read from a HISTORY file  ****/
/****************************************************************/

#define MAXFRAMES  2500
#define MAX_OUT_LIST 20000
#define MAXRDF    5000
#define MAX_CORR_VECS 8000
/****************************************************************/
/***** Maximum array dimension for users of ends in searches ****/
/****************************************************************/

#define MAX_ENDS 200

/****************************************************************/
/***** Maximum number of chosen pairs ***************************/
/****************************************************************/

#define MAX_CHOOSE 10

/****************************************************************/
/***** Maximum number of bins for histograms ********************/
/****************************************************************/

#define MAX_BINS 500


#define MAX_MONIT 20

/****************************************************************/
/***** Maximum number of atoms used to define a window **********/
/****************************************************************/

#define MAX_WINDOW_ATOMS 5      /* Max atoms to use in a windows definition  */
#define MAX_WINATOMS_PER_MOL 50 /* Max atoms in windows per molecule allowed */
#define MAX_WINDOWS_PER_MOL 10 /* Max atoms in windows per molecule allowed */
#define MAX_WIN_TOT 1000       /* Max windows allowed per frame *************/


/*** Ellipse and line drawing definitions *****************************/
/*** We draw three ellipses/lines so need three times the storage *****/
#define NUM_ELLIPSE_DOTS 50
#define NUM_ELLIPSE_DOTS_3 150
#define NUM_LINE_DOTS 50
#define NUM_LINE_DOTS_3 150
