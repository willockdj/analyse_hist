enum { PRIME_DIRECTIVE= 1, SECONDARY_DIRECTIVE, TERTIARY_DIRECTIVE, REMAINING_DIRECTIVE, ATOM_TYPE };


enum {TITLE = 1, 
      HISTORY_FILE, FIELD_FILE, RDF_FILE, OUT_FILE,  NUM_MOLS, 
      NUM_ATOMS, USE_TYPES,
      ANALYSE, MOI, LARGE, MIDDLE, SMALLEST, CONTACT, CONFIG_FILE, READ_PDB_FILE,
      COORDS, ANGLES, DIPOLES, FORSTER, ADD_ATOM, SUM_COMPS, CORRELATION, 
      TIME_STEP, ENG_TOT, TEMP_TOT, ENG_CFG, ENG_VDW, ENG_COU,
      ENG_BND, ENG_ANG, ENG_DIH, ENG_TET, ENG_PV, TEMP_ROT, VIR_CFG,
      VIR_VDW, VIR_COU, VIR_BND, VIR_ANG, VIR_CON, VIR_TET, VOLUME,
      TEMP_SHL, ENG_SHL, VIR_SHL, ALPHA, BETA, GAMMA, VIR_PMF,
      PRESS, PDB_FILE, MSD, STATIS_FILE, RAW_TRAJ, SKIP_OUT, BONDS, MONITOR, 
      LIMITS, BINS, XDATCAR_FILE, POSCAR, POTCAR, TRAJ, RGYR, WINDOWS, DRAW, ELLI_CALC, };

enum {RADIUS= 101, 
      TYPE, CHARGE, BISECTOR, ANY, SPEC, BOND2, ANGLE, TORSION, GYRATION, ATOM,
      AVERAGE, LINE, ELLIPSE, ELLI_START, ELLI_ATYPE, ELLI_VOL };

enum {STERIC= 201, ELECTROSTATIC };

enum {DEFAULT = 666};    
enum {UNIT = 998};    
enum {PARSE = 1001};

enum {NODIRECT = 999};
enum {BLANK_DIRECT = 999};

#define PRIME_DIRECTIVE_LIST \
	"titl", TITLE,"hist", HISTORY_FILE , "fiel", FIELD_FILE,  \
        "rdf_", RDF_FILE, "out_", OUT_FILE, \
        "num_", NUM_MOLS,"atom", NUM_ATOMS, "use_", USE_TYPES, \
        "anal", ANALYSE, "mome",MOI, "larg",LARGE,"midd",MIDDLE, \
        "smal",SMALLEST, "cont", CONTACT, "conf", CONFIG_FILE, "read", READ_PDB_FILE, \
        "coor", COORDS, "angl", ANGLES, "dipo", DIPOLES,\
        "fors", FORSTER, "add_", ADD_ATOM, "sum_", SUM_COMPS,\
        "corr", CORRELATION, "time_frame", TIME_STEP,\
        "eng_tot", ENG_TOT, "temp_tot", TEMP_TOT, "eng_cfg", ENG_CFG,\
        "eng_vdw",ENG_VDW, "eng_cou", ENG_COU, "eng_bnd", ENG_BND,\
        "eng_ang", ENG_ANG, "eng_dih", ENG_DIH, "eng_tet", ENG_TET,\
        "eng_pv", ENG_PV, "temp_rot", TEMP_ROT, "vir_cfg", VIR_CFG,\
        "vir_vdw", VIR_VDW, "vir_cou", VIR_COU, "vir_bnd", VIR_BND,\
        "vir_ang", VIR_ANG, "vir_con", VIR_CON, "vir_tet",VIR_TET,\
        "volume", VOLUME, "temp_shl",TEMP_SHL, "eng_shl", ENG_SHL,\
        "vir_shl", VIR_SHL, "alpha", ALPHA, "beta", BETA,\
        "gamma", GAMMA, "vir_pmf", VIR_PMF, "press", PRESS,\
        "pdb_",PDB_FILE,"msd",MSD,"statis",STATIS_FILE,\
        "raw_traj",RAW_TRAJ,"skip", SKIP_OUT, "bonds", BONDS, \
        "limit", LIMITS, "bins", BINS, "moni", MONITOR, \
        "xdat", XDATCAR_FILE, "poscar", POSCAR, "potcar", POTCAR, "traj", TRAJ, \
        "radi", RADIUS, "wind", WINDOWS, "draw", DRAW, "elli", ELLI_CALC, "",NODIRECT 
                        

#define SECOND_DIRECTIVE_LIST \
	"radi", RADIUS, "type", TYPE, "char", CHARGE, \
        "bise", BISECTOR,"any", ANY, "spec", SPEC, \
        "bond", BOND2, "angl", ANGLE, "dihe", TORSION, \
        "gyra", GYRATION, "atom", ATOM, "aver", AVERAGE, \
        "line", LINE, "elli", ELLIPSE, "atyp", ELLI_ATYPE, \
        "star", ELLI_START, "vol", ELLI_VOL, "", NODIRECT

#define THIRD_DIRECTIVE_LIST \
	 "ster", STERIC, "char", ELECTROSTATIC, "", NODIRECT

#define NULL_DIRECTIVE_LIST \
       "an", UNIT,             "(a", UNIT,           "of", UNIT, \
       "(c", UNIT,             "k ", UNIT,           "(k",  UNIT, \
       "on", UNIT,             "#", PARSE, "", NODIRECT

typedef struct 
{
  char *directive;
  int token_index;
}list;

static char  target[BUFFER];
static char  buf[BUFFER];
static char  *line, *last_tok;
extern int   read_new_line, line_no;

