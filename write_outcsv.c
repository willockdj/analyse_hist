#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void write_csv(FILE *fp, char *p_title_x, char *p_title_y, 
               char *p_title_z,
               double *p_x, double *p_y, double *p_z,
               int have_z, int num);

void write_outcsv(double *p_eng_tot, 
                  double *p_temp_tot,  double *p_eng_cfg, 
                  double *p_eng_vdw,   double *p_eng_cou, 
                  double *p_eng_bnd,   double *p_eng_ang, 
                  double *p_eng_dih,   double *p_eng_tet, 
                  double *p_time,      double *p_eng_pv,   
                  double *p_temp_rot,  double *p_vir_cfg,  
                  double *p_vir_vdw,   double *p_vir_cou,   
                  double *p_vir_bnd,   double *p_vir_ang,   
                  double *p_vir_con,   double *p_vir_tet,   
                  double *p_volume,    double *p_temp_shl,
                  double *p_eng_shl,   double *p_vir_shl,
                  double *p_alpha,     double *p_beta,
                  double *p_gamma,     double *p_vir_pmf,
                  double *p_press,     double *p_enthalpy,
                  double *p_stat_time, int num_stat_list,
                  int num_in_list, 
                  int want_rdf,        rdf *p_rdf,
                  int num_rdfs,
                  analysis *p_anal_flags )
{
int is_start, i, done_time;
int last_start, place, num, occurances, iloop;
int first_time, done, found;

double *p_this_time;

char *p_key;
char dummy[100];

FILE *file_fp;

printf("In write_outcsv time list length : %d\n", num_in_list);

if (num_stat_list > 0 )
  {
    printf("Calling enthalpy write %10.6f %10.6f\n", *p_stat_time, *p_enthalpy);
    file_fp=fopen("enthalpy.csv","w");
    write_csv(file_fp, "time/ps", "enthalpy", " ",
              p_stat_time, p_enthalpy, p_enthalpy, FALSE, 
              num_stat_list);
    fclose(file_fp);
  }

if ( p_anal_flags->eng_tot )
  {
    file_fp=fopen("eng_tot.csv","w");
    write_csv(file_fp, "time/ps", "eng_tot", " ",
              p_time, p_eng_tot, p_eng_tot, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->temp_tot )
  {
    file_fp=fopen("temp_tot.csv","w");
    write_csv(file_fp, "time/ps", "temp_tot", " ",
              p_time, p_temp_tot, p_temp_tot, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_cfg )
  {
    file_fp=fopen("eng_cfg.csv","w");
    write_csv(file_fp, "time/ps", "eng_cfg", " ",
              p_time, p_eng_cfg, p_eng_cfg, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_vdw )
  {
    file_fp=fopen("eng_vdw.csv","w");
    write_csv(file_fp, "time/ps", "eng_vdw", " ",
              p_time, p_eng_vdw, p_eng_vdw, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_cou )
  {
    file_fp=fopen("eng_cou.csv","w");
    write_csv(file_fp, "time/ps", "eng_cou", " ",
              p_time, p_eng_cou, p_eng_cou, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_bnd )
  {
    file_fp=fopen("eng_bnd.csv","w");
    write_csv(file_fp, "time/ps", "eng_bnd", " ",
              p_time, p_eng_bnd, p_eng_bnd, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_dih )
  {
    file_fp=fopen("eng_dih.csv","w");
    write_csv(file_fp, "time/ps", "eng_dih", " ",
              p_time, p_eng_dih, p_eng_dih, FALSE, 
              num_in_list);
    fclose(file_fp);
  }   

if ( p_anal_flags->eng_tet )
  {
    file_fp=fopen("eng_tet.csv","w");
    write_csv(file_fp, "time/ps", "eng_tet", " ",
              p_time, p_eng_tet, p_eng_tet, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->eng_pv )
  {
    file_fp=fopen("eng_pv.csv","w");
    write_csv(file_fp, "time/ps", "eng_pv", " ",
              p_time, p_eng_pv, p_eng_pv, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->temp_rot )
  {
    file_fp=fopen("temp_rot.csv","w");
    write_csv(file_fp, "time/ps", "temp_rot", " ",
              p_time, p_temp_rot, p_temp_rot, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_cfg )
  {
    file_fp=fopen("vir_cfg.csv","w");
    write_csv(file_fp, "time/ps", "vir_cfg", " ",
              p_time, p_vir_cfg, p_vir_cfg, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_vdw )
  {
    file_fp=fopen("vir_vdw.csv","w");
    write_csv(file_fp, "time/ps", "vir_vdw", " ",
              p_time, p_vir_vdw, p_vir_vdw, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_cou )
  {
    file_fp=fopen("vir_cou.csv","w");
    write_csv(file_fp, "time/ps", "vir_cou", " ",
              p_time, p_vir_cou, p_vir_cou, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_bnd )
  {
    file_fp=fopen("vir_bnd.csv","w");
    write_csv(file_fp, "time/ps", "vir_bnd", " ",
              p_time, p_vir_bnd, p_vir_bnd, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_ang )
  {
    file_fp=fopen("vir_ang.csv","w");
    write_csv(file_fp, "time/ps", "vir_ang", " ",
              p_time, p_vir_ang, p_vir_ang, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_con )
  {
    file_fp=fopen("vir_con.csv","w");
    write_csv(file_fp, "time/ps", "vir_con", " ",
              p_time, p_vir_con, p_vir_con, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->vir_tet )
  {
    file_fp=fopen("vir_tet.csv","w");
    write_csv(file_fp, "time/ps", "vir_tet", " ",
              p_time, p_vir_tet, p_vir_tet, FALSE, 
              num_in_list);
    fclose(file_fp);
  }

if ( p_anal_flags->volume )
  {
    file_fp=fopen("volume.csv","w");
    write_csv(file_fp, "time/ps", "volume", " ",
              p_time, p_volume, p_volume, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->temp_shl )
  {
    file_fp=fopen("temp_shl.csv","w");
    write_csv(file_fp, "time/ps", "temp_shl", " ",
              p_time, p_temp_shl, p_temp_shl, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->eng_shl )
  {
    file_fp=fopen("eng_shl.csv","w");
    write_csv(file_fp, "time/ps", "eng_shl", " ",
              p_time, p_eng_shl, p_eng_shl, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->vir_shl )
  {
    file_fp=fopen("vir_shl.csv","w");
    write_csv(file_fp, "time/ps", "vir_shl", " ",
              p_time, p_vir_shl, p_vir_shl, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->alpha )
  {
    file_fp=fopen("alpha.csv","w");
    write_csv(file_fp, "time/ps", "alpha", " ",
              p_time, p_alpha, p_alpha, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->beta )
  {
    file_fp=fopen("beta.csv","w");
    write_csv(file_fp, "time/ps", "beta", " ",
              p_time, p_beta, p_beta, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->gamma )
  {
    file_fp=fopen("gamma.csv","w");
    write_csv(file_fp, "time/ps", "gamma", " ",
              p_time, p_gamma, p_gamma, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->vir_pmf )
  {
    file_fp=fopen("vir_pmf.csv","w");
    write_csv(file_fp, "time/ps", "vir_pmf", " ",
              p_time, p_vir_pmf, p_vir_pmf, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if ( p_anal_flags->press )
  {
    file_fp=fopen("press.csv","w");
    write_csv(file_fp, "time/ps", "press", " ",
              p_time, p_press, p_press, FALSE, 
              num_in_list);
    fclose(file_fp);
  } 

if (want_rdf)
  {
    printf("writing %d rdf lists\n",num_rdfs);
    for ( i=0; i <= num_rdfs; i++)
      {
        if (p_rdf->len > 0)
          {
            sprintf(dummy,"%s%s_rdf.csv", p_rdf->atom1, p_rdf->atom2);        
            printf("writing rdf file %s\n",dummy);
            file_fp=fopen(dummy,"w");

            write_csv(file_fp, "r", "g(r)", "n(r)",
                      &(p_rdf->r[0]), &(p_rdf->g[0]), &(p_rdf->n[0]), TRUE, 
                      p_rdf->len);
            fclose(file_fp);
          }
        p_rdf++;
      } 
  }

return;
}

