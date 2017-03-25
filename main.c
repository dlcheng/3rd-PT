#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  if(argc == 3)
    {
    int flag = atoi(argv[1]);
    double z = atof(argv[2]);
    jing_simulation_comparison(flag, z);
    }
  else
    {
     printf("Need two arguments.\n");
     printf("For example: 3rd_PT 8600 1.0\n");
     printf("where the first is the simulation flag, and the sencond is the redshift.\n");
    }    

  return 0;
}   /* end main */	

void Our_alphas()
{

  double kmin = 1;
  double kmax = 10;
  double k;
  int N_bin  = 10;
  double k_dis = log10(kmax/kmin) / (double) N_bin;
  int i;

  double avg_alpha_dd = 0.0;
  double avg_alpha_tt = 0.0;
  double avg_alpha_dt = 0.0;


  for(i=0; i<N_bin; i++)
    {
    k = log10(kmin) + (i+0.5) * k_dis;
    k = pow(10, k);

    double pik = pow(k, 3) / 2.0 / PI / PI;
    double pi_factor = 1.0 / pow(2.0 * PI, 3);

    double f3pl = f3_pl(k);
    double g3pl = g3_pl(k);
    double D11 = p_linear(k) * pik;

    double D13_dd = 6.0 * p_linear(k) * f3pl * pi_factor * pik;
    double D22_dd = 2.0 * two_two_pl_pl(k, 3) * pi_factor * pik;
    double D13_tt = 6.0 * p_linear(k) * g3pl * pi_factor * pik;
    double D22_tt = 2.0 * two_two_pl_pl(k, 2) * pi_factor  * pik; 
    double D13_dt = 3.0 * p_linear(k) * (g3pl + f3pl) * pi_factor * pik;
    double D22_dt = 2.0 * two_two_pl_pl(k, 1) * pi_factor * pik;   

    double alpha_dd_n = (D13_dd+D22_dd)/D11/D11;
    double alpha_tt_n = (D13_tt+D22_tt)/D11/D11;
    double alpha_dt_n = (D13_dt+D22_dt)/D11/D11;
    
    avg_alpha_dd += alpha_dd_n;
    avg_alpha_tt += alpha_tt_n;
    avg_alpha_dt += alpha_dt_n;
//    printf("%.6e\t%.6e\n", k, alpha_dd_n);
    }

    avg_alpha_dd /= N_bin;
    avg_alpha_tt /= N_bin;
    avg_alpha_dt /= N_bin;

    printf("%.6e\t%.6e\t%.6e\t%.6e\n", ns, avg_alpha_dd, avg_alpha_tt, avg_alpha_dt);
}  /* end Bernardeau_alpha_dd */    

void jing_simulation_comparison(int flag, double z)
{
  double kmin = 0.;
  double kmax = 0.40;;
  double k;
  int N_bin  = 20;
  double k_dis = (kmax-kmin) / (double) N_bin;
  int i;
  double a = 1.0 / (1.0 + z);
  double local_z = z;

  if (flag == 8601 || flag == 8600)
    set_parameters(1.0, 0.83, 2, 1e-5, 10000);     /* this is for simulation 8601, 8600 of SCDM power
                                                      flag == 2
                                                      Gamma = 0.5
                                                   */
  if (flag == 8610 || flag == 8411)                                                
    set_parameters(0.968, 0.83, 4, 1e-5, 10000);   /* for simulation 8610 and 8411, reading from WMAP TK */   

  if (flag == 8411)
    { 
     init_growth();
     a = growth_factor(z);                         /* find equivelance z */  
     local_z = 1.0 / a - 1.0;  
    }

  printf("%s\t%s\t%s\t%s\t%s\n","#(1)k", "(2)P_dd(k,z)", "(3)P_dt(k,z)", "(4)P_tt(k,z)", "(5)P_lin(k,z)");                                                

  for(i=0; i<N_bin; i++)
    {
     k = kmin + (i + 0.5) * k_dis;
     double p_dd, p_dt, p_tt;
     p_dd = p_delta_delta(k, local_z);
     p_tt = p_theta_theta(k, local_z);
     p_dt = p_delta_theta(k, local_z);
     printf("%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", k, p_dd, p_dt, p_tt, p_linear(k)*a*a);
    }

}  /* end jing_simulation_comparison */    
