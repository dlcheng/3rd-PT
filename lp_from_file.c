#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"

/* Return the linear power of simulation 8610 at z=3.48 between Kl and Kt */

static gsl_interp_accel *sp_lp_acc;
static gsl_spline *sp_lp;              /* spline of tk */
static double *k_bin;                  /* the k bins */
static double *p_bin;                  /* the power bins */

static double k_max;                   /* the local k max */
static double k_min;                   /* the local k min */
static int Num_l = 20;                 /* the line number of the file */
static double z_l = 3.48;              /* the redshift of the file */

void init_power_from_file()
{
   char filename[500] = "/Users/dalongcheng/Dropbox/Documents/Python-workspace/RSD-window/Data/Bin-20/wf8610_0255_z_3.48_2048.txt";
   FILE * fp;
   fp = fopen(filename, "r");
   int i;
   double temp_k, temp_wk, temp_dk;
   int other_1;   
   double a = 1.0 / (1.0 + z_l);
   printf("Lp file open successfully!\n");

   k_bin  = (double *)malloc(Num_l * sizeof(double));
   p_bin  = (double *)malloc(Num_l * sizeof(double));
   sp_lp_acc = gsl_interp_accel_alloc();
   sp_lp     = gsl_spline_alloc(gsl_interp_cspline, Num_l);

   rewind(fp);

   for(i=0; i<Num_l; i++)
   	 {
      fscanf(fp, "%lf %lf %lf %d\n", &temp_k, &temp_wk, &temp_dk, &other_1);
      k_bin[i] = temp_k;
      p_bin[i] = temp_dk * 2. * PI * PI / pow(temp_k, 3.) / a / a;  /* the power at z=0 */

      if(i == 0)
      	k_min = temp_k;

      if(i == Num_l - 1)
      	k_max = temp_k;
   	 } 

   gsl_spline_init(sp_lp, k_bin, p_bin, Num_l);

   fclose(fp);
}  /* end init_tk_from_file */

double lp_from_file(double k)
{
  if(k < k_min)
  	return low_k_extrapolation(k);

  return gsl_spline_eval(sp_lp, k, sp_lp_acc);
}  /* end tk_from_file */

/* k < kmin, using power law lp=A*k^a for estimation */
double low_k_extrapolation(double k)
{
  double a = log(p_bin[1]/p_bin[0]) / log(k_bin[1]/k_bin[0]);

  double result = p_bin[0] * pow(k/k_min, a);

  return result;
}  /* end low_k_extrapolation */
