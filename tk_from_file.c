#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "allvars.h"
#include "proto.h"


static gsl_interp_accel *sp_tk_acc;
static gsl_spline *sp_tk;              /* spline of tk */
static double *kk;                     /* the k bins */
static double *tkk;                    /* the transfer function of all matter */

static double kf_max;                  /* the maximum k in the transfer function file */
static double kf_min;                  /* the minimum k */
static int Num_l = 472;                /* the line number of the file 
                                        * WMAP:   transf0.dat        ->  472
                                        * PLANCK: planck_transf.dat  ->  470
                                        */
static double omegab0 = 0.166;          /* the baryon fraction among all matter
                                        * WMAP:   0.166  
                                        * PLANCK: 0.155
                                        */

void init_tk_from_file()
{
   char filename[500] = "/Users/dalongcheng/Dropbox/Projects/Code/Codes/RSD-window/JV-1.3/TK/transf0.dat";
   FILE * fp;
   fp = fopen(filename, "r");
   int i;
   double temp_k, temp_tkc, temp_tkb, other_1, other_2;   
//   printf("Tk file open successfully!\n");

   kk  = (double *)malloc(Num_l * sizeof(double));
   tkk = (double *)malloc(Num_l * sizeof(double));
   sp_tk_acc = gsl_interp_accel_alloc();
   sp_tk     = gsl_spline_alloc(gsl_interp_cspline, Num_l);

   rewind(fp);

   for(i=0; i<Num_l; i++)
   	 {
      fscanf(fp, " %lf %lf %lf %lf %lf\n", &temp_k, &temp_tkc, &temp_tkb, &other_1, &other_2);
      kk[i] = temp_k;
      tkk[i] = temp_tkc * (1. - omegab0) + temp_tkb * omegab0;

      if(i == 0)
      	kf_min = temp_k;

      if(i == Num_l - 1)
      	kf_max = temp_k;
   	 } 

   gsl_spline_init(sp_tk, kk, tkk, Num_l);

   fclose(fp);
}  /* end init_tk_from_file */

double tk_from_file(double k)
{
  if(k < kf_min)
  	return 1.;

  if(k > kf_max)
  	return large_k_extrapolation(k);        /* or neglect all the modes larger than kf_max ~ 100 h/Mpc */

  return gsl_spline_eval(sp_tk, k, sp_tk_acc);
}  /* end tk_from_file */


/* k > kmax, using power law tk=A*k^a for estimation */
double large_k_extrapolation(double k)
{
  double a = log(tkk[Num_l-1]/tkk[Num_l-2]) / log(kk[Num_l-1]/kk[Num_l-2]);

  double result = tkk[Num_l-1] * pow(k/kf_max, a);

  return result;
}  /* end large_k_extrapolation */
