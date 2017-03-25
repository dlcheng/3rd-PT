#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

/* 

*/                             
double p_delta_delta(double k, double z)
{
   double a = 1.0 / (1.0 + z);
   double p11 = a * a * p_linear(k);
   double p13 = 6.0 * pow(a, 4) * p_linear(k) * f3_pl(k);
   double p22 = 2.0 * pow(a, 4) * two_two_pl_pl(k, 3);

   double p = p11 + (p13 + p22)/pow(2*PI, 3);

   return p;
}        /* end p_delta_delta */

double p_theta_theta(double k, double z)
{
  double a = 1.0 / (1.0 + z);
  double p11 = a * p_linear(k);
  double p13 = 6.0 * pow(a, 3) * p_linear(k) * g3_pl(k);
  double p22 = 2.0 * pow(a, 3) * two_two_pl_pl(k, 2);

  double p = p11 + (p13 + p22)/pow(2*PI, 3);

  return p * a;
}        /* end p_theta_theta */             

double p_delta_theta(double k , double z)
{
  double a = 1.0 / (1.0 + z);
  double p11 =  -1.0 * pow(a, 1.5) * p_linear(k);
  double p13 =  -3.0 * pow(a, 3.5) * p_linear(k) * (g3_pl(k) + f3_pl(k));
  double p22 =  -2.0 * pow(a, 3.5) * two_two_pl_pl(k, 1);

  double p = p11 + (p13 + p22)/pow(2*PI, 3);

  return -1.0 * p  * sqrt(a);
}         /* end p_deta_theta */
