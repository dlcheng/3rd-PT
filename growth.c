#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/* for simulation 8411 */
static double Omega_m = 0.268;
static double Omega_v = 0.732;
static double Agr;

double g_factor(double z)
{ 
  
  return Omega_m * pow(1.0+z, 3) + (1.0 - Omega_m - Omega_v) * pow(1.0+z, 2) + Omega_v;	
}          /* end g_factor */

double unnorm_growth(double z)
{
  double g = g_factor(z);
  double O1 = Omega_m * pow(1.0+z, 3) / g;
  double O2 = Omega_v / g;
  
  double result;
  
  result = 1.0 / (1.0 + z) * 2.5 * O1;
  result /= pow(O1, 4.0/7.0) - O2 + (1.0 + 0.5 * O1) * (1.0 + O2/70.0);
  
  return result;
}         /* end unnorm_growth */

void init_growth()
{
  Agr = 1.0 / unnorm_growth(0);
}          /* end init_growth */

double growth_factor(double z)
{	
  return Agr * unnorm_growth(z);	
}        /* end growth_factor */  
