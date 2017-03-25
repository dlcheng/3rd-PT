#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

double win_top_hat(double x)
{
  double result = 1;
  
  if(x < 1e-8)
    result = 1.0 - 0.1 * x * x;
  else  
    result = 3.0 / pow(x, 3) * (sin(x) - x*cos(x));
      
    return result;    
}           /* end win_top_hat */

double sigma_int_element(double k, void *param)
{	
   double R = *(double *) param;		
   return 0.5 / PI / PI * pow(k, ns + 2.0) * prim_tk(k) * prim_tk(k) * win_top_hat(k * R) * win_top_hat(k * R);	
}      /* end sigma_init_element*/

double unnorm_sigma(double R)
{
   gsl_integration_workspace * w = gsl_integration_workspace_alloc(10000);             	
   double result, error;
   		
   gsl_function F;
   F.function = &sigma_int_element;
   F.params = &R;

   gsl_integration_qagiu(&F, 0, 0, 1e-6, 10000, w, &result, &error);   
   
   gsl_integration_workspace_free(w);   
   return result;
}     /* end unnorm_sigma */	

void init_power_amplitude()
{
  
if(tk_flag == 1)
  { 
   double zi = 16.38;          /* ns = -2, simulation 2730 */
   double N = 2048;
   double L = 1200;
   double C = 0.4241;          /* ns = -2, simulation 2730 */
     
   Aps = pow(1.0 + zi, 2);
   Aps *= (C * C / N / N / N);
   Aps *= pow(L / N / PI,  ns);
   Aps *= pow(L, 3);
  }
 else
  Aps = sigma_8 * sigma_8 / unnorm_sigma(8.0);

}  /* end init_sigma_m */	

