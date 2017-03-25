#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

static double relative_error = 1e-4;

/* The kernel integration of \int dq  F_3 P_L (q) , r is the integration parameter*/
double f3_pl_element(double r, void *param)
{	
   double k = *(double *) param;		
   double result;

  if(r > 200)
    result = -97.6;                   /* using the limit value to prevent round-off error */
  else
    {      
     result  =  12.0 / r / r;
     result += -158 + 100 * r *r - 42 * r * r * r * r;
     result +=  3.0 / pow(r, 3) * pow((r * r - 1), 3) * ( 2.0 + 7 * r * r) *log(fabs((1+r)/(1-r)));
    }

   result *= p_linear(k * r);

   return result;
}      /* end f3_pl_element */

double f3_pl(double k)
{  
   gsl_integration_workspace * gsl_w = gsl_integration_workspace_alloc(10000);
              	
   double result, error;
   double w = k/Kc;                              /* r cuts off at 1/w */
   		
   gsl_function F;
   F.function = &f3_pl_element;
   F.params = &k;

   gsl_integration_qag(&F, Epsilon, 1.0/w, 0, relative_error,  10000, GSL_INTEG_GAUSS15, gsl_w, &result, &error);   
   
   result *= PI / 756.0 * pow(k, 3);
   gsl_integration_workspace_free(gsl_w);   
   return result;
}     /* end sigma_m */	

/* The kernel integration of \int dq  G_3 P_L (q) , r is the integration parameter
  *
  *
  *
  */
double g3_pl_element(double r, void *param)
{	
   double k = *(double *) param;		
   double result;
   
   if(r > 200)
        result = -100.8;           /* using the limit value to prevent round-off error */
   else
        {
        result  =  12.0 / r / r;
        result += -82 + 4 * r *r - 6 * r * r * r * r;
        result +=  3.0 / pow(r, 3) * pow((r * r - 1), 3) * ( 2.0 + r * r) *log(fabs((1+r)/(1-r)));
        }

   result *= p_linear(k * r);

   return result;
}      /* end sigma_m_before_int_kernel*/

double g3_pl(double k)
{  
   gsl_integration_workspace * gsl_w = gsl_integration_workspace_alloc(10000);
              	
   double result, error;
   double w = k/Kc;                              /* r cuts off at 1/w */
   		
   gsl_function F;
   F.function = &g3_pl_element;
   F.params = &k;

  gsl_integration_qag(&F, Epsilon, 1.0/w, 0, relative_error, 10000, GSL_INTEG_GAUSS15, gsl_w, &result, &error);   
   
  result *= PI / 252.0 * pow(k, 3);
  gsl_integration_workspace_free(gsl_w);   
  return result;
}     /* end g3_pl */	

 /* x is the integration parameter, r and k are passed in as addtional information 
   *
   *
   *
   */   
 double two_two_pl_pl_element_element(double x, void * param)
 {
   double *pm = (double *) param;
   double r = pm[0];
   double k = pm[1];
   int flag = (int) pm[2];
   /* 
       flag ==1    G2F2
       flag ==2    G2G2
       flag ==3    F2F2
  */
   double result;
   double ep = sqrt(1 + r*r - 2*r*x);

  if(ep < Epsilon || k * ep >= Kc)
    return 0.0;                  /* mode is cut if it is smaller than the infrared limit or larger than the ultravoilet limit */

  if(flag == 1)	
    result  = (6 * r * x * x - 7 * x + r)* (10 * r * x  * x - 7 * x - 3 * r);

  if(flag == 2)
    result  = pow(6 * r * x * x - 7 * x + r, 2);

  if(flag == 3)
    result  = pow(10 * r * x  * x - 7 * x - 3 * r, 2);

   result /= pow(1 + r*r - 2*r*x , 2.0);
   result *= p_linear(k * ep);

   return result;
 }    /* end two_two_pl_pl_element */

/* here r is the parameter for the next integration  with the inetergation already done for x */
double two_two_pl_pl_element(double r, void * param)
{
   if(r < Epsilon)
   	return 0.0;

   double *pm   = (double *) param;
   double k = pm[0];
   int flag = (int) pm[1];

   double pm_pass[3];      /* the parameters pass to two_two_pl_pl_element_element */
   pm_pass[0] = r;
   pm_pass[1] = k;
   pm_pass[2] = flag;

   gsl_integration_workspace * gsl_w = gsl_integration_workspace_alloc(10000);             	
   double result, error;
   		
   gsl_function F;
   F.function = &two_two_pl_pl_element_element;
   F.params = pm_pass;

  /* GSL integration dealing with possible sigularities */
   gsl_integration_qag(&F, -1, 1, 0, relative_error, 10000, GSL_INTEG_GAUSS15, gsl_w, &result, &error);   
   
   result *= p_linear(k * r);
   gsl_integration_workspace_free(gsl_w);   

   return result;
}  /* end two_two_pl_pl_element */

/*       
    flag == 1   G2F2
    flag ==2    G2G2
    flag ==3    F2F2
*/
double two_two_pl_pl(double k, int flag)
{
  double w = k/Kc;
  double pm_pass[2];
  pm_pass[0] = k;
  pm_pass[1] = flag;

  gsl_integration_workspace * gsl_w = gsl_integration_workspace_alloc(10000);             	
  double result, error;
   		
  gsl_function F;
  F.function = &two_two_pl_pl_element;
  F.params = pm_pass;

  /* GSL integration dealing with possible sigularities */
  gsl_integration_qag(&F, 0, 1.0/w, 0, relative_error, 10000, GSL_INTEG_GAUSS15, gsl_w, &result, &error);   
   
  result *= (PI/98.0*pow(k, 3));
  gsl_integration_workspace_free(gsl_w);   

  return result;  
}   /* end two_two_pl_pl */

