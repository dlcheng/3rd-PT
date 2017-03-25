#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

double prim_tk(double k)
{
 double tk = 1;
 
if(tk_flag == 1)
    tk = 1;       /* untouched for scale-free simulations */

if(tk_flag == 2)
	tk = prim_tk_bbks(k);

if(tk_flag == 3)
	tk = 1.0;     /* for E&H */

if(tk_flag == 4 || tk_flag == 5)
	tk = tk_from_file(k);

 return tk;
}    /* end prim_tk */

/* the Linear power spectrum at a =1 */
 double p_linear(double k)
 {
   double result, tk, Kt, Kl;
   
   if(tk_flag == 5)
   	 {
      Kt = 0.2;          /* switching to the power of N-body simualtion */
   	Kl = 2.0 * PI / 1200.0; 

      tk = prim_tk(k);
   
      if(k < Kc && k >= Kt)
        result = Aps * pow(k, ns) * tk * tk;
   
      if(k >= Kc)
   	  result = 0.;

      if(k <= Kl)
   	  result = 0.;

      if(k > Kl && k < Kt)
   	  result = lp_from_file(k);
   	  }
   	else
   	  {
   	   tk = prim_tk(k);
   	   result = Aps * pow(k, ns) * tk * tk;	
   	  } 
   
   return result;
 }  /* P_linear */
