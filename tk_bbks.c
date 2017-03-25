#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

double prim_tk_bbks(double k)
{
 double tk = 1.0;

 double b1, b2;
 double T_sig = 0.5;            /* 
                                    Shape parameter 
                                    0.25 for LCDM
                                    0.5  for SCDM
                                 */
 double q = k / T_sig;

 double a0 = 2.34;
 double a1 = 3.89;
 double a2 = 16.19;
 double a3 = 5.46;
 double a4 = 6.71;
    
 if(q <= 1e-8)
    tk = 1.0;                    /* for k = 0 return 1 */
 else 
   {
     b1 = 1.0 + a1 * q + a2 * a2 * q * q + a3 * a3 * a3 * q * q * q + a4 * a4 * a4 * a4 * q * q * q * q;
     b1 = 1.0 / pow (b1 , 0.25);
     b2 = log(1.0 + a0 * q)/ a0 / q;
     tk = b1 * b2;
   }
   
 return tk;		
}	                            /* end prim_tk_bbks */
