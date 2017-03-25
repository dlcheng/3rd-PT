#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"

/* Cosmological parameters */
 double ns;         
 double sigma_8;

/* Power amplitude */
 double Aps;
 int tk_flag;                      /* the types of the transfer function
                                          * tk_flag = 1 --> scale free simulations
                                          * tk_flag = 2 --> BBKS
                                          * tk_flag = 3 --> E&H
                                          * tk_flag = 4 --> Read from file
                                          */

/* Define the infrared and ultravoilet cutoff */
 double  Epsilon;                  /* for given k, the modes smaller than EPSILON * k are considered with no contribution */
 double  Kc;                       /* in unit of h/Mpc, modes larger than KC are considered with no contribution */ 

