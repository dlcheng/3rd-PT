#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

void set_parameters(double x1, double x2, int flag, double x3, double x4)
{
  ns      = x1;             
  sigma_8 = x2;
  
  tk_flag = flag;
  
  Epsilon = x3;
  Kc = x4;

  if(tk_flag == 4 || tk_flag == 5)
  	init_tk_from_file();

  if(tk_flag == 5)
  	init_power_from_file();
  
  init_power_amplitude();
}    /* end set_cosmological_parameters */