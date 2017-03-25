#ifndef ALLVAR_H
 #include "allvars.h"
#endif

double f3_pl_element(double r, void *param);
double f3_pl(double k);
double g3_pl_element(double r, void *param);
double g3_pl(double k);
double two_two_pl_pl_element_element(double x, void * param);
double two_two_pl_pl_element(double r, void * param);
double two_two_pl_pl(double k, int flag);

double p_delta_delta(double k, double z);
double p_theta_theta(double k, double z);
double p_delta_theta(double k, double z);

void set_parameters(double x1, double x2, int flag, double x3, double x4);

double prim_tk_bbks(double k);

double prim_tk(double k);
double p_linear(double k);

double win_top_hat(double x);
double sigma_int_element(double k, void *param);
double unnorm_sigma(double R);
void init_power_amplitude();

void Our_alphas();
void jing_simulation_comparison(int flag, double z);

void init_tk_from_file();
double tk_from_file(double k);
double large_k_extrapolation(double k);

void init_power_from_file();
double lp_from_file(double k);
double low_k_extrapolation(double k);

double g_factor(double z);
double unnorm_growth(double z);
void init_growth();
double growth_factor(double z);