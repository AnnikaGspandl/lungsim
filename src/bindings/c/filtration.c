
#include "filtration.h"

#include "string.h"

void calculate_flow_c(double *L_p, double *sigma, double *pi_c, double *pi_alv, double *c_L);

void calculate_flow(double L_p, double sigma, double pi_c, double pi_alv, double c_L)
{
	calculate_flow_c(&L_p, &sigma, &pi_c, &pi_alv, &c_L);
}
