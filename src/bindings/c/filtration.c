
#include "filtration.h"

#include "string.h"

void calculate_flow_c(const char *mesh_type, int *mesh_type_len, int *grav_dirn, double *grav_factor, const char *bc_type, int *bc_type_len, double *inlet_bc, double *outlet_bc, double *L_p, double *sigma, double *pi_c, double *pi_alv, double *c_L);

void calculate_flow(const char *mesh_type, int grav_dirn, double grav_factor, const char *bc_type, double inlet_bc, double outlet_bc, double L_p, double sigma, double pi_c, double pi_alv, double c_L)
{
	int mesh_type_len = strlen(mesh_type);
	int bc_type_len = strlen(bc_type);
	
	calculate_flow_c(mesh_type, &mesh_type_len, &grav_dirn, &grav_factor, bc_type, &bc_type_len, &inlet_bc, &outlet_bc, &L_p, &sigma, &pi_c, &pi_alv, &c_L);
}
