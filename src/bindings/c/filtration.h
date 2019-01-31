
#ifndef AETHER_FILTRATION_H
#define AETHER_FILTRATION_H

#include "symbol_export.h"

SHO_PUBLIC void calculate_flow(const char *mesh_type, int grav_dirn, double grav_factor, const char *bc_type, double inlet_bc, double outlet_bc, double L_p, double sigma, double pi_c, double pi_alv, double c_L);

#endif /* AETHER_FILTRATION_H */
