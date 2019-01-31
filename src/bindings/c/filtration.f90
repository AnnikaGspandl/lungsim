module filtration_c
  implicit none
  private
  
contains
! 
!#########################################################################
! 
subroutine calculate_flow_c(mesh_type, mesh_type_len, grav_dirn, grav_factor, bc_type, bc_type_len, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L) bind(C, name="calculate_flow_c")

use arrays,only: dp
use iso_c_binding, only: c_ptr
use utils_c, only: strncpy
use filtration
use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
implicit none

integer,intent(in) :: mesh_type_len, bc_type_len
type(c_ptr), value, intent(in) :: mesh_type, bc_type
character(len=MAX_STRING_LEN) :: mesh_type_f, bc_type_f
    
integer, intent(in) :: grav_dirn
real(dp), intent(in) :: grav_factor, inlet_bc, outlet_bc

real(dp), intent(in) :: L_p
real(dp), intent(in) :: sigma
real(dp), intent(in) :: pi_c
real(dp), intent(in) :: pi_alv
real(dp), intent(in) :: c_L

call strncpy(mesh_type_f, mesh_type, mesh_type_len)
call strncpy(bc_type_f, bc_type, bc_type_len)

#if defined _WIN32 && defined __INTEL_COMPILER
call so_calculate_flow(mesh_type_f, grav_dirn, grav_factor, bc_type_f, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L)
#else
call calculate_flow(mesh_type_f, grav_dirn, grav_factor, bc_type_f, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L)
#endif

end subroutine calculate_flow_c

! 
!#########################################################################
! 

end module filtration_c
