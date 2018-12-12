module filtration_c
  implicit none
  private
  
contains
! 
!#########################################################################
! 
  subroutine calculate_flow_c(L_p, sigma, pi_c, pi_alv, c_L) bind(C, name="calculate_flow_c")
  
    use arrays,only: dp
    use filtration
    implicit none

    real(dp), intent(in) :: L_p
    real(dp), intent(in) :: sigma
    real(dp), intent(in) :: pi_c
    real(dp), intent(in) :: pi_alv
    real(dp), intent(in) :: c_L

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_calculate_flow(L_p, sigma, pi_c, pi_alv, c_L)
#else
    call calculate_flow(L_p, sigma, pi_c, pi_alv, c_L)
#endif

  end subroutine calculate_flow_c

! 
!#########################################################################
! 

end module filtration_c
