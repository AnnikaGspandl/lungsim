module filtration
  !*Brief Description:* This module handles all code specific to simulating filtration in the lungs
  !
  !*LICENSE:*
  !
  !
  !
  !*Full Description:*
  !
  ! This module handles all code specific to simulating filtration in the lungs

  use arrays,only: dp
  implicit none
  
  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private
  public calculate_flow
  !public evaluate_vent_edema

    contains

  !!!###################################################################################
  !*calculate_flow:* --- Main Subroutine, gets called by Python Script ---
  ! Calculates transcapillary flow (filtration) J_v according to Starling equation: 
  ! J_v = L_p * S * ((P_c - P_alv) - sigma * (pi_c - pi_alv))
  ! 
  ! Assumptions:    - L_p, sigma, pi_c, pi_alv = constant over time & over whole lung
  !                 - Fluid flow in alveoli (no interstitial space)
  !                 - Thickness of fluid layer in alveoli increases -> ball sheet  
  subroutine calculate_flow(mesh_type, grav_dirn, grav_factor, bc_type, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L)  
    !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_CALCULATE_FLOW" :: CALCULATE_FLOW
  
    use diagnostics, only: enter_exit,get_diagnostics_on
    use arrays, only: num_elems,num_elems_aw,elem_field
    use pressure_resistance_flow
    use indices!, only: num_ne
    implicit none
    
    ! In-/Output
    real(dp),intent(in) :: L_p     ! Hydraulic conductivity in mm/(s*Pa)
    real(dp),intent(in) :: sigma   ! Reflection coefficient
    real(dp),intent(in) :: pi_c    ! Capillary osmotic pressure in mmH2O
    real(dp),intent(in) :: pi_alv  ! Alveolar fluid osmotic pressure in mmH2O
    real(dp),intent(in) :: c_L     ! Lymphatic clearance/ flow in ml/(s*g)
    
    ! Inputs for Perfusion Model
    real(dp) :: grav_factor,inlet_bc,outlet_bc
    integer :: grav_dirn
    character(len=60) :: mesh_type,bc_type
    
    ! local variables
    integer :: ne
    character(len=60) :: sub_name
    logical :: arg
    logical :: diags
    
    ! ###########################################################################
    
    !sub_name = 'calculate_flow'
    call enter_exit(sub_name,1)
    call get_diagnostics_on(diags)
    num_elems_aw = 0
    do ne=1,num_elems
      if(elem_field(ne_group,ne).eq.0.0_dp)then !Artery or airway
        num_elems_aw = num_elems_aw + 1
      endif
    enddo


    if(diags)write(*,*) "Call SR evaluate_vent_edema"
    call evaluate_vent_edema(mesh_type, grav_dirn, grav_factor, bc_type, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L)
       
    call enter_exit(sub_name,2)
    
  end subroutine calculate_flow
  
  
!!!###################################################################################
!!!################################ventilation########################################
!*evaluate_vent_edema:* Sets up and solves venilation model
  subroutine evaluate_vent_edema(mesh_type, grav_dirn, grav_factor, bc_type, inlet_bc, outlet_bc, L_p, sigma, pi_c, pi_alv, c_L)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_EVALUATE_VENT_EDEMA" :: EVALUATE_VENT_EDEMA
    use arrays,only: dp,elem_field,elem_units_below,node_field,num_elems,num_elems_aw,num_units,units_vent,unit_field
    use indices!,only: ne_Vdot,ne_Vdot0,ne_t_resist,nj_aw_press,nu_air_press,nu_comp,nu_dpdt,nu_pe,nu_rad,nu_SA,nu_vt
    use exports,only: export_1d_elem_field,export_terminal_solution, export_starling_variables
    use other_consts !PI
    use diagnostics, only: enter_exit,get_diagnostics_on
    use geometry,only: set_initial_volume,volume_of_mesh
    use pressure_resistance_flow
    implicit none

    ! Inputs for Perfusion Model
    real(dp) :: grav_factor,inlet_bc,outlet_bc
    integer :: grav_dirn
    integer :: perf_call_number
    character(len=60) :: mesh_type,bc_type
    
    ! In-/Output for Filtration
    real(dp),intent(in) :: L_p     ! Hydraulic conductivity in mm/(s*Pa)
    real(dp),intent(in) :: sigma   ! Reflection coefficient
    real(dp),intent(in) :: pi_c    ! Capillary osmotic pressure in mmH2O
    real(dp),intent(in) :: pi_alv  ! Alveolar fluid osmotic pressure in mmH2O
    real(dp),intent(in) :: c_L     ! Lymphatic clearance/ flow in ml/(s*g)
    
    ! Local variables
    integer :: Gdirn,iter_step,n,ne,num_brths,num_itns,nunit,count_aw, timestep_int
    real(dp) :: pleural,ChestWallRestVol,chest_wall_compliance,constrict,COV,&
         dpmus,dt,endtime,err_est,err_tol,FRC,i_to_e_ratio,init_vol,&
         last_vol,now_vol,Pcw,p_mus,pmus_factor_in,pmus_factor_ex,pmus_step,&
         ppl_current,pptrans,press_in,prev_flow,ptrans_frc,refvol,RMaxMean,&
         RMinMean,sum_dpmus,sum_dpmus_ei,sum_expid,sum_tidal,Texpn,T_interval,&
         time,Tinsp,totalc,Tpass,ttime,undef,volume_target,volume_tree,WOBe,&
         WOBr,WOBe_insp,WOBr_insp,WOB_insp
    real(dp) :: ppl_edema(num_units)
    real(dp) :: Jv_sum              ! Summarized Filtration in whole lung per time step
    real(dp), allocatable :: check_vol(:)
    character :: expiration_type*(10)
    logical :: CONTINUE,converged
    logical :: diags
    character(len=36) :: filenametime
    character(len=26) :: filename
    character(len=10) :: timestep_char

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'evaluate_vent_edema'
    call enter_exit(sub_name,1)
    call get_diagnostics_on(diags)

    init_vol = 1.0_dp
    
!!! -------------  DESCRIPTION OF IMPORTANT VARIABLES ---------------
!!! pmus_factor (_in and _ex) are used to scale the driving pressureS, to converge
!!! the tidal volume and expired volume to the target volume. Note that this does
!!! not represent passive recoil of the lung. It assumes a sinusoidal pressure
!!! acting external to the lung, approximately the opposite of the inspiration
!!! muscle pressure.

!!! T_interval !the total length of the breath

!!! Gdirn !direction of gravitational gradient for tissue volumes
!!! 'Gdirn' is 1(x), 2(y), 3(z); upright lung (for our models) is z, supine is y.

!!! press_in !constant value of pressure at the entry to the model

!!! COV !coefficient of variation of tissue compliance

!!! RMaxMean !ratio max to mean volume

!!! RMinMean !ratio min to mean volume

!!! i_to_e_ratio ! ratio of inspiration time to expiration time

!!! refvol !sets the 'zero stress' state at half of our lung volume (OK for FRC)

!!! volume_target ! the target tidal volume, in mm^3

!!! pmus_step !the change in Ppl that is driving the flow, in Pa

!!! Texpn ! time for expiration

!!! Tinsp ! time for inspiration

!!! undef !!!! undef is the volume at which our stress-strain relationship has zero
!!! stress. i.e. the tissue is 'undeformed'. Lung tissue is never stress-free,
!!! however our continuum model approach requires this to be defined (a theoretical
!!! state: should be less than our minimum simulation volume, e.g. RV).

!!! sum_tidal !tidal volume for a breath

!!! sum_expid ! sum of expired volume for a breath. Should equal the sum inspired.
!!! -----------------------------------------------------------------------------

!!! set default values for the parameters that control the breathing simulation
!!! these should be controlled by user input (showing hard-coded for now)

    perf_call_number = 0

    ! Initialize summarizd filtration per unit for whole breath
        unit_field(nu_filt_cleared_sum,nunit) = 0.0_dp
    
    call read_params_evaluate_flow_edema(Gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)
    call read_params_main_edema(num_brths, num_itns, dt, err_tol)

    allocate(check_vol(num_brths))
    check_vol = 0.0_dp
    
!!! calculate key variables from the boundary conditions/problem parameters
    Texpn = T_interval / (1.0_dp+i_to_e_ratio)
    Tinsp = T_interval - Texpn

!!! store initial branch lengths, radii, resistance etc. in array 'elem_field'
    call update_elem_field_edema

    call volume_of_mesh(init_vol,volume_tree,1) ! to get deadspace volume

!!! distribute the initial tissue unit volumes along the gravitational axis.
    call set_initial_volume(Gdirn,COV,FRC*1.0e+6_dp,RMaxMean,RMinMean,1)
    undef = refvol * (FRC*1.0e+6_dp-volume_tree)/DBLE(elem_units_below(1))

!!! calculate the total model volume
    call volume_of_mesh(init_vol,volume_tree,1)

    write(*,'('' Anatomical deadspace = '',F8.3,'' ml'')') volume_tree/1.0e+3_dp ! in mL
    write(*,'('' Respiratory volume   = '',F8.3,'' L'')') (init_vol-volume_tree)/1.0e+6_dp !in L
    write(*,'('' Total lung volume    = '',F8.3,'' L'')') init_vol/1.0e+6_dp !in L
    unit_field(nu_dpdt,1:num_units) = 0.0_dp

!!! calculate the compliance of each tissue unit
    call tissue_compliance_edema(chest_wall_compliance,undef)
    totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance

    call update_mean_pleural_pressure_edema(ppl_current) !calculate new pleural pressure
    pptrans=SUM(unit_field(nu_pe,1:num_units))/num_units

    
    
    ChestWallRestVol = init_vol + 0.2e+6_dp/98.0665_dp * (-ppl_current)
    Pcw = (ChestWallRestVol - init_vol)/(0.2e+6_dp/98.0665_dp)

!!! write out the header information for run-time output
    write(*,'(2X,''Time'',3X,''Inflow'',4X,''V_t'',5X,''Raw'',5X,&
         &''Comp'',4X,''Ppl'',5X,''Ptp'',5X,''VolL'',4X,''Pmus'',&
         &4X,''Pcw'',2X,''Pmus-Pcw'')')
    write(*,'(3X,''(s)'',4X,''(mL/s)'',3X,''(mL)'',1X,''(cmH/L.s)'',&
         &1X,''(L/cmH)'',1X,''(...cmH2O...)'',&
         &4X,''(L)'',5X,''(......cmH2O.......)'')')

    write(*,'(F7.3,2(F8.1),8(F8.2))') &
         0.0,0.0,0.0, &  !time, flow, tidal
         elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !airway resistance (cmH2O/L.s)
         totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
         ppl_current/98.0665_dp, & !Ppl (cmH2O)
         -ppl_current/98.0665_dp, & !mean Ptp (cmH2O)
         init_vol/1.0e+6_dp, & !total model volume (L)
         0.0, & !Pmuscle (cmH2O)
         Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
         (-Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)

!!! Initialise variables:
    pmus_factor_in = 1.0_dp
    pmus_factor_ex = 1.0_dp
    time = 0.0_dp !initialise the simulation time.
    n = 0 !initialise the 'breath number'. This is incremented at the start of each breath.
    sum_tidal = 0.0_dp ! initialise the inspired and expired volumes
    sum_expid = 0.0_dp
    last_vol = init_vol

    CONTINUE=.TRUE.
    do while(CONTINUE)
       n=n+1 !increment the breath number
       ttime=0.0_dp !each breath starts with ttime=0
       endtime = T_interval * n - 0.5_dp * dt !the end time of the breath
       p_mus = 0.0_dp !initialise the muscle pressure to zero
       ptrans_frc=SUM(unit_field(nu_pe,1:num_units))/num_units !ptrans at frc
       if(n.gt.1)then !write out 'end of breath' information
          write(*,'('' End of breath, inspired = '',F10.2,'' L'')') sum_tidal/1.0e+6_dp
          write(*,'('' End of breath, expired  = '',F10.2,'' L'')') sum_expid/1.0e+6_dp
          write(*,'('' Peak muscle pressure    = '',F10.2,'' cmH2O'')') &
               pmus_step*pmus_factor_in/98.0665_dp
          write(*,'('' Drift in FRC from start = '',F10.2,'' %'')') &
               100*(now_vol-init_vol)/init_vol
          write(*,'('' Difference from target Vt = '',F8.2,'' %'')') &
               100*(volume_target-sum_tidal)/volume_target
          write(*,'('' Total Work of Breathing ='',F7.3,''J/min'')')WOB_insp
          write(*,'('' elastic WOB ='',F7.3,''J/min'')')WOBe_insp
          write(*,'('' resistive WOB='',F7.3,''J/min'')')WOBr_insp

          if(DABS(volume_target).gt.1.0e-5_dp)THEN
             ! modify driving muscle pressure by volume_target/sum_tidal
             ! this increases p_mus for volume_target>sum_tidal, and
             ! decreases p_mus for volume_target<sum_tidal
             pmus_factor_in=pmus_factor_in*DABS(volume_target/sum_tidal)
             pmus_factor_ex=pmus_factor_ex*DABS(volume_target/sum_expid)
          endif
          sum_tidal=0.0_dp !reset the tidal volume
          sum_expid=0.0_dp !reset the expired volume
          unit_field(nu_vt,1:num_units) = 0.0_dp !reset acinar tidal volume
          sum_dpmus=0.0_dp
          sum_dpmus_ei=0.0_dp
       endif

       timestep_int = 0
       
!!! Do the simulation of each breath
       do while (time.LT.endtime)
!          time = time + dt
          ttime = ttime + dt
          timestep_int = timestep_int + 1
          ! set the increment in driving (muscle) pressure
          if(expiration_type(1:6).eq.'active')then
             if(ttime.lt.Tinsp)then
                dpmus=pmus_step*pmus_factor_in*PI* &
                     sin(2.0_dp*pi/(2.0_dp*Tinsp)*ttime)/(2.0_dp*Tinsp)*dt
             elseif(ttime.LE.Tinsp+Texpn)then
                dpmus=pmus_step*pmus_factor_ex*PI* &
                     sin(2.0_dp*pi*(0.5d0+(ttime-Tinsp)/(2.0_dp*Texpn)))/ &
                     (2.0_dp*Texpn)*dt
             endif
          elseif(expiration_type(1:7).eq.'passive')then
             if(ttime.le.Tinsp+0.5d0*dt)then
                dpmus=pmus_step*pmus_factor_in*PI*dt* &
                     sin(pi*ttime/Tinsp)/(2.0_dp*Tinsp)
                sum_dpmus=sum_dpmus+dpmus
                sum_dpmus_ei=sum_dpmus
             else
                Tpass=0.1d0
                dpmus=MIN(-sum_dpmus_ei/(Tpass*Texpn)*dt,-sum_dpmus)
                sum_dpmus=sum_dpmus+dpmus
             endif
          endif

          p_mus = p_mus + dpmus !current value for muscle pressure
          prev_flow=elem_field(ne_Vdot,1)

!!! Solve for a new flow and pressure field
!!! We will estimate the flow into each terminal lumped
!!! parameter unit (assumed to be an acinus), so we can calculate flow
!!! throughout the rest of the tree simply by summation. After summing
!!! the flows we can use the resistance equation (P0-P1=R1*Q1) to update
!!! the pressures throughout the tree.

          !initialise Qinit to the previous flow
          elem_field(ne_Vdot0,1:num_elems_aw) = elem_field(ne_Vdot,1:num_elems_aw)
          converged = .FALSE.
          iter_step=0
          do while (.not.converged)
             iter_step=iter_step+1 !count the iterative steps
             call estimate_flow_edema(dpmus,dt,err_est) !analytic solution for Q
             if(iter_step.gt.1.and.err_est.lt.err_tol)then
                converged=.TRUE.
             else if(iter_step.gt.num_itns)then
                converged=.TRUE.
                write(*,'('' Warning: lower convergence '// &
                     'tolerance and time step - check values, Error='',D10.3)') &
                     err_est
             endif
             call sum_elem_field_from_periphery_edema(ne_Vdot) !sum the flows recursively UP the tree
             call update_elem_field_edema !updates resistances
             call update_node_pressures_edema(press_in) !updates the pressures at nodes
             call update_unit_dpdt_edema(dt) ! update dP/dt at the terminal units
          enddo !converged

          call update_unit_volume_edema(dt,Tinsp,Texpn) ! Update tissue unit volumes, and unit tidal volumes
          call volume_of_mesh(now_vol,volume_tree,1) !calculate the mesh volume, store in 'now_vol'
          call update_elem_field_edema  !update element lengths, volumes, resistances
          call tissue_compliance_edema(chest_wall_compliance,undef) !update the unit compliances, uses 'undef' as input
          totalc = SUM(unit_field(nu_comp,1:num_units)) !the total model compliance
          call update_mean_pleural_pressure_edema(ppl_current) !calculate new pleural pressure
          call update_proximal_pressure_edema !updates values of pressure at proximal nodes of end branches
          call calculate_work_edema(now_vol-init_vol,now_vol-last_vol,WOBe,WOBr,pptrans)!calculate work of breathing
          last_vol=now_vol
          Pcw = (now_vol - ChestWallRestVol)/(0.2d6/98.0665_dp)
                 
          
          ! Calculate Filatration in last breath (after convergence achieved)

          ! Check: Volume converged in last breath? If yes -> calculate filtration
          if(n.GT.1)then
            if(check_vol(n-1).LT.0.1_dp)then
                
           !if(diags)
           write(*,*) "Filtration calculated in this breath :"
                    
          ! Calculate Ppl = -Pel + Palv
            call update_pleural_pressure_edema
    
          ! Get filtration surface area SA in mm^2
            call update_surface_area_edema
            perf_call_number = perf_call_number + 1
            
          ! Get capillary hydrostatic pressure P_c in mmH2O (-> perfusion model)
            if(diags)write(*,*) "Call SR evaluate_prq from perfusion model", time,n,perf_call_number
            call evaluate_prq(mesh_type, grav_dirn, grav_factor, bc_type, inlet_bc, outlet_bc)
        
          ! Calculate transcapillary fluid flow J_v per time step according to Starling Equation in m^3
            
            call update_filtration_edema(L_p, sigma, pi_c, pi_alv, c_L, dt, Jv_sum)
            
          ! Calculate average values for terminal export file
            call calculate_terminal_exports(timestep_int)
            
          ! Define name of export file in this time step
            filename = '_starling_variables.exnode'
            write (timestep_char, "(A8,I2)")  "Timestep", timestep_int
            timestep_char = trim(timestep_char)
            filenametime =  timestep_char//filename
            
          ! Export all the values that go into the Starling equation
            call export_starling_variables(filenametime,'filt_model')!(time, ttime,Jv_sum)
          
            endif   
          endif
         
          
          time = time + dt
          write(*,'(F7.3,2(F8.1),8(F8.2))') &
               time, & !time through breath (s)
               elem_field(ne_Vdot,1)/1.0e+3_dp, & !flow at the inlet (mL/s)
               (now_vol - init_vol)/1.0e+3_dp, & !current tidal volume (mL)
               elem_field(ne_t_resist,1)*1.0e+6_dp/98.0665_dp, & !airway resistance (cmH2O/L.s)
               totalC*98.0665_dp/1.0e+6_dp, & !total model compliance
               ppl_current/98.0665_dp, & !Ppl (cmH2O)
               pptrans/98.0665_dp, & !mean Ptp (cmH2O)
               now_vol/1.0e+6_dp, & !total model volume (L)
               p_mus/98.0665_dp, & !Pmuscle (cmH2O)
               -Pcw/98.0665_dp, & !Pchest_wall (cmH2O)
               (p_mus+Pcw)/98.0665_dp !Pmuscle - Pchest_wall (cmH2O)

          ! increment the tidal volume, or the volume expired
          if(elem_field(ne_Vdot,1).gt.0.0_dp)then
             sum_tidal=sum_tidal+elem_field(ne_Vdot,1)*dt
          else
             sum_expid=sum_expid-elem_field(ne_Vdot,1)*dt
             if(prev_flow.gt.0.0_dp)then
                WOBe_insp=(WOBe+sum_tidal*ptrans_frc*1.0e-9_dp)*(30.0_dp/Tinsp)
                WOBr_insp=WOBr*(30.0_dp/Tinsp)
                WOB_insp=WOBe_insp+WOBr_insp
                WOBe=0.0_dp
                WOBr=0.0_dp
             endif
          endif

       ENDDO !while ttime<endtime

       !...  CHECK WHETHER TO CONTINUE     
       check_vol(n) = DABS(100.0_dp*(volume_target-sum_tidal)/volume_target)
       !write(*,*) "check_vol for this breath calculated"
       if(n.ge.num_brths)then
          CONTINUE=.FALSE.
       elseif(DABS(volume_target).GT.1.0e-3_dp)THEN
          if(check_vol(n).GT.0.1_dp.OR.(n.LT.2))then
             CONTINUE=.TRUE.
          ! after convergence achieved -> do one last loop for calculation of filtration
          elseif((check_vol(n).LE.0.1_dp).AND.(check_vol(n-1).GT.0.1_dp))then
              CONTINUE=.TRUE.
          else
              CONTINUE=.FALSE.
          endif
       endif

    enddo !...WHILE(CONTINUE)
    deallocate(check_vol)
    
    write(*,'(''End of breath, inspired = '',F10.2,'' L'')') sum_tidal/1.0e+6_dp
    write(*,'(''End of breath, expired  = '',F10.2,'' L'')') sum_expid/1.0e+6_dp
    write(*,'(''Peak muscle pressure    = '',F10.2,'' cmH2O'')') &
         pmus_step*pmus_factor_in/98.0665_dp
    write(*,'(''Drift in FRC from start = '',F10.2,'' %'')') &
         100.0_dp*(now_vol-init_vol)/init_vol
    write(*,'(''Difference from target Vt = '',F8.2,'' %'')') &
         100*(volume_target-sum_tidal)/volume_target
    write(*,'(''Total Work of Breathing ='',F7.3,''J/min'')')WOB_insp
    write(*,'(''elastic WOB ='',F7.3,''J/min'')')WOBe_insp
    write(*,'(''resistive WOB='',F7.3,''J/min'')')WOBr_insp

!!! Transfer the tidal volume for each elastic unit to the terminal branches, and sum up the tree.
!!! Divide by inlet flow. This gives the time-averaged and normalised flow field for the tree.
    do nunit=1,num_units !for each terminal element only (with tissue units attached)
       ne=units_vent(nunit) !local element number
       elem_field(ne_Vdot,ne) = unit_field(nu_vt,nunit)
    enddo
    call sum_elem_field_from_periphery_edema(ne_Vdot)
    elem_field(ne_Vdot,1:num_elems_aw) = elem_field(ne_Vdot,1:num_elems_aw)/elem_field(ne_Vdot,1)

!    call export_terminal_solution(TERMINAL_EXNODEFILE,'terminals')

    call enter_exit(sub_name,2)

  end subroutine evaluate_vent_edema

!!!#############################################################################

  subroutine calculate_work_edema(breath_vol,dt_vol,WOBe,WOBr,pptrans)
    use arrays,only: dp,elem_nodes,node_field,&
         num_units,units_vent,unit_field
    use diagnostics, only: enter_exit
    use indices,only: nj_aw_press,nu_pe
    implicit none
    real(dp) :: breath_vol,dt_vol,WOBe,WOBr,pptrans
    ! Local variables
    integer :: ne,np1,nunit
    real(dp) :: p_resis,p_trans

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'calculate_work_edema'
    call enter_exit(sub_name,1)

    p_resis=0.0d0
    !estimate elastic and resistive WOB for each dt (sum dP.V)
    p_trans=SUM(unit_field(nu_pe,1:num_units))/num_units
    do nunit=1,num_units
       ne=units_vent(nunit)
       np1=elem_nodes(2,ne)
       p_resis=p_resis+node_field(nj_aw_press,1)-node_field(nj_aw_press,np1)
    enddo
    p_resis=p_resis/num_units
    ! vol in mm3 *1e-9=m3, pressure in Pa, hence *1d-9 = P.m3 (Joules)
    WOBe=WOBe+(p_trans-pptrans)*breath_vol*1d-9
    WOBr=WOBr+p_resis*dt_vol*1d-9

    pptrans=p_trans

    call enter_exit(sub_name,2)

  end subroutine calculate_work_edema
  
!!!####################################################################

  subroutine estimate_flow_edema(dpmus,dt,err_est)
    use arrays,only: dp,elem_field,num_units,&
         units_vent,unit_field
    use diagnostics, only: enter_exit
    use indices,only: ne_Vdot,ne_Vdot0,ne_aw_resist,nu_comp,&
         nu_dpdt,nu_Vdot0,nu_Vdot1,nu_Vdot2
    implicit none

    real(dp),intent(in) :: dpmus,dt
    real(dp),intent(out) :: err_est

    integer :: ne,nunit
    real(dp) :: alpha,beta,flow_diff,flow_sum,Q,Qinit

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'estimate_flow_edema'
    call enter_exit(sub_name,1)

    err_est = 0.0_dp
    flow_sum = 0.0_dp
  !!! For each elastic unit, calculate Qbar (equation 4.13 from Swan thesis)
    do nunit=1,num_units !for each terminal element only (with tissue units attached)
       ne=units_vent(nunit) !local element number
       ! Calculate the mean flow into the unit in the time step
       ! alpha is the rate of change of pressure at start node of terminal element
       alpha=unit_field(nu_dpdt,nunit) !dPaw/dt, initialised to zero and updated each iter
       Qinit=elem_field(ne_Vdot0,ne) !terminal element flow, updated each dt
       ! beta is the rate of change of pleural pressure
       beta = dpmus/dt ! == dPmus/dt (-ve for inspiration), updated each dt

       !      Q = C*(alpha-beta)+(Qinit-C*(alpha-beta))*DEXP(-dt/(C*R))
       Q = unit_field(nu_comp,nunit)*(alpha-beta)+ &
            (Qinit-unit_field(nu_comp,nunit)*(alpha-beta))* &
            DEXP(-dt/(unit_field(nu_comp,nunit)*elem_field(ne_aw_resist,ne)))

       unit_field(nu_Vdot2,nunit)=unit_field(nu_Vdot1,nunit) !flow at iter-2
       unit_field(nu_Vdot1,nunit)=unit_field(nu_Vdot0,nunit) !flow at iter-1

       ! flow estimate for current iter includes flow estimates at previous two iters
       !       unit_field(nu_Vdot0,nunit) = 0.75d0*unit_field(nu_Vdot2,nunit)+ &
       !            0.25d0*(Q+unit_field(nu_Vdot1,nunit))*0.5d0
       ! flow estimate for current iter includes flow estimate at previous iter
  !       unit_field(nu_Vdot0,nunit) = (Q + unit_field(nu_Vdot1,nunit))*0.5d0
  !!! from original code:
       unit_field(nu_Vdot0,nunit) = 0.75d0*unit_field(nu_Vdot2,nunit)+ &
            0.25d0*(Q+unit_field(nu_Vdot1,nunit))*0.5d0

       flow_diff=unit_field(nu_Vdot0,nunit) - elem_field(ne_Vdot,ne)
       err_est=err_est+DABS(flow_diff)**2 !sum up the error for all elements
       flow_sum=flow_sum+unit_field(nu_Vdot0,nunit)**2

  !!! ARC: DO NOT CHANGE BELOW. THIS IS NEEDED FOR THE ITERATIVE STEP
  !!! - SIMPLER OPTIONS JUST FORCE IT TO CONVERGE WHEN ITS NOT
       elem_field(ne_Vdot,ne) = (unit_field(nu_Vdot0,nunit)&
            +unit_field(nu_Vdot1,nunit))/2.0_dp
       unit_field(nu_Vdot0,nunit) = elem_field(ne_Vdot,ne)
    enddo !nunit

    ! the estimate of error for the iterative solution
    err_est=err_est/(flow_sum*DBLE(num_units))

    call enter_exit(sub_name,2)

  end subroutine estimate_flow_edema  
  
!!!###########################################################################
  subroutine read_params_evaluate_flow_edema (Gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(out) :: Gdirn
    real(dp),intent(out) :: chest_wall_compliance, constrict, COV,&
       FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target
    character,intent(out) :: expiration_type*(*)

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_evaluate_flow_edema'
    call enter_exit(sub_name,1)

    ! following values are examples from control.txt
    !    T_interval = 4.0_dp !s
    !    Gdirn = 3
    !    press_in = 0.0_dp !Pa
    !    COV = 0.2d0
    !    RMaxMean = 1.29d0
    !    RMinMean = 0.78d0
    !    i_to_e_ratio = 0.5d0 !dimensionless
    !    refvol = 0.6d0 !dimensionless
    !    volume_target = 8.d5 !mm^3  800 ml
    !    pmus_step = -5.4d0 * 98.0665_dp !-5.4 cmH2O converted to Pa
    !    expiration_type = 'passive' ! or 'active'
    !    chest_wall_compliance = 0.2d6/98.0665_dp !(0.2 L/cmH2O --> mm^3/Pa)

    open(fh, file='C:\Users\agsp886\functional-models\filtration\Parameters\params_evaluate_flow.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('FRC')
             read(buffer, *, iostat=ios) FRC
             print *, 'Read FRC: ', FRC
          case ('constrict')
             read(buffer, *, iostat=ios) constrict
             print *, 'Read constrict: ', constrict
          case ('T_interval')
             read(buffer, *, iostat=ios) T_interval
             print *, 'Read T_interval: ', T_interval
          case ('Gdirn')
             read(buffer, *, iostat=ios) Gdirn
             print *, 'Read Gdirn: ', Gdirn
          case ('press_in')
             read(buffer, *, iostat=ios) press_in
             print *, 'Read press_in: ', press_in
          case ('COV')
             read(buffer, *, iostat=ios) COV
             print *, 'Read COV: ', COV
          case ('RMaxMean')
             read(buffer, *, iostat=ios) RMaxMean
             print *, 'Read RMaxMean: ', RMaxMean
          case ('RMinMean')
             read(buffer, *, iostat=ios) RMinMean
             print *, 'Read RMinMean: ', RMinMean
          case ('i_to_e_ratio')
             read(buffer, *, iostat=ios) i_to_e_ratio
             print *, 'Read i_to_e_ratio: ', i_to_e_ratio
          case ('refvol')
             read(buffer, *, iostat=ios) refvol
             print *, 'Read refvol: ', refvol
          case ('volume_target')
             read(buffer, *, iostat=ios) volume_target
             print *, 'Read volume_target: ', volume_target
          case ('pmus_step')
             read(buffer, *, iostat=ios) pmus_step
             print *, 'Read pmus_step_coeff: ', pmus_step
          case ('expiration_type')
             read(buffer, *, iostat=ios) expiration_type
             print *, 'Read expiration_type: ', expiration_type
          case ('chest_wall_compliance')
             read(buffer, *, iostat=ios) chest_wall_compliance
             print *, 'Read chest_wall_compliance: ', chest_wall_compliance
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_evaluate_flow_edema


!!!###########################################################################

  subroutine read_params_main_edema(num_brths, num_itns, dt, err_tol)
    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(out) :: num_brths, num_itns
    real(dp) :: dt,err_tol

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'read_params_main_edema'
    call enter_exit(sub_name,1)

    ios = 0
    line = 0
    open(fh, file='C:\Users\agsp886\functional-models\filtration\Parameters\params_main.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('num_brths')
             read(buffer, *, iostat=ios) num_brths
             print *, 'Read num_brths: ', num_brths
          case ('num_itns')
             read(buffer, *, iostat=ios) num_itns
             print *, 'Read num_itns: ', num_itns
          case ('dt')
             read(buffer, *, iostat=ios) dt
             print *, 'Read dt: ', dt
          case ('err_tol')
             read(buffer, *, iostat=ios) err_tol
             print *, 'Read err_tol: ', err_tol
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_main_edema

!!!###################################################################################

  subroutine sum_elem_field_from_periphery_edema(ne_field)
    use arrays,only: dp,elem_cnct,elem_field,elem_symmetry,num_elems_aw
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: ne_field

    !Local parameters
    real(dp) :: field_value
    integer :: i,ne,ne2

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'sum_elem_field_from_periphery_edema'
    call enter_exit(sub_name,1)

    do ne=num_elems_aw,1,-1
       if(elem_cnct(1,0,ne).ne.1.or.ne.lt.70)then !not terminal
          field_value=0.d0
          do i=1,elem_cnct(1,0,ne) !for each possible daughter branch (maximum 2)
             ne2=elem_cnct(1,i,ne) !the daughter element number
             field_value=field_value+elem_symmetry(ne2)*elem_field(ne_field,ne2) !sum daughter fields
          enddo !noelem2
          elem_field(ne_field,ne)=field_value
       endif
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine sum_elem_field_from_periphery_edema  

!!!###################################################################################

  subroutine tissue_compliance_edema(chest_wall_compliance,undef)
    use arrays,only: dp,num_units,units_vent,unit_field
    use indices,only: nu_comp,nu_pe,nu_vol
    use diagnostics, only: enter_exit
    implicit none

    !     Parameter List
    real(dp), intent(in) :: chest_wall_compliance,undef

    integer :: ne,nunit
    real(dp),parameter :: a = 0.433d0, b = -0.611d0, cc = 2500.0_dp
    real(dp) :: exp_term,lambda,ratio

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_tissue_compliance_edema'
    call enter_exit(sub_name,1)

    !.....dV/dP=1/[(1/2h^2).c/2.(3a+b)exp().(4h(h^2-1)^2)+(h^2+1)/h^2)]

    do nunit=1,num_units
       ne=units_vent(nunit)
       !calculate a compliance for the tissue unit
       ratio=unit_field(nu_vol,nunit)/undef
       lambda = ratio**(1.0_dp/3.0_dp) !uniform extension ratio
       exp_term=DEXP(0.75d0*(3.0_dp*a+b)*(lambda**2-1.0_dp)**2)

       unit_field(nu_comp,nunit)=cc*exp_term/6.0_dp*(3.0_dp*(3.0_dp*a+b)**2 &
            *(lambda**2-1.0_dp)**2/lambda**2+(3.0_dp*a+b) &
            *(lambda**2+1.0_dp)/lambda**4)
       unit_field(nu_comp,nunit)=undef/unit_field(nu_comp,nunit) !in units of volume/pressure
       ! add the chest wall (proportionately) in parallel
       unit_field(nu_comp,nunit) = 1.0_dp/(1.0_dp/unit_field(nu_comp,nunit)&
            +1.0_dp/(chest_wall_compliance/DBLE(num_units)))

       !estimate an elastic recoil pressure for the unit
       unit_field(nu_pe,nunit) = cc/2.0_dp*(3.0_dp*a+b)*(lambda**2.0_dp &
            -1.0_dp)*exp_term/lambda
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine tissue_compliance_edema  
  
!!!####################################################################

  subroutine update_elem_field_edema
    use arrays,only: dp,elem_field,elem_nodes,node_xyz,num_elems_aw
    use diagnostics, only: enter_exit
    use indices,only: ne_Vdot,ne_length, &
         ne_aw_radius,ne_aw_resist,ne_t_resist,ne_vol
    use other_consts
    implicit none

    ! Local variables
    integer :: ne,np1,np2
    real(dp),parameter :: gas_density = 0.1146d-5 ! g.mm^-3
    real(dp),parameter :: gas_viscosity = 0.18d-4 ! Pa.s
    real(dp) :: gamma,resistance,reynolds,zeta
    real(dp) :: rad,le

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_elem_volume'
    call enter_exit(sub_name,1)

    do ne=1,num_elems_aw
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)

       ! element length
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)

       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_aw_radius,ne)**2 * &
            elem_field(ne_length,ne)

       le=elem_field(ne_length,ne)
       rad=elem_field(ne_aw_radius,ne)

       ! element Poiseuille (laminar) resistance in units of Pa.s.mm-3
       resistance = 8.0_dp*GAS_VISCOSITY*elem_field(ne_length,ne)/ &
            (PI*elem_field(ne_aw_radius,ne)**4) !laminar resistance

       ! element turbulent resistance (flow in bifurcating tubes)
       gamma = 0.357_dp !inspiration
       if(elem_field(ne_Vdot,ne).lt.0.0_dp) gamma = 0.46_dp !expiration

       reynolds=DABS(elem_field(ne_Vdot,ne)*2.0_dp*GAS_DENSITY/ &
            (PI*elem_field(ne_aw_radius,ne)*GAS_VISCOSITY))
       zeta = MAX(1.0_dp,dsqrt(2.0_dp*elem_field(ne_aw_radius,ne)* &
            reynolds/elem_field(ne_length,ne))*gamma)
       elem_field(ne_aw_resist,ne) = resistance * zeta

       elem_field(ne_t_resist,ne) = elem_field(ne_aw_resist,ne)

    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine update_elem_field_edema

!!!###################################################################################

  subroutine update_filtration_edema(L_p, sigma, pi_c, pi_alv, c_L, dt, Jv_sum)
  !!!   Update filtration for each elastic unit and stores it in unit_field(nu_filt,nunit)
  !     Calculate transcapillary fluid flow J_v according to Starling Equation in m^3/s
  !     J_v = L_p * S * ((P_c - P_alv) - sigma * (pi_c - pi_alv))
          
    use arrays,only: dp,num_units,unit_field, elem_nodes,node_field,units_vent
    use diagnostics, only: enter_exit,get_diagnostics_on
    use indices!,only: nu_SA,nu_rad,nu_filt,nu_filt_cleared
    use other_consts,only: PI
    implicit none

    ! In-/Output
    real(dp),intent(in) :: L_p          ! Hydraulic conductivity in mm/(s*Pa)
    real(dp),intent(in) :: sigma        ! Reflection coefficient
    real(dp),intent(in) :: pi_c         ! Capillary osmotic pressure in mmH2O
    real(dp),intent(in) :: pi_alv       ! Alveolar fluid osmotic pressure in mmH2O
    real(dp),intent(in) :: c_L          ! Lymphatic clearance/ flow in ml/(min*100g)
    real(dp),intent(in) :: dt           ! Time Step
    real(dp),intent(out) :: Jv_sum      ! Summarized Filtration in whole lung per time step
    
    ! Local variables
    !real(dp) :: J_v                    ! Filtration per time step (without lymphatic clearance)
    real(dp) :: c_L_t                   ! Lymphatic clearance/flow per time step in ml
    real(dp) :: const                   ! Constant term of Starling equation (osmotic pressures)
    
    integer :: nunit, ne, np2
    logical :: diags=.True.

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_filtration_edema'
    call enter_exit(sub_name,1)
    call get_diagnostics_on(diags)

    ! Calculate sigma * (pi_c - pi_alv) = const. in Pa
        const = sigma * (pi_c - pi_alv) * 9.80665
    
    ! Calculate Lymphatic Clearance per unit per time step c_L_t in ml
        c_L_t = c_L*dt
    
    ! Initialise summarized filtration Jv_sum over whole lung per time step
        Jv_sum = 0.0_dp
        
    ! calculate filtration J_v for all elastic units per time step (= unit_field(nu_filt,nunit))
        do  nunit=1,num_units
            ne=units_vent(nunit)
            np2=elem_nodes(2,ne)

            if(diags)write(*,*) 'blood press', nunit,unit_field(nu_blood_press,nunit),nu_blood_press
            if(unit_field(nu_blood_press,nunit).eq.0.0_dp)then
            endif
            ! calculate Filtration J_v per time step (without lymphatic clearance)
            unit_field(nu_filt,nunit) = L_p * unit_field(nu_SA,nunit) * &
                ((unit_field(nu_blood_press,nunit) - node_field(nj_aw_press,np2)) &
                - const) * dt
            
            ! Save Lymphatic Clearance in unit_field (later: calculate clearance for each unit)
            unit_field(nu_clearance,nunit) = c_L_t
            
            ! Check: Filtration > 0?
            if (unit_field(nu_filt,nunit).gt.0.0_dp)then
                ! Check: Lymphatic Clearance c_L_t > Filtration?
                if(c_L_t.gt.unit_field(nu_filt,nunit)) then
                 unit_field(nu_filt_cleared,nunit) = 0.0_dp         ! Filtration gets cleared by lymph flow
                 write(*,*) "Filtration cleared by lymphatic clearance"
                else
                 unit_field(nu_filt_cleared,nunit) = unit_field(nu_filt,nunit) - c_L_t ! Lymph flow too low to clear filtration -> edema
                 write(*,*) "Filtration not cleared by lymphatic clearance -> edema"
                endif
            else
                unit_field(nu_filt,nunit) = 0.0_dp
                unit_field(nu_filt_cleared,nunit) = 0.0_dp  
                write(*,*) 'Calculated filtration negative -> set to zero'
            endif
        
            if(diags)write(*,*) "Filtration without clearance: ", unit_field(nu_filt,nunit),unit_field(nu_blood_press,nunit),&
              node_field(nj_aw_press,np2),dt

            write(*,*) "Lymphatic Clearance: ", c_L_t
            write(*,*) "Filtration", unit_field(nu_filt,nunit)
            write(*,*) "Pcap, Palv, const, Lp*SA", unit_field(nu_blood_press,nunit), node_field(nj_aw_press,np2), const, L_p * unit_field(nu_SA,nunit)        
          
            ! Summarized Filtration for whole lung (all units) per time step (after lymphatic clearance)
                Jv_sum = unit_field(nu_filt_cleared,nunit) + Jv_sum
           
        enddo !nounit
        
    call enter_exit(sub_name,2)

  end subroutine update_filtration_edema    
  
 !!!###################################################################################
 
   subroutine calculate_terminal_exports(timestep_int)!, Jv_sum)
  !!!   Calculates the average values Pcap, Ppl (perfusion model and ventilation model), Palv,&
  !     Filtration (with/without clearance) for one whole breath and the summarized filtration &
  !     after clearance for whole lung
          
    use arrays,only: dp,num_units,unit_field, elem_nodes,node_field,units_vent
    use diagnostics, only: enter_exit,get_diagnostics_on
    use indices!,only: nu_SA,nu_rad,nu_filt,nu_filt_cleared
    implicit none

    ! In-/Output
    integer,intent(in) :: timestep_int  ! Current timestep
!    real(dp),intent(out) :: Jv_sum      ! Summarized Filtration in whole lung per time step
        
    integer :: nunit, ne, np2
    logical :: diags=.True.

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'calculate_terminal_exports'
    call enter_exit(sub_name,1)
    call get_diagnostics_on(diags)
   
    ! Initialise summarized filtration Jv_sum over whole lung per time step
!        Jv_sum = 0.0_dp
        
    ! calculate average values for all elastic units per time step (= unit_field(nu_...,nunit))
        do  nunit=1,num_units
            ne=units_vent(nunit)
            np2=elem_nodes(2,ne)
           
            ! Summarize unit values per unit of all until now calculated timesteps
                ! Sum Filtration without clearance
                unit_field(nu_filt_av,nunit) = unit_field(nu_filt_av,nunit) + unit_field(nu_filt,nunit)
                ! Sum Filtration after clearance
                unit_field(nu_filt_cleared_sum,nunit) = unit_field(nu_filt_cleared_sum,nunit) + unit_field(nu_filt_cleared,nunit)
                ! Sum capillary pressure
                unit_field(nu_blood_press_av,nunit) = unit_field(nu_blood_press_av,nunit) + unit_field(nu_blood_press,nunit)
                ! Sum pleural pressure (gradient from perfusion model)
                unit_field(nu_perfppl_av,nunit) = unit_field(nu_perfppl_av,nunit) + unit_field(nu_perfppl,nunit)
                ! Sum pleural pressure (ventilation model)
                unit_field(nu_ppl_av,nunit) = unit_field(nu_ppl_av,nunit) + unit_field(nu_ppl,nunit)
                ! Sum alveolar pressure
                node_field(nj_aw_press_av,np2) = node_field(nj_aw_press_av,np2) + node_field(nj_aw_press,np2)
                            
        enddo !nounit
        
        ! Calculate average values for terminal export file of all until now calculated timesteps
            ! Average Filtration without clearance
            unit_field(nu_filt_av,:) = unit_field(nu_filt_av,:)/timestep_int
            ! Average Filtration after clearance
            unit_field(nu_filt_cleared_av,:) = unit_field(nu_filt_cleared_sum,:)/timestep_int
            ! Sum capillary pressure
            unit_field(nu_blood_press_av,:) = unit_field(nu_blood_press_av,:)/timestep_int
            ! Average pleural pressure (gradient from perfusion model)
            unit_field(nu_perfppl_av,:) = unit_field(nu_perfppl_av,:)/timestep_int
            ! Sum pleural pressure (ventilation model)
            unit_field(nu_ppl_av,:) = unit_field(nu_ppl_av,:)/timestep_int
            ! Sum alveolar pressure
            node_field(nj_aw_press_av,:) = node_field(nj_aw_press_av,:)/timestep_int
            
            
            
        
    call enter_exit(sub_name,2)

  end subroutine calculate_terminal_exports  
  
 !!!################################################################################### 
  
  
  subroutine update_node_pressures_edema(press_in)
    use arrays,only: dp,elem_field,elem_nodes,node_field,num_elems_aw
    use indices,only: ne_Vdot,ne_aw_resist,nj_aw_press
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(in) :: press_in
    !Local parameters
    integer :: ne,np1,np2

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_node_pressures_edema'
    call enter_exit(sub_name,1)

!!! Use the known resistances and flows to calculate nodal pressures through whole tree

    ! set the initial node pressure to be the input pressure (usually zero)
    ne=1 !element number at top of tree, usually = 1
    np1=elem_nodes(1,ne) !first node in element
    node_field(nj_aw_press,np1)=press_in !set pressure at top of tree

    do ne=1,num_elems_aw !for each element
       np1=elem_nodes(1,ne) !start node number
       np2=elem_nodes(2,ne) !end node number
       !P(np2) = P(np1) - Resistance(ne)*Flow(ne)
       node_field(nj_aw_press,np2)=node_field(nj_aw_press,np1) &
            -elem_field(ne_aw_resist,ne)*elem_field(ne_Vdot,ne)
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine update_node_pressures_edema  

!!!###################################################################################

  subroutine update_mean_pleural_pressure_edema(ppl_current)
  !!! Update the mean pleural pressure based on current Pel (=Ptp) and Palv,
  !!! i.e. Ppl(unit) = -Pel(unit)+Palv(unit)
    use arrays,only: dp,elem_nodes,node_field,num_units,units_vent,unit_field
    use diagnostics, only: enter_exit
    use indices,only: nj_aw_press,nu_pe
    implicit none

    real(dp),intent(out) :: ppl_current
    integer :: ne,np2,nunit

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_mean_pleural_pressure_edema'
    call enter_exit(sub_name,1)

    ppl_current = 0.0_dp
    do nunit=1,num_units
       ne=units_vent(nunit)
       np2=elem_nodes(2,ne)
       ppl_current = ppl_current - unit_field(nu_pe,nunit) + &
            node_field(nj_aw_press,np2)
    enddo !noelem
    ppl_current = ppl_current/num_units

    call enter_exit(sub_name,2)

  end subroutine update_mean_pleural_pressure_edema  


!!!###################################################################################

  subroutine update_pleural_pressure_edema!(ppl_edema)
  !!! Update the pleural pressure based on current Pel (=Ptp) and Palv,
  !!! i.e. Ppl(unit) = -Pel(unit) + Palv(unit)
    use arrays,only: dp,elem_nodes,node_field,num_units,units_vent,unit_field
    use diagnostics, only: enter_exit
    use indices,only: nj_aw_press,nu_pe,nu_ppl
    implicit none

    !real(dp),intent(out) :: ppl_edema(:)
    integer :: ne,np2,nunit
    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_pleural_pressure_edema'
    call enter_exit(sub_name,1)

    do nunit=1,num_units
       ne=units_vent(nunit)
       np2=elem_nodes(2,ne)
       unit_field(nu_ppl,nunit) = abs(- unit_field(nu_pe,nunit) + &
            node_field(nj_aw_press,np2))
    enddo !noelem
    call enter_exit(sub_name,2)

  end subroutine update_pleural_pressure_edema 

!!!###################################################################################
  
  subroutine update_proximal_pressure_edema
    use arrays,only: elem_nodes,node_field,num_units,units_vent,unit_field
    use indices,only: nj_aw_press,nu_air_press
    use diagnostics, only: enter_exit
    implicit none

    integer :: ne,np1,nunit

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_proximal_pressure_edema'
    call enter_exit(sub_name,1)

    ! update the pressure at the proximal node of the element that feeds the elastic unit
    do nunit=1,num_units
       ne=units_vent(nunit)
       np1=elem_nodes(1,ne)
       unit_field(nu_air_press,nunit)=node_field(nj_aw_press,np1) !store the pressure at entry node
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine update_proximal_pressure_edema  

!!!###################################################################################

  subroutine update_surface_area_edema
  !!! Update surface area for each elastic unit and stores it in unit_field(nu_SA,nunit) in mikroliter
    use arrays,only: dp,num_units,unit_field
    use diagnostics, only: enter_exit
    use indices!,only: nu_SA,nu_rad
    use other_consts,only: PI
    implicit none

    integer :: nunit

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_surface_area_edema'
    call enter_exit(sub_name,1)
    
    ! calculate radii and surface areas for all elastic units depending on unit volume (= unit_field(nu_vol,nunit))
    do nunit=1,num_units
        
       ! calculate inner radius for each unit in millimeter
       unit_field(nu_rad,nunit) = (3/(4*PI) * unit_field(nu_vol,nunit))**(1.0/3.0)

       !!write(*,*) "Radius unit nu"
       !!write(*,*) unit_field(nu_rad,nunit)
       
       ! calculate surface area for each unit in millimeter^2
       unit_field(nu_SA,nunit) = 4.0 * PI * unit_field(nu_rad,nunit)**2.0
       
       !!write(*,*) "Surface area unit nu"
       !!write(*,*) unit_field(nu_SA,nunit)
       
    enddo !nounit

    call enter_exit(sub_name,2)

  end subroutine update_surface_area_edema    
  
 !!!###################################################################################

  subroutine update_unit_dpdt_edema(dt)
    use arrays,only: dp,elem_nodes,node_field,&
         num_units,units_vent,unit_field
    use indices,only: nj_aw_press,nu_dpdt,nu_air_press
    use diagnostics, only: enter_exit
    implicit none

    real(dp), intent(in) :: dt
    integer :: ne,np1,nunit
    real(dp) :: est

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_dpdt_edema'
    call enter_exit(sub_name,1)

! this is the rate of change of pressure at the proximal end of the element
! that supplies the tissue unit, i.e. not the rate of change of pressure within the unit.
    do nunit=1,num_units
       ne=units_vent(nunit)
       np1=elem_nodes(1,ne)
       ! linear estimate
       est=(node_field(nj_aw_press,np1) &
            -unit_field(nu_air_press,nunit))/dt
       ! weight new estimate with the previous dP/dt
       unit_field(nu_dpdt,nunit)=0.5d0*(est+unit_field(nu_dpdt,nunit))
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine update_unit_dpdt_edema  
  
!!!###################################################################################

  subroutine update_unit_volume_edema(dt,Tinsp, Texpn)
    use arrays,only: dp,elem_field,elem_nodes,num_units,&
         units_vent,unit_field
    use diagnostics, only: enter_exit
    use indices,only: ne_Vdot,nu_vol,nu_vt,nu_vent
    implicit none

    real(dp),intent(in) :: dt,Tinsp,Texpn
    integer :: ne,np,nunit

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'update_unit_volume_edema'
    call enter_exit(sub_name,1)

    do nunit=1,num_units
       ne=units_vent(nunit)
       np=elem_nodes(2,ne)
       ! update the volume of the lumped tissue unit
       unit_field(nu_vol,nunit)=unit_field(nu_vol,nunit)+dt* &
            elem_field(ne_Vdot,ne) !in mm^3
       if(elem_field(ne_Vdot,1).gt.0.0_dp)then  !only store inspired volume
          unit_field(nu_vt,nunit)=unit_field(nu_vt,nunit)+dt* &
               elem_field(ne_Vdot,ne)
        unit_field(nu_vent,nunit)=unit_field(nu_vt,nunit)/(TInsp+Texpn) ! volume  per secondto alvolus
       endif
    enddo !nunit

    call enter_exit(sub_name,2)

  end subroutine update_unit_volume_edema

!!!#################################################################################
  
end module filtration
