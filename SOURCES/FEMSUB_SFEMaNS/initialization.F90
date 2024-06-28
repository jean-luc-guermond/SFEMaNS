!
!Authors Jean-Luc Guermond, Caroline Nore, Copyrights 2005
!Revised June 2008, Jean-Luc Guermond
!Revised for PETSC, Jean-Luc Guermond, Francky Luddens, January 2011
!Revised July 9th 2013, JLG, Loic Cappanera, Remi Menard, Daniel Castanon
!Revised July 25th 2016, JLG, Loic Cappanera, Raphael Zanella
!Revised July 20th 2019, JLG, Hughes Faller (p3/p2 FE + Taylor algorithm for NS)

MODULE initialization
  USE def_type_mesh
  USE symmetric_field
  USE input_data
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
  USE fourier_to_real_for_vtu
  IMPLICIT NONE
  PUBLIC:: initial, save_run, run_SFEMaNS
  PUBLIC:: prodmat_maxwell_int_by_parts
  PRIVATE

  !Logicals for equations-----------------------------------------------------
  LOGICAL                                         :: if_momentum, if_mass, if_induction, if_energy
  LOGICAL                                         :: if_concentration

  !Fields for Navier-Stokes---------------------------------------------------
  TYPE(mesh_type), TARGET                         :: pp_mesh, vv_mesh
  TYPE(petsc_csr_LA)                              :: vv_1_LA, pp_1_LA
  TYPE(petsc_csr_LA)                              :: vv_3_LA   ! for stress bc
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: un, un_m1
  ! CN-HF 16/01/2020
  TYPE(dyn_real_array_three), TARGET, ALLOCATABLE, DIMENSION(:)   :: der_un
  ! (noeuds,type,mode) composante du champ de vitesse a deux instants sur vv_mesh
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: pn, pn_m1
  TYPE(dyn_real_array_three), ALLOCATABLE, DIMENSION(:)   :: der_pn
  ! (noeuds,type,mode) composante du champ de pression a deux instants sur pp_mesh
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)          :: incpn, incpn_m1
  !---------------------------------------------------------------------------

  !Fields for level sets in Navier-Stokes-------------------------------------
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:,:):: level_set, level_set_m1
  !Fields for density---------------------------------------------------------
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: density, density_m1, density_m2
  !Entropy viscosity for level-set--------------------------------------------
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:)            :: visc_entro_level
  !Maximum of velocity--------------------------------------------------------
  REAL(KIND=8)                                         :: max_vel
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:,:):: visc_LES
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:,:):: visc_LES_level

  !Fields for temperature-----------------------------------------------------
  ! (noeuds,type,mode) composante du champ de phase a deux instants sur vv_mesh
  TYPE(mesh_type), TARGET                              :: temp_mesh
  TYPE(petsc_csr_LA)                                   :: temp_1_LA
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: tempn, tempn_m1
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:)      :: vol_heat_capacity_field
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:)      :: temperature_diffusivity_field
  !---------------------------------------------------------------------------
  !Fields for concentration-----------------------------------------------------
  TYPE(mesh_type), TARGET                              :: conc_mesh
  TYPE(petsc_csr_LA)                                   :: conc_1_LA
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:)  :: concn, concn_m1
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:)      :: concentration_diffusivity_field


  !Fields for Maxwell---------------------------------------------------------
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: Hn, Hn1, Hext, phin, phin1
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: Bn, Bn1, Bext
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:)     :: sigma_field, mu_H_field
  TYPE(mesh_type), TARGET                             :: H_mesh, phi_mesh, pmag_mesh
  TYPE(petsc_csr_LA)                                  :: LA_H, LA_pmag, LA_phi, LA_mhd
  TYPE(interface_type), TARGET                        :: interface_H_mu, interface_H_phi
  !---------------------------------------------------------------------------

  !Periodic structures--------------------------------------------------------
  TYPE(periodic_type)                                 :: H_phi_per
  TYPE(periodic_type)                                 :: vvrt_per
  TYPE(periodic_type)                                 :: vvrtz_per
  TYPE(periodic_type)                                 :: vvz_per
  TYPE(periodic_type)                                 :: pp_per
  TYPE(periodic_type)                                 :: temp_per
  TYPE(periodic_type)                                 :: level_set_per
  TYPE(periodic_type)                                 :: conc_per
  !---------------------------------------------------------------------------

  !Coupling variables---------------------------------------------------------
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: v_to_Max
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: H_to_NS
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: B_to_NS
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: T_to_H
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: j_H_to_T
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: j_H_to_conc
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: conc_to_v
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: conc_to_H
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: v_to_conc
  !October 7, 2008, JLG
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:)         :: T_to_NS
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: v_to_energy
  REAL(KIND=8), TARGET, ALLOCATABLE, DIMENSION(:,:,:) :: H_to_energy, pdt_H_to_energy
  !---------------------------------------------------------------------------

  !Connectivity structures----------------------------------------------------
  !October 7, 2008, JLG
  INTEGER, ALLOCATABLE, DIMENSION(:)                  :: jj_v_to_H
  INTEGER, ALLOCATABLE, DIMENSION(:)                  :: jj_v_to_temp
  INTEGER, ALLOCATABLE, DIMENSION(:)                  :: jj_T_to_H
  INTEGER, ALLOCATABLE, DIMENSION(:)                  :: jj_c_to_v
  INTEGER, ALLOCATABLE, DIMENSION(:)                  :: jj_c_to_H
  !---------------------------------------------------------------------------

  !Modes list-----------------------------------------------------------------
  INTEGER,      TARGET, ALLOCATABLE, DIMENSION(:) :: list_mode
  INTEGER                                         :: m_max_c
  !---------------------------------------------------------------------------

  !Names already used---------------------------------------------------------
  REAL(KIND=8)                                    :: time
  REAL(KIND=8)                                    :: R_fourier
  INTEGER                                         :: index_fourier
  !---------------------------------------------------------------------------

  !Communicators for Petsc, in space and Fourier------------------------------
  MPI_Comm                                        :: comm_cart
  MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d, comm_one_d_ns
  MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d_temp, coord_cart
  MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d_conc
  !---------------------------------------------------------------------------

  !-------------END OF DECLARATIONS------------------------------------------
CONTAINS
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE initial(vv_mesh_out, pp_mesh_out, H_mesh_out, phi_mesh_out, temp_mesh_out, &
       conc_mesh_out,interface_H_phi_out, interface_H_mu_out, list_mode_out, &
       un_out, pn_out, Hn_out, Bn_out, phin_out, v_to_Max_out, &
       vol_heat_capacity_field_out, temperature_diffusivity_field_out, concentration_diffusivity_field_out, &
       mu_H_field_out, sigma_field_out, &
       time_out, m_max_c_out, comm_one_d_out, comm_one_d_ns_out, comm_one_d_temp_out, comm_one_d_conc_out, &
       tempn_out, concn_out, level_set_out, density_out, der_un_out, visc_LES_out, visc_LES_level_out)
    USE fourier_to_real_for_vtu
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                  :: pp_mesh_out, vv_mesh_out
    TYPE(mesh_type), POINTER                  :: H_mesh_out, phi_mesh_out
    TYPE(mesh_type), POINTER                  :: temp_mesh_out
    TYPE(mesh_type), POINTER                  :: conc_mesh_out
    TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un_out
    TYPE(interface_type), POINTER             :: interface_H_mu_out, interface_H_phi_out
    INTEGER,      POINTER,  DIMENSION(:)      :: list_mode_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)  :: un_out, pn_out, Hn_out, Bn_out, phin_out, v_to_Max_out, tempn_out, density_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)  :: concn_out
    REAL(KIND=8), POINTER,  DIMENSION(:)      :: concentration_diffusivity_field_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:):: level_set_out
    REAL(KIND=8), POINTER,  DIMENSION(:)      :: sigma_field_out, mu_H_field_out
    REAL(KIND=8), POINTER,  DIMENSION(:)      :: vol_heat_capacity_field_out, temperature_diffusivity_field_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:):: visc_LES_out
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:):: visc_LES_level_out
    REAL(KIND=8)                              :: time_out
    INTEGER                                   :: m_max_c_out
    MPI_Comm, DIMENSION(:), POINTER           :: comm_one_d_out, comm_one_d_ns_out
    MPI_Comm, DIMENSION(:), POINTER           :: comm_one_d_temp_out, comm_one_d_conc_out

    CALL INIT

    !===Initialize meshes for vtu post processing
    CALL sfemans_initialize_postprocessing(comm_one_d, vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, &
         conc_mesh, list_mode, inputs%number_of_planes_in_real_space)
    vv_mesh_out => vv_mesh
    pp_mesh_out => pp_mesh
    H_mesh_out => H_mesh
    phi_mesh_out => phi_mesh
    temp_mesh_out => temp_mesh
    interface_H_mu_out => interface_H_mu
    interface_H_phi_out => interface_H_phi
    list_mode_out => list_mode
    un_out   => un
    der_un_out   => der_un
    pn_out   => pn
    tempn_out   => tempn
    level_set_out => level_set
    density_out => density
    visc_LES_out => visc_LES
    visc_LES_level_out => visc_LES_level
    Hn_out   => Hn
    Bn_out   => Bn
    phin_out => phin
    v_to_Max_out => v_to_Max
    vol_heat_capacity_field_out => vol_heat_capacity_field
    temperature_diffusivity_field_out => temperature_diffusivity_field
    mu_H_field_out => mu_H_field
    sigma_field_out => sigma_field
    time_out = time
    m_max_c_out = m_max_c
    comm_one_d_out => comm_one_d
    comm_one_d_ns_out => comm_one_d_ns
    comm_one_d_temp_out => comm_one_d_temp
    conc_mesh_out => conc_mesh
    concn_out => concn
    concentration_diffusivity_field_out => concentration_diffusivity_field
    comm_one_d_conc_out => comm_one_d_conc

  END SUBROUTINE initial
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE run_SFEMaNS(time_in, it)
    USE subroutine_mass
    USE update_temperature
    USE subroutine_concentration
    USE update_navier_stokes
    USE update_maxwell
    USE update_taylor_navier_stokes
    USE input_data
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN)                                  :: time_in
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: visco_dyn
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: one_over_sigma_ns_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: heat_density_ns_m1
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: heat_density_ns
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: heat_density_ns_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))     :: heat_diffusivity_ns
    INTEGER                                                   :: it
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode))      :: j_Hn
    REAL(KIND=8), DIMENSION(SIZE(level_set,1),SIZE(level_set,2),2,SIZE(list_mode)) :: level_set_reg

    CALL zero_out_modes

    time = time_in

    IF (if_mass) THEN
       IF (inputs%variation_temp_param_fluid) THEN
          CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set_m1, &
               inputs%heat_capacity_fluid*inputs%density_fluid, heat_density_ns_m1)
       ELSE
          heat_density_ns_m1  = 0.d0
       END IF

       CALL three_level_mass(comm_one_d_ns, time, pp_1_LA, vv_1_LA, list_mode, pp_mesh, vv_mesh, &
            2*un-un_m1, max_vel, level_set_per, density_m2, density_m1, density, level_set_m1, level_set,&
            visc_entro_level, level_set_reg, visc_LES_level)
       CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set_m1, &
            inputs%dyna_visc_fluid, visco_dyn)

       IF (inputs%variation_sigma_fluid) THEN
          CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set, &
               1.d0/inputs%sigma_fluid, one_over_sigma_ns_p1)
       ELSE
          one_over_sigma_ns_p1 = 0.d0
       END IF

       IF (inputs%variation_temp_param_fluid) THEN
          CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set_m1, &
               inputs%heat_capacity_fluid*inputs%density_fluid, heat_density_ns)
          CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set, &
               inputs%heat_capacity_fluid*inputs%density_fluid, heat_density_ns_p1)
          CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set, &
               inputs%heat_diffu_fluid, heat_diffusivity_ns)
       ELSE
          heat_density_ns     = 0.d0
          heat_density_ns_p1  = 0.d0
          heat_diffusivity_ns = 0.d0
       END IF
    END IF

    IF (if_concentration) THEN
       IF (if_momentum) THEN
          CALL projection_velocity(conc_mesh, 2*un-un_m1, jj_c_to_v, .TRUE., v_to_conc)
       END IF

       IF (if_induction) THEN
          CALL compute_rot_h(2*Hn-Hn1, j_Hn)
          !j_Hn = Hn
          CALL projection_mag_field(conc_mesh, j_Hn, jj_c_to_H, .TRUE., j_H_to_conc)
       END IF
       CALL three_level_concentration(comm_one_d_conc, time, conc_1_LA, inputs%dt, list_mode, &
            conc_mesh, concn_m1, concn, v_to_conc,& ! H_to_energy,&
            concentration_diffusivity_field, inputs%my_par_concentration,&
            inputs%concentration_list_dirichlet_sides, inputs%concentration_list_robin_sides, &
            inputs%convection_coeff_conc_lhs, inputs%convection_coeff_conc_rhs,inputs%exterior_concentration, &
            conc_per, j_H_to_conc)
    END IF

    IF (if_energy) THEN
       IF (if_momentum) THEN
          CALL projection_velocity(temp_mesh, 2*un-un_m1, jj_v_to_temp, .FALSE., v_to_energy)
       END IF
       IF (inputs%type_pb=='fhd') THEN
          CALL projection_mag_field(vv_mesh, 2*Hn-Hn1, jj_v_to_H, .TRUE., H_to_NS)
          CALL projection_mag_field(temp_mesh, H_to_NS, jj_v_to_temp, .FALSE., H_to_energy)
          IF (.NOT. inputs%if_steady_current_fhd) THEN
             CALL projection_mag_field(vv_mesh, (Hn-Hn1)/inputs%dt, jj_v_to_H, .TRUE., H_to_NS)
             CALL projection_mag_field(temp_mesh, H_to_NS, jj_v_to_temp, .FALSE., pdt_H_to_energy)
          END IF
       END IF
       CALL temperature_decouple(comm_one_d_temp, time, temp_1_LA, list_mode, &
            temp_mesh, tempn_m1, tempn, v_to_energy, H_to_energy, pdt_H_to_energy, &
            vol_heat_capacity_field, temperature_diffusivity_field, temp_per, &
            heat_density_ns_m1, heat_density_ns, heat_density_ns_p1, heat_diffusivity_ns, jj_v_to_temp)
    END IF

    IF (if_momentum) THEN
       IF (if_energy) THEN
          CALL projection_temperature(vv_mesh, tempn, jj_v_to_temp, .TRUE., T_to_NS)
       END IF
       IF (if_concentration) THEN
          CALL projection_concentration(vv_mesh, concn, jj_c_to_v, conc_to_v)
       END IF
       IF (if_induction) THEN
          CALL projection_mag_field(vv_mesh, 2*Hn-Hn1, jj_v_to_H, .TRUE., H_to_NS)
          CALL projection_mag_field(vv_mesh, 2*Bn-Bn1, jj_v_to_H, .TRUE., B_to_NS)
       END IF
       !===JLG July 20, 2019, p3 mesh
       !===HF April 2019
       IF (inputs%if_navier_stokes_with_taylor) THEN
          CALL navier_stokes_taylor(comm_one_d_ns, time, vv_3_LA, pp_1_LA, &
               list_mode, pp_mesh, vv_mesh, pn, der_pn, un, der_un, vvz_per, &
               pp_per, density, tempn, concn)
       ELSE
          CALL navier_stokes_decouple(comm_one_d_ns,time, vv_3_LA, pp_1_LA, &
               list_mode, pp_mesh, vv_mesh, incpn_m1, incpn, &
               pn_m1, pn, un_m1, un, vvz_per, pp_per, H_to_NS, B_to_NS, &
               density_m2, density_m1, density, visco_dyn, T_to_NS, conc_to_v, &
               level_set_m1, level_set, visc_entro_level, level_set_reg, visc_LES)
       END IF
       !===HF April 2019
       !===JLG July 20, 2019, p3 mesh
    END IF

    IF (if_induction) THEN
       IF (inputs%type_pb == 'fhd' .AND. inputs%if_steady_current_fhd) THEN !===If steady current, computation of H only once
          IF (it == 1) THEN
             IF (if_momentum) THEN
                CALL projection_velocity(H_mesh, un, jj_v_to_H, .FALSE., v_to_Max)
             END IF
             CALL maxwell_decouple(comm_one_d, H_mesh, pmag_mesh, phi_mesh, &
                  interface_H_phi, interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, v_to_Max, &
                  inputs%stab, inputs%stab_jump_h,sigma_field, R_fourier, index_fourier, mu_H_field, inputs%mu_phi, &
                  time, inputs%dt, inputs%Rem, list_mode, H_phi_per, LA_H, LA_pmag, LA_phi, &
                  LA_mhd, one_over_sigma_ns_p1, jj_v_to_H,conc_to_H)
             Hn1 = Hn
             Bn1 = Bn
             phin1 = phin
          END IF
       ELSE
          IF (if_momentum) THEN
             CALL projection_velocity(H_mesh, un, jj_v_to_H, .FALSE., v_to_Max)
          END IF
          IF (if_energy) THEN
             CALL projection_temperature(H_mesh, 2*tempn-tempn_m1, jj_T_to_H, .FALSE.,&
                  T_to_H)
          END IF
          IF (if_concentration) THEN
             CALL projection_concentration(H_mesh, 2*concn-concn_m1, jj_c_to_H,conc_to_H)
          END IF
          CALL maxwell_decouple(comm_one_d, H_mesh, pmag_mesh, phi_mesh, &
               interface_H_phi, interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, v_to_Max, &
               inputs%stab, inputs%stab_jump_h, sigma_field, R_fourier, index_fourier, mu_H_field, inputs%mu_phi, &
               time, inputs%dt, inputs%Rem, list_mode, H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns_p1, &
               jj_v_to_H, conc_to_H)
       END IF
    END IF

  END SUBROUTINE run_SFEMaNS
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE zero_out_modes
    USE input_data
    IMPLICIT NONE
    LOGICAL,                                  SAVE            :: once_zero_out_mode=.TRUE.
    INTEGER,      POINTER,  DIMENSION(:),     SAVE            :: select_mode_ns, select_mode_mxw
    INTEGER                                                   :: kp
    IF (once_zero_out_mode) THEN
       once_zero_out_mode=.FALSE.
       IF (inputs%if_zero_out_modes) THEN
          CALL prepare_zero_out_modes(list_mode, inputs%list_select_mode_ns, select_mode_ns)
          CALL prepare_zero_out_modes(list_mode, inputs%list_select_mode_mxw, select_mode_mxw)
       END IF
    END IF
    IF (inputs%if_zero_out_modes) THEN
       IF (H_mesh%me /=0 .AND. SIZE(select_mode_mxw)>0) THEN
          Hn(:,:,select_mode_mxw) = 0.d0
          Hn1(:,:,select_mode_mxw) = 0.d0
       END IF
       IF (phi_mesh%me /= 0 .AND. SIZE(select_mode_mxw)>0) THEN
          phin(:,:,select_mode_mxw) = 0.d0
          phin1(:,:,select_mode_mxw) = 0.d0
       END IF
       IF (vv_mesh%me /=0 .AND. SIZE(select_mode_ns)>0) THEN
          un(:,:,select_mode_ns) = 0.d0
          pn(:,:,select_mode_ns) = 0.d0
          IF(inputs%if_navier_stokes_with_taylor) THEN
             DO kp = 1, inputs%taylor_order-1
                der_un(kp)%DRT( :,:,select_mode_ns) = 0.d0
                der_pn(kp)%DRT( :,:,select_mode_ns) = 0.d0
             END DO
          ELSE
             incpn(:,:,select_mode_ns) = 0.d0
             un_m1(:,:,select_mode_ns) = 0.d0
             pn_m1(:,:,select_mode_ns) = 0.d0
             incpn_m1(:,:,select_mode_ns) = 0.d0
          END IF
       END IF
    END IF
  END SUBROUTINE zero_out_modes
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE prepare_zero_out_modes(list_mode, list_mode_to_zero_out, select_mode)
    USE input_data
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: list_mode_to_zero_out
    INTEGER,          DIMENSION(:) :: list_mode
    INTEGER, POINTER, DIMENSION(:) :: select_mode
    INTEGER               :: i, inc
    INTEGER, DIMENSION(1) :: kloc
    inc = 0
    DO i = 1, SIZE(list_mode_to_zero_out)
       IF (MINVAL(ABS(list_mode-list_mode_to_zero_out(i)))==0) THEN
          inc = inc + 1
       END IF
    END DO
    ALLOCATE(select_mode(inc))
    inc = 0
    DO i = 1, SIZE(list_mode_to_zero_out)
       IF (MINVAL(ABS(list_mode-list_mode_to_zero_out(i)))==0) THEN
          inc = inc + 1
          kloc = MINLOC(ABS(list_mode-list_mode_to_zero_out(i)))
          select_mode(inc) = kloc(1)
       END IF
    END DO
  END SUBROUTINE prepare_zero_out_modes
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE projection_velocity(mesh, vn, connectivity_structure, if_restriction, coupling_variable)
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                               :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
    INTEGER, DIMENSION(:), INTENT(IN)                         :: connectivity_structure
    LOGICAL, INTENT(IN)                                       :: if_restriction
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
    INTEGER                                                   :: i, j, k

    IF (if_restriction) THEN
       IF (mesh%me /=0 ) THEN !===construction of the restricted velocity field
          DO j = 1, SIZE(connectivity_structure)
             IF (connectivity_structure(j) == -1) CYCLE
             coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
          END DO
       END IF
    ELSE
       IF (mesh%np>0) THEN !===Check that ns is a subset of temp before constructing the extended vv field
          DO i = 1, m_max_c
             DO k= 1, 6 !===The user has to code extension_vel
                coupling_variable(:,k,i) = extension_velocity(k, mesh, list_mode(i), time, 1)
             END DO
          END DO
       END IF
       IF (mesh%me /=0) THEN !===construction of the extended field
          DO j = 1, SIZE(connectivity_structure)
             IF (connectivity_structure(j) == -1) CYCLE
             coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
          END DO
       END IF
    END IF
  END SUBROUTINE projection_velocity
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  SUBROUTINE projection_concentration(mesh, vn, connectivity_structure, coupling_variable)
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                               :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
    INTEGER, DIMENSION(:), INTENT(IN)                         :: connectivity_structure
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
    INTEGER                                                   :: i,j,k

    !IF (mesh%me /=0 ) THEN !===construction of the restricted  concentration field
    !   DO j = 1, SIZE(connectivity_structure)
    !      IF (connectivity_structure(j) == -1) CYCLE
    !      coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
    !  END DO
    !END IF
    IF (mesh%np>0) THEN !===Check that temp is a subset of H before constructing the extended vv field
       DO i = 1, m_max_c
          DO k= 1, 2 !===The user has to code extension_concentration
             coupling_variable(:,k,i) = extension_concentration(k, mesh, list_mode(i), time, 1)
          END DO
       END DO
    END IF
    IF (mesh%me /=0) THEN !===construction of the extended field
       DO j = 1, SIZE(connectivity_structure)
          IF (connectivity_structure(j) == -1) CYCLE
          coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
       END DO
    END IF
  END SUBROUTINE projection_concentration
  !---------------------------------------------------------------------------

  SUBROUTINE projection_temperature(mesh, vn, connectivity_structure, if_restriction, coupling_variable)
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                               :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
    INTEGER, DIMENSION(:), INTENT(IN)                         :: connectivity_structure
    LOGICAL, INTENT(IN)                                       :: if_restriction
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
    INTEGER                                                   :: i,j,k

    IF (if_restriction) THEN
       IF (mesh%me /=0 ) THEN !===construction of the restricted temperature field
          DO j = 1, SIZE(connectivity_structure)
             IF (connectivity_structure(j) == -1) CYCLE
             coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
          END DO
       END IF
    ELSE
       IF (mesh%np>0) THEN !===Check that temp is a subset of H before constructing the extended vv field
          DO i = 1, m_max_c
             DO k= 1, 2 !===The user has to code extension_temperature
                coupling_variable(:,k,i) = extension_temperature(k, mesh, list_mode(i), time, 1)
             END DO
          END DO
       END IF
       IF (mesh%me /=0) THEN !===construction of the extended field
          DO j = 1, SIZE(connectivity_structure)
             IF (connectivity_structure(j) == -1) CYCLE
             coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
          END DO
       END IF
    END IF

  END SUBROUTINE projection_temperature
  !---------------------------------------------------------------------------
  SUBROUTINE projection_mag_field(mesh, vn, connectivity_structure, if_restriction, coupling_variable)
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                               :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
    INTEGER, DIMENSION(:), INTENT(IN)                         :: connectivity_structure
    LOGICAL, INTENT(IN)                                       :: if_restriction
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
    INTEGER                                                   :: j, m

    DO j = 1, SIZE(connectivity_structure)
       IF (connectivity_structure(j) == -1) CYCLE
       IF (if_restriction) THEN
          coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
       ELSE
          coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
       END IF
    END DO
    IF (if_restriction) THEN
       IF (H_mesh%gauss%n_w/=mesh%gauss%n_w) THEN
          DO m = 1, mesh%me
             coupling_variable(mesh%jj(4,m),:,:) = (coupling_variable(mesh%jj(2,m),:,:) &
                  + coupling_variable(mesh%jj(3,m),:,:))/2
             coupling_variable(mesh%jj(5,m),:,:) = (coupling_variable(mesh%jj(3,m),:,:) &
                  + coupling_variable(mesh%jj(1,m),:,:))/2
             coupling_variable(mesh%jj(6,m),:,:) = (coupling_variable(mesh%jj(1,m),:,:) &
                  + coupling_variable(mesh%jj(2,m),:,:))/2
          END DO
       END IF
    END IF

  END SUBROUTINE projection_mag_field
  !---------------------------------------------------------------------------

  !----------------SAVE RUN---------------------------------------------------
  SUBROUTINE save_run(it, freq_restart)
    USE restart
!TEST LC LES_SUITE 2024/06
    USE subroutine_compute_visc_LES_level
!TEST LC LES_SUITE 2024/06
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: it, freq_restart
!TEST LC LES_SUITE 2024/06
    INTEGER :: n
!TEST LC LES_SUITE 2024/06

    IF (if_momentum) THEN
       IF (pp_mesh%me /= 0) THEN
          IF(inputs%if_navier_stokes_with_taylor) THEN
             CALL write_restart_ns_taylor(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                  list_mode, un, der_un, pn, der_pn, inputs%file_name, it, freq_restart)
          ELSE
             IF (.NOT. if_mass) THEN
                CALL write_restart_ns(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                     list_mode, un, un_m1, pn, pn_m1, &
                     incpn, incpn_m1, inputs%file_name, it, freq_restart)
                IF (inputs%LES) THEN
                   CALL write_restart_LES(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                        list_mode, inputs%file_name, it, freq_restart, &
                        opt_LES_NS=visc_LES)
                END IF
             ELSE
                CALL write_restart_ns(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                     list_mode, un, un_m1, pn, pn_m1, &
                     incpn, incpn_m1, inputs%file_name, it, freq_restart, &
                     opt_level_set=level_set, opt_level_set_m1=level_set_m1,opt_max_vel=max_vel)

!TEST LC LES_SUITE 2024/06
                IF (inputs%if_level_set_P2) THEN
                   DO n = 1, inputs%nb_fluid-1
                      CALL compute_visc_LES_level(comm_one_d_ns, time, vv_1_LA, list_mode, vv_mesh, &
                           level_set_m1(n,:,:,:), level_set(n,:,:,:), inputs%my_par_level_set, &
                           inputs%level_set_list_dirichlet_sides, level_set_per, n, &
                           visc_entro_level, visc_LES_level(n,:,:,:))
                   END DO
                ELSE
                   DO n = 1, inputs%nb_fluid-1
                      CALL compute_visc_LES_level(comm_one_d_ns, time, pp_1_LA, list_mode, pp_mesh, &
                           level_set_m1(n,:,:,:), level_set(n,:,:,:), inputs%my_par_level_set, &
                           inputs%level_set_list_dirichlet_sides, level_set_per, n, &
                           visc_entro_level, visc_LES_level(n,:,:,:))
                   END DO
                END IF

                IF (inputs%LES .AND. inputs%if_LES_in_momentum) THEN
                   CALL write_restart_LES(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                        list_mode, inputs%file_name, it, freq_restart, &
                        opt_LES_NS=visc_LES, opt_LES_level=visc_LES_level)
                ELSE
                   CALL write_restart_LES(comm_one_d_ns, vv_mesh, pp_mesh, time, &
                        list_mode, inputs%file_name, it, freq_restart, &
                        opt_LES_level=visc_LES_level)
                END IF
!TEST LC LES_SUITE 2024/06
             END IF
          END IF
       END IF
    END IF

    IF (if_induction) THEN
       CALL write_restart_maxwell(comm_one_d, H_mesh, phi_mesh, &
            time, list_mode, Hn, Hn1, Bn, Bn1, &
            phin, phin1, inputs%file_name, it, freq_restart)
    END IF

    IF (if_energy) THEN
       CALL write_restart_temp(comm_one_d_temp, temp_mesh, time, &
            list_mode, tempn, tempn_m1, inputs%file_name, it, freq_restart)
    END IF
    IF (if_concentration) THEN
       CALL write_restart_conc(comm_one_d_conc, conc_mesh, time, &
            list_mode, concn, concn_m1, inputs%file_name, it, freq_restart)
    END IF

  END SUBROUTINE save_run
  !---------------------------------------------------------------------------

  SUBROUTINE INIT
    !==================
    USE my_util
    USE chaine_caractere
    USE periodic
    USE prep_maill
    USE prep_mesh_interface
    USE restart
    USE boundary
    USE sub_plot
    USE def_type_mesh
    USE create_comm
    USE metis_sfemans
    !===JLG July 20, 2019, p3 mesh
    USE mod_gauss_points_2d
    !===JLG July 20, 2019, p3 mesh
    USE st_matrix
    USE st_csr_mhd
    USE matrix_type
    USE symmetric_field
    USE tn_axi
    USE fourier_to_real_for_vtu
    USE subroutine_mass
    USE sft_parallele
    USE update_taylor_navier_stokes
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                         :: vv_mesh_glob, pp_mesh_glob
    TYPE(mesh_type)                         :: H_mesh_glob, phi_mesh_glob, pmag_mesh_glob, temp_mesh_glob
    TYPE(mesh_type)                         :: conc_mesh_glob
    TYPE(mesh_type)                         :: p1_mesh_glob, p2_mesh_glob, p1_c0_mesh_glob, p2_c0_mesh_glob_temp
    TYPE(mesh_type)                         :: p2_c0_mesh_glob_conc
    TYPE(mesh_type)                         :: p3_mesh_glob !===JLG july 20, 2019, p3 mesh
    TYPE(interface_type)                    :: interface_H_phi_glob, interface_H_mu_glob
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_H, list_dom_H_ref
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_temp, list_dom_temp_ref
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom_ns
    INTEGER, DIMENSION(:), ALLOCATABLE      :: list_dom, list_inter, part, list_dummy, list_inter_temp, list_inter_conc
    INTEGER, DIMENSION(:), ALLOCATABLE      :: H_in_to_new, H_in_to_new_ref
    INTEGER, DIMENSION(:), ALLOCATABLE      :: temp_in_to_new, temp_in_to_new_ref
    INTEGER, DIMENSION(:), ALLOCATABLE      :: vv_in_to_new
!!$    CHARACTER(len=200)                      :: data_file
    CHARACTER(len=200)                      :: data_directory
    CHARACTER(len=200)                      :: tit_part, mesh_part_name
    CHARACTER(len=200)                      :: data_fichier
    INTEGER                                 :: nsize
    INTEGER                                 :: k, kp, m, n, i, j
    INTEGER                                 :: code, rank, rank_S, nb_procs, petsc_rank, bloc_size, m_max_pad
    REAL(KIND=8)                            :: time_u, time_h, time_T, error, max_vel_S
    REAL(KIND=8)                            :: time_conc
    LOGICAL                                 :: ns_periodic, mxw_periodic, temp_periodic
    LOGICAL                                 :: conc_periodic
!!$    CHARACTER(len=2)                        :: tit
    PetscMPIInt                             :: nb_procs_F, nb_procs_S

    !===Get numbers of processors===================================================
    CALL MPI_COMM_SIZE(PETSC_COMM_WORLD,nb_procs,code)
    CALL MPI_COMM_RANK(PETSC_COMM_WORLD,petsc_rank,code)

    !===Check if regression=========================================================
    CALL regression_initialize

!!$    !===Decide whether debugging or not=============================================
!!$    CALL sfemansinitialize
!!$    IF (inputs%test_de_convergence) THEN
!!$       IF (inputs%numero_du_test_debug<1 .OR. inputs%numero_du_test_debug>40) THEN
!!$          CALL error_Petsc('BUG in INIT: debug_test_number is not in the correct range')
!!$       END IF
!!$       WRITE(tit,'(i2)') inputs%numero_du_test_debug
!!$       data_file = 'data_'//TRIM(ADJUSTL(tit))
!!$       data_directory = inputs%data_directory_debug
!!$       data_fichier = TRIM(ADJUSTL(data_directory))//'/debug_'//TRIM(ADJUSTL(data_file))
!!$    ELSE
!!$       data_directory = '.'
!!$       data_file='data'
!!$       data_fichier = TRIM(ADJUSTL(data_directory))//'/'//TRIM(ADJUSTL(data_file))
!!$    END IF

    !===Assign correct user functions in boundary module
    call assign_boundary

    !===Read data file==============================================================
    data_fichier = TRIM(ADJUSTL('data'))
    CALL read_my_data(data_fichier)

    !===Debugging===================================================================
    IF (inputs%test_de_convergence) THEN
       inputs%directory = data_directory
    END IF

    !===Initialization for empty vacuum=============================================
    IF (inputs%nb_dom_phi==0) THEN
       inputs%phi_nb_dirichlet_sides = 0
       inputs%nb_inter = 0
       inputs%mu_phi = 1.d0
       inputs%type_fe_phi = -1
       inputs%stab(2) = 0.d0
       IF (ASSOCIATED(inputs%phi_list_dirichlet_sides)) DEALLOCATE(inputs%phi_list_dirichlet_sides)
       ALLOCATE(inputs%phi_list_dirichlet_sides(0))
       IF (ASSOCIATED(inputs%list_inter_H_phi)) DEALLOCATE(inputs%list_inter_H_phi)
       ALLOCATE(inputs%list_inter_H_phi(0))
    END IF

    !===Control inputs==============================================================
    nb_procs_S = inputs%ndim(1)
    nb_procs_F = inputs%ndim(2)
    IF (inputs%ndim(1)*inputs%ndim(2)/=nb_procs) THEN
       CALL error_Petsc('BUG in INIT, nb_proc_space*nb_proc_fourier/=nb_procs')
    END IF

    !===Create communicators========================================================
    CALL create_cart_comm(inputs%ndim,comm_cart,comm_one_d,coord_cart)
    CALL MPI_COMM_SIZE(comm_one_d(2),nb_procs,code)
    CALL MPI_COMM_RANK(comm_one_d(2),rank,code)
    IF (nb_procs_F/=nb_procs) THEN
       CALL error_Petsc('BUG in INIT, nb_procs_F/=nb_procs')
    END IF

    !===Sort out fourier modes======================================================
    m_max_c = inputs%m_max/nb_procs_F
    ALLOCATE(list_mode(m_max_c))
    IF (m_max_c==0) THEN
       CALL error_Petsc('BUG in INIT, m_max_c==0')
    END IF
    IF (inputs%select_mode) THEN
       DO i = 1, m_max_c
          list_mode(i) = inputs%list_mode_lect(i + rank*m_max_c)
       END DO
    ELSE
       DO i = 1, m_max_c
          list_mode(i) = i + rank*m_max_c - 1
       END DO
    END IF

    !===Check periodicity===========================================================
    IF (inputs%my_periodic%nb_periodic_pairs< 1) THEN
       ns_periodic=.FALSE.; mxw_periodic=.FALSE.; temp_periodic = .FALSE.; conc_periodic = .FALSE.
       vvrtz_per%n_bord = 0; vvrt_per%n_bord = 0; vvz_per%n_bord = 0; pp_per%n_bord = 0
       H_phi_per%n_bord = 0 ; temp_per%n_bord = 0 ; conc_per%n_bord = 0
       level_set_per%n_bord = 0
    ELSE
       ns_periodic=.TRUE.; mxw_periodic=.TRUE.; temp_periodic=.TRUE. ; conc_periodic=.TRUE.
    END IF

    !===Creation of logicals for equations==========================================
    if_mass = inputs%if_level_set
    if_momentum = inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
         .OR. inputs%type_pb=='mhs'
    if_induction = inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='mxx' &
         .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs'
    if_energy = inputs%if_temperature
    if_concentration = inputs%if_concentration

    !===JLG july 20, 2019, p3 mesh
    IF (inputs%if_navier_stokes_with_taylor) THEN
       IF (petsc_rank==0) WRITE(*,*) 'INIT: Everything that is not Navier-Stokes is disabled, for Taylor Method'
       if_mass = .FALSE. !===Disable everything that is not NS, for the time being
       if_momentum = inputs%type_pb=='nst'
       if_induction = .FALSE.
       if_energy = .FALSE.
       if_concentration = .FALSE.
    END IF

    !===Check mesh that conc_mesh is a subset of H_mesh===============================
    IF (if_induction) THEN
       ALLOCATE(list_dom_H(inputs%nb_dom_H), H_in_to_new(inputs%nb_dom_H))
       ALLOCATE(list_dom_H_ref(inputs%nb_dom_H), H_in_to_new_ref(inputs%nb_dom_H))
       IF (if_concentration) THEN
          IF (SIZE(list_dom_H) < SIZE(inputs%list_dom_conc)) THEN
             CALL error_Petsc(' BUG: conc must be a subset of Maxwell ')
          END IF
          DO k = 1, inputs%nb_dom_conc
             IF (MINVAL(ABS(inputs%list_dom_H - inputs%list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc(' BUG: conc must be a subset of Maxwell ')
             END IF
             DO kp = 1, inputs%nb_dom_H
                IF (inputs%list_dom_H(kp) == inputs%list_dom_conc(k)) EXIT
             END DO
             H_in_to_new(k) = kp
             list_dom_H(k) = inputs%list_dom_conc(k)
          END DO
          m = inputs%nb_dom_conc
          DO k = 1, inputs%nb_dom_H
             IF (MINVAL(ABS(inputs%list_dom_H(k) - inputs%list_dom_conc)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = k
             list_dom_H(m) = inputs%list_dom_H(k)
          END DO
          IF (m/=inputs%nb_dom_H) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_H ')
          END IF
       ELSE
          DO k = 1, inputs%nb_dom_H
             H_in_to_new(k) = k
          END DO
          list_dom_H = inputs%list_dom_H
       END IF
    ELSE
       ALLOCATE(H_in_to_new(0))
       ALLOCATE(H_in_to_new_ref(0))
    END IF

    !===Check mesh that conc_mesh is a subset of temp_mesh===============================
    IF (if_energy) THEN
       ALLOCATE(list_dom_temp(inputs%nb_dom_temp), temp_in_to_new(inputs%nb_dom_temp))
       ALLOCATE(list_dom_temp_ref(inputs%nb_dom_temp), temp_in_to_new_ref(inputs%nb_dom_temp))
       IF (if_concentration) THEN
          IF (SIZE(list_dom_temp) < SIZE(inputs%list_dom_conc)) THEN
             CALL error_Petsc(' BUG: conc must be a subset of temp ')
          END IF
          DO k = 1, inputs%nb_dom_conc
             IF (MINVAL(ABS(inputs%list_dom_temp - inputs%list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc(' BUG: conc must be a subset of temp ')
             END IF
             DO kp = 1, inputs%nb_dom_temp
                IF (inputs%list_dom_temp(kp) == inputs%list_dom_conc(k)) EXIT
             END DO
             temp_in_to_new(k) = kp
             list_dom_temp(k) = inputs%list_dom_conc(k)
          END DO
          m = inputs%nb_dom_conc
          DO k = 1, inputs%nb_dom_temp
             IF (MINVAL(ABS(inputs%list_dom_temp(k) - inputs%list_dom_conc)) == 0) CYCLE
             m = m + 1
             temp_in_to_new(m) = k
             list_dom_temp(m) = inputs%list_dom_temp(k)
          END DO
          IF (m/=inputs%nb_dom_temp) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_temp ')
          END IF
       ELSE
          DO k = 1, inputs%nb_dom_temp
             temp_in_to_new(k) = k
          END DO
          list_dom_temp = inputs%list_dom_temp
       END IF
    ELSE
       ALLOCATE(temp_in_to_new(0))
       ALLOCATE(temp_in_to_new_ref(0))
    END IF

    !===Check mesh that conc_mesh is a subset of vv_mesh===============================
    IF (if_momentum) THEN
       ALLOCATE(list_dom_ns(inputs%nb_dom_ns), vv_in_to_new(inputs%nb_dom_ns))
       IF (if_concentration) THEN
          IF (SIZE(list_dom_ns) < SIZE(inputs%list_dom_conc)) THEN
             CALL error_Petsc(' BUG: conc must be a subset of NS ')
          END IF
          DO k = 1, inputs%nb_dom_conc
             IF (MINVAL(ABS(inputs%list_dom_ns - inputs%list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc(' BUG: conc must be a subset of NS ')
             END IF
             DO kp = 1, inputs%nb_dom_ns
                IF (inputs%list_dom_ns(kp) == inputs%list_dom_conc(k)) EXIT
             END DO
             vv_in_to_new(k) = kp
             list_dom_ns(k) = inputs%list_dom_conc(k)
          END DO
          m = inputs%nb_dom_conc
          DO k = 1, inputs%nb_dom_ns
             IF (MINVAL(ABS(inputs%list_dom_ns(k) - inputs%list_dom_conc)) == 0) CYCLE
             m = m + 1
             vv_in_to_new(m) = k
             list_dom_ns(m) = inputs%list_dom_ns(k)
          END DO
          IF (m/=inputs%nb_dom_ns) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_ns ')
          END IF
       ELSE
          DO k = 1, inputs%nb_dom_ns
             vv_in_to_new(k) = k
          END DO
          list_dom_ns = inputs%list_dom_ns
       END IF
    ELSE
       ALLOCATE(vv_in_to_new(0))
    END IF

    !===Check mesh that vv_mesh is a subset of H_mesh===============================
    IF (if_induction) THEN
       H_in_to_new_ref=H_in_to_new
       list_dom_H_ref=list_dom_H
       IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_ns)) THEN
             CALL error_Petsc(' BUG: NS must be a subset of Maxwell ')
          END IF
          DO k = 1+inputs%nb_dom_conc, inputs%nb_dom_ns
             IF (MINVAL(ABS(list_dom_H - list_dom_ns(k))) /= 0) THEN
                CALL error_Petsc(' BUG: NS must be a subset of Maxwell ')
             END IF
             DO kp = 1+inputs%nb_dom_conc, inputs%nb_dom_H
                IF (list_dom_H_ref(kp) == list_dom_ns(k)) EXIT
             END DO
             H_in_to_new(k) = H_in_to_new_ref(kp)
             list_dom_H(k) = list_dom_ns(k)
          END DO
          m = inputs%nb_dom_ns
          DO k = 1+inputs%nb_dom_conc, inputs%nb_dom_H
             IF (MINVAL(ABS(list_dom_H_ref(k) - list_dom_ns)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = H_in_to_new_ref(k)
             list_dom_H(m) = list_dom_H_ref(k)
          END DO
          IF (m/=inputs%nb_dom_H) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_H ')
          END IF
       END IF
    END IF

    !===Check mesh that vv_mesh is a subset of temp_mesh============================
    IF (if_energy) THEN
       temp_in_to_new_ref=temp_in_to_new
       list_dom_temp_ref=list_dom_temp
       IF (if_momentum) THEN
          IF (SIZE(list_dom_temp) < SIZE(list_dom_ns)) THEN
             CALL error_Petsc(' BUG: NS must be a subset of temp ')
          END IF
          DO k = 1+inputs%nb_dom_conc, inputs%nb_dom_ns
             IF (MINVAL(ABS(list_dom_temp - list_dom_ns(k))) /= 0) THEN
                CALL error_Petsc(' BUG: NS must be a subset of temp ')
             END IF
             DO kp = 1+inputs%nb_dom_conc, inputs%nb_dom_temp
                IF (list_dom_temp_ref(kp) == list_dom_ns(k)) EXIT
             END DO
             temp_in_to_new(k) = temp_in_to_new_ref(kp)
             list_dom_temp(k) = list_dom_ns(k)
          END DO
          m = inputs%nb_dom_ns
          DO k = 1+inputs%nb_dom_conc, inputs%nb_dom_temp
             IF (MINVAL(ABS(list_dom_temp_ref(k) - list_dom_ns)) == 0) CYCLE
             m = m + 1
             temp_in_to_new(m) = temp_in_to_new_ref(k)
             list_dom_temp(m) = list_dom_temp_ref(k)
          END DO
          IF (m/=inputs%nb_dom_temp) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_temp ')
          END IF
       END IF
    END IF

    !===Check mesh that temp_mesh is a subset of H_mesh=============================
    IF (if_induction) THEN
       H_in_to_new_ref=H_in_to_new
       list_dom_H_ref=list_dom_H
       IF (if_energy) THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_temp)) THEN
             CALL error_Petsc(' BUG: temp must be a subset of H ')
          END IF
          DO k = 1+inputs%nb_dom_ns, inputs%nb_dom_temp
             IF (MINVAL(ABS(list_dom_H - list_dom_temp(k))) /= 0) THEN
                CALL error_Petsc(' BUG: temp must be a subset of H ')
             END IF
             DO kp = 1+inputs%nb_dom_ns, inputs%nb_dom_H
                IF (list_dom_H_ref(kp) == list_dom_temp(k)) EXIT
             END DO
             H_in_to_new(k) = H_in_to_new_ref(kp)
             list_dom_H(k) = list_dom_temp(k)
          END DO
          m = inputs%nb_dom_temp
          DO k = 1+inputs%nb_dom_ns, inputs%nb_dom_H
             IF (MINVAL(ABS(list_dom_H_ref(k) - list_dom_temp)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = H_in_to_new_ref(k)
             list_dom_H(m) = list_dom_H_ref(k)
          END DO
          IF (m/=inputs%nb_dom_H) THEN
             CALL error_Petsc(' BUG: m/=inputs%nb_dom_H ')
          END IF
       END IF
    END IF

    !===Create interfaces in meshes=================================================
    IF (if_momentum .AND. (.NOT. if_induction)) THEN
       IF (if_energy) THEN
          nsize = SIZE(list_dom_temp)
          ALLOCATE(list_dom(nsize))
          list_dom = list_dom_temp
          ALLOCATE(list_inter(SIZE(inputs%list_inter_v_T)))
          list_inter = inputs%list_inter_v_T
       ELSE
          nsize = SIZE(list_dom_ns)
          ALLOCATE(list_dom(nsize))
          list_dom = list_dom_ns
          IF (if_concentration) THEN
             ALLOCATE(list_inter(SIZE(inputs%list_inter_c_v)))
             list_inter = inputs%list_inter_c_v
          ELSE
             ALLOCATE(list_inter(0))
          END IF
       END IF
    ELSE
       nsize = SIZE(list_dom_H)+SIZE(inputs%list_dom_phi)
       ALLOCATE(list_dom(nsize))
       list_dom(1:SIZE(list_dom_H)) = list_dom_H
       IF (SIZE(inputs%list_dom_phi)>0) THEN
          list_dom(SIZE(list_dom_H)+1:) = inputs%list_dom_phi
       END IF
       nsize = SIZE(inputs%list_inter_mu)+SIZE(inputs%list_inter_H_phi)
       ALLOCATE(list_inter(nsize))
       IF (SIZE(inputs%list_inter_mu)>0) THEN
          list_inter(1:SIZE(inputs%list_inter_mu)) = inputs%list_inter_mu
       END IF
       IF (SIZE(inputs%list_inter_H_phi)>0) THEN
          list_inter(SIZE(inputs%list_inter_mu)+1:) = inputs%list_inter_H_phi
       END IF
    END IF
    IF (if_energy) THEN
       ALLOCATE(list_inter_temp(0))
    END IF
    IF (if_concentration) THEN
       ALLOCATE(list_inter_conc(0))
    END IF

    !===Create meshes===============================================================
    CALL load_dg_mesh_free_format(inputs%directory, inputs%file_name, list_dom, &
         list_inter, 1, p1_mesh_glob, inputs%iformatted)
    CALL load_dg_mesh_free_format(inputs%directory, inputs%file_name, list_dom, &
         list_inter, 2, p2_mesh_glob, inputs%iformatted)

    !===JLG july 20, 2019, p3 mesh
    IF (inputs%type_fe_velocity==3) THEN
       CALL create_p3_mesh(p1_mesh_glob, p2_mesh_glob, p3_mesh_glob, 3)
    END IF
    !===JLG july 20, 2019, p3 mesh
    IF (if_concentration) THEN
       !===MODIFICATION: Dirichlet nodes in temp_mesh not created if list_dom > list_dom_conc
       CALL load_dg_mesh_free_format(inputs%directory, inputs%file_name, inputs%list_dom_conc, &
            list_inter_conc, 2, p2_c0_mesh_glob_conc, inputs%iformatted)
    END IF
    IF (if_energy) THEN
       !===MODIFICATION: Dirichlet nodes in temp_mesh not created if list_dom > list_dom_temp
       CALL load_dg_mesh_free_format(inputs%directory, inputs%file_name, list_dom_temp, &
            list_inter_temp, 2, p2_c0_mesh_glob_temp, inputs%iformatted)
    END IF
    IF (if_induction) THEN
       ALLOCATE(list_dummy(0))
       CALL load_dg_mesh_free_format(inputs%directory, inputs%file_name, list_dom, &
            inputs%list_inter_H_phi, 1, p1_c0_mesh_glob, inputs%iformatted)
    END IF

    !===Start Metis mesh generation=================================================
    ALLOCATE(part(p1_mesh_glob%me))
    WRITE(tit_part,'(i4)') inputs%ndim(1)
    mesh_part_name='mesh_part_S'//TRIM(ADJUSTL(tit_part))//'.'//TRIM(ADJUSTL(inputs%file_name))
    IF (inputs%if_read_partition) THEN
       OPEN(UNIT=51, FILE=mesh_part_name, STATUS='unknown', FORM='formatted')
       READ(51,*) part
       CLOSE(51)
       WRITE(*,*) 'read partition'
    ELSE
       WRITE(*,*) 'create partition'
       CALL part_mesh_M_T_H_phi(nb_procs_S, inputs%list_dom_conc, inputs%list_dom_ns, inputs%list_dom_temp, &
            inputs%list_dom_H, inputs%list_dom_phi, p1_mesh_glob, list_inter, part, inputs%my_periodic)
       IF (petsc_rank==0) THEN
          OPEN(UNIT=51, FILE=mesh_part_name, STATUS='replace', FORM='formatted')
          WRITE(51,*) part
          CLOSE(51)
       END IF
    END IF

    !===Extract local meshes from global meshes=====================================
    !===Specific to momentum (velocity)
    IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
       IF (inputs%type_fe_velocity==2) THEN !===JLG july 20, 2019, p3 mesh
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p1_mesh_glob,part,list_dom_ns,pp_mesh_glob,pp_mesh)
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_mesh_glob,part,list_dom_ns,vv_mesh_glob,vv_mesh)
       ELSE IF (inputs%type_fe_velocity==3) THEN
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_mesh_glob,part,list_dom_ns,pp_mesh_glob,pp_mesh)
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p3_mesh_glob,part,list_dom_ns,vv_mesh_glob,vv_mesh)
       ELSE
          CALL error_PETSC('Bug in INIT, inputs%type_fe_velocity not correct')
       END IF
       !===JLG july 20, 2019, p3 mesh
       !===Use increasing vertex index enumeration
       !CALL incr_vrtx_indx_enumeration(vv_mesh,inputs%type_fe_velocity)
       !CALL incr_vrtx_indx_enumeration(pp_mesh,inputs%type_fe_velocity-1)
       !===JLG july 20, 2019, p3 mesh
       ALLOCATE(comm_one_d_ns(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_ns(2), code)
       !Correction FL-ID, 4/4/13
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (pp_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_ns(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_ns(1),code)
       END IF

       !===Test whether pp_mesh and vv_mesh meshes coincide=========================
       DO m = 1, vv_mesh%me
          DO n = 1, 3
             IF (MINVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(n,m)) &
                  -pp_mesh%rr(1,pp_mesh%jj(:,m))))>1.d-16 &
                  .OR. MINVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(n,m))&
                  -pp_mesh%rr(2,pp_mesh%jj(:,m))))>1.d-16) THEN
                CALL error_Petsc('BUG in INIT, vv and pp global meshes are different')
             END IF
          END DO
       END DO

       !===Create periodic structures===============================================
       IF (ns_periodic) THEN
          CALL prep_periodic(inputs%my_periodic, pp_mesh, pp_per)
          CALL prep_periodic(inputs%my_periodic, vv_mesh, vvz_per)
          CALL prep_periodic_bloc(inputs%my_periodic, vv_mesh, vvrt_per, 2)
          CALL prep_periodic_bloc(inputs%my_periodic, vv_mesh, vvrtz_per, 3)
          IF (inputs%if_level_set_P2) THEN
             CALL prep_periodic(inputs%my_periodic, vv_mesh, level_set_per)
          ELSE
             CALL prep_periodic(inputs%my_periodic, pp_mesh, level_set_per)
          END IF
       END IF

       !===Create global csr structure==============================================
       IF (pp_mesh%me/=0) THEN
          CALL st_aij_csr_glob_block(comm_one_d_ns(1),1,vv_mesh_glob,vv_mesh,vv_1_LA, opt_per=vvz_per)
          CALL st_aij_csr_glob_block(comm_one_d_ns(1),3,vv_mesh_glob,vv_mesh,vv_3_LA, opt_per=vvrtz_per)
          CALL st_aij_csr_glob_block(comm_one_d_ns(1),1,pp_mesh_glob,pp_mesh,pp_1_LA, opt_per=pp_per)
          !===Prepare csr structure for post processing curl_u==========================
          CALL st_aij_csr_glob_block(comm_one_d_ns(1),1,vv_mesh_glob,vv_mesh,vizu_rot_u_LA)
       END IF

       !===Create symmetric points==================================================
       IF (inputs%is_mesh_symmetric) THEN
          ALLOCATE(vv_mz_LA(vv_mesh%np))
          CALL symmetric_points(vv_mesh, vv_mesh_glob, vv_mz_LA)
       END IF

       !===Deallocate global meshes=================================================
       CALL free_mesh(vv_mesh_glob)
       CALL free_mesh(pp_mesh_glob)

       !===Start Gauss points generation============================================
       !===JLG july 20, 2019, p3 mesh
       vv_mesh%edge_stab=.FALSE.
       pp_mesh%edge_stab=.FALSE.
       CALL gauss_points_2d(vv_mesh,inputs%type_fe_velocity)
       CALL gauss_points_2d(pp_mesh,inputs%type_fe_velocity-1)
       !===JLG july 20, 2019, p3 mesh
    END IF !=== (if_momentum .OR. inputs%type_pb=='mxx')

    !===Extract local meshes from global meshes for Maxwell=========================
    IF (if_induction) THEN
       CALL  extract_mesh(comm_one_d(1),nb_procs_S,p1_c0_mesh_glob,part,list_dom_H,pmag_mesh_glob,pmag_mesh)
       IF (inputs%type_fe_H==1) THEN
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p1_mesh_glob,part,list_dom_H,H_mesh_glob,H_mesh)
       ELSE
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_mesh_glob,part,list_dom_H,H_mesh_glob,H_mesh)
       END IF
       IF (inputs%type_fe_phi==1) THEN
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p1_mesh_glob,part,inputs%list_dom_phi,phi_mesh_glob,phi_mesh)
       ELSE
          CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_mesh_glob,part,inputs%list_dom_phi,phi_mesh_glob,phi_mesh)
       END IF
       !===JLG july 20, 2019, p3 mesh
       !===Use increasing vertex index enumeration
       !CALL incr_vrtx_indx_enumeration(pmag_mesh,1)
       !CALL incr_vrtx_indx_enumeration(H_mesh,inputs%type_fe_H)
       !CALL incr_vrtx_indx_enumeration(phi_mesh,inputs%type_fe_phi)
       !===JLG july 20, 2019, p3 mesh
    END IF

    !===Extract local meshes from global meshes for temperature=====================
    IF (if_energy) THEN
       CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_c0_mesh_glob_temp,part,list_dom_temp,temp_mesh_glob,temp_mesh)
       ALLOCATE(comm_one_d_temp(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_temp(2), code)
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (temp_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_temp(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_temp(1),code)
       END IF
    END IF

    !===Extract local meshes from global meshes for concentration=====================
    IF (if_concentration) THEN
       CALL extract_mesh(comm_one_d(1),nb_procs_S,p2_c0_mesh_glob_conc,part,inputs%list_dom_conc,conc_mesh_glob,conc_mesh)
       ALLOCATE(comm_one_d_conc(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_conc(2), code)
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (conc_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_conc(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_conc(1),code)
       END IF
    END IF

    !===Deallocate global meshes====================================================
    CALL free_mesh(p1_mesh_glob)
    CALL free_mesh(p2_mesh_glob)
    !===JLG July 20, 2019, p3 mesh
    IF (inputs%type_fe_velocity==3) THEN
       CALL free_mesh(p3_mesh_glob)
    END IF
    !===JLG July 20, 2019, p3 mesh
    IF (if_induction) THEN
       DEALLOCATE(list_dummy)
       CALL free_mesh(p1_c0_mesh_glob)
    END IF
    IF (if_energy) THEN
       DEALLOCATE(list_inter_temp)
       CALL free_mesh(p2_c0_mesh_glob_temp)
    END IF
    IF (if_concentration) THEN
       DEALLOCATE(list_inter_conc)
       CALL free_mesh(p2_c0_mesh_glob_conc)
    END IF

    !===Specific to induction equation==============================================
    IF (if_induction) THEN

       !===Verify that pmag_mesh and H_mesh coincide================================
       IF (pmag_mesh%me/=0) THEN
          error = 0.d0
          DO k = 1, 2
             DO n = 1, SIZE(pmag_mesh%jj,1)
                error = error + MAXVAL(ABS(pmag_mesh%rr(k,pmag_mesh%jj(n,:))-H_mesh%rr(k,H_mesh%jj(n,1:pmag_mesh%me))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
             CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14')
          END IF
       END IF

       !===Verify if H and temp meshes coincide on the temp domain=======================
       IF (if_energy) THEN
          IF (temp_mesh%me/=0) THEN
             error = 0.d0
             DO k = 1, 2
                DO n = 1, SIZE(H_mesh%jj,1)
                   error = error + MAXVAL(ABS(temp_mesh%rr(k,temp_mesh%jj(n,:))-H_mesh%rr(k,H_mesh%jj(n,1:temp_mesh%me))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14')
             END IF

             error = error + MAXVAL(ABS(temp_mesh%rr(1,temp_mesh%jj(4,1:temp_mesh%me)) &
                  -(H_mesh%rr(1,H_mesh%jj(2,1:temp_mesh%me))+H_mesh%rr(1,H_mesh%jj(3,1:temp_mesh%me)))/2))&
                  + MAXVAL(ABS(temp_mesh%rr(1,temp_mesh%jj(5,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(3,1:temp_mesh%me))+H_mesh%rr(1,H_mesh%jj(1,1:temp_mesh%me)))/2))&
                  + MAXVAL(ABS(temp_mesh%rr(1,temp_mesh%jj(6,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(1,1:temp_mesh%me))+H_mesh%rr(1,H_mesh%jj(2,1:temp_mesh%me)))/2))&
                  + MAXVAL(ABS(temp_mesh%rr(2,temp_mesh%jj(4,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(2,1:temp_mesh%me))+H_mesh%rr(2,H_mesh%jj(3,1:temp_mesh%me)))/2))&
                  + MAXVAL(ABS(temp_mesh%rr(2,temp_mesh%jj(5,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(3,1:temp_mesh%me))+H_mesh%rr(2,H_mesh%jj(1,1:temp_mesh%me)))/2))&
                  + MAXVAL(ABS(temp_mesh%rr(2,temp_mesh%jj(6,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(1,1:temp_mesh%me))+H_mesh%rr(2,H_mesh%jj(2,1:temp_mesh%me)))/2))
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                WRITE(*,*) ' WARNING: temp_mesh and H_mesh do not coincide on the temp domain.'
                WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                     error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:)))
             END IF

             error=0.d0
             DO k = 1, temp_mesh%me
                DO n = 1, 2
                   error = error+ MAXVAL(ABS(temp_mesh%rr(n,temp_mesh%jj(1:3,k))-H_mesh%rr(n,H_mesh%jj(1:3,k))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('temp_mesh and H_mesh do NOT have the same P1 nodes')
             END IF

          END IF
       END IF

       !===Verify if H and NS meshes coincide on the NS domain=======================
       IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
          IF (vv_mesh%me/=0) THEN
             error = 0.d0
             DO k = 1, 2
                DO n = 1, SIZE(H_mesh%jj,1)
                   error = error + MAXVAL(ABS(vv_mesh%rr(k,vv_mesh%jj(n,:))-H_mesh%rr(k,H_mesh%jj(n,1:vv_mesh%me))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14')
             END IF

             error = error + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(4,1:vv_mesh%me)) &
                  -(H_mesh%rr(1,H_mesh%jj(2,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(3,1:vv_mesh%me)))/2))&
                  + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(5,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(3,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(1,1:vv_mesh%me)))/2))&
                  + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(6,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(1,1:vv_mesh%me))+H_mesh%rr(1,H_mesh%jj(2,1:vv_mesh%me)))/2))&
                  + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(4,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(2,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(3,1:vv_mesh%me)))/2))&
                  + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(5,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(3,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(1,1:vv_mesh%me)))/2))&
                  + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(6,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(1,1:vv_mesh%me))+H_mesh%rr(2,H_mesh%jj(2,1:vv_mesh%me)))/2))
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                WRITE(*,*) ' WARNING: vv_mesh and H_mesh do not coincide on the NS domain.'
                WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                     error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:)))
             END IF

             error=0.d0
             DO k = 1, vv_mesh%me
                DO n = 1, 2
                   error = error+ MAXVAL(ABS(vv_mesh%rr(n,vv_mesh%jj(1:3,k))-H_mesh%rr(n,H_mesh%jj(1:3,k))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('vv_mesh and H_mesh do NOT have the same P1 nodes')
             END IF

          END IF
       END IF

       !===Verify if H and conc meshes coincide on the conc domain=======================
       IF (if_concentration) THEN
          IF (conc_mesh%me/=0) THEN
             error = 0.d0
             DO k = 1, 2
                DO n = 1, SIZE(H_mesh%jj,1)
                   error = error + MAXVAL(ABS(conc_mesh%rr(k,conc_mesh%jj(n,:))-H_mesh%rr(k,H_mesh%jj(n,1:conc_mesh%me))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14')
             END IF

             error = error + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(4,1:conc_mesh%me)) &
                  -(H_mesh%rr(1,H_mesh%jj(2,1:conc_mesh%me))+H_mesh%rr(1,H_mesh%jj(3,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(5,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(3,1:conc_mesh%me))+H_mesh%rr(1,H_mesh%jj(1,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(6,:)) &
                  -(H_mesh%rr(1,H_mesh%jj(1,1:conc_mesh%me))+H_mesh%rr(1,H_mesh%jj(2,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(4,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(2,1:conc_mesh%me))+H_mesh%rr(2,H_mesh%jj(3,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(5,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(3,1:conc_mesh%me))+H_mesh%rr(2,H_mesh%jj(1,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(6,:)) &
                  -(H_mesh%rr(2,H_mesh%jj(1,1:conc_mesh%me))+H_mesh%rr(2,H_mesh%jj(2,1:conc_mesh%me)))/2))
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                WRITE(*,*) ' WARNING: conc_mesh and H_mesh do not coincide on the conc domain.'
                WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                     error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:)))
             END IF

             error=0.d0
             DO k = 1, conc_mesh%me
                DO n = 1, 2
                   error = error+ MAXVAL(ABS(conc_mesh%rr(n,conc_mesh%jj(1:3,k))-H_mesh%rr(n,H_mesh%jj(1:3,k))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(H_mesh%rr(1,1) -H_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('conc_mesh and H_mesh do NOT have the same P1 nodes')
             END IF

          END IF
       END IF

       !===Create interface structures==============================================
       CALL load_interface(H_mesh_glob, H_mesh_glob, inputs%list_inter_mu, interface_H_mu_glob, .FALSE.)
       CALL load_interface(H_mesh_glob, phi_mesh_glob, inputs%list_inter_H_phi, interface_H_phi_glob, .TRUE.)

       IF (H_mesh%me /=0) THEN
          CALL load_interface(H_mesh, H_mesh, inputs%list_inter_mu, interface_H_mu, .FALSE.)
       ELSE
          interface_H_mu%mes = 0
       END IF

       IF (H_mesh%me * phi_mesh%me /=0) THEN
          CALL load_interface(H_mesh, phi_mesh, inputs%list_inter_H_phi, interface_H_phi, .TRUE.)
       ELSE
          interface_H_phi%mes = 0
       END IF
       !===JLG july 20, 2019, p3 mesh
       !===Use increasing vertex index enumeration
       !CALL incr_vrtx_indx_enumeration_for_interfaces(interface_H_phi,inputs%type_fe_H+1,inputs%type_fe_phi+1)
       !===JLG july 20, 2019, p3 mesh

       !===Create periodic structures for Maxwell===================================
       IF (mxw_periodic) THEN
          CALL prep_periodic_H_p_phi_bc(inputs%my_periodic, H_mesh, pmag_mesh, phi_mesh, H_phi_per)
          WRITE(*,*) 'periodic MHD done'
       END IF

       !===Create global csr structure==============================================
       CALL st_scr_maxwell_mu_H_p_phi(comm_one_d(1), H_mesh_glob, H_mesh, pmag_mesh_glob, pmag_mesh, &
            phi_mesh_glob, phi_mesh, interface_H_phi_glob, interface_H_mu_glob, &
            LA_H, LA_pmag, LA_phi, LA_mhd, opt_per=H_phi_per)

       !===Prepare csr structure for post processing grad phi=======================+++
       CALL st_aij_csr_glob_block(comm_one_d(1),1,phi_mesh_glob,phi_mesh, vizu_grad_phi_LA)

       !===Prepare csr structure for post processing rot !h==========================+++
       CALL st_aij_csr_glob_block(comm_one_d(1),1,H_mesh_glob,H_mesh, vizu_rot_h_LA)

       IF (inputs%is_mesh_symmetric) THEN
          ALLOCATE(H_mz_LA(H_mesh%np))
          CALL symmetric_points(H_mesh, H_mesh_glob, H_mz_LA)
       END IF

       !===Deallocate global meshes=================================================
       CALL free_mesh(H_mesh_glob)
       CALL free_mesh(pmag_mesh_glob)
       CALL free_mesh(phi_mesh_glob)
       CALL free_interface(interface_H_mu_glob)
       CALL free_interface(interface_H_phi_glob)

       !===Start Gauss points generation============================================
       !===JLG july 20, 2019, p3 mesh
       H_mesh%edge_stab    = .FALSE.
       pmag_mesh%edge_stab = .FALSE.
       phi_mesh%edge_stab  = .FALSE.
       IF (1.LE.inputs%type_fe_H .AND. inputs%type_fe_H.LE.2) THEN
          CALL gauss_points_2d(H_mesh,inputs%type_fe_H)
       END IF
       CALL gauss_points_2d(pmag_mesh,1)
       IF (1.LE.inputs%type_fe_phi .AND. inputs%type_fe_phi.LE.2) THEN
          CALL gauss_points_2d(phi_mesh,inputs%type_fe_phi)
       END IF
       !===JLG july 20, 2019, p3 mesh

       !===Create sigma_field=======================================================
       ALLOCATE(sigma_field(H_mesh%me)) !===sigma field is P0 and defined on H_mesh
       DO m = 1, H_mesh%me
          DO k=1, inputs%nb_dom_H
             IF (H_mesh%i_d(m) == list_dom_H(k)) THEN
                sigma_field(m) = inputs%sigma(H_in_to_new(k))
             ENDIF
          ENDDO
       END DO

       !===Create mu_H_field========================================================
       ALLOCATE(mu_H_field(H_mesh%np)) !===mu_H_field is defined at nodes of H_mesh
       IF (inputs%analytical_permeability) THEN !===JLG + DCQ, June 26 2013
          mu_H_field = mu_bar_in_fourier_space(H_mesh,1,H_mesh%np)
       ELSE
          DO m = 1, H_mesh%me
             DO k=1, inputs%nb_dom_H
                IF (H_mesh%i_d(m) == list_dom_H(k)) THEN
                   mu_H_field(H_mesh%jj(:,m)) = inputs%mu_H(H_in_to_new(k))
                ENDIF
             ENDDO
          END DO
       END IF
       !===Create mu_H_field========================================================
       !===Artificial boundary condition on phi on sphere of radius R
       !===d(phi)/dR + (1/R)*phi = 0. Assumes that phi=C/r at infinity
       !===Feature is currently disabled.
       R_fourier=-1.d0 !Negative radius disables the boundary condition
       index_fourier=0 !Index of spherical boundary where Fourier BC enforced

    END IF

    !===Specific to temperature=====================================================
    IF (if_energy) THEN
       !===Verify if temp and NS coincide on the NS domain=======================
       IF (vv_mesh%me/=0) THEN
          error = 0.d0
          DO k = 1, 2
             DO n = 1, SIZE(temp_mesh%jj,1)
                error = error + MAXVAL(ABS(vv_mesh%rr(k,vv_mesh%jj(n,:))-temp_mesh%rr(k,temp_mesh%jj(n,1:vv_mesh%me))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
             CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14')
          END IF

          error = error + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(4,1:vv_mesh%me)) &
               -(temp_mesh%rr(1,temp_mesh%jj(2,1:vv_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(3,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(5,:)) &
               -(temp_mesh%rr(1,temp_mesh%jj(3,1:vv_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(1,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(1,vv_mesh%jj(6,:)) &
               -(temp_mesh%rr(1,temp_mesh%jj(1,1:vv_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(2,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(4,:)) &
               -(temp_mesh%rr(2,temp_mesh%jj(2,1:vv_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(3,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(5,:)) &
               -(temp_mesh%rr(2,temp_mesh%jj(3,1:vv_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(1,1:vv_mesh%me)))/2))&
               + MAXVAL(ABS(vv_mesh%rr(2,vv_mesh%jj(6,:)) &
               -(temp_mesh%rr(2,temp_mesh%jj(1,1:vv_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(2,1:vv_mesh%me)))/2))
          IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
             WRITE(*,*) ' WARNING: vv_mesh and temp_mesh do not coincide on the NS domain.'
             WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                  error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:)))
          END IF

          error=0.d0
          DO k = 1, vv_mesh%me
             DO n = 1, 2
                error = error+ MAXVAL(ABS(vv_mesh%rr(n,vv_mesh%jj(1:3,k))-temp_mesh%rr(n,temp_mesh%jj(1:3,k))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
             CALL error_Petsc('vv_mesh and temp_mesh do NOT have the same P1 nodes')
          END IF

       END IF

       !===Verify if temp and conc coincide on the conc domain=======================
       IF (if_concentration) THEN
          IF (conc_mesh%me/=0) THEN
             error = 0.d0
             DO k = 1, 2
                DO n = 1, SIZE(temp_mesh%jj,1)
                   error = error + MAXVAL(ABS(conc_mesh%rr(k,conc_mesh%jj(n,:))-temp_mesh%rr(k,temp_mesh%jj(n,1:conc_mesh%me))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14')
             END IF

             error = error + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(4,1:conc_mesh%me)) &
                  -(temp_mesh%rr(1,temp_mesh%jj(2,1:conc_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(3,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(5,:)) &
                  -(temp_mesh%rr(1,temp_mesh%jj(3,1:conc_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(1,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(6,:)) &
                  -(temp_mesh%rr(1,temp_mesh%jj(1,1:conc_mesh%me))+temp_mesh%rr(1,temp_mesh%jj(2,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(4,:)) &
                  -(temp_mesh%rr(2,temp_mesh%jj(2,1:conc_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(3,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(5,:)) &
                  -(temp_mesh%rr(2,temp_mesh%jj(3,1:conc_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(1,1:conc_mesh%me)))/2))&
                  + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(6,:)) &
                  -(temp_mesh%rr(2,temp_mesh%jj(1,1:conc_mesh%me))+temp_mesh%rr(2,temp_mesh%jj(2,1:conc_mesh%me)))/2))
             IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
                WRITE(*,*) ' WARNING: conc_mesh and temp_mesh do not coincide on the conc domain.'
                WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                     error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:)))
             END IF

             error=0.d0
             DO k = 1, conc_mesh%me
                DO n = 1, 2
                   error = error+ MAXVAL(ABS(conc_mesh%rr(n,conc_mesh%jj(1:3,k))-temp_mesh%rr(n,temp_mesh%jj(1:3,k))))
                END DO
             END DO
             IF (error/MAXVAL(ABS(temp_mesh%rr(1,1) -temp_mesh%rr(1,:))) .GE. 5.d-14) THEN
                CALL error_Petsc('conc_mesh and temp_mesh do NOT have the same P1 nodes')
             END IF

          END IF
       END IF
       !===Create periodic structures for temperature===============================
       IF (temp_periodic) THEN
          CALL prep_periodic(inputs%my_periodic, temp_mesh, temp_per)
       END IF

       !===Create csr structure for temperature=====================================
       CALL st_aij_csr_glob_block(comm_one_d_temp(1),1,temp_mesh_glob,temp_mesh,temp_1_LA, opt_per=temp_per)

       !===Deallocate global meshes=================================================
       CALL free_mesh(temp_mesh_glob)

       !===Start Gauss points generation============================================
       !===JLG July 20, 2019, p3 mesh
       temp_mesh%edge_stab    = .FALSE.
       CALL gauss_points_2d(temp_mesh,2)
       !===JLG July 20, 2019, p3 mesh


       !===Create temperature_diffusivity_field=====================================
       ALLOCATE(temperature_diffusivity_field(temp_mesh%me))
       DO m = 1, temp_mesh%me
          DO k=1, inputs%nb_dom_temp
             IF (temp_mesh%i_d(m) == list_dom_temp(k)) THEN
                temperature_diffusivity_field(m) = inputs%temperature_diffusivity(temp_in_to_new(k))
             END IF
          END DO
       END DO

       !===Create vol_heat_capacity_field===========================================
       ALLOCATE(vol_heat_capacity_field(temp_mesh%me))
       DO m = 1, temp_mesh%me
          DO k=1, inputs%nb_dom_temp
             IF (temp_mesh%i_d(m) == list_dom_temp(k)) THEN
                vol_heat_capacity_field(m) = inputs%vol_heat_capacity(temp_in_to_new(k))
             END IF
          END DO
       END DO
    END IF


    !===Specific to concentration=====================================================
    IF (if_concentration) THEN
       !===Verify if vv_mesh and conc_mesh coincide on the conc domain=======================
       IF (vv_mesh%me/=0) THEN
          error = 0.d0
          DO k = 1, 2
             DO n = 1, SIZE(vv_mesh%jj,1)
                error = error + MAXVAL(ABS(conc_mesh%rr(k,conc_mesh%jj(n,:))-vv_mesh%rr(k,vv_mesh%jj(n,1:conc_mesh%me))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(vv_mesh%rr(1,1) -vv_mesh%rr(1,:))) .GE. 5.d-14) THEN
             CALL error_Petsc('BUG in INIT, (error/MAXVAL(ABS(conc_mesh%rr(1,1) -conc_mesh%rr(1,:))) .GE. 5.d-14')
          END IF

          error = error + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(4,1:conc_mesh%me)) &
               -(vv_mesh%rr(1,vv_mesh%jj(2,1:conc_mesh%me))+vv_mesh%rr(1,vv_mesh%jj(3,1:conc_mesh%me)))/2))&
               + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(5,:)) &
               -(vv_mesh%rr(1,vv_mesh%jj(3,1:conc_mesh%me))+vv_mesh%rr(1,vv_mesh%jj(1,1:conc_mesh%me)))/2))&
               + MAXVAL(ABS(conc_mesh%rr(1,conc_mesh%jj(6,:)) &
               -(vv_mesh%rr(1,vv_mesh%jj(1,1:conc_mesh%me))+vv_mesh%rr(1,vv_mesh%jj(2,1:conc_mesh%me)))/2))&
               + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(4,:)) &
               -(vv_mesh%rr(2,vv_mesh%jj(2,1:conc_mesh%me))+vv_mesh%rr(2,vv_mesh%jj(3,1:conc_mesh%me)))/2))&
               + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(5,:)) &
               -(vv_mesh%rr(2,vv_mesh%jj(3,1:conc_mesh%me))+vv_mesh%rr(2,vv_mesh%jj(1,1:conc_mesh%me)))/2))&
               + MAXVAL(ABS(conc_mesh%rr(2,conc_mesh%jj(6,:)) &
               -(vv_mesh%rr(2,vv_mesh%jj(1,1:conc_mesh%me))+vv_mesh%rr(2,vv_mesh%jj(2,1:conc_mesh%me)))/2))
          IF (error/MAXVAL(ABS(vv_mesh%rr(1,1) -vv_mesh%rr(1,:))) .GE. 5.d-14) THEN
             WRITE(*,*) ' WARNING: vv_mesh and conc_mesh do not coincide on the conc domain.'
             WRITE(*,*) ' WARNING: Either you use curved elements P2 elements or BUG, ', &
                  error/MAXVAL(ABS(vv_mesh%rr(1,1) -vv_mesh%rr(1,:)))
          END IF

          error=0.d0
          DO k = 1, conc_mesh%me
             DO n = 1, 2
                error = error+ MAXVAL(ABS(vv_mesh%rr(n,vv_mesh%jj(1:3,k))-conc_mesh%rr(n,conc_mesh%jj(1:3,k))))
             END DO
          END DO
          IF (error/MAXVAL(ABS(vv_mesh%rr(1,1) -vv_mesh%rr(1,:))) .GE. 5.d-14) THEN
             CALL error_Petsc('vv_mesh and conc_mesh do NOT have the same P1 nodes')
          END IF

       END IF

       !===Create periodic structures for concentration===============================
       IF (conc_periodic) THEN
          CALL prep_periodic(inputs%my_periodic, conc_mesh, conc_per)
       END IF

       !===Create csr structure for concentration=====================================
       CALL st_aij_csr_glob_block(comm_one_d_conc(1),1,conc_mesh_glob,conc_mesh,conc_1_LA, opt_per=conc_per)

       !===Deallocate global meshes=================================================
       CALL free_mesh(conc_mesh_glob)

       !===Start Gauss points generation============================================
       !===JLG July 20, 2019, p3 mesh
       conc_mesh%edge_stab    = .FALSE.
       CALL gauss_points_2d(conc_mesh,2)
       !===JLG July 20, 2019, p3 mesh


       !===Create concentration_diffusivity_field=====================================
       ALLOCATE(concentration_diffusivity_field(conc_mesh%me))
       DO m = 1, conc_mesh%me
          DO k=1, inputs%nb_dom_conc
             IF (conc_mesh%i_d(m) == inputs%list_dom_conc(k)) THEN
                concentration_diffusivity_field(m) = inputs%concentration_diffusivity(k)
             END IF
          END DO
       END DO

    END IF
    !===Check coherence of vv_mesh and H_mesh=======================================
    IF ((if_induction .AND. if_momentum) .OR. inputs%type_pb=='mxx') THEN
       IF (vv_mesh%me /=0) THEN
          DO m = 1, vv_mesh%me
             IF (MAXVAL(ABS(H_mesh%rr(:,H_mesh%jj(1:3,m))-vv_mesh%rr(:,vv_mesh%jj(1:3,m))))/=0.d0) THEN
                CALL error_Petsc( ' BUG in init: H_mesh and vv_mesh do not coincide ')
             END IF
          END DO
       END IF
    END IF
    IF (ALLOCATED(list_dom_H)) DEALLOCATE(list_dom_H)
    IF (ALLOCATED(list_dom_H_ref)) DEALLOCATE(list_dom_H_ref)

    !===Check coherence of vv_mesh and temp_mesh====================================
    IF (if_energy) THEN
       IF (vv_mesh%me /=0) THEN
          DO m = 1, vv_mesh%me
             IF (MAXVAL(ABS(temp_mesh%rr(:,temp_mesh%jj(1:3,m))-vv_mesh%rr(:,vv_mesh%jj(1:3,m))))/=0.d0) THEN
                CALL error_Petsc( ' BUG in init: temp_mesh and vv_mesh do not coincide ')
             END IF
          END DO
       END IF
    END IF
    IF (ALLOCATED(list_dom_temp)) DEALLOCATE(list_dom_temp)
    IF (ALLOCATED(list_dom_temp_ref)) DEALLOCATE(list_dom_temp_ref)

    !===Check coherence of vv_mesh and conc_mesh====================================
    IF (if_concentration) THEN
       IF (vv_mesh%me /=0) THEN
          DO m = 1, conc_mesh%me
             IF (MAXVAL(ABS(conc_mesh%rr(:,conc_mesh%jj(1:3,m))-vv_mesh%rr(:,vv_mesh%jj(1:3,m))))/=0.d0) THEN
                CALL error_Petsc( ' BUG in init: conc_mesh and vv_mesh do not coincide ')
             END IF
          END DO
       END IF
    END IF
    IF (ALLOCATED(list_dom_ns)) DEALLOCATE(list_dom_ns)



    !===Compute local mesh size for stabilization================================
    IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
       CALL compute_local_mesh_size(pp_mesh)
       CALL compute_local_mesh_size(vv_mesh)
       CALL compute_local_mesh_size_level_set
    END IF
    IF (if_induction) THEN
       CALL compute_local_mesh_size(pmag_mesh)
       CALL compute_local_mesh_size(phi_mesh)
       CALL compute_local_mesh_size(H_mesh)
    END IF
    IF (if_energy) THEN
       CALL compute_local_mesh_size(temp_mesh)
    END IF
    IF (if_concentration) THEN
       CALL compute_local_mesh_size(conc_mesh)
    END IF

    !===Allocate arrays for Navier-Stokes===========================================
    IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
       IF(inputs%if_navier_stokes_with_taylor) THEN
          ALLOCATE(der_un(inputs%taylor_order-1))
          ALLOCATE(der_pn(inputs%taylor_order-1))
          DO kp = 1, inputs%taylor_order-1
             ALLOCATE(der_un(kp)%DRT(vv_mesh%np, 6, m_max_c))
             ALLOCATE(der_pn(kp)%DRT(pp_mesh%np, 2, m_max_c))
          END DO
       ELSE
          ALLOCATE(un_m1      (vv_mesh%np, 6, m_max_c))
          ALLOCATE(pn_m1      (pp_mesh%np, 2, m_max_c))
          ALLOCATE(incpn_m1   (pp_mesh%np, 2, m_max_c))
          ALLOCATE(incpn      (pp_mesh%np, 2, m_max_c))
          ALLOCATE(density_m2 (vv_mesh%np, 2, m_max_c))
          ALLOCATE(density_m1 (vv_mesh%np, 2, m_max_c))
          ALLOCATE(density    (vv_mesh%np, 2, m_max_c))
       END IF
       ALLOCATE(un         (vv_mesh%np, 6, m_max_c))
       ALLOCATE(pn         (pp_mesh%np, 2, m_max_c))
       IF (inputs%LES) THEN
          ALLOCATE(visc_LES(3, vv_mesh%gauss%l_G*vv_mesh%dom_me, 6, m_max_c))
       END IF
       !===Allocate arrays for Level sets===========================================
       IF (if_mass) THEN
          IF (inputs%if_level_set_P2) THEN
             ALLOCATE(level_set_m1   (inputs%nb_fluid-1, vv_mesh%np, 2, m_max_c))
             ALLOCATE(level_set      (inputs%nb_fluid-1, vv_mesh%np, 2, m_max_c))
             CALL MPI_COMM_SIZE(comm_one_d_ns(2), nb_procs, code)
             bloc_size = vv_mesh%gauss%l_G*vv_mesh%dom_me/nb_procs+1
             bloc_size = vv_mesh%gauss%l_G*(bloc_size/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
             m_max_pad = 3*SIZE(list_mode)*nb_procs/2
             ALLOCATE(visc_entro_level(2*m_max_pad-1,bloc_size))
!TEST LC LES_SUITE 2024/06
             ALLOCATE(visc_LES_level(inputs%nb_fluid-1, vv_mesh%gauss%l_G*vv_mesh%dom_me, 6, m_max_c))
!TEST LC LES_SUITE 2024/06
          ELSE
             ALLOCATE(level_set_m1   (inputs%nb_fluid-1, pp_mesh%np, 2, m_max_c))
             ALLOCATE(level_set      (inputs%nb_fluid-1, pp_mesh%np, 2, m_max_c))
             CALL MPI_COMM_SIZE(comm_one_d_ns(2), nb_procs, code)
             bloc_size = pp_mesh%gauss%l_G*pp_mesh%dom_me/nb_procs+1
             bloc_size = pp_mesh%gauss%l_G*(bloc_size/pp_mesh%gauss%l_G)+pp_mesh%gauss%l_G
             m_max_pad = 3*SIZE(list_mode)*nb_procs/2
             ALLOCATE(visc_entro_level(2*m_max_pad-1,bloc_size))
!TEST LC LES_SUITE 2024/06
             ALLOCATE(visc_LES_level(inputs%nb_fluid-1, pp_mesh%gauss%l_G*pp_mesh%dom_me, 6, m_max_c))
!TEST LC LES_SUITE 2024/06
          END IF
       END IF
    END IF

    !===Allocate arrays for Maxwell=================================================
    IF (if_induction) THEN
       ALLOCATE(Hn1  (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Hn   (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Bn1  (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Bn   (H_mesh%np,  6,  m_max_c))
       ALLOCATE(phin1(phi_mesh%np,2,  m_max_c))
       ALLOCATE(phin (phi_mesh%np,2,  m_max_c))
    END IF

    !===Allocate arrays for temperature=============================================
    IF (if_energy) THEN
       ALLOCATE(tempn_m1   (temp_mesh%np, 2, m_max_c))
       ALLOCATE(tempn      (temp_mesh%np, 2, m_max_c))
    END IF

    !===Allocate arrays for concentration=============================================
    IF (if_concentration) THEN
       ALLOCATE(concn_m1   (conc_mesh%np, 2, m_max_c))
       ALLOCATE(concn      (conc_mesh%np, 2, m_max_c))
    END IF


    !===Create data structure jj_v_to_H=============================================
    IF (if_induction) THEN
       ALLOCATE(v_to_Max(H_mesh%np, 6, m_max_c))
       IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
          ALLOCATE(jj_v_to_H(H_mesh%np))
          jj_v_to_H = -1
          DO m = 1, vv_mesh%me
             jj_v_to_H(H_mesh%jj(:,m)) = vv_mesh%jj(1:H_mesh%gauss%n_w,m)
          END DO
       ELSE
          ALLOCATE(jj_v_to_H(H_mesh%np))
          jj_v_to_H = -1
       END IF
    END IF
    IF (if_momentum) THEN
       ALLOCATE(H_to_NS(vv_mesh%np, 6, m_max_c))
       ALLOCATE(B_to_NS(vv_mesh%np, 6, m_max_c))
       ALLOCATE(Hext(H_mesh%np, 6, m_max_c))
       ALLOCATE(Bext(H_mesh%np, 6, m_max_c))
    ELSE
       ALLOCATE(H_to_NS(1, 1, 1))
       ALLOCATE(B_to_NS(1, 1, 1))
    END IF

    !===Create data structure jj_v_to_temp==========================================
    IF (if_energy) THEN
       ALLOCATE(v_to_energy(temp_mesh%np, 6, m_max_c))
       ALLOCATE(jj_v_to_temp(temp_mesh%np))
       jj_v_to_temp = -1
       IF (vv_mesh%me/=0) THEN
          DO m = 1, vv_mesh%me
             jj_v_to_temp(temp_mesh%jj(:,m)) = vv_mesh%jj(1:temp_mesh%gauss%n_w,m)
          END DO
          ALLOCATE(T_to_NS(vv_mesh%np, 2, m_max_c))
       ELSE
          ALLOCATE(T_to_NS(1, 1, 1))
       END IF
    END IF

    !===Create data structure jj_c_to_v==========================================
    IF (if_concentration) THEN
       ALLOCATE(conc_to_v(vv_mesh%np, 2, m_max_c))
       ALLOCATE(jj_c_to_v(vv_mesh%np))
       jj_c_to_v = -1
       IF (conc_mesh%me/=0) THEN
          DO m = 1, conc_mesh%me
             jj_c_to_v(vv_mesh%jj(:,m)) = conc_mesh%jj(1:vv_mesh%gauss%n_w,m)
          END DO
          ALLOCATE(v_to_conc(conc_mesh%np, 6, m_max_c))
       ELSE
          ALLOCATE(v_to_conc(1, 1, 1))
       END IF
    END IF

    !===Create data structure jj_T_to_H==========================================
    IF (if_induction) THEN
       ALLOCATE(T_to_H(H_mesh%np, 2, m_max_c))
       IF (if_energy) THEN
          ALLOCATE(jj_T_to_H(H_mesh%np))
          jj_T_to_H = -1
          DO m = 1, temp_mesh%me
             jj_T_to_H(H_mesh%jj(:,m)) = temp_mesh%jj(1:H_mesh%gauss%n_w,m)
          END DO
          ALLOCATE(j_H_to_T(temp_mesh%np, 6, m_max_c))
       ELSE
          ALLOCATE(jj_T_to_H(H_mesh%np))
          jj_T_to_H = -1
          ALLOCATE(j_H_to_T(1, 1, 1))
       END IF
    END IF

    !===Create data structure jj_c_to_H==========================================
    IF (if_induction) THEN
       ALLOCATE(conc_to_H(H_mesh%np, 2, m_max_c))
       IF (if_concentration) THEN
          ALLOCATE(jj_c_to_H(H_mesh%np))
          jj_c_to_H = -1
          DO m = 1, conc_mesh%me
             jj_c_to_H(H_mesh%jj(:,m)) = conc_mesh%jj(1:H_mesh%gauss%n_w,m)
          END DO
          ALLOCATE(j_H_to_conc(conc_mesh%np, 6, m_max_c))
       ELSE
          ALLOCATE(jj_c_to_H(H_mesh%np))
          jj_c_to_H = -1
          ALLOCATE(j_H_to_conc(1, 1, 1))
       END IF
    END IF

    !===Create coupling variables H and pdt H to temp===============================
    IF (if_energy .AND. if_induction) THEN
       ALLOCATE(H_to_energy(temp_mesh%np, 6, m_max_c))
       H_to_energy = 0.d0
       ALLOCATE(pdt_H_to_energy(temp_mesh%np, 6, m_max_c))
       pdt_H_to_energy = 0.d0
    END IF

    !===Initialize Navier-Stokes====================================================
    time_u = 0.d0
    IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
       IF (vv_mesh%me/=0) THEN
          IF (inputs%irestart_u) THEN
             IF(inputs%if_navier_stokes_with_taylor) THEN
                CALL read_restart_ns_taylor(comm_one_d_ns, time_u, list_mode, un, der_un, pn, der_pn, inputs%file_name)
             ELSE
                IF (.NOT. if_mass) THEN
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un, un_m1, pn, pn_m1, &
                        incpn, incpn_m1, inputs%file_name)
!TEST LC LES_SUITE 2024/06
                   IF (inputs%LES) THEN
                      IF (inputs%irestart_LES) THEN
                         CALL read_restart_LES(comm_one_d_ns, time_u, list_mode, inputs%file_name, opt_LES_NS=visc_LES)
                      ELSE
                         visc_LES = 0.d0
                      END IF
                   END IF
!TEST LC LES_SUITE 2024/06
                ELSE
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un, un_m1, pn, pn_m1, &
                        incpn, incpn_m1, inputs%file_name, opt_level_set=level_set, &
                        opt_level_set_m1=level_set_m1, opt_max_vel=max_vel)
!TEST LC LES_SUITE 2024/06
                   IF (inputs%irestart_LES) THEN !inputs%LES=.t. true when inputs%if_mass=.t.
                      IF (inputs%if_LES_in_momentum) THEN
                         CALL read_restart_LES(comm_one_d_ns, time_u, list_mode, inputs%file_name, &
                             opt_LES_NS=visc_LES, opt_LES_level=visc_LES_level)
                      ELSE
                         CALL read_restart_LES(comm_one_d_ns, time_u, list_mode, inputs%file_name, &
                             opt_LES_level=visc_LES_level)
                         visc_LES = 0.d0
                      END IF
                   ELSE
                      visc_LES = 0.d0
                      visc_LES_level = 0.d0
                   END IF
!TEST LC LES_SUITE 2024/06
                END IF
             END IF
          ELSE
             IF(inputs%if_navier_stokes_with_taylor) THEN
                time_u = 0.d0 !===ATTENTION: Fixed initialization time for Taylor method
                CALL init_velocity_pressure_taylor(vv_mesh, pp_mesh, list_mode, time_u, pn, der_pn, un, der_un)
             ELSE
                CALL init_velocity_pressure(vv_mesh, pp_mesh, time_u, &
                     inputs%dt, list_mode, un_m1, un, pn_m1, pn, incpn_m1, incpn)
             END IF
             IF (inputs%LES) THEN
                visc_LES = 0.d0
             END IF
             IF (if_mass) THEN
                IF (inputs%if_level_set_P2) THEN
                   CALL init_level_set(vv_mesh, time_u, &
                        inputs%dt, list_mode, level_set_m1, level_set)
                ELSE
                   CALL init_level_set(pp_mesh, time_u, &
                        inputs%dt, list_mode, level_set_m1, level_set)
                END IF
                bloc_size = SIZE(un,1)/inputs%ndim(2) + 1
                m_max_pad = 3*SIZE(list_mode)*inputs%ndim(2)/2
                CALL FFT_MAX_VEL_DCL(comm_one_d_ns(2), 2*un-un_m1, max_vel_S, inputs%ndim(2), bloc_size, m_max_pad)
                CALL MPI_ALLREDUCE(max_vel_S, max_vel, 1, MPI_DOUBLE_PRECISION, &
                     MPI_MAX, comm_one_d_ns(1), code)
                max_vel = MAX(1.1d0*max_vel, 0.1d0)
!TEST LC LES_SUITE 2024/06
                IF (inputs%LES) THEN
                   visc_LES_level = 0.d0
                END IF
!TEST LC LES_SUITE 2024/06
             END IF
          END IF

          IF (if_mass) THEN
             CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set_m1, &
                  inputs%density_fluid, density_m1)
             CALL reconstruct_variable(comm_one_d_ns, list_mode, pp_mesh, vv_mesh, level_set, &
                  inputs%density_fluid, density)
             visc_entro_level = 0.d0
          END IF

          !===Force sine coefficients of zero mode to zero==========================
          DO i = 1, m_max_c
             IF (list_mode(i) == 0) THEN
                DO k= 1, 3
                   un(:,2*k,i)       = 0.d0
                   IF(inputs%if_navier_stokes_with_taylor) THEN
                      DO kp = 1, inputs%taylor_order-1
                         der_un(kp)%DRT( :,2*k,i)= 0.d0
                      END DO
                   ELSE
                      un_m1(:,2*k,i) = 0.d0
                   END IF
                END DO
                pn(:,2,i)         = 0.d0
                IF(inputs%if_navier_stokes_with_taylor) THEN
                   DO kp = 1, inputs%taylor_order-1
                      der_pn(kp)%DRT( :,2,i)  = 0.d0
                   END DO
                ELSE
                   incpn(:,2,i)      = 0.d0
                   density(:,2,i)    = 0.d0
                   pn_m1(:,2,i)      = 0.d0
                   incpn_m1(:,2,i)   = 0.d0
                   density_m1(:,2,i) = 0.d0
                END IF
                IF (if_mass) THEN
                   level_set(:,:,2,i)       = 0.d0
                   level_set_m1(:,:,2,i)    = 0.d0
                END IF
             END IF
          END DO
          !===End force sine coefficients of zero mode to zero======================
          IF(.NOT.inputs%if_navier_stokes_with_taylor) THEN
             density_m2 = density_m1
          END IF
       END IF
    END IF

    !===Initialize velocity (time-independent) if mxw===============================
    IF ( (inputs%type_pb=='mxw') .AND. (H_mesh%me/=0) ) THEN
       DO i = 1, m_max_c       !===Initialization of vel
          v_to_Max(:,:,i) = Vexact(list_mode(i), H_mesh)
       END DO
    ENDIF

    !===Initialize velocity using un (time-independent) if mxx======================
    IF (inputs%type_pb=='mxx') THEN
       !===Use extension_velocity===================================================
       IF (H_mesh%np>0) THEN !We extend v_to_Max
          DO i = 1, m_max_c
             DO k= 1, 6 !===The user has to code extension_vel
                v_to_Max(:,k,i) = extension_velocity(k, H_mesh, list_mode(i), time_u, 1)
             END DO
          END DO
       END IF

       !===Use extension_velocity===================================================
       IF (vv_mesh%me /=0) THEN
          DO j = 1, SIZE(jj_v_to_H)
             IF (jj_v_to_H(j) == -1) CYCLE
             v_to_Max(j,:,:) = un(jj_v_to_H(j),:,:)
          END DO
       END IF

       !===Cleanup==================================================================
       IF(inputs%if_navier_stokes_with_taylor) THEN
          DO kp = 1, inputs%taylor_order-1
             DEALLOCATE(der_un(kp)%DRT)
             DEALLOCATE(der_pn(kp)%DRT)
          END DO
          DEALLOCATE(der_un)
          DEALLOCATE(der_pn)
       ELSE
          DEALLOCATE(un_m1)
          DEALLOCATE(pn_m1)
          DEALLOCATE(incpn_m1)
          DEALLOCATE(incpn)
          DEALLOCATE(density_m2)
          DEALLOCATE(density_m1)
          DEALLOCATE(density)
       END IF
       DEALLOCATE(un)
       DEALLOCATE(pn)
       IF (if_mass) THEN
          IF (ALLOCATED(level_set)) DEALLOCATE(level_set,level_set_m1)
          IF (ALLOCATED(visc_LES_level)) DEALLOCATE(visc_LES_level)
       END IF
       IF (inputs%LES) THEN
          IF (ALLOCATED(visc_LES)) DEALLOCATE(visc_LES)
       END IF
    ENDIF

    !===Initialize Maxwell==========================================================
    time_h = 0.d0
    IF (if_induction) THEN
       IF (inputs%irestart_h) THEN
          CALL read_restart_maxwell(comm_one_d, H_mesh, phi_mesh, time_h, list_mode, Hn, Hn1, Bn, Bn1, &
               phin, phin1, inputs%file_name)
       ELSE
          CALL init_maxwell(H_mesh,phi_mesh,time_h,inputs%dt,mu_H_field,inputs%mu_phi,list_mode,&
               Hn1,Hn,phin1,phin)
          !===Initialize Bn and Bn1
          IF (H_mesh%me/=0) THEN
             IF (inputs%if_permeability_variable_in_theta) THEN
                CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
                m_max_pad = 3*SIZE(list_mode)*nb_procs/2
                bloc_size = SIZE(Bn,1)/nb_procs+1
                CALL FFT_PAR_VAR_ETA_PROD_T_DCL(comm_one_d(2), mu_in_real_space, &
                     H_mesh, Hn, Bn, nb_procs, bloc_size, m_max_pad, time)
                CALL FFT_PAR_VAR_ETA_PROD_T_DCL(comm_one_d(2), mu_in_real_space, &
                     H_mesh, Hn1, Bn1, nb_procs, bloc_size, m_max_pad, time)
             ELSE
                DO i = 1, m_max_c
                   DO k = 1, 6
                      Bn(:,k,i)  =  mu_H_field*Hn(:,k,i)
                      Bn1(:,k,i) =  mu_H_field*Hn1(:,k,i)
                   END DO
                END DO
             END IF
          END IF
       END IF
       !===Force sine coefficients of zero mode to zero=============================
       DO i = 1, m_max_c
          IF (list_mode(i) == 0) THEN
             IF (H_mesh%me/=0) THEN
                DO k = 1, 3
                   Hn(:,2*k,i)  = 0.d0
                   Hn1(:,2*k,i) = 0.d0
                END DO
             END IF
             IF (phi_mesh%me/=0) THEN
                phin(:,2,i)  = 0.d0
                phin1(:,2,i) = 0.d0
             END IF
          END IF
       END DO
    END IF

    !===Initialize temperature======================================================
    time_T = 0.d0
    IF (if_energy) THEN
       IF (temp_mesh%me/=0) THEN
          IF (inputs%irestart_T) THEN
             CALL read_restart_temp(comm_one_d_temp, time_T, list_mode, tempn, tempn_m1, &
                  inputs%file_name)
          ELSE
             CALL init_temperature(temp_mesh, time_T, inputs%dt, list_mode, tempn_m1, tempn)
          END IF
          !===Force sine coefficients of zero mode to zero==========================
          DO i = 1, m_max_c
             IF (list_mode(i) == 0) THEN
                tempn(:,2,i)       = 0.d0
                tempn_m1(:,2,i)    = 0.d0
             END IF
          END DO

       END IF
    END IF


    !===Initialize concentration======================================================
    time_conc = 0.d0
    IF (if_concentration) THEN
       IF (conc_mesh%me/=0) THEN
          IF (inputs%irestart_conc) THEN
             CALL read_restart_conc(comm_one_d_conc, time_conc, list_mode, concn, concn_m1, &
                  inputs%file_name)
          ELSE
             CALL init_concentration(conc_mesh, time_conc, inputs%dt, list_mode, concn_m1, concn)
          END IF
          !===Force sine coefficients of zero mode to zero==========================
          DO i = 1, m_max_c
             IF (list_mode(i) == 0) THEN
                concn(:,2,i)       = 0.d0
                concn_m1(:,2,i)    = 0.d0
             END IF
          END DO

       END IF
    END IF

    !===Initialize time=============================================================
    IF (inputs%irestart_h .OR. inputs%irestart_u .OR. inputs%irestart_T .OR. inputs%irestart_conc) THEN
       time = MAX(time_u,time_h,time_T,time_conc)
    ELSE
       time = 0.d0
    END IF

    !===Guardrail===================================================================
    IF (if_mass.AND.inputs%variation_sigma_fluid) THEN
       IF (inputs%analytical_permeability.OR.inputs%nb_dom_phi>0) THEN
          CALL error_Petsc(' BUG in INIT : sigma reconstruct via level set not implemented with vacuum or variable mu')
       END IF
    END IF

  END SUBROUTINE INIT

  SUBROUTINE prodmat_maxwell_int_by_parts(vect_in, vect_out, ndim, i)
    USE update_maxwell
    IMPLICIT NONE
    INTEGER :: ndim
    REAL(KIND=8), DIMENSION(ndim) :: vect_in
    REAL(KIND=8), DIMENSION(ndim) :: vect_out
    INTEGER                       :: i
    INTEGER :: TYPE, i_deb, i_fin
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode))  :: sigma_ns

    time = 0.d0
    DO TYPE = 1, 6
       i_deb = (TYPE-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(TYPE,2)==0 .AND. list_mode(i)==0) THEN
          Hn(:,TYPE,i)   = 0.d0
       ELSE
          Hn(:,TYPE,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
    DO TYPE = 1, 2
       phin(:,TYPE,i) =0.d0
    END DO

    DO TYPE = 1, 6
       i_deb = 6*H_mesh%np + (TYPE-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       IF (MODULO(TYPE,2)==0 .AND. list_mode(i)==0) THEN
          Hn1(:,TYPE,i)   = 0.d0
       ELSE
          Hn1(:,TYPE,i)   = vect_in(i_deb:i_fin)
       END IF
    END DO
    DO TYPE = 1, 2
       phin1(:,TYPE,i) =0.d0
    END DO

    CALL maxwell_decouple(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, interface_H_mu, &
         Hn, Bn, phin, Hn1, Bn1, phin1, v_to_Max, inputs%stab, inputs%stab_jump_h,sigma_field, R_fourier, index_fourier, &
         mu_H_field, inputs%mu_phi, time, inputs%dt, inputs%Rem, list_mode, H_phi_per, LA_H, LA_pmag, &
         LA_phi, LA_mhd, sigma_ns, jj_v_to_H,conc_to_H)

    DO TYPE = 1, 6
       i_deb = (TYPE-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn(:,TYPE,i)
    END DO

    DO TYPE = 1, 6
       i_deb = 6*H_mesh%np + (TYPE-1)*H_mesh%np+1
       i_fin = i_deb + H_mesh%np -1
       vect_out(i_deb:i_fin) = Hn1(:,TYPE,i)
    END DO

  END SUBROUTINE prodmat_maxwell_int_by_parts

!!$  SUBROUTINE sfemansinitialize
!!$    IMPLICIT NONE
!!$    INTEGER                                          :: narg, i
!!$    CHARACTER(len=200),DIMENSION(:,:), ALLOCATABLE   :: inline
!!$    LOGICAL                                          :: ok
!!$    CHARACTER(len=3)                                 :: tit
!!$
!!$    narg = 0
!!$    ok = .TRUE.
!!$
!!$    DO WHILE (ok)
!!$       CALL getarg(narg+1,tit)
!!$       IF (tit == '   ') THEN
!!$          ok = .FALSE.
!!$       ELSE
!!$          narg = narg+1
!!$       END IF
!!$    END DO
!!$
!!$    narg = narg/2
!!$    ALLOCATE(inline(2,narg))
!!$
!!$    DO i = 1, narg
!!$       CALL getarg(2*(i-1)+1,inline(1,i))
!!$       CALL getarg(2*i      ,inline(2,i))
!!$    END DO
!!$
!!$    inputs%test_de_convergence = .FALSE.
!!$    inputs%numero_du_test_debug = 0
!!$    inputs%data_directory_debug = '.'
!!$    DO i = 1, narg
!!$       IF (TRIM(ADJUSTL(inline(1,i))) == 'debug') THEN
!!$          inputs%test_de_convergence = .TRUE.
!!$          READ(inline(2,i),*) inputs%numero_du_test_debug
!!$       ELSE IF (TRIM(ADJUSTL(inline(1,i))) == 'debug_dir') THEN
!!$          inputs%data_directory_debug=inline(2,i)
!!$       END IF
!!$    END DO
!!$  END SUBROUTINE sfemansinitialize

  SUBROUTINE regression_initialize
    IMPLICIT NONE
    INTEGER                                          :: narg
    CHARACTER(len=200)                               :: inline
    LOGICAL                                          :: ok
    CHARACTER(len=3)                                 :: tit

    narg = 0
    ok = .TRUE.

    DO WHILE (ok)
       CALL getarg(narg+1,tit)
       IF (tit == '   ') THEN
          ok = .FALSE.
       ELSE
          narg = narg+1
       END IF
    END DO

    inputs%if_regression = .FALSE.
    IF (narg==1) THEN
       CALL getarg(1,inline)
       IF (TRIM(ADJUSTL(inline)) == 'regression') THEN
          inputs%if_regression = .TRUE.
       END IF
    END IF
  END SUBROUTINE regression_initialize

  SUBROUTINE compute_local_mesh_size(mesh)
    USE def_type_mesh
    TYPE(mesh_type) :: mesh
    INTEGER         :: m, type_fe, index, l, ierr, i
    REAL(KIND=8)    :: diam_loc, diam

    IF (mesh%gauss%n_w==3) THEN
       type_fe = 1
    ELSE IF (mesh%gauss%n_w==6) THEN
       type_fe = 2
    ELSE IF (mesh%gauss%n_w==10) THEN
       type_fe = 3
    ELSE
       CALL error_petsc('BUG in compute_local_mesh_size')
    END IF
    ALLOCATE(mesh%hloc(mesh%dom_me))
    ALLOCATE(mesh%hloc_gauss(mesh%gauss%l_G*mesh%dom_me))

    index = 0
    DO m = 1, mesh%dom_me
       mesh%hloc(m) = SQRT(SUM(mesh%gauss%rj(:,m)))/type_fe
       DO l = 1, mesh%gauss%l_G
          index = index + 1
          mesh%hloc_gauss(index) = mesh%hloc(m)
       END DO
    END DO

    !===JLG+CN Dec 14 2016
    !===Compute characteristic diameter of meridian section of domain associated with mesh
    !===BUG diam_loc = SUM(mesh%gauss%rj) !===BUG JLG+HF (not defined if mesh%dom_me=0)
    if(mesh%dom_me==0) THEN !===JLG+HF Dec 9 2019
       diam_loc =0.d0
    else
       diam_loc = SUM(mesh%gauss%rj)
    end if
    CALL MPI_ALLREDUCE(diam_loc,diam,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm_one_d(1),ierr)
    mesh%global_diameter = SQRT(diam)
    !===End JLG+CN Dec 14 2016
    !===hm (JLG April 7, 2017)
    diam_loc = MAXVAL(mesh%rr(1,:))
    CALL MPI_ALLREDUCE(diam_loc,diam,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_one_d(1),ierr)
    ALLOCATE(mesh%hm(m_max_c))
    DO i = 1, m_max_c
       mesh%hm(i) = 0.5d0*diam/inputs%m_max
    END DO
    !===end hm (JLG April 7, 2017)
  END SUBROUTINE compute_local_mesh_size

  SUBROUTINE compute_local_mesh_size_level_set
    REAL(KIND=8)    :: h_min, h_min_F
    INTEGER         :: code
    !===Compute h_min
    IF (inputs%if_level_set_P2) THEN
       h_min_F=MINVAL(vv_mesh%hloc_gauss)
    ELSE
       h_min_F=MINVAL(pp_mesh%hloc_gauss)
    END IF
    CALL MPI_ALLREDUCE(h_min_F, h_min, 1, MPI_DOUBLE_PRECISION, &
         MPI_MIN, comm_one_d_ns(1), code)
    inputs%h_min_distance = inputs%multiplier_for_h_min_distance*h_min
    !===End compute h_min
  END SUBROUTINE compute_local_mesh_size_level_set

END MODULE initialization
