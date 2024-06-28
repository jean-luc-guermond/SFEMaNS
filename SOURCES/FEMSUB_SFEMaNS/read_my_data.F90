MODULE my_data_module
  USE def_type_mesh
  USE solve_petsc
  TYPE my_data
     LOGICAL                                 :: iformatted
     CHARACTER(LEN=200)                      :: directory, file_name
     LOGICAL                                 :: if_read_partition
     CHARACTER(LEN=3)                        :: type_pb
     INTEGER, DIMENSION(2)                   :: ndim
     LOGICAL                                 :: select_mode
     INTEGER                                 :: m_max
     INTEGER, DIMENSION(:), POINTER          :: list_mode_lect
     !===Data time iterations========================================================
     INTEGER                                 :: nb_iteration
     REAL(KIND=8)                            :: dt
     !===Data for LES================================================================
     LOGICAL                                 :: LES
     REAL(KIND=8)                            :: LES_coeff1, LES_coeff2, LES_coeff3, LES_coeff4
     REAL(KIND=8)                            :: LES_coeff1_mom
     LOGICAL                                 :: if_LES_in_momentum
     !===Data for precession=========================================================
     REAL(KIND=8)                            :: taux_precession, angle_s_pi
     LOGICAL                                 :: precession
     !===Data for NS penalization====================================================
     LOGICAL                                 :: if_ns_penalty
     LOGICAL                                 :: if_impose_vel_in_solids
     LOGICAL                                 :: if_compute_momentum_pseudo_force
     REAL(KIND=8)                            :: div_stab_in_ns
     !===Data for Navier-Stokes======================================================
     LOGICAL                                 :: if_navier_stokes_with_u
     LOGICAL                                 :: if_navier_stokes_art_comp
     REAL(KIND=8)                            :: penal_coeff_art_comp
     LOGICAL                                 :: if_tensor_sym
     LOGICAL                                 :: if_moment_bdf2
     LOGICAL                                 :: if_temp_bdf2
     LOGICAL                                 :: if_level_bdf2
     LOGICAL                                 :: irestart_u
     LOGICAL                                 :: irestart_LES
     LOGICAL                                 :: if_variable_visco
     LOGICAL                                 :: if_navier_stokes_with_taylor !===JLG july 20, 2019, p3 mesh
     INTEGER                                 :: taylor_order !===JLG Dec 2020
     REAL(KIND=8)                            :: taylor_lambda!===JLG Dec 2020
     INTEGER                                 :: type_fe_velocity !===JLG july 20, 2019, p3 mesh
     INTEGER                                 :: nb_dom_ns
     INTEGER, DIMENSION(:), POINTER          :: list_dom_ns
     REAL(KIND=8)                            :: Re
     REAL(KIND=8)                            :: coeff_lorentz
     TYPE(solver_param)                      :: my_par_vv, my_par_pp, my_par_mass
     INTEGER                                 :: pp_nb_dirichlet_sides
     INTEGER, DIMENSION(:), POINTER          :: pp_list_dirichlet_sides
     INTEGER, DIMENSION(3)                   :: vv_nb_dirichlet_sides
     TYPE(dyn_int_line), DIMENSION(3)        :: vv_list_dirichlet_sides
     INTEGER                                 :: vv_nb_dirichlet
     INTEGER                                 :: vv_nb_dirichlet_normal_velocity
     INTEGER, DIMENSION(:), POINTER          :: vv_list_dirichlet_normal_velocity_sides
     !===Data for Maxwell============================================================
     LOGICAL                                 :: if_maxwell_with_H
     LOGICAL                                 :: irestart_h
     INTEGER                                 :: nb_dom_H, nb_dom_phi, nb_inter, nb_inter_mu
     INTEGER                                 :: type_fe_H, type_fe_phi, nb_dirichlet_sides_H
     INTEGER, DIMENSION(:), POINTER          :: list_dom_H, list_dom_phi, list_inter_H_phi
     INTEGER, DIMENSION(:), POINTER          :: list_inter_mu, list_dirichlet_sides_H, list_inter_rot_h_jump
     REAL(KIND=8)                            :: Rem, mu_phi
     LOGICAL                                 :: analytical_permeability
     LOGICAL                                 :: if_use_fem_integration_for_mu_bar
     LOGICAL                                 :: if_permeability_variable_in_theta
     REAL(KIND=8), DIMENSION(:), POINTER     :: mu_H, sigma
     REAL(KIND=8), DIMENSION(3)              :: stab
     TYPE(solver_param)                      :: my_par_H_p_phi
     INTEGER                                 :: phi_nb_dirichlet_sides
     INTEGER, DIMENSION(:), POINTER          :: phi_list_dirichlet_sides
     LOGICAL                                 :: if_quasi_static_approx
     LOGICAL                                 :: if_steady_current_fhd
     REAL(KIND=8)                            :: stab_jump_h
     INTEGER                                 :: rot_h_nb_jump_sides
     !===Data for temperature========================================================
     LOGICAL                                 :: if_temperature
     LOGICAL                                 :: if_temperature_with_T
     LOGICAL                                 :: irestart_T
     LOGICAL                                 :: if_helmholtz_force
     REAL(KIND=8)                            :: gravity_coefficient
     REAL(KIND=8)                            :: mag_force_coefficient
     REAL(KIND=8), DIMENSION(:), POINTER     :: temperature_diffusivity
     REAL(KIND=8), DIMENSION(:), POINTER     :: density, heat_capacity, vol_heat_capacity
     TYPE(solver_param)                      :: my_par_temperature
     INTEGER                                 :: temperature_nb_dirichlet_sides
     INTEGER, DIMENSION(:), POINTER          :: temperature_list_dirichlet_sides
     INTEGER                                 :: temperature_nb_robin_sides
     INTEGER, DIMENSION(:), POINTER          :: temperature_list_robin_sides
     INTEGER                                 :: temperature_nb_neumann_sides
     INTEGER, DIMENSION(:), POINTER          :: temperature_list_neumann_sides
     REAL(KIND=8), DIMENSION(:), POINTER     :: convection_coeff
     REAL(KIND=8), DIMENSION(:), POINTER     :: exterior_temperature
     INTEGER                                 :: nb_dom_temp
     INTEGER, DIMENSION(:), POINTER          :: list_dom_temp
     INTEGER                                 :: nb_inter_v_T
     INTEGER, DIMENSION(:), POINTER          :: list_inter_v_T
     !===Data for concentration======================================================
     LOGICAL                                 :: if_concentration
     LOGICAL                                 :: irestart_conc
     REAL(KIND=8), DIMENSION(:), POINTER     :: concentration_diffusivity
     TYPE(solver_param)                      :: my_par_concentration
     INTEGER                                 :: concentration_nb_dirichlet_sides
     INTEGER, DIMENSION(:), POINTER          :: concentration_list_dirichlet_sides
     INTEGER                                 :: concentration_nb_robin_sides
     INTEGER, DIMENSION(:), POINTER          :: concentration_list_robin_sides
     INTEGER                                 :: concentration_nb_neumann_sides
     INTEGER, DIMENSION(:), POINTER          :: concentration_list_neumann_sides
     REAL(KIND=8), DIMENSION(:), POINTER     :: convection_coeff_conc_lhs, convection_coeff_conc_rhs
     REAL(KIND=8), DIMENSION(:), POINTER     :: exterior_concentration
     INTEGER                                 :: nb_dom_conc
     INTEGER, DIMENSION(:), POINTER          :: list_dom_conc
     INTEGER                                 :: nb_inter_c_v
     INTEGER, DIMENSION(:), POINTER          :: list_inter_c_v
     !===Data for mhs problems=======================================================
     LOGICAL                                 :: if_coupling_H_x
     LOGICAL                                 :: if_coupling_analytical
     REAL(KIND=8)                            :: MA
     REAL(KIND=8)                            :: MB
     REAL(KIND=8)                            :: rho_0percent
     REAL(KIND=8)                            :: xi
     REAL(KIND=8)                            :: pot_slope
     REAL(KIND=8)                            :: faraday_cst
     !===Data for level set==========================================================
     LOGICAL                                 :: if_level_set
     LOGICAL                                 :: if_level_set_fixed
     TYPE(solver_param)                      :: my_par_level_set
     INTEGER                                 :: level_set_nb_dirichlet_sides
     INTEGER, DIMENSION(:), POINTER          :: level_set_list_dirichlet_sides
     REAL(KIND=8)                            :: level_set_cmax
     REAL(KIND=8)                            :: level_set_comp_factor
     CHARACTER(LEN=3)                        :: level_set_reconstruction_type
     REAL(KIND=8)                            :: level_set_tanh_coeff_reconstruction
     INTEGER                                 :: nb_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: density_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: dyna_visc_fluid
     LOGICAL                                 :: variation_sigma_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: sigma_fluid
     LOGICAL                                 :: variation_temp_param_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: heat_capacity_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: heat_diffu_fluid
     REAL(KIND=8), DIMENSION(:), POINTER     :: heat_grav_fluid
     LOGICAL                                 :: if_surface_tension
     REAL(KIND=8), DIMENSION(:), POINTER     :: coeff_surface
     LOGICAL                                 :: if_mass_correction
     LOGICAL                                 :: if_kill_overshoot
     REAL(KIND=8)                            :: multiplier_for_h_min_distance
     REAL(KIND=8)                            :: sigma_min, mu_min
     !===Computed data for level set
     REAL(KIND=8)                            :: h_min_distance
     LOGICAL                                 :: if_level_set_P2
     LOGICAL                                 :: if_compression_mthd_JLG
     !===Data for periodicity========================================================
     TYPE(periodic_data)                     :: my_periodic
     !===Data for convergence tests==================================================
     LOGICAL                                 :: if_regression
     LOGICAL                                 :: test_de_convergence
     CHARACTER(len=200)                      :: data_directory_debug
     INTEGER                                 :: numero_du_test_debug
     REAL(KIND=8), DIMENSION(4)              :: norm_ref
     !===Data for Arpack=============================================================
     LOGICAL                                 :: if_arpack
     CHARACTER(len=2)                        :: arpack_which
     INTEGER                                 :: nb_eigenvalues, arpack_iter_max
     REAL(KIND=8)                            :: arpack_tolerance
     LOGICAL                                 :: if_arpack_vtu_2d
     !===Data for stress bc==========================================================
     REAL(KIND=8)                            :: stab_bdy_ns
     !===Data for postprocessing=====================================================
     INTEGER                                 :: number_of_planes_in_real_space
     LOGICAL                                 :: check_numerical_stability
     LOGICAL                                 :: is_mesh_symmetric
     LOGICAL                                 :: if_zero_out_modes
     INTEGER                                 :: nb_select_mode_ns, nb_select_mode_mxw
     INTEGER, DIMENSION(:), POINTER          :: list_select_mode_ns, list_select_mode_mxw
     INTEGER                                 :: freq_restart, freq_en, freq_plot
     LOGICAL                                 :: verbose_timing, verbose_divergence, verbose_CFL
     LOGICAL                                 :: if_just_processing
     LOGICAL                                 :: if_post_proc_init
     LOGICAL                                 :: if_xml
     LOGICAL                                 :: if_plot_2D
     LOGICAL                                 :: if_compute_error
     !===Data for anemometer (postprocessing)========================================
     LOGICAL                                 :: if_anemo_conc, if_anemo_v, if_anemo_T, if_anemo_h
     INTEGER                                 :: nb_anemo_r_conc, nb_anemo_z_conc
     INTEGER                                 :: nb_anemo_r_v, nb_anemo_z_v
     INTEGER                                 :: nb_anemo_r_T, nb_anemo_z_T
     INTEGER                                 :: nb_anemo_r_h, nb_anemo_z_h
     REAL(KIND=8), DIMENSION(:), POINTER     :: r_anemo_conc, z_anemo_conc
     REAL(KIND=8), DIMENSION(:), POINTER     :: r_anemo_v, z_anemo_v
     REAL(KIND=8), DIMENSION(:), POINTER     :: r_anemo_T, z_anemo_T
     REAL(KIND=8), DIMENSION(:), POINTER     :: r_anemo_h, z_anemo_h

   CONTAINS
     PROCEDURE, PUBLIC                       :: init
  END TYPE my_data
CONTAINS
  SUBROUTINE init(a)
    CLASS(my_data), INTENT(INOUT) :: a
    !===Logicals
    a%iformatted=.FALSE.
    a%if_read_partition=.FALSE.
    a%select_mode=.FALSE.
    a%LES=.FALSE.
    a%if_LES_in_momentum=.FALSE.
    a%precession=.FALSE.
    a%if_ns_penalty=.FALSE.
    a%if_impose_vel_in_solids=.FALSE.
    a%if_compute_momentum_pseudo_force=.FALSE.
    a%if_navier_stokes_with_u=.FALSE.
    a%if_navier_stokes_art_comp=.FALSE.
    a%if_tensor_sym=.FALSE.
    a%if_moment_bdf2=.FALSE.
    a%if_temp_bdf2=.FALSE.
    a%if_level_bdf2=.FALSE.
    a%irestart_u=.FALSE.
    a%irestart_LES=.FALSE.
    a%if_maxwell_with_H=.FALSE.
    a%irestart_h=.FALSE.
    a%irestart_T=.FALSE.
    a%irestart_conc=.FALSE.
    a%if_helmholtz_force=.FALSE.
    a%analytical_permeability=.FALSE.
    a%if_use_fem_integration_for_mu_bar=.FALSE.
    a%if_permeability_variable_in_theta=.FALSE.
    a%if_quasi_static_approx=.FALSE.
    a%if_steady_current_fhd=.FALSE.
    a%if_temperature=.FALSE.
    a%if_temperature_with_T=.FALSE.
    a%if_concentration=.FALSE.
    a%if_level_set=.FALSE.
    a%if_level_set_fixed=.FALSE.
    a%variation_sigma_fluid=.FALSE.
    a%if_surface_tension=.FALSE.
    a%if_mass_correction=.FALSE.
    a%if_kill_overshoot=.FALSE.
    a%if_level_set_P2=.FALSE.
    a%if_compression_mthd_JLG=.FALSE.
    a%if_variable_visco=.FALSE.
    a%if_xml = .TRUE.
    a%if_plot_2D=.FALSE.
    a%if_compute_error=.FALSE.
    a%if_anemo_conc=.FALSE.
    a%if_anemo_v=.FALSE.
    a%if_anemo_T=.FALSE.
    a%if_anemo_h=.FALSE.
    a%if_navier_stokes_with_taylor =.FALSE.
    a%if_coupling_H_x=.FALSE.
    a%if_coupling_analytical=.FALSE.
    !===Done in sfemansinitialize. Do not touch. a%test_de_convergence
    a%if_arpack=.FALSE.
    a%if_arpack_vtu_2d=.FALSE.
    a%check_numerical_stability=.FALSE.
    a%is_mesh_symmetric=.FALSE.
    a%if_zero_out_modes=.FALSE.
    a%verbose_timing =.FALSE.
    a%verbose_divergence=.FALSE.
    a%verbose_CFL=.FALSE.
    a%if_just_processing=.FALSE.
    a%if_post_proc_init=.FALSE.
    a%my_par_vv%verbose=.FALSE.
    a%my_par_pp%verbose=.FALSE.
    a%my_par_mass%verbose=.FALSE.
    a%my_par_H_p_phi%verbose=.FALSE.
    a%my_par_temperature%verbose=.FALSE.
    a%my_par_concentration%verbose=.FALSE.
    a%my_par_level_set%verbose=.FALSE.
    !===Reals
    a%LES_coeff1=0.d0
    a%LES_coeff2=0.d0
    a%LES_coeff3=0.d0
    a%LES_coeff4=0.d0
    a%LES_coeff1_mom=0.d0
    a%taux_precession=0.d0
    a%angle_s_pi=0.d0
    a%div_stab_in_ns=0.d0
    a%stab=0.d0
    a%gravity_coefficient=0.d0
    a%mag_force_coefficient=0.d0
    a%level_set_cmax=0.d0
    a%level_set_comp_factor=0.d0
    a%level_set_tanh_coeff_reconstruction=0.d0
    a%norm_ref=0.d0
    a%arpack_tolerance=0.d0
    a%stab_bdy_ns=0.d0
    a%h_min_distance=0.d0
    a%multiplier_for_h_min_distance=0.d0
    a%taylor_lambda = 1.d0
    a%penal_coeff_art_comp=1.d0
    a%MA=0.d0
    a%MB=1.d0
    a%rho_0percent=0.d0
    a%xi=0.d0
    a%pot_slope=1.d0
    a%faraday_cst=1.d0
    !===Integers
    a%vv_nb_dirichlet=0
    a%vv_nb_dirichlet_normal_velocity=0
    a%freq_restart = 10000000
    a%freq_en =  10000000
    a%freq_plot = 10000000
    a%type_fe_velocity = 2
    a%taylor_order = -1
    a%nb_dom_conc=0
    a%nb_dom_ns=0
    a%nb_dom_temp=0
    a%nb_dom_H=0
    a%nb_dom_phi=0
  END SUBROUTINE init
END MODULE my_data_module

MODULE input_data
  USE my_data_module
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE

CONTAINS

  SUBROUTINE read_my_data(data_fichier)
    USE chaine_caractere
    USE my_util
    IMPLICIT NONE
    LOGICAL                        :: test, test_sigma_min
    CHARACTER(len=200), INTENT(IN) :: data_fichier
    INTEGER                        :: k
    CHARACTER(len=8), DIMENSION(3) :: vel_component=(/'uradial','utheta ','uzaxis '/)

    !===Initialize data to zero and false by default===============================
    CALL inputs%init

    !===Open data file=============================================================
    OPEN(UNIT = 21, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')

    !===Location of mesh============================================================
    CALL read_until(21, '===Is mesh file formatted (true/false)?')
    READ(21,*) inputs%iformatted
    CALL read_until(21, '===Directory and name of mesh file')
    READ(21,*) inputs%directory, inputs%file_name

    !===Processor distribution======================================================
    CALL read_until(21, '===Number of processors in meridian section')
    READ (21,*) inputs%ndim(1)
    CALL read_until(21, '===Number of processors in Fourier space')
    READ (21,*) inputs%ndim(2)

    !===Fourier modes===============================================================
    CALL read_until(21, '===Number of Fourier modes')
    READ(21,*) inputs%m_max
    IF (MOD(inputs%m_max,inputs%ndim(2))/= 0) THEN
       CALL error_Petsc('BUG in read_my_data, MOD(nb_select_mode,nb_procs_F)/= 0')
    END IF
    CALL find_string(21,'===Select Fourier modes? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%select_mode
    ELSE
       inputs%select_mode = .FALSE.
    END IF
    IF (inputs%select_mode) THEN
       ALLOCATE(inputs%list_mode_lect(inputs%m_max))
       CALL read_until(21, '===List of Fourier modes (if select_mode=.TRUE.)')
       READ(21,*) inputs%list_mode_lect
    END IF

    !===Type of problem to be solved================================================
    CALL find_string(21,'===Problem type: (nst, mxw, mhd, fhd, mhs)', test)
    IF (test) THEN
       READ (21,*) inputs%type_pb
    ELSE
       CALL read_until(21, '===Problem type: (nst, mxw, mhd, fhd)')
       READ(21,*) inputs%type_pb
    END IF
    !   CALL read_until(21, '===Problem type: (nst, mxw, mhd, fhd, mhs)')
    !   READ(21,*) inputs%type_pb
    IF (inputs%type_pb/='nst' .AND. inputs%type_pb/='mhd' .AND. inputs%type_pb/='mxw'&
         .AND. inputs%type_pb/='fhd' .AND. inputs%type_pb/='mhs') THEN
       CALL error_Petsc('BUG in read_my_data, type_pb of probleme not yet defined')
    END IF

    !===Restarts====================================================================
    CALL read_until(21, '===Restart on velocity (true/false)')
    READ (21 ,*) inputs%irestart_u
    IF (inputs%type_pb=='mhd' .OR. inputs%type_pb=='mxw' .OR. inputs%type_pb=='fhd' &
         .OR. inputs%type_pb=='mhs') THEN
       CALL read_until(21, '===Restart on magnetic field (true/false)')
       READ (21 ,*) inputs%irestart_h
    ELSE
       inputs%irestart_h = .FALSE.
    END IF
    CALL find_string(21, '===Restart on temperature (true/false)', test)
    IF (test) THEN
       READ (21 ,*) inputs%irestart_T
    ELSE
       inputs%irestart_T = .FALSE.
    END IF
    CALL find_string(21, '===Restart on concentration (true/false)', test)
    IF (test) THEN
       READ (21 ,*) inputs%irestart_conc
    ELSE
       inputs%irestart_conc = .FALSE.
    END IF

    !===Mesh partitioning===========================================================
    CALL find_string(21, '===Do we read metis partition? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%if_read_partition
    ELSE
       inputs%if_read_partition = .FALSE.
    END IF
    IF (.NOT.inputs%if_read_partition .AND. (inputs%irestart_u .OR. inputs%irestart_h &
         .OR. inputs%irestart_T .OR. inputs%irestart_conc)) THEN
       call error_petsc('Possibility of bug: set "===Do we read metis partition? (true/false) "' // &
            'parameter to .true. when restarting a computation since number of procs, ' // &
            'or machine, or type_pb may have changed. Make sure your mesh_part file is the correct one.')
    END IF

    !===Time integration============================================================
    CALL read_until(21, '===Time step and number of time iterations')
    READ(21,*) inputs%dt, inputs%nb_iteration

    !===Check numerical stability===================================================
    CALL find_string(21, '===Check numerical stability (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%check_numerical_stability
    ELSE
       inputs%check_numerical_stability = .FALSE.
    END IF

    !===Periodicity=================================================================
    CALL find_string(21, '===How many pieces of periodic boundary?', test)
    IF (test) THEN
       READ (21,*) inputs%my_periodic%nb_periodic_pairs
    ELSE
       inputs%my_periodic%nb_periodic_pairs = 0
    END IF
    IF (inputs%my_periodic%nb_periodic_pairs.GE.1) THEN
       ALLOCATE(inputs%my_periodic%list_periodic(2,inputs%my_periodic%nb_periodic_pairs))
       ALLOCATE(inputs%my_periodic%vect_e(2,inputs%my_periodic%nb_periodic_pairs))
       CALL read_until(21, '===Indices of periodic boundaries and corresponding vectors')
       DO k = 1, inputs%my_periodic%nb_periodic_pairs
          READ(21,*) inputs%my_periodic%list_periodic(:,k), inputs%my_periodic%vect_e(:,k)
       END DO
    END IF

    !===Mesh symmetry===============================================================
    CALL find_string(21, '===Is the mesh symmetric (true/false)?', test)
    IF (test) THEN
       READ (21,*) inputs%is_mesh_symmetric
    ELSE
       inputs%is_mesh_symmetric = .FALSE.
    END IF


    !===Navier Stokes data==========================================================
    IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
         .OR. inputs%type_pb=='mhs' .OR. inputs%irestart_u) THEN

       !==========Navier Stokes with artifical compression or projection-correction schemes===============!
       CALL find_string(21, '===Solve Navier-Stokes with art comp scheme (true) or (false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_navier_stokes_art_comp
       ELSE
          inputs%if_navier_stokes_art_comp = .FALSE.
       END IF
       IF (inputs%if_navier_stokes_art_comp) THEN
          CALL find_string(21, '===Penalty coefficient for artifical compression', test)
          IF (test) THEN
             READ(21,*) inputs%penal_coeff_art_comp
          ELSE
             inputs%penal_coeff_art_comp=1.d0
          END IF
       END IF

       !==========Navier Stokes with u or m===============!
       CALL find_string(21, '===Solve Navier-Stokes with u (true) or m (false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_navier_stokes_with_u
       ELSE
          inputs%if_navier_stokes_with_u = .TRUE.
       END IF
       !==========Finite element type for velocity========!
       CALL find_string(21, '===Type of finite element for velocity field', test)
       IF (test) THEN
          READ(21,*) inputs%type_fe_velocity
       ELSE
          inputs%type_fe_velocity=2
       END IF
!!$       CALL find_string(21, '===Is tensor symmetric (true/false)?', test)
!!$       IF (test) THEN
!!$          READ(21,*) inputs%if_tensor_sym
!!$       ELSE
!!$          IF(inputs%if_navier_stokes_with_u) THEN
!!$             inputs%if_tensor_sym = .FALSE.
!!$          ELSE
!!$             inputs%if_tensor_sym = .TRUE.
!!$          END IF
!!$       END IF
       !===Tensor symmetric always used.
       inputs%if_tensor_sym = .TRUE.

       IF (.NOT.inputs%if_navier_stokes_with_u) THEN
          CALL find_string(21, '===Do we solve momentum with bdf2 (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_moment_bdf2
          ELSE
             inputs%if_moment_bdf2=.FALSE.
          END IF
          IF (inputs%if_moment_bdf2) THEN
             inputs%if_level_bdf2=.TRUE.
          ELSE
             CALL find_string(21, '===Do we solve level set with bdf2 (true/false)?', test)
             IF (test) THEN
                READ(21,*) inputs%if_level_bdf2
             ELSE
                inputs%if_level_bdf2=.FALSE.
             END IF
          END IF
       ELSE
          !===JLG July 20, 2019, p3 mesh
          CALL find_string(21, '===Do we use Taylor method to solve momentum (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_navier_stokes_with_taylor
             CALL read_until(21, '===Time accuracy of Taylor method (3 or 4)?')
             READ(21,*) inputs%taylor_order
             CALL read_until(21, '===Lambda parameter for Taylor method?')
             READ(21,*) inputs%taylor_lambda
          ELSE
             inputs%if_navier_stokes_with_taylor = .FALSE.
             inputs%taylor_order = -1
          END IF
          !===JLG July 20, 2019, p3 mesh
       END IF

       CALL read_until(21, '===Number of subdomains in Navier-Stokes mesh')
       READ(21,*) inputs%nb_dom_ns
       ALLOCATE(inputs%list_dom_ns(inputs%nb_dom_ns))
       CALL read_until(21, '===List of subdomains for Navier-Stokes mesh')
       READ(21,*) inputs%list_dom_ns

       !==========Dirichlet BCs for velocity==============!
       DO k = 1, 3
          CALL find_string(21, '===How many boundary pieces for Dirichlet BCs on '//TRIM(vel_component(k))//'?', test)
          IF (test) THEN
             CALL error_petsc( '===How many boundary pieces for Dirichlet BCs on '//TRIM(vel_component(k))//'? is disabled')
             !READ(21,*) inputs%vv_nb_dirichlet_sides(k)
          ELSE
             inputs%vv_nb_dirichlet_sides(k) = 0
          END IF
          IF (inputs%vv_nb_dirichlet_sides(k)>0) THEN
             ALLOCATE(inputs%vv_list_dirichlet_sides(k)%DIL(inputs%vv_nb_dirichlet_sides(k)))
             CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on '//TRIM(vel_component(k)))
             READ(21,*) inputs%vv_list_dirichlet_sides(k)%DIL
          ELSE
             ALLOCATE(inputs%vv_list_dirichlet_sides(k)%DIL(0))
          END IF
       END DO
       CALL find_string(21, '===How many boundary pieces for full Dirichlet BCs on velocity?', test)
       IF (test) THEN
          READ(21,*) inputs%vv_nb_dirichlet
       ELSE
          inputs%vv_nb_dirichlet=0
       END IF
       IF (inputs%vv_nb_dirichlet>0) THEN
          DO k = 1, 3
             ALLOCATE(inputs%vv_list_dirichlet_sides(k)%DIL(inputs%vv_nb_dirichlet))
             CALL read_until(21, '===List of boundary pieces for full Dirichlet BCs on velocity')
             READ(21,*) inputs%vv_list_dirichlet_sides(k)%DIL
          END DO
       ELSE
          DO k = 1, 3
             ALLOCATE(inputs%vv_list_dirichlet_sides(k)%DIL(0))
          END DO
       END IF

       !===Homogeneous normal velocity====================!
       CALL find_string(21, '===How many boundary pieces for homogeneous normal velocity?', test)
       IF (test) THEN
          READ(21,*) inputs%vv_nb_dirichlet_normal_velocity
       ELSE
          inputs%vv_nb_dirichlet_normal_velocity=0
       END IF
       IF (inputs%vv_nb_dirichlet_normal_velocity>0) THEN
          ALLOCATE(inputs%vv_list_dirichlet_normal_velocity_sides(inputs%vv_nb_dirichlet_normal_velocity))
          CALL read_until(21, '===List of boundary pieces for homogeneous normal velocity')
          READ(21,*) inputs%vv_list_dirichlet_normal_velocity_sides
          CALL find_string(21, '===stab_bdy_ns', test)
          IF (test) THEN
             READ (21,*) inputs%stab_bdy_ns
          ELSE
             inputs%stab_bdy_ns=1.d0
          END IF
          DO k = 1, inputs%vv_nb_dirichlet_normal_velocity
             IF (MINVAL(ABS(inputs%vv_list_dirichlet_normal_velocity_sides(k) &
                  - inputs%vv_list_dirichlet_sides(1)%DIL))==0) THEN
                CALL error_petsc('Boundary conditions are mixed up')
             END IF
          END DO
       END IF

       !===Stress boundary conditions=====================! !Disabled
       CALL find_string(21, '===Stress boundary conditions? (true/false)', test)
       IF (test) THEN
          CALL error_petsc('===Stress boundary conditions Does not work any more')
       END IF

       !==========Dirichlet BCs for pressure==============! !Disabled
       !CALL find_string(21, '===How many boundary pieces for Dirichlet BCs on pressure?', test)
       test=.FALSE.
       IF (test) THEN
          READ(21,*) inputs%pp_nb_dirichlet_sides
       ELSE
          inputs%pp_nb_dirichlet_sides = 0
       END IF
       IF (inputs%pp_nb_dirichlet_sides>0) THEN
          ALLOCATE(inputs%pp_list_dirichlet_sides(inputs%pp_nb_dirichlet_sides))
          CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on pressure')
          READ(21,*) inputs%pp_list_dirichlet_sides
       ELSE
          ALLOCATE(inputs%pp_list_dirichlet_sides(0))
       END IF

       !==========Reynolds number=========================!
       CALL find_string(21, '===Reynolds number',test)
       IF (test) THEN
          READ(21,*) inputs%Re
       ELSE
          CALL read_until(21, '===Kinematic viscosity')
          READ(21,*) inputs%Re
          inputs%Re = 1.d0 / inputs%Re
       END IF


       !==========Lorentz Force Coefficient===============!
       CALL find_string(21, '===Coefficient for Lorentz force', test)
       IF (test) THEN
          READ(21,*) inputs%coeff_lorentz
       ELSE
          inputs%coeff_lorentz=1.d0
       END IF

       !==========Variable viscosity======================!
       CALL find_string(21, '===Variable viscosity (true/false)?',test)
       IF (test) THEN
          READ(21,*) inputs%if_variable_visco
       ELSE
          inputs%if_variable_visco = .FALSE.
       END IF

       !==========DIV penalty=============================!
       IF (inputs%if_navier_stokes_with_u) THEN
          CALL find_string(21, '===Coefficient for penalty of divergence in NS?', test)
          IF (test) THEN
             READ(21,*) inputs%div_stab_in_ns
          ELSE
             inputs%div_stab_in_ns = 0.d0
          END IF
       END IF

       !==========Precession==============================!
       CALL find_string(21, '===Is there a precession term (true/false)?', test)
       IF (test) THEN
          READ(21,*) inputs%precession
       ELSE
          inputs%precession = .FALSE.
       END IF
       IF (inputs%precession) THEN
          CALL read_until(21, '===Precession rate')
          READ(21,*) inputs%taux_precession
          CALL read_until(21, '===Precession angle over pi')
          READ(21,*) inputs%angle_s_pi
       ELSE
          inputs%taux_precession = 0.d0
          inputs%angle_s_pi = 0.d0
       END IF

       !==========NS penalty==============================!
       CALL find_string(21, '===Use penalty in NS domain (true/false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_ns_penalty
       ELSE
          inputs%if_ns_penalty = .FALSE.
       END IF
       IF (inputs%if_ns_penalty) THEN
          CALL find_string(21, '===Use nonzero velocity in solids (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_impose_vel_in_solids
          ELSE
             inputs%if_impose_vel_in_solids = .FALSE.
          END IF
          CALL find_string(21, '===Compute z momentum (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_compute_momentum_pseudo_force
          ELSE
             inputs%if_compute_momentum_pseudo_force = .FALSE.
          END IF
       ELSE
          inputs%if_impose_vel_in_solids = .FALSE.
          inputs%if_compute_momentum_pseudo_force = .FALSE.
       END IF

       !==========Solver parameters for velocity==========!
       CALL read_until(21, '===Maximum number of iterations for velocity solver')
       READ(21,*) inputs%my_par_vv%it_max
       CALL read_until(21, '===Relative tolerance for velocity solver')
       READ(21,*) inputs%my_par_vv%rel_tol
       CALL read_until(21, '===Absolute tolerance for velocity solver')
       READ(21,*) inputs%my_par_vv%abs_tol
       CALL find_string(21, '===Velocity solver verbose? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%my_par_vv%verbose
       END IF
       CALL read_until(21, '===Solver type for velocity (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_vv%solver
       CALL read_until(21, '===Preconditionner type for velocity solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_vv%precond

       !==========Solver parameters for pressure==========!
       CALL read_until(21, '===Maximum number of iterations for pressure solver')
       READ(21,*) inputs%my_par_pp%it_max
       CALL read_until(21, '===Relative tolerance for pressure solver')
       READ(21,*) inputs%my_par_pp%rel_tol
       CALL read_until(21, '===Absolute tolerance for pressure solver')
       READ(21,*) inputs%my_par_pp%abs_tol
       CALL find_string(21, '===Pressure solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_pp%verbose
       END IF
       CALL read_until(21, '===Solver type for pressure (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_pp%solver
       CALL read_until(21, '===Preconditionner type for pressure solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_pp%precond

       !==========Solver parameters for mass matrix=======!
       CALL read_until(21, '===Maximum number of iterations for mass matrix solver')
       READ(21,*) inputs%my_par_mass%it_max
       CALL read_until(21, '===Relative tolerance for mass matrix solver')
       READ(21,*) inputs%my_par_mass%rel_tol
       CALL read_until(21, '===Absolute tolerance for mass matrix solver')
       READ(21,*) inputs%my_par_mass%abs_tol
       CALL find_string(21, '===Mass matrix solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_mass%verbose
       END IF
       CALL read_until(21, '===Solver type for mass matrix (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_mass%solver
       CALL read_until(21, '===Preconditionner type for mass matrix solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_mass%precond

       !==========LES coefficients========================!
       CALL find_string(21, '===Use LES? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%LES
       ELSE
          inputs%LES = .FALSE.
       END IF

       IF (inputs%LES) THEN
          CALL find_string(21, '===Coefficients for LES', test)
          IF (test) THEN
             CALL error_petsc('===Coefficients for LES is disabled')
          END IF
          !===LES_coeff4 is fixed and coeff3 is disabled
          inputs%LES_coeff4 = 0.125d0
          inputs%LES_coeff3 = 0.d0
          !===LES_coeff2 in [0.15,1.0] for viscosity entropy
          !===LES_coeff2 is equal to 1.d10 for test with first order viscosity
          CALL read_until(21, '===Coefficient multiplying residual')
          READ(21,*) inputs%LES_coeff2
          !===LES_coeff1 is equal to LES_coeff4*MAX(velocity)=0.125d0*MAX(velocity)
          CALL find_string(21, '===Coefficient for explicit LES', test)
          IF (test) THEN
             READ(21,*) inputs%LES_coeff1
          ELSE
             inputs%LES_coeff1 = 0.d0
          END IF
          !===Check if LES in restart file
          IF (inputs%irestart_u) THEN
             CALL find_string(21, '===Restart on LES (true/false)',test)
             IF (test) THEN
                READ(21,*) inputs%irestart_LES
             ELSE
                inputs%irestart_LES = .FALSE.
             END IF
          ELSE
             inputs%irestart_LES = .FALSE.
          END IF
       ELSE
          inputs%LES_coeff1 = 0.d0
          inputs%LES_coeff2 = 0.d0
          inputs%LES_coeff3 = 0.d0
          inputs%LES_coeff4 = 0.d0
          inputs%irestart_LES = .FALSE.
       END IF
    END IF
    IF (.NOT.inputs%if_navier_stokes_with_u) THEN
       CALL find_string(21, '===Use LES in momentum? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%if_LES_in_momentum
       ELSE
          inputs%if_LES_in_momentum=.TRUE.
       END IF
       IF (inputs%if_LES_in_momentum) THEN
          inputs%LES_coeff1_mom=inputs%LES_coeff1
       ELSE
          inputs%LES_coeff1_mom=0.d0
       END IF
    END IF

    !===Maxwell data================================================================
    IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd'&
         .OR. inputs%type_pb=='mhs') THEN
       !==========Maxwell with H or B=====================!
       CALL find_string(21, '===Solve Maxwell with H (true) or B (false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_maxwell_with_H
       ELSE
          inputs%if_maxwell_with_H = .FALSE.
       END IF

       !==========Quasi-static approximation==============!
       CALL find_string(21, '===Quasi-static approximation (true) or (false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_quasi_static_approx
       ELSE
          inputs%if_quasi_static_approx = .FALSE.
       END IF

       !==========Steady current in fhd===================!
       IF (inputs%type_pb=='fhd') THEN
          CALL find_string(21, '===Steady current for fhd (true or false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_steady_current_fhd
          ELSE
             inputs%if_steady_current_fhd = .FALSE.
          END IF
       END IF

       !==========Subdomains for H========================!
       CALL read_until(21, '===Number of subdomains in magnetic field (H) mesh')
       READ(21,*) inputs%nb_dom_H  ! number of sub_domains for H
       ALLOCATE(inputs%list_dom_H(inputs%nb_dom_H))
       CALL read_until(21, '===List of subdomains for magnetic field (H) mesh')
       READ(21,*) inputs%list_dom_H

       !==========Interfaces H/H==========================!
       CALL read_until(21, '===Number of interfaces in H mesh')
       READ(21,*) inputs%nb_inter_mu
       ALLOCATE(inputs%list_inter_mu(inputs%nb_inter_mu))
       IF (inputs%nb_inter_mu>0) THEN
          CALL read_until(21, '===List of interfaces in H mesh')
          READ(21,*) inputs%list_inter_mu
       END IF

       !==========Interfaces jump rot h==========================!
       CALL find_string(21, '===How many boundary pieces for jump BCs on rot H?', test)
       IF (test) THEN
          READ(21,*) inputs%rot_h_nb_jump_sides
       ELSE
          inputs%rot_h_nb_jump_sides = 0
       END IF
       IF (test .AND. (inputs%rot_h_nb_jump_sides>0)) THEN
          ALLOCATE(inputs%list_inter_rot_h_jump(inputs%rot_h_nb_jump_sides))
          CALL read_until(21, '===List of boundary pieces for jump BCs on rot H')
          READ(21,*) inputs%list_inter_rot_h_jump
       ELSE
          ALLOCATE(inputs%list_inter_rot_h_jump(0))
       END IF

       !==========Dirichlet BCs for H=====================!
       CALL read_until(21, '===Number of Dirichlet sides for Hxn')
       READ(21,*) inputs%nb_dirichlet_sides_H
       ALLOCATE(inputs%list_dirichlet_sides_H(inputs%nb_dirichlet_sides_H))
       IF (inputs%nb_dirichlet_sides_H>0) THEN
          CALL read_until(21, '===List of Dirichlet sides for Hxn')
          READ(21,*) inputs%list_dirichlet_sides_H
       END IF

       !==========Permeability for H======================!
       CALL find_string(21, '===Is permeability defined analytically (true/false)?', test)
       IF (test) THEN
          READ(21,*) inputs%analytical_permeability
       ELSE
          inputs%analytical_permeability = .FALSE.
       END IF
       IF (.NOT.inputs%analytical_permeability) THEN
          CALL read_until(21, '===Permeability in the conductive part (1:nb_dom_H)')
          ALLOCATE(inputs%mu_H(inputs%nb_dom_H))
          READ(21,*) inputs%mu_H
       END IF
       inputs%if_permeability_variable_in_theta = .FALSE.
       IF (inputs%analytical_permeability) THEN
          CALL find_string(21, '===Is permeability variable in theta (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_permeability_variable_in_theta
          END IF
       END IF

       !==========FEM or Gaussian integration for mu_bar==!
       inputs%if_use_fem_integration_for_mu_bar = .TRUE.
       IF (inputs%analytical_permeability) THEN
          CALL find_string(21, '===Use FEM Interpolation for magnetic permeability (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_use_fem_integration_for_mu_bar
          END IF
       END IF

       !==========Conductivity for H======================!
       CALL read_until(21, '===Conductivity in the conductive part (1:nb_dom_H)')
       ALLOCATE(inputs%sigma(inputs%nb_dom_H))
       READ(21,*) inputs%sigma

       !==========Minimum of conductivity=================!
       CALL find_string(21,'===Minimum of conductivity in the whole domain', test)
       IF (test) THEN
          READ(21,*) inputs%sigma_min
       ELSE
          test_sigma_min = .TRUE.
          inputs%sigma_min=MINVAL(inputs%sigma) !JLG (Feb 23/2017)
       END IF

       !==========Minimum of permeability=================!
       CALL find_string(21,'===Minimum of permeability in the whole domain', test)
       IF (test) THEN
          READ(21,*) inputs%mu_min
       ELSE
          IF (.NOT.inputs%analytical_permeability) THEN ! JLG (FEB 23, 2017), begin
             inputs%mu_min = MINVAL(inputs%mu_H)
          ELSE
             inputs%mu_min=1.d0
          END IF ! JLG (FEB 23, 2017), end
       END IF

       !==========Finite element type=====================!
       CALL read_until(21, '===Type of finite element for magnetic field')
       READ(21,*) inputs%type_fe_H

       !==========Magnetic Reynolds number================!
       CALL read_until(21, '===Magnetic Reynolds number')
       READ(21,*) inputs%Rem

       !==========Stabilization parameters================!
       CALL read_until(21, '===Stabilization coefficient (divergence)')
       READ(21,*) inputs%stab(1)
       IF (inputs%nb_inter_mu>0 .OR. inputs%nb_dirichlet_sides_H>0) THEN
          CALL read_until(21, '===Stabilization coefficient for Dirichlet H and/or interface H/H')
          READ(21,*) inputs%stab(3)
       ELSE
          inputs%stab(3) = 0.d0
       END IF
       !==========Stab coefficient for interface jump on H========================!
       IF (inputs%nb_inter_mu>0) THEN
          !CALL read_until(21, '===Stabilization coefficient for interface H/H')
          CALL find_string(21, '===Stabilization coefficient for interface H/H', test)
          IF (test) THEN
             READ(21,*) inputs%stab_jump_h
          ELSE
             inputs%stab_jump_h=1.d0
          END IF
       ELSE
          inputs%stab_jump_h = 0.d0
       END IF
       !==========Subdomains for phi======================!
       CALL find_string(21, '===Number of subdomains in magnetic potential (phi) mesh', test)
       IF (test) THEN
          READ (21,*) inputs%nb_dom_phi
       ELSE
          inputs%nb_dom_phi = 0
       END IF
       ALLOCATE(inputs%list_dom_phi(inputs%nb_dom_phi))
       IF (inputs%nb_dom_phi>0) THEN
          !==========List of subdomains for phi======================!
          CALL read_until(21, '===List of subdomains for magnetic potential (phi) mesh')
          READ(21,*) inputs%list_dom_phi

          !==========Dirichlet BCs on phi====================!
          CALL read_until(21, '===How many boundary pieces for Dirichlet BCs on phi?')
          READ(21,*) inputs%phi_nb_dirichlet_sides
          IF (inputs%phi_nb_dirichlet_sides>0) THEN
             ALLOCATE(inputs%phi_list_dirichlet_sides(inputs%phi_nb_dirichlet_sides))
             CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on phi')
             READ(21,*) inputs%phi_list_dirichlet_sides
          ELSE
             ALLOCATE(inputs%phi_list_dirichlet_sides(0))
          END IF

          !==========H/phi interface=========================!
          CALL read_until(21, '===Number of interfaces between H and phi')
          READ(21,*) inputs%nb_inter  ! number of interfaces between H and phi
          ALLOCATE(inputs%list_inter_H_phi(inputs%nb_inter))
          CALL read_until(21, '===List of interfaces between H and phi')
          IF (inputs%nb_inter>0) THEN
             READ(21,*) inputs%list_inter_H_phi
          END IF

          !==========Permeability in vacuum==================!
          CALL read_until(21, '===Permeability in vacuum')
          READ(21,*) inputs%mu_phi

          !==========Finite element type=====================!
          CALL read_until(21, '===Type of finite element for scalar potential')
          READ(21,*) inputs%type_fe_phi

          !==========Stabilization parameters================!
          CALL read_until(21, '===Stabilization coefficient (interface H/phi)')
          READ(21,*) inputs%stab(2)
       END IF

       !==========Solver parameters=======================!
       CALL read_until(21, '===Maximum number of iterations for Maxwell solver')
       READ(21,*) inputs%my_par_H_p_phi%it_max
       CALL read_until(21, '===Relative tolerance for Maxwell solver')
       READ(21,*) inputs%my_par_H_p_phi%rel_tol
       CALL read_until(21, '===Absolute tolerance for Maxwell solver')
       READ(21,*) inputs%my_par_H_p_phi%abs_tol
       CALL find_string(21, '===Maxwell solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_H_p_phi%verbose
       END IF
       CALL read_until(21, '===Solver type for Maxwell (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_H_p_phi%solver
       CALL read_until(21, '===Preconditionner type for Maxwell solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_H_p_phi%precond
    END IF


    !===Data temperature==============================================================
    CALL find_string(21, '===Is there a temperature field?', test)
    IF (test) THEN
       READ (21,*) inputs%if_temperature
    ELSE
       inputs%if_temperature = .FALSE.
    END IF
    IF (inputs%if_temperature) THEN
       !==========Temperature with T or e=================!
       CALL find_string(21, '===Solve Temperature with T (true) or e (false)?', test)
       IF (test) THEN
          READ(21,*) inputs%if_temperature_with_T
       ELSE
          inputs%if_temperature_with_T = .TRUE.
       END IF

       IF (.NOT. inputs%if_temperature_with_T) THEN
          CALL find_string(21, '===Do we solve heat equation with bdf2 (true/false)?', test)
          IF (test) THEN
             READ(21,*) inputs%if_temp_bdf2
          ELSE
             inputs%if_temp_bdf2=.FALSE.
          END IF
          IF (inputs%if_temp_bdf2) THEN
             inputs%if_level_bdf2=.TRUE.
          END IF
       END IF
       !==========Gravity coefficient for temperature=====!
       CALL find_string(21, '===Non-dimensional gravity coefficient', test)
       IF (test) THEN
          READ (21,*) inputs%gravity_coefficient
       ELSE
          inputs%gravity_coefficient = 0.d0
       END IF
       !==========Magnetic force coefficient==============!
       CALL find_string(21, '===Helmholtz force for ferrohydrodynamics? (true/false)', test)
       IF (test) THEN
          READ (21,*) inputs%if_helmholtz_force
       END IF
       CALL find_string(21, '===Non-dimensional magnetic force coefficient for ferrohydrodynamics', test)
       IF (test) THEN
          READ (21,*) inputs%mag_force_coefficient
       END IF
       !==========Subdomains for temp=====================!
       CALL read_until(21, '===Number of subdomains in temperature mesh')
       READ(21,*) inputs%nb_dom_temp  ! number of sub_domains for temp
       ALLOCATE(inputs%list_dom_temp(inputs%nb_dom_temp))
       CALL read_until(21, '===List of subdomains for temperature mesh')
       READ(21,*) inputs%list_dom_temp
       ALLOCATE(inputs%vol_heat_capacity(inputs%nb_dom_temp)) ! Two choices for the user, indicate 1) the heat capacity and the conductivity (necessary if non uniform heat capacity) or 2) the diffusivity only
       ALLOCATE(inputs%temperature_diffusivity(inputs%nb_dom_temp)) ! In both cases, temperature_diffusivity is used, it contains the conductivities in 1) and the diffusivities in 2)
       CALL find_string(21, '===Density (1:nb_dom_temp)', test)
       IF (test) THEN ! Case 1)
          ALLOCATE(inputs%density(inputs%nb_dom_temp))
          ALLOCATE(inputs%heat_capacity(inputs%nb_dom_temp))
          READ(21,*) inputs%density
          CALL read_until(21, '===Heat capacity (1:nb_dom_temp)')
          READ(21,*) inputs%heat_capacity
          inputs%vol_heat_capacity = inputs%density * inputs%heat_capacity
          DEALLOCATE (inputs%density)
          DEALLOCATE (inputs%heat_capacity)
          CALL read_until(21, '===Thermal conductivity (1:nb_dom_temp)')
          READ(21,*) inputs%temperature_diffusivity
       ELSE ! Case 1) bis
          CALL find_string(21, '===Volumetric heat capacity (1:nb_dom_temp)', test)
          IF (test) THEN
             READ(21,*) inputs%vol_heat_capacity
             CALL read_until(21, '===Thermal conductivity (1:nb_dom_temp)')
             READ(21,*) inputs%temperature_diffusivity
          ELSE ! Case 2)
             inputs%vol_heat_capacity = 1.d0 ! Heat capacity is equalized to one so that it does not impact the calculus
             CALL read_until(21, '===Diffusivity coefficient for temperature (1:nb_dom_temp)')
             READ(21,*) inputs%temperature_diffusivity
          END IF
       END IF
       !==========Dirichlet BCs on temperature============!
       CALL read_until(21, '===How many boundary pieces for Dirichlet BCs on temperature?')
       READ(21,*) inputs%temperature_nb_dirichlet_sides
       IF (inputs%temperature_nb_dirichlet_sides>0) THEN
          ALLOCATE(inputs%temperature_list_dirichlet_sides(inputs%temperature_nb_dirichlet_sides))
          CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on temperature')
          READ(21,*) inputs%temperature_list_dirichlet_sides
       ELSE
          ALLOCATE(inputs%temperature_list_dirichlet_sides(0))
       END IF
       !==========Robin BCs on temperature================!
       CALL find_string(21, '===How many boundary pieces for Robin BCs on temperature?', test)
       IF (test) THEN
          READ(21,*) inputs%temperature_nb_robin_sides
       ELSE
          inputs%temperature_nb_robin_sides = 0
       END IF
       IF (test .AND. (inputs%temperature_nb_robin_sides>0)) THEN
          ALLOCATE(inputs%temperature_list_robin_sides(inputs%temperature_nb_robin_sides))
          CALL read_until(21, '===List of boundary pieces for Robin BCs on temperature')
          READ(21,*) inputs%temperature_list_robin_sides
          ALLOCATE(inputs%convection_coeff(inputs%temperature_nb_robin_sides))
          CALL read_until(21, '===Convection heat transfert coefficient (1:temperature_nb_robin_sides)')
          READ(21,*) inputs%convection_coeff
          ALLOCATE(inputs%exterior_temperature(inputs%temperature_nb_robin_sides))
          CALL read_until(21, '===Exterior temperature (1:temperature_nb_robin_sides)')
          READ(21,*) inputs%exterior_temperature
       ELSE
          ALLOCATE(inputs%temperature_list_robin_sides(0))
       END IF
       !==========Inhomogeneous Neumann BCs on temperature================!
       ! SB 29/04/2022 Additional parameter for inhomogeneous Neumann BC for temp related to H
       IF ((inputs%type_pb=='mhs') .AND. (inputs%rot_h_nb_jump_sides>0)) THEN
          ALLOCATE(inputs%temperature_list_neumann_sides(inputs%rot_h_nb_jump_sides))
          inputs%temperature_list_neumann_sides = inputs%list_inter_rot_h_jump
       END IF

       !==========Interfaces between vel and temp=========!
       CALL find_string(21, '===Number of interfaces between velocity and temperature only domains (for nst applications)', test)
       IF (test) THEN
          READ(21,*) inputs%nb_inter_v_T
          ALLOCATE(inputs%list_inter_v_T(inputs%nb_inter_v_T))
          CALL read_until(21, '===List of interfaces between velocity and temperature only domains (for nst applications)')
          READ(21,*) inputs%list_inter_v_T
       ELSE
          ALLOCATE(inputs%list_inter_v_T(0))
       END IF
       !==========Solver parameters=======================!
       CALL read_until(21, '===Maximum number of iterations for temperature solver')
       READ(21,*) inputs%my_par_temperature%it_max
       CALL read_until(21, '===Relative tolerance for temperature solver')
       READ(21,*) inputs%my_par_temperature%rel_tol
       CALL read_until(21, '===Absolute tolerance for temperature solver')
       READ(21,*) inputs%my_par_temperature%abs_tol
       CALL find_string(21, '===Temperature solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_temperature%verbose
       END IF
       CALL read_until(21, '===Solver type for temperature (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_temperature%solver
       CALL read_until(21, '===Preconditionner type for temperature solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_temperature%precond
    END IF



    !===Data concentration==============================================================
    CALL find_string(21, '===Is there a concentration field?', test)
    IF (test) THEN
       READ (21,*) inputs%if_concentration
    ELSE
       inputs%if_concentration = .FALSE.
    END IF
    IF (inputs%if_concentration) THEN
       !==========Subdomains for conc=====================!
       CALL read_until(21, '===Number of subdomains in concentration mesh')
       READ(21,*) inputs%nb_dom_conc  ! number of sub_domains for conc
       ALLOCATE(inputs%list_dom_conc(inputs%nb_dom_conc))
       CALL read_until(21, '===List of subdomains for concentration mesh')
       READ(21,*) inputs%list_dom_conc
       ALLOCATE(inputs%concentration_diffusivity(inputs%nb_dom_conc)) ! In both cases, concentration_diffusivity is used, it contains the conductivities in 1) and the diffusivities in 2)
       CALL read_until(21, '===Diffusivity coefficient for concentration (1:nb_dom_conc)')
       READ(21,*) inputs%concentration_diffusivity
       !==========Dirichlet BCs on concentration============!
       CALL read_until(21, '===How many boundary pieces for Dirichlet BCs on concentration?')
       READ(21,*) inputs%concentration_nb_dirichlet_sides
       IF (inputs%concentration_nb_dirichlet_sides>0) THEN
          ALLOCATE(inputs%concentration_list_dirichlet_sides(inputs%concentration_nb_dirichlet_sides))
          CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on concentration')
          READ(21,*) inputs%concentration_list_dirichlet_sides
       ELSE
          ALLOCATE(inputs%concentration_list_dirichlet_sides(0))
       END IF
       !==========Robin BCs on concentration================!
       CALL find_string(21, '===How many boundary pieces for Robin BCs on concentration?', test)
       IF (test) THEN
          READ(21,*) inputs%concentration_nb_robin_sides
       ELSE
          inputs%concentration_nb_robin_sides = 0
       END IF
       IF (test .AND. (inputs%concentration_nb_robin_sides>0)) THEN
          ALLOCATE(inputs%concentration_list_robin_sides(inputs%concentration_nb_robin_sides))
          CALL read_until(21, '===List of boundary pieces for Robin BCs on concentration')
          READ(21,*) inputs%concentration_list_robin_sides
          ALLOCATE(inputs%convection_coeff_conc_lhs(inputs%concentration_nb_robin_sides))
          CALL read_until(21, '===Convection heat transfert coefficient of lhs (1:concentration_nb_robin_sides)')
          READ(21,*) inputs%convection_coeff_conc_lhs
          ALLOCATE(inputs%convection_coeff_conc_rhs(inputs%concentration_nb_robin_sides))
          CALL read_until(21, '===Convection heat transfert coefficient of rhs (1:concentration_nb_robin_sides)')
          READ(21,*) inputs%convection_coeff_conc_rhs
          ALLOCATE(inputs%exterior_concentration(inputs%concentration_nb_robin_sides))
          CALL read_until(21, '===Exterior concentration (1:concentration_nb_robin_sides)')
          READ(21,*) inputs%exterior_concentration
       ELSE
          ALLOCATE(inputs%concentration_list_robin_sides(0))
       END IF

       !==========Neumann BCs on concentration================!
       CALL find_string(21, '===How many boundary pieces for Neumann BCs on concentration?', test)
       IF (test) THEN
          READ(21,*) inputs%concentration_nb_neumann_sides
       ELSE
          inputs%concentration_nb_neumann_sides = 0
       END IF
       IF (test .AND. (inputs%concentration_nb_neumann_sides>0)) THEN
          ALLOCATE(inputs%concentration_list_neumann_sides(inputs%concentration_nb_neumann_sides))
          CALL read_until(21, '===List of boundary pieces for Neumann BCs on concentration')
          READ(21,*) inputs%concentration_list_neumann_sides
       ELSE
          ALLOCATE(inputs%concentration_list_neumann_sides(0))
       END IF
       !==========Interfaces between vel and conc=========!
       CALL find_string(21, '===Number of interfaces between velocity and concentration only domains (for nst applications)', test)
       IF (test) THEN
          READ(21,*) inputs%nb_inter_c_v
          ALLOCATE(inputs%list_inter_c_v(inputs%nb_inter_c_v))
          CALL read_until(21, '===List of interfaces between velocity and concentration only domains (for nst applications)')
          READ(21,*) inputs%list_inter_c_v
       ELSE
          ALLOCATE(inputs%list_inter_c_v(0))
       END IF
       !==========Solver parameters=======================!
       CALL read_until(21, '===Maximum number of iterations for concentration solver')
       READ(21,*) inputs%my_par_concentration%it_max
       CALL read_until(21, '===Relative tolerance for concentration solver')
       READ(21,*) inputs%my_par_concentration%rel_tol
       CALL read_until(21, '===Absolute tolerance for concentration solver')
       READ(21,*) inputs%my_par_concentration%abs_tol
       CALL find_string(21, '===Temperature solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_concentration%verbose
       END IF
       CALL read_until(21, '===Solver type for concentration (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_concentration%solver
       CALL read_until(21, '===Preconditionner type for concentration solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_concentration%precond
    END IF

    !===Data For mhs problem or with concentration======================================
    CALL find_string(21, '===Is there a molar fraction law as a function of density? (true/false)', test)
    IF (test) THEN
       CALL read_until(21, '===Mass of A for molar fraction law')
       READ(21,*) inputs%MA
       CALL read_until(21, '===Mass of B for molar fraction law')
       READ(21,*) inputs%MB
       CALL read_until(21, '===Reference density of pure B')
       READ(21,*) inputs%rho_0percent
       CALL read_until(21, '===Slope of density law')
       READ(21,*) inputs%xi
       CALL read_until(21, '===Slope of potential law')
       READ(21,*) inputs%pot_slope
       CALL read_until(21, '===Faraday constant')
       READ(21,*) inputs%faraday_cst
       CALL find_string(21, '===Is there a coupling between H and molar fraction? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%if_coupling_H_x
       ELSE
          inputs%if_coupling_H_x=.FALSE.
       END IF
    ELSE
       inputs%if_coupling_H_x = .FALSE.
       inputs%MA = 0.d0
       inputs%MB = 1.d0
       inputs%rho_0percent = 0.d0
       inputs%xi = 0.d0
       inputs%pot_slope = 1.d0
       inputs%faraday_cst = 1.d0
    END IF
    CALL find_string(21, '===Is there an analytical coupling for mhs problems? (true/false)', test)
    IF (test) THEN
       READ(21,*) inputs%if_coupling_analytical
    ELSE
       inputs%if_coupling_analytical = .FALSE.
    END IF

    !===Data Level set================================================================
    CALL find_string(21, '===Is there a level set?', test)
    IF (test) THEN
       READ (21,*) inputs%if_level_set
    ELSE
       inputs%if_level_set = .FALSE.
    END IF
    IF (inputs%if_level_set) THEN
       !==========Fix level set to test code==============!
       CALL find_string(21, '===Do we fix level set? (true/false)', test)
       IF (test) THEN
          READ (21, *) inputs%if_level_set_fixed
       ELSE
          inputs%if_level_set_fixed = .FALSE.
       END IF
       !==========Number of fluids========================!
       CALL read_until(21, '===How many fluids?')
       READ(21,*) inputs%nb_fluid
       !==========Level_set multiplier for h_min==========!
       CALL read_until(21, '===multiplier for h_min for level set')
       READ(21,*) inputs%multiplier_for_h_min_distance
       !==========Level_set c_max=========================!
       !===Disabled functionality
       inputs%level_set_cmax=0.d0
       !==========Level_set compression factor============!
       CALL read_until(21, '===Compression factor for level set')
       READ(21,*) inputs%level_set_comp_factor
       !==========Densities of fluids=====================!
       ALLOCATE(inputs%density_fluid( inputs%nb_fluid))
       CALL read_until(21, '===Density of fluid 0, fluid 1, ...')
       READ(21,*) inputs%density_fluid
       !==========Dynamic viscosities of fluids===========!
       ALLOCATE(inputs%dyna_visc_fluid( inputs%nb_fluid))
       CALL read_until(21, '===Dynamic viscosity of fluid 0, fluid 1, ...')
       READ(21,*) inputs%dyna_visc_fluid
       !==========Conductivities of fluids================!
       IF (inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' .OR. inputs%type_pb=='mhs') THEN
          CALL find_string(21, '===Conductivity of fluid 0, fluid 1, ...', test)
          IF (test) THEN
             inputs%variation_sigma_fluid =.TRUE.
             ALLOCATE(inputs%sigma_fluid(inputs%nb_fluid))
             READ(21,*) inputs%sigma_fluid
             !===sigma_min has not been been read=========!
             IF(test_sigma_min) THEN
                inputs%sigma_min=MINVAL(inputs%sigma_fluid) !===JLG May 5,
             END IF
          ELSE
             inputs%variation_sigma_fluid =.FALSE.
          END IF
       END IF
       !==========Temperature parameters of fluids========!
       IF(.NOT.inputs%if_temperature_with_T) THEN
          ALLOCATE(inputs%heat_capacity_fluid(inputs%nb_fluid))
          ALLOCATE(inputs%heat_diffu_fluid(inputs%nb_fluid))
          ALLOCATE(inputs%heat_grav_fluid(inputs%nb_fluid))
          CALL find_string(21, '===Heat capacity of fluid 0, fluid 1, ...', test)
          IF (test) THEN
             inputs%variation_temp_param_fluid=.TRUE.
             READ(21,*) inputs%heat_capacity_fluid
          ELSE
             inputs%variation_temp_param_fluid=.FALSE.
             inputs%heat_capacity_fluid=1.d0
          END IF
          CALL find_string(21, '===Thermal conductivity of fluid 0, fluid 1, ...', test)
          IF (test) THEN
             IF(.NOT. inputs%variation_temp_param_fluid) THEN
                CALL error_petsc('BUG in read_my_data: if define Thermal Conductivity in fluid,' // &
                     ' one also needs to define Heat capacity in fluid')
             ELSE
                READ(21,*) inputs%heat_diffu_fluid
             END IF
          ELSE
             inputs%heat_diffu_fluid=1.d0
          END IF
          CALL find_string(21, '===Thermal expansion coefficient of fluid 0, fluid 1, ...', test)
          IF (test) THEN
             IF(.NOT. inputs%variation_temp_param_fluid) THEN
                CALL error_petsc('BUG in read_my_data: if define Thermal Expansion in fluid,' // &
                     ' one also needs to define Heat capacity in fluid')
             ELSE
                READ(21,*) inputs%heat_grav_fluid
             END IF
          ELSE
             inputs%heat_grav_fluid=0.d0
          END IF
       ELSE
          ALLOCATE(inputs%heat_capacity_fluid(0))
          ALLOCATE(inputs%heat_diffu_fluid(0))
          ALLOCATE(inputs%heat_grav_fluid(0))
       END IF
       !==========Surface tension=========================!
       CALL find_string(21, '===Is there a surface tension?', test)
       IF (test) THEN
          READ (21,*) inputs%if_surface_tension
       ELSE
          inputs%if_surface_tension = .FALSE.
       END IF
       IF (inputs%if_surface_tension) THEN
          !==========Coefficient for surface tension======!
          CALL read_until(21, '===Coefficients of surface tension for level set 0, level set 1, ...')
          ALLOCATE(inputs%coeff_surface(inputs%nb_fluid-1))
          READ(21,*) inputs%coeff_surface
       END IF
       !==========Mass correction=========================!
       CALL find_string(21, '===Do we apply mass correction? (true/false)', test)
       IF (test) THEN
          READ (21,*) inputs%if_mass_correction
       ELSE
          inputs%if_mass_correction=.TRUE.
       END IF
       !==========Dirichlet BCs on level_set==============!
       CALL read_until(21, '===How many boundary pieces for Dirichlet BCs on level set?')
       READ(21,*) inputs%level_set_nb_dirichlet_sides
       IF (inputs%level_set_nb_dirichlet_sides>0) THEN
          ALLOCATE(inputs%level_set_list_dirichlet_sides(inputs%level_set_nb_dirichlet_sides))
          CALL read_until(21, '===List of boundary pieces for Dirichlet BCs on level set')
          READ(21,*) inputs%level_set_list_dirichlet_sides
       ELSE
          ALLOCATE(inputs%level_set_list_dirichlet_sides(0))
       END IF
       !==========Solver parameters=======================!
       CALL read_until(21, '===Maximum number of iterations for level set solver')
       READ(21,*) inputs%my_par_level_set%it_max
       CALL read_until(21, '===Relative tolerance for level set solver')
       READ(21,*) inputs%my_par_level_set%rel_tol
       CALL read_until(21, '===Absolute tolerance for level set solver')
       READ(21,*) inputs%my_par_level_set%abs_tol
       CALL find_string(21, '===Level set solver verbose? (true/false)',test)
       IF (test) THEN
          READ(21,*) inputs%my_par_level_set%verbose
       END IF
       CALL read_until(21, '===Solver type for level set (FGMRES, CG, ...)')
       READ(21,*) inputs%my_par_level_set%solver
       CALL read_until(21, '===Preconditionner type for level set solver (HYPRE, JACOBI, MUMPS...)')
       READ(21,*) inputs%my_par_level_set%precond

       !==========Reconstruction parameters===============!
       CALL find_string(21, '===How are the variables reconstructed from the level set function? (lin, reg)', test)
       IF (test) THEN
          READ (21,*) inputs%level_set_reconstruction_type
       ELSE
          inputs%level_set_reconstruction_type = 'reg'
       END IF
       IF (inputs%level_set_reconstruction_type=='reg') THEN
          CALL find_string(21, '===Value of the regularization coefficient in (0,0.5]', test)
          IF (test) THEN
             READ (21,*) inputs%level_set_tanh_coeff_reconstruction
             IF (inputs%level_set_tanh_coeff_reconstruction.LE.1d-2 .OR. &
                  .5d0 .LE. inputs%level_set_tanh_coeff_reconstruction) THEN
                inputs%level_set_tanh_coeff_reconstruction = 0.45d0
             END IF
          ELSE
             inputs%level_set_tanh_coeff_reconstruction = 0.45d0
          END IF
       ELSE IF (inputs%level_set_reconstruction_type/='lin') THEN
          CALL error_petsc('BUG in read_my_data: variable reconstruction type not correct')
       END IF
       !==========Level set overshoot=====================!
       IF(inputs%level_set_reconstruction_type/='reg') THEN
          CALL find_string(21, '===Do we kill level set overshoot? (true/false)', test)
          IF (test) THEN
             READ (21,*) inputs%if_kill_overshoot
          ELSE
             inputs%if_kill_overshoot=.FALSE.
          END IF
       ELSE
          inputs%if_kill_overshoot=.FALSE.
       END IF

       !==========Level set type finite element===========!
       CALL find_string(21, '===Do we use P2 finite element for level_set? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%if_level_set_P2
       ELSE
          inputs%if_level_set_P2=.TRUE.
       END IF

       !==========Level set compression method============!
       CALL find_string(21, '===Do we use JLG compression method for level_set? (true/false)', test)
       IF (test) THEN
          READ(21,*) inputs%if_compression_mthd_JLG
       ELSE
          inputs%if_compression_mthd_JLG=.TRUE.
       END IF

       !==========Check coherence=========================!
       IF (inputs%type_pb=='mxw') THEN
          CALL error_petsc('BUG in read_my_data: Level_set without Navier-Stokes')
       END IF
    END IF

    !===Data for arpack (eigenvalue problems)==========!
    !==========Frequency parameters====================!
    IF (inputs%type_pb=='mxw') THEN
       CALL find_string(21, '===Do we use Arpack?', test)
       IF (test) THEN
          READ (21,*) inputs%if_arpack
       ELSE
          inputs%if_arpack = .FALSE.
       END IF
       IF (inputs%if_arpack)THEN
          CALL read_until(21, '===Number of eigenvalues to compute')
          READ (21,*) inputs%nb_eigenvalues
          CALL read_until(21, '===Maximum number of Arpack iterations')
          READ (21,*) inputs%arpack_iter_max
          CALL read_until(21, '===Tolerance for Arpack')
          READ (21,*) inputs%arpack_tolerance
          CALL read_until(21, &
               '===Which eigenvalues (''LM'', ''SM'', ''SR'', ''LR'' ''LI'', ''SI'')')
          READ (21,*) inputs%arpack_which
          CALL find_string(21, '===Create 2D vtu files for Arpack? (true/false)', test)
          IF (test) THEN
             READ (21,*) inputs%if_arpack_vtu_2d
          ELSE
             inputs%if_arpack_vtu_2d = .FALSE.
          END IF
       END IF
    ELSE
       !==========Arpack currently works only for Maxwell=!
       inputs%if_arpack = .FALSE.
    END IF

    !===Format for paraview============================!
    CALL find_string(21, '===Vtu files in xml format? (true=xml/false=ascii)', test)
    IF (test) THEN
       READ (21,*) inputs%if_xml
    ELSE
       inputs%if_xml=.TRUE.
    END IF

    !===Data for post processing=======================!
    CALL find_string(21, '===Number of planes in real space for Visualization', test)
    IF (test) THEN
       READ (21,*) inputs%number_of_planes_in_real_space
    ELSE
       inputs%number_of_planes_in_real_space = 10
    END IF

    !==========Frequency parameters====================!
    CALL find_string(21, '===Frequency to write restart file', test)
    IF (test) THEN
       READ (21,*) inputs%freq_restart
    ELSE
       inputs%freq_restart = 100000000
    END IF
    CALL find_string(21, '===Frequency to write energies', test)
    IF (test) THEN
       READ (21,*) inputs%freq_en
    ELSE
       inputs%freq_en = 100000000
    END IF
    CALL find_string(21, '===Frequency to create plots', test)
    IF (test) THEN
       READ (21,*) inputs%freq_plot
    ELSE
       inputs%freq_plot = 100000000
    END IF
    CALL find_string(21, '===Just postprocessing without computing? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%if_just_processing
    ELSE
       inputs%if_just_processing = .FALSE.
    END IF
    CALL find_string(21, '===Should I do post proc init? (true/false)', test)
    IF (test) THEN
       READ (21, *) inputs%if_post_proc_init
    ELSE
       inputs%if_post_proc_init=.FALSE.
    END IF
    CALL find_string(21, '===Create 2D plots (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%if_plot_2D
    ELSE
       inputs%if_plot_2D = .FALSE.
    END IF
    CALL find_string(21, '===Compute L2 and H1 relative errors (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%if_compute_error
    ELSE
       inputs%if_compute_error = .FALSE.
    END IF

    !==========Anemometer parameters===================!
    !===Anemometers for concentration
    CALL find_string(21, '===Anemometers (conc) ? (true/false)', test)
    IF (test) THEN
       READ (21, *) inputs%if_anemo_conc
       IF (inputs%if_anemo_conc) THEN
          CALL read_until(21, '===Number of anemo_conc (r,z)')
          READ (21, *) inputs%nb_anemo_r_conc, inputs%nb_anemo_z_conc
          ALLOCATE(inputs%r_anemo_conc(inputs%nb_anemo_r_conc))
          ALLOCATE(inputs%z_anemo_conc(inputs%nb_anemo_z_conc))
          CALL read_until(21, '===List of r anemo_conc')
          READ (21, *) inputs%r_anemo_conc
          CALL read_until(21, '===List of z anemo_conc')
          READ (21, *) inputs%z_anemo_conc
       ELSE
          inputs%if_anemo_conc = .FALSE.
          inputs%nb_anemo_r_conc=0
          inputs%nb_anemo_z_conc=0
       END IF
    ELSE
       inputs%if_anemo_conc = .FALSE.
       inputs%nb_anemo_r_conc=0
       inputs%nb_anemo_z_conc=0
    END IF
    !===Anemometers for velocity
    CALL find_string(21, '===Anemometers (v) ? (true/false)', test)
    IF (test) THEN
       READ (21, *) inputs%if_anemo_v
       IF (inputs%if_anemo_v) THEN
          CALL read_until(21, '===Number of anemo_v (r,z)')
          READ (21, *) inputs%nb_anemo_r_v, inputs%nb_anemo_z_v
          ALLOCATE(inputs%r_anemo_v(inputs%nb_anemo_r_v))
          ALLOCATE(inputs%z_anemo_v(inputs%nb_anemo_z_v))
          CALL read_until(21, '===List of r anemo_v')
          READ (21, *) inputs%r_anemo_v
          CALL read_until(21, '===List of z anemo_v')
          READ (21, *) inputs%z_anemo_v
       ELSE
          inputs%if_anemo_v = .FALSE.
          inputs%nb_anemo_r_v=0
          inputs%nb_anemo_z_v=0
       END IF
    ELSE
       inputs%if_anemo_v = .FALSE.
       inputs%nb_anemo_r_v=0
       inputs%nb_anemo_z_v=0
    END IF
    !===Anemometers for temperature
    CALL find_string(21, '===Anemometers (T) ? (true/false)', test)
    IF (test) THEN
       READ (21, *) inputs%if_anemo_T
       IF (inputs%if_anemo_T) THEN
          CALL read_until(21, '===Number of anemo_T (r,z)')
          READ (21, *) inputs%nb_anemo_r_T, inputs%nb_anemo_z_T
          ALLOCATE(inputs%r_anemo_T(inputs%nb_anemo_r_T))
          ALLOCATE(inputs%z_anemo_T(inputs%nb_anemo_z_T))
          CALL read_until(21, '===List of r anemo_T')
          READ (21, *) inputs%r_anemo_T
          CALL read_until(21, '===List of z anemo_T')
          READ (21, *) inputs%z_anemo_T
       ELSE
          inputs%if_anemo_T = .FALSE.
          inputs%nb_anemo_r_T=0
          inputs%nb_anemo_z_T=0
       END IF
    ELSE
       inputs%if_anemo_T = .FALSE.
       inputs%nb_anemo_r_T=0
       inputs%nb_anemo_z_T=0
    END IF
    !===Anemometers for magnetic field
    CALL find_string(21, '===Anemometers (H) ? (true/false)', test)
    IF (test) THEN
       READ (21, *) inputs%if_anemo_h
       IF (inputs%if_anemo_h) THEN
          CALL read_until(21, '===Number of anemo_h (r,z)')
          READ (21, *) inputs%nb_anemo_r_h, inputs%nb_anemo_z_h
          ALLOCATE(inputs%r_anemo_h(inputs%nb_anemo_r_h))
          ALLOCATE(inputs%z_anemo_h(inputs%nb_anemo_z_h))
          CALL read_until(21, '===List of r anemo_h')
          READ (21, *) inputs%r_anemo_h
          CALL read_until(21, '===List of z anemo_h')
          READ (21, *) inputs%z_anemo_h
       ELSE
          inputs%if_anemo_h = .FALSE.
          inputs%nb_anemo_r_h=0
          inputs%nb_anemo_z_h=0
       END IF
    ELSE
       inputs%if_anemo_h = .FALSE.
       inputs%nb_anemo_r_h=0
       inputs%nb_anemo_z_h=0
    END IF


    !==========Modes to be zeroed out==================!
    IF (inputs%type_pb=='mhd' .OR. inputs%type_pb=='mhs') THEN
       CALL find_string(21, '===Should some modes be zeroed out?', test)
       IF (test) THEN
          READ (21,*) inputs%if_zero_out_modes
       ELSE
          inputs%if_zero_out_modes = .FALSE.
       END IF
       IF (inputs%if_zero_out_modes) THEN
          CALL read_until(21, '===How many Navier-Stokes modes to zero out?')
          READ(21,*) inputs%nb_select_mode_ns
          ALLOCATE(inputs%list_select_mode_ns(inputs%nb_select_mode_ns))
          CALL read_until(21, '===List of Navier-Stokes modes to zero out?')
          READ(21,*) inputs%list_select_mode_ns
          CALL read_until(21, '===How Maxwell modes to zero out?')
          READ(21,*) inputs%nb_select_mode_mxw
          ALLOCATE(inputs%list_select_mode_mxw(inputs%nb_select_mode_mxw))
          CALL read_until(21, '===List of Maxwell modes to zero out?')
          READ(21,*) inputs%list_select_mode_mxw
       END IF
    END IF

    !==========Verbose=================================!
    CALL find_string(21, '===Verbose timing? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%verbose_timing
    ELSE
       inputs%verbose_timing = .FALSE.
    END IF
    CALL find_string(21, '===Verbose divergence? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%verbose_divergence
    ELSE
       inputs%verbose_divergence = .FALSE.
    END IF
    CALL find_string(21, '===Verbose CFL? (true/false)', test)
    IF (test) THEN
       READ (21,*) inputs%verbose_CFL
    ELSE
       inputs%verbose_CFL = .FALSE.
    END IF

    !===Norms for reference tests===================================================
    IF (inputs%test_de_convergence) THEN
       CALL read_until(21, '===Reference results')
       DO k = 1, 4
          READ(21,*) inputs%norm_ref(k)
       END DO
    ELSE
       inputs%norm_ref = 1.d0
    END IF

    !===mxw becomes mxx if restart_u and mxw========================================
    IF (inputs%type_pb == 'mxw' .AND. inputs%irestart_u) THEN
       inputs%type_pb = 'mxx'
    END IF

    CLOSE(21)

    !===Check coherence of data=====================================================
    CALL check_coherence_of_data
    RETURN

  END SUBROUTINE read_my_data

  SUBROUTINE check_coherence_of_data
    USE my_util
    IMPLICIT NONE
    INTEGER :: k, mode_min, mode_max, Delta_mode
    LOGICAL :: test

    !===Dirichlet BCs for Navier-Stokes=============================================
    IF (inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
            .OR. inputs%type_pb=='mhs'  &
            .OR. inputs%irestart_u) THEN
          !==========Velocity================================!
          DO k = 1, 3
             IF (inputs%vv_nb_dirichlet_sides(k)<1) CYCLE
             test = check_coherence_with_periodic_bcs(inputs%vv_list_dirichlet_sides(k)%DIL)
             IF (test) THEN
                CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                     ' and periodic BCs on velocity')
             END IF
          END DO
          !==========Pressure================================!
          IF (inputs%pp_nb_dirichlet_sides>0) THEN
             test = check_coherence_with_periodic_bcs(inputs%pp_list_dirichlet_sides)
             IF (test) THEN
                CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                     ' and periodic BCs on pressure')
             END IF
          END IF
       END IF
    END IF

    !===No penalty with momentum formulation========================================
    IF (.NOT.inputs%if_navier_stokes_with_u .AND. inputs%if_ns_penalty) THEN
       CALL error_petsc('BUG in read_my_data: No penalty with momentum formulation')
    END IF

    !===Dirichlet BCs for temperature for Navier-Stokes=============================
    IF (inputs%if_temperature.AND. inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%temperature_nb_dirichlet_sides>0) THEN
          test = check_coherence_with_periodic_bcs(inputs%temperature_list_dirichlet_sides)
          IF (test) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and periodic BCs on temperature')
          END IF
       END IF
    END IF
    !===Dirichlet BCs for concentration for Navier-Stokes=============================
    IF (inputs%if_concentration.AND. inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%concentration_nb_dirichlet_sides>0) THEN
          test = check_coherence_with_periodic_bcs(inputs%concentration_list_dirichlet_sides)
          IF (test) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and periodic BCs on concentration')
          END IF
       END IF
    END IF
    !===Dirichlet BCs for Level_set for Navier-Stokes=============================
    IF (inputs%if_level_set.AND. inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%level_set_nb_dirichlet_sides>0) THEN
          test = check_coherence_with_periodic_bcs(inputs%level_set_list_dirichlet_sides)
          IF (test) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and periodic BCs on level_set')
          END IF
       END IF
    END IF

    !===Robin and Dirichlet BCs===================================================== ! MODIFICATION/ The user should not specify sides that are Robin and Dirichlet
    IF (inputs%temperature_nb_robin_sides>0) THEN
       DO k = 1, inputs%temperature_nb_robin_sides
          IF (MINVAL(ABS(inputs%temperature_list_dirichlet_sides - inputs%temperature_list_robin_sides(k))) == 0) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and Robin BCs for temperature')
          END IF
       END DO
    END IF

    !===Dirichlet BCs for concentration for Navier-Stokes=============================
    IF (inputs%if_concentration.AND. inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%concentration_nb_dirichlet_sides>0) THEN
          test = check_coherence_with_periodic_bcs(inputs%concentration_list_dirichlet_sides)
          IF (test) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and periodic BCs on concentration')
          END IF
       END IF
    END IF

    !===Robin and Dirichlet BCs===================================================== ! MODIFICATION/ The user should not specify sides that are Robin and Dirichlet
    IF (inputs%concentration_nb_robin_sides>0) THEN
       DO k = 1, inputs%concentration_nb_robin_sides
          IF (MINVAL(ABS(inputs%concentration_list_dirichlet_sides - inputs%concentration_list_robin_sides(k))) == 0) THEN
             CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                  ' and Robin BCs for concentration')
          END IF
       END DO
    END IF

    !===Dirichlet BCs magnetic field for Maxwell====================================
    IF (inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%type_pb=='mxw' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd' &
            .OR. inputs%type_pb=='mhs') THEN
          !==========Magnetic field==========================!
          IF (inputs%nb_dirichlet_sides_H>0) THEN
             test = check_coherence_with_periodic_bcs(inputs%list_dirichlet_sides_H)
             IF (test) THEN
                CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                     ' and periodic BCs on magnetic field')
             END IF
          END IF

          !==========Scalar potential========================!
          IF (inputs%nb_dom_phi>0 .AND. inputs%phi_nb_dirichlet_sides>0) THEN
             test = check_coherence_with_periodic_bcs(inputs%phi_list_dirichlet_sides)
             IF (test) THEN
                CALL error_petsc(' BUG in read_my_data: Incompatible Dirichlet'// &
                     ' and periodic BCs on scalar potential')
             END IF
          END IF
       END IF
    END IF

    !===Check temperature with Maxwell==============================================
    IF (inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%if_temperature .AND. inputs%type_pb=='mxw') THEN
          CALL error_petsc('Bug in read_my_data: incompatible temperature with maxwell')
       END IF
    END IF

    !===Check concentration with Maxwell==============================================
    IF (inputs%my_periodic%nb_periodic_pairs>0) THEN
       IF (inputs%if_concentration .AND. inputs%type_pb=='mxw') THEN
          CALL error_petsc('Bug in read_my_data: incompatible concentration with maxwell')
       END IF
    END IF

    !===Check temperature with Ferrohydrodynamics===================================
    IF ((.NOT. inputs%if_temperature) .AND. inputs%type_pb=='fhd') THEN
       CALL error_petsc('Bug in read_my_data: ferrohydrodynamics but no temperature')
    END IF

    !===Check Arpack================================================================
    IF (inputs%if_arpack) THEN
       IF (inputs%ndim(2) /= inputs%m_max) THEN
          CALL error_petsc('Bug in read_my_data: #Fourier modes'// &
               ' not equal to #processors in Fourier direction')
       END IF
    END IF

    !===Check Fourier modes=========================================================
    IF (inputs%select_mode .AND. .NOT.inputs%if_arpack) THEN
       IF (SIZE(inputs%list_mode_lect)/=1) THEN
          mode_max = MAXVAL(inputs%list_mode_lect)
          mode_min = MINVAL(inputs%list_mode_lect)
          IF (MOD(mode_max-mode_min,SIZE(inputs%list_mode_lect)-1)/=0) THEN
             CALL error_petsc('Bug in read_my_data: Fourier modes not equally spaced ')
          END IF
          Delta_mode = (mode_max-mode_min)/(SIZE(inputs%list_mode_lect)-1)
          DO k = 0, SIZE(inputs%list_mode_lect)-1
             IF (MINVAL(ABS(inputs%list_mode_lect-(Delta_mode*k+mode_min)))/=0) THEN
                CALL error_petsc('Bug in read_my_data: Fourier modes not equally spaced ')
             END IF
          END DO
          DO k = 1, SIZE(inputs%list_mode_lect)-1
             IF (inputs%list_mode_lect(k+1).LE.inputs%list_mode_lect(k)) THEN
                CALL error_petsc('Bug in read_my_data: Fourier modes not in increasing order ')
             END IF
          END DO
       END IF
    END IF

    !===Irestart for postprocessing=================================================
    IF (inputs%if_just_processing) THEN
       inputs%irestart_u = .FALSE.
       inputs%irestart_h = .FALSE.
       IF (inputs%type_pb/='mxw') inputs%irestart_u = .TRUE.
       IF (inputs%type_pb/='nst') inputs%irestart_h = .TRUE.
    END IF

    !===Allocate dummy list_dom_* for new partitioning==============================
    IF (.NOT. ASSOCIATED(inputs%list_dom_ns))  ALLOCATE(inputs%list_dom_ns(0))
    IF (.NOT. ASSOCIATED(inputs%list_dom_H))   ALLOCATE(inputs%list_dom_H(0))
    IF (.NOT. ASSOCIATED(inputs%list_dom_phi)) ALLOCATE(inputs%list_dom_phi(0))
    IF (.NOT. ASSOCIATED(inputs%list_dom_temp))   ALLOCATE(inputs%list_dom_temp(0))
    IF (.NOT. ASSOCIATED(inputs%list_dom_conc))   ALLOCATE(inputs%list_dom_conc(0))

  END SUBROUTINE check_coherence_of_data

  FUNCTION check_coherence_with_periodic_bcs(list) RESULT(test)
    IMPLICIT NONE
    INTEGER, DIMENSION(:) :: list
    LOGICAL               :: test
    INTEGER               :: k, n

    test = .FALSE.
    DO k = 1, inputs%my_periodic%nb_periodic_pairs
       DO n = 1, SIZE(list)
          IF (MINVAL(ABS(list(n)-inputs%my_periodic%list_periodic(:,k)))==0) THEN
             test = .TRUE.
             RETURN
          END IF
       END DO
    END DO
  END FUNCTION check_coherence_with_periodic_bcs

END MODULE input_data
