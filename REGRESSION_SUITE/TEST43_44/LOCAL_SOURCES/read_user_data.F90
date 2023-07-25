MODULE user_data_module
  TYPE personalized_data
     !===I declare my own data here==================================================
     LOGICAL                                 :: if_my_stuff
     LOGICAL                                 :: if_post_proc_init
     LOGICAL                                 :: if_anemo_h, if_anemo_v, if_anemo_T, if_anemo_conc
     LOGICAL                                 :: if_confine
     REAL(KIND=8)                            :: J0
     LOGICAL                                 :: if_Hz
     REAL(KIND=8)                            :: Hz
     LOGICAL                                 :: if_Rcc, if_Rfoam
     REAL(KIND=8)                            :: Rcc, Rmax, Rfoam
     INTEGER                                 :: nb_surf_temp
     INTEGER, DIMENSION(:), POINTER          :: list_surf_sides
     LOGICAL                                 :: if_grav_tilt
     REAL(KIND=8)                            :: time_grav_tilt
     LOGICAL                                 :: if_nonlin
     REAL(KIND=8)                            :: amp_nl
     LOGICAL                                 :: if_LMW_LS
     REAL(KIND=8)                            :: HsR, Fr
     REAL(KIND=8)                            :: height_sigma_bar
     INTEGER                                 :: nb_anemo_r_h, nb_anemo_z_h
     INTEGER                                 :: nb_anemo_r_v, nb_anemo_z_v
     INTEGER                                 :: nb_anemo_r_T, nb_anemo_z_T
     INTEGER                                 :: nb_anemo_r_conc, nb_anemo_z_conc
     REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_anemo_h, z_anemo_h
     REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_anemo_v, z_anemo_v
     REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_anemo_T, z_anemo_T
     REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: r_anemo_conc, z_anemo_conc
     LOGICAL                                 :: if_density_law, if_sigma_law
     REAL(KIND=8)                            :: hmin_mult
     REAL(KIND=8)                            :: rho_lin_coeff, rho_ref, rho_min, rho_max
     REAL(KIND=8), DIMENSION(2)              :: rho_range
     REAL(KIND=8)                            :: one_over_sigma_lin_coeff, sigma_ref, sigma_min, sigma_max
     REAL(KIND=8), DIMENSION(2)              :: sigma_range
     REAL(KIND=8)                            :: h_half
     INTEGER                                 :: nb_surf_conc
     INTEGER, DIMENSION(:), POINTER          :: list_surf_sides_conc
     LOGICAL                                 :: if_evf, if_buoyancy

     !===Validation tests===========================================================
     LOGICAL                                 :: validation_test
     INTEGER                                 :: nb_validation_test
     !.......Continue here ................................
  END TYPE personalized_data
END MODULE user_data_module

MODULE user_data
  USE user_data_module
  IMPLICIT NONE
  PUBLIC :: read_user_data
  TYPE(personalized_data), PUBLIC  :: user
  PRIVATE

CONTAINS

  SUBROUTINE read_user_data(data_file)
    USE my_util
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(*),       INTENT(IN) :: data_file
    INTEGER                        :: unit_file=22
    LOGICAL                        :: test

    OPEN(UNIT=unit_file, FILE = data_file, FORM = 'formatted', STATUS = 'unknown')

    !===Template====================================================================
    CALL find_string(unit_file, '===Should I read my stuff? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_my_stuff
    ELSE
       user%if_my_stuff = .FALSE.
    END IF
    !.......Continue here ................................

    !===Post proc init==============================================================
    CALL find_string(unit_file, '===Should I do post proc init? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_post_proc_init
    ELSE
       user%if_post_proc_init=.FALSE.
    END IF

    !===Nonlinear restart================================================
    CALL find_string(unit_file, '===Is it a nonlinear restart on modes/=0 ? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_nonlin
       IF (user%if_nonlin) THEN
          CALL find_string(unit_file, '===Amplitude of amp_nl', test)
          READ(unit_file,*) user%amp_nl
       ELSE
          user%if_nonlin = .FALSE.
          user%amp_nl = 1.d0
       END IF
    ELSE
       user%if_nonlin = .FALSE.
       user%amp_nl = 1.d0
    END IF

    !===Confinement for RFP==============================================
    CALL find_string(unit_file, '===Is there a vertical current ? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_confine
       IF (user%if_confine) THEN
          CALL find_string(unit_file, '===Amplitude of J0', test)
          READ(unit_file,*) user%J0
       ELSE
          user%if_confine = .FALSE.
          user%J0 = 0.d0
       END IF
    ELSE
       user%if_confine = .FALSE.
       user%J0 = 0.d0
    END IF

    !===Vertical magnetic field for Metal Pad Roll=========================
    CALL find_string(unit_file, '===Is there a vertical magnetic field ? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_Hz
       IF (user%if_Hz) THEN
          CALL find_string(unit_file, '===Amplitude of Hz', test)
          READ(unit_file,*) user%Hz
       ELSE
          user%if_Hz = .FALSE.
          user%Hz = 0.d0
       END IF
    ELSE
       user%if_Hz = .FALSE.
       user%Hz = 0.d0
    END IF

    !===Solutal Convection================================
    CALL find_string(unit_file, '===Is there a density law for solutal convection? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_density_law
       IF (user%if_density_law) THEN
          CALL find_string(unit_file, '===Linear coefficient for density law', test)
          READ(unit_file,*) user%rho_lin_coeff
          CALL find_string(unit_file, '===Reference density for density law', test)
          READ(unit_file,*) user%rho_ref
          CALL find_string(unit_file, '===Minimum density', test)
          READ(unit_file,*) user%rho_min
          CALL find_string(unit_file, '===Maximum density', test)
          READ(unit_file,*) user%rho_max
          user%rho_range(1)=user%rho_min
          user%rho_range(2)=user%rho_max
          CALL find_string(unit_file, '===Multiplier coefficient hmin', test)
          READ(unit_file,*) user%hmin_mult
       ELSE
          user%if_density_law = .FALSE.
          user%rho_lin_coeff = 0.d0
          user%rho_ref = 1.d0
          user%rho_min=0.d0
          user%rho_max=0.d0
          user%rho_range=0.d0
          user%hmin_mult=1.d0
       END IF
    ELSE
       user%if_density_law = .FALSE.
       user%rho_lin_coeff = 0.d0
       user%rho_ref = 1.d0
       user%rho_min=0.d0
       user%rho_max=0.d0
       user%rho_range=0.d0
       user%hmin_mult=1.d0
    END IF

    CALL find_string(unit_file, '===Is there a conductivity law for solutal convection? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_sigma_law
       IF (user%if_density_law) THEN
          CALL find_string(unit_file, '===Linear coefficient for one over conductivity law', test)
          READ(unit_file,*) user%one_over_sigma_lin_coeff
          CALL find_string(unit_file, '===Reference conductivity for conductivity law', test)
          READ(unit_file,*) user%sigma_ref
          CALL find_string(unit_file, '===Minimum sigma', test)
          READ(unit_file,*) user%sigma_min
          CALL find_string(unit_file, '===Maximum sigma', test)
          READ(unit_file,*) user%sigma_max
          user%sigma_range(1)=user%sigma_min
          user%sigma_range(2)=user%sigma_max
       ELSE
          user%if_sigma_law = .FALSE.
          user%one_over_sigma_lin_coeff = 0.d0
          user%sigma_ref = 1.d0
          user%sigma_min=0.d0
          user%sigma_max=0.d0
          user%sigma_range=0.d0
       END IF
    ELSE
       user%if_sigma_law = .FALSE.
       user%one_over_sigma_lin_coeff = 0.d0
       user%sigma_ref = 1.d0
       user%sigma_min=0.d0
       user%sigma_max=0.d0
       user%sigma_range=0.d0
    END IF

    !===Current collector for EVF=========================
    CALL find_string(unit_file, '===Is there a current collector ?  (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_Rcc
       IF (user%if_Rcc) THEN
          CALL find_string(unit_file, '===Length of Rcc', test)
          READ(unit_file,*) user%Rcc
          CALL find_string(unit_file, '===Maximum radius Rmax', test)
          READ(unit_file,*) user%Rmax
       ELSE
          user%if_Rcc = .FALSE.
          user%Rcc = 0.d0
          user%Rmax = 0.d0
       END IF
    ELSE
       user%if_Rcc = .FALSE.
       user%Rcc = 0.d0
       user%Rmax = 0.d0
    END IF
    CALL find_string(unit_file, '===Is there a foam ?  (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_Rfoam
       IF (user%if_Rfoam) THEN
          CALL find_string(unit_file, '===Length of Rfoam', test)
          READ(unit_file,*) user%Rfoam
       ELSE
          user%if_Rfoam = .FALSE.
          user%Rfoam = 0.d0
       END IF
    ELSE
       user%if_Rfoam = .FALSE.
       user%Rfoam = 0.d0
    END IF

    ! === EVF
    CALL find_string(unit_file, '===Is there EVF ?  (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_evf
    ELSE
       user%if_evf = .FALSE.
    END IF

    ! === Buoyancy
    CALL find_string(unit_file, '===Is there buoyancy ?  (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_buoyancy
    ELSE
       user%if_buoyancy = .FALSE.
    END IF

    !===Surface integral of temperature===========================
    CALL find_string(unit_file, '===How many boundary pieces for surface integral', test)
    IF (test) THEN
       READ(unit_file,*) user%nb_surf_temp
    ELSE
       user%nb_surf_temp = 0
    END IF
    IF (user%nb_surf_temp > 0) THEN
       ALLOCATE(user%list_surf_sides(user%nb_surf_temp))
       CALL find_string(unit_file, '===List of boundary pieces for surface integral',test)
       READ(unit_file,*) user%list_surf_sides
    ELSE
       ALLOCATE(user%list_surf_sides(0))
    END IF

    !===Surface integral of concentration===========================
    CALL find_string(unit_file, '===How many boundary pieces for surface integral of concentration', test)
    IF (test) THEN
       READ(unit_file,*) user%nb_surf_conc
    ELSE
       user%nb_surf_conc = 0
    END IF
    IF (user%nb_surf_conc > 0) THEN
       ALLOCATE(user%list_surf_sides_conc(user%nb_surf_conc))
       CALL find_string(unit_file, '===List of boundary pieces for surface integral of concentration',test)
       READ(unit_file,*) user%list_surf_sides_conc
    ELSE
       ALLOCATE(user%list_surf_sides_conc(0))
    END IF

    !===Tilted gravity for Metal Pad Roll=========================
    CALL find_string(unit_file, '===Is there a tilted gravity ? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_grav_tilt
       IF (user%if_grav_tilt) THEN
          CALL find_string(unit_file, '===Time of grav_tilt', test)
          READ(unit_file,*) user%time_grav_tilt
       ELSE
          user%if_grav_tilt = .FALSE.
          user%time_grav_tilt = 0.d0
       END IF
    ELSE
       user%if_grav_tilt = .FALSE.
       user%time_grav_tilt = 0.d0
    END IF

    !===Level Set for LMW cases=============================================
    CALL find_string(unit_file, '===Is it a LMW LSet ? (true/false)', test)
    IF (test) THEN
       READ(unit_file,*) user%if_LMW_LS
       IF (user%if_LMW_LS) THEN
          CALL find_string(unit_file, '===HfluidsR_init, Fr', test)
          READ(unit_file,*) user%HsR, user%Fr !, user%ep
       ELSE
          user%if_LMW_LS = .FALSE.
          user%HsR = 0.d0
          user%Fr = 0.d0
       END IF
    ELSE
       user%if_LMW_LS = .FALSE.
       user%HsR = 0.d0
       user%Fr = 0.d0
    END IF

    !===Height interfaces for 2 layers ======================================
    CALL find_string(unit_file, '===Height interfaces for layers: h_half', test)
    IF (test) THEN
       READ(unit_file,*) user%h_half
    ELSE
       user%h_half=0.d0
    END IF

    !===Height interface for sigma bar======================================
    CALL find_string(unit_file, '===Height interface for sigma bar', test)
    IF (test) THEN
       READ(unit_file,*) user%height_sigma_bar
    ELSE
       user%height_sigma_bar=0.d0
    END IF

    !===Anemometers for velocity==============================================
    CALL find_string(unit_file, '===Anemometers (v) ? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_anemo_v
       IF (user%if_anemo_v) THEN
          CALL read_until(unit_file, '===Number of anemo_v (r,z)')
          READ (unit_file, *) user%nb_anemo_r_v, user%nb_anemo_z_v
          ALLOCATE(user%r_anemo_v(user%nb_anemo_r_v))
          ALLOCATE(user%z_anemo_v(user%nb_anemo_z_v))
          CALL read_until(unit_file, '===List of r anemo_v')
          READ (unit_file, *) user%r_anemo_v
          CALL read_until(unit_file, '===List of z anemo_v')
          READ (unit_file, *) user%z_anemo_v
       ELSE
          user%if_anemo_v = .false.
          user%nb_anemo_r_v=0
          user%nb_anemo_z_v=0
       END IF
    ELSE
       user%if_anemo_v = .false.
       user%nb_anemo_r_v=0
       user%nb_anemo_z_v=0
    END IF

    !===Anemometers for H======================================================
    CALL find_string(unit_file, '===Anemometers (H) ? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_anemo_h
       IF (user%if_anemo_h) THEN
          CALL read_until(unit_file, '===Number of anemo_h (r,z)')
          READ (unit_file, *) user%nb_anemo_r_h, user%nb_anemo_z_h
          ALLOCATE(user%r_anemo_h(user%nb_anemo_r_h))
          ALLOCATE(user%z_anemo_h(user%nb_anemo_z_h))
          CALL read_until(unit_file, '===List of r anemo_h')
          READ (unit_file, *) user%r_anemo_h
          CALL read_until(unit_file, '===List of z anemo_h')
          READ (unit_file, *) user%z_anemo_h
       ELSE
          user%if_anemo_h = .false.
          user%nb_anemo_r_h=0
          user%nb_anemo_z_h=0
       END IF
    ELSE
       user%if_anemo_h = .false.
       user%nb_anemo_r_h=0
       user%nb_anemo_z_h=0
    END IF

    !===Anemometers for T======================================================
    CALL find_string(unit_file, '===Anemometers (T) ? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_anemo_T
       IF (user%if_anemo_T) THEN
          CALL read_until(unit_file, '===Number of anemo_T (r,z)')
          READ (unit_file, *) user%nb_anemo_r_T, user%nb_anemo_z_T
          ALLOCATE(user%r_anemo_T(user%nb_anemo_r_T))
          ALLOCATE(user%z_anemo_T(user%nb_anemo_z_T))
          CALL read_until(unit_file, '===List of r anemo_T')
          READ (unit_file, *) user%r_anemo_T
          CALL read_until(unit_file, '===List of z anemo_T')
          READ (unit_file, *) user%z_anemo_T
       ELSE
          user%if_anemo_T = .false.
          user%nb_anemo_r_T=0
          user%nb_anemo_z_T=0
       END IF
    ELSE
       user%if_anemo_T = .false.
       user%nb_anemo_r_T=0
       user%nb_anemo_z_T=0
    END IF
    !===Anemometers for conc======================================================
    CALL find_string(unit_file, '===Anemometers (conc) ? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_anemo_conc
       IF (user%if_anemo_conc) THEN
          CALL read_until(unit_file, '===Number of anemo_conc (r,z)')
          READ (unit_file, *) user%nb_anemo_r_conc, user%nb_anemo_z_conc
          ALLOCATE(user%r_anemo_conc(user%nb_anemo_r_conc))
          ALLOCATE(user%z_anemo_conc(user%nb_anemo_z_conc))
          CALL read_until(unit_file, '===List of r anemo_conc')
          READ (unit_file, *) user%r_anemo_conc
          CALL read_until(unit_file, '===List of z anemo_conc')
          READ (unit_file, *) user%z_anemo_conc
       ELSE
          user%if_anemo_conc = .false.
          user%nb_anemo_r_conc=0
          user%nb_anemo_z_conc=0
       END IF
    ELSE
       user%if_anemo_conc = .false.
       user%nb_anemo_r_conc=0
       user%nb_anemo_z_conc=0
    END IF

    !===End template=================================================================

    CLOSE(unit_file)
  END SUBROUTINE read_user_data

END MODULE user_data
