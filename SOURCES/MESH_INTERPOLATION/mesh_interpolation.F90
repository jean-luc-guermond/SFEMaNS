MODULE mesh_interpolation

  PUBLIC :: mesh_interpol
  PRIVATE
CONTAINS
  SUBROUTINE mesh_interpol

    USE chaine_caractere
    USE def_type_mesh
    USE my_util
    USE create_comm
    USE interpolation_tools
    USE metis_sfemans
    USE prep_maill
    USE sub_plot
    USE restart
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE

    TYPE(mesh_type)                       :: p1_mesh_glob, p2_mesh_glob
    TYPE(mesh_type)                       :: p2_c0_mesh_glob_temp, p2_c0_mesh_glob_conc

    TYPE(mesh_type), TARGET               :: H_mesh, phi_mesh
    TYPE(mesh_type), TARGET               :: H_mesh_glob, phi_mesh_glob
    TYPE(mesh_type), POINTER              :: H_mesh_in, H_mesh_out, phi_mesh_in, phi_mesh_out

    TYPE(mesh_type), TARGET               :: vv_mesh, pp_mesh
    TYPE(mesh_type), TARGET               :: vv_mesh_glob, pp_mesh_glob
    TYPE(mesh_type), POINTER              :: vv_mesh_in, vv_mesh_out, pp_mesh_in, pp_mesh_out

    TYPE(mesh_type), TARGET               :: temp_mesh
    TYPE(mesh_type), TARGET               :: temp_mesh_glob
    TYPE(mesh_type), POINTER              :: temp_mesh_in, temp_mesh_out

    TYPE(mesh_type), TARGET               :: conc_mesh
    TYPE(mesh_type), TARGET               :: conc_mesh_glob
    TYPE(mesh_type), POINTER              :: conc_mesh_in, conc_mesh_out

    CHARACTER(len=200)                    :: old_filename, old_directory
    CHARACTER(len=200)                    :: new_filename, new_directory
    CHARACTER(len=200)                    :: file_name, directory
    CHARACTER(len=200)                    :: file_name_m, directory_m
    LOGICAL                               :: if_momentum, if_mass, if_induction, if_energy, if_conc, &
         old_is_form, new_is_form, iformatted, is_form_m, mono_in, mono_out, is_in, &
         if_concentration_in, if_concentration_out, if_concentration, &
         if_temperature_in, if_temperature_out, if_temperature, if_level_set, inter_mesh, &
         rw_conc, rw_ns, rw_temp, rw_mxw, check_plt, test, if_select

    INTEGER                               :: nb_S, nb_F, nb_procs, rank, type_fe_H, type_fe_phi, &
         nb_dom_conc, nb_dom_ns, nb_dom_temp, nb_dom_H, nb_dom_phi, nsize, code, m, &
         m_max_c, rank_S, nb_inter, nb_inter_mu, nb_inter_c_v, nb_inter_v_T, &
         k, kp, i, nb_mode, rang_conc_S, rang_ns_S, rang_temp_S, rang_S, l, lblank, nb_fluid
    INTEGER                               :: nb_fic, index_start
    REAL(KIND=8)                          :: time_h, time_u, time_T, time_conc, max_vel
    TYPE(periodic_data)                   :: my_periodic

    INTEGER, DIMENSION(:), ALLOCATABLE    :: list_inter, list_inter_temp, list_inter_conc, list_inter_v_T, &
         list_inter_c_v, list_inter_mu, list_inter_H_phi, list_dom_H, list_dom_conc, list_dom_ns, &
         list_dom_ns_in, list_dom_temp, list_dom_temp_in, list_dom_phi, list_dom_H_in, list_dom, part, &
         list_mode, controle_H, controle_phi, &
         controle_vv, controle_pp, controle_temp, controle_conc, vel_in_to_new, temp_in_to_new, H_in_to_new, &
         l_t_g_vv, l_t_g_pp, l_t_g_H, l_t_g_phi, l_t_g_temp, l_t_g_conc
    INTEGER, DIMENSION(:), ALLOCATABLE    :: list_dom_H_ref, H_in_to_new_ref, list_dom_temp_ref, temp_in_to_new_ref
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: Hn, Hn1, Bn, Bn1, phin, phin1
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: Hn_glob, Hn1_glob, Bn_glob, Bn1_glob, phin_glob, phin1_glob
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: un, un_m1, pn, pn_m1
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: incpn, incpn_m1, tempn, tempn_m1
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: concn, concn_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER     :: level_setn, level_setn_m1
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: un_glob, un_m1_glob, pn_glob, pn_m1_glob
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: incpn_glob, incpn_m1_glob, tempn_glob, tempn_m1_glob
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: concn_glob, concn_m1_glob
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER     :: level_setn_glob, level_setn_m1_glob

    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: Hn_in, Hn1_in, Bn_in, Bn1_in, phin_in, phin1_in
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: Hn_out, Hn1_out, Bn_out, Bn1_out, phin_out, phin1_out
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: un_in, un_m1_in, pn_in, pn_m1_in
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: incpn_in, incpn_m1_in, tempn_in, tempn_m1_in
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: concn_in, concn_m1_in
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER     :: level_setn_in, level_setn_m1_in
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: un_out, un_m1_out, pn_out, pn_m1_out
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: incpn_out, incpn_m1_out, tempn_out, tempn_m1_out
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER     :: concn_out, concn_m1_out
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER     :: level_setn_out, level_setn_m1_out

    CHARACTER(len=3)             :: type_pb, tit_m, tit_s, tit
    CHARACTER(len=4)             :: tit_part
    CHARACTER(len=200)           :: mesh_part_name
    INTEGER                      :: petsc_rank
    LOGICAL                      :: if_read_partition
    LOGICAL                      :: if_level_set_P2
    MPI_Comm                                        :: comm_cart
    MPI_Comm, DIMENSION(:), POINTER                 :: comm_one_d, comm_one_d_ns, &
         comm_one_d_temp, comm_one_d_conc, coord_cart
    PetscErrorCode                                  :: ierr

    CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, petsc_rank, code)

    OPEN(UNIT=22, FILE='data_interpol', STATUS='unknown', ACCESS='sequential')

    !===Compute is_in (is_in=.TRUE. if Petsc -> NonPetsc; is_in=.FALSE. if NonPetsc -> Petsc)
    CALL which_pb(is_in)

    IF (is_in) THEN
       CALL read_until(22, '===Number of processors in meridian section (Input)')
       READ(22,*) nb_S
       CALL read_until(22, '===Problem type (Input): (nst, mxw, mhd, fhd, mhs)')
       READ(22,*) type_pb
       CALL find_string(22, '===Is there an input concentration field?', test)
       IF (test) THEN
          READ (22, *) if_concentration_in
       ELSE
          if_concentration_in = .FALSE.
       END IF
       CALL find_string(22, '===Is there an input temperature field?', test)
       IF (test) THEN
          READ (22, *) if_temperature_in
       ELSE
          if_temperature_in = .FALSE.
       END IF
       if_temperature = if_temperature_in
       if_concentration = if_concentration_in
       inter_mesh = .FALSE.
       if_read_partition = .TRUE.
    ELSE
       CALL read_until(22, '===Number of processors in meridian section (Output)')
       READ(22,*) nb_S
       CALL read_until(22, '===Problem type (Output): (nst, mxw, mhd, fhd, mhs)')
       READ(22,*) type_pb
       CALL find_string(22, '===Is there an output concentration field?', test)
       IF (test) THEN
          READ (22, *) if_concentration_out
       ELSE
          if_concentration_out = .FALSE.
       END IF
       CALL find_string(22, '===Is there an output temperature field?', test)
       IF (test) THEN
          READ (22, *) if_temperature_out
       ELSE
          if_temperature_out = .FALSE.
       END IF
       if_concentration = if_concentration_out
       if_temperature = if_temperature_out
       CALL read_until(22, '===Should data be interpolated on new mesh? (True/False)')
       READ(22,*) inter_mesh
       if_read_partition = .FALSE.
    END IF

    !===Set rw_ns, rw_mxw and rw_temp (function of Problem type (Input))
    CALL read_until(22, '===Problem type (Input): (nst, mxw, mhd, fhd, mhs)')
    READ(22,*) tit
    IF (tit=='nst') THEN
       rw_ns = .TRUE.
       rw_mxw = .FALSE.
    ELSE IF ((tit=='mhd').OR.(tit=='fhd').OR.(tit=='mhs')) THEN
       rw_ns = .TRUE.
       rw_mxw = .TRUE.
    ELSE IF (tit=='mxw') THEN
       rw_ns = .FALSE.
       rw_mxw = .TRUE.
    ELSE
       CALL error_Petsc('BUG in interpol, type_pb (Input) not correct')
    END IF
    CALL find_string(22, '===Is there an input concentration field?', test)
    IF (test) THEN
       READ (22, *) rw_conc
    ELSE
       rw_conc = .FALSE.
    END IF
    CALL find_string(22, '===Is there an input temperature field?', test)
    IF (test) THEN
       READ (22, *) rw_temp
    ELSE
       rw_temp = .FALSE.
    END IF

    !===Default parameters
    nb_F = 1

    !===Check construction with plotmtv
    CALL find_string(22, '===Check construction with plotmtv? (True/False)', test)
    IF (test) THEN
       READ (22, *) check_plt
    ELSE
       check_plt  = .FALSE.
    END IF

    !===Number of time steps to process
    CALL read_until(22, '===How many files (times steps) should be converted?')
    READ(22, *) nb_fic

    !===Starting index
    CALL find_string(22, '===What is the starting index? (suite*_I*index.mesh*), index is an integer in [1,999]', test)
    IF (test) THEN
       READ (22, *) index_start
    ELSE
       index_start  = 1
    END IF

    !===Number of Fourier modes
    CALL read_until(22, '===Number of Fourier modes')
    READ(22, *) nb_mode
    ALLOCATE(list_mode(nb_mode))
    CALL find_string(22, '===Select Fourier modes? (true/false)',test)
    IF (test) THEN
       READ(22,*) if_select
       IF (if_select) THEN
          CALL read_until(22, '===List of Fourier modes (if select_mode=.TRUE.)')
          READ(22,*) list_mode
       ELSE
          list_mode = (/ (m, m=0, nb_mode-1) /)
       END IF
    ELSE
       list_mode = (/ (m, m=0, nb_mode-1) /)
    END IF

    !===Input and Output mesh files
    CALL read_until(22, '===Directory and name of input mesh file')
    READ(22,*) old_directory, old_filename
    CALL read_until(22, '===Is input mesh file formatted (true/false)?')
    READ(22,*) old_is_form
    IF (is_in) THEN !Interpolation is done at second step only
       new_directory = old_directory
       new_filename = old_filename
       new_is_form = old_is_form
    ELSE
       CALL read_until(22, '===Directory and name of output mesh file')
       READ(22,*) new_directory, new_filename
       CALL read_until(22, '===Is output mesh file formatted (true/false)?')
       READ(22,*) new_is_form
    END IF

     !===Data for concentration
    IF (if_concentration) THEN
       CALL read_until(22, '===Number of subdomains in concentration mesh')
       READ(22,*) nb_dom_conc
       ALLOCATE(list_dom_conc(nb_dom_conc))
       CALL read_until(22, '===List of subdomains for concentration mesh')
       READ (22, *)  list_dom_conc
    ELSE
       ALLOCATE(list_dom_conc(0))
    END IF

    !===Data for NS
    IF ( (type_pb == 'nst') .OR. (type_pb == 'mhd') .OR. (type_pb == 'fhd') &
         .OR. (type_pb == 'mhs')) THEN
       CALL read_until(22, '===Number of subdomains in Navier-Stokes mesh')
       READ(22,*) nb_dom_ns
       ALLOCATE(list_dom_ns_in(nb_dom_ns), list_dom_ns(nb_dom_ns), vel_in_to_new(nb_dom_ns))
       CALL read_until(22, '===List of subdomains for Navier-Stokes mesh')
       READ (22, *)  list_dom_ns_in
       CALL find_string(22, '===Is there a level set?', test)
       IF (test) THEN
          READ (22, *) if_level_set
          IF (if_level_set) THEN
             CALL read_until(22, '===How many fluids?')
             READ(22,*) nb_fluid

             CALL find_string(22, '===Do we use P2 finite element for level_set? (true/false)', test)
             IF (test) THEN
                READ (22, *) if_level_set_P2
             ELSE
                if_level_set_P2=.TRUE.
             END IF
          END IF
       ELSE
          if_level_set = .FALSE.
       END IF
       IF (type_pb == 'nst') ALLOCATE(list_dom_H(0), list_dom_phi(0))

       CALL find_string(22, '===Number of interfaces between velocity and concentration only domains (for nst applications)', test)
       IF (test) THEN
          READ(22,*) nb_inter_c_v
          ALLOCATE(list_inter_c_v(nb_inter_c_v))
          CALL read_until(22, '===List of interfaces between velocity and concentration only domains (for nst applications)')
          READ(22,*) list_inter_c_v
       ELSE
          ALLOCATE(list_inter_c_v(0))
       END IF
    ELSE
       ALLOCATE(list_dom_ns(0))
    END IF

    !===Data for temperature
    IF (if_temperature) THEN
       CALL read_until(22, '===Number of subdomains in temperature mesh')
       READ(22,*) nb_dom_temp
       ALLOCATE(list_dom_temp_in(nb_dom_temp), list_dom_temp(nb_dom_temp), temp_in_to_new(nb_dom_temp))
       CALL read_until(22, '===List of subdomains for temperature mesh')
       READ (22, *)  list_dom_temp_in
       CALL find_string(22, '===Number of interfaces between velocity and temperature only domains (for nst applications)', test)
       IF (test) THEN
          READ(22,*) nb_inter_v_T
          ALLOCATE(list_inter_v_T(nb_inter_v_T))
          CALL read_until(22, '===List of interfaces between velocity and temperature only domains (for nst applications)')
          READ(22,*) list_inter_v_T
       ELSE
          ALLOCATE(list_inter_v_T(0))
       END IF
    ELSE
       ALLOCATE(list_dom_temp(0))
    END IF

    !===Data for Maxwell
    IF ( (type_pb == 'mhd') .OR. (type_pb == 'mxw') .OR. (type_pb == 'fhd') &
         .OR. (type_pb == 'mhs')) THEN
       CALL read_until(22, '===Type of finite element for magnetic field')
       READ (22, *) type_fe_H
       CALL read_until(22, '===Number of subdomains in magnetic field (H) mesh')
       READ(22,*) nb_dom_H  ! number of sub_domains for H
       ALLOCATE(list_dom_H_in(nb_dom_H), list_dom_H(nb_dom_H), H_in_to_new(nb_dom_H)) ! JLG/AR Nov 17 2008
       CALL read_until(22, '===List of subdomains for magnetic field (H) mesh')
       READ(22,*) list_dom_H_in

       !==========Interfaces H/H==========================!
       CALL read_until(22, '===Number of interfaces in H mesh')
       READ(22,*) nb_inter_mu
       ALLOCATE(list_inter_mu(nb_inter_mu))
       IF (nb_inter_mu>0) THEN
          CALL read_until(22, '===List of interfaces in H mesh')
          READ(22,*) list_inter_mu
       END IF

       !==========Subdomains for phi======================!
       CALL find_string(22, '===Number of subdomains in magnetic potential (phi) mesh', test)
       IF (test) THEN
          READ (22, *) nb_dom_phi
       ELSE
          nb_dom_phi = 0
       END IF
       ALLOCATE(list_dom_phi(nb_dom_phi))
       IF (nb_dom_phi>0) THEN
          !==========List of subdomains for phi======================!
          CALL read_until(22, '===List of subdomains for magnetic potential (phi) mesh')
          READ(22,*) list_dom_phi

          !==========H/phi interface=========================!
          CALL read_until(22, '===Number of interfaces between H and phi')
          READ(22,*) nb_inter  ! number of interfaces between H and phi
          ALLOCATE(list_inter_H_phi(nb_inter))
          IF (nb_inter>0) THEN
             CALL read_until(22, '===List of interfaces between H and phi')
             READ(22,*) list_inter_H_phi
          END IF

          !==========Finite element type=====================!
          CALL read_until(22, '===Type of finite element for scalar potential')
          READ(22,*) type_fe_phi
       ELSE
          ALLOCATE(list_inter_H_phi(0))
          type_fe_phi = -1
       END IF
       IF (type_pb == 'mxw') THEN
          ALLOCATE(list_dom_ns(0))
       END IF
    END IF

    CALL find_string(22, '===How many pieces of periodic boundary?', test)
    IF (test) THEN
       READ (22, *) my_periodic%nb_periodic_pairs
    ELSE
       my_periodic%nb_periodic_pairs = 0
    END IF

    IF (my_periodic%nb_periodic_pairs.GE.1) THEN
       ALLOCATE(my_periodic%list_periodic(2,my_periodic%nb_periodic_pairs))
       ALLOCATE(my_periodic%vect_e(2,my_periodic%nb_periodic_pairs))
       CALL read_until(22, '===Indices of periodic boundaries and corresponding vectors')
       DO k = 1, my_periodic%nb_periodic_pairs
          READ(22,*) my_periodic%list_periodic(:,k), my_periodic%vect_e(:,k)
       END DO
    END IF
    CLOSE(22)

    !===Creation of logicals for equations
    if_conc = if_concentration
    if_mass = if_level_set
    if_momentum = (type_pb=='nst' .OR. type_pb=='mhd' .OR. type_pb=='fhd' &
                 .OR. type_pb=='mhs')
    if_induction = (type_pb=='mxw' .OR. type_pb=='mhd' .OR. type_pb=='fhd' &
                 .OR. type_pb=='mhs')
    if_energy = if_temperature

    !===Check mesh that conc_mesh is a subset of H_mesh===============================
    IF (if_induction) THEN
       ALLOCATE(list_dom_H_ref(nb_dom_H), H_in_to_new_ref(nb_dom_H))
       IF (if_conc) THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_conc)) THEN
             CALL error_Petsc('BUG mesh_interpol: conc must be a subset of Maxwell ')
          END IF
          DO k = 1, nb_dom_conc
             IF (MINVAL(ABS(list_dom_H_in - list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: conc must be a subset of Maxwell ')
             END IF
             DO kp = 1, nb_dom_H
                IF (list_dom_H_in(kp) == list_dom_conc(k)) EXIT
             END DO
             H_in_to_new(k) = kp
             list_dom_H(k) = list_dom_conc(k)
          END DO
          m = nb_dom_conc
          DO k = 1, nb_dom_H
             IF (MINVAL(ABS(list_dom_H_in(k) - list_dom_conc)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = k
             list_dom_H(m) = list_dom_H_in(k)
          END DO
          IF (m/=nb_dom_H) THEN
             CALL error_Petsc('BUG in mesh_interpol conc subset H: m/=nb_dom_H ')
          END IF
       ELSE
          DO k = 1, nb_dom_H
             H_in_to_new(k) = k
          END DO
          list_dom_H = list_dom_H_in
       END IF
    ELSE
       ALLOCATE(H_in_to_new_ref(0))
    END IF

    !===Check mesh that conc_mesh is a subset of temp_mesh===============================
    IF (if_energy) THEN
       ALLOCATE(list_dom_temp_ref(nb_dom_temp), temp_in_to_new_ref(nb_dom_temp))
       IF (if_conc) THEN
          IF (SIZE(list_dom_temp) < SIZE(list_dom_conc)) THEN
             CALL error_Petsc('BUG in mesh_interpol: conc must be a subset of temp ')
          END IF
          DO k = 1, nb_dom_conc
             IF (MINVAL(ABS(list_dom_temp_in - list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: conc must be a subset of temp ')
             END IF
             DO kp = 1, nb_dom_temp
                IF (list_dom_temp_in(kp) == list_dom_conc(k)) EXIT
             END DO
             temp_in_to_new(k) = kp
             list_dom_temp(k) = list_dom_conc(k)
          END DO
          m = nb_dom_conc
          DO k = 1, nb_dom_temp
             IF (MINVAL(ABS(list_dom_temp_in(k) - list_dom_conc)) == 0) CYCLE
             m = m + 1
             temp_in_to_new(m) = k
             list_dom_temp(m) = list_dom_temp(k)
          END DO
          IF (m/=nb_dom_temp) THEN
             CALL error_Petsc('BUG in mesh_interpol conc subset T: m/=nb_dom_temp ')
          END IF
       ELSE
          DO k = 1, nb_dom_temp
             temp_in_to_new(k) = k
          END DO
          list_dom_temp = list_dom_temp_in
       END IF
    ELSE
       ALLOCATE(temp_in_to_new(0))
       ALLOCATE(temp_in_to_new_ref(0))
    END IF

    !===Check mesh that conc_mesh is a subset of vv_mesh===============================
    IF (if_momentum) THEN
       IF (if_conc) THEN
          IF (SIZE(list_dom_ns) < SIZE(list_dom_conc)) THEN
             CALL error_Petsc('BUG in mesh_interpol: conc must be a subset of NS ')
          END IF
          DO k = 1, nb_dom_conc
             IF (MINVAL(ABS(list_dom_ns_in - list_dom_conc(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: conc must be a subset of NS ')
             END IF
             DO kp = 1, nb_dom_ns
                IF (list_dom_ns_in(kp) == list_dom_conc(k)) EXIT
             END DO
             vel_in_to_new(k) = kp
             list_dom_ns(k) = list_dom_conc(k)
          END DO
          m = nb_dom_conc
          DO k = 1, nb_dom_ns
             IF (MINVAL(ABS(list_dom_ns_in(k) - list_dom_conc)) == 0) CYCLE
             m = m + 1
             vel_in_to_new(m) = k
             list_dom_ns(m) = list_dom_ns(k)
          END DO
          IF (m/=nb_dom_ns) THEN
             CALL error_Petsc('BUG in mesh_interpol conc subset vel: m/=nb_dom_ns ')
          END IF
       ELSE
          DO k = 1, nb_dom_ns
             vel_in_to_new(k) = k
          END DO
          list_dom_ns = list_dom_ns_in
       END IF
    END IF

    !===Check that vv_mesh is a subset of H_mesh
    IF (if_induction) THEN
       H_in_to_new_ref=H_in_to_new
       list_dom_H_ref=list_dom_H
       IF (if_momentum) THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_ns)) THEN
             CALL error_Petsc('BUG in mesh_interpol: NS must be a subset of Maxwell ')
          END IF
          DO k = 1+nb_dom_conc, nb_dom_ns
             IF (MINVAL(ABS(list_dom_H - list_dom_ns(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: NS must be a subset of Maxwell ')
             END IF
             DO kp = 1+nb_dom_conc, nb_dom_H
                IF (list_dom_H_ref(kp) == list_dom_ns(k)) EXIT
             END DO
             H_in_to_new(k) = H_in_to_new_ref(kp)
             list_dom_H(k) = list_dom_ns(k)
          END DO
          m = nb_dom_ns
          DO k = 1+nb_dom_conc, nb_dom_H
             IF (MINVAL(ABS(list_dom_H_ref(k) - list_dom_ns)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = H_in_to_new_ref(k)
             list_dom_H(m) = list_dom_H_ref(k)
          END DO
          IF (m/=nb_dom_H) THEN
             CALL error_Petsc('BUG in mesh_interpol vel subset H: m/=nb_dom_H ')
          END IF
       END IF
    END IF

    !===Check mesh that vv_mesh is a subset of temp_mesh============================
    IF (if_energy) THEN
       temp_in_to_new_ref=temp_in_to_new
       list_dom_temp_ref=list_dom_temp
       IF (if_momentum) THEN
          IF (SIZE(list_dom_temp) < SIZE(list_dom_ns)) THEN
             CALL error_Petsc('BUG in mesh_interpol: NS must be a subset of temp ')
          END IF
          DO k = 1+nb_dom_conc, nb_dom_ns
             IF (MINVAL(ABS(list_dom_temp - list_dom_ns(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: NS must be a subset of temp ')
             END IF
             DO kp = 1+nb_dom_conc, nb_dom_temp
                IF (list_dom_temp_ref(kp) == list_dom_ns(k)) EXIT
             END DO
             temp_in_to_new(k) = temp_in_to_new_ref(kp)
             list_dom_temp(k) = list_dom_ns(k)
          END DO
          m = nb_dom_ns
          DO k = 1+nb_dom_conc, nb_dom_temp
             IF (MINVAL(ABS(list_dom_temp_ref(k) - list_dom_ns)) == 0) CYCLE
             m = m + 1
             temp_in_to_new(m) = temp_in_to_new_ref(k)
             list_dom_temp(m) = list_dom_temp_ref(k)
          END DO
          IF (m/=nb_dom_temp) THEN
             CALL error_Petsc('BUG in mesh_interpol vel subset temp: m/=nb_dom_temp ')
          END IF
       END IF
    END IF

    !===Check that temp_mesh is a subset of H_mesh
    IF (if_induction) THEN
       H_in_to_new_ref=H_in_to_new
       list_dom_H_ref=list_dom_H
       IF (if_energy) THEN
          IF (SIZE(list_dom_H) < SIZE(list_dom_temp)) THEN
             CALL error_Petsc('BUG in mesh_interpol: temp must be a subset of H ')
          END IF
          DO k = 1+nb_dom_ns, nb_dom_temp
             IF (MINVAL(ABS(list_dom_H - list_dom_temp(k))) /= 0) THEN
                CALL error_Petsc('BUG in mesh_interpol: temp must be a subset of H ')
             END IF
             DO kp = 1+nb_dom_ns, nb_dom_H
                IF (list_dom_H_ref(kp) == list_dom_temp(k)) EXIT
             END DO
             H_in_to_new(k) = H_in_to_new_ref(kp)
             list_dom_H(k) = list_dom_temp(k)
          END DO
          m = nb_dom_temp
          DO k = 1+nb_dom_ns, nb_dom_H
             IF (MINVAL(ABS(list_dom_H_ref(k) - list_dom_temp)) == 0) CYCLE
             m = m + 1
             H_in_to_new(m) = H_in_to_new_ref(k)
             list_dom_H(m) = list_dom_H_ref(k)
          END DO
          IF (m/=nb_dom_H) THEN
             CALL error_Petsc('BUG in mesh_interpol temp subset H: m/=nb_dom_H ')
          END IF
       END IF
    END IF

    !===Create interfaces in meshes
    IF (if_momentum .AND. (.NOT. if_induction)) THEN
       IF (if_energy) THEN
          nsize = SIZE(list_dom_temp)
          ALLOCATE(list_dom(nsize))
          list_dom = list_dom_temp
          ALLOCATE(list_inter(SIZE(list_inter_v_T)))
          list_inter = list_inter_v_T
       ELSE
          nsize = SIZE(list_dom_ns)
          ALLOCATE(list_dom(nsize))
          list_dom = list_dom_ns
          IF (if_conc) THEN
             ALLOCATE(list_inter(SIZE(list_inter_c_v)))
             list_inter = list_inter_c_v
          ELSE
             ALLOCATE(list_inter(0))
          END IF
       END IF
    ELSE
       nsize = SIZE(list_dom_H)+SIZE(list_dom_phi)
       ALLOCATE(list_dom(nsize))
       list_dom(1:SIZE(list_dom_H)) = list_dom_H
       IF (SIZE(list_dom_phi)>0) THEN
          list_dom(SIZE(list_dom_H)+1:) = list_dom_phi
       END IF
       nsize = SIZE(list_inter_mu)+SIZE(list_inter_H_phi)
       ALLOCATE(list_inter(nsize))
       IF (SIZE(list_inter_mu)>0) THEN
          list_inter(1:SIZE(list_inter_mu)) = list_inter_mu
       END IF
       IF (SIZE(list_inter_H_phi)>0) THEN
          list_inter(SIZE(list_inter_mu)+1:) = list_inter_H_phi
       END IF
    END IF
    IF (if_energy) THEN
       ALLOCATE(list_inter_temp(0))
    END IF
    IF (if_conc) THEN
       ALLOCATE(list_inter_conc(0))
    END IF

    !===Directory, file name and format
    IF (is_in) THEN
       directory  = old_directory
       file_name  = old_filename
       iformatted = old_is_form
    ELSE
       directory  = new_directory
       file_name  = new_filename
       iformatted = new_is_form
    END IF

    IF (is_in) THEN
       directory_m = new_directory
       file_name_m = new_filename
       is_form_m   = new_is_form
    ELSE
       directory_m = old_directory
       file_name_m = old_filename
       is_form_m   = old_is_form
    END IF

    !===Check coherence
    CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, nb_procs, ierr)
    IF ( nb_S*nb_F /= nb_procs ) THEN
       CALL error_Petsc('BUG: nb_S * nb_F /= nb_procs')
    END IF

    CALL create_cart_comm((/nb_S, nb_F/),comm_cart,comm_one_d,coord_cart)
    CALL MPI_COMM_SIZE(comm_one_d(2),nb_procs,code)
    CALL MPI_COMM_RANK(comm_one_d(2),rank,code)
    IF (nb_F/=nb_procs) THEN
       CALL error_Petsc('BUG in INIT, nb_procs_F/=nb_procs')
    END IF

    IF ( (rw_ns) .AND. (type_pb=='mxw')) THEN
       CALL error_Petsc('WARNING: type_pb and/or rw_ns might be wrong')
    END IF
    IF ( (rw_mxw) .AND. ((type_pb=='nst')) ) THEN
       CALL error_Petsc('WARNING: type_pb and/or rw_mxw might be wrong')
    END IF

    !===Prepare meshes and pointers
    CALL load_dg_mesh_free_format(directory, file_name, list_dom, list_inter, 1, p1_mesh_glob, iformatted)
    CALL load_dg_mesh_free_format(directory, file_name, list_dom, list_inter, 2, p2_mesh_glob, iformatted)
    IF (if_conc) THEN
!       CALL load_dg_mesh_free_format(directory, file_name, list_dom, &
       CALL load_dg_mesh_free_format(directory, file_name, list_dom_conc, & !TEST LC
            list_inter_conc, 2, p2_c0_mesh_glob_conc, iformatted)
    END IF
    IF (if_energy) THEN
!       CALL load_dg_mesh_free_format(directory, file_name, list_dom, &
       CALL load_dg_mesh_free_format(directory, file_name, list_dom_temp, & !TEST LC
            list_inter_temp, 2, p2_c0_mesh_glob_temp, iformatted)
    END IF

    !===Start Metis mesh generation=================================================
    ALLOCATE(part(p1_mesh_glob%me))
    WRITE(tit_part,'(i4)') nb_S
    mesh_part_name='mesh_part_S'//TRIM(ADJUSTL(tit_part))//'.'//TRIM(ADJUSTL(file_name))
    IF (if_read_partition) THEN
       OPEN(UNIT=51, FILE=mesh_part_name, STATUS='unknown', FORM='formatted')
       READ(51,*) part
       CLOSE(51)
       WRITE(*,*) 'read partition'
    ELSE
       WRITE(*,*) 'create partition'
       CALL part_mesh_M_T_H_phi(nb_S, list_dom_conc, list_dom_ns, list_dom_temp, list_dom_H, &
            list_dom_phi, p1_mesh_glob, list_inter, part, my_periodic)
       IF (petsc_rank==0) THEN
          OPEN(UNIT=51, FILE=mesh_part_name, STATUS='replace', FORM='formatted')
          WRITE(51,*) part
          CLOSE(51)
       END IF
    END IF

    !===Extract local meshes from global meshes
    IF (if_conc) THEN
       CALL extract_mesh(comm_one_d(1),nb_S,p2_c0_mesh_glob_conc,part,list_dom_conc,conc_mesh_glob,conc_mesh)
       ALLOCATE(comm_one_d_conc(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_conc(2), code)
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (conc_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_conc(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_conc(1),code)
       END IF
    END IF

    IF (if_momentum) THEN
       CALL extract_mesh(comm_one_d(1),nb_S,p1_mesh_glob,part,list_dom_ns,pp_mesh_glob,pp_mesh)
       CALL extract_mesh(comm_one_d(1),nb_S,p2_mesh_glob,part,list_dom_ns,vv_mesh_glob,vv_mesh)
       ALLOCATE(comm_one_d_ns(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_ns(2), code)
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (pp_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_ns(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_ns(1),code)
       END IF
    END IF

    IF (if_energy) THEN
       CALL extract_mesh(comm_one_d(1),nb_S,p2_c0_mesh_glob_temp,part,list_dom_temp,temp_mesh_glob,temp_mesh)
       ALLOCATE(comm_one_d_temp(2))
       CALL MPI_COMM_DUP(comm_one_d(2), comm_one_d_temp(2), code)
       CALL MPI_COMM_RANK(comm_one_d(1),rank_S,code)
       IF (temp_mesh%me/=0) THEN
          CALL MPI_COMM_SPLIT (comm_one_d(1),1,rank_S,comm_one_d_temp(1),code)
       ELSE
          CALL MPI_COMM_SPLIT (comm_one_d(1),MPI_UNDEFINED,rank_S,comm_one_d_temp(1),code)
       END IF
    END IF

    IF (if_induction) THEN
       IF (type_fe_H==1) THEN
          CALL extract_mesh(comm_one_d(1),nb_S,p1_mesh_glob,part,list_dom_H,H_mesh_glob,H_mesh)
       ELSE
          CALL extract_mesh(comm_one_d(1),nb_S,p2_mesh_glob,part,list_dom_H,H_mesh_glob,H_mesh)
       END IF
       IF (type_fe_phi==1) THEN
          CALL extract_mesh(comm_one_d(1),nb_S,p1_mesh_glob,part,list_dom_phi,phi_mesh_glob,phi_mesh)
       ELSE
          CALL extract_mesh(comm_one_d(1),nb_S,p2_mesh_glob,part,list_dom_phi,phi_mesh_glob,phi_mesh)
       END IF
    END IF

    !===Cleanup
    CALL free_mesh(p1_mesh_glob)
    CALL free_mesh(p2_mesh_glob)
    IF (if_conc) THEN
       DEALLOCATE(list_inter_conc)
       CALL free_mesh(p2_c0_mesh_glob_conc)
    END IF
    IF (if_energy) THEN
       DEALLOCATE(list_inter_temp)
       CALL free_mesh(p2_c0_mesh_glob_temp)
    END IF

    m_max_c = nb_mode/nb_F

    !===Load meshes for monoproc
    IF (if_conc) THEN
       CALL free_mesh(conc_mesh_glob)
       CALL load_mesh_free_format(directory_m, file_name_m, list_dom_conc,  2, conc_mesh_glob, is_form_m)
       IF (check_plt) THEN
          CALL plot_const_p1_label(conc_mesh_glob%jj, conc_mesh_glob%rr, 1.d0*conc_mesh_glob%i_d, 'conc.plt')
       END IF
    END IF
    IF (if_momentum) THEN
       CALL free_mesh(vv_mesh_glob)
       CALL free_mesh(pp_mesh_glob)
       CALL load_mesh_free_format(directory_m, file_name_m, list_dom_ns,  2, vv_mesh_glob, is_form_m)
       CALL load_mesh_free_format(directory_m, file_name_m, list_dom_ns,  1, pp_mesh_glob, is_form_m)
       IF (check_plt) THEN
          CALL plot_const_p1_label(vv_mesh_glob%jj, vv_mesh_glob%rr, 1.d0*vv_mesh_glob%i_d, 'vv.plt')
       END IF
    END IF
    IF (if_energy) THEN
       CALL free_mesh(temp_mesh_glob)
       CALL load_mesh_free_format(directory_m, file_name_m, list_dom_temp,  2, temp_mesh_glob, is_form_m)
       IF (check_plt) THEN
          CALL plot_const_p1_label(temp_mesh_glob%jj, temp_mesh_glob%rr, 1.d0*temp_mesh_glob%i_d, 'temp.plt')
       END IF
    END IF
    IF (if_induction) THEN
       CALL free_mesh(H_mesh_glob)
       CALL free_mesh(phi_mesh_glob)
       CALL load_dg_mesh_free_format(directory_m, file_name_m, list_dom_H, list_inter_mu,  type_fe_H, H_mesh_glob, is_form_m)
       CALL load_mesh_free_format(directory_m, file_name_m, list_dom_phi, type_fe_phi, phi_mesh_glob, is_form_m)
       IF (check_plt) THEN
          CALL plot_const_p1_label(H_mesh_glob%jj, H_mesh_glob%rr, 1.d0*H_mesh_glob%i_d, 'HH.plt')
          CALL plot_const_p1_label(phi_mesh_glob%jj, phi_mesh_glob%rr, 1.d0*phi_mesh_glob%i_d, 'phi.plt')
       END IF
    END IF

    !===Array allocation
    IF (if_conc) THEN
       ALLOCATE(concn_m1_glob(conc_mesh_glob%np, 2, m_max_c))
       ALLOCATE(concn_glob   (conc_mesh_glob%np, 2, m_max_c))
       ALLOCATE(concn_m1     (conc_mesh%np, 2, m_max_c))
       ALLOCATE(concn        (conc_mesh%np, 2, m_max_c))
       concn_m1_glob = 0.d0
       concn_glob    = 0.d0
       concn_m1      = 0.d0
       concn         = 0.d0
    END IF

    IF (if_momentum) THEN
       ALLOCATE(un_glob      (vv_mesh_glob%np, 6, m_max_c))
       ALLOCATE(un_m1_glob   (vv_mesh_glob%np, 6, m_max_c))
       ALLOCATE(pn_glob      (pp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(pn_m1_glob   (pp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(incpn_m1_glob(pp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(incpn_glob   (pp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(un      (vv_mesh%np, 6, m_max_c))
       ALLOCATE(un_m1   (vv_mesh%np, 6, m_max_c))
       ALLOCATE(pn      (pp_mesh%np, 2, m_max_c))
       ALLOCATE(pn_m1   (pp_mesh%np, 2, m_max_c))
       ALLOCATE(incpn_m1(pp_mesh%np, 2, m_max_c))
       ALLOCATE(incpn   (pp_mesh%np, 2, m_max_c))
       un_glob       = 0.d0
       un_m1_glob    = 0.d0
       pn_glob       = 0.d0
       pn_m1_glob    = 0.d0
       incpn_m1_glob = 0.d0
       incpn_glob    = 0.d0
       un      = 0.d0
       un_m1   = 0.d0
       pn      = 0.d0
       pn_m1   = 0.d0
       incpn_m1= 0.d0
       incpn   = 0.d0
       IF (if_mass) THEN
          IF (if_level_set_P2) THEN
             ALLOCATE(level_setn_m1_glob(nb_fluid-1,vv_mesh_glob%np, 2, m_max_c))
             ALLOCATE(level_setn_glob   (nb_fluid-1,vv_mesh_glob%np, 2, m_max_c))
             ALLOCATE(level_setn_m1     (nb_fluid-1,vv_mesh%np, 2, m_max_c))
             ALLOCATE(level_setn        (nb_fluid-1,vv_mesh%np, 2, m_max_c))
          ELSE
             ALLOCATE(level_setn_m1_glob(nb_fluid-1,pp_mesh_glob%np, 2, m_max_c))
             ALLOCATE(level_setn_glob   (nb_fluid-1,pp_mesh_glob%np, 2, m_max_c))
             ALLOCATE(level_setn_m1     (nb_fluid-1,pp_mesh%np, 2, m_max_c))
             ALLOCATE(level_setn        (nb_fluid-1,pp_mesh%np, 2, m_max_c))
          END IF
          level_setn_m1_glob = 0.d0
          level_setn_glob    = 0.d0
          level_setn_m1      = 0.d0
          level_setn         = 0.d0
       END IF
    END IF

    IF (if_energy) THEN
       ALLOCATE(tempn_m1_glob(temp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(tempn_glob   (temp_mesh_glob%np, 2, m_max_c))
       ALLOCATE(tempn_m1     (temp_mesh%np, 2, m_max_c))
       ALLOCATE(tempn        (temp_mesh%np, 2, m_max_c))
       tempn_m1_glob = 0.d0
       tempn_glob    = 0.d0
       tempn_m1      = 0.d0
       tempn         = 0.d0
    END IF

    IF (if_induction) THEN
       ALLOCATE(Hn1  (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Hn   (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Bn1  (H_mesh%np,  6,  m_max_c))
       ALLOCATE(Bn   (H_mesh%np,  6,  m_max_c))
       ALLOCATE(phin1(phi_mesh%np,2,  m_max_c))
       ALLOCATE(phin (phi_mesh%np,2,  m_max_c))
       ALLOCATE(Hn1_glob   (H_mesh_glob%np,  6,  m_max_c))
       ALLOCATE(Hn_glob    (H_mesh_glob%np,  6,  m_max_c))
       ALLOCATE(Bn1_glob   (H_mesh_glob%np,  6,  m_max_c))
       ALLOCATE(Bn_glob    (H_mesh_glob%np,  6,  m_max_c))
       ALLOCATE(phin1_glob (phi_mesh_glob%np,2,  m_max_c))
       ALLOCATE(phin_glob  (phi_mesh_glob%np,2,  m_max_c))
       Hn1        = 0.d0
       Hn         = 0.d0
       Bn1        = 0.d0
       Bn         = 0.d0
       phin1      = 0.d0
       phin       = 0.d0
       Hn1_glob   = 0.d0
       Hn_glob    = 0.d0
       Bn1_glob   = 0.d0
       Bn_glob    = 0.d0
       phin1_glob = 0.d0
       phin_glob  = 0.d0
    END IF

    !===Pointers
    IF (is_in) THEN
       mono_in = .FALSE.
       mono_out = .TRUE.
       IF (if_induction) THEN
          Hn_in     => Hn
          Hn1_in    => Hn1
          Bn_in     => Bn
          Bn1_in    => Bn1
          phin_in   => phin
          phin1_in  => phin1
          H_mesh_in => H_mesh
          phi_mesh_in => phi_mesh
          Hn_out     => Hn_glob
          Hn1_out    => Hn1_glob
          Bn_out     => Bn_glob
          Bn1_out    => Bn1_glob
          phin_out   => phin_glob
          phin1_out  => phin1_glob
          H_mesh_out => H_mesh_glob
          phi_mesh_out => phi_mesh_glob
       END IF
       IF (if_momentum) THEN
          un_in      => un
          un_m1_in   => un_m1
          pn_in      => pn
          pn_m1_in   => pn_m1
          incpn_in   => incpn
          incpn_m1_in=> incpn_m1
          vv_mesh_in => vv_mesh
          pp_mesh_in => pp_mesh
          un_out      => un_glob
          un_m1_out   => un_m1_glob
          pn_out      => pn_glob
          pn_m1_out   => pn_m1_glob
          incpn_out   => incpn_glob
          incpn_m1_out=> incpn_m1_glob
          vv_mesh_out => vv_mesh_glob
          pp_mesh_out => pp_mesh_glob
          IF (if_mass) THEN
             level_setn_in    => level_setn
             level_setn_m1_in => level_setn_m1
             level_setn_out   => level_setn_glob
             level_setn_m1_out=> level_setn_m1_glob
          END IF
       END IF
       IF (if_energy) THEN
          tempn_in    => tempn
          tempn_m1_in => tempn_m1
          tempn_out   => tempn_glob
          tempn_m1_out => tempn_m1_glob
          temp_mesh_in => temp_mesh
          temp_mesh_out => temp_mesh_glob
       END IF
       IF (if_conc) THEN
          concn_in    => concn
          concn_m1_in => concn_m1
          concn_out   => concn_glob
          concn_m1_out => concn_m1_glob
          conc_mesh_in => conc_mesh
          conc_mesh_out => conc_mesh_glob
       END IF

    ELSE
       mono_in = .TRUE.
       mono_out = .FALSE.
       IF (if_induction) THEN
          Hn_in     => Hn_glob
          Hn1_in    => Hn1_glob
          Bn_in     => Bn_glob
          Bn1_in    => Bn1_glob
          phin_in   => phin_glob
          phin1_in  => phin1_glob
          H_mesh_in => H_mesh_glob
          phi_mesh_in => phi_mesh_glob
          Hn_out     => Hn
          Hn1_out    => Hn1
          Bn_out     => Bn
          Bn1_out    => Bn1
          phin_out   => phin
          phin1_out  => phin1
          H_mesh_out => H_mesh
          phi_mesh_out => phi_mesh
       END IF
       IF (if_momentum) THEN
          un_in      => un_glob
          un_m1_in   => un_m1_glob
          pn_in      => pn_glob
          pn_m1_in   => pn_m1_glob
          incpn_in   => incpn_glob
          incpn_m1_in=> incpn_m1_glob
          vv_mesh_in => vv_mesh_glob
          pp_mesh_in => pp_mesh_glob
          un_out      => un
          un_m1_out   => un_m1
          pn_out      => pn
          pn_m1_out   => pn_m1
          incpn_out   => incpn
          incpn_m1_out=> incpn_m1
          vv_mesh_out => vv_mesh
          pp_mesh_out => pp_mesh
          IF (if_mass) THEN
             level_setn_in    => level_setn_glob
             level_setn_m1_in => level_setn_m1_glob
             level_setn_out   => level_setn
             level_setn_m1_out=> level_setn_m1
          END IF
       END IF
       IF (if_energy) THEN
          tempn_in    => tempn_glob
          tempn_m1_in => tempn_m1_glob
          tempn_out   => tempn
          tempn_m1_out=> tempn_m1
          temp_mesh_in => temp_mesh_glob
          temp_mesh_out => temp_mesh
       END IF
       IF (if_conc) THEN
          concn_in    => concn_glob
          concn_m1_in => concn_m1_glob
          concn_out   => concn
          concn_m1_out=> concn_m1
          conc_mesh_in => conc_mesh_glob
          conc_mesh_out => conc_mesh
       END IF
    END IF

    !===Interpolation for Maxwell
    IF (rw_mxw) THEN
       IF (rank==0) WRITE(*,*) 'Start interpolation Maxwell'
       IF (inter_mesh) THEN
          ALLOCATE(controle_H(H_mesh_out%np), controle_phi(phi_mesh_out%np))
          DO m = index_start, index_start+nb_fic-1
             Hn1        = 0.d0
             Hn         = 0.d0
             Bn1        = 0.d0
             Bn         = 0.d0
             phin1      = 0.d0
             phin       = 0.d0
             Hn1_glob   = 0.d0
             Hn_glob    = 0.d0
             Bn1_glob   = 0.d0
             Bn_glob    = 0.d0
             phin1_glob = 0.d0
             phin_glob  = 0.d0
             WRITE(tit, '(i3)') m
             lblank = eval_blank(3,tit)
             DO l = 1, lblank - 1
                tit(l:l) = '0'
             END DO

             IF (petsc_rank==0) CALL system('mv suite_maxwell_I'//tit//'.'//old_filename//'suite_maxwell.'//old_filename)
             CALL MPI_Barrier( MPI_Comm_WORLD, code)
             CALL read_restart_maxwell(comm_one_d, H_mesh_in, phi_mesh_in, time_h, list_mode, Hn_in, Hn1_in, &
                  Bn_in, Bn1_in, phin_in, phin1_in, old_filename, interpol=.FALSE., opt_mono=mono_in)

             !===Controle_H and controle_phi are initialized to zero in interp_mesh
             CALL interp_mesh(H_mesh_in, H_mesh_out, Hn_in, Hn_out, controle_H, type_fe_H)
             CALL interp_mesh(H_mesh_in, H_mesh_out, Hn1_in, Hn1_out, controle_H, type_fe_H)
             CALL interp_mesh(H_mesh_in, H_mesh_out, Bn_in, Bn_out, controle_H, type_fe_H)
             CALL interp_mesh(H_mesh_in, H_mesh_out, Bn1_in, Bn1_out, controle_H, type_fe_H)
             CALL interp_mesh(phi_mesh_in, phi_mesh_out, phin_in, phin_out, controle_phi, type_fe_phi)
             CALL interp_mesh(phi_mesh_in, phi_mesh_out, phin1_in, phin1_out, controle_phi, type_fe_phi)

             IF ( MIN(MINVAL(controle_H),MINVAL(controle_phi)) == 0) THEN
                CALL error_Petsc('certains points non trouve H/phi 2')
             END IF

             CALL write_restart_maxwell(comm_one_d, H_mesh_out, phi_mesh_out, time_h, list_mode, Hn_out, Hn1_out, &
                  Bn_out, Bn1_out, phin_out, phin1_out, new_filename, m, 1, opt_mono=mono_out)
             CALL MPI_COMM_RANK(comm_one_d(1), rank_S, ierr)
          END DO
          DEALLOCATE(controle_H, controle_phi)
       ELSE
          ALLOCATE(l_t_g_H(H_mesh%np), l_t_g_phi(phi_mesh%np))
          l_t_g_H = 0
          l_t_g_phi = 0
          CALL loc_to_glob(H_mesh, H_mesh_glob, l_t_g_H)
          CALL loc_to_glob(phi_mesh, phi_mesh_glob, l_t_g_phi)
          DO m = index_start, index_start+nb_fic-1
             Hn1        = 0.d0
             Hn         = 0.d0
             Bn1        = 0.d0
             Bn         = 0.d0
             phin1      = 0.d0
             phin       = 0.d0
             Hn1_glob   = 0.d0
             Hn_glob    = 0.d0
             Bn1_glob   = 0.d0
             Bn_glob    = 0.d0
             phin1_glob = 0.d0
             phin_glob  = 0.d0
             WRITE(tit, '(i3)') m
             lblank = eval_blank(3,tit)
             DO l = 1, lblank - 1
                tit(l:l) = '0'
             END DO

             IF (is_in) THEN
                CALL MPI_COMM_RANK(comm_one_d(1), rank_S, ierr)
                WRITE(tit_s,'(i3)') rank_S
                lblank = eval_blank(3,tit_s)
                DO l = 1, lblank - 1
                   tit_s(l:l) = '0'
                END DO
                CALL system('mv suite_maxwell_S'//tit_s//'_I'//tit//'.'//old_filename//'suite_maxwell_S'//tit_s//'.'//old_filename)
             ELSE
                IF (petsc_rank==0) CALL system('mv suite_maxwell_I'//tit//'.'//old_filename//'suite_maxwell.'//old_filename)
                CALL MPI_Barrier( MPI_Comm_WORLD, code)
             END IF
             CALL read_restart_maxwell(comm_one_d, H_mesh_in, phi_mesh_in, time_h, list_mode, Hn_in, Hn1_in, &
                  Bn_in, Bn1_in, phin_in, phin1_in, old_filename, interpol=.FALSE., opt_mono=mono_in)

             CALL inter_mesh_loc_to_glob(H_mesh_in, H_mesh_out, Hn_in, Hn_out, l_t_g_H, is_in, comm_one_d(1))
             CALL inter_mesh_loc_to_glob(H_mesh_in, H_mesh_out, Hn1_in, Hn1_out, l_t_g_H, is_in, comm_one_d(1))
             CALL inter_mesh_loc_to_glob(H_mesh_in, H_mesh_out, Bn_in, Bn_out, l_t_g_H, is_in, comm_one_d(1))
             CALL inter_mesh_loc_to_glob(H_mesh_in, H_mesh_out, Bn1_in, Bn1_out, l_t_g_H, is_in, comm_one_d(1))
             CALL inter_mesh_loc_to_glob(phi_mesh_in, phi_mesh_out, phin_in, phin_out, l_t_g_phi, is_in, comm_one_d(1))
             CALL inter_mesh_loc_to_glob(phi_mesh_in, phi_mesh_out, phin1_in, phin1_out, l_t_g_phi, is_in, comm_one_d(1))

             CALL write_restart_maxwell(comm_one_d, H_mesh_out, phi_mesh_out, time_h, list_mode, Hn_out, Hn1_out, &
                  Bn_out, Bn1_out, phin_out, phin1_out, new_filename, m, 1, opt_mono=mono_out)
             CALL MPI_COMM_RANK(comm_one_d(1), rank_S, ierr)
          END DO

       END IF

       IF (check_plt) THEN
          IF (H_mesh%me /= 0) THEN
             CALL MPI_COMM_RANK(comm_one_d(1), rang_S, ierr)
             WRITE(tit_s,'(i3)') rang_S
             lblank = eval_blank(3,tit_s)
             DO l = 1, lblank - 1
                tit_s(l:l) = '0'
             END DO
             DO i = 1, SIZE(list_mode)
                WRITE(tit_m,'(i3)') list_mode(i)
                lblank = eval_blank(3,tit_m)
                DO l = 1, lblank - 1
                   tit_m(l:l) = '0'
                END DO
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,1,i), 'H_r_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,2,i), 'H_r_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,3,i), 'H_t_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,4,i), 'H_t_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,5,i), 'H_z_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(H_mesh%jj, H_mesh%rr, Hn(:,6,i), 'H_z_sin_m='//tit_m//'_'//tit_s//'_999.plt' )

             END DO

             IF (rang_S == 0) THEN
                DO i = 1, SIZE(list_mode)
                   WRITE(tit_m,'(i3)') list_mode(i)
                   lblank = eval_blank(3,tit_m)
                   DO l = 1, lblank - 1
                      tit_m(l:l) = '0'
                   END DO
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,1,i), 'gH_r_cos_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,2,i), 'gH_r_sin_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,3,i), 'gH_t_cos_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,4,i), 'gH_t_sin_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,5,i), 'gH_z_cos_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(H_mesh_glob%jj, H_mesh_glob%rr, Hn_glob(:,6,i), 'gH_z_sin_m='//tit_m//'_999.plt' )
                END DO
             END IF
          END IF
       END IF
       IF (rank==0) WRITE(*,*) 'End interpolation Maxwell'

    END IF

    !===Interpolation for NS
    IF (rw_ns) THEN
       IF (rank==0) WRITE(*,*) 'Start interpolation NS'

       IF (inter_mesh) THEN

          ALLOCATE(controle_vv(vv_mesh_out%np), controle_pp(pp_mesh_out%np))
          DO m = index_start, index_start+nb_fic-1
             un_glob       = 0.d0
             un_m1_glob    = 0.d0
             pn_glob       = 0.d0
             pn_m1_glob    = 0.d0
             incpn_m1_glob = 0.d0
             incpn_glob    = 0.d0
             un      = 0.d0
             un_m1   = 0.d0
             pn      = 0.d0
             pn_m1   = 0.d0
             incpn_m1= 0.d0
             incpn   = 0.d0
             IF (if_level_set) THEN
                level_setn_m1_glob = 0.d0
                level_setn_glob    = 0.d0
                level_setn_m1      = 0.d0
                level_setn         = 0.d0
             END IF

             WRITE(tit, '(i3)') m
             lblank = eval_blank(3,tit)
             DO l = 1, lblank - 1
                tit(l:l) = '0'
             END DO

             IF (petsc_rank==0) THEN
                CALL system('mv suite_ns_I'//tit//'.'//old_filename//'suite_ns.'//old_filename)
             END IF
             CALL MPI_Barrier( MPI_Comm_WORLD, code)

             IF (vv_mesh%me/=0) THEN
                IF (if_level_set)THEN
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un_in, un_m1_in, &
                        pn_in, pn_m1_in, incpn_in, incpn_m1_in, old_filename,opt_level_set=level_setn_in, &
                        opt_level_set_m1=level_setn_m1_in, opt_max_vel=max_vel, opt_mono = mono_in)

                ELSE
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un_in, un_m1_in, &
                        pn_in, pn_m1_in, incpn_in, incpn_m1_in, old_filename, opt_mono = mono_in)
                END IF
             END IF
             controle_vv = 0
             controle_pp = 0
             CALL interp_mesh(vv_mesh_in, vv_mesh_out, un_in, un_out, controle_vv, 2)
             CALL interp_mesh(vv_mesh_in, vv_mesh_out, un_m1_in, un_m1_out, controle_vv, 2)
             CALL interp_mesh(pp_mesh_in, pp_mesh_out, pn_in, pn_out, controle_pp, 1)
             CALL interp_mesh(pp_mesh_in, pp_mesh_out, pn_m1_in, pn_m1_out, controle_pp, 1)
             CALL interp_mesh(pp_mesh_in, pp_mesh_out, incpn_in, incpn_out, controle_pp, 1)
             CALL interp_mesh(pp_mesh_in, pp_mesh_out, incpn_m1_in, incpn_m1_out, controle_pp, 1)
             IF (if_level_set) THEN
                IF (if_level_set_P2) THEN
                   DO k = 1, nb_fluid-1
                      CALL interp_mesh(vv_mesh_in, vv_mesh_out, level_setn_in(k,:,:,:),&
                           level_setn_out(k,:,:,:), controle_vv, 2)
                      CALL interp_mesh(vv_mesh_in, vv_mesh_out, level_setn_m1_in(k,:,:,:),&
                           level_setn_m1_out(k,:,:,:), controle_vv, 2)
                   END DO
                ELSE
                   DO k = 1, nb_fluid-1
                      CALL interp_mesh(pp_mesh_in, pp_mesh_out, level_setn_in(k,:,:,:),&
                           level_setn_out(k,:,:,:), controle_pp, 1)
                      CALL interp_mesh(pp_mesh_in, pp_mesh_out, level_setn_m1_in(k,:,:,:),&
                           level_setn_m1_out(k,:,:,:), controle_pp, 1)
                   END DO
                END IF
             END IF

             IF ( MIN(MINVAL(controle_vv),MINVAL(controle_pp)) == 0) THEN
                CALL error_Petsc('certains points non trouve')
             END IF

             IF (pp_mesh%me /=0) THEN
                IF (if_level_set) THEN
                   CALL write_restart_ns(comm_one_d_ns, vv_mesh_out, pp_mesh_out, time_u, list_mode, un_out, un_m1_out, &
                        pn_out, pn_m1_out, incpn_out, incpn_m1_out, new_filename, m, 1, opt_level_set=level_setn_out, &
                        opt_level_set_m1=level_setn_m1_out, opt_max_vel=max_vel, opt_mono=mono_out)
                ELSE
                   CALL write_restart_ns(comm_one_d_ns, vv_mesh_out, pp_mesh_out, time_u, list_mode, un_out, un_m1_out, &
                        pn_out, pn_m1_out, incpn_out, incpn_m1_out, new_filename, m, 1, opt_mono=mono_out)
                END IF
                CALL MPI_COMM_RANK(comm_one_d_ns(1), rang_ns_S, ierr)
             END IF
          END DO
          DEALLOCATE(controle_vv, controle_pp)

       ELSE

          ALLOCATE(l_t_g_vv(vv_mesh%np), l_t_g_pp(pp_mesh%np))
          l_t_g_vv = 0
          l_t_g_pp = 0
          CALL loc_to_glob(vv_mesh, vv_mesh_glob, l_t_g_vv)
          CALL loc_to_glob(pp_mesh, pp_mesh_glob, l_t_g_pp)
          IF (pp_mesh%me /=0) THEN
             DO m = index_start, index_start+nb_fic-1
                un_glob       = 0.d0
                un_m1_glob    = 0.d0
                pn_glob       = 0.d0
                pn_m1_glob    = 0.d0
                incpn_m1_glob = 0.d0
                incpn_glob    = 0.d0
                un      = 0.d0
                un_m1   = 0.d0
                pn      = 0.d0
                pn_m1   = 0.d0
                incpn_m1= 0.d0
                incpn   = 0.d0
                IF (if_level_set) THEN
                   level_setn_m1_glob = 0.d0
                   level_setn_glob    = 0.d0
                   level_setn_m1      = 0.d0
                   level_setn         = 0.d0
                END IF

                WRITE(tit, '(i3)') m
                lblank = eval_blank(3,tit)
                DO l = 1, lblank - 1
                   tit(l:l) = '0'
                END DO
                IF (is_in) THEN
                   CALL MPI_COMM_RANK(comm_one_d_ns(1), rang_ns_S, ierr)
                   WRITE(tit_s,'(i3)') rang_ns_S
                   lblank = eval_blank(3,tit_s)
                   DO l = 1, lblank - 1
                      tit_s(l:l) = '0'
                   END DO

                   CALL system('mv suite_ns_S'//tit_s//'_I'//tit//'.'//old_filename//'suite_ns_S'//tit_s//'.'//old_filename)

                ELSE
                   IF (petsc_rank==0) CALL system('mv suite_ns_I'//tit//'.'//old_filename//'suite_ns.'//old_filename)
                   CALL MPI_Barrier( MPI_Comm_WORLD, code)
                END IF

                IF (if_level_set) THEN
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un_in, un_m1_in, &
                        pn_in, pn_m1_in, incpn_in, incpn_m1_in, old_filename,  opt_level_set=level_setn_in, &
                        opt_level_set_m1=level_setn_m1_in, opt_max_vel=max_vel, opt_mono = mono_in)
                ELSE
                   CALL read_restart_ns(comm_one_d_ns, time_u, list_mode, un_in, un_m1_in, &
                        pn_in, pn_m1_in, incpn_in, incpn_m1_in, old_filename, opt_mono = mono_in)
                END IF

                CALL inter_mesh_loc_to_glob(vv_mesh_in, vv_mesh_out, un_in, un_out, l_t_g_vv, is_in, &
                     comm_one_d_ns(1))
                CALL inter_mesh_loc_to_glob(vv_mesh_in, vv_mesh_out, un_m1_in, un_m1_out, l_t_g_vv, is_in, &
                     comm_one_d_ns(1))
                CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, pn_in, pn_out, l_t_g_pp, is_in, &
                     comm_one_d_ns(1))
                CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, pn_m1_in, pn_m1_out, l_t_g_pp, is_in, &
                     comm_one_d_ns(1))
                CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, incpn_in, incpn_out, l_t_g_pp, is_in, &
                     comm_one_d_ns(1))
                CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, incpn_m1_in, incpn_m1_out, l_t_g_pp, is_in, &
                     comm_one_d_ns(1))
                IF (if_level_set) THEN
                   IF (if_level_set_P2) THEN
                      DO k = 1, nb_fluid-1
                         CALL inter_mesh_loc_to_glob(vv_mesh_in, vv_mesh_out, level_setn_in(k,:,:,:),&
                              level_setn_out(k,:,:,:), l_t_g_vv, is_in, comm_one_d_ns(1))
                         CALL inter_mesh_loc_to_glob(vv_mesh_in, vv_mesh_out, level_setn_m1_in(k,:,:,:),&
                              level_setn_m1_out(k,:,:,:), l_t_g_vv, is_in, comm_one_d_ns(1))
                      END DO
                   ELSE
                      DO k = 1, nb_fluid-1
                         CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, level_setn_in(k,:,:,:),&
                              level_setn_out(k,:,:,:), l_t_g_pp, is_in, comm_one_d_ns(1))
                         CALL inter_mesh_loc_to_glob(pp_mesh_in, pp_mesh_out, level_setn_m1_in(k,:,:,:),&
                              level_setn_m1_out(k,:,:,:), l_t_g_pp, is_in, comm_one_d_ns(1))
                      END DO
                   END IF
                END IF

                IF (if_level_set) THEN
                   CALL write_restart_ns(comm_one_d_ns, vv_mesh_out, pp_mesh_out, time_u, list_mode, un_out, un_m1_out, &
                        pn_out, pn_m1_out, incpn_out, incpn_m1_out, new_filename, m, 1, opt_level_set=level_setn_out, &
                        opt_level_set_m1=level_setn_m1_out, opt_max_vel=max_vel, opt_mono=mono_out)
                ELSE
                   CALL write_restart_ns(comm_one_d_ns, vv_mesh_out, pp_mesh_out, time_u, list_mode, un_out, un_m1_out, &
                        pn_out, pn_m1_out, incpn_out, incpn_m1_out, new_filename, m, 1, opt_mono=mono_out)
                END IF
                CALL MPI_COMM_RANK(comm_one_d_ns(1), rang_ns_S, ierr)
             END DO
          END IF
       END IF

       IF (check_plt) THEN
          IF (pp_mesh%me /= 0) THEN
             CALL MPI_COMM_RANK(comm_one_d_ns(1), rang_ns_S, ierr)
             WRITE(tit_s,'(i3)') rang_ns_S
             lblank = eval_blank(3,tit_s)
             DO l = 1, lblank - 1
                tit_s(l:l) = '0'
             END DO
             DO i = 1, SIZE(list_mode)
                WRITE(tit_m,'(i3)') list_mode(i)
                lblank = eval_blank(3,tit_m)
                DO l = 1, lblank - 1
                   tit_m(l:l) = '0'
                END DO
                CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,1,i), 'u_r_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,3,i), 'u_t_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,4,i), 'u_t_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(vv_mesh%jj, vv_mesh%rr, un(:,5,i), 'u_z_cos_m='//tit_m//'_'//tit_s//'_999.plt' )

             END DO

             IF (rang_ns_S == 0) THEN
                DO i = 1, SIZE(list_mode)
                   WRITE(tit_m,'(i3)') list_mode(i)
                   lblank = eval_blank(3,tit_m)
                   DO l = 1, lblank - 1
                      tit_m(l:l) = '0'
                   END DO
                   CALL plot_scalar_field(vv_mesh_glob%jj, vv_mesh_glob%rr, un_glob(:,1,i), 'gu_r_cos_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(vv_mesh_glob%jj, vv_mesh_glob%rr, un_glob(:,3,i), 'gu_t_cos_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(vv_mesh_glob%jj, vv_mesh_glob%rr, un_glob(:,4,i), 'gu_t_sin_m='//tit_m//'_999.plt' )
                   CALL plot_scalar_field(vv_mesh_glob%jj, vv_mesh_glob%rr, un_glob(:,5,i), 'gu_z_cos_m='//tit_m//'_999.plt' )
                END DO
             END IF
          END IF
       END IF
       IF (rank==0) WRITE(*,*) 'End interpolation NS'
    END IF

    !===Interpolation for temperature
    IF (rw_temp) THEN
       IF (rank==0) WRITE(*,*) 'Start interpolation temperature'

       IF (inter_mesh) THEN
          ALLOCATE(controle_temp(temp_mesh_out%np))
          DO m = index_start, index_start+nb_fic-1
             tempn_m1_glob = 0.d0
             tempn_glob    = 0.d0
             tempn_m1      = 0.d0
             tempn         = 0.d0
             WRITE(tit, '(i3)') m
             lblank = eval_blank(3,tit)
             DO l = 1, lblank - 1
                tit(l:l) = '0'
             END DO

             IF (petsc_rank==0) THEN
                CALL system('mv suite_temp_I'//tit//'.'//old_filename//'suite_temp.'//old_filename)
             END IF
             CALL MPI_Barrier( MPI_Comm_WORLD, code)

             IF (temp_mesh%me/=0) THEN
                CALL read_restart_temp(comm_one_d_temp, time_T, list_mode, tempn_in, tempn_m1_in, old_filename, &
                     opt_mono = mono_in)
             END IF
             controle_temp = 0
             CALL interp_mesh(temp_mesh_in, temp_mesh_out, tempn_in, tempn_out, controle_temp, 2)
             CALL interp_mesh(temp_mesh_in, temp_mesh_out, tempn_m1_in, tempn_m1_out, controle_temp, 2)

             IF (temp_mesh%me /= 0) THEN
                CALL write_restart_temp(comm_one_d_temp, temp_mesh_out, time_T, &
                     list_mode, tempn_out, tempn_m1_out, new_filename, m, 1, opt_mono = mono_out)
                CALL MPI_COMM_RANK(comm_one_d_temp(1), rang_temp_S, ierr)
             END IF
          END DO
          DEALLOCATE(controle_temp)

       ELSE

          ALLOCATE(l_t_g_temp(temp_mesh%np))
          l_t_g_temp = 0
          CALL loc_to_glob(temp_mesh, temp_mesh_glob, l_t_g_temp)
          IF (temp_mesh%me /=0) THEN
             DO m = index_start, index_start+nb_fic-1
                tempn_m1_glob = 0.d0
                tempn_glob    = 0.d0
                tempn_m1      = 0.d0
                tempn         = 0.d0

                WRITE(tit, '(i3)') m
                lblank = eval_blank(3,tit)
                DO l = 1, lblank - 1
                   tit(l:l) = '0'
                END DO
                IF (is_in) THEN
                   CALL MPI_COMM_RANK(comm_one_d_temp(1), rang_temp_S, ierr)
                   WRITE(tit_s,'(i3)') rang_temp_S
                   lblank = eval_blank(3,tit_s)
                   DO l = 1, lblank - 1
                      tit_s(l:l) = '0'
                   END DO

                   CALL system('mv suite_temp_S'//tit_s//'_I'//tit//'.'//old_filename//'suite_temp_S'//tit_s//'.'//old_filename)
                ELSE
                   IF (petsc_rank==0) CALL system('mv suite_temp_I'//tit//'.'//old_filename//'suite_temp.'//old_filename)
                   CALL MPI_Barrier( MPI_Comm_WORLD, code)
                END IF

                CALL read_restart_temp(comm_one_d_temp, time_T, list_mode, tempn_in, tempn_m1_in, old_filename, &
                     opt_mono = mono_in)

                CALL inter_mesh_loc_to_glob(temp_mesh_in, temp_mesh_out, tempn_in, tempn_out, l_t_g_temp, is_in, &
                     comm_one_d_temp(1))
                CALL inter_mesh_loc_to_glob(temp_mesh_in, temp_mesh_out, tempn_m1_in, tempn_m1_out, l_t_g_temp, is_in, &
                     comm_one_d_temp(1))

                CALL write_restart_temp(comm_one_d_temp, temp_mesh_out, time_T, list_mode, tempn_out, tempn_m1_out, &
                     new_filename, m, 1, opt_mono = mono_out)
                CALL MPI_COMM_RANK(comm_one_d_temp(1), rang_temp_S, ierr)
             END DO
          END IF
       END IF

       IF (check_plt) THEN
          IF (temp_mesh%me /= 0) THEN
             CALL MPI_COMM_RANK(comm_one_d_temp(1), rang_temp_S, ierr)
             WRITE(tit_s,'(i3)') rang_temp_S
             lblank = eval_blank(3,tit_s)
             DO l = 1, lblank - 1
                tit_s(l:l) = '0'
             END DO
             DO i = 1, SIZE(list_mode)
                WRITE(tit_m,'(i3)') list_mode(i)
                lblank = eval_blank(3,tit_m)
                DO l = 1, lblank - 1
                   tit_m(l:l) = '0'
                END DO
                CALL plot_scalar_field(temp_mesh%jj, temp_mesh%rr, tempn(:,1,i), 'temp_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(temp_mesh%jj, temp_mesh%rr, tempn(:,2,i), 'temp_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
             END DO
             IF (rang_temp_S == 0) THEN
                DO i = 1, SIZE(list_mode)
                   WRITE(tit_m,'(i3)') list_mode(i)
                   lblank = eval_blank(3,tit_m)
                   DO l = 1, lblank - 1
                      tit_m(l:l) = '0'
                   END DO
                   CALL plot_scalar_field(temp_mesh_glob%jj, temp_mesh_glob%rr, tempn(:,1,i), &
                        'gtemp_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                   CALL plot_scalar_field(temp_mesh_glob%jj, temp_mesh_glob%rr, tempn(:,2,i), &
                        'gtemp_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
                END DO
             END IF
          END IF
       END IF
       IF (rank==0) WRITE(*,*) 'End interpolation temperature'
    END IF

    !===Interpolation for concentration
    IF (rw_conc) THEN
       IF (rank==0) WRITE(*,*) 'Start interpolation concentration'

       IF (inter_mesh) THEN
          ALLOCATE(controle_conc(conc_mesh_out%np))
          DO m = index_start, index_start+nb_fic-1
             concn_m1_glob = 0.d0
             concn_glob    = 0.d0
             concn_m1      = 0.d0
             concn         = 0.d0
             WRITE(tit, '(i3)') m
             lblank = eval_blank(3,tit)
             DO l = 1, lblank - 1
                tit(l:l) = '0'
             END DO

             IF (petsc_rank==0) THEN
                CALL system('mv suite_conc_I'//tit//'.'//old_filename//'suite_conc.'//old_filename)
             END IF
             CALL MPI_Barrier( MPI_Comm_WORLD, code)

             IF (conc_mesh%me/=0) THEN
                CALL read_restart_conc(comm_one_d_conc, time_conc, list_mode, concn_in, concn_m1_in, old_filename, &
                     opt_mono = mono_in)
             END IF
             controle_conc = 0
             CALL interp_mesh(conc_mesh_in, conc_mesh_out, concn_in, concn_out, controle_conc, 2)
             CALL interp_mesh(conc_mesh_in, conc_mesh_out, concn_m1_in, concn_m1_out, controle_conc, 2)

             IF (conc_mesh%me /= 0) THEN
                CALL write_restart_conc(comm_one_d_conc, conc_mesh_out, time_conc, &
                     list_mode, concn_out, concn_m1_out, new_filename, m, 1, opt_mono = mono_out)
                CALL MPI_COMM_RANK(comm_one_d_conc(1), rang_conc_S, ierr)
             END IF
          END DO
          DEALLOCATE(controle_conc)

       ELSE

          ALLOCATE(l_t_g_conc(conc_mesh%np))
          l_t_g_conc = 0
          CALL loc_to_glob(conc_mesh, conc_mesh_glob, l_t_g_conc)
          IF (conc_mesh%me /=0) THEN
             DO m = index_start, index_start+nb_fic-1
                concn_m1_glob = 0.d0
                concn_glob    = 0.d0
                concn_m1      = 0.d0
                concn         = 0.d0

                WRITE(tit, '(i3)') m
                lblank = eval_blank(3,tit)
                DO l = 1, lblank - 1
                   tit(l:l) = '0'
                END DO
                IF (is_in) THEN
                   CALL MPI_COMM_RANK(comm_one_d_conc(1), rang_conc_S, ierr)
                   WRITE(tit_s,'(i3)') rang_conc_S
                   lblank = eval_blank(3,tit_s)
                   DO l = 1, lblank - 1
                      tit_s(l:l) = '0'
                   END DO

                   CALL system('mv suite_conc_S'//tit_s//'_I'//tit//'.'//old_filename//'suite_conc_S'//tit_s//'.'//old_filename)
                ELSE
                   IF (petsc_rank==0) CALL system('mv suite_conc_I'//tit//'.'//old_filename//'suite_conc.'//old_filename)
                   CALL MPI_Barrier( MPI_Comm_WORLD, code)
                END IF

                CALL read_restart_conc(comm_one_d_conc, time_conc, list_mode, concn_in, concn_m1_in, old_filename, &
                     opt_mono = mono_in)

                CALL inter_mesh_loc_to_glob(conc_mesh_in, conc_mesh_out, concn_in, concn_out, l_t_g_conc, is_in, &
                     comm_one_d_conc(1))
                CALL inter_mesh_loc_to_glob(conc_mesh_in, conc_mesh_out, concn_m1_in, concn_m1_out, l_t_g_conc, is_in, &
                     comm_one_d_conc(1))

                CALL write_restart_conc(comm_one_d_conc, conc_mesh_out, time_conc, list_mode, concn_out, concn_m1_out, &
                     new_filename, m, 1, opt_mono = mono_out)
                CALL MPI_COMM_RANK(comm_one_d_conc(1), rang_conc_S, ierr)
             END DO
          END IF
       END IF

       IF (check_plt) THEN
          IF (conc_mesh%me /= 0) THEN
             CALL MPI_COMM_RANK(comm_one_d_conc(1), rang_conc_S, ierr)
             WRITE(tit_s,'(i3)') rang_conc_S
             lblank = eval_blank(3,tit_s)
             DO l = 1, lblank - 1
                tit_s(l:l) = '0'
             END DO
             DO i = 1, SIZE(list_mode)
                WRITE(tit_m,'(i3)') list_mode(i)
                lblank = eval_blank(3,tit_m)
                DO l = 1, lblank - 1
                   tit_m(l:l) = '0'
                END DO
                CALL plot_scalar_field(conc_mesh%jj, conc_mesh%rr, concn(:,1,i), 'conc_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                CALL plot_scalar_field(conc_mesh%jj, conc_mesh%rr, concn(:,2,i), 'conc_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
             END DO
             IF (rang_conc_S == 0) THEN
                DO i = 1, SIZE(list_mode)
                   WRITE(tit_m,'(i3)') list_mode(i)
                   lblank = eval_blank(3,tit_m)
                   DO l = 1, lblank - 1
                      tit_m(l:l) = '0'
                   END DO
                   CALL plot_scalar_field(conc_mesh_glob%jj, conc_mesh_glob%rr, concn(:,1,i), &
                        'gconc_cos_m='//tit_m//'_'//tit_s//'_999.plt' )
                   CALL plot_scalar_field(conc_mesh_glob%jj, conc_mesh_glob%rr, concn(:,2,i), &
                        'gconc_sin_m='//tit_m//'_'//tit_s//'_999.plt' )
                END DO
             END IF
          END IF
       END IF
       IF (rank==0) WRITE(*,*) 'End interpolation concentration'
    END IF

    IF (is_in) THEN
       IF (petsc_rank==0 .AND. rw_mxw) CALL system('rm -rf suite_maxwell_S*')
       IF (petsc_rank==0 .AND. rw_ns) CALL system('rm -rf suite_ns_S*')
       IF (petsc_rank==0 .AND. rw_temp) CALL system('rm -rf suite_temp_S*')
       IF (petsc_rank==0 .AND. rw_conc) CALL system('rm -rf suite_conc_S*')
    ELSE
       IF (petsc_rank==0 .AND. rw_mxw) CALL system('rm -rf suite_maxwell.*')
       IF (petsc_rank==0 .AND. rw_ns) CALL system('rm -rf suite_ns.*')
       IF (petsc_rank==0 .AND. rw_temp) CALL system('rm -rf suite_temp.*')
       IF (petsc_rank==0 .AND. rw_conc) CALL system('rm -rf suite_conc.*')
    END IF

    IF (check_plt .AND. petsc_rank==0) THEN
       CALL system('mkdir -p PLT')
       CALL system('mv *.plt PLT/')
    END IF
    CALL MPI_Barrier( MPI_Comm_WORLD, code)

    CALL PetscFinalize(ierr)

  END SUBROUTINE mesh_interpol

  SUBROUTINE which_pb(is_in)
    USE my_util
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: is_in
    CHARACTER(len=200)   :: inline
    CALL getarg(1,inline)
    IF (TRIM(ADJUSTL(inline)) == 'petsc_nonpetsc') THEN
       is_in = .TRUE.
    ELSE IF (TRIM(ADJUSTL(inline)) == 'nonpetsc_petsc') THEN
       is_in = .FALSE.
    ELSE
       CALL error_Petsc('BUG in which_pb')
    END IF
  END SUBROUTINE which_pb

END MODULE mesh_interpolation
