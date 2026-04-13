PROGRAM prog_write_ns
  USE def_type_mesh
  USE initialization
  USE my_util
  USE input_data
  !USE arpack_mhd
  USE fourier_to_real_for_vtu
  USE user_data
  USE post_processing_debug
  USE verbose
  USE sub_plot
  USE chaine_caractere
#include "petsc/finclude/petsc.h"
  USE petsc
  USE restart
  USE def_type_field 
  USE user_data_module
  USE symmetric_field
  USE tn_axi
  
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE         :: ek, ek_sym, ek_anti
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE         :: em, em_sym, em_anti
  REAL(KIND=8)                                    :: norm
  REAL(KIND=8) :: ur_th_rpi_s, ut_th_rpi_s, uz_th_rpi_s, ur_th_rpi_a, ut_th_rpi_a, uz_th_rpi_a
  REAL(KIND=8) :: hr_th_rpi_s, ht_th_rpi_s, hz_th_rpi_s, hr_th_rpi_a, ht_th_rpi_a, hz_th_rpi_a
  REAL(KIND=8), PARAMETER :: PI = ACOS(-1.d0), tor_over_pol = 0.73d0,  alpha = 0.5d0
!   REAL(KIND=8), PARAMETER :: heigth_conc = 1.6d0, heigth_vel = 1.6d0, heigth_temp = 2.0d0, heigth_mag = 2.4d0
  LOGICAL :: test_loc, test_glob=.TRUE., test_glob_comm
  REAL(KIND=8), PARAMETER :: heigth_conc=1.2d0, heigth_vel=1.2d0, heigth_temp=1.6d0, heigth_mag=2.0d0, heigth_phi=2.4d0

  TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un
  !===Maxwell fields==============================================================
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: vel
  TYPE(mag_field_type), POINTER                   :: mag_field
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: sigma_field, mu_H_field
  !===Temperature field===========================================================
  TYPE(mesh_type), POINTER                        :: temp_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: temperature
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: vol_heat_capacity_field
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: temperature_diffusivity_field
  !===Concentration field===========================================================
  TYPE(mesh_type), POINTER                        :: conc_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: concentration
  REAL(KIND=8), POINTER,      DIMENSION(:)        :: concentration_diffusivity_field
  !===Level_set===================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: level_set
  !===Density=====================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: density
  !===LES=========================================================================
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: visc_LES
  REAL(KIND=8), POINTER, DIMENSION(:,:,:,:)       :: visc_LES_level
  !===Fourier modes===============================================================
  INTEGER                                         :: m_max_c
  INTEGER,      POINTER,      DIMENSION(:)        :: list_mode
  !===Time iterations=============================================================
  REAL(KIND=8)                                    :: time
  INTEGER                                         :: i
  !===Timing======================================================================
  REAL(KIND=8)                                    :: tps, tploc, tploc_max=0.d0
  !===Declare PETSC===============================================================
  PetscErrorCode :: ierr
  PetscMPIInt    :: rank
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d, comm_one_d_ns, comm_one_d_temp
  MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_conc

  !===Start PETSC and MPI (mandatory)=============================================
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  !===User reads his/her own data=================================================
  CALL read_user_data('data')

  !===Initialize SFEMANS (mandatory)==============================================
  CALL initial(vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, conc_mesh,&
       interface_H_phi, interface_H_mu, list_mode, &
       un, pn, mag_field, vel, &
       vol_heat_capacity_field, temperature_diffusivity_field, &
       concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
       comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
       concentration, level_set, density, &
       der_un, visc_LES, visc_LES_level)

!=== TESTING VELOCITY FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(heigth_vel, comm_one_d_ns(1), vv_mesh, un, 'u', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in velocity field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in velocity field modes ", list_mode
   END IF

!=== TESTING MAGNETIC FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(heigth_mag, comm_one_d(1), H_mesh, mag_field%Hn, 'h', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in magnetic field modes ", list_mode
   ELSE 
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in magnetic field modes ", list_mode
   END IF

!=== TESTING TEMPERATURE FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(heigth_temp, comm_one_d_temp(1), temp_mesh, temperature, 'T', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in temperature field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in temperature field modes ", list_mode
   END IF

!=== TESTING CONCENTRATION FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(heigth_conc, comm_one_d_conc(1), conc_mesh, concentration, 'c', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in concentration field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in concentration field modes ", list_mode
   END IF

!=== TESTING CONCENTRATION FIELD SYMMETRIES 

   CALL test_symmetry_phi(heigth_phi, heigth_mag, comm_one_d(1), phi_mesh, mag_field%phin, test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in phi field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in phi field modes ", list_mode
   END IF

!=== SUMMARY OF TESTS
   CALL MPI_Allreduce(test_glob, test_glob_comm, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
   IF (test_glob_comm) THEN
       WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS"
       WRITE(*,*) "1234567891"
   ELSE
       WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE"
   END IF

   CALL error_petsc("END OF ctest RPI/CENTRO SYM")
  
  CONTAINS 
  !===End of code=================================================================

  SUBROUTINE test_symmetry_u_H_temp_conc(heigth, communicator, mesh, field, if_u_h_T_c, test_loc)
       IMPLICIT NONE
       REAL(KIND=8),                         INTENT(IN)   :: heigth
       TYPE(mesh_type),                      INTENT(IN)   :: mesh
       REAL(KIND=8), DIMENSION(:,:,:),       INTENT(IN)   :: field
       CHARACTER(LEN=1),                     INTENT(IN)   :: if_u_h_T_c
       LOGICAL,                              INTENT(OUT)  :: test_loc
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ek, ek_sym, ek_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ref_rpi_sym, ref_rpi_anti, ref_centro_sym, ref_centro_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: err_rpi_sym, err_rpi_anti, err_centro_sym, err_centro_anti
       INTEGER                                            :: i
       REAL(KIND=8) :: ur_th_rpi_s, ut_th_rpi_s, uz_th_rpi_s, ur_th_rpi_a, ut_th_rpi_a, uz_th_rpi_a
       REAL(KIND=8) :: relat_tol = 1d-6
       MPI_Comm     :: communicator

!========== ANALYTICAL SOLUTIONS
      ALLOCATE(ref_rpi_sym(SIZE(list_mode)), ref_rpi_anti(SIZE(list_mode)), SOURCE=0.d0)
      ALLOCATE(ref_centro_sym(SIZE(list_mode)), ref_centro_anti(SIZE(list_mode)), SOURCE=0.d0)

      ur_th_rpi_s = 2*PI * heigth/2.d0 * (PI/heigth)**2 * 1/60
      IF (SIZE(field, 2) == 6) THEN
         ut_th_rpi_s = 2*PI * heigth/2.d0 * (4*tor_over_pol)**2 * 1/60
         uz_th_rpi_s = 2*PI * heigth/2.d0 * 1/4
      ELSE
         ut_th_rpi_s = 0.d0
         uz_th_rpi_s = 0.d0
      END IF

      ur_th_rpi_a = (alpha/2)**2 * ur_th_rpi_s
      IF (SIZE(field, 2) == 6) THEN
         ut_th_rpi_a = 2.d0/heigth*heigth * alpha**2 * ut_th_rpi_s
         uz_th_rpi_a = alpha**2 * uz_th_rpi_s
      ELSE
         ut_th_rpi_a = 0.d0
         uz_th_rpi_a = 0.d0
      END IF

   DO i=1, SIZE(list_mode)
   !=== mF = 0
      IF (list_mode(i) == 0) THEN
         ref_rpi_sym(i) = (ur_th_rpi_s + ut_th_rpi_s + uz_th_rpi_s)/2
         ref_rpi_anti(i) = (ur_th_rpi_a + ut_th_rpi_a + uz_th_rpi_a)/2

         ref_centro_sym(i) = (ur_th_rpi_s + ut_th_rpi_a + uz_th_rpi_s)/2
         ref_centro_anti(i) = (ur_th_rpi_a + ut_th_rpi_s + uz_th_rpi_a)/2
   !=== mF = 1
      ELSE IF (list_mode(i) == 1) THEN
         ref_rpi_sym(i) = (ur_th_rpi_s/2 + ut_th_rpi_s/2 + uz_th_rpi_s/2)/2
         ref_rpi_anti(i) = (ur_th_rpi_a/2 + ut_th_rpi_a/2 + uz_th_rpi_a/2)/2

         ref_centro_sym(i) = (ur_th_rpi_a/2 + ut_th_rpi_s/2 + uz_th_rpi_a/2)/2
         ref_centro_anti(i) = (ur_th_rpi_s/2 + ut_th_rpi_a/2 + uz_th_rpi_s/2)/2
      END IF
   END DO
!========== ANALYTICAL SOLUTIONS

      ALLOCATE(ek(SIZE(list_mode)), ek_sym(SIZE(list_mode)), ek_anti(SIZE(list_mode)))
      ALLOCATE(err_rpi_sym(SIZE(list_mode)), err_rpi_anti(SIZE(list_mode)))
      ALLOCATE(err_centro_sym(SIZE(list_mode)), err_centro_anti(SIZE(list_mode)))
      test_loc = .TRUE.

!========== TEST OF RPI SYMMETRY
      CALL val_ener_sym_rpi(communicator, mesh, list_mode, field, ek, ek_sym, ek_anti, if_u_h_T_c)

      err_rpi_sym = ABS(ek_sym - ref_rpi_sym)/ABS(ek_sym)
      err_rpi_anti = ABS(ek_anti - ref_rpi_anti)/ABS(ek_anti)
      DO i=1,SIZE(list_mode)
         IF ((list_mode(i)<2)) THEN
            IF ((MAXVAL(err_rpi_sym) > relat_tol) .OR. MAXVAL(err_rpi_anti) > relat_tol) THEN
               WRITE(*,*) 'FAILURE RPI symmetry test for '//trim(adjustl(if_u_h_T_c))//' mF=', list_mode(i)
               test_loc = .FALSE.
               WRITE(*,*) if_u_h_T_c,' ==> Rpi for mF=', list_mode(i), ': ', ek(i)
               
               WRITE(*,*) if_u_h_T_c,' ==> Rpi Sym for mF=', list_mode(i), ': ', ek_sym(i),&
               "ref=",ref_rpi_sym(i), 'err=', err_rpi_sym(i)
               
               WRITE(*,*) if_u_h_T_c,' ==> Rpi Anti-sym for mF=', list_mode(i), ': ', ek_anti(i),&
               "ref=",ref_rpi_anti(i), 'err=', err_rpi_anti(i)
            END IF
         END IF
      END DO

!========== TEST OF CENTRO SYMMETRY
      CALL val_ener_sym_centrale(communicator, mesh, list_mode, field, ek, ek_sym, ek_anti, if_u_h_T_c)

      err_centro_sym = ABS(ek_sym - ref_centro_sym)/ABS(ek_sym)
      err_centro_anti = ABS(ek_anti - ref_centro_anti)/ABS(ek_anti)
      DO i=1,SIZE(list_mode)
         IF ((list_mode(i)<2)) THEN
            IF ((MAXVAL(err_centro_sym) > relat_tol) .OR. MAXVAL(err_centro_anti) > relat_tol) THEN
               WRITE(*,*) 'FAILURE CENTRO symmetry test for '//trim(adjustl(if_u_h_T_c))//' mF=', list_mode(i)
               test_loc = .FALSE.
               
               WRITE(*,*) if_u_h_T_c,' ==> centro for mF=', list_mode(i), ': ', ek(i)
               
               WRITE(*,*) if_u_h_T_c,' ==> centro Sym for mF=', list_mode(i), ': ', ek_sym(i),&
               "ref=",ref_centro_sym(i), 'err=', err_centro_sym(i)
               
               WRITE(*,*) if_u_h_T_c,' ==> centro Anti-sym for mF=', list_mode(i), ': ', ek_anti(i),&
               "ref=",ref_centro_anti(i), 'err=', err_centro_anti(i)
            END IF
         END IF
      END DO
  END SUBROUTINE test_symmetry_u_H_temp_conc

  SUBROUTINE test_symmetry_phi(heigth_phi, heigth_mag, communicator, mesh, field, test_loc)
       IMPLICIT NONE
       REAL(KIND=8),                         INTENT(IN)   :: heigth_phi, heigth_mag
       TYPE(mesh_type),                      INTENT(IN)   :: mesh
       REAL(KIND=8), DIMENSION(:,:,:),       INTENT(IN)   :: field
       LOGICAL,                              INTENT(OUT)  :: test_loc
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ek, ek_sym, ek_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ref_rpi_sym, ref_rpi_anti, ref_centro_sym, ref_centro_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: err_rpi_sym, err_rpi_anti, err_centro_sym, err_centro_anti
       INTEGER                                            :: i
       REAL(KIND=8) :: ur_th_rpi_s, ut_th_rpi_s, uz_th_rpi_s, ur_th_rpi_a, ut_th_rpi_a, uz_th_rpi_a
       REAL(KIND=8) :: int_z, relat_tol = 1d-6
       MPI_Comm     :: communicator

!========== ANALYTICAL SOLUTIONS
      ALLOCATE(ref_rpi_sym(SIZE(list_mode)), ref_rpi_anti(SIZE(list_mode)), SOURCE=0.d0)
      ALLOCATE(ref_centro_sym(SIZE(list_mode)), ref_centro_anti(SIZE(list_mode)), SOURCE=0.d0)

      int_z = (heigth_phi-heigth_mag)/4.d0 - heigth_phi/(8*PI)*SIN(2*PI*heigth_mag/heigth_phi)
      ur_th_rpi_s = 2*PI * 2*(int_z) * (PI/heigth_phi)**2 * 1/60

      int_z = (heigth_phi-heigth_mag)/4.d0 + heigth_phi/(4*PI)*SIN(PI*heigth_mag/heigth_phi)
      ur_th_rpi_a = (alpha/2)**2 * (2*PI * 2*(int_z) * (PI/heigth_phi)**2 * 1/60)

      ut_th_rpi_s = 0.d0; ut_th_rpi_a = 0.d0; uz_th_rpi_s = 0.d0; uz_th_rpi_a = 0.d0

   DO i=1, SIZE(list_mode)
   !=== mF = 0
      IF (list_mode(i) == 0) THEN
         ref_rpi_sym(i) = (ur_th_rpi_s + ut_th_rpi_s + uz_th_rpi_s)/2
         ref_rpi_anti(i) = (ur_th_rpi_a + ut_th_rpi_a + uz_th_rpi_a)/2

         ref_centro_sym(i) = (ur_th_rpi_s + ut_th_rpi_a + uz_th_rpi_s)/2
         ref_centro_anti(i) = (ur_th_rpi_a + ut_th_rpi_s + uz_th_rpi_a)/2
   !=== mF = 1
      ELSE IF (list_mode(i) == 1) THEN
         ref_rpi_sym(i) = (ur_th_rpi_s/2 + ut_th_rpi_s/2 + uz_th_rpi_s/2)/2
         ref_rpi_anti(i) = (ur_th_rpi_a/2 + ut_th_rpi_a/2 + uz_th_rpi_a/2)/2

         ref_centro_sym(i) = (ur_th_rpi_a/2 + ut_th_rpi_s/2 + uz_th_rpi_a/2)/2
         ref_centro_anti(i) = (ur_th_rpi_s/2 + ut_th_rpi_a/2 + uz_th_rpi_s/2)/2
      END IF
   END DO
!========== ANALYTICAL SOLUTIONS

      ALLOCATE(ek(SIZE(list_mode)), ek_sym(SIZE(list_mode)), ek_anti(SIZE(list_mode)))
      ALLOCATE(err_rpi_sym(SIZE(list_mode)), err_rpi_anti(SIZE(list_mode)))
      ALLOCATE(err_centro_sym(SIZE(list_mode)), err_centro_anti(SIZE(list_mode)))
      test_loc = .TRUE.

!========== TEST OF RPI SYMMETRY
      CALL val_ener_sym_rpi(communicator, mesh, list_mode, field, ek, ek_sym, ek_anti, 'phi')

      err_rpi_sym = ABS(ek_sym - ref_rpi_sym)/ABS(ek_sym)
      err_rpi_anti = ABS(ek_anti - ref_rpi_anti)/ABS(ek_anti)
      DO i=1,SIZE(list_mode)
         IF ((list_mode(i)<2)) THEN
            IF ((MAXVAL(err_rpi_sym) > relat_tol) .OR. MAXVAL(err_rpi_anti) > relat_tol) THEN
               WRITE(*,*) 'FAILURE RPI symmetry test for '//trim(adjustl('phi'))//' mF=', list_mode(i)
               test_loc = .FALSE.
               WRITE(*,*) 'phi',' ==> Rpi for mF=', list_mode(i), ': ', ek(i)
               
               WRITE(*,*) 'phi',' ==> Rpi Sym for mF=', list_mode(i), ': ', ek_sym(i),&
               "ref=",ref_rpi_sym(i), 'err=', err_rpi_sym(i)
               
               WRITE(*,*) 'phi',' ==> Rpi Anti-sym for mF=', list_mode(i), ': ', ek_anti(i),&
               "ref=",ref_rpi_anti(i), 'err=', err_rpi_anti(i)
            END IF
         END IF
      END DO

!========== TEST OF CENTRO SYMMETRY
      CALL val_ener_sym_centrale(communicator, mesh, list_mode, field, ek, ek_sym, ek_anti, 'phi')

      err_centro_sym = ABS(ek_sym - ref_centro_sym)/ABS(ek_sym)
      err_centro_anti = ABS(ek_anti - ref_centro_anti)/ABS(ek_anti)
      DO i=1,SIZE(list_mode)
         IF ((list_mode(i)<2)) THEN
            IF ((MAXVAL(err_centro_sym) > relat_tol) .OR. MAXVAL(err_centro_anti) > relat_tol) THEN
               WRITE(*,*) 'FAILURE CENTRO symmetry test for '//trim(adjustl('phi'))//' mF=', list_mode(i)
               test_loc = .FALSE.
               
               WRITE(*,*) 'phi' ,' ==> centro for mF=', list_mode(i), ': ', ek(i)
               
               WRITE(*,*) 'phi',' ==> centro Sym for mF=', list_mode(i), ': ', ek_sym(i),&
               "ref=",ref_centro_sym(i), 'err=', err_centro_sym(i)
               
               WRITE(*,*) 'phi',' ==> centro Anti-sym for mF=', list_mode(i), ': ', ek_anti(i),&
               "ref=",ref_centro_anti(i), 'err=', err_centro_anti(i)
            END IF
         END IF
      END DO
  END SUBROUTINE test_symmetry_phi

  
  SUBROUTINE my_post_processing(it)
    USE sub_plot
    USE chaine_caractere
    USE tn_axi
    USE boundary
    USE fft_parallele
    USE verbose
    USE plot_vtk
    USE sfemans_tools
    USE subroutine_mass
    USE user_data
    IMPLICIT NONE
    INTEGER,                             INTENT(IN) :: it
    REAL(KIND=8)                                    :: err, norm
    INTEGER                                         :: i, it_plot
    CHARACTER(LEN=3)                                :: what
    INTEGER                                         :: rank_S, rank_F
    INTEGER                                         :: rank_ns_S, rank_ns_F
    REAL(KIND=8), DIMENSION(vv_mesh%np, 2, SIZE(list_mode)) :: level_1_P2
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: level_1_P1
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode)) :: chi, one_chi
    INTEGER :: nb_procs, m_max_pad, bloc_size
    LOGICAL, SAVE :: once_plot=.TRUE.
    !===VTU 2d======================================================================
    CHARACTER(LEN=3)   :: st_mode
    CHARACTER(LEN=200) :: header
    CHARACTER(LEN=3)   :: name_of_field

    !===Check ranks
    IF (vv_mesh%me /=0) THEN
       CALL MPI_Comm_rank(comm_one_d_ns(1), rank_ns_S, ierr)
       CALL MPI_Comm_rank(comm_one_d_ns(2), rank_ns_F, ierr)
    ELSE
       rank_ns_S = -1
       rank_ns_F = -1
    END IF
    CALL MPI_Comm_rank(comm_one_d(1), rank_S, ierr)
    CALL MPI_Comm_rank(comm_one_d(2), rank_F, ierr)

    !===Check divergence of fields
    IF (inputs%check_numerical_stability) THEN
       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
          norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un)
       ELSE
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, mag_field%Hn)
       END IF
       IF (norm>1.d8 .OR. isnan(norm)) THEN
          CALL error_petsc('From my_post_processing: numerical unstability')
       END IF
    END IF

    !===Put your postprocessing stuff here
    IF (MOD(it,inputs%freq_en) == 0) THEN

       !===Verbose
       CALL write_verbose(rank)

       IF (inputs%type_pb=='nst' .OR. inputs%type_pb=='mhd' .OR. inputs%type_pb=='fhd') THEN
          IF (inputs%if_compute_momentum_pseudo_force) THEN
             !===Compute the term -2/pi*integral((1-chi)*(u-u_solid).e_z/dt)
             !===chi is the penalty function, u_solid the velocity of the solid, dt the time step
             !==Output written in the file fort.12
             CALL FORCES_AND_MOMENTS(time,vv_mesh,comm_one_d_ns,list_mode,un)
          END IF

          err = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)
          norm = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             !===Divergence of velocity field
             WRITE(31,*) time, err/norm
          END IF
          DO i=1,SIZE(list_mode)
             norm = norm_S(comm_one_d, 'L2', vv_mesh, list_mode(i:i), un(:,:,i:i))
             IF (rank_ns_S == 0) THEN
                !===L2 norm of Fourier mode list_mode(i) of velocity field un
                WRITE(100+list_mode(i),*) time, norm
             END IF
          END DO

          err = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
          norm = norm_SF(comm_one_d, 'sH1', vv_mesh, list_mode, un)
          IF (rank == 0) THEN
             WRITE(98,*) time, err
             WRITE(*,*) 'norm L2 of velocity', time, err
             WRITE(*,*) 'semi norm H1 of velocity', time, norm
          END IF

          err = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn)
          IF (rank == 0) THEN
             WRITE(*,*) 'norm L2 of pressure', time, err
          END IF

          IF (inputs%if_level_set) THEN
             !===Compute the term integral(level_set-level_set_t=0)/integral(level_set_t=0)
             !===Output written in file fort.97
             IF (inputs%if_level_set_P2) THEN
                CALL compute_level_set_conservation(time,vv_mesh,comm_one_d_ns,list_mode,level_set)
             ELSE
                CALL compute_level_set_conservation(time,pp_mesh,comm_one_d_ns,list_mode,level_set)
             END IF
          END IF

          IF (inputs%if_temperature) THEN
             norm = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, temperature)
             IF (rank == 0) THEN
                WRITE(99,*) 'norm L2 of temperature', time, norm
                WRITE(*,*) 'norm L2 of temperature', time, norm
             END IF
          END IF
       END IF ! end nst or mhd or fhd

       IF (inputs%type_pb/='nst') THEN
          err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, mag_field%Hn)
          IF (rank == 0) THEN
             !===L2 norm of magnetic field
             WRITE(41,*) time, err
          END IF
          err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, mag_field%Bn)
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, mag_field%Bn)
          IF (rank == 0) THEN
             !===L2 norm of div(Bn)
             WRITE(51,*) time, err, err/norm
             WRITE(52,*) time, err, norm
             WRITE(*,*) 'norm L2 of magnetic field', time, norm
          END IF
          DO i=1,SIZE(list_mode)
             norm = norm_S(comm_one_d, 'L2', H_mesh, list_mode(i:i), mag_field%Hn(:,:,i:i))
             IF (rank_S == 0) THEN
                !===L2 norm of Fourier mode list_mode(i) of magnetic field Hn
                WRITE(200+list_mode(i),*) time, norm
             END IF
          END DO
       END IF ! end /=nst
    END IF ! end freq_en

    IF (MOD(it,inputs%freq_plot) == 0) THEN
       !===Plot whatever you want here
       !IF (it==inputs%freq_plot) THEN
       IF (once_plot) THEN
          once_plot=.FALSE.
          what = 'new'
       ELSE
          what = 'old'
       END IF
       it_plot = it/inputs%freq_plot
       !===Generation of 3D plots
       IF (inputs%if_level_set) THEN
          IF (inputs%if_level_set_P2) THEN
             level_1_P2=level_set(1,:,:,:)
             CALL vtu_3d(comm_one_d, level_1_P2, 'vv_mesh', 'Level_1', 'level_1', what, opt_it=it_plot)
          ELSE
             level_1_P1=level_set(1,:,:,:)
             CALL vtu_3d(comm_one_d, level_1_P1, 'pp_mesh', 'Level_1', 'level_1', what, opt_it=it_plot)
          END IF
          !CALL vtu_3d(density, 'vv_mesh', 'Density', 'density', what, opt_it=it_plot)
       END IF
       IF (inputs%type_pb/='mxw' .AND. inputs%type_pb/='mxx') THEN
          IF (inputs%type_fe_velocity==3) THEN
             CALL vtu_3d(comm_one_d, un, 'vv_mesh', 'Velocity', 'vel', what, opt_it=it_plot, opt_mesh_in=vv_mesh)
          ELSE
             CALL vtu_3d(comm_one_d, un, 'vv_mesh', 'Velocity', 'vel', what, opt_it=it_plot)
          ENDIF
          CALL vtu_3d(comm_one_d, pn, 'pp_mesh', 'Pressure', 'pre', what, opt_it=it_plot)
          CALL vtu_3d(comm_one_d, un, 'vv_mesh', 'Vorticity', 'vor', what, opt_it=it_plot, &
               opt_grad_curl='curl_u', opt_mesh_in=vv_mesh)
          !===Visualization of the penalty function
          IF (inputs%if_ns_penalty) THEN
             CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, ierr)
             m_max_pad = 3*SIZE(list_mode)*nb_procs/2
             bloc_size = SIZE(vv_mesh%rr,2)/nb_procs+1
             DO i=1,SIZE(list_mode)
                IF (list_mode(i)==0) THEN
                   one_chi(:,1,i)=1.d0
                   one_chi(:,2,i)=0.d0
                ELSE
                   one_chi(:,:,i)=0.d0
                END IF
             END DO
             CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(comm_one_d(2), penal_in_real_space, vv_mesh, &
                  one_chi, chi, nb_procs, bloc_size, m_max_pad, vv_mesh%rr, time)
             CALL vtu_3d(comm_one_d, chi, 'vv_mesh', 'Chi', 'chi', what, opt_it=it_plot)
          END IF
          !===Visualization of the penalty function
          !===2D plots for each mode of the velocity field and the pressure
          DO i = 1, m_max_c
             WRITE(st_mode,'(I3)') list_mode(i)
             header = 'Vn_'//'mode_'//trim(adjustl(st_mode))
             name_of_field = 'Vn'
             CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, un(:,:,i), name_of_field, what, opt_it=it_plot)
             WRITE(st_mode,'(I3)') list_mode(i)
             header = 'Pn_'//'mode_'//trim(adjustl(st_mode))
             name_of_field = 'Pn'
             CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, pn(:,:,i), name_of_field, what, opt_it=it_plot)
          END DO
       END IF
       IF (inputs%if_temperature) THEN
          !CALL vtu_3d(comm_one_d, temperature, 'temp_mesh', 'Temperature', 'temp', what, opt_it=it_plot)
       END IF
       IF (inputs%type_pb/='nst') THEN
          !CALL vtu_3d(comm_one_d, Hn, 'H_mesh', 'MagField', 'mag', what, opt_it=it_plot)
          !CALL vtu_3d(comm_one_d, Hn, 'H_mesh', 'Current', 'cur', what, opt_it=it_plot, &
          !     opt_grad_curl='curl_h', opt_2D=.FALSE.)
          !IF (inputs%nb_dom_phi>0) THEN
          !   CALL vtu_3d(comm_one_d, phin, 'phi_mesh', 'ScalPot', 'phi', what, opt_it=it_plot)
          !END IF
       END IF
       !==End generation of 3D plots

       !===Generation of 2D plots for each Fourier mode (uncomment if wanted)
       !===Proceed as follows to make 2D plots in the Fourier space (using Hn for instance)
       !===what = 'new' if first image, what= 'old' for later images  (CHARACTER(LEN=3)   :: what)
       !===WRITE(st_mode,'(I3)') list_mode(i)                         (CHARACTER(LEN=3)   :: st_mode)
       !===header = 'Hn_'//'_mode_'//trim(adjustl(st_mode))           (CHARACTER(LEN=200) :: header)
       !===name_of_field = 'Hn' (for instance)                        (CHARACTER(LEN=3)   :: name_of_field)
       !===CALL make_vtu_file_2D(comm_one_(1), H_mesh, header, Hn(:,:,i), name_of_field, what, opt_it=1)
!!$       IF (inputs%if_level_set) THEN
!!$          !===2D plots for each mode of the first level set
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Ln_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Ln'
!!$             IF (inputs%if_level_set_P2) THEN
!!$                level_1_P2(:,:,i)=level_set(1,:,:,i)
!!$                CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, level_1_P2(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             ELSE
!!$                level_1_P1(:,:,i)=level_set(1,:,:,i)
!!$                CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, level_1_P1(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             END IF
!!$             header = 'Dn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Dn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, density(:,:,i), name_of_field, what, opt_it=it_plot)
!!$          END DO
!!$       END IF
!!$       IF (inputs%type_pb/='mxw' .AND. inputs%type_pb/='mxx') THEN
!!$          !===2D plots for each mode of the velocity field and the pressure
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Vn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Vn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), vv_mesh, header, un(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Pn_'//'mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Pn'
!!$             CALL make_vtu_file_2D(comm_one_d_ns(1), pp_mesh, header, pn(:,:,i), name_of_field, what, opt_it=it_plot)
!!$          END DO
!!$       END IF
!!$       IF (inputs%if_temperature) THEN
!!$          !===2D plots for each mode of the temperature
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Tn_'//'_mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Tn'
!!$             CALL make_vtu_file_2D(comm_one_d_temp(1), temp_mesh, header, temperature(:,:,i), name_of_field, what, opt_it=it_plot) ! MODIFICATION: comm_one_d_temp instead of comm_one_d
!!$          END DO
!!$       END IF
!!$       IF (inputs%type_pb/='nst') THEN
!!$          !===2D plots for each mode of the magnetic field and the scalar potential
!!$          DO i = 1, m_max_c
!!$             WRITE(st_mode,'(I3)') list_mode(i)
!!$             header = 'Hn_'//'_mode_'//trim(adjustl(st_mode))
!!$             name_of_field = 'Hn'
!!$             CALL make_vtu_file_2D(comm_one_d(1), H_mesh, header, Hn(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             IF (inputs%nb_dom_phi>0) THEN
!!$                WRITE(st_mode,'(I3)') list_mode(i)
!!$                header = 'Phin_'//'_mode_'//trim(adjustl(st_mode))
!!$                name_of_field = 'Phin'
!!$                CALL make_vtu_file_2D(comm_one_d(1), phi_mesh, header, phin(:,:,i), name_of_field, what, opt_it=it_plot)
!!$             END IF
!!$          END DO
!!$       END IF
       !===End Generation of 2D plots for each Fourier mode (uncomment if wanted)

    END IF ! end freq_plot

  END SUBROUTINE my_post_processing

  SUBROUTINE FORCES_AND_MOMENTS(time,vv_mesh,communicator,list_mode,un)
    USE def_type_mesh
    USE input_data
    USE boundary
    USE fft_parallele
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                                        INTENT(IN)          :: vv_mesh
    INTEGER,      DIMENSION(:),                             INTENT(IN)          :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),                         INTENT(IN)          :: un
    REAL(KIND=8),                                           INTENT(IN)          :: time
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: vel_gauss, vel_gauss_penal
    REAL(KIND=8), DIMENSION(2,vv_mesh%gauss%l_G*vv_mesh%dom_me)                 :: rr_gauss
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                                  :: j_loc
    REAL(KIND=8)                                                                :: vel_torque, vel_torque_tot
    INTEGER ::  m, l , i, mode, index, type, nb_procs, m_max_pad, bloc_size
    PetscErrorCode                   :: ierr
    MPI_Comm,DIMENSION(2)            :: communicator

    index = 0
    DO m = 1, vv_mesh%dom_me
       j_loc = vv_mesh%jj(:,m)
       DO l = 1, vv_mesh%gauss%l_G
          index = index + 1
          rr_gauss(1,index) = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
          rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))
       END DO
    END DO

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO type = 1, 6
                vel_gauss(index,type,i) = SUM(un(j_loc,type,i)*vv_mesh%gauss%ww(:,l))*(3/(2*inputs%dt))
             END DO
          END DO
       END DO
       IF(inputs%if_impose_vel_in_solids) THEN
          IF (mode==0) THEN
             vel_gauss(:,:,i)  =  vel_gauss(:,:,i) - imposed_velocity_by_penalty(rr_gauss(:,:),time)
          ENDIF
       END IF
    END DO

    CALL MPI_COMM_SIZE(communicator(2), nb_procs, ierr)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(vel_gauss,1)/nb_procs+1
    CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator(2), penal_in_real_space, vv_mesh, &
         vel_gauss, vel_gauss_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)

    vel_torque   = 0.d0
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       IF (mode/=0) THEN
          CYCLE
       ELSE
          index = 0
          DO m = 1, vv_mesh%dom_me
             j_loc = vv_mesh%jj(:,m)
             DO l = 1, vv_mesh%gauss%l_G
                index = index + 1
                !===Force is int_domain ((1-chi)*(u-u_solid)/dt )dx
                vel_torque =  vel_torque + (vel_gauss_penal(index,5,i) - vel_gauss(index,5,i)) &
                     *rr_gauss(1,index)*vv_mesh%gauss%rj(l,m)
             END DO
          END DO
          CALL MPI_ALLREDUCE(vel_torque, vel_torque_tot,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator(1), ierr)
          WRITE(*,*) ' FORCES_AND_MOMENTS ', time, 2*ACOS(-1.d0)*vel_torque_tot/(0.5d0*ACOS(-1.d0))
          WRITE(12,*) time, 2*ACOS(-1.d0)*vel_torque_tot/(0.5d0*ACOS(-1.d0))
       END IF
    END DO
  END SUBROUTINE FORCES_AND_MOMENTS

  SUBROUTINE compute_level_set_conservation(time, mesh, communicator, list_mode, level_set)
    USE def_type_mesh
    USE input_data
    USE boundary
    USE fft_parallele
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                                        INTENT(IN)          :: mesh
    INTEGER,      DIMENSION(:),                             INTENT(IN)          :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:),                       INTENT(IN)          :: level_set
    REAL(KIND=8),                                           INTENT(IN)          :: time
    LOGICAL,                                    SAVE        :: once_compute=.TRUE.
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:),    SAVE        :: volum_init
    REAL(KIND=8)                                            :: volum_init_loc, volum_init_F
    REAL(KIND=8)                                            :: inte_fft_loc, inte_fft_tot_F
    REAL(KIND=8), DIMENSION(inputs%nb_fluid-1)              :: inte_fft_tot
    REAL(KIND=8), DIMENSION(mesh%np, 2, SIZE(list_mode))    :: level_posi_fft
    REAL(KIND=8)                                            :: ray
    INTEGER,      DIMENSION(mesh%gauss%n_w)                 :: j_loc
    INTEGER                                                 :: m, l , i, nb_inter
    INTEGER                                                 :: my_petscworld_rank, code
    MPI_Comm,DIMENSION(2)            :: communicator

    CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)

103 FORMAT(1500(e22.9,2x))

    !===Computation of initial integral of level set
    IF (once_compute) THEN
       once_compute = .FALSE.

       ALLOCATE(volum_init(SIZE(level_set,1)))

       DO nb_inter=1, SIZE(level_set,1)
          volum_init_loc = 0.d0
          DO i = 1, SIZE(list_mode)
             IF (list_mode(i)==0) THEN
                DO m = 1, mesh%me
                   j_loc = mesh%jj(:,m)
                   DO l = 1, mesh%gauss%l_G
                      !===Compute radius of Gauss point
                      ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
                      volum_init_loc = volum_init_loc + SUM(level_set_exact(nb_inter,1,mesh%rr(:,j_loc),list_mode(i),0.d0)* &
                           mesh%gauss%ww(:,l))*ray*mesh%gauss%rj(l,m)
                   END DO
                END DO
             END IF
          END DO
          volum_init_loc = volum_init_loc*2*ACOS(-1.d0)
          CALL  MPI_ALLREDUCE(volum_init_loc, volum_init_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
               communicator(2), code)
          CALL MPI_ALLREDUCE(volum_init_F, volum_init(nb_inter), 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
               communicator(1), code)
       END DO
       IF (my_petscworld_rank==0) THEN
          WRITE(*,*) 'mass initial = ', time, volum_init
       END IF
    END IF !end once_compute

    !===Computation of level set conservation relative error
    DO nb_inter=1, SIZE(level_set,1)
       level_posi_fft = level_set(nb_inter,:,:,:)
       inte_fft_loc = 0.d0
       DO i = 1, SIZE(list_mode)
          IF (list_mode(i)==0) THEN
             DO m = 1, mesh%me
                j_loc = mesh%jj(:,m)
                DO l = 1, mesh%gauss%l_G
                   !===Compute radius of Gauss point
                   ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
                   inte_fft_loc = inte_fft_loc + SUM(level_posi_fft(j_loc,1,i)*mesh%gauss%ww(:,l))* &
                        ray*mesh%gauss%rj(l,m)
                END DO
             END DO
          END IF
       END DO
       inte_fft_loc = inte_fft_loc*2*ACOS(-1.d0)
       CALL  MPI_ALLREDUCE(inte_fft_loc, inte_fft_tot_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            communicator(2), code)
       CALL MPI_ALLREDUCE(inte_fft_tot_F, inte_fft_tot(nb_inter), 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            communicator(1), code)
    END DO
    IF (my_petscworld_rank==0) THEN
       WRITE(*,*) 'relative mass error of level set at t = ', &
            time, ABS(1-inte_fft_tot/(volum_init+1.d-14))
       WRITE(97,103)  time, ABS(1-inte_fft_tot/(volum_init+1.d-14))
    END IF
  END SUBROUTINE compute_level_set_conservation

END PROGRAM prog_write_ns
