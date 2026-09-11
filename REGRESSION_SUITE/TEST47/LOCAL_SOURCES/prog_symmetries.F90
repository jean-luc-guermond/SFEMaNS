PROGRAM prog_write_ns
   !----------------
   USE petscdmda, ONLY: PETSC_NULL_CHARACTER, PetscInitialize
   USE petscmpi
   !----------------
  USE def_type_mesh
  USE initialization
  USE my_util
  USE input_data
  USE fourier_to_real_for_vtu
  USE user_data
  USE post_processing_debug
  USE verbose
  USE sub_plot
  USE chaine_caractere
  USE restart
  USE def_type_field_mag
  USE def_type_field_vel
  USE user_data_module
  USE symmetric_field
  USE tn_axi
  
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  TYPE(vel_field_type), POINTER                   :: vel_n, vel_n1
  REAL(KIND=8), PARAMETER :: PI = ACOS(-1.d0), tor_over_pol = 0.73d0,  alpha = 0.5d0
!   REAL(KIND=8), PARAMETER :: height_conc = 1.6d0, height_vel = 1.6d0, height_temp = 2.0d0, height_mag = 2.4d0
  LOGICAL :: test_loc, test_glob=.TRUE., test_glob_comm
  REAL(KIND=8), PARAMETER :: height_conc=1.2d0, height_vel=1.2d0, height_temp=1.6d0, height_mag=2.0d0, height_phi=2.4d0

  TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un
  !===Maxwell fields==============================================================
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: vel
  TYPE(mag_field_type), POINTER                   :: mag_n, mag_n1
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
  !===Declare PETSC===============================================================
  INTEGER :: ierr
  INTEGER    :: rank
  INTEGER, DIMENSION(:), POINTER  :: comm_one_d, comm_one_d_ns, comm_one_d_temp
  INTEGER, DIMENSION(:), POINTER  :: comm_one_d_conc

  !===Start PETSC and MPI (mandatory)=============================================
  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)

  !===User reads his/her own data=================================================
  CALL read_user_data('data')

  !===Initialize SFEMANS (mandatory)==============================================
  CALL initial(vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, conc_mesh,&
       interface_H_phi, interface_H_mu, list_mode, &
       vel_n, vel_n1, mag_n, mag_n1, vel, &
       vol_heat_capacity_field, temperature_diffusivity_field, &
       concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
       comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
       concentration, level_set, density, &
       der_un, visc_LES, visc_LES_level)

!=== TESTING VELOCITY FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(height_vel, comm_one_d_ns(1), vv_mesh, vel_n%u, 'u', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in velocity field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in velocity field modes ", list_mode
   END IF

!=== TESTING MAGNETIC FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(height_mag, comm_one_d(1), H_mesh, mag_n%H, 'h', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in magnetic field modes ", list_mode
   ELSE 
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in magnetic field modes ", list_mode
   END IF

!=== TESTING TEMPERATURE FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(height_temp, comm_one_d_temp(1), temp_mesh, temperature, 'T', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in temperature field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in temperature field modes ", list_mode
   END IF

!=== TESTING CONCENTRATION FIELD SYMMETRIES 

   CALL test_symmetry_u_H_temp_conc(height_conc, comm_one_d_conc(1), conc_mesh, concentration, 'c', test_loc)
   test_glob = test_glob .AND. test_loc
   IF (.NOT. test_loc) THEN
      WRITE(*,*) "ctest RPI/CENTRO SYM: FAILURE in concentration field modes ", list_mode
   ELSE
      WRITE(*,*) "ctest RPI/CENTRO SYM: SUCCESS in concentration field modes ", list_mode
   END IF

!=== TESTING CONCENTRATION FIELD SYMMETRIES 

   CALL test_symmetry_phi(height_phi, height_mag, comm_one_d(1), phi_mesh, mag_n%phi, test_loc)
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

  SUBROUTINE test_symmetry_u_H_temp_conc(height, communicator, mesh, field, if_u_h_T_c, test_loc)
       IMPLICIT NONE
       REAL(KIND=8),                         INTENT(IN)   :: height
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
       INTEGER     :: communicator

!========== ANALYTICAL SOLUTIONS
      ALLOCATE(ref_rpi_sym(SIZE(list_mode)), ref_rpi_anti(SIZE(list_mode)), SOURCE=0.d0)
      ALLOCATE(ref_centro_sym(SIZE(list_mode)), ref_centro_anti(SIZE(list_mode)), SOURCE=0.d0)

      ur_th_rpi_s = 2*PI * height/2.d0 * (PI/height)**2 * 1/60
      IF (SIZE(field, 2) == 6) THEN
         ut_th_rpi_s = 2*PI * height/2.d0 * (4*tor_over_pol)**2 * 1/60
         uz_th_rpi_s = 2*PI * height/2.d0 * 1/4
      ELSE
         ut_th_rpi_s = 0.d0
         uz_th_rpi_s = 0.d0
      END IF

      ur_th_rpi_a = (alpha/2)**2 * ur_th_rpi_s
      IF (SIZE(field, 2) == 6) THEN
         ut_th_rpi_a = 2.d0/height*height * alpha**2 * ut_th_rpi_s
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

  SUBROUTINE test_symmetry_phi(height_phi, height_mag, communicator, mesh, field, test_loc)
       IMPLICIT NONE
       REAL(KIND=8),                         INTENT(IN)   :: height_phi, height_mag
       TYPE(mesh_type),                      INTENT(IN)   :: mesh
       REAL(KIND=8), DIMENSION(:,:,:),       INTENT(IN)   :: field
       LOGICAL,                              INTENT(OUT)  :: test_loc
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ek, ek_sym, ek_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: ref_rpi_sym, ref_rpi_anti, ref_centro_sym, ref_centro_anti
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE            :: err_rpi_sym, err_rpi_anti, err_centro_sym, err_centro_anti
       INTEGER                                            :: i
       REAL(KIND=8) :: ur_th_rpi_s, ut_th_rpi_s, uz_th_rpi_s, ur_th_rpi_a, ut_th_rpi_a, uz_th_rpi_a
       REAL(KIND=8) :: int_z, relat_tol = 1d-6
       INTEGER     :: communicator

!========== ANALYTICAL SOLUTIONS
      ALLOCATE(ref_rpi_sym(SIZE(list_mode)), ref_rpi_anti(SIZE(list_mode)), SOURCE=0.d0)
      ALLOCATE(ref_centro_sym(SIZE(list_mode)), ref_centro_anti(SIZE(list_mode)), SOURCE=0.d0)

      int_z = (height_phi-height_mag)/4.d0 - height_phi/(8*PI)*SIN(2*PI*height_mag/height_phi)
      ur_th_rpi_s = 2*PI * 2*(int_z) * (PI/height_phi)**2 * 1/60

      int_z = (height_phi-height_mag)/4.d0 + height_phi/(4*PI)*SIN(PI*height_mag/height_phi)
      ur_th_rpi_a = (alpha/2)**2 * (2*PI * 2*(int_z) * (PI/height_phi)**2 * 1/60)

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

END PROGRAM prog_write_ns
