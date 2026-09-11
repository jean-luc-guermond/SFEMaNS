PROGRAM prog_write_ns
   !----------------
   USE petscdmda, ONLY: PETSC_NULL_CHARACTER, PetscInitialize
   !----------------
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
  USE petscmpi
  USE restart
  USE def_type_field_mag
  USE def_type_field_vel
  USE user_data_module
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  TYPE(vel_field_type), POINTER                   :: vel_n, vel_n1
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
  !===Timing======================================================================
!  REAL(KIND=8)                                    :: tps, tploc, tploc_max=0.d0
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


  CALL write_restart_ns(comm_one_d_ns, vv_mesh, pp_mesh, time, list_mode, &
       vel_n%u, vel_n%u, vel_n%p, vel_n%p, vel_n%p, vel_n%p, inputs%file_name, 1, 1)

  CALL error_petsc("Finished writing NS restart file. Stop here for testing purposes.")
  
  CONTAINS 
  !===End of code=================================================================

END PROGRAM prog_write_ns
