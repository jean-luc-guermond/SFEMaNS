MODULE exponential_propagator
#ifdef USE_LIGHTKRYLOV
   USE LightKrylov, only: abstract_exptA_linop_rdp, abstract_vector_rdp
#endif
   USE def_type_mesh
   USE def_type_field
#include "petsc/finclude/petsc.h"
   USE petsc
   USE COMM_LK_SFEM
   USE def_type_field

 
  IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: this_module = 'exponential_propagator'
      PRIVATE 

       !===Navier-Stokes fields========================================================
       TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
       REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn
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
       !===Declare PETSC===============================================================
       PetscErrorCode :: ierr
       PetscMPIInt    :: rank
       MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d, comm_one_d_ns, comm_one_d_temp
       MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_conc

#ifdef USE_LIGHTKRYLOV
   TYPE, EXTENDS(abstract_exptA_linop_rdp), PUBLIC :: exptA_linop
#else
   TYPE,                                    PUBLIC :: exptA_linop
#endif
   CONTAINS
       PRIVATE
!       PROCEDURE, PASS(self),                   PUBLIC :: init => init_exptA
       PROCEDURE, PASS(self),                   PUBLIC :: matvec => exptA_matvec
       PROCEDURE, PASS(self),                   PUBLIC :: rmatvec => exptA_matvec
   END TYPE exptA_linop


   CONTAINS

!===========================
!===========================
!=== Procedures for exptA
!===========================
!===========================
   
   SUBROUTINE exptA_matvec(self, vec_in, vec_out)

       USE COMM_LK_SFEM
       USE def_type_field
       USE def_type_mesh
       USE initialization
       USE my_util
       USE input_data
       USE fourier_to_real_for_vtu
       USE post_processing_debug
       USE verbose
       USE sub_plot
       USE chaine_caractere
#include "petsc/finclude/petsc.h"
       USE petsc

     
       CLASS(exptA_linop),              INTENT(INOUT) :: self
#ifdef USE_LIGHTKRYLOV
       CLASS(abstract_vector_rdp),      INTENT(IN)    :: vec_in
       CLASS(abstract_vector_rdp),     INTENT(OUT)    :: vec_out
#else
       TYPE(mag_field_type),      INTENT(IN)    :: vec_in
       TYPE(mag_field_type),     INTENT(OUT)    :: vec_out
#endif
       !===Time iterations=============================================================
       REAL(KIND=8)                                    :: time
       INTEGER                                         :: it
       INTEGER, SAVE                                   :: it_tot=0
       !===Timing======================================================================
       REAL(KIND=8)                                    :: tps, tploc, tploc_max=0.d0
       
       !===Start PETSC and MPI (mandatory)=============================================
      !  CALL PetscInitialize(PETSC_NULL_CHARACTER,ierr)
       CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
     
      !===User reads his/her own data=================================================
     
       !===Initialize SFEMANS (mandatory)==============================================
       CALL build_pointers(vv_mesh, pp_mesh, H_mesh, phi_mesh, temp_mesh, conc_mesh,&
            interface_H_phi, interface_H_mu, list_mode, &
            un, pn, mag_field, vel, &
            vol_heat_capacity_field, temperature_diffusivity_field, &
            concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
            comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
            concentration, level_set, density, &
            der_un, visc_LES, visc_LES_level)
#ifdef USE_LIGHTKRYLOV
   SELECT TYPE(vec_in)
       TYPE is(mag_field_type)
           SELECT TYPE(vec_out)
            TYPE is(mag_field_type)
#endif
               CALL mag_field%zero()
               CALL LK_2_SFEM(mag_field, vec_in)
               time = mag_field%time
               !===============================================================================
               !                        TIME INTEGRATION                                      !
               !===============================================================================
               !===Start time loop
               tps = user_time()
               tploc_max = 0.d0

               DO it = 1, inputs%nb_iteration
                  it_tot = it_tot + 1
                  tploc =  user_time()
                  time = time + inputs%dt
                  CALL run_SFEMaNS(time, it)
                  
                  ! !===My postprocessing
                  ! IF (.NOT.inputs%test_de_convergence) THEN
                  !    CALL my_post_processing(it)
                  ! END IF
                  !===Write restart file
                  IF ((MOD(it_tot, inputs%freq_restart) == 0) .AND. (.NOT. inputs%LK%eigs)) THEN
                     CALL save_run(it_tot,inputs%freq_restart)
                  ENDIF
                  ! !===Write snapshot file
                  ! IF (MOD(it_tot, inputs%freq_snapshot) == 0) THEN
                  !    CALL save_snapshot(it_tot,inputs%freq_snapshot)
                  ! ENDIF
                  !===Timing
                  tploc = user_time() - tploc
                  IF (it>1) tploc_max = tploc_max + tploc
               ENDDO

               mag_field%time = time 
               CALL SFEM_2_LK(vec_out, mag_field)
             
               !===Timing======================================================================
               tps = user_time() - tps
               CALL write_verbose(rank,opt_tps=tps,opt_tploc_max=tploc_max)
             
               !===End of code=================================================================
#ifdef USE_LIGHTKRYLOV
            CLASS DEFAULT
               CALL type_error('vec_out','mag_field_type','OUT',this_module,'exptA_matvec')
            END SELECT
       CLASS DEFAULT
           CALL type_error('vec_in','mag_field_type','IN',this_module,'exptA_matvec')
       END SELECT
#endif
   END SUBROUTINE exptA_matvec

END MODULE exponential_propagator
