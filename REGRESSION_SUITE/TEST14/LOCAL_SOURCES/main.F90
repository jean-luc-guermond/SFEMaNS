PROGRAM mhd_prog
  USE def_type_mesh
  USE initialization
  USE my_util
  USE input_data
  !USE arpack_mhd
  USE fourier_to_real_for_vtu
  USE user_data
  USE post_processing_debug
  USE verbose
#include "petsc/finclude/petsc.h"
  USE petsc
  IMPLICIT NONE
  !===Navier-Stokes fields========================================================
  TYPE(mesh_type), POINTER                        :: pp_mesh, vv_mesh
  REAL(KIND=8), POINTER, DIMENSION(:,:,:)         :: un, pn
  TYPE(dyn_real_array_three), POINTER, DIMENSION(:):: der_un
  !===Maxwell fields==============================================================
  TYPE(mesh_type), POINTER                        :: H_mesh, phi_mesh
  TYPE(interface_type), POINTER                   :: interface_H_mu, interface_H_phi
  REAL(KIND=8), POINTER,      DIMENSION(:,:,:)    :: Hn, Bn, phin, vel
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
  INTEGER                                         :: it
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
       un, pn, Hn, Bn, phin, vel, &
       vol_heat_capacity_field, temperature_diffusivity_field, &
       concentration_diffusivity_field,mu_H_field, sigma_field, time, m_max_c, &
       comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc,temperature, &
       concentration, level_set, density, &
       der_un, visc_LES, visc_LES_level)

  !===============================================================================
  !                        VISUALIZATION WITHOUT COMPUTING                       !
  !===============================================================================
  IF (inputs%if_just_processing) THEN
     inputs%freq_plot=1
     CALL my_post_processing(1)
     CALL error_petsc('End post_processing')
  END IF

  !===============================================================================
  !                        EIGENVALUE PROBLEMS/ARPACK                            !
  !===============================================================================
  !IF (inputs%if_arpack) THEN
  !   !ATTENTION: m_max_c should be equal to 1, meaning each processors is dealing with 1 Fourier mode
  !   CALL solver_arpack_mhd(comm_one_d,H_mesh,phi_mesh,&
  !        inputs%dt,list_mode,mu_H_field)
  !   !===Postprocessing to check convergence
  !   IF (inputs%test_de_convergence) THEN
  !      CALL post_proc_test(vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
  !           un, pn, Hn, Bn, phin, temperature, level_set, mu_H_field, &
  !           time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp)
  !      CALL error_Petsc('End of convergence test')
  !      !IF (rank==0) WRITE(*,*) 'End of convergence test'
  !      !RETURN
  !   END IF
  !   !=== Put your postprocessing here
  !
  !   !===End of code for ARPACK problem
  !   CALL error_Petsc('END OF ARPACK, EXITING PRGM')
  !IF (rank==0) WRITE(*,*) 'END OF ARPACK, EXITING PRGM'
  !   !RETURN
  !END IF
  !===Test 14 disabled so user don't need to install arpack
  IF (inputs%if_regression) THEN
     !CALL regression(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
     !     un, pn, Hn, Bn, phin, temperature, level_set, concentration, mu_H_field, &
     !     time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)
     it= 123
     CALL PetscFinalize(ierr)
     CALL EXIT(it)
     CALL error_Petsc('End of convergence test')
  END IF

  !===============================================================================
  !                        TIME INTEGRATION                                      !
  !===============================================================================
  !===Start time loop
  tps = user_time()
  DO it = 1, inputs%nb_iteration
     tploc =  user_time()
     time = time + inputs%dt

     CALL run_SFEMaNS(time, it)

     !===My postprocessing
     IF (.NOT.inputs%test_de_convergence) THEN
        CALL my_post_processing(it)
     END IF

     !===Write restart file
     IF (MOD(it, inputs%freq_restart) == 0) THEN
        CALL  save_run(it,inputs%freq_restart)
     ENDIF

     !===Timing
     tploc = user_time() - tploc
     IF (it>1) tploc_max = tploc_max + tploc
  ENDDO

  !===Timing======================================================================
  tps = user_time() - tps
  CALL write_verbose(rank,opt_tps=tps,opt_tploc_max=tploc_max)

  !===Postprocessing to check convergence=========================================
  !IF (inputs%test_de_convergence) THEN
  IF (inputs%if_regression) THEN
     CALL regression(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
          un, pn, Hn, Bn, phin, temperature, level_set, concentration, mu_H_field, &
          time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)
     CALL error_Petsc('End of convergence test')
  END IF

  !===End of code=================================================================
  CALL error_Petsc('End of SFEMaNS')
CONTAINS

  SUBROUTINE my_post_processing(it)
    USE sub_plot
    USE chaine_caractere
    USE tn_axi
    USE boundary
    USE sft_parallele
    USE verbose
    USE vtk_viz
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
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
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
          err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
          IF (rank == 0) THEN
             !===L2 norm of magnetic field
             WRITE(41,*) time, err
          END IF
          err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
          norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
          IF (rank == 0) THEN
             !===L2 norm of div(Bn)
             WRITE(51,*) time, err, err/norm
             WRITE(52,*) time, err, norm
             WRITE(*,*) 'norm L2 of magnetic field', time, norm
          END IF
          DO i=1,SIZE(list_mode)
             norm = norm_S(comm_one_d, 'L2', H_mesh, list_mode(i:i), Hn(:,:,i:i))
             IF (rank_S == 0) THEN
                !===L2 norm of Fourier mode list_mode(i) of magnetic field Hn
                WRITE(200+list_mode(i),*) time, norm
             END IF
          END DO
       END IF ! end /=nst
    END IF ! end freq_en

    IF (MOD(it,inputs%freq_plot) == 0) THEN
       !===Plot whatever you want here
       IF (it==inputs%freq_plot) THEN
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
    USE sft_parallele
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
    USE sft_parallele
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

END PROGRAM mhd_prog
