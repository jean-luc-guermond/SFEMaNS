!Authors Jean-Luc Guermond, Hugues Faller, 2019
!cf https://doi.org/10.1137/18M1209301
MODULE update_taylor_navier_stokes
  USE my_util
  USE input_data
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC :: navier_stokes_taylor, init_velocity_pressure_taylor
  PRIVATE

CONTAINS

  SUBROUTINE navier_stokes_taylor(comm_one_d_ns,time, vv_3_LA, pp_1_LA, &
       list_mode, pp_mesh, vv_mesh, pn, der_pn, un, der_un, vvz_per, pp_per, density, tempn, concn)
    USE def_type_mesh
    USE periodic

    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: vv_mesh, pp_mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: un
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT)  :: pn
    TYPE(dyn_real_array_three), DIMENSION(:)       :: der_un
    TYPE(dyn_real_array_three), DIMENSION(:)       :: der_pn
    REAL(KIND=8),                   INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: density
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: tempn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: concn
    TYPE(periodic_type),            INTENT(IN)     :: vvz_per, pp_per
    TYPE(petsc_csr_LA)                             :: vv_3_LA, pp_1_LA
    LOGICAL                                        :: if_navier_stokes_with_taylor = .TRUE. !HF We will nead a 'read_data'
#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_ns

    IF (if_navier_stokes_with_taylor) THEN
       IF (inputs%taylor_order==3) THEN
          CALL update_ns_with_taylor(comm_one_d_ns,time,vv_3_LA, pp_1_LA, vvz_per, pp_per, &
               inputs%dt, inputs%Re, inputs%taylor_lambda, list_mode, pp_mesh, vv_mesh, pn, &
               der_pn, un, der_un, density, tempn, concn)
       ELSE IF (inputs%taylor_order==4) THEN
          CALL update_ns_with_taylor_fourth(comm_one_d_ns,time,vv_3_LA, pp_1_LA, vvz_per, pp_per, &
               inputs%dt, inputs%Re, inputs%taylor_lambda, list_mode, pp_mesh, vv_mesh, pn, der_pn,&
               un, der_un, density, tempn, concn)
       ELSE
          CALL error_petsc('BUG in navier_stokes_taylor: inputs%taylor_order not well defined')
       END IF
    END IF
  END SUBROUTINE navier_stokes_taylor

  SUBROUTINE init_velocity_pressure_taylor(vv_mesh, pp_mesh, list_mode, time, pn, der_pn, un, der_un)
    USE def_type_mesh
    !USE boundary_generic
    !JLG+MC 2022
    USE boundary
    !JLG+MC 2022
    IMPLICIT NONE
    INTEGER,      DIMENSION(:),                            INTENT(IN)    :: list_mode
    TYPE(mesh_type),                                       INTENT(IN)    :: vv_mesh, pp_mesh
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)), INTENT(INOUT) :: un
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT) :: pn
    TYPE(dyn_real_array_three), DIMENSION(inputs%taylor_order-1)         :: der_un
    TYPE(dyn_real_array_three), DIMENSION(inputs%taylor_order-1)         :: der_pn
    REAL(KIND=8), INTENT(IN) :: time
    REAL(KIND=8) :: dt
    INTEGER :: i, mode, k, kp, m_max_c

    dt = inputs%dt
    m_max_c = SIZE(list_mode)
    DO i= 1, m_max_c
       mode = list_mode(i)
       DO k = 1, 6
          !===initialize velocity at time=t0
          un(:,k,i) = vv_exact(k,vv_mesh%rr, mode,time)
          der_un(1)%DRT(:,k,i)= &
               (-2*vv_exact(k,vv_mesh%rr, mode,time-3*dt)&
               +9*vv_exact(k,vv_mesh%rr, mode,time-2*dt)&
               -18*vv_exact(k,vv_mesh%rr, mode,time-dt)&
               +11*vv_exact(k,vv_mesh%rr, mode,time))/(6*dt)
          der_un(2)%DRT(:,k,i)= &
               (-1*vv_exact(k,vv_mesh%rr, mode,time-3*dt)&
               +4*vv_exact(k,vv_mesh%rr, mode,time-2*dt)&
               -5*vv_exact(k,vv_mesh%rr, mode,time-dt)&
               +2*vv_exact(k,vv_mesh%rr, mode,time))/(dt**2)
       END DO
       DO k = 1, 2
          !===pressure at time=t0
          pn(:,k,i) = pp_exact(k,pp_mesh%rr, mode,time)
          der_pn(1)%DRT(:,k,i)= &
               (-2*pp_exact(k,pp_mesh%rr, mode,time-3*dt)&
               +9*pp_exact(k,pp_mesh%rr, mode,time-2*dt)&
               -18*pp_exact(k,pp_mesh%rr, mode,time-dt)&
               +11*pp_exact(k,pp_mesh%rr, mode,time))/(6*dt)
          der_pn(2)%DRT(:,k,i)= &
               (-1*pp_exact(k,pp_mesh%rr, mode,time-3*dt)&
               +4*pp_exact(k,pp_mesh%rr, mode,time-2*dt)&
               -5*pp_exact(k,pp_mesh%rr, mode,time-dt)&
               +2*pp_exact(k,pp_mesh%rr, mode,time))/(dt**2)
       ENDDO
       IF (inputs%taylor_order==4) THEN
          DO k = 1, 6
             der_un(3)%DRT(:,k,i)= &
                  (-1*vv_exact(k,vv_mesh%rr, mode,time-3*dt)&
                  +3*vv_exact(k,vv_mesh%rr, mode,time-2*dt)&
                  -3*vv_exact(k,vv_mesh%rr, mode,time-dt)&
                  +1*vv_exact(k,vv_mesh%rr, mode,time))/(dt**3)
          END DO
          DO k = 1, 2
             der_pn(3)%DRT(:,k,i)= &
                  (-1*pp_exact(k,pp_mesh%rr, mode,time-3*dt)&
                  +3*pp_exact(k,pp_mesh%rr, mode,time-2*dt)&
                  -3*pp_exact(k,pp_mesh%rr, mode,time-dt)&
                  +1*pp_exact(k,pp_mesh%rr, mode,time))/(dt**3)
          END DO
       END IF
    ENDDO

    !===Set to zero
    IF (vv_mesh%me/=0) THEN
       DO kp = 1, inputs%taylor_order-1
          DO i = 1, m_max_c
             IF (list_mode(i) == 0) THEN
                DO k= 1, 3!===velocity
                   der_un(kp)%DRT(:,2*k,i) = 0.d0
                END DO
                der_pn(kp)%DRT(:,2,i) = 0.d0!===pressure
             END IF
          END DO
       END DO
    END IF
  END SUBROUTINE init_velocity_pressure_taylor

  SUBROUTINE update_ns_with_taylor(comm_one_d, time, vv_3_LA, pp_1_LA, vvz_per, pp_per, &
       dt, Re, lambda, list_mode, pp_mesh, vv_mesh, pn, der_pn, un, der_un, density, tempn, concn)
    !==============================
    USE def_type_mesh
    USE fem_M_axi
    USE fem_rhs_axi
    USE Dir_nodes_petsc
    USE periodic
    USE st_matrix
    USE solve_petsc
    USE boundary
    USE st_matrix
    USE input_data
    USE rhs_para_assembling
    USE tn_axi
    USE verbose
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8)                                        :: time, dt, Re, lambda
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: list_mode
    REAL(KIND=8),      DIMENSION(:,:,:),     INTENT(IN) :: density
    REAL(KIND=8),      DIMENSION(:,:,:),     INTENT(IN) :: tempn
    REAL(KIND=8),      DIMENSION(:,:,:),     INTENT(IN) :: concn
    TYPE(mesh_type),                INTENT(IN)          :: pp_mesh, vv_mesh
    TYPE(petsc_csr_LA)                                  :: vv_3_LA, pp_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: vvz_per, pp_per
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: un
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT)  :: pn
    TYPE(dyn_real_array_three), DIMENSION(:)            :: der_un
    TYPE(dyn_real_array_three), DIMENSION(:)            :: der_pn
    !===Saved variables
    INTEGER,                                       SAVE :: m_max_c
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: pp_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: pp_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3),              SAVE :: vv_js_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: vv_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: vel_global_D
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: my_petscworld_rank
    !===End saved variables
    INTEGER,          POINTER, DIMENSION(:)  :: pp_1_ifrom, vv_3_ifrom
    INTEGER                                  :: i, m, n, n1, n2, n3, n123
    INTEGER                                  :: nb_procs, code, nu_mat, mode
    REAL(KIND=8)                             :: moyenne
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: div
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)                  :: un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rotv_v1, rotv_v2, rotv_v3
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))   :: uext, dtt_un_p1, dt_un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np)                     :: vel_loc, vel_tot
    REAL(KIND=8) :: tps, tps_tot, tps_cumul, coeff, vloc, cfl, cfl_max, norm
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/
    !===Communicators for Petsc, in space and Fourier
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: vel_mat
    Mat,                        SAVE :: mass_mat
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Vec,                        SAVE :: pb_1,    pb_2,    px_1, px_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: vel_ksp
    KSP,                        SAVE :: mass_ksp
    !===END OF DECLARATIONS

    IF (once) THEN
       once = .FALSE.

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)
       !===CREATE PETSC VECTORS AND GHOSTS
       CALL create_my_ghost(vv_mesh,vv_3_LA,vv_3_ifrom)
       n = 3*vv_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
       CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
       CALL VecDuplicate(vx_3, vb_3_145, ierr)
       CALL VecDuplicate(vx_3, vb_3_236, ierr)

       CALL create_my_ghost(pp_mesh,pp_1_LA,pp_1_ifrom)
       n = pp_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(pp_1_ifrom), pp_1_ifrom, px_1, ierr)
       CALL VecGhostGetLocalForm(px_1, px_1_ghost, ierr)
       CALL VecDuplicate(px_1, pb_1, ierr)
       CALL VecDuplicate(px_1, pb_2, ierr)
       !===End CREATE PETSC VECTORS AND GHOSTS

       !===Number of modes on proc
       m_max_c = SIZE(list_mode)
       !===End umber of modes on proc


       !===PREPARE pp_mode_global_js_D ARRAY FOR PRESSURE
       CALL scalar_without_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
       ALLOCATE(pp_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(pp_global_D(i)%DRL(SIZE(pp_mode_global_js_D(i)%DIL)))
       END DO
       !===End PREPARE pp_mode_global_js_D  ARRAY FOR PRESSURE

       !===PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY
       CALL vector_glob_js_D(vv_mesh, list_mode, vv_3_LA, inputs%vv_list_dirichlet_sides, &
            vv_js_D, vv_mode_global_js_D)

       ALLOCATE(vel_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(vel_global_D(i)%DRL(SIZE(vv_mode_global_js_D(i)%DIL)))
       END DO
       !===END PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY

       !===ASSEMBLE MASS MATRIX
       CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, mass_mat, CLEAN=.FALSE.)
       CALL qs_diff_mass_scal_M (pp_mesh, pp_1_LA, 0.d0, 1.d0, 0.d0, 0, mass_mat)
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_matrix_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, mass_mat, pp_1_LA)
       END IF
       DO i = 1, m_max_c
          CALL Dirichlet_M_parallel(mass_mat,pp_mode_global_js_D(i)%DIL)
       END DO
       CALL init_solver(inputs%my_par_mass,mass_ksp,mass_mat,comm_one_d(1),&
            solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)
       !===END ASSEMBLE MASS MATRIX

       !===ASSEMBLING VELOCITY MATRICES
       ALLOCATE(vel_mat(2*m_max_c),vel_ksp(2*m_max_c))
       DO i = 1, m_max_c
          mode = list_mode(i)
          !===VELOCITY
          nu_mat = 2*i-1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_taylor (1, vv_3_LA, vv_mesh, one/Re, one/dt, & !HF new coeffs
               lambda, i, mode, vel_mat(nu_mat))

          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          nu_mat = nu_mat+1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_taylor (2, vv_3_LA, vv_mesh, one/Re, one/dt, & !HF new coeffs
               lambda, i, mode, vel_mat(nu_mat))

          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          !===End VELOCITY
       ENDDO
       !===End ASSEMBLING VELOCITY MATRICES
    ENDIF !===End of once

    !===Compute NL for dtt_un by FFT at Gauss points
    tps_tot = user_time()
    tps_cumul = 0
    tps = user_time()

    uext = un + 2.d0*dt*der_un(1)%DRT + 2.d0*dt**2 *der_un(2)%DRT !HF 1st NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v1,inputs%precession)

    uext = un + dt*der_un(1)%DRT + 0.5d0*dt**2*der_un(2)%DRT !HF 2nd NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v2,inputs%precession)

    uext = un
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v3,inputs%precession)
    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !===End Compute NL by FFT at Gauss points

    !===Computation of CFL
    IF (inputs%verbose_CFL) THEN
       vel_loc = 0.d0
       DO i = 1, m_max_c
          IF (list_mode(i)==0) THEN
             coeff = 1.d0
          ELSE
             coeff = .5d0
          END IF
          vel_loc = vel_loc + coeff*(un(:,1,i)**2+un(:,2,i)**2+un(:,3,i)**2 &
               +un(:,4,i)**2+un(:,5,i)**2+un(:,6,i)**2)
       END DO
       CALL MPI_COMM_SIZE(comm_one_d(2),nb_procs,code)
       CALL MPI_ALLREDUCE(vel_loc,vel_tot,vv_mesh%np,MPI_DOUBLE_PRECISION, MPI_SUM, comm_one_d(2), code)
       vel_tot = sqrt(vel_tot)
       cfl = 0.d0
       DO m = 1, vv_mesh%dom_me
          vloc = MAXVAL(vel_tot(vv_mesh%jj(:,m)))
          cfl = MAX(vloc*dt/MIN(vv_mesh%hloc(m),MAXVAL(vv_mesh%hm)),cfl)
       END DO
       CALL MPI_ALLREDUCE(cfl,cfl_max,1,MPI_DOUBLE_PRECISION, MPI_MAX, comm_one_d(1), code)
       talk_to_me%CFL=cfl_max
       talk_to_me%time=time
    END IF
    !===End Computation of CFL

    !===Computation of second order derivatives
    !===Computation of rhs for dtt_un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         der_un(2)%DRT/dt, der_pn(2)%DRT, (rotv_v1-2.d0*rotv_v2+rotv_v3)/(dt**2), rhs_gauss,2, density, tempn, concn)
    !===End Computation of rhs for dtt_un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        =(vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time+dt)&
            -2.d0*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)   &
            +vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =(vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time+dt)&
            -2.d0*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)   &
            +vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n1+n2+1:n123)=(vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time+dt)&
            -2.d0*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)   &
            +vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)

       vel_global_D(i)%DRL(1:n1)        =(vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time+dt)&
            -2.d0*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)   &
            +vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =(vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time+dt)&
            -2.d0*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)   &
            +vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n1+n2+1:n123)=(vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time+dt)&
            -2.d0*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)   &
            +vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt))/(dt**2)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, dtt_ur_c, dtt_ut_s, dtt_uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, dtt_ur_s, dtt_ut_c, dtt_uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of dtt_velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for dtt_pressure
       pp_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for dtt_pressure

       !===Solve mass matrix for dtt_pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for dtt_pressure correction

       !===Compute 2nd order derrivative pressure for this mode
       der_pn(2)%DRT(:,:,i) =  der_pn(2)%DRT(:,:,i) - lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End 2nd order derrivative Pressure computation we will then take mean out

       !===UPDATES
       tps = user_time()
       !===Handling 2nd order derrivative of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, der_pn(2)%DRT(:,1,i),moyenne)
          der_pn(2)%DRT(:,1,i) = der_pn(2)%DRT(:,1,i)-moyenne
       ENDIF
       !===End of handling of 2nd order derrivative mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)    = 0.d0
          un_p1 (:,4)    = 0.d0
          un_p1 (:,6)    = 0.d0
          der_pn(2)%DRT (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATES 2nd order derrivative
       dtt_un_p1(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES 2nd order derrivative

    ENDDO
    !===END Computation of second order derivatives

    !===Computation of first order derivatives
    !===Computation of rhs for dt_un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         (der_un(1)%DRT-dt/2.d0*(dtt_un_p1-der_un(2)%DRT))/dt, der_pn(1)%DRT+ dt*der_pn(2)%DRT, &
         (rotv_v1-rotv_v3)/(2.d0*dt), rhs_gauss,1, density, tempn, concn)
    !===End Computation of rhs for dt_un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        =(vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time+dt)&
            -vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =(vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time+dt)&
            -vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n1+n2+1:n123)=(vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time+dt)&
            -vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        =(vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time+dt)&
            -vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =(vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time+dt)&
            -vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n1+n2+1:n123)=(vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time+dt)&
            -vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt))/(2.d0*dt)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, dt_ur_c, dt_ut_s, dt_uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, dt_ur_s, dt_ut_c, dt_uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of dt_velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for dt_pressure
       pp_global_D(i)%DRL = 0.d0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for dt_pressure

       !===Solve mass matrix for dt_pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for dt_pressure correction

       !===Compute first order derrivative pressure for this mode
       !===Pressure first order derrivative computation
       der_pn(1)%DRT(:,:,i) = der_pn(1)%DRT(:,:,i) +dt*der_pn(2)%DRT(:,:,i)- lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End dt_Pressure computation

       !===UPDATES
       tps = user_time()
       !===Handling of mean first order derrivative pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, der_pn(1)%DRT(:,1,i),moyenne)
          der_pn(1)%DRT(:,1,i) = der_pn(1)%DRT(:,1,i)-moyenne
       ENDIF
       !===End of handling of mean first order derrivative pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)   = 0.d0
          un_p1 (:,4)   = 0.d0
          un_p1 (:,6)   = 0.d0
          der_pn(1)%DRT (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATE Velocity first order derrivative
       dt_un_p1(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATE Velocity first order derrivative

    ENDDO ! i = 1, m_max_c
    !===END Computation of first order derivatives

    !===Computation of 0th order derivatives
    !===Computation of rhs for un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         (un-dt/2.d0*(dt_un_p1-der_un(1)%DRT)-dt**2/12.d0*(dtt_un_p1-der_un(2)%DRT))/dt, &
         pn+dt*der_pn(1)%DRT- (dt**2)/2.d0*der_pn(2)%DRT, rotv_v2, rhs_gauss,0, density, tempn, concn)
    !===End Computation of rhs for un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        = vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+n2+1:n123)= vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        = vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+n2+1:n123)= vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, ur_c, ut_s, uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, ur_s, ut_c, uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for pressure
       pp_global_D(i)%DRL = 0.d0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for pressure

       !===Solve mass matrix for pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for pressure correction

       !===Pressure computation
       pn(:,:,i) = pn(:,:,i) +dt*der_pn(1)%DRT(:,:,i)-(dt**2)/2.d0*der_pn(2)%DRT(:,:,i)- lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Pressure computation

       !===UPDATES
       tps = user_time()
       !===Handling of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, pn(:,1,i),moyenne)
          pn(:,1,i) = pn(:,1,i)-moyenne
       ENDIF
       !===End of handling of mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)   = 0.d0
          un_p1 (:,4)   = 0.d0
          un_p1 (:,6)   = 0.d0
          pn    (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATES Velocity 0th order derrivative
       un(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES Velocity 0th order derrivative

    ENDDO ! i = 1, m_max_c
    !===END Computation of 0th order derivatives

    !===Update dtt_un un and dt_un once everyone has been compuded
    DO i = 1, m_max_c
       der_un(2)%DRT(:,:,i)=dtt_un_p1(:,:,i)
       der_un(1)%DRT(:,:,i)=dt_un_p1(:,:,i)
    END DO

    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'H1',  vv_mesh, list_mode, un)
       talk_to_me%div_L2  = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)/norm
       talk_to_me%weak_div_L2  = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, div)/norm
    END IF

    tps_tot = user_time() - tps_tot
  END SUBROUTINE update_ns_with_taylor
  !============================================

  SUBROUTINE update_ns_with_taylor_fourth(comm_one_d, time, vv_3_LA, pp_1_LA, vvz_per, pp_per, &
       dt, Re, lambda, list_mode, pp_mesh, vv_mesh, pn, der_pn, un, der_un, density, tempn, concn)
    !==============================
    USE def_type_mesh
    USE fem_M_axi
    USE fem_rhs_axi
    USE Dir_nodes_petsc
    USE periodic
    USE st_matrix
    USE solve_petsc
    USE boundary
    USE st_matrix
    USE input_data
    USE rhs_para_assembling
    USE tn_axi
    USE verbose
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8)                                        :: time, dt, Re, lambda
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: list_mode
    TYPE(mesh_type),                INTENT(IN)          :: pp_mesh, vv_mesh
    TYPE(petsc_csr_LA)                                  :: vv_3_LA, pp_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: vvz_per, pp_per
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)), INTENT(INOUT) :: un
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT) :: pn
    TYPE(dyn_real_array_three), DIMENSION(inputs%taylor_order-1)         :: der_un
    TYPE(dyn_real_array_three), DIMENSION(inputs%taylor_order-1)         :: der_pn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: density
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: tempn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: concn
    !===Saved variables
    INTEGER,                                       SAVE :: m_max_c
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: pp_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: pp_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3),              SAVE :: vv_js_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: vv_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: vel_global_D
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: my_petscworld_rank
    !===End saved variables
    INTEGER,          POINTER, DIMENSION(:)  :: pp_1_ifrom, vv_3_ifrom
    INTEGER                                  :: i, m, n, n1, n2, n3, n123
    INTEGER                                  :: nb_procs, code, nu_mat, mode
    REAL(KIND=8)                             :: moyenne
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: div
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)                  :: un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rotv_v1, rotv_v2, rotv_v3, rotv_v4
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))   :: uext, dttt_un_p1, dtt_un_p1, dt_un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np)                     :: vel_loc, vel_tot
    REAL(KIND=8) :: tps, tps_tot, tps_cumul, coeff, vloc, cfl, cfl_max, norm
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/
    !===Communicators for Petsc, in space and Fourier
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: vel_mat
    Mat,                        SAVE :: mass_mat
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Vec,                        SAVE :: pb_1,    pb_2,    px_1, px_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: vel_ksp
    KSP,                        SAVE :: mass_ksp
    !===END OF DECLARATION

    IF (once) THEN
       once = .FALSE.

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)
       !===CREATE PETSC VECTORS AND GHOSTS
       CALL create_my_ghost(vv_mesh,vv_3_LA,vv_3_ifrom)
       n = 3*vv_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
       CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
       CALL VecDuplicate(vx_3, vb_3_145, ierr)
       CALL VecDuplicate(vx_3, vb_3_236, ierr)

       CALL create_my_ghost(pp_mesh,pp_1_LA,pp_1_ifrom)
       n = pp_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(pp_1_ifrom), pp_1_ifrom, px_1, ierr)
       CALL VecGhostGetLocalForm(px_1, px_1_ghost, ierr)
       CALL VecDuplicate(px_1, pb_1, ierr)
       CALL VecDuplicate(px_1, pb_2, ierr)
       !===End CREATE PETSC VECTORS AND GHOSTS

       !===Number of modes on proc
       m_max_c = SIZE(list_mode)
       !===End umber of modes on proc

       !===PREPARE pp_mode_global_js_D ARRAY FOR PRESSURE
       CALL scalar_without_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
       ALLOCATE(pp_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(pp_global_D(i)%DRL(SIZE(pp_mode_global_js_D(i)%DIL)))
       END DO
       !===End PREPARE pp_mode_global_js_D  ARRAY FOR PRESSURE

       !===PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY
       CALL vector_glob_js_D(vv_mesh, list_mode, vv_3_LA, inputs%vv_list_dirichlet_sides, &
            vv_js_D, vv_mode_global_js_D)

       ALLOCATE(vel_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(vel_global_D(i)%DRL(SIZE(vv_mode_global_js_D(i)%DIL)))
       END DO
       !===END PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY

       !===ASSEMBLE MASS MATRIX
       CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, mass_mat, CLEAN=.FALSE.)
       CALL qs_diff_mass_scal_M (pp_mesh, pp_1_LA, 0.d0, 1.d0, 0.d0, 0, mass_mat)
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_matrix_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, mass_mat, pp_1_LA)
       END IF
       DO i = 1, m_max_c
          CALL Dirichlet_M_parallel(mass_mat,pp_mode_global_js_D(i)%DIL)
       END DO
       CALL init_solver(inputs%my_par_mass,mass_ksp,mass_mat,comm_one_d(1),&
            solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)
       !===END ASSEMBLE MASS MATRIX

       !===ASSEMBLING VELOCITY MATRICES
       ALLOCATE(vel_mat(2*m_max_c),vel_ksp(2*m_max_c))
       DO i = 1, m_max_c
          mode = list_mode(i)
          !===VELOCITY
          nu_mat = 2*i-1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_taylor (1, vv_3_LA, vv_mesh, one/Re, one/dt, & !HF new coeffs
               lambda, i, mode, vel_mat(nu_mat))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          nu_mat = nu_mat+1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_taylor (2, vv_3_LA, vv_mesh, one/Re, one/dt, & !HF new coeffs
               lambda, i, mode, vel_mat(nu_mat))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          !===End VELOCITY
       ENDDO
       !===End ASSEMBLING VELOCITY MATRICES

    ENDIF !===End of once

    !===Compute NL for dtt_un by FFT at Gauss points
    tps_tot = user_time()
    tps_cumul = 0
    tps = user_time()

    uext = un +(4*dt/3)*der_un(1)%DRT + (4*dt/3)**2/2*der_un(2)%DRT + (4*dt/3)**3/6*der_un(3)%DRT !===1st NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v1,inputs%precession)
    !===rotv_v1 is now a non linear term

    uext = un + dt*der_un(1)%DRT + dt**2/2*der_un(2)%DRT + dt**3/6*der_un(3)%DRT !===2nd NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v2,inputs%precession)
    !===rotv_v2 is now a non linear term

    uext = un +(dt/2)*der_un(1)%DRT + (dt/2)**2/2*der_un(2)%DRT + (dt/2)**3/6*der_un(3)%DRT  !===3rd NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v3,inputs%precession)
    !===rotv_v3 is now a non linear term

    uext = un  !===4th NL term
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v4,inputs%precession)
    !===rotv_v4 is now a non linear term
    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !===End Compute NL by FFT at Gauss points

    !===Computation of CFL
    IF (inputs%verbose_CFL) THEN
       vel_loc = 0.d0
       DO i = 1, m_max_c
          IF (list_mode(i)==0) THEN
             coeff = 1.d0
          ELSE
             coeff = 0.5d0
          END IF
          vel_loc = vel_loc + coeff*(un(:,1,i)**2+un(:,2,i)**2+un(:,3,i)**2 &
               +un(:,4,i)**2+un(:,5,i)**2+un(:,6,i)**2)
       END DO
       CALL MPI_COMM_SIZE(comm_one_d(2),nb_procs,code)
       CALL MPI_ALLREDUCE(vel_loc,vel_tot,vv_mesh%np,MPI_DOUBLE_PRECISION, MPI_SUM, comm_one_d(2), code)
       vel_tot = sqrt(vel_tot)
       cfl = 0.d0
       DO m = 1, vv_mesh%dom_me
          vloc = MAXVAL(vel_tot(vv_mesh%jj(:,m)))
          cfl = MAX(vloc*dt/MIN(vv_mesh%hloc(m),MAXVAL(vv_mesh%hm)),cfl)
       END DO
       CALL MPI_ALLREDUCE(cfl,cfl_max,1,MPI_DOUBLE_PRECISION, MPI_MAX, comm_one_d(1), code)
       talk_to_me%CFL=cfl_max
       talk_to_me%time=time
    END IF
    !===End Computation of CFL

    !===Computation of third order derivatives
    !===Computation of rhs for dttt_un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         der_un(3)%DRT/dt, der_pn(3)%DRT, 9*(9*rotv_v1-20*rotv_v2+16*rotv_v3-5*rotv_v4)/(5*dt**3), &
         rhs_gauss,3, density, tempn, concn)
    !===End Computation of rhs for dttt_un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)       = &
            (  vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)      &
            -3*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)   &
            +3*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt) &
            -  vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = &
            (  vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)      &
            -3*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)   &
            +3*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt) &
            -  vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n1+n2+1:n123)= &
            (  vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)      &
            -3*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)   &
            +3*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt) &
            -  vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)       = &
            (  vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)      &
            -3*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)   &
            +3*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt) &
            -  vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = &
            (  vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)      &
            -3*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)   &
            +3*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt) &
            -  vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n1+n2+1:n123)= &
            (  vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)      &
            -3*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)   &
            +3*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt) &
            -  vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(dt**3)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, dtt_ur_c, dtt_ut_s, dtt_uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, dtt_ur_s, dtt_ut_c, dtt_uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))
       !===dttt_un_p1 is temporarily in un_p1

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of dttt_velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for dttt_pressure
       pp_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for dttt_pressure

       !===Solve mass matrix for dttt_pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for dttt_pressure correction

       !===Compute 3rd order derrivative pressure for this mode
       der_pn(3)%DRT(:,:,i) =  der_pn(3)%DRT(:,:,i) - lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End 3rd order derrivative pressure

       !===UPDATES
       tps = user_time()
       !===Handling 2nd order derrivative of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, der_pn(3)%DRT(:,1,i),moyenne)
          der_pn(3)%DRT(:,1,i) = der_pn(3)%DRT(:,1,i)-moyenne
       ENDIF
       !===End of handling of 3rd order derrivative mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)     = 0.d0
          un_p1 (:,4)     = 0.d0
          un_p1 (:,6)     = 0.d0
          der_pn(3)%DRT (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATES 3rd order derrivative
       dttt_un_p1(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES 3rd order derrivative

    ENDDO
    !===END Computation of third order derivatives

    !===Computation of second order derivatives
    !===Computation of rhs for dtt_un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         der_un(2)%DRT/dt-dt/2*(dttt_un_p1-der_un(3)%DRT)/dt, der_pn(2)%DRT+dt*der_pn(3)%DRT, &
         (81*rotv_v1-140*rotv_v2+64*rotv_v3-5*rotv_v4)/(10*dt**2), rhs_gauss,2, density, tempn, concn)
    !===End Computation of rhs for dtt_un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        =                              &
            (2*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)     &
            -5*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)  &
            +4*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt)&
            -  vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =                              &
            (2*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)     &
            -5*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)  &
            +4*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt)&
            -  vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n1+n2+1:n123)=                              &
            (2*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)     &
            -5*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)  &
            +4*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt)&
            -  vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)

       vel_global_D(i)%DRL(1:n1)        =                              &
            (2*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)     &
            -5*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)  &
            +4*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt)&
            -  vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =                              &
            (2*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)     &
            -5*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)  &
            +4*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt)&
            -  vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n1+n2+1:n123)=                              &
            (2*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)     &
            -5*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)  &
            +4*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt)&
            -  vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(dt**2)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0 !===BUG JLG DEc 11 2019
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, dtt_ur_c, dtt_ut_s, dtt_uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, dtt_ur_s, dtt_ut_c, dtt_uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))
       !rk: dtt_un_p1 is temporarly in un_p1

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of dtt_velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for dtt_pressure
       pp_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for dtt_pressure

       !===Solve mass matrix for dtt_pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for dtt_pressure correction

       !===Compute 2nd order derrivative pressure for this mode
       der_pn(2)%DRT(:,:,i) =  der_pn(2)%DRT(:,:,i) + dt*der_pn(3)%DRT(:,:,i)- lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End 2nd order derrivative Pressure computation we will then take mean out

       !===UPDATES
       tps = user_time()
       !===Handling 2nd order derrivative of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, der_pn(2)%DRT(:,1,i),moyenne)
          der_pn(2)%DRT(:,1,i) = der_pn(2)%DRT(:,1,i)-moyenne
       ENDIF
       !===End of handling of 2nd order derrivative mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)    = 0.d0
          un_p1 (:,4)    = 0.d0
          un_p1 (:,6)    = 0.d0
          der_pn(2)%DRT (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATES 2nd order derrivative
       dtt_un_p1(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES 2nd order derrivative

    ENDDO
    !===END Computation of second order derivatives

    !===Computation of first order derivatives
    !===Computation of rhs for dt_un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         der_un(1)%DRT/dt-dt/2*(dtt_un_p1-der_un(2)%DRT)/dt-dt**2/12*(dttt_un_p1-der_un(3)%DRT)/dt, &
         der_pn(1)%DRT+dt*der_pn(2)%DRT-dt**2/2*der_pn(3)%DRT, (27*rotv_v1-32*rotv_v3+5*rotv_v4)/(20*dt), rhs_gauss,1, &
         density, tempn, concn)
    !===End Computation of rhs for dt_un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        =                               &
            (11*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)     &
            -18*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)  &
            + 9*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt)&
            - 2*vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =                               &
            (11*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)     &
            -18*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)  &
            + 9*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt)&
            - 2*vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n1+n2+1:n123)=                               &
            (11*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)     &
            -18*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)  &
            + 9*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt)&
            - 2*vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        =                               &
            (11*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time)     &
            -18*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-dt)  &
            + 9*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-2*dt)&
            - 2*vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n1+1:n1+n2)  =                               &
            (11*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time)     &
            -18*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-dt)  &
            + 9*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-2*dt)&
            - 2*vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n1+n2+1:n123)=                               &
            (11*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time)     &
            -18*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-dt)  &
            + 9*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-2*dt)&
            - 2*vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode,time-3*dt))/(6*dt)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, dt_ur_c, dt_ut_s, dt_uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))
       !===dt_un_p1 is temporarily in un_p1

       !===Solve system 2, dt_ur_s, dt_ut_c, dt_uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of dt_velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for dt_pressure
       pp_global_D(i)%DRL = 0.d0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for dt_pressure

       !===Solve mass matrix for dt_pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for dt_pressure correction

       !===Compute first order derrivative pressure for this mode
       !===Pressure first order derrivative computation
       der_pn(1)%DRT(:,:,i) = der_pn(1)%DRT(:,:,i) + dt*der_pn(2)%DRT(:,:,i) - dt**2/2*der_pn(3)%DRT(:,:,i) - lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End dt_Pressure computation

       !===UPDATES
       tps = user_time()
       !===Handling of mean first order derrivative pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, der_pn(1)%DRT(:,1,i),moyenne)
          der_pn(1)%DRT(:,1,i) = der_pn(1)%DRT(:,1,i)-moyenne
       ENDIF
       !===End of handling of mean first order derrivative pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)   = 0.d0
          un_p1 (:,4)   = 0.d0
          un_p1 (:,6)   = 0.d0
          der_pn(1)%DRT (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATE Velocity first order derrivative
       dt_un_p1(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATE Velocity first order derrivative

    ENDDO !=== i = 1, m_max_c
    !===END Computation of first order derivatives

    !===Computation of 0th order derivatives
    !===Computation of rhs for un at Gauss points for every mode
    CALL rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         (un-dt/2*(dt_un_p1-der_un(1)%DRT)-dt**2/12*(dtt_un_p1-der_un(2)%DRT))/dt, &
         pn+dt*der_pn(1)%DRT-dt**2/2*der_pn(2)%DRT+dt**3/6*der_pn(3)%DRT,rotv_v2,rhs_gauss,0, density, tempn, concn)
    !===End Computation of rhs for un

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       !===Periodic boundary conditions
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        = vv_exact(1,vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = vv_exact(4,vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+n2+1:n123)= vv_exact(5,vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        = vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+n2+1:n123)= vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
       vel_global_D(i)%DRL(n123+1:)     = 0.d0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
       tps = user_time()
       !===Solve system 1, ur_c, ut_s, uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !===Solve system 2, ur_s, ut_c, uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve linear system for momentum equation

       !===Assemble divergence of velocity in arrays pb_1, pb_2
       tps = user_time()
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===Boundary condition on axis for pressure
       pp_global_D(i)%DRL = 0.d0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for pressure

       !===Solve mass matrix for pressure correction
       tps = user_time()
       CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for pressure correction

       !===Pressure computation
       pn(:,:,i) = pn(:,:,i)+dt*der_pn(1)%DRT(:,:,i)-dt**2/2*der_pn(2)%DRT(:,:,i)+dt**3/6*der_pn(3)%DRT(:,:,i)-lambda*div(:,:,i)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Pressure computation

       !===UPDATES
       tps = user_time()
       !===Handling of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, pn(:,1,i),moyenne)
          pn(:,1,i) = pn(:,1,i)-moyenne
       ENDIF
       !===End of handling of mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2)   = 0.d0
          un_p1 (:,4)   = 0.d0
          un_p1 (:,6)   = 0.d0
          pn    (:,2,i) = 0.d0
       END IF
       !===Correction of zero mode

       !===UPDATES Velocity 0th order derrivative
       un(:,:,i)  = un_p1
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES Velocity 0th order derrivative

    ENDDO !=== i = 1, m_max_c
    !===END Computation of 0th order derivatives

    !===update dttt_un, dtt_un, and dt_un once everyone has been compuded
    DO i = 1, m_max_c
       der_un(3)%DRT(:,:,i)=dttt_un_p1(:,:,i)
       der_un(2)%DRT(:,:,i)= dtt_un_p1(:,:,i)
       der_un(1)%DRT(:,:,i)=  dt_un_p1(:,:,i)
    END DO

    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'H1',  vv_mesh, list_mode, un)
       talk_to_me%div_L2  = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)/norm
       talk_to_me%weak_div_L2  = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, div)/norm
    END IF

    tps_tot = user_time() - tps_tot
  END SUBROUTINE update_ns_with_taylor_fourth
  !============================================


  !===PRECESSION 28/07/09
  SUBROUTINE smb_cross_prod_gauss_sft_par(communicator,mesh,list_mode,V_in,V_out,precession_in)
    !=================================
    USE sft_parallele
    USE chaine_caractere
    USE input_data
    USE def_type_mesh
    IMPLICIT NONE

    TYPE(mesh_type),                INTENT(IN)  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: V_out
    LOGICAL,                        INTENT(IN)  :: precession_in
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: RotV, W, Om
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs, Omega_s
    REAL(KIND=8)   :: ray
    REAL(KIND=8)   :: tps, PI
    REAL(KIND=8), DIMENSION(3)                        :: temps
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE, SAVE ::  Omega
    LOGICAL, SAVE :: once=.TRUE.
    !===FOR FFT_PAR_CROSS_PROD_DCL
    INTEGER       :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm       :: communicator

    IF (once) THEN
       once = .FALSE.
       ALLOCATE(omega(mesh%np,6,SIZE(list_mode)))
       omega = 0.d0
       PI = ACOS(-1.d0)
       DO i=1, SIZE(list_mode)
          IF (list_mode(i) == 1) THEN
             !precession selon un axe penche d un angle angle_s_pi*PI
             omega(:,1,i) = inputs%taux_precession * SIN(inputs%angle_s_pi*PI)
             omega(:,4,i) = -inputs%taux_precession * SIN(inputs%angle_s_pi*PI)
          ELSE IF (list_mode(i) == 0) THEN
             !precession selon un axe penche d un angle angle_s_pi*PI
             omega(:,5,i) = inputs%taux_precession * COS(inputs%angle_s_pi*PI)
          ENDIF
       ENDDO
    ENDIF ! fin du once

    tps = user_time()
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
             Omega_s(:,k) = Omega(mesh%jj(:,m),k,i)
          END DO

          DO l = 1, mesh%gauss%l_G
             index = index + 1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Vs(:,1)*mesh%gauss%ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*mesh%gauss%ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*mesh%gauss%ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*mesh%gauss%ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*mesh%gauss%ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*mesh%gauss%ww(:,l))
             !-----------------vecteur rotation sur les points de Gauss---------------------------
             Om(index,1,i) = SUM(Omega_s(:,1)*mesh%gauss%ww(:,l))
             Om(index,3,i) = SUM(Omega_s(:,3)*mesh%gauss%ww(:,l))
             Om(index,5,i) = SUM(Omega_s(:,5)*mesh%gauss%ww(:,l))

             Om(index,2,i) = SUM(Omega_s(:,2)*mesh%gauss%ww(:,l))
             Om(index,4,i) = SUM(Omega_s(:,4)*mesh%gauss%ww(:,l))
             Om(index,6,i) = SUM(Omega_s(:,6)*mesh%gauss%ww(:,l))
             !-----------------rotational sur les points de Gauss---------------------------
             !coeff sur les cosinus
             RotV(index,1,i) = mode/ray*W(index,6,i) &
                  -SUM(Vs(:,3)*dw_loc(2,:))
             RotV(index,4,i) =          SUM(Vs(:,2)*dw_loc(2,:)) &
                  -SUM(Vs(:,6)*dw_loc(1,:))
             RotV(index,5,i) =    1/ray*W(index,3,i) &
                  +SUM(Vs(:,3)*dw_loc(1,:)) &
                  -mode/ray*W(index,2,i)
             !coeff sur les sinus
             RotV(index,2,i) =-mode/ray*W(index,5,i) &
                  -SUM(Vs(:,4)*dw_loc(2,:))
             RotV(index,3,i) =         SUM(Vs(:,1)*dw_loc(2,:)) &
                  -SUM(Vs(:,5)*dw_loc(1,:))
             RotV(index,6,i) =    1/ray*W(index,4,i) &
                  +SUM(Vs(:,4)*dw_loc(1,:))&
                  +mode/ray*W(index,1,i)
          ENDDO
       ENDDO
    END DO

    !JLG, FL, July 23 2010, There was a Bug here
    IF (.NOT.precession_in) THEN
       Om=0.d0
    END IF
    !
    ! cas PRECESSION 28/07/09
    !take care to the -F= 2 Om x U term in LHS
    RotV = RotV + 2.d0 * Om
    temps = 0


    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(RotV,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    CALL FFT_PAR_CROSS_PROD_DCL(communicator, RotV, W, V_out, nb_procs, bloc_size, m_max_pad, temps)
    tps = user_time() - tps

  END SUBROUTINE smb_cross_prod_gauss_sft_par

  SUBROUTINE Moy(communicator,mesh,p,RESLT)
    !===========================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(:)  ,   INTENT(IN)  :: p
    REAL(KIND=8)                ,   INTENT(OUT) :: RESLT
    REAL(KIND=8)                                :: vol_loc, vol_out, r_loc, r_out
    INTEGER ::  m, l , i , ni, code
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: j_loc
    REAL(KIND=8)   :: ray
    MPI_Comm                                    :: communicator
    r_loc = 0.d0
    vol_loc = 0.d0

    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)
       DO l = 1, mesh%gauss%l_G
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = j_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          r_loc = r_loc +  SUM(p(j_loc(:))*mesh%gauss%ww(:,l))*ray*mesh%gauss%rj(l,m)
          vol_loc = vol_loc + ray*mesh%gauss%rj(l,m)

       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(r_loc,r_out,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    CALL MPI_ALLREDUCE(vol_loc,vol_out,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    RESLT = r_out / vol_out

  END SUBROUTINE Moy

  SUBROUTINE rhs_ns_gauss_3x3_taylor(vv_mesh, pp_mesh, communicator, list_mode, time, V1m, pn, rotv_v, &
       rhs_gauss, der_ord, density, tempn, concn)
    !=================================
    !RHS for Navier-Stokes_taylor : RHS for (d/dt)**der_ord un
    USE def_type_mesh
    USE my_util
    USE input_data
    USE fem_rhs_axi
    USE sft_parallele
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type)                                        :: vv_mesh, pp_mesh
    REAL(KIND=8),                               INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: rotv_v
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: V1m
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: pn
    INTEGER,      DIMENSION(:),                 INTENT(IN) :: list_mode
    INTEGER,                                    INTENT(IN) :: der_ord
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)                                   :: fs, ft
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8), DIMENSION(vv_mesh%np,6)                        :: ff, imposed_vel
    REAL(KIND=8), DIMENSION(vv_mesh%np,2)                        :: P
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(vv_mesh%dom_me*vv_mesh%gauss%l_G,6)  :: imposed_vel_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%dom_me*vv_mesh%gauss%l_G,6,SIZE(list_mode)) :: fp, rhs_gauss_penal
    REAL(KIND=8), DIMENSION(2,vv_mesh%gauss%l_G*vv_mesh%dom_me)  :: rr_gauss
    REAL(KIND=8) :: ray
    INTEGER :: m, l , i, k, index, TYPE
    INTEGER :: nb_procs, m_max_pad, bloc_size
#include "petsc/finclude/petsc.h"
    PetscErrorCode                   :: ierr
    MPI_Comm                         :: communicator

    IF ((der_ord-0)*(der_ord-1)*(der_ord-2)*(der_ord-3)/=0)THEN
       CALL error_Petsc('BUG in rhs_ns_gauss_3x3_taylor, der_ord not correct')
    END IF
    DO i = 1, SIZE(list_mode)
       DO k= 1, 6 !===external forces
          IF (der_ord==0) THEN
             ff(:,k) = source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time, &
                  inputs%Re,'ns', density, tempn, concn)
          ELSEIF (der_ord==1) THEN
             ff(:,k) = &
                  (11*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time, &
                  inputs%Re,'ns', density, tempn, concn) &
                  -18*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  + 9*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-2*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  - 2*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-3*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn))/(6*inputs%dt)
          ELSEIF (der_ord==2) THEN
             ff(:,k) = &
                  (2*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time, &
                  inputs%Re,'ns', density, tempn, concn) &
                  -5*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  +4*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-2*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  -  source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-3*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn))/(inputs%dt**2)
          ELSEIF (der_ord==3) THEN
             ff(:,k) = &
                  (  source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time, &
                  inputs%Re,'ns', density, tempn, concn) &
                  -3*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  +3*source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-2*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn) &
                  -  source_in_NS_momentum(k, vv_mesh%rr, list_mode(i),i, time-3*inputs%dt, &
                  inputs%Re,'ns', density, tempn, concn))/(inputs%dt**3)
          END IF
       END DO
       DO k = 1, 2 !===pressure
          !===JLG+HF July 24th, 2019
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
       ENDDO

       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)

          DO l = 1, vv_mesh%gauss%l_G !===compute values on gauss pts
             index  = index +1
             dw_loc = vv_mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
             rr_gauss(1,index) = ray
             rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))

             !=== u0(1,:) <--> f(r,m,c)
             fs(1) = SUM(ff(j_loc,1) * vv_mesh%gauss%ww(:,l))    !external forcing
             ft(1) = SUM(V1m(j_loc,1,i) * vv_mesh%gauss%ww(:,l)) !inertia
             fp(index,1,i) = -SUM(P(j_loc,1)*dw_loc(1,:))        !pressure term
             !=== u0(2,:) <--> f(r,m,s)
             fs(2) = SUM(ff(j_loc,2) * vv_mesh%gauss%ww(:,l))
             ft(2) = SUM(V1m(j_loc,2,i) * vv_mesh%gauss%ww(:,l))
             fp(index,2,i) = -SUM(P(j_loc,2)*dw_loc(1,:))
             !=== u0(3,:) <--> f(th,m,c)
             fs(3) = SUM(ff(j_loc,3) * vv_mesh%gauss%ww(:,l))
             ft(3) = SUM(V1m(j_loc,3,i) * vv_mesh%gauss%ww(:,l))
             fp(index,3,i) = -SUM(P(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !=== u0(4,:) <--> f(th,m,s)
             fs(4) = SUM(ff(j_loc,4) * vv_mesh%gauss%ww(:,l))
             ft(4) = SUM(V1m(j_loc,4,i) * vv_mesh%gauss%ww(:,l))
             fp(index,4,i) = SUM(P(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !=== u0(5,:) <--> f(z,m,c)
             fs(5) = SUM(ff(j_loc,5) * vv_mesh%gauss%ww(:,l))
             ft(5) = SUM(V1m(j_loc,5,i) * vv_mesh%gauss%ww(:,l))
             fp(index,5,i) = -SUM(P(j_loc,1)*dw_loc(2,:))
             !=== u0(6,:) <--> f(z,m,s)
             fs(6) = SUM(ff(j_loc,6) * vv_mesh%gauss%ww(:,l))
             ft(6) = SUM(V1m(j_loc,6,i) * vv_mesh%gauss%ww(:,l))
             fp(index,6,i) = -SUM(P(j_loc,2)*dw_loc(2,:))

             rhs_gauss(index,:,i) =  (ft+fs-rotv_v(index,:,i))

          ENDDO
       ENDDO
       IF (inputs%if_ns_penalty) THEN
          IF(inputs%if_impose_vel_in_solids) THEN
             IF (list_mode(i)==0) THEN
                IF (der_ord==0) THEN
                   imposed_vel(:,:) = imposed_velocity_by_penalty(vv_mesh%rr,time)/(inputs%dt)
                   index = 0
                   DO m = 1, vv_mesh%dom_me
                      j_loc = vv_mesh%jj(:,m)
                      DO l = 1, vv_mesh%gauss%l_G
                         index  = index +1
                         DO TYPE = 1, 6
                            imposed_vel_gauss(index,TYPE) = SUM(imposed_vel(j_loc,TYPE) &
                                 * vv_mesh%gauss%ww(:,l))
                         END DO
                      END DO
                   END DO
                   rhs_gauss(:,:,i) = rhs_gauss(:,:,i) - imposed_vel_gauss(:,:)
                ELSE IF (der_ord==1) THEN
                   !===CN HF adapt to 4th order
                   !   imposed_vel(:,:) = (imposed_velocity_by_penalty(vv_mesh%rr,time+inputs%dt) &
                   !        - imposed_velocity_by_penalty(vv_mesh%rr,time-inputs%dt))/(2.d0*inputs%dt**2)
                   imposed_vel(:,:) = &
                        (-2*imposed_velocity_by_penalty(vv_mesh%rr,time-3*inputs%dt)&
                        +9*imposed_velocity_by_penalty(vv_mesh%rr,time-2*inputs%dt)&
                        -18*imposed_velocity_by_penalty(vv_mesh%rr,time-1*inputs%dt)&
                        +11*imposed_velocity_by_penalty(vv_mesh%rr,time))/(6*inputs%dt**2)
                   !===CN HF adapt to 4th order
                   index = 0
                   DO m = 1, vv_mesh%dom_me
                      j_loc = vv_mesh%jj(:,m)
                      DO l = 1, vv_mesh%gauss%l_G
                         index  = index +1
                         DO TYPE = 1, 6
                            imposed_vel_gauss(index,TYPE) = SUM(imposed_vel(j_loc,TYPE) &
                                 * vv_mesh%gauss%ww(:,l))
                         END DO
                      END DO
                   END DO
                   rhs_gauss(:,:,i) = rhs_gauss(:,:,i) - imposed_vel_gauss(:,:)
                ELSE IF (der_ord==2) THEN
                   !===CN HF adapt to 4th order: nothing to do
                   imposed_vel(:,:) = &
                        (-1*imposed_velocity_by_penalty(vv_mesh%rr,time-3*inputs%dt)&
                        +4*imposed_velocity_by_penalty(vv_mesh%rr,time-2*inputs%dt)&
                        -5*imposed_velocity_by_penalty(vv_mesh%rr,time-1*inputs%dt)&
                        +2*imposed_velocity_by_penalty(vv_mesh%rr,time))/(inputs%dt**3)
                   !===CN HF adapt to 4th order: nothing to do
                   index = 0
                   DO m = 1, vv_mesh%dom_me
                      j_loc = vv_mesh%jj(:,m)
                      DO l = 1, vv_mesh%gauss%l_G
                         index  = index +1
                         DO TYPE = 1, 6
                            imposed_vel_gauss(index,TYPE) = SUM(imposed_vel(j_loc,TYPE) &
                                 * vv_mesh%gauss%ww(:,l))
                         END DO
                      END DO
                   END DO
                   rhs_gauss(:,:,i) = rhs_gauss(:,:,i) - imposed_vel_gauss(:,:)
                   !===CN HF adapt to 4th order
                ELSE IF (der_ord==3) THEN
                   imposed_vel(:,:) = &
                        (-1*imposed_velocity_by_penalty(vv_mesh%rr,time-3*inputs%dt)&
                        +3*imposed_velocity_by_penalty(vv_mesh%rr,time-2*inputs%dt)&
                        -3*imposed_velocity_by_penalty(vv_mesh%rr,time-1*inputs%dt)&
                        +1*imposed_velocity_by_penalty(vv_mesh%rr,time))/(inputs%dt**4)
                   index = 0
                   DO m = 1, vv_mesh%dom_me
                      j_loc = vv_mesh%jj(:,m)
                      DO l = 1, vv_mesh%gauss%l_G
                         index  = index +1
                         DO TYPE = 1, 6
                            imposed_vel_gauss(index,TYPE) = SUM(imposed_vel(j_loc,TYPE) &
                                 * vv_mesh%gauss%ww(:,l))
                         END DO
                      END DO
                   END DO
                   rhs_gauss(:,:,i) = rhs_gauss(:,:,i) - imposed_vel_gauss(:,:)
                   !===CN HF adapt to 4th order
                ENDIF
             END IF
          END IF
       END IF
    END DO

    IF (inputs%if_ns_penalty) THEN
       IF(inputs%if_impose_vel_in_solids) THEN
          CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
          m_max_pad = 3*SIZE(list_mode)*nb_procs/2
          bloc_size = SIZE(rhs_gauss,1)/nb_procs+1

          CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator, penal_in_real_space, vv_mesh, &
               rhs_gauss, rhs_gauss_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
          DO i = 1, SIZE(list_mode)
             IF (list_mode(i)==0) THEN
                rhs_gauss(:,:,i) = rhs_gauss_penal(:,:,i) + imposed_vel_gauss(:,:)
             ELSE
                rhs_gauss(:,:,i) = rhs_gauss_penal(:,:,i)
             END IF
          END DO
       END IF
    END IF

    rhs_gauss = rhs_gauss + fp

  END SUBROUTINE rhs_ns_gauss_3x3_taylor

  SUBROUTINE qs_diff_mass_vect_3x3_taylor(type_op, LA, mesh, visco, mass, lambda, i_mode, mode, matrix)
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    USE my_util
    USE input_data
    IMPLICIT NONE

    INTEGER     ,                 INTENT(IN)    :: type_op, mode, i_mode
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass, lambda
    TYPE(petsc_csr_la)                          :: LA
    TYPE(mesh_type), TARGET                     :: mesh

    INTEGER :: k, l, m, ni, nj, i, j, np, ki, kj, k_max, ls, ms, n_w, n_ws
    REAL(KIND=8) :: xij, viscolm, div_penal
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij, cij, dij, eij, fij
    REAL(KIND=8) :: ray, eps1, eps2, z
    REAL(KIND=8) :: two = 2.d0
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w,3*mesh%gauss%n_w)   :: mat_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                    :: idxn, jdxn
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_ws,2*mesh%gauss%n_ws) :: mat_loc_s
    INTEGER,      DIMENSION(2*mesh%gauss%n_ws)                   :: idxn_s, jdxn_s
    INTEGER                                                      :: ix, jx, iglob, jglob
    INTEGER,      DIMENSION(mesh%gauss%n_w)                      :: jj_loc
    REAL(KIND=8) :: viscomode, hm
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    np   = SIZE(mesh%rr,2)
    n_w  = mesh%gauss%n_w
    n_ws = mesh%gauss%n_ws

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 3) THEN
       !cas du laplacien scalaire
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 !Structure scalaire
    ELSE
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1
       CALL error_petsc('probleme de type d''operateur')
    ENDIF


    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, mesh%me
       jj_loc = mesh%jj(:,m)

       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm  = visco*mesh%gauss%rj(l,m)
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)
          hm=MIN(mesh%hm(i_mode),mesh%hloc(m))!WRONG choice
          viscomode = visco*mesh%gauss%rj(l,m)

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)

             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.d0
                DO k = 1, mesh%gauss%k_d
                   xij =  xij + mesh%gauss%dw(k,nj,l,m) * mesh%gauss%dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                z = ray * viscolm* xij    &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                cij(ni,nj) =  cij(ni,nj) + z
                aij(ni,nj) =  aij(ni,nj) + z + viscomode*eps1*wwprod(ni,nj,l)/ray
                !blocs couplant
                bij(ni,nj) = bij(ni,nj) + eps2*viscomode*2*mode*wwprod(ni,nj,l)/ray
             ENDDO
          ENDDO

       ENDDO

       mat_loc = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxn(ix) = iglob-1
             DO kj = 1, 3
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   jdxn(jx) = jglob-1
                   IF ((ki .LT. 3) .AND. (kj .LT. 3)) THEN
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
                      ELSE
                         mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
                      END IF
                   ELSE ! ki=3 OR kj=3
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

       !+++---------------------
!!$       IF (type_op==3) CYCLE
       !==Calcul de visco (grad u)T . (grad v)
       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       dij = 0.d0
       eij = 0.d0
       fij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm  = visco*mesh%gauss%rj(l,m)*ray
          !HF 2019/04/11
          div_penal = lambda*mesh%gauss%rj(l,m)*ray
          !HF 2019/04/11

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)
             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)
                aij(ni,nj) = aij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(1,ni,l,m)*mesh%gauss%dw(1,nj,l,m) + wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(mesh%gauss%dw(1,nj,l,m) &
                     + mesh%gauss%ww(nj,l)/ray)
                bij(ni,nj) = bij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%ww(ni,l)*mesh%gauss%dw(1,nj,l,m)/ray+eps2*mode*wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(eps2*(mode/ray)*mesh%gauss%ww(nj,l))
                cij(ni,nj) = cij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(1,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*mesh%gauss%dw(2,nj,l,m)
                dij(ni,nj) = dij(ni,nj) &
                     + viscolm*(-mesh%gauss%dw(1,ni,l,m)*mesh%gauss%ww(nj,l)/ray+(mode/ray)**2*wwprod(ni,nj,l) &
                     -mesh%gauss%dw(1,nj,l,m)*mesh%gauss%ww(ni,l)/ray) &
                     + div_penal*wwprod(ni,nj,l)*(mode/ray)**2
                eij(ni,nj) = eij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%dw(2,ni,l,m)*mesh%gauss%ww(nj,l)/ray) &
                     + div_penal*eps2*(mode/ray)*mesh%gauss%ww(ni,l)*mesh%gauss%dw(2,nj,l,m)
                fij(ni,nj) = fij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m))
             END DO
          END DO
       END DO
       !++++++++++++++++
       mat_loc=0.d0
       idxn=0
       jdxn=0
       DO ni = 1, n_w
          DO ki = 1, 3
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w + ni
             idxn(ix) = iglob-1
          END DO
       END DO
       jdxn=idxn

       DO ni = 1, n_w
          DO nj = 1, n_w
             !=== Line i 1 (Vr)
             ix = ni
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)

             !=== Line i 2 (Vt)
             ix = ni+n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + dij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(ni,nj)

             !=== Line i 3 (Vz)
             ix = ni+2*n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(nj,ni)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + fij(ni,nj)
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

    ENDDO
    !== Fin du Calcul de visco (grad u)T . (grad v)
    !++++++++++++++++------

    IF (inputs%vv_nb_dirichlet_normal_velocity>0) THEN
       !===Surface terms
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)- inputs%vv_list_dirichlet_normal_velocity_sides)).NE.0) CYCLE
          aij = 0.d0
          bij = 0.d0
          cij = 0.d0
          dij = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
             IF (ray.LT.1.d-10) CYCLE
             z = two*mesh%gauss%rjs(ls,ms)*ray*visco

             DO ni = 1, mesh%gauss%n_ws
                DO nj = 1, mesh%gauss%n_ws
                   aij(ni,nj) = aij(ni,nj) - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   bij(ni,nj) = bij(ni,nj) - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                   cij(ni,nj) = cij(ni,nj) - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   dij(ni,nj) = dij(ni,nj) - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                END DO
             END DO
          END DO

          !++++++++++++++++
          !=== In the following loops, ki=1 for Vr and ki=2 for Vz
          !===   ==> 2ki-1 = 1 for Vr, 2ki-1 = 3 for Vz   (LA%loc_to_glob(2*ki-1,i))
          mat_loc_s = 0.d0
          idxn_s    = 0
          jdxn_s    = 0
          DO ki = 1, 2
             DO ni = 1, n_ws
                i = mesh%jjs(ni,ms)
                iglob = LA%loc_to_glob(2*ki-1,i)
                ix = (ki-1)*n_ws+ni
                idxn_s(ix) = iglob-1
                DO kj = 1, 2
                   DO nj = 1, n_ws
                      j = mesh%jjs(nj,ms)
                      jglob = LA%loc_to_glob(2*kj-1,j)
                      jx = (kj-1)*n_ws+nj
                      jdxn_s(jx) = jglob-1
                      IF ((ki == 1) .AND. (kj == 1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + aij(ni,nj) + aij(nj,ni)
                      ELSE IF ((ki == 1) .AND. (kj==2)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + bij(ni,nj) + cij(nj,ni)
                      ELSE IF ((ki == 2) .AND. (kj==1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + cij(ni,nj) + bij(nj,ni)
                      ELSE
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + dij(ni,nj) + dij(nj,ni)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          CALL MatSetValues(matrix, 2*n_ws, idxn_s, 2*n_ws, jdxn_s, mat_loc_s, ADD_VALUES, ierr)
       END DO
    END IF

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)


  END SUBROUTINE qs_diff_mass_vect_3x3_taylor
END MODULE update_taylor_navier_stokes
