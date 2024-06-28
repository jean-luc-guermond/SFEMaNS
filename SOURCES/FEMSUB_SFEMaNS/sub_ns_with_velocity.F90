!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_ns_with_u
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
  USE subroutine_ns_with_m ! MODIFICATION: for smb_explicit_diffu_sym
  PUBLIC :: moy, BDF2_ns_stress_bc_with_u
  PRIVATE
CONTAINS

  SUBROUTINE BDF2_ns_stress_bc_with_u(comm_one_d, time, vv_3_LA, pp_1_LA, vvz_per, pp_per, dt, Re, list_mode, pp_mesh, &
       vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un,  &
!       chmp_mag, Bn_p2, density, tempn, concn)
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
       chmp_mag, Bn_p2, density, tempn, concn, visco_entro_sym_grad_u)
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
    !==============================
    USE def_type_mesh
    USE fem_M_axi
    USE fem_rhs_axi
    USE fem_tn_axi
    USE Dir_nodes_petsc
    USE periodic
    USE st_matrix
    USE solve_petsc
    USE dyn_line
    USE boundary
    USE chaine_caractere
    USE sub_plot
    USE st_matrix
    USE input_data
    USE rhs_gauss_computing
    USE rhs_para_assembling
    USE tn_axi
    USE verbose
    USE entropy_viscosity
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8)                                   :: time, dt, Re
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    TYPE(mesh_type),                INTENT(IN)     :: pp_mesh, vv_mesh
    TYPE(petsc_csr_LA)                             :: vv_3_LA, pp_1_LA
    TYPE(periodic_type),            INTENT(IN)     :: vvz_per, pp_per
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: incpn_m1, incpn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: pn_m1, pn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: density
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: concn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: tempn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: chmp_mag, Bn_p2
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(INOUT):: visco_entro_sym_grad_u
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)):: chmp_mag_aux
    !===Saved variables
    INTEGER,                                       SAVE :: m_max_c
    !INTEGER,     DIMENSION(:),   POINTER,          SAVE :: pp_js_D
    !INTEGER,     DIMENSION(:),   POINTER,          SAVE :: pp_js_axis_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: pp_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: pp_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3),              SAVE :: vv_js_D
    !INTEGER,     DIMENSION(:),   POINTER,          SAVE :: vv_js_axis_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: vv_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: vel_global_D
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: my_petscworld_rank
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
!    REAL(KIND=8), DIMENSION(:,:,:,:),ALLOCATABLE,  SAVE :: visco_entro_sym_grad_u
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
    !===End saved variables

    !----------Declaration sans save-------------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: pp_1_ifrom, vv_3_ifrom

    INTEGER                                  :: i, k, m, n, n1, n2, n3, n123
    INTEGER                                  :: nb_procs, code, nu_mat, mode
    REAL(KIND=8)                             :: moyenne
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: div
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2)   :: pn_p1, phi
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)   :: un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6, SIZE(list_mode))   :: un_m2
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rotv_v, rotb_b, rotb_b_aux
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: mag_force
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)) :: uext
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode)) :: nu_tilde ! MODIFICATION: variable part of the kinematic viscosity
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: two_nu_tilde_grads_u_ext ! MODIFICATION: 2 times variable part of nu times symmetric gradient of u
    REAL(KIND=8), DIMENSION(vv_mesh%np) :: vel_loc, vel_tot
    REAL(KIND=8)   :: tps, tps_tot, tps_cumul, coeff, vloc, cfl, cfl_max, norm
    !April 17th 2008, JLG
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/

    !Communicators for Petsc, in space and Fourier------------------------------
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: vel_mat
    Mat, DIMENSION(:), POINTER, SAVE :: press_mat
    Mat,                        SAVE :: mass_mat, mass_mat0
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Vec,                        SAVE :: pb_1,    pb_2,    px_1, px_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: vel_ksp, press_ksp
    KSP,                        SAVE :: mass_ksp, mass_ksp0
    !------------------------------END OF DECLARATION--------------------------------------

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
       !===ATTENTION pressure BCs are no longer implemented
       !===JLG JUne 9 2017
       !CALL scalar_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
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
          IF (list_mode(i)==0) CYCLE
          CALL Dirichlet_M_parallel(mass_mat,pp_mode_global_js_D(i)%DIL)
       END DO
       CALL init_solver(inputs%my_par_mass,mass_ksp,mass_mat,comm_one_d(1),&
            solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)

       IF (MINVAL(list_mode)==0) THEN
          CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, mass_mat0, CLEAN=.FALSE.)
          CALL qs_diff_mass_scal_M (pp_mesh, pp_1_LA, 0.d0, 1.d0, 0.d0, 0, mass_mat0)
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, mass_mat0, pp_1_LA)
          END IF
          DO i = 1, m_max_c
             IF (list_mode(i).NE.0) CYCLE
             CALL Dirichlet_M_parallel(mass_mat0,pp_mode_global_js_D(i)%DIL)
             CALL init_solver(inputs%my_par_mass,mass_ksp0,mass_mat0,comm_one_d(1),&
                  solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)
          END DO
       END IF
       !===END ASSEMBLE MASS MATRIX

       !===ASSEMBLING VELOCITY MATRICES
!TEST LC 2024/02/26 (now initialize or read in initialization.F90)
!       ALLOCATE(visco_entro_sym_grad_u(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)))
!       visco_entro_sym_grad_u = 0.d0
!TEST LC 2024/02/26
       ALLOCATE(vel_mat(2*m_max_c),vel_ksp(2*m_max_c))
       ALLOCATE(press_mat(m_max_c),press_ksp(m_max_c))
       DO i = 1, m_max_c
          mode = list_mode(i)
          !===VELOCITY
          nu_mat = 2*i-1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_divpenal (1, vv_3_LA, vv_mesh, one/Re, three/(2*dt), &
               inputs%LES_coeff1, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          nu_mat = nu_mat+1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          CALL qs_diff_mass_vect_3x3_divpenal (2, vv_3_LA, vv_mesh, one/Re, three/(2*dt), &
               inputs%LES_coeff1, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          !===End VELOCITY

          !===PRESSURE
          CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, press_mat(i), clean=.FALSE.)
          !JLG Jan 2014 (regularize pressure matrix)
          !CALL qs_diff_mass_scal_M(pp_mesh, pp_1_LA, one, zero, zero, mode, press_mat(i))
          CALL qs_diff_mass_scal_M(pp_mesh, pp_1_LA, one, 1.d-10, zero, mode, press_mat(i))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, press_mat(i), pp_1_LA)
          END IF
          CALL Dirichlet_M_parallel(press_mat(i),pp_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_pp,press_ksp(i),press_mat(i),comm_one_d(1),&
               solver=inputs%my_par_pp%solver,precond=inputs%my_par_pp%precond)
          !===End PRESSURE
       ENDDO
       !===End ASSEMBLING VELOCITY MATRICES
    ENDIF !===End of once

    !===Compute NL by FFT at Gauss points
    tps_tot = user_time()
    tps_cumul = 0
    tps = user_time()
    uext = 2*un-un_m1
    CALL smb_cross_prod_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,uext,rotv_v,inputs%precession)

    IF (inputs%type_pb=='mhd' .OR. inputs%type_pb=='mhs') THEN
       !===Compute Lorentz force if mhd in quasi-static limit
       IF (inputs%if_quasi_static_approx) THEN
          DO i = 1, m_max_c
             mode = list_mode(i)
             chmp_mag_aux(:,:,i) = H_B_quasi_static('H', vv_mesh%rr, mode)
          END DO
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,chmp_mag_aux,Bn_p2,rotb_b_aux)
          DO i = 1, m_max_c
             mode = list_mode(i)
             chmp_mag_aux(:,:,i)  = H_B_quasi_static('B', vv_mesh%rr, mode)
          END DO
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,chmp_mag,chmp_mag_aux,rotb_b)
          rotb_b = rotb_b + rotb_b_aux
       ELSE !===Compute Lorentz force if mhd
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,chmp_mag,Bn_p2,rotb_b)
       END IF
       rotv_v = rotv_v - inputs%coeff_lorentz*rotb_b
    END IF

    IF (inputs%type_pb=='fhd') THEN
       !===Computate magnetic force if fhd
       IF (inputs%if_helmholtz_force) THEN
          CALL smb_helmholtz_force_gauss_fft_par(comm_one_d(2),vv_mesh,list_mode,tempn,chmp_mag,mag_force)
       ELSE
          CALL smb_kelvin_force_gauss_fft_par(comm_one_d(2),vv_mesh,list_mode,tempn,chmp_mag,mag_force)
       END IF
       rotv_v = rotv_v - mag_force
    END IF

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
          !cfl = MAX(vloc*dt/MIN(vv_mesh%hloc(m),0.5d0/inputs%m_max),cfl)
          cfl = MAX(vloc*dt/MIN(vv_mesh%hloc(m),MAXVAL(vv_mesh%hm)),cfl) !(JLG April 7 2017)
       END DO
       CALL MPI_ALLREDUCE(cfl,cfl_max,1,MPI_DOUBLE_PRECISION, MPI_MAX, comm_one_d(1), code)
       talk_to_me%CFL=cfl_max
       talk_to_me%time=time
    END IF
    !===End Computation of CFL

    !===Computation of rhs at Gauss points for every mode
    CALL rhs_ns_gauss_3x3(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
         (4*un-un_m1)/(2*inputs%dt), pn, (4.d0*incpn-incpn_m1)/3.d0, &
         rotv_v, rhs_gauss, density, tempn, concn)
    !===End Computation of rhs

    !===Computation of 2 * nu_tilde(T) * grad^s u
    IF (inputs%if_variable_visco) THEN
       nu_tilde = tempn
       CALL smb_nu_tilde_fft_par(comm_one_d(2), list_mode, nu_tilde)
       IF (.NOT. (inputs%verbose_CFL)) THEN
          CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       END IF
       CALL smb_explicit_diffu_sym(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            nu_tilde, 2 * un - un_m1, two_nu_tilde_grads_u_ext)
    END IF
    !===End Computation of 2 * nu_tilde(T) * grad^s u

    DO i = 1, m_max_c
       tps = user_time()
       mode = list_mode(i)

       !===Compute phi
       pn_p1(:,:) = pn(:,:,i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       IF (inputs%LES) THEN
          CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236, &
               opt_tensor=-visco_entro_sym_grad_u(:,:,:,i))
       ELSE IF (inputs%if_variable_visco) THEN
          CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236, &
               opt_tensor=-two_nu_tilde_grads_u_ext(:,:,:,i)) ! MODIFICATION: variable part of nu term in the RHS, constant part of nu (in data) in the LHS
       ELSE
          CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236)
       END IF
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
       vel_global_D(i)%DRL(n123+1:)     = 0.D0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        = vv_exact(2,vv_mesh%rr(:,vv_js_D(1)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = vv_exact(3,vv_mesh%rr(:,vv_js_D(2)%DIL), mode, time)
       vel_global_D(i)%DRL(n1+n2+1:n123)= vv_exact(6,vv_mesh%rr(:,vv_js_D(3)%DIL), mode, time)
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
       CALL qs_01_div_hybrid_generic(inputs%type_fe_velocity,vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2) !===JLG Dec 9 2019
       !CALL qs_01_div_hybrid_2006(vv_mesh, pp_mesh, pp_1_LA, mode, un_p1, pb_1, pb_2)
       !===End assembling; pb_1, and pb_2 are petsc vectors for the rhs divergence

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_1, pp_1_LA)
          CALL periodic_rhs_petsc(pp_per%n_bord, pp_per%list, pp_per%perlist, pb_2, pp_1_LA)
       END IF

       !===ATENTION BCs are no longer implemented for pressure
       !===Boundary condition on axis for pressure
       pp_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for pressure

       !===Solve -LAP(PH3) = -3/(2*dt)*DIV(un_p1)
       CALL solver(press_ksp(i),pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_pp%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,phi(:,1))

       CALL solver(press_ksp(i),pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_pp%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,phi(:,2))

       phi = -phi*(1.5d0/dt)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Solve -LAP(PH3) = -3/(2*dt)*DIV(un_p1)

       !===Solve mass matrix for pressure correction
       tps = user_time()
       IF (mode==0) THEN
          CALL solver(mass_ksp0,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       ELSE
          CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       END IF
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))
       IF (mode==0) THEN
          CALL solver(mass_ksp0,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       ELSE
          CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       END IF
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))
       !===End solve mass matrix for pressure correction

       !===Pressure computation
       DO k=1, 2
          !pn_p1(:,k) = pn_p1(:,k) + phi(:,k) - div(:,k,i)/Re
          !pn_p1(:,k) = pn_p1(:,k) + phi(:,k)  !- 0.25*div(:,k,i)/Re - inputs%div_stab_in_ns*div(:,k,i)
          !LC 2017/01/27
          pn_p1(:,k) = pn_p1(:,k) + phi(:,k) - 2*div(:,k,i)/Re - inputs%div_stab_in_ns/Re*div(:,k,i)
       END DO
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Pressure computation

       !===UPDATES
       tps = user_time()
       !===Handling of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, pn_p1(:,1),moyenne)
          pn_p1(:,1) = pn_p1(:,1)-moyenne
       ENDIF
       !===End of handling of mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2) = 0.d0
          un_p1 (:,4) = 0.d0
          un_p1 (:,6) = 0.d0
          pn_p1 (:,2) = 0.d0
       END IF
       !===Correction of zero mode

       IF (inputs%LES) THEN
          un_m2(:,:,i) = un_m1(:,:,i)
       END IF
       !===UPDATES
       un_m1(:,:,i)    = un (:,:,i)
       un   (:,:,i)    = un_p1
       pn_m1(:,:,i)    = pn(:,:,i)
       pn   (:,:,i)    = pn_p1
       incpn_m1(:,:,i) = incpn(:,:,i)
       incpn   (:,:,i) = phi
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End UPDATES

    ENDDO

    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'H1',  vv_mesh, list_mode, un)
       talk_to_me%div_L2  = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)/norm
       talk_to_me%weak_div_L2  = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, div)/norm
    END IF

    !===Compute entropy viscosity
    IF (inputs%LES) THEN
       CALL compute_entropy_viscosity(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, vvz_per, &
            un, un_m1, un_m2, pn_m1, rotv_v, visco_entro_sym_grad_u, density, tempn, concn)
    END IF
    !===End Compute entropy viscosity

    tps_tot = user_time() - tps_tot

  END SUBROUTINE BDF2_ns_stress_bc_with_u
  !============================================

  ! cas PRECESSION 28/07/09
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

  SUBROUTINE smb_CurlH_cross_B_gauss_sft_par(communicator,mesh,list_mode,V_in,W_in,V_out)
    !=================================
    USE sft_parallele
    USE chaine_caractere
    USE input_data
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in, W_in
    REAL(KIND=8), DIMENSION(:,:,:) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: RotV, V, W
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs, Ws
    REAL(KIND=8)   :: ray
    REAL(KIND=8)   :: tps
    REAL(KIND=8), DIMENSION(3)                  :: temps
    !===FOR FFT_PAR_CROSS_PROD_DCL
    INTEGER       :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm       :: communicator

    tps = user_time()
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
             Ws(:,k) = W_in(j_loc,k,i)
          END DO

          DO l = 1, mesh%gauss%l_G
             index = index + 1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Ws(:,1)*mesh%gauss%ww(:,l))
             W(index,3,i) = SUM(Ws(:,3)*mesh%gauss%ww(:,l))
             W(index,5,i) = SUM(Ws(:,5)*mesh%gauss%ww(:,l))

             W(index,2,i) = SUM(Ws(:,2)*mesh%gauss%ww(:,l))
             W(index,4,i) = SUM(Ws(:,4)*mesh%gauss%ww(:,l))
             W(index,6,i) = SUM(Ws(:,6)*mesh%gauss%ww(:,l))
             V(index,1,i) = SUM(Vs(:,1)*mesh%gauss%ww(:,l))
             V(index,3,i) = SUM(Vs(:,3)*mesh%gauss%ww(:,l))
             V(index,5,i) = SUM(Vs(:,5)*mesh%gauss%ww(:,l))

             V(index,2,i) = SUM(Vs(:,2)*mesh%gauss%ww(:,l))
             V(index,4,i) = SUM(Vs(:,4)*mesh%gauss%ww(:,l))
             V(index,6,i) = SUM(Vs(:,6)*mesh%gauss%ww(:,l))
             !-----------------rotational sur les points de Gauss---------------------------
             !coeff sur les cosinus
             RotV(index,1,i) = mode/ray*V(index,6,i) &
                  -SUM(Vs(:,3)*dw_loc(2,:))
             RotV(index,4,i) =          SUM(Vs(:,2)*dw_loc(2,:)) &
                  -SUM(Vs(:,6)*dw_loc(1,:))
             RotV(index,5,i) =    1/ray*V(index,3,i) &
                  +SUM(Vs(:,3)*dw_loc(1,:)) &
                  -mode/ray*V(index,2,i)
             !coeff sur les sinus
             RotV(index,2,i) =-mode/ray*V(index,5,i) &
                  -SUM(Vs(:,4)*dw_loc(2,:))
             RotV(index,3,i) =         SUM(Vs(:,1)*dw_loc(2,:)) &
                  -SUM(Vs(:,5)*dw_loc(1,:))
             RotV(index,6,i) =    1/ray*V(index,4,i) &
                  +SUM(Vs(:,4)*dw_loc(1,:))&
                  +mode/ray*V(index,1,i)
          ENDDO
       ENDDO
    END DO

    !tps = user_time() - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    !tps = user_time()
    temps = 0


    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(RotV,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    CALL FFT_PAR_CROSS_PROD_DCL(communicator, RotV, W, V_out, nb_procs, bloc_size, m_max_pad, temps)
    tps = user_time() - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)
    !DEALLOCATE(Om, W, RotV)

  END SUBROUTINE smb_CurlH_cross_B_gauss_sft_par

  SUBROUTINE smb_kelvin_force_gauss_fft_par(communicator,mesh,list_mode,scal_in,vect_in,vect_out)
    USE Gauss_points
    USE sft_parallele
    USE def_type_mesh
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: scal_in, vect_in
    REAL(KIND=8), DIMENSION(:,:,:)              :: vect_out
    REAL(KIND=8), DIMENSION(mesh%np,2,SIZE(list_mode)) :: chi_coeff, vect_in_square
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: chi_coeff_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: grad_vect_in_square
    INTEGER                                                           :: i, mode, index, m, k, l
    INTEGER,      DIMENSION(:,:), POINTER                             :: jj
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc
    REAL(KIND=8)                                                      :: rad
    !===FOR FFT
    INTEGER                                     :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm                                    :: communicator

    ! We compute the Kelvin force: Kelvin force = chi_coeff(T) * grad(H**2/2)

    !===nb_procs and m_max_pad calculus for FFT
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===chi_coeff(T) on nodes computation
    chi_coeff = scal_in
    bloc_size = SIZE(chi_coeff)/nb_procs+1
    CALL FFT_PAR_SCAL_FUNCT(communicator, chi_coeff, chi_coeff_law, nb_procs, bloc_size, m_max_pad)

    !===chi_coeff(T) on Gauss nodes computation
    CALL gauss(mesh)
    jj => mesh%jj
    DO i = 1, SIZE(list_mode)
       index = 0 ! global index of Gauss node
       DO m = 1, mesh%dom_me
          DO l=1, l_G
             index = index + 1
             DO k = 1, 2
                chi_coeff_gauss(index,k,i) = SUM(chi_coeff(jj(:,m),k,i)*ww(:,l))
             END DO
          END DO
       END DO
    END DO

    !===H**2 on nodes computation
    bloc_size = SIZE(vect_in,1)/nb_procs+1
    CALL FFT_PAR_DOT_PROD_DCL(communicator, vect_in, vect_in, vect_in_square, nb_procs, bloc_size, m_max_pad)

    !===grad(H**2) on Gauss nodes computation
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          DO l=1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)
             rad = SUM(mesh%rr(1,jj(:,m))*ww(:,l)) ! radius of Gauss node
             grad_vect_in_square(index,1,i) = SUM(vect_in_square(jj(:,m),1,i)*dw_loc(1,:))
             grad_vect_in_square(index,2,i) = SUM(vect_in_square(jj(:,m),2,i)*dw_loc(1,:))
             grad_vect_in_square(index,3,i) = mode/rad * SUM(vect_in_square(jj(:,m),2,i)*ww(:,l))
             grad_vect_in_square(index,4,i) = - mode/rad * SUM(vect_in_square(jj(:,m),1,i)*ww(:,l))
             grad_vect_in_square(index,5,i) = SUM(vect_in_square(jj(:,m),1,i)*dw_loc(2,:))
             grad_vect_in_square(index,6,i) = SUM(vect_in_square(jj(:,m),2,i)*dw_loc(2,:))
          END DO
       END DO
    END DO

    !===Kelvin force = chi_coeff(T) * grad(H**2/2) on Gauss nodes computation
    bloc_size = SIZE(chi_coeff_gauss,1)/nb_procs+1
    CALL FFT_SCALAR_VECT_DCL(communicator, 0.5*grad_vect_in_square, chi_coeff_gauss, vect_out, 1, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_kelvin_force_gauss_fft_par

  SUBROUTINE smb_helmholtz_force_gauss_fft_par(communicator,mesh,list_mode,scal_in,vect_in,vect_out)
    USE Gauss_points
    USE sft_parallele
    USE def_type_mesh
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)  :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: scal_in, vect_in
    REAL(KIND=8), DIMENSION(:,:,:)              :: vect_out
    REAL(KIND=8), DIMENSION(mesh%np,2,SIZE(list_mode)) :: vect_in_square, chi_coeff
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: vect_in_square_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: grad_chi_coeff
    INTEGER                                                           :: i, mode, index, m, k, l
    INTEGER,      DIMENSION(:,:), POINTER                             :: jj
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc
    REAL(KIND=8)                                                      :: rad
    !===FOR FFT
    INTEGER                                     :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm                                    :: communicator

    ! We compute the Helmholtz force: Helmholtz force = - H**2/2 * grad(chi_coeff(T))

    !===nb_procs and m_max_pad calculus for FFT
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===H**2 on nodes computation
    bloc_size = SIZE(vect_in,1)/nb_procs+1
    CALL FFT_PAR_DOT_PROD_DCL(communicator, vect_in, vect_in, vect_in_square, nb_procs, bloc_size, m_max_pad)

    !===H**2 on Gauss nodes computation
    CALL gauss(mesh)
    jj => mesh%jj
    DO i = 1, SIZE(list_mode)
       index = 0 ! global index of Gauss node
       DO m = 1, mesh%dom_me
          DO l=1, l_G
             index = index + 1
             DO k = 1, 2
                vect_in_square_gauss(index,k,i) = SUM(vect_in_square(jj(:,m),k,i)*ww(:,l))
             END DO
          END DO
       END DO
    END DO

    !===chi_coeff(T) on nodes computation
    chi_coeff = scal_in
    bloc_size = SIZE(chi_coeff)/nb_procs+1
    CALL FFT_PAR_SCAL_FUNCT(communicator, chi_coeff, chi_coeff_law, nb_procs, bloc_size, m_max_pad)

    !===grad(chi_coeff(T)) on Gauss nodes computation
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          DO l=1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)
             rad = SUM(mesh%rr(1,jj(:,m))*ww(:,l)) ! radius of Gauss node
             grad_chi_coeff(index,1,i) = SUM(chi_coeff(jj(:,m),1,i)*dw_loc(1,:))
             grad_chi_coeff(index,2,i) = SUM(chi_coeff(jj(:,m),2,i)*dw_loc(1,:))
             grad_chi_coeff(index,3,i) = mode/rad * SUM(chi_coeff(jj(:,m),2,i)*ww(:,l))
             grad_chi_coeff(index,4,i) = - mode/rad * SUM(chi_coeff(jj(:,m),1,i)*ww(:,l))
             grad_chi_coeff(index,5,i) = SUM(chi_coeff(jj(:,m),1,i)*dw_loc(2,:))
             grad_chi_coeff(index,6,i) = SUM(chi_coeff(jj(:,m),2,i)*dw_loc(2,:))
          END DO
       END DO
    END DO

    !===Helmholtz force = - H**2/2 * grad(chi_coeff(T)) on Gauss nodes computation
    bloc_size = SIZE(grad_chi_coeff,1)/nb_procs+1
    CALL FFT_SCALAR_VECT_DCL(communicator, grad_chi_coeff, -0.5*vect_in_square_gauss, vect_out, 1, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_helmholtz_force_gauss_fft_par

  SUBROUTINE smb_nu_tilde_fft_par(communicator,list_mode,scal_inout) ! MODIFICATION: computation of nu tilde function of temperature
    USE Gauss_points
    USE sft_parallele
    USE def_type_mesh
    USE boundary
    IMPLICIT NONE
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: scal_inout
    !===FOR FFT
    INTEGER                                     :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm                                    :: communicator

    !===nb_procs and m_max_pad calculus for FFT
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    bloc_size = SIZE(scal_inout)/nb_procs+1
    CALL FFT_PAR_SCAL_FUNCT(communicator, scal_inout, nu_tilde_law, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_nu_tilde_fft_par

  SUBROUTINE Moy(communicator,mesh,p,RESLT)
    !===========================
    !moyenne
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

END MODULE subroutine_ns_with_u
