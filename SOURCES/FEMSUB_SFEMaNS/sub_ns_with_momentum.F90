!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_ns_with_m
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC :: three_level_ns_tensor_sym_with_m, smb_explicit_diffu_sym, smb_surface_tension, &
       smb_explicit_LES, momentum_dirichlet, smb_buoyancy
  PRIVATE
CONTAINS

  SUBROUTINE three_level_ns_tensor_sym_with_m(comm_one_d, time, vv_3_LA, pp_1_LA, &
       dt, Re, list_mode, pp_mesh, vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un,  &
       Hn_p2, Bn_p2, density_m1, density, density_p1, visco_dyn, tempn, concn, level_set, level_set_p1, &
       visc_entro_level, level_set_reg, visc_entro_grad_mom)

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
    USE sft_parallele
    USE tn_axi
    USE verbose
    USE entropy_viscosity
    USE subroutine_mass
    IMPLICIT NONE
    REAL(KIND=8)                                          :: time, dt, Re
    INTEGER,      DIMENSION(:),       INTENT(IN)          :: list_mode
    TYPE(mesh_type),                  INTENT(IN)          :: pp_mesh, vv_mesh
    TYPE(petsc_csr_LA)                                    :: vv_3_LA, pp_1_LA
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(INOUT)       :: incpn_m1, incpn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(INOUT)       :: pn_m1, pn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(INOUT)       :: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: density_m1, density, density_p1
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)          :: level_set, level_set_p1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: visco_dyn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: concn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: Hn_p2, Bn_p2
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)) :: Hn_p2_aux
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)         :: visc_entro_level
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)          :: level_set_reg
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(INOUT)       :: visc_entro_grad_mom
!TEST LC LES_SUITE 2024/06
    INTEGER,                                         SAVE :: m_max_c
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE,   SAVE :: pp_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,       SAVE :: pp_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3),                SAVE :: vv_js_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,       SAVE :: vv_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE,   SAVE :: vel_global_D
    LOGICAL,                                         SAVE :: once = .TRUE.
!TEST LC LES_SUITE 2024/06
    LOGICAL,                                         SAVE :: once_LES = .TRUE.
!TEST LC LES_SUITE 2024/06
    INTEGER,                                         SAVE :: my_petscworld_rank
    REAL(KIND=8),                                    SAVE :: mu_bar, nu_bar, rho_bar, sqrt_2d_vol
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE,     SAVE :: momentum, momentum_m1, momentum_m2
    INTEGER,                                         SAVE :: bloc_size, m_max_pad, nb_procs
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:),     SAVE :: rotb_b_m1
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:),   SAVE :: visc_grad_vel_m1
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:),   SAVE :: tensor_m1
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:),       SAVE :: visc_entro_real
    !----------FIN SAVE--------------------------------------------------------------------

    !----------Declaration sans save-------------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: pp_1_ifrom, vv_3_ifrom
    INTEGER                                  :: i, k, m, n
    INTEGER                                  :: code,nu_mat, mode
    REAL(KIND=8)                             :: moyenne
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: div
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2)   :: pn_p1, phi
    REAL(KIND=8), DIMENSION(vv_mesh%np, 2)   :: p_p2
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)   :: un_p1, src
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: rotb_b, rotb_b_aux
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: visc_grad_vel
!TEST LC LES_SUITE 2024/06
!    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: visc_entro_grad_mom
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: tensor_surface_gauss
    REAL(KIND=8), DIMENSION(3,vv_mesh%np,6,SIZE(list_mode))                       :: tensor
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6)                 :: visc_grad_vel_ext
    REAL(KIND=8), DIMENSION(3,vv_mesh%np,6)                                       :: tensor_ext
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                         :: uext, momentumext, momentum_exact
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                         :: buoyancy
    REAL(KIND=8), DIMENSION(inputs%nb_fluid-1, vv_mesh%np, 2, SIZE(list_mode))    :: level_set_FEM_P2
    REAL(KIND=8), DIMENSION(vv_mesh%np)                      :: vel_loc, vel_tot
    REAL(KIND=8), DIMENSION(SIZE(list_mode))                 :: normalization_mt
    REAL(KINd=8)   :: vel_tot_max_S,vel_tot_max
    INTEGER        :: n1, n2, n3, n123, nb_inter
    REAL(KIND=8)   :: tps, tps_tot, tps_cumul, coeff, cfl, cfl_max, norm
    INTEGER        :: nb_procs_LES, bloc_size_LES, m_max_pad_LES
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

       !-----CREATE PETSC VECTORS AND GHOSTS-----------------------------------------
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
       !------------------------------------------------------------------------------

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------

       !-------------MOMENTUM INITIALIZATION------------------------------------------
       ALLOCATE(momentum(SIZE(un,1),SIZE(un,2), SIZE(un,3)),&
            momentum_m1(SIZE(un,1),SIZE(un,2), SIZE(un,3)), &
            momentum_m2(SIZE(un,1),SIZE(un,2), SIZE(un,3)))
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       bloc_size = SIZE(un,1)/nb_procs+1
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2

       IF (inputs%if_level_set) THEN
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), un_m1, density_m1, momentum_m1, 1, nb_procs, &
               bloc_size, m_max_pad)
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), un, density, momentum, 1, nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentum_m1 = un_m1
          momentum = un
       END IF

       ALLOCATE(rotb_b_m1(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)))
       rotb_b_m1 = 0.d0

       ALLOCATE(tensor_m1(3,vv_mesh%np,6,SIZE(list_mode)))
       bloc_size = SIZE(un,1)/nb_procs+1
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       CALL FFT_TENSOR_DCL(comm_one_d(2), momentum_m1, un_m1, tensor_m1, nb_procs, bloc_size, m_max_pad)

       ALLOCATE(visc_grad_vel_m1(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)))
       CALL smb_explicit_diffu_sym(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            visco_dyn/Re, un_m1, visc_grad_vel_m1)

       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs_LES, code)
       bloc_size_LES = vv_mesh%gauss%l_G*vv_mesh%dom_me/nb_procs_LES+1
       bloc_size_LES = vv_mesh%gauss%l_G*(bloc_size_LES/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
       m_max_pad_LES = 3*SIZE(list_mode)*nb_procs_LES/2
       ALLOCATE(visc_entro_real(2*m_max_pad_LES-1,bloc_size_LES))
       visc_entro_real = 0.d0

       !---------PREPARE pp_js_D ARRAY FOR PRESSURE-----------------------------------
       !CALL dirichlet_nodes_parallel(pp_mesh, inputs%pp_list_dirichlet_sides, pp_js_D)
       !------------------------------------------------------------------------------
       !===PREPARE pp_mode_global_js_D ARRAY FOR PRESSURE
       !===ATTENTION pressure BCs are no longer implemented
       !===JLG June 9 2017
       !CALL scalar_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
       CALL scalar_without_glob_js_D(pp_mesh, list_mode, pp_1_LA, pp_mode_global_js_D)
       ALLOCATE(pp_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(pp_global_D(i)%DRL(SIZE(pp_mode_global_js_D(i)%DIL)))
       END DO

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
       DO i = 1, m_max_c
          IF (list_mode(i)==0) CYCLE
          CALL Dirichlet_M_parallel(mass_mat,pp_mode_global_js_D(i)%DIL)
       END DO
       CALL init_solver(inputs%my_par_mass,mass_ksp,mass_mat,comm_one_d(1),&
            solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)

       IF (MINVAL(list_mode)==0) THEN
          CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, mass_mat0, CLEAN=.FALSE.)
          CALL qs_diff_mass_scal_M (pp_mesh, pp_1_LA, 0.d0, 1.d0, 0.d0, 0, mass_mat0)
          DO i = 1, m_max_c
             IF (list_mode(i).NE.0) CYCLE
             CALL Dirichlet_M_parallel(mass_mat0,pp_mode_global_js_D(i)%DIL)
             CALL init_solver(inputs%my_par_mass,mass_ksp0,mass_mat0,comm_one_d(1),&
                  solver=inputs%my_par_mass%solver,precond=inputs%my_par_mass%precond)
          END DO
       END IF
       !===END ASSEMBLE MASS MATRIX

       !-------------ASSEMBLE VELOCITY MATRICES----------------------
       ALLOCATE(vel_mat(2*m_max_c),vel_ksp(2*m_max_c))
       ALLOCATE(press_mat(m_max_c),press_ksp(m_max_c))

       ! Definition of nu_bar
       IF (inputs%if_level_set) THEN
          !To be tested thoroughly
          mu_bar  = MINVAL(inputs%dyna_visc_fluid)
          rho_bar = MINVAL(inputs%density_fluid)
          nu_bar = 0.d0
          DO n = 1, inputs%nb_fluid
             nu_bar = MAX(nu_bar,inputs%dyna_visc_fluid(n)/inputs%density_fluid(n))
          END DO
          !nu_bar = 2.0d0*nu_bar
          nu_bar = 1.1d0*nu_bar
       ELSE
          mu_bar  = 1.d0
          nu_bar  = 1.d0
          rho_bar = 1.d0
       END IF

       DO i = 1, m_max_c
          mode = list_mode(i)
          !---VELOCITY
          nu_mat = 2*i-1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          IF (inputs%if_moment_bdf2) THEN
             CALL qs_diff_mass_vect_3x3 (1, vv_3_LA, vv_mesh, nu_bar/Re, three/(2*dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          ELSE
             CALL qs_diff_mass_vect_3x3 (1, vv_3_LA, vv_mesh, nu_bar/Re, one/dt, &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)

          nu_mat = nu_mat+1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          IF (inputs%if_moment_bdf2) THEN
             CALL qs_diff_mass_vect_3x3 (2, vv_3_LA, vv_mesh, nu_bar/Re, three/(2*dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          ELSE
             CALL qs_diff_mass_vect_3x3 (2, vv_3_LA, vv_mesh, nu_bar/Re, one/dt, &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, i, mode, vel_mat(nu_mat))
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)

          !---PRESSURE
          CALL create_local_petsc_matrix(comm_one_d(1), pp_1_LA, press_mat(i), clean=.FALSE.)
          CALL qs_diff_mass_scal_M(pp_mesh, pp_1_LA, one, 1.d-10, zero, mode, press_mat(i))
          CALL Dirichlet_M_parallel(press_mat(i),pp_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_pp,press_ksp(i),press_mat(i),comm_one_d(1),&
               solver=inputs%my_par_pp%solver,precond=inputs%my_par_pp%precond)

       ENDDO

       CALL twoD_volume(comm_one_d(1),vv_mesh,sqrt_2d_vol)
       sqrt_2d_vol =  SQRT(sqrt_2d_vol)

    ENDIF ! Fin du once
    tps_tot = user_time()
    tps_cumul = 0


    !===Compute rhs by FFT at Gauss points
    tps = user_time()

    IF (inputs%if_moment_bdf2) THEN
       uext = 2*un-un_m1
       IF (inputs%if_level_set) THEN
          ! BDF2: momentumext = rho_np1*(2*un-un_m1)
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), uext, density_p1, momentumext, 1,nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentumext = 2*momentum - momentum_m1
       END IF
    ELSE
       uext = un
       IF (inputs%if_level_set) THEN
          ! BDF1: momentumext = rho_np1*un
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), un, density_p1, momentumext, 1, nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentumext = momentum
       END IF
    END IF

    !===Precession should be in condlim
    IF (inputs%precession) THEN
       CALL error_petsc('for momentum ns : precession should be in condlim')
    END IF

    !===Compute Lorentz force if mhd
    IF (inputs%type_pb=='mhd') THEN
       !===Compute Lorentz force if mhd in quasi-static limit
       IF (inputs%if_quasi_static_approx) THEN
          DO i = 1, m_max_c
             mode = list_mode(i)
             Hn_p2_aux(:,:,i) = H_B_quasi_static('H', vv_mesh%rr, mode)
          END DO
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,Hn_p2_aux,Bn_p2,rotb_b_aux)
          DO i = 1, m_max_c
             mode = list_mode(i)
             Hn_p2_aux(:,:,i)  = H_B_quasi_static('B', vv_mesh%rr, mode)
          END DO
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,Hn_p2,Hn_p2_aux,rotb_b)
          rotb_b = rotb_b + rotb_b_aux
       ELSE !===Compute Lorentz force if mhd
          CALL smb_CurlH_cross_B_gauss_sft_par(comm_one_d(2),vv_mesh,list_mode,Hn_p2,Bn_p2,rotb_b)
       END IF
    ELSE
       rotb_b = 0.d0
    END IF
    !===End compute Lorentz force if mhd

    !===Compute diffusion/corrections and surface tension
    IF (inputs%if_level_set) THEN
       !===Compute visco_dyn/Re*Grad(un)
       CALL smb_explicit_diffu_sym(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            visco_dyn/Re, un, visc_grad_vel)
       !===End compute visco_dyn/Re*Grad(un)

       IF (inputs%if_surface_tension) THEN
          !===Compute coeff_surface*Grad(level_set):Grad(level_set)
          IF (inputs%if_level_set_P2) THEN
             CALL smb_surface_tension(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
                  level_set_reg, tensor_surface_gauss)
          ELSE
             DO nb_inter = 1, inputs%nb_fluid-1
                DO i = 1, SIZE(list_mode)
                   DO k = 1, 2
                      CALL inject_P1_P2(pp_mesh%jj, vv_mesh%jj, level_set_reg(nb_inter,:,k,i), &
                           level_set_FEM_P2(nb_inter,:,k,i))
                   END DO
                END DO
             END DO
             CALL smb_surface_tension(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
                  level_set_FEM_P2, tensor_surface_gauss)
          END IF
       ELSE
          tensor_surface_gauss = 0.d0
       END IF

       !===Compute buoyancy force heat_grav*density*T e_z
       IF (inputs%if_temperature) THEN
          CALL smb_buoyancy(comm_one_d, vv_mesh, pp_mesh, list_mode, level_set_p1, tempn, buoyancy)
       ELSE
          buoyancy = 0.d0
       END IF
    ELSE
       visc_grad_vel        = 0.d0
       tensor_surface_gauss = 0.d0
       buoyancy = 0.d0 !Defined in condlim.F90
    END IF
    !===End Compute diffusion/corrections and surface tension

!TEST LC LES_SUITE 2024/06
    !===First iteration without restart: compute -LES_coeff1*Grad(momentumext)
    IF (once_LES) THEN
       once_LES = .FALSE.
       IF (.NOT.inputs%irestart_LES) THEN
          !visc_entro_grad_mom initialized to zero in initialization.F90
          !visc_entro_real is set to zero in once
          IF (inputs%LES.AND.inputs%if_LES_in_momentum) THEN
             CALL smb_explicit_LES(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
               visc_entro_real, momentumext, visc_entro_grad_mom)
          END IF
       END IF
    END IF
    !===End First iteration without restart: compute -LES_coeff1*Grad(momentumext)
!TEST LC LES_SUITE 2024/06

    !===Compute tensor product of momentum by velocity
    CALL FFT_TENSOR_DCL(comm_one_d(2), momentum, un, tensor, nb_procs, bloc_size, m_max_pad)

    !===PREPARE BOUNDARY CONDITION FOR MOMENTUM
    CALL momentum_dirichlet(comm_one_d(2), vv_mesh, list_mode, time, nb_procs, density_p1, &
         momentum_exact, vv_js_D)

    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Tps fft vitesse', tps

    !------------BEGINING LOOP ON FOURIER MODES------------------
    DO i = 1, m_max_c
       mode = list_mode(i)
       !===Compute phi
       tps = user_time()
       !jan 29 2007
       pn_p1(:,:) = pn(:,:,i)
       IF (inputs%if_moment_bdf2) THEN
          phi = pn_p1(:,:) + (4.d0 * incpn(:,:,i) - incpn_m1(:,:,i))/3.d0
       ELSE
          phi = pn_p1(:,:) + incpn(:,:,i)
       END IF
       !jan 29 2007

       !===Inject pressure P1 -> P2
       DO k = 1, 2
          CALL inject(pp_mesh%jj, vv_mesh%jj, phi(:,k), p_p2(:,k))
       ENDDO

       !===Prediction step
       DO k = 1, 6
          src(:,k) = source_in_NS_momentum(k, vv_mesh%rr, mode, i, time, Re, 'ns', &
               density_p1, tempn, concn) + buoyancy(:,k,i)
       END DO

       IF (inputs%if_moment_bdf2) THEN
          tensor_ext               = 2*tensor(:,:,:,i) - tensor_m1(:,:,:,i)
          visc_grad_vel_ext        = 2*visc_grad_vel(:,:,:,i) - visc_grad_vel_m1(:,:,:,i)

          !===Update terms at m1
          rotb_b_m1(:,:,i) = rotb_b(:,:,i)
          tensor_m1(:,:,:,i) = tensor(:,:,:,i)
          visc_grad_vel_m1(:,:,:,i) = visc_grad_vel(:,:,:,i)
          CALL qs_ns_momentum_stab_3x3(vv_mesh, vv_3_LA, mode, src, &
               (4*momentum(:,:,i)-momentum_m1(:,:,i))/(2*dt), p_p2(:,:), &
               vb_3_145, vb_3_236, rotb_b(:,:,i), tensor_ext,&
               tensor_surface_gauss(:,:,:,i), nu_bar/Re, momentumext(:,:,i), &
               -visc_grad_vel_ext-visc_entro_grad_mom(:,:,:,i))
       ELSE
          CALL qs_ns_momentum_stab_3x3(vv_mesh, vv_3_LA, mode, src, &
               momentum(:,:,i)/dt, p_p2(:,:), &
               vb_3_145, vb_3_236, rotb_b(:,:,i), tensor(:,:,:,i),&
               tensor_surface_gauss(:,:,:,i), nu_bar/Re, momentumext(:,:,i), &
               -visc_grad_vel(:,:,:,i)-visc_entro_grad_mom(:,:,:,i))
       END IF

       !===Axis boundary conditions
       n1 = SIZE(vv_js_D(1)%DIL)
       n2 = SIZE(vv_js_D(2)%DIL)
       n3 = SIZE(vv_js_D(3)%DIL)
       n123 = n1+n2+n3
       vel_global_D(i)%DRL(1:n1)        = momentum_exact(vv_js_D(1)%DIL,1,i)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = momentum_exact(vv_js_D(2)%DIL,4,i)
       vel_global_D(i)%DRL(n1+n2+1:n123)= momentum_exact(vv_js_D(3)%DIL,5,i)
       vel_global_D(i)%DRL(n123+1:)     = 0.D0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        = momentum_exact(vv_js_D(1)%DIL,2,i)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = momentum_exact(vv_js_D(2)%DIL,3,i)
       vel_global_D(i)%DRL(n1+n2+1:n123)= momentum_exact(vv_js_D(3)%DIL,6,i)
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss


       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps second membre vitesse', tps
       !-------------------------------------------------------------------------------------

       !--------------------INVERSION DE L'OPERATEUR 1 --------------
       tps = user_time()
       !Solve system 1, ur_c, ut_s, uz_c
       nu_mat  =2*i-1
       CALL solver(vel_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,1))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,4))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,5))

       !Solve system 2, ur_s, ut_c, uz_s
       nu_mat = nu_mat + 1
       CALL solver(vel_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,un_p1(:,2))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,un_p1(:,3))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,un_p1(:,6))

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de vitesse', tps, 'for mode ', mode
       !-------------------------------------------------------------------------------------

       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010
       IF (mode==0) THEN
          un_p1 (:,2) = 0.d0
          un_p1 (:,4) = 0.d0
          un_p1 (:,6) = 0.d0
          pn_p1 (:,2) = 0.d0
       END IF
       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010

       momentum_m2(:,:,i) = momentum_m1 (:,:,i)
       momentum_m1(:,:,i) = momentum (:,:,i)
       momentum (:,:,i) = un_p1

    END DO

    !---------------UPDATE VELOCITY---------------------
    un_m1 = un
    IF (inputs%if_level_set) THEN
       CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), momentum, density_p1, un, 2, nb_procs, &
            bloc_size, m_max_pad)
    ELSE
       un = momentum
    END IF


    DO i = 1, m_max_c
       mode = list_mode(i)
       !===Compute phi
       tps = user_time()
       !jan 29 2007
       pn_p1(:,:) = pn(:,:,i)
       !jan 29 2007

       !---------------SOLUTION OF THE POISSON EQUATION--------------
       ! BDF2 : solve -LAP(PH3) = -3*rho_bar/(2*dt)*DIV(un_p1)
       ! BDF1 : solve -LAP(PH3) = -rho_bar/dt*DIV(un_p1)
       tps = user_time()
       CALL qs_01_div_hybrid_2006(vv_mesh, pp_mesh, pp_1_LA, mode, un(:,:,i), pb_1, pb_2)
       !pb_1, and pb_2 are petsc vectors for the rhs divergence

       !===ATENTION BCs are no longer implemented for pressure
       !===Boundary condition on axis for pressure
       pp_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_1)
       CALL dirichlet_rhs(pp_mode_global_js_D(i)%DIL-1,pp_global_D(i)%DRL,pb_2)
       !===End boundary condition on axis for pressure

       CALL solver(press_ksp(i),pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_pp%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,phi(:,1))

       CALL solver(press_ksp(i),pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_pp%verbose)
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,phi(:,2))

       !Don't forget the -(1.5d0/dt)*min(rh0) coefficient
       IF (inputs%if_moment_bdf2) THEN
          phi = -phi*(1.5d0/dt)*rho_bar
       ELSE
          phi = -phi/dt*rho_bar
       END IF
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de pression', tps, 'for mode ', mode
       !-------------------------------------------------------------------------------------


       !---------------CORRECTION DE LA PRESSION-----------------------
       tps = user_time()
       !CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       IF (mode==0) THEN
          CALL solver(mass_ksp0,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       ELSE
          CALL solver(mass_ksp,pb_1,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       END IF
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,1,i))

       !CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       IF (mode==0) THEN
          CALL solver(mass_ksp0,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       ELSE
          CALL solver(mass_ksp,pb_2,px_1,reinit=.FALSE.,verbose=inputs%my_par_mass%verbose)
       END IF
       CALL VecGhostUpdateBegin(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(px_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(px_1_ghost,1,1,pp_1_LA,div(:,2,i))

       !Pressure computation
       !jan 29 2007
       DO k=1, 2
          pn_p1(:,k) = pn_p1(:,k) + phi(:,k) - div(:,k,i)*(mu_bar/Re)
          !pn_p1(:,k) = pn_p1(:,k) + phi(:,k) !- div(:,k,i)*(mu_bar/Re)
       END DO
       !jan 29 2007
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps correction de divergence', tps
       !-------------------------------------------------------------------------------------

       !---------------UPDATE PRESSURE---------------------
       tps = user_time()
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, pn_p1(:,1),moyenne)
          pn_p1(:,1) = pn_p1(:,1)-moyenne
       ENDIF

       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010
       IF (mode==0) THEN
          pn_p1 (:,2) = 0.d0
       END IF
       !JLG AR, Dec 18 2008/JLG Bug corrige Jan 23 2010

       pn_m1(:,:,i)    = pn(:,:,i)
       pn   (:,:,i)    = pn_p1

       incpn_m1(:,:,i) = incpn(:,:,i)
       incpn   (:,:,i) = phi

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps  des updates', tps
       !-------------------------------------------------------------------------------------

    ENDDO
    tps_tot = user_time() - tps_tot

    !===Verbose divergence of velocity
    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'H1',  vv_mesh, list_mode, un)
       talk_to_me%div_L2  = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)/norm
       talk_to_me%weak_div_L2  = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, div)/norm
    END IF
    !===End verbose divergence of velocity

    !===Compute vel max and CFL
    IF (inputs%verbose_CFL.OR.inputs%LES) THEN
       vel_loc = 0.d0
       DO i = 1, m_max_c
          IF (list_mode(i)==0) THEN
             coeff = 1.d0
          ELSE
             coeff = .5d0
          END IF
          vel_loc = vel_loc + coeff*(un(:,1,i)**2+un(:,2,i)**2+un(:,3,i)**2+&
               un(:,4,i)**2+un(:,5,i)**2+un(:,6,i)**2)
       END DO
       CALL MPI_COMM_SIZE(comm_one_d(2),nb_procs,code)
       CALL MPI_ALLREDUCE(vel_loc,vel_tot,vv_mesh%np,MPI_DOUBLE_PRECISION, MPI_SUM, comm_one_d(2), code)

       vel_tot = sqrt(vel_tot)
       vel_tot_max_S = MAXVAL(vel_tot)
       CALL MPI_ALLREDUCE(vel_tot_max_S,vel_tot_max,1,MPI_DOUBLE_PRECISION, MPI_MAX, comm_one_d(1), code)

       !===New normalization (JLG; April 24, 2015)
       DO i = 1, m_max_c
          normalization_mt(i) = norm_S(comm_one_d, 'L2', vv_mesh, list_mode, momentum)/(sqrt(2.d0)*sqrt_2d_vol)
       END DO
       !===New normalization (JLG; April 24, 2015)

       !===Computation of CFL
       IF (inputs%verbose_CFL) THEN
          cfl = 0
          DO m = 1, vv_mesh%dom_me
             cfl = MAX(vel_tot_max*dt/vv_mesh%hloc(m),cfl)
          END DO
          CALL MPI_ALLREDUCE(cfl,cfl_max,1,MPI_DOUBLE_PRECISION, MPI_MAX, comm_one_d(1), code)
          talk_to_me%CFL=cfl_max
          talk_to_me%time=time
          !IF (my_petscworld_rank==0) THEN
          !   WRITE(*,'(3(A,e10.3))') ' Time = ', time, ', CFL = ', cfl_max, 'vel_tot_max', vel_tot_max
          !END IF
       END IF
    END IF

    !===Compute entropy viscosity
    IF (inputs%if_level_set) THEN
       IF (inputs%LES) THEN
          CALL compute_entropy_viscosity_mom(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, &
               momentum, momentum_m1, momentum_m2, pn_m1, un_m1, tensor, visc_grad_vel, tensor_surface_gauss, &
               rotb_b, visco_dyn, density_m1, density, density_p1, tempn, concn, visc_entro_real, visc_entro_level)
       ELSE
          visc_entro_real = 0.d0
          visc_entro_level= 0.d0
       END IF
    ELSE
       IF (inputs%if_LES_in_momentum) THEN
          CALL compute_entropy_viscosity_mom_no_level_set(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, &
               momentum, momentum_m1, momentum_m2, pn_m1, un_m1, tensor, rotb_b, density, tempn, concn, visc_entro_real)
       ELSE
          visc_entro_real=0.d0
       END IF
       !visc_entro_level not allocated so nothing to do
    END IF
    !===End compute entropy viscosity

!TEST LC LES_SUITE 2024/06
       !===Compute entropy viscosity stabilization: (-LES_coeff1+visc_entro)*Grad(momentumext)
       IF (inputs%LES.AND.inputs%if_LES_in_momentum) THEN
          CALL smb_explicit_LES(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
               visc_entro_real, momentumext, visc_entro_grad_mom)
       ELSE
          visc_entro_grad_mom=0.d0
       END IF
       !===End compute entropy viscosity stabilization: (-LES_coeff1+visc_entro)*Grad(momentumext)
!TEST LC LES_SUITE 2024/06

    !===Dummies variables to avoid warning
    i=SIZE(level_set,1)
    !===Dummies variables to avoid warning
  END SUBROUTINE three_level_ns_tensor_sym_with_m
  !============================================

  SUBROUTINE inject(jj_c, jj_f, pp_c, pp_f)
    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
    REAL(KIND=8) :: half = 0.5
    INTEGER:: m
    IF (SIZE(jj_c,1)==3) THEN
       DO m = 1, SIZE(jj_f,2)
          pp_f(jj_f(1:3,m)) =  pp_c(jj_c(:,m))
          pp_f(jj_f(4,m)) = (pp_c(jj_c(2,m)) + pp_c(jj_c(3,m)))*half
          pp_f(jj_f(5,m)) = (pp_c(jj_c(3,m)) + pp_c(jj_c(1,m)))*half
          pp_f(jj_f(6,m)) = (pp_c(jj_c(1,m)) + pp_c(jj_c(2,m)))*half
       END DO

    ELSE
       DO m = 1, SIZE(jj_f,2)
          pp_f(jj_f(1:4,m)) =  pp_c(jj_c(:,m))
       END DO
       pp_f(jj_f(5,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(4,:)))*half
       pp_f(jj_f(6,:)) = (pp_c(jj_c(4,:)) + pp_c(jj_c(2,:)))*half
       pp_f(jj_f(7,:)) = (pp_c(jj_c(2,:)) + pp_c(jj_c(3,:)))*half
       pp_f(jj_f(8,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(4,:)))*half
       pp_f(jj_f(9,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(1,:)))*half
       pp_f(jj_f(10,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(2,:)))*half

    END IF

  END SUBROUTINE inject

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

  SUBROUTINE smb_explicit_diffu_sym(communicator, mesh, list_mode, nb_procs, visc_dyn, vel, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: visc_dyn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: grad1_vel, grad2_vel, grad3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: gradT1_vel, gradT2_vel, gradT3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_tensor_1, part_tensor_2
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_visc_sym_grad_1
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_visc_sym_grad_2
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: visc_dyn_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, l , i, mode, index, k
    INTEGER                                     :: m_max_pad, bloc_size
    MPI_Comm          :: communicator

    CALL gauss(mesh)

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vel_loc(:,k) = vel(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------Dynamic viscosity on Gauss points---------------------------
             visc_dyn_gauss(index,1,i) = SUM(visc_dyn(j_loc,1,i)*ww(:,l))
             visc_dyn_gauss(index,2,i) = SUM(visc_dyn(j_loc,2,i)*ww(:,l))

             !-----------------Grad u_r on Gauss points------------------------------------
             grad1_vel(index,1,i) = SUM(vel_loc(:,1)*dw_loc(1,:))
             grad1_vel(index,2,i) = SUM(vel_loc(:,2)*dw_loc(1,:))
             grad1_vel(index,3,i) = (mode*SUM(vel_loc(:,2)*ww(:,l)) - SUM(vel_loc(:,3)*ww(:,l)))/ray
             grad1_vel(index,4,i) = (-mode*SUM(vel_loc(:,1)*ww(:,l)) - SUM(vel_loc(:,4)*ww(:,l)))/ray
             grad1_vel(index,5,i) =  SUM(vel_loc(:,1)*dw_loc(2,:))
             grad1_vel(index,6,i) =  SUM(vel_loc(:,2)*dw_loc(2,:))

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(index,1,i) = SUM(vel_loc(:,3)*dw_loc(1,:))
             grad2_vel(index,2,i) = SUM(vel_loc(:,4)*dw_loc(1,:))
             grad2_vel(index,3,i) = (mode*SUM(vel_loc(:,4)*ww(:,l)) + SUM(vel_loc(:,1)*ww(:,l)))/ray
             grad2_vel(index,4,i) = (-mode*SUM(vel_loc(:,3)*ww(:,l)) + SUM(vel_loc(:,2)*ww(:,l)))/ray
             grad2_vel(index,5,i) = SUM(vel_loc(:,3)*dw_loc(2,:))
             grad2_vel(index,6,i) = SUM(vel_loc(:,4)*dw_loc(2,:))

             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(index,1,i) = SUM(vel_loc(:,5)*dw_loc(1,:))
             grad3_vel(index,2,i) = SUM(vel_loc(:,6)*dw_loc(1,:))
             grad3_vel(index,3,i) = mode*SUM(vel_loc(:,6)*ww(:,l))/ray
             grad3_vel(index,4,i) = -mode*SUM(vel_loc(:,5)*ww(:,l))/ray
             grad3_vel(index,5,i) = SUM(vel_loc(:,5)*dw_loc(2,:))
             grad3_vel(index,6,i) = SUM(vel_loc(:,6)*dw_loc(2,:))
          ENDDO
       ENDDO
    END DO

    IF (inputs%if_tensor_sym) THEN
       !-----------------GradT u_r on Gauss points------------------------------------
       gradT1_vel(:,1,:) = grad1_vel(:,1,:)
       gradT1_vel(:,2,:) = grad1_vel(:,2,:)
       gradT1_vel(:,3,:) = grad2_vel(:,1,:)
       gradT1_vel(:,4,:) = grad2_vel(:,2,:)
       gradT1_vel(:,5,:) = grad3_vel(:,1,:)
       gradT1_vel(:,6,:) = grad3_vel(:,2,:)
       !-----------------GradT u_th on Gauss points-----------------------------------
       gradT2_vel(:,1,:) = grad1_vel(:,3,:)
       gradT2_vel(:,2,:) = grad1_vel(:,4,:)
       gradT2_vel(:,3,:) = grad2_vel(:,3,:)
       gradT2_vel(:,4,:) = grad2_vel(:,4,:)
       gradT2_vel(:,5,:) = grad3_vel(:,3,:)
       gradT2_vel(:,6,:) = grad3_vel(:,4,:)
       !-----------------GradT u_z on Gauss points------------------------------------
       gradT3_vel(:,1,:) = grad1_vel(:,5,:)
       gradT3_vel(:,2,:) = grad1_vel(:,6,:)
       gradT3_vel(:,3,:) = grad2_vel(:,5,:)
       gradT3_vel(:,4,:) = grad2_vel(:,6,:)
       gradT3_vel(:,5,:) = grad3_vel(:,5,:)
       gradT3_vel(:,6,:) = grad3_vel(:,6,:)

       !===Grad = Grad + Grad^T
       grad1_vel = grad1_vel + gradT1_vel
       grad2_vel = grad2_vel + gradT2_vel
       grad3_vel = grad3_vel + gradT3_vel
    END IF

    !===Compute (visc_dyn + nu_entro)*Grad^s(vel)
    part_tensor_1          = grad1_vel(:,:,:)
    part_tensor_2(:,1:4,:) = grad2_vel(:,3:6,:)
    part_tensor_2(:,5:6,:) = grad3_vel(:,5:6,:)

    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(grad1_vel,1)/nb_procs+1
    CALL FFT_COMPUTE_DIFFU_MOM(communicator,part_tensor_1, part_tensor_2, visc_dyn_gauss, &
         part_visc_sym_grad_1, part_visc_sym_grad_2, nb_procs, bloc_size, m_max_pad)

    !===Using strain rate tensor symmetry to get V1_out= visc_dyn*Grad^s(vel)
    V_out(1,:,:,:)   = part_visc_sym_grad_1
    V_out(2,:,1:2,:) = part_visc_sym_grad_1(:,3:4,:)
    V_out(2,:,3:6,:) = part_visc_sym_grad_2(:,1:4,:)
    V_out(3,:,1:2,:) = part_visc_sym_grad_1(:,5:6,:)
    V_out(3,:,3:4,:) = part_visc_sym_grad_2(:,3:4,:)
    V_out(3,:,5:6,:) = part_visc_sym_grad_2(:,5:6,:)

  END SUBROUTINE smb_explicit_diffu_sym

  SUBROUTINE smb_explicit_LES(communicator, mesh, list_mode, nb_procs, visc_entro_real, mom, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: visc_entro_real
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: mom
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: grad1_mom, grad2_mom, grad3_mom
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: mom_loc
    REAL(KIND=8)                                :: ray, hh, hm
    INTEGER                                     :: m, l , i, mode, index, k
    INTEGER                                     :: m_max_pad, bloc_size
    MPI_Comm          :: communicator

    CALL gauss(mesh)

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             mom_loc(:,k) = mom(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute local mesh sizes
             hh=mesh%hloc_gauss(index)
             hm=MIN(mesh%hm(i),hh) !WRONG choice
             !hm=0.5d0/inputs%m_max
             !hm=mesh%hm(i) !(JLG April 7 2017)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------Grad mom_r on Gauss points------------------------------------
             grad1_mom(index,1,i) = SUM(mom_loc(:,1)*dw_loc(1,:))*hh
             grad1_mom(index,2,i) = SUM(mom_loc(:,2)*dw_loc(1,:))*hh
             grad1_mom(index,3,i) = (mode*SUM(mom_loc(:,2)*ww(:,l))  - &
                  SUM(mom_loc(:,3)*ww(:,l)))/ray*hm
             grad1_mom(index,4,i) = (-mode*SUM(mom_loc(:,1)*ww(:,l)) - &
                  SUM(mom_loc(:,4)*ww(:,l)))/ray*hm
             grad1_mom(index,5,i) =  SUM(mom_loc(:,1)*dw_loc(2,:))*hh
             grad1_mom(index,6,i) =  SUM(mom_loc(:,2)*dw_loc(2,:))*hh

             !-----------------Grad mom_th on Gauss points-----------------------------------
             grad2_mom(index,1,i) = SUM(mom_loc(:,3)*dw_loc(1,:))*hh
             grad2_mom(index,2,i) = SUM(mom_loc(:,4)*dw_loc(1,:))*hh
             grad2_mom(index,3,i) = (mode*SUM(mom_loc(:,4)*ww(:,l))  + &
                  SUM(mom_loc(:,1)*ww(:,l)))/ray*hm
             grad2_mom(index,4,i) = (-mode*SUM(mom_loc(:,3)*ww(:,l)) + &
                  SUM(mom_loc(:,2)*ww(:,l)))/ray*hm
             grad2_mom(index,5,i) = SUM(mom_loc(:,3)*dw_loc(2,:))*hh
             grad2_mom(index,6,i) = SUM(mom_loc(:,4)*dw_loc(2,:))*hh

             !-----------------Grad mom_z on Gauss points------------------------------------
             grad3_mom(index,1,i) = SUM(mom_loc(:,5)*dw_loc(1,:))*hh
             grad3_mom(index,2,i) = SUM(mom_loc(:,6)*dw_loc(1,:))*hh
             grad3_mom(index,3,i) = mode*SUM(mom_loc(:,6)*ww(:,l))/ray*hm
             grad3_mom(index,4,i) = -mode*SUM(mom_loc(:,5)*ww(:,l))/ray*hm
             grad3_mom(index,5,i) = SUM(mom_loc(:,5)*dw_loc(2,:))*hh
             grad3_mom(index,6,i) = SUM(mom_loc(:,6)*dw_loc(2,:))*hh
          ENDDO
       ENDDO
    END DO

    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(grad1_mom,1)/nb_procs+1
    bloc_size = mesh%gauss%l_G*(bloc_size/mesh%gauss%l_G)+mesh%gauss%l_G
    CALL FFT_COMPUTE_ENTROPY_VISC_GRAD_MOM(communicator,grad1_mom, grad2_mom, grad3_mom, visc_entro_real, &
         V_out, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_explicit_LES

  SUBROUTINE smb_surface_tension(communicator, mesh, list_mode, nb_procs, level_set, tensor_surface_gauss)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                       :: mesh
    INTEGER,                          INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),       INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)  :: level_set
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: tensor_surface_gauss
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))              :: tensor
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))                :: grad_level
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                :: level_set_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, l , i, mode, index, k, n
    INTEGER                                     :: m_max_pad, bloc_size
    MPI_Comm          :: communicator

    CALL gauss(mesh)

    tensor_surface_gauss = 0.d0
    DO n = 1, inputs%nb_fluid-1
       DO i = 1, SIZE(list_mode)
          mode = list_mode(i)
          index = 0
          DO m = 1, mesh%dom_me
             j_loc = mesh%jj(:,m)
             DO k = 1, 2
                level_set_loc(:,k) = level_set(n,j_loc,k,i)
             END DO
             DO l = 1, l_G
                index = index + 1
                dw_loc = dw(:,:,l,m)

                !===Compute radius of Gauss point
                ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

                !-----------------Grad u_r on Gauss points------------------------------------
                grad_level(index,1,i) = SUM(level_set_loc(:,1)*dw_loc(1,:))
                grad_level(index,2,i) = SUM(level_set_loc(:,2)*dw_loc(1,:))
                grad_level(index,3,i) = mode/ray*SUM(level_set_loc(:,2)*ww(:,l))
                grad_level(index,4,i) = -mode/ray*SUM(level_set_loc(:,1)*ww(:,l))
                grad_level(index,5,i) =  SUM(level_set_loc(:,1)*dw_loc(2,:))
                grad_level(index,6,i) =  SUM(level_set_loc(:,2)*dw_loc(2,:))
             END DO
          ENDDO
       ENDDO
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       bloc_size = SIZE(grad_level,1)/nb_procs+1
       CALL FFT_TENSOR_DCL(communicator, grad_level, grad_level, tensor, nb_procs, bloc_size, m_max_pad, opt_tension=.TRUE.)
       tensor_surface_gauss = tensor_surface_gauss + inputs%coeff_surface(n)*tensor
    END DO

  END SUBROUTINE smb_surface_tension

  SUBROUTINE smb_buoyancy(communicator, mesh_P2, mesh_P1, list_mode, level_set, temperature, V_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE subroutine_mass
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                       :: mesh_P2, mesh_P1
    INTEGER,      DIMENSION(:),       INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)  :: level_set
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)  :: temperature
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh_P2%np,2,SIZE(list_mode)) :: grav_density
    REAL(KIND=8)                                :: tps
    REAL(KIND=8), DIMENSION(3)                  :: temps
    INTEGER                                     :: code, m_max_pad, bloc_size, nb_procs
    !MPI_Comm       :: communicator
    MPI_Comm, DIMENSION(:), POINTER  :: communicator

    !===Set buoyancy to zero (for radial and azimuthal components)
    V_out = 0.d0

    !===Compute gravity_temp*density using level_set
    CALL reconstruct_variable(communicator, list_mode, mesh_P1, mesh_P2, level_set, &
         inputs%heat_grav_fluid*inputs%density_fluid, grav_density)

    !===Compute vertical component of buoyancy using FFT
    temps=0.d0
    CALL MPI_COMM_SIZE(communicator(2), nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(temperature,1)/nb_procs+1
    CALL FFT_PAR_PROD_DCL(communicator(2), grav_density, temperature, V_out(:,5:6,:), nb_procs, bloc_size, m_max_pad, temps)

  END SUBROUTINE smb_buoyancy

  SUBROUTINE momentum_dirichlet(communicator, mesh, list_mode, t, nb_procs, density, momentum_exact, vv_js_D)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    REAL(KIND=8),                   INTENT(IN)  :: t
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: density
    TYPE(dyn_int_line), DIMENSION(3),INTENT(IN) :: vv_js_D
    REAL(KIND=8), DIMENSION(mesh%np,6,SIZE(list_mode)), INTENT(OUT) :: momentum_exact
    REAL(KIND=8), DIMENSION(SIZE(vv_js_D(1)%DIL),2,SIZE(list_mode))   :: vel_exact_r, mr
    REAL(KIND=8), DIMENSION(SIZE(vv_js_D(2)%DIL),2,SIZE(list_mode))   :: vel_exact_t, mt
    REAL(KIND=8), DIMENSION(SIZE(vv_js_D(3)%DIL),2,SIZE(list_mode))   :: vel_exact_z, mz
    INTEGER                                     :: i, k, kk
    INTEGER                                     :: m_max_pad, bloc_size
    MPI_Comm          :: communicator

    IF (inputs%if_level_set) THEN
       DO i = 1, SIZE(list_mode)
          DO k = 1, 2
             vel_exact_r(:,k,i) = vv_exact(k,  mesh%rr(:,vv_js_D(1)%DIL),list_mode(i),t)
             vel_exact_t(:,k,i) = vv_exact(k+2,mesh%rr(:,vv_js_D(2)%DIL),list_mode(i),t)
             vel_exact_z(:,k,i) = vv_exact(k+4,mesh%rr(:,vv_js_D(3)%DIL),list_mode(i),t)
          END DO
       END DO

       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       bloc_size = SIZE(vv_js_D(1)%DIL)/nb_procs+1
       CALL FFT_PAR_PROD_DCL(communicator, density(vv_js_D(1)%DIL,:,:), vel_exact_r, &
            mr, nb_procs, bloc_size, m_max_pad)
       momentum_exact(vv_js_D(1)%DIL,1:2,:) = mr

       bloc_size = SIZE(vv_js_D(2)%DIL)/nb_procs+1
       CALL FFT_PAR_PROD_DCL(communicator, density(vv_js_D(2)%DIL,:,:), vel_exact_t, &
            mt, nb_procs, bloc_size, m_max_pad)
       momentum_exact(vv_js_D(2)%DIL,3:4,:) = mt

       bloc_size = SIZE(vv_js_D(3)%DIL)/nb_procs+1
       CALL FFT_PAR_PROD_DCL(communicator, density(vv_js_D(3)%DIL,:,:), vel_exact_z, &
            mz, nb_procs, bloc_size, m_max_pad)
       momentum_exact(vv_js_D(3)%DIL,5:6,:) = mz
    ELSE
       DO i = 1, SIZE(list_mode)
          DO k = 1, 6
             kk = (k-1)/2 + 1
             momentum_exact(vv_js_D(kk)%DIL,k,i) = vv_exact(k,mesh%rr(:,vv_js_D(kk)%DIL),list_mode(i),t)
          END DO
       END DO
    END IF

  END SUBROUTINE momentum_dirichlet

  SUBROUTINE twoD_volume(communicator,mesh,RESLT)
    !===========================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    REAL(KIND=8),                INTENT(OUT) :: RESLT
    REAL(KIND=8)                             :: vol_loc, vol_out
    INTEGER,      DIMENSION(mesh%gauss%n_w)  :: j_loc
    INTEGER                                  ::  m, l , i , ni, code
    REAL(KIND=8)                             :: ray
    MPI_Comm                                 :: communicator
    vol_loc = 0.d0
    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)
       DO l = 1, mesh%gauss%l_G
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = j_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          vol_loc = vol_loc + ray*mesh%gauss%rj(l,m)
       ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(vol_loc,vol_out,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    RESLT = vol_out
  END SUBROUTINE TwoD_volume

  SUBROUTINE Moy(communicator,mesh,p,RESLT)
    !===========================
    !===average of p
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

END MODULE subroutine_ns_with_m
