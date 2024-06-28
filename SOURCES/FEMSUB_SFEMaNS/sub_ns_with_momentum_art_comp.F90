!
!Authors LC 2019/11 based on sub_ns_with_vel
!
MODULE subroutine_ns_with_m_art_comp
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
  USE subroutine_ns_with_m
  PUBLIC :: BDF1_art_comp_with_m, smb_compute_tensor_gauss
  PRIVATE
CONTAINS

  SUBROUTINE BDF1_art_comp_with_m(comm_one_d, time, vv_3_LA, pp_1_LA, vvz_per, pp_per, &
       dt, Re, list_mode, pp_mesh, vv_mesh, pn_m1, pn, un_m1, un, Hn_p2, Bn_p2, tempn, concn, &
       density_m1, density, density_p1, visco_dyn, level_set_p1, visc_entro_level, level_set_reg, &
       visc_entro_grad_mom)
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
    USE rhs_gauss_computing
    USE rhs_para_assembling
    USE tn_axi
    USE verbose
    USE entropy_viscosity
    USE subroutine_mass
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8)                                           :: time, dt, Re
    INTEGER,      DIMENSION(:),              INTENT(IN)    :: list_mode
    TYPE(mesh_type),                         INTENT(IN)    :: pp_mesh, vv_mesh
    TYPE(petsc_csr_LA)                                     :: vv_3_LA, pp_1_LA
    TYPE(periodic_type),                     INTENT(IN)    :: vvz_per, pp_per
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(INOUT) :: pn_m1, pn
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(INOUT) :: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN)    :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN)    :: concn
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN)    :: Hn_p2, Bn_p2
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN)    :: density_m1, density, density_p1
    REAL(KIND=8), DIMENSION(:,:,:,:),        INTENT(IN)    :: level_set_p1
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN)    :: visco_dyn
    REAL(KIND=8), DIMENSION(:,:),            INTENT(OUT)   :: visc_entro_level
    REAL(KIND=8), DIMENSION(:,:,:,:),        INTENT(IN)    :: level_set_reg
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(:,:,:,:),        INTENT(INOUT) :: visc_entro_grad_mom
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))  :: Hn_p2_aux
    !===Saved variables
    INTEGER,                                       SAVE :: m_max_c
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: pp_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: pp_mode_global_js_D
    TYPE(dyn_int_line), DIMENSION(3),              SAVE :: vv_js_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: vv_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: vel_global_D
    LOGICAL,                                       SAVE :: once = .TRUE.
!TEST LC LES_SUITE 2024/06
    LOGICAL,                                       SAVE :: once_LES = .TRUE.
!TEST LC LES_SUITE 2024/06
    INTEGER,                                       SAVE :: my_petscworld_rank
    REAL(KIND=8),                                  SAVE :: nu_bar
    REAL(KIND=8),                                  SAVE :: stab_bar
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE,   SAVE :: momentum, momentum_m1, momentum_m2
    INTEGER,                                       SAVE :: bloc_size, m_max_pad, nb_procs
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: visc_grad_vel_m1
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: tensor_m1_gauss
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE,   SAVE :: visc_entro_real
    !===End saved variables

    !----------Declaration sans save-------------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: pp_1_ifrom, vv_3_ifrom

    INTEGER                                  :: i, k, m, n, n1, n2, n3, n123, nb_inter
    INTEGER                                  :: code, nu_mat, mode
    REAL(KIND=8)                             :: moyenne
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: div
!!$    REAL(KIND=8), DIMENSION(pp_mesh%np, 2, SIZE(list_mode)) :: visco_dyn_P1, visc_div
!!$    INTEGER :: bloc_size_P1
    REAL(KIND=8), DIMENSION(pp_mesh%np, 2)   :: pn_p1
    REAL(KIND=8), DIMENSION(vv_mesh%np, 6)   :: un_p1
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: rotb_b, rotb_b_aux
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: visc_grad_vel
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: stab_grad_mom
!TEST LC LES_SUITE 2024/06
!    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: visc_entro_grad_mom
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,2,SIZE(list_mode))   :: stab_div_vel
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: tensor_surface_gauss
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)) :: tensor_gauss
    REAL(KIND=8), DIMENSION(3,vv_mesh%np,6,SIZE(list_mode))                       :: tensor
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: rhs_gauss
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6)                 :: visc_grad_vel_ext
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6)                 :: tensor_gauss_ext
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)) :: uext, momentumext, momentum_exact, momentumext_div
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                         :: buoyancy
    REAL(KIND=8), DIMENSION(inputs%nb_fluid-1, vv_mesh%np, 2, SIZE(list_mode))    :: level_set_FEM_P2

    REAL(KIND=8), DIMENSION(vv_mesh%np) :: vel_loc, vel_tot
    REAL(KIND=8)   :: coeff, vloc, cfl, cfl_max, norm
    INTEGER        :: nb_procs_LES, bloc_size_LES, m_max_pad_LES
    !April 17th 2008, JLG
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/

    !Communicators for Petsc, in space and Fourier------------------------------
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: vel_mat
    Mat,                        SAVE :: mass_mat, mass_mat0
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Vec,                        SAVE :: pb_1,    pb_2,    px_1, px_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: vel_ksp
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

       !===Momentum Initialization
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
       !===End Momentum Initialization

       !===Tensors_m1 allocation and initialization
       ALLOCATE(tensor_m1_gauss(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)))
       CALL smb_compute_tensor_gauss(comm_one_d(2), vv_mesh, list_mode, un_m1, momentum_m1, &
            nb_procs, bloc_size, m_max_pad, tensor, tensor_m1_gauss)

       ALLOCATE(visc_grad_vel_m1(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)))
       CALL smb_explicit_diffu_sym(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            visco_dyn/Re, un_m1, visc_grad_vel_m1)
       !===End Tensors_m1 allocation and initialization

       !===Allocation entropy visc
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs_LES, code)
       bloc_size_LES = vv_mesh%gauss%l_G*vv_mesh%dom_me/nb_procs_LES+1
       bloc_size_LES = vv_mesh%gauss%l_G*(bloc_size_LES/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
       m_max_pad_LES = 3*SIZE(list_mode)*nb_procs_LES/2
       ALLOCATE(visc_entro_real(2*m_max_pad_LES-1,bloc_size_LES))
       visc_entro_real = 0.d0
       !===End Allocation entropy visc

       !===PREPARE pp_mode_global_js_D ARRAY FOR PRESSURE
       !===ATTENTION pressure BCs are no longer implemented
       !===JLG June 9 2017
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

       !===ASSEMBLING VELOCITY MATRICES
       ALLOCATE(vel_mat(2*m_max_c),vel_ksp(2*m_max_c))

       ! Definition of nu_bar
       IF (inputs%if_level_set) THEN
          nu_bar = 0.d0
          stab_bar = 0.d0
          DO n = 1, inputs%nb_fluid
             nu_bar = MAX(nu_bar,inputs%dyna_visc_fluid(n)/inputs%density_fluid(n))
             stab_bar = MAX(stab_bar,1.d0/inputs%density_fluid(n))
          END DO
          !LC: To be tested thoroughly
          !nu_bar    = 2.0d0*nu_bar
          nu_bar    = 1.1d0*nu_bar
          stab_bar  = 1.1d0*stab_bar

          ! Normalization penal_coeff_art_comp
          inputs%penal_coeff_art_comp=inputs%penal_coeff_art_comp*MAX(1.d0,nu_bar/stab_bar)
       ELSE
          nu_bar    = 1.d0
          stab_bar  = 1.d0
       END IF

       DO i = 1, m_max_c
          mode = list_mode(i)
          !===VELOCITY
          nu_mat = 2*i-1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          IF (inputs%if_moment_bdf2) THEN
             CALL qs_diff_mass_vect_3x3_divpenal_art_comp (1, vv_3_LA, vv_mesh, nu_bar/Re, three/(2*dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, stab_bar*inputs%penal_coeff_art_comp, &
                  i, mode, vel_mat(nu_mat))
          ELSE
             CALL qs_diff_mass_vect_3x3_divpenal_art_comp (1, vv_3_LA, vv_mesh, nu_bar/Re, one/(dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, stab_bar*inputs%penal_coeff_art_comp, &
                  i, mode, vel_mat(nu_mat))
          END IF
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  vel_mat(nu_mat), vv_3_LA)
          END IF
          CALL Dirichlet_M_parallel(vel_mat(nu_mat),vv_mode_global_js_D(i)%DIL)
          CALL init_solver(inputs%my_par_vv,vel_ksp(nu_mat),vel_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
          nu_mat = nu_mat+1
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA, vel_mat(nu_mat), clean=.FALSE.)
          IF (inputs%if_moment_bdf2) THEN
             CALL qs_diff_mass_vect_3x3_divpenal_art_comp (2, vv_3_LA, vv_mesh, nu_bar/Re, three/(2*dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, stab_bar*inputs%penal_coeff_art_comp, &
                  i, mode, vel_mat(nu_mat))
          ELSE
             CALL qs_diff_mass_vect_3x3_divpenal_art_comp (2, vv_3_LA, vv_mesh, nu_bar/Re, one/(dt), &
                  inputs%LES_coeff1_mom, inputs%stab_bdy_ns, stab_bar*inputs%penal_coeff_art_comp, &
                  i, mode, vel_mat(nu_mat))
          END IF
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

    !===Compute rhs by FFT at Gauss points

    !===Compute extrapolation velocity and momentum
    IF (inputs%if_moment_bdf2) THEN
       uext = 2*un-un_m1
       IF (inputs%if_level_set) THEN
          !Second time order: momentumext = rho_np1*(2*un-un_m1)
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), uext, density_p1, momentumext, 1,nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentumext = 2*momentum - momentum_m1
       END IF
       !TEST LC (above 2nd order extrapolation not stable)
       IF (inputs%if_level_set) THEN
          !First time order: momentumext = rho_np1*un
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), un, density_p1, momentumext_div, 1, nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentumext_div = momentum
       END IF
       !END TEST LC
    ELSE
       uext = un
       IF (inputs%if_level_set) THEN
          !First time order: momentumext = rho_np1*un
          CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), un, density_p1, momentumext, 1, nb_procs, &
               bloc_size, m_max_pad)
       ELSE
          momentumext = momentum
       END IF
       momentumext_div=momentumext
    END IF
    !===End Compute extrapolation velocity and momentum

    !===Precession should be in condlim
    IF (inputs%precession) THEN
       CALL error_petsc('for momentum ns: precession should be in condlim')
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

    !===Compute diffusion/artificial compression corrections and surface tension
    IF (inputs%if_level_set) THEN
       !===Compute visco_dyn/Re*Grad(un) and nu_bar*Grad(momentum)
       CALL smb_explicit_diffu_correction(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            visco_dyn/Re, nu_bar/Re, un, momentumext, visc_grad_vel, stab_grad_mom)
       !===End compute visco_dyn/Re*Grad(un) and nu_bar/Re*Grad(momentum)

       !===Compute Div(un) - stab_bar*Div(momentum)
       CALL smb_explicit_div_correction(comm_one_d(2), vv_mesh, list_mode, stab_bar, un, momentumext_div, stab_div_vel)
       stab_div_vel=inputs%penal_coeff_art_comp*stab_div_vel
       !===End Compute Div(un) - stab_bar*Div(momentum)

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
       visc_grad_vel = 0.d0
       stab_grad_mom = 0.d0
       stab_div_vel  = 0.d0
       tensor_surface_gauss = 0.d0
       buoyancy = 0.d0 !Defined in condlim.F90
    END IF
    !===End Compute diffusion/ artificial compression corrections and surface tension

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
    CALL smb_compute_tensor_gauss(comm_one_d(2), vv_mesh, list_mode, un, momentum, &
         nb_procs, bloc_size, m_max_pad, tensor, tensor_gauss)
    !===End Compute tensor product of momentum by velocity

    !===PREPARE BOUNDARY CONDITION FOR MOMENTUM
    CALL momentum_dirichlet(comm_one_d(2), vv_mesh, list_mode, time, nb_procs, density_p1, &
         momentum_exact, vv_js_D)
    !===End PREPARE BOUNDARY CONDITION FOR MOMENTUM

    !===End Compute rhs by FFT at Gauss points

    !===Computation of rhs at Gauss points for every mode
    IF (inputs%if_moment_bdf2) THEN
       CALL rhs_ns_gauss_3x3_art_comp_mom(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
            (4*momentum-momentum_m1)/(2*inputs%dt), pn, -rotb_b, rhs_gauss, tempn, concn, density_p1, &
            buoyancy)
    ELSE
       CALL rhs_ns_gauss_3x3_art_comp_mom(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time, &
            (momentum)/(inputs%dt), pn, -rotb_b, rhs_gauss, tempn, concn, density_p1, &
            buoyancy)
    END IF
    !===End Computation of rhs at Gauss points for every mode

    !------------BEGIN LOOP ON FOURIER MODES FOR MOMENTUM---------------
    DO i = 1, m_max_c
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       IF (inputs%if_moment_bdf2) THEN

          !===Compute 2nd time order explicit tensors
          tensor_gauss_ext         = 2*tensor_gauss(:,:,:,i) - tensor_m1_gauss(:,:,:,i)
          visc_grad_vel_ext        = 2*visc_grad_vel(:,:,:,i) - visc_grad_vel_m1(:,:,:,i)
          !===Update terms at m1
          tensor_m1_gauss(:,:,:,i)  = tensor_gauss(:,:,:,i)
          visc_grad_vel_m1(:,:,:,i) = visc_grad_vel(:,:,:,i)

          CALL rhs_3x3_art_comp(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236, &
               opt_tensor=tensor_gauss_ext+tensor_surface_gauss(:,:,:,i) &
               -visc_grad_vel_ext+stab_grad_mom(:,:,:,i)-visc_entro_grad_mom(:,:,:,i), &
               opt_grad_div=-stab_div_vel(:,:,i))
       ELSE
          CALL rhs_3x3_art_comp(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236, &
               opt_tensor=tensor_gauss(:,:,:,i)+tensor_surface_gauss(:,:,:,i) &
               -visc_grad_vel(:,:,:,i)+stab_grad_mom(:,:,:,i)-visc_entro_grad_mom(:,:,:,i), &
               opt_grad_div=-stab_div_vel(:,:,i))
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
       vel_global_D(i)%DRL(1:n1)        = momentum_exact(vv_js_D(1)%DIL,1,i)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = momentum_exact(vv_js_D(2)%DIL,4,i)
       vel_global_D(i)%DRL(n1+n2+1:n123)= momentum_exact(vv_js_D(3)%DIL,5,i)
       vel_global_D(i)%DRL(n123+1:)     = 0.D0
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_145)
       vel_global_D(i)%DRL(1:n1)        = momentum_exact(vv_js_D(1)%DIL,2,i)
       vel_global_D(i)%DRL(n1+1:n1+n2)  = momentum_exact(vv_js_D(2)%DIL,3,i)
       vel_global_D(i)%DRL(n1+n2+1:n123)= momentum_exact(vv_js_D(3)%DIL,6,i)
       CALL dirichlet_rhs(vv_mode_global_js_D(i)%DIL-1,vel_global_D(i)%DRL,vb_3_236)
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss

       !===Solve linear system for momentum equation
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
       !===End Solve linear system for momentum equation

       !===Correction of zero mode
       IF (mode==0) THEN
          un_p1 (:,2) = 0.d0
          un_p1 (:,4) = 0.d0
          un_p1 (:,6) = 0.d0
       END IF
       !===End Correction of zero mode

       !===Update momentum
       momentum_m2(:,:,i) = momentum_m1 (:,:,i)
       momentum_m1(:,:,i) = momentum (:,:,i)
       momentum (:,:,i) = un_p1
       !===End Update momentum
    END DO
    !------------END LOOP ON FOURIER MODES FOR MOMENTUM---------------

    !===Update Velocity
    un_m1 = un
    IF (inputs%if_level_set) THEN
       CALL FFT_SCALAR_VECT_DCL(comm_one_d(2), momentum, density_p1, un, 2, nb_procs, &
            bloc_size, m_max_pad)
    ELSE
       un = momentum
    END IF
    !===End Update Velocity

    !------------BEGIN LOOP ON FOURIER MODES FOR DIVERGENCE P1-----------
    DO i = 1, m_max_c
       mode = list_mode(i)

       !===Assemble divergence of velocity in arrays pb_1, pb_2
       CALL qs_01_div_hybrid_2006(vv_mesh, pp_mesh, pp_1_LA, mode, un(:,:,i), pb_1, pb_2)
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

       !===Solve mass matrix for pressure correction
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
    END DO
    !------------END LOOP ON FOURIER MODES FOR DIVERGENCE P1-------------

!!$    !===Compute visco_dyn*div on P1 elements
!!$    DO i = 1, m_max_c
!!$       DO k=1, 2
!!$          CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, visco_dyn(:,k,i), visco_dyn_P1(:,k,i))
!!$       END DO
!!$    END DO
!!$    bloc_size_P1 = SIZE(div,1)/nb_procs+1
!!$    CALL FFT_PAR_PROD_DCL(comm_one_d(2), visco_dyn_P1, div, visc_div, nb_procs, bloc_size_P1, m_max_pad)
!!$    !===End Compute visco_dyn*div on P1 elements

    !------------BEGIN LOOP ON FOURIER MODES FOR PRESSURE----------------
    DO i = 1, m_max_c
       mode = list_mode(i)
       !===Pressure computation
       pn_p1(:,:) = pn(:,:,i)
       DO k=1, 2
          pn_p1(:,k) = pn_p1(:,k) - inputs%penal_coeff_art_comp*div(:,k,i)
       END DO
       !===End Pressure computation

       !===Handling of mean pressure
       IF (mode == 0)  THEN
          CALL Moy(comm_one_d(1),pp_mesh, pn_p1(:,1),moyenne)
          pn_p1(:,1) = pn_p1(:,1)-moyenne
       ENDIF
       !===End of handling of mean pressure

       !===Correction of zero mode
       IF (mode==0) THEN
          pn_p1 (:,2) = 0.d0
       END IF
       !===Correction of zero mode

       !===Update pressures
       pn_m1(:,:,i)    = pn(:,:,i)
       pn   (:,:,i)    = pn_p1
       !===End Update pressures
    ENDDO
    !------------END LOOP ON FOURIER MODES FOR PRESSURE------------------

    !===Verbose divergence of velocity
    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'H1',  vv_mesh, list_mode, un)
       talk_to_me%div_L2  = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)/norm
       talk_to_me%weak_div_L2  = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, div)/norm
    END IF
    !===End verbose divergence of velocity

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
    !===End Compute entropy viscosity

!TEST LC LES_SUITE 2024/06
    !===Compute entropy viscosity stablization: (-LES_coeff1+visc_entro)*Grad(momentumext)
    IF (inputs%LES.AND.inputs%if_LES_in_momentum) THEN
       CALL smb_explicit_LES(comm_one_d(2), vv_mesh, list_mode, nb_procs, &
            visc_entro_real, momentumext, visc_entro_grad_mom)
    ELSE
       visc_entro_grad_mom=0.d0
    END IF
    !===End compute entropy viscosity stabilization: (-LES_coeff1+visc_entro)*Grad(momentumext)
!TEST LC LES_SUITE 2024/06

  END SUBROUTINE BDF1_art_comp_with_m
  !============================================

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
    temps = 0

    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(RotV,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    CALL FFT_PAR_CROSS_PROD_DCL(communicator, RotV, W, V_out, nb_procs, bloc_size, m_max_pad, temps)
    tps = user_time() - tps

  END SUBROUTINE smb_CurlH_cross_B_gauss_sft_par

  SUBROUTINE smb_explicit_diffu_correction(communicator, mesh, list_mode, nb_procs, visc_dyn, stab_mom, &
       vel, mom, V_out, V2_out)
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
    REAL(KIND=8)                                :: stab_mom
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: mom
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V2_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: grad1_vel, grad2_vel, grad3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: gradT1_vel, gradT2_vel, gradT3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: grad1_mom, grad2_mom, grad3_mom
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: gradT1_mom, gradT2_mom, gradT3_mom
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_tensor_1, part_tensor_2
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_visc_sym_grad_1
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))   :: part_visc_sym_grad_2
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: visc_dyn_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc, mom_loc
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
             mom_loc(:,k) = mom(j_loc,k,i)
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
             !-----------------Grad mom_r on Gauss points----------------------------------
             grad1_mom(index,1,i) = stab_mom*SUM(mom_loc(:,1)*dw_loc(1,:))
             grad1_mom(index,2,i) = stab_mom*SUM(mom_loc(:,2)*dw_loc(1,:))
             grad1_mom(index,3,i) = stab_mom*(mode*SUM(mom_loc(:,2)*ww(:,l)) - SUM(mom_loc(:,3)*ww(:,l)))/ray
             grad1_mom(index,4,i) = stab_mom*(-mode*SUM(mom_loc(:,1)*ww(:,l)) - SUM(mom_loc(:,4)*ww(:,l)))/ray
             grad1_mom(index,5,i) = stab_mom*SUM(mom_loc(:,1)*dw_loc(2,:))
             grad1_mom(index,6,i) = stab_mom*SUM(mom_loc(:,2)*dw_loc(2,:))

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(index,1,i) = SUM(vel_loc(:,3)*dw_loc(1,:))
             grad2_vel(index,2,i) = SUM(vel_loc(:,4)*dw_loc(1,:))
             grad2_vel(index,3,i) = (mode*SUM(vel_loc(:,4)*ww(:,l)) + SUM(vel_loc(:,1)*ww(:,l)))/ray
             grad2_vel(index,4,i) = (-mode*SUM(vel_loc(:,3)*ww(:,l)) + SUM(vel_loc(:,2)*ww(:,l)))/ray
             grad2_vel(index,5,i) = SUM(vel_loc(:,3)*dw_loc(2,:))
             grad2_vel(index,6,i) = SUM(vel_loc(:,4)*dw_loc(2,:))
             !-----------------Grad mom_th on Gauss points---------------------------------
             grad2_mom(index,1,i) = stab_mom*SUM(mom_loc(:,3)*dw_loc(1,:))
             grad2_mom(index,2,i) = stab_mom*SUM(mom_loc(:,4)*dw_loc(1,:))
             grad2_mom(index,3,i) = stab_mom*(mode*SUM(mom_loc(:,4)*ww(:,l)) + SUM(mom_loc(:,1)*ww(:,l)))/ray
             grad2_mom(index,4,i) = stab_mom*(-mode*SUM(mom_loc(:,3)*ww(:,l)) + SUM(mom_loc(:,2)*ww(:,l)))/ray
             grad2_mom(index,5,i) = stab_mom*SUM(mom_loc(:,3)*dw_loc(2,:))
             grad2_mom(index,6,i) = stab_mom*SUM(mom_loc(:,4)*dw_loc(2,:))


             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(index,1,i) = SUM(vel_loc(:,5)*dw_loc(1,:))
             grad3_vel(index,2,i) = SUM(vel_loc(:,6)*dw_loc(1,:))
             grad3_vel(index,3,i) = mode*SUM(vel_loc(:,6)*ww(:,l))/ray
             grad3_vel(index,4,i) = -mode*SUM(vel_loc(:,5)*ww(:,l))/ray
             grad3_vel(index,5,i) = SUM(vel_loc(:,5)*dw_loc(2,:))
             grad3_vel(index,6,i) = SUM(vel_loc(:,6)*dw_loc(2,:))
             !-----------------Grad mom_z on Gauss points----------------------------------
             grad3_mom(index,1,i) = stab_mom*SUM(mom_loc(:,5)*dw_loc(1,:))
             grad3_mom(index,2,i) = stab_mom*SUM(mom_loc(:,6)*dw_loc(1,:))
             grad3_mom(index,3,i) = stab_mom*mode*SUM(mom_loc(:,6)*ww(:,l))/ray
             grad3_mom(index,4,i) = stab_mom*(-mode*SUM(mom_loc(:,5)*ww(:,l)))/ray
             grad3_mom(index,5,i) = stab_mom*SUM(mom_loc(:,5)*dw_loc(2,:))
             grad3_mom(index,6,i) = stab_mom*SUM(mom_loc(:,6)*dw_loc(2,:))
          ENDDO
       ENDDO
    END DO

    !-----------------GradT u_r on Gauss points------------------------------------
    gradT1_vel(:,1,:) = grad1_vel(:,1,:)
    gradT1_vel(:,2,:) = grad1_vel(:,2,:)
    gradT1_vel(:,3,:) = grad2_vel(:,1,:)
    gradT1_vel(:,4,:) = grad2_vel(:,2,:)
    gradT1_vel(:,5,:) = grad3_vel(:,1,:)
    gradT1_vel(:,6,:) = grad3_vel(:,2,:)
    !-----------------GradT mom_r on Gauss points----------------------------------
    gradT1_mom(:,1,:) = grad1_mom(:,1,:)
    gradT1_mom(:,2,:) = grad1_mom(:,2,:)
    gradT1_mom(:,3,:) = grad2_mom(:,1,:)
    gradT1_mom(:,4,:) = grad2_mom(:,2,:)
    gradT1_mom(:,5,:) = grad3_mom(:,1,:)
    gradT1_mom(:,6,:) = grad3_mom(:,2,:)

    !-----------------GradT u_th on Gauss points-----------------------------------
    gradT2_vel(:,1,:) = grad1_vel(:,3,:)
    gradT2_vel(:,2,:) = grad1_vel(:,4,:)
    gradT2_vel(:,3,:) = grad2_vel(:,3,:)
    gradT2_vel(:,4,:) = grad2_vel(:,4,:)
    gradT2_vel(:,5,:) = grad3_vel(:,3,:)
    gradT2_vel(:,6,:) = grad3_vel(:,4,:)
    !-----------------GradT mom_th on Gauss points---------------------------------
    gradT2_mom(:,1,:) = grad1_mom(:,3,:)
    gradT2_mom(:,2,:) = grad1_mom(:,4,:)
    gradT2_mom(:,3,:) = grad2_mom(:,3,:)
    gradT2_mom(:,4,:) = grad2_mom(:,4,:)
    gradT2_mom(:,5,:) = grad3_mom(:,3,:)
    gradT2_mom(:,6,:) = grad3_mom(:,4,:)

    !-----------------GradT u_z on Gauss points------------------------------------
    gradT3_vel(:,1,:) = grad1_vel(:,5,:)
    gradT3_vel(:,2,:) = grad1_vel(:,6,:)
    gradT3_vel(:,3,:) = grad2_vel(:,5,:)
    gradT3_vel(:,4,:) = grad2_vel(:,6,:)
    gradT3_vel(:,5,:) = grad3_vel(:,5,:)
    gradT3_vel(:,6,:) = grad3_vel(:,6,:)
    !-----------------GradT mom_z on Gauss points----------------------------------
    gradT3_mom(:,1,:) = grad1_mom(:,5,:)
    gradT3_mom(:,2,:) = grad1_mom(:,6,:)
    gradT3_mom(:,3,:) = grad2_mom(:,5,:)
    gradT3_mom(:,4,:) = grad2_mom(:,6,:)
    gradT3_mom(:,5,:) = grad3_mom(:,5,:)
    gradT3_mom(:,6,:) = grad3_mom(:,6,:)

    !===Grad_vel = Grad_vel + Grad_vel^T
    grad1_vel = grad1_vel + gradT1_vel
    grad2_vel = grad2_vel + gradT2_vel
    grad3_vel = grad3_vel + gradT3_vel
    !===Grad_mom = Grad_mom + Grad_mom^T
    grad1_mom = grad1_mom + gradT1_mom
    grad2_mom = grad2_mom + gradT2_mom
    grad3_mom = grad3_mom + gradT3_mom

    !===Compute (visc_dyn)*Grad^s(vel)
    part_tensor_1          = grad1_vel(:,:,:)
    part_tensor_2(:,1:4,:) = grad2_vel(:,3:6,:)
    part_tensor_2(:,5:6,:) = grad3_vel(:,5:6,:)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(grad1_vel,1)/nb_procs+1
    CALL FFT_COMPUTE_DIFFU_MOM(communicator,part_tensor_1, part_tensor_2, visc_dyn_gauss, &
         part_visc_sym_grad_1, part_visc_sym_grad_2, nb_procs, bloc_size, m_max_pad)

    !===Using strain rate tensor symmetry to get V_out = visc_dyn*Grad^s(vel)
    V_out(1,:,:,:)   = part_visc_sym_grad_1
    V_out(2,:,1:2,:) = part_visc_sym_grad_1(:,3:4,:)
    V_out(2,:,3:6,:) = part_visc_sym_grad_2(:,1:4,:)
    V_out(3,:,1:2,:) = part_visc_sym_grad_1(:,5:6,:)
    V_out(3,:,3:4,:) = part_visc_sym_grad_2(:,3:4,:)
    V_out(3,:,5:6,:) = part_visc_sym_grad_2(:,5:6,:)

    !===Compute V2_out = stab_mom*Grad^s(mom)
    V2_out(1,:,:,:) = grad1_mom
    V2_out(2,:,:,:) = grad2_mom
    V2_out(3,:,:,:) = grad3_mom

  END SUBROUTINE smb_explicit_diffu_correction

  SUBROUTINE smb_explicit_div_correction(communicator, mesh, list_mode, stab, vel, mom, c_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8),                   INTENT(IN)  :: stab
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: mom
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: div_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: div_mom
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: mom_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, l , i, mode, index, k
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    MPI_Comm       :: communicator
    CALL MPI_Comm_rank(communicator,rank,ierr)
    CALL gauss(mesh)

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vel_loc(:,k) = vel(j_loc,k,i)
             mom_loc(:,k) = mom(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------Divergence vel on Gauss points---------------------------------
             div_vel(index,1,i) = SUM(vel_loc(:,1)*dw_loc(1,:)) &
                  + SUM(vel_loc(:,1)*ww(:,l))/ray &
                  + mode*SUM(vel_loc(:,4)*ww(:,l))/ray &
                  + SUM(vel_loc(:,5)*dw_loc(2,:))
             div_vel(index,2,i) = SUM(vel_loc(:,2)*dw_loc(1,:)) &
                  + SUM(vel_loc(:,2)*ww(:,l))/ray &
                  - mode*SUM(vel_loc(:,3)*ww(:,l))/ray &
                  + SUM(vel_loc(:,6)*dw_loc(2,:))

             !-----------------Divergence mom on Gauss points---------------------------------
             div_mom(index,1,i) = SUM(mom_loc(:,1)*dw_loc(1,:)) &
                  + SUM(mom_loc(:,1)*ww(:,l))/ray &
                  + mode*SUM(mom_loc(:,4)*ww(:,l))/ray &
                  + SUM(mom_loc(:,5)*dw_loc(2,:))
             div_mom(index,2,i) = SUM(mom_loc(:,2)*dw_loc(1,:)) &
                  + SUM(mom_loc(:,2)*ww(:,l))/ray &
                  - mode*SUM(mom_loc(:,3)*ww(:,l))/ray &
                  + SUM(mom_loc(:,6)*dw_loc(2,:))
          ENDDO
       ENDDO
    END DO

    !===Compute c_out
    c_out = div_vel - stab*div_mom

  END SUBROUTINE smb_explicit_div_correction

!!$  SUBROUTINE smb_explicit_div_correction_new(communicator, mesh, list_mode, nb_procs, visc_dyn, stab, vel, mom, c_out)
!!$    USE Gauss_points
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE boundary
!!$    USE input_data
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type), TARGET                     :: mesh
!!$    INTEGER,                        INTENT(IN)  :: nb_procs
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: visc_dyn
!!$    REAL(KIND=8),                   INTENT(IN)  :: stab
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: mom
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)), INTENT(OUT) :: c_out
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: div_vel
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: div_mom
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: visc_div_vel
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode))   :: visc_dyn_gauss
!!$    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: mom_loc
!!$    REAL(KIND=8)                                :: ray
!!$    INTEGER                                     :: m, l , i, mode, index, k
!!$    INTEGER                                     :: m_max_pad, bloc_size
!!$
!!$    PetscErrorCode :: ierr
!!$    PetscMPIInt    :: rank
!!$    MPI_Comm       :: communicator
!!$    CALL MPI_Comm_rank(communicator,rank,ierr)
!!$    CALL gauss(mesh)
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO m = 1, mesh%dom_me
!!$          j_loc = mesh%jj(:,m)
!!$          DO k = 1, 6
!!$             vel_loc(:,k) = vel(j_loc,k,i)
!!$             mom_loc(:,k) = mom(j_loc,k,i)
!!$          END DO
!!$          DO l = 1, l_G
!!$             index = index + 1
!!$             dw_loc = dw(:,:,l,m)
!!$
!!$             !===Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))
!!$
!!$             !-----------------Dynamic viscosity on Gauss points---------------------------
!!$             visc_dyn_gauss(index,1,i) = SUM(visc_dyn(j_loc,1,i)*ww(:,l))
!!$             visc_dyn_gauss(index,2,i) = SUM(visc_dyn(j_loc,2,i)*ww(:,l))
!!$
!!$             !-----------------Divergence vel on Gauss points---------------------------------
!!$             div_vel(index,1,i) = SUM(vel_loc(:,1)*dw_loc(1,:)) &
!!$                  + SUM(vel_loc(:,1)*ww(:,l))/ray &
!!$                  + mode*SUM(vel_loc(:,4)*ww(:,l))/ray &
!!$                  + SUM(vel_loc(:,5)*dw_loc(2,:))
!!$             div_vel(index,2,i) = SUM(vel_loc(:,2)*dw_loc(1,:)) &
!!$                  + SUM(vel_loc(:,2)*ww(:,l))/ray &
!!$                  - mode*SUM(vel_loc(:,3)*ww(:,l))/ray &
!!$                  + SUM(vel_loc(:,6)*dw_loc(2,:))
!!$
!!$             !-----------------Divergence mom on Gauss points---------------------------------
!!$             div_mom(index,1,i) = SUM(mom_loc(:,1)*dw_loc(1,:)) &
!!$                  + SUM(mom_loc(:,1)*ww(:,l))/ray &
!!$                  + mode*SUM(mom_loc(:,4)*ww(:,l))/ray &
!!$                  + SUM(mom_loc(:,5)*dw_loc(2,:))
!!$             div_mom(index,2,i) = SUM(mom_loc(:,2)*dw_loc(1,:)) &
!!$                  + SUM(mom_loc(:,2)*ww(:,l))/ray &
!!$                  - mode*SUM(mom_loc(:,3)*ww(:,l))/ray &
!!$                  + SUM(mom_loc(:,6)*dw_loc(2,:))
!!$          ENDDO
!!$       ENDDO
!!$    END DO
!!$
!!$    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$    bloc_size = SIZE(div_vel,1)/nb_procs+1
!!$    CALL FFT_PAR_PROD_DCL(communicator, visc_dyn_gauss, div_vel, visc_div_vel, nb_procs, bloc_size, m_max_pad)
!!$
!!$    !===Compute c_out
!!$    c_out = visc_div_vel - stab*div_mom
!!$
!!$  END SUBROUTINE smb_explicit_div_correction_new

  SUBROUTINE smb_compute_tensor_gauss(communicator, mesh, list_mode, vel, mom, nb_procs, &
       bloc_size, m_max_pad, tensor, tensor_gauss)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                       :: mesh
    INTEGER,                          INTENT(IN)  :: nb_procs
    INTEGER,                          INTENT(IN)  :: bloc_size
    INTEGER,                          INTENT(IN)  :: m_max_pad
    INTEGER,      DIMENSION(:),       INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)  :: mom
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: tensor_gauss
    REAL(KIND=8), DIMENSION(3,mesh%np,6,SIZE(list_mode))                   , INTENT(OUT) :: tensor
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: tensor_loc
    INTEGER                                                  :: m, l , i, mode, index, k, n
    MPI_Comm          :: communicator

    CALL gauss(mesh)

    CALL FFT_TENSOR_DCL(communicator, mom, vel, tensor, nb_procs, bloc_size, m_max_pad)

    DO n = 1, 3
       DO i = 1, SIZE(list_mode)
          mode = list_mode(i)
          index = 0
          DO m = 1, mesh%dom_me
             j_loc = mesh%jj(:,m)
             DO k = 1, 6
                tensor_loc(:,k) = tensor(n,j_loc,k,i)
             END DO
             DO l = 1, l_G
                index = index + 1
                !-----------------Tensor on Gauss points-----------------
                tensor_gauss(n,index,1,i) = SUM(tensor_loc(:,1)*ww(:,l))
                tensor_gauss(n,index,2,i) = SUM(tensor_loc(:,2)*ww(:,l))
                tensor_gauss(n,index,3,i) = SUM(tensor_loc(:,3)*ww(:,l))
                tensor_gauss(n,index,4,i) = SUM(tensor_loc(:,4)*ww(:,l))
                tensor_gauss(n,index,5,i) = SUM(tensor_loc(:,5)*ww(:,l))
                tensor_gauss(n,index,6,i) = SUM(tensor_loc(:,6)*ww(:,l))
             END DO
          ENDDO
       ENDDO
    END DO

  END SUBROUTINE smb_compute_tensor_gauss

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

END MODULE subroutine_ns_with_m_art_comp
