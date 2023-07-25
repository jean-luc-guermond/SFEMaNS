MODULE entropy_viscosity
  USE my_util
  PUBLIC :: compute_entropy_viscosity, compute_entropy_viscosity_mom, &
       compute_entropy_viscosity_mom_no_level_set
  PRIVATE
CONTAINS
  SUBROUTINE compute_entropy_viscosity(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, vvz_per, &
       un, un_m1, un_m2, pn_m1, rotv_v_m1, visco_entro_grad_u, density, tempn, concn)
    USE def_type_mesh
    USE fem_M_axi
    USE solve_petsc
    USE periodic
    USE rhs_gauss_computing
    USE rhs_para_assembling
    USE Dir_nodes_petsc
    USE st_matrix
    USE sft_parallele
    USE tn_axi
    USE input_data
    USE my_util
    USE boundary
    USE fourier_to_real_for_vtu
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                    :: vv_3_LA
    TYPE(mesh_type),                INTENT(IN)            :: pp_mesh, vv_mesh
    TYPE(periodic_type),            INTENT(IN)            :: vvz_per
    REAL(KIND=8),                   INTENT(IN)            :: time
    INTEGER,      DIMENSION(:),     INTENT(IN)            :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)         :: un, un_m1, un_m2
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)         :: pn_m1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)            :: rotv_v_m1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)            :: tempn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)            :: concn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)            :: density
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT)         :: visco_entro_grad_u
    TYPE(dyn_int_line), DIMENSION(3),                SAVE :: vv_js_D
    LOGICAL,                                         SAVE :: once = .TRUE.
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:),         SAVE :: zero_dir1
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:),         SAVE :: zero_dir2
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:),         SAVE :: zero_dir3
    REAL(KIND=8),                                    SAVE :: Volume_3D
    INTEGER,                                         SAVE :: iteration

    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: strain_rate_tensor
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_Gs*vv_mesh%dom_mes,6,SIZE(list_mode))   :: strain_rate_tensor_scal_n_bdy
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: un_m1_gauss, res_ns_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                           :: res_ns
    INTEGER,          POINTER, DIMENSION(:)                     :: vv_3_ifrom
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8)    :: norm_vel_L2
    !REAL(KIND=8)    :: norm_res_ns_L2
    !INTEGER         :: rank
    INTEGER         :: i, k, nu_mat, mode, m, l, TYPE, index, n
    INTEGER         :: nb_procs, code, bloc_size, m_max_pad
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Mat, DIMENSION(:), POINTER, SAVE :: LES_mat
    KSP, DIMENSION(:), POINTER, SAVE :: LES_ksp

    IF (once) THEN

       once = .FALSE.

       !===CREATE PETSC VECTORS AND GHOSTS
       CALL create_my_ghost(vv_mesh,vv_3_LA,vv_3_ifrom)
       n = 3*vv_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
       CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
       CALL VecDuplicate(vx_3, vb_3_145, ierr)
       CALL VecDuplicate(vx_3, vb_3_236, ierr)

       !===Compute Volume
       CALL twoD_volume(comm_one_d(1),vv_mesh,Volume_3D)
       Volume_3D = Volume_3D*2*ACOS(-1.d0)

       !===PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY
       DO k = 1, 3
          CALL dirichlet_nodes_parallel(vv_mesh, inputs%vv_list_dirichlet_sides(k)%DIL, vv_js_D(k)%DIL)
       END DO
       !===End PREPARE TYPE OF BOUNDARY CONDITIONS AND js_D ARRAY FOR VELOCITY

       !===ASSEMBLING NS RESIDUAL MATRICES
       ALLOCATE(LES_mat(2),LES_ksp(2))
       ALLOCATE(zero_dir1(SIZE(vv_js_D(1)%DIL)))
       ALLOCATE(zero_dir2(SIZE(vv_js_D(2)%DIL)))
       ALLOCATE(zero_dir3(SIZE(vv_js_D(3)%DIL)))

       DO k = 1, 2
          nu_mat = k
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA,LES_mat(nu_mat), clean=.FALSE.)
          CALL qs_mass_vect_3x3(vv_3_LA, vv_mesh, 1.d0, LES_mat(nu_mat))
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
                  LES_mat(nu_mat), vv_3_LA)
          END IF
          CALL init_solver(inputs%my_par_vv,LES_ksp(nu_mat),LES_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
       END DO
       iteration = 0
    END IF !end of once

    iteration = iteration + 1

    !===Compute Strain rate tensor
    CALL smb_explicit_strain_rate_tensor(vv_mesh, list_mode, un_m1, strain_rate_tensor)
    strain_rate_tensor = 1/inputs%Re*strain_rate_tensor
    CALL smb_explicit_strain_rate_tensor_bdy(vv_mesh, list_mode, un_m1, strain_rate_tensor_scal_n_bdy)
    strain_rate_tensor_scal_n_bdy = -1/inputs%Re*strain_rate_tensor_scal_n_bdy

    !===Computation of rhs at Gauss points for every mode without the strain rate tensor
    CALL rhs_residual_ns_gauss_3x3(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time-inputs%dt, &
         (un-un_m2)/(2*inputs%dt), pn_m1, rotv_v_m1, rhs_gauss, density, tempn, concn)
    !===End Computation of rhs

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236,&
            opt_tensor=strain_rate_tensor(:,:,:,i), opt_tensor_scaln_bdy=strain_rate_tensor_scal_n_bdy(:,:,i))

       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
       END IF
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss


       !===Solve linear system for momentum equation
       !Solve system 1, ur_c, ut_s, uz_c
       nu_mat  = 1
       CALL solver(LES_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,1,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,4,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,5,i))

       !Solve system 2, ur_s, ut_c, uz_s
       nu_mat = 2
       CALL solver(LES_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,2,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,3,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,6,i))
       !===End Solve linear system for momentum equation
    END DO

    !norm_res_ns_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, res_ns)
    !CALL vtu_3d(res_ns, 'vv_mesh', 'Resns', 'resns', 'new')
    !WRITE(*,*) 'norm_res_ns_L2 = ', norm_res_ns_L2
!!$    IF (MOD(iteration,inputs%freq_en) == 0) THEN
!!$       CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!!$       norm_res_ns_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, res_ns)
!!$       !CALL vtu_3d(res_ns, 'vv_mesh', 'Resns', 'resns', 'new')
!!$       IF (rank==0) THEN
!!$          WRITE(*,*) 'norm_res_ns_L2 = ', norm_res_ns_L2
!!$       END IF
!!$    END IF

    !===Compute un_m1 and res_ns on gauss points
    DO i = 1, SIZE(list_mode)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO TYPE = 1, 6
                un_m1_gauss(index,TYPE,i)  = SUM(un_m1(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                res_ns_gauss(index,TYPE,i) = SUM(res_ns(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
             END DO
          END DO
       END DO
    END DO
    !===End compute un_m1 and res_ns on gauss points

    !===Compute entropy viscosity times gradient of 2*un-un_m1 in real space by FFT
    CALL smb_explicit_grad_vel_LES(vv_mesh, list_mode, 2*un-un_m1, strain_rate_tensor)

    norm_vel_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_m1)

    CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
    bloc_size = SIZE(un_m1_gauss,1)/nb_procs+1
    bloc_size = vv_mesh%gauss%l_G*(bloc_size/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    CALL FFT_COMPUTE_ENTROPY_VISC(comm_one_d(2), comm_one_d(1), un_m1_gauss, res_ns_gauss, &
         strain_rate_tensor(1,:,:,:), strain_rate_tensor(2,:,:,:), strain_rate_tensor(3,:,:,:), &
         vv_mesh%hloc_gauss, visco_entro_grad_u(1,:,:,:), visco_entro_grad_u(2,:,:,:), &
         visco_entro_grad_u(3,:,:,:), nb_procs, bloc_size, m_max_pad, norm_vel_L2**2/Volume_3D, vv_mesh%gauss%l_G)
    !===End compute entropy viscosity times gradient of 2*un-un_m1 in real space by FFT

  END SUBROUTINE compute_entropy_viscosity

  SUBROUTINE compute_entropy_viscosity_mom(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, &
       momentum, momentum_m1, momentum_m2, pn_m1, un_m1, tensor_m1, visc_grad_vel_m1, tensor_surface_gauss, &
       rotb_b_m1, visco_dyn_m1, density_m2, density_m1, density, tempn, concn, visc_entro_real, visc_entro_level_real)
    USE def_type_mesh
    USE fem_M_axi
    USE solve_petsc
    USE periodic
    USE rhs_gauss_computing
    USE rhs_para_assembling
    USE Dir_nodes_petsc
    USE st_matrix
    USE sft_parallele
    USE tn_axi
    USE input_data
    USE my_util
    USE subroutine_mass
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                    :: vv_3_LA
    TYPE(mesh_type),                  INTENT(IN)          :: pp_mesh, vv_mesh
    REAL(KIND=8),                     INTENT(IN)          :: time
    INTEGER,      DIMENSION(:),       INTENT(IN)          :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: momentum, momentum_m1, momentum_m2
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: un_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: pn_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)          :: tensor_m1, visc_grad_vel_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)          :: tensor_surface_gauss !t=tn_m1 in bdf1 / tn if bdf2)
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: rotb_b_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: visco_dyn_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: density_m2, density_m1, density
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: concn
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)         :: visc_entro_real
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)         :: visc_entro_level_real
    !TYPE(dyn_int_line), DIMENSION(3),                SAVE :: vv_js_D
    LOGICAL,                                         SAVE :: once = .TRUE.
    REAL(KIND=8),                                    SAVE :: Volume_3D
    INTEGER,                                         SAVE :: iteration
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: tensor_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_Gs*vv_mesh%dom_mes,6,SIZE(list_mode))   :: visc_grad_vel_m1_scal_n_bdy
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_Gs*vv_mesh%dom_mes,6,SIZE(list_mode))   :: tensor_m1_scal_n_bdy
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: un_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: momentum_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: res_ns_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                           :: res_ns
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,2,SIZE(list_mode))     :: res_mass_gauss
    REAL(KIND=8), DIMENSION(pp_mesh%np,6,SIZE(list_mode))                           :: un_m1_P1, momentum_m1_P1, res_ns_P1
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode))                           :: density_P1, density_m2_P1
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%l_G*pp_mesh%dom_me,6,SIZE(list_mode))     :: un_m1_P1_gauss
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%l_G*pp_mesh%dom_me,6,SIZE(list_mode))     :: momentum_m1_P1_gauss
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%l_G*pp_mesh%dom_me,6,SIZE(list_mode))     :: res_ns_P1_gauss
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%l_G*pp_mesh%dom_me,2,SIZE(list_mode))     :: res_mass_P1_gauss
    INTEGER,          POINTER, DIMENSION(:)                     :: vv_3_ifrom
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                  :: j_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_w)                  :: j_loc_P1
    REAL(KIND=8)    :: norm_vel_L2, norm_mom_L2
    !REAL(KIND=8)    :: norm_res_ns_L2
    !INTEGER         :: rank
    INTEGER         :: i, k, nu_mat, mode, m, l, TYPE, index, n
    INTEGER         :: nb_procs, code, bloc_size, m_max_pad
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Mat, DIMENSION(:), POINTER, SAVE :: LES_mat
    KSP, DIMENSION(:), POINTER, SAVE :: LES_ksp

    IF (once) THEN

       once = .FALSE.

       !===CREATE PETSC VECTORS AND GHOSTS
       CALL create_my_ghost(vv_mesh,vv_3_LA,vv_3_ifrom)
       n = 3*vv_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
       CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
       CALL VecDuplicate(vx_3, vb_3_145, ierr)
       CALL VecDuplicate(vx_3, vb_3_236, ierr)

       !===Compute Volume
       CALL twoD_volume(comm_one_d(1),vv_mesh,Volume_3D)
       Volume_3D = Volume_3D*2*ACOS(-1.d0)

       !===ASSEMBLING NS RESIDUAL MATRICES
       ALLOCATE(LES_mat(2),LES_ksp(2))
       DO k = 1, 2
          nu_mat = k
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA,LES_mat(nu_mat), clean=.FALSE.)
          CALL qs_mass_vect_3x3(vv_3_LA, vv_mesh, 1.d0, LES_mat(nu_mat))
!!$          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
!!$             CALL periodic_matrix_petsc(vvz_per%n_bord, vvz_per%list,vvz_per%perlist, &
!!$                  LES_mat(nu_mat), vv_3_LA)
!!$          END IF
          CALL init_solver(inputs%my_par_vv,LES_ksp(nu_mat),LES_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
       END DO

       iteration = 0
    END IF !end of once

    iteration = iteration + 1

    !===Compute Strain rate tensor scalar normal on bdy
    CALL smb_explicit_strain_rate_tensor_bdy_mom(comm_one_d(2), vv_mesh, list_mode, visco_dyn_m1, un_m1, &
         visc_grad_vel_m1_scal_n_bdy)
    visc_grad_vel_m1_scal_n_bdy = -1/inputs%Re*visc_grad_vel_m1_scal_n_bdy

    !===Compute tensor scalar normal on bdy
    CALL smb_explicit_tensor_bdy(vv_mesh, list_mode, tensor_m1, tensor_m1_scal_n_bdy)

    !===Computation of rhs at Gauss points for every mode without tensors
    CALL rhs_residual_ns_gauss_3x3_mom(vv_mesh, pp_mesh, list_mode, time-inputs%dt, &
         (momentum-momentum_m2)/(2*inputs%dt), pn_m1, density_m1, rotb_b_m1, rhs_gauss, tempn, concn)
    !===End Computation of rhs

    !===Computation of tensor_m1(=m x u) on gauss points
    DO i = 1, SIZE(list_mode)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO TYPE = 1, 6
                DO k = 1, 3
                   tensor_m1_gauss(k,index,TYPE,i)=SUM(tensor_m1(k,j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                END DO
             END DO
          END DO
       END DO
    END DO
    !===End computation of tensor_m1(=m x u) on gauss points

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236,&
            opt_tensor=-tensor_m1_gauss(:,:,:,i)+visc_grad_vel_m1(:,:,:,i)+tensor_surface_gauss(:,:,:,i), &
            opt_tensor_scaln_bdy=visc_grad_vel_m1_scal_n_bdy(:,:,i)+tensor_m1_scal_n_bdy(:,:,i))

!!$       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
!!$          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_236, vv_3_LA)
!!$          CALL periodic_rhs_petsc(vvz_per%n_bord, vvz_per%list, vvz_per%perlist, vb_3_145, vv_3_LA)
!!$       END IF
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss


       !===Solve linear system for momentum equation
       !Solve system 1, ur_c, ut_s, uz_c
       nu_mat  = 1
       CALL solver(LES_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,1,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,4,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,5,i))

       !Solve system 2, ur_s, ut_c, uz_s
       nu_mat = 2
       CALL solver(LES_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,2,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,3,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,6,i))
       !===End Solve linear system for momentum equation
    END DO

!!$    IF (MOD(iteration,inputs%freq_en) == 0) THEN
!!$       CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!!$       norm_res_ns_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, res_ns)
!!$       !CALL vtu_3d(res_ns, 'vv_mesh', 'Resns', 'resns', 'new')
!!$       IF (rank==0) THEN
!!$          WRITE(*,*) 'norm_res_ns_L2 = ', norm_res_ns_L2
!!$       END IF
!!$    END IF

    !===Compute un_m1 and res_ns on gauss points
    DO i = 1, SIZE(list_mode)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO TYPE = 1, 6
                un_m1_gauss(index,TYPE,i)  = SUM(un_m1(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                momentum_m1_gauss(index,TYPE,i)  = SUM(momentum_m1(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                res_ns_gauss(index,TYPE,i) = SUM(res_ns(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
             END DO
          END DO
       END DO
    END DO
    !===End compute un_m1 and res_ns on gauss points

    IF (.NOT.inputs%if_LES_in_momentum) THEN
       res_ns_gauss=0.d0
    END IF

    !===Compute res_mass on gauss points
    CALL  compute_res_mass_gauss(vv_mesh, list_mode, density_m2, density, momentum_m1, res_mass_gauss)
    !===End compute res_mass on gauss points

    !===Compute entropy viscosity for momentum and level set equations on gauss points
    norm_vel_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_m1)
    norm_mom_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, momentum_m1)

    CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
    bloc_size = SIZE(un_m1_gauss,1)/nb_procs+1
    bloc_size = vv_mesh%gauss%l_G*(bloc_size/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    IF (inputs%if_level_set_P2) THEN
       !===Entropy viscosity for momentum
       CALL FFT_COMPUTE_ENTROPY_VISC_MOM(comm_one_d(2), comm_one_d(1), un_m1_gauss, momentum_m1_gauss, res_ns_gauss, &
            res_mass_gauss, vv_mesh%hloc_gauss, visc_entro_real, nb_procs, bloc_size, m_max_pad,&
            vv_mesh%gauss%l_G, opt_c2_real_out=visc_entro_level_real)
    ELSE
       !===Entropy viscosity for momentum
       CALL FFT_COMPUTE_ENTROPY_VISC_MOM(comm_one_d(2), comm_one_d(1), un_m1_gauss, momentum_m1_gauss, res_ns_gauss, &
            res_mass_gauss, vv_mesh%hloc_gauss, visc_entro_real, nb_procs, bloc_size, m_max_pad,&
            vv_mesh%gauss%l_G)

       !===Compute un_m1, momentum_m1, res_ns, density and density_m2 on P1 nodes
       DO i = 1, SIZE(list_mode)
          DO k = 1, 6
             CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, un_m1(:,k,i), un_m1_P1(:,k,i))
             CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, momentum_m1(:,k,i), momentum_m1_P1(:,k,i))
             CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, res_ns(:,k,i), res_ns_P1(:,k,i))
          END DO
          DO k = 1, 2
             CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, density(:,k,i), density_P1(:,k,i))
             CALL project_P2_P1(vv_mesh%jj, pp_mesh%jj, density_m2(:,k,i), density_m2_P1(:,k,i))
          END DO
       END DO
       !===End compute un_m1, momentum_m1, res_ns, density and density_m2 on P1 nodes

       !===Compute un_m1 and res_ns on P1 gauss points
       DO i = 1, SIZE(list_mode)
          index = 0
          DO m = 1, pp_mesh%dom_me
             j_loc_P1 = pp_mesh%jj(:,m)
             DO l = 1, pp_mesh%gauss%l_G
                index = index + 1
                DO TYPE = 1, 6
                   un_m1_P1_gauss(index,TYPE,i)  = SUM(un_m1_P1(j_loc_P1,TYPE,i)*pp_mesh%gauss%ww(:,l))
                   momentum_m1_P1_gauss(index,TYPE,i)  = SUM(momentum_m1_P1(j_loc_P1,TYPE,i)*pp_mesh%gauss%ww(:,l))
                   res_ns_P1_gauss(index,TYPE,i) = SUM(res_ns_P1(j_loc_P1,TYPE,i)*pp_mesh%gauss%ww(:,l))
                END DO
             END DO
          END DO
       END DO
       !===End compute un_m1 and res_ns on P1 gauss points

       IF (.NOT.inputs%if_LES_in_momentum) THEN
          res_ns_P1_gauss=0.d0
       END IF

       !===Compute res_mass on P1 gauss points
       CALL  compute_res_mass_gauss(pp_mesh, list_mode, density_m2_P1, density_P1, momentum_m1_P1, res_mass_P1_gauss)
       !===End compute res_mass on P1 gauss points

       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       bloc_size = SIZE(un_m1_P1_gauss,1)/nb_procs+1
       bloc_size = pp_mesh%gauss%l_G*(bloc_size/pp_mesh%gauss%l_G)+pp_mesh%gauss%l_G
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2

       !===Entropy viscosity for level set
       CALL FFT_COMPUTE_ENTROPY_VISC_MOM(comm_one_d(2), comm_one_d(1), un_m1_P1_gauss, momentum_m1_P1_gauss,&
            res_ns_P1_gauss, res_mass_P1_gauss, pp_mesh%hloc_gauss, visc_entro_level_real, nb_procs, bloc_size, m_max_pad,&
            pp_mesh%gauss%l_G)
    END IF

!!$    IF (MOD(iteration,inputs%freq_en) == 0) THEN
!!$       CALL MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr)
!!$       norm_res_ns_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, visc_entro)
!!$       IF (rank==0) THEN
!!$          WRITE(*,*) 'norm_visc_entro_L2 = ', norm_res_ns_L2
!!$       END IF
!!$       norm_res_ns_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, visc_entro_level)
!!$       IF (rank==0) THEN
!!$          WRITE(*,*) 'norm_visc_entro_level_L2 = ', norm_res_ns_L2
!!$       END IF
!!$    END IF

  END SUBROUTINE compute_entropy_viscosity_mom

  SUBROUTINE compute_entropy_viscosity_mom_no_level_set(comm_one_d, vv_3_LA, vv_mesh, pp_mesh, time, list_mode, &
       momentum, momentum_m1, momentum_m2, pn_m1, un_m1, tensor_m1, &
       rotb_b_m1, density, tempn, concn, visc_entro_real)
    USE def_type_mesh
    USE fem_M_axi
    USE solve_petsc
    USE periodic
    USE rhs_gauss_computing
    USE rhs_para_assembling
    USE Dir_nodes_petsc
    USE st_matrix
    USE sft_parallele
    USE tn_axi
    USE input_data
    USE my_util
    USE subroutine_mass
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                    :: vv_3_LA
    TYPE(mesh_type),                  INTENT(IN)          :: pp_mesh, vv_mesh
    REAL(KIND=8),                     INTENT(IN)          :: time
    INTEGER,      DIMENSION(:),       INTENT(IN)          :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: momentum, momentum_m1, momentum_m2
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: un_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: pn_m1
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)          :: tensor_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: rotb_b_m1
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: concn
    REAL(KIND=8), DIMENSION(:,:,:),   INTENT(IN)          :: density
    REAL(KIND=8), DIMENSION(:,:),     INTENT(OUT)         :: visc_entro_real
    !TYPE(dyn_int_line), DIMENSION(3),                SAVE :: vv_js_D
    LOGICAL,                                         SAVE :: once = .TRUE.
    REAL(KIND=8),                                    SAVE :: Volume_3D
    INTEGER,                                         SAVE :: iteration
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: tensor_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_Gs*vv_mesh%dom_mes,6,SIZE(list_mode))   :: tensor_m1_scal_n_bdy
    REAL(KIND=8), DIMENSION(3,vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))   :: visc_grad_vel_m1
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_Gs*vv_mesh%dom_mes,6,SIZE(list_mode))   :: visc_grad_vel_m1_scal_n_bdy
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: rhs_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,2,SIZE(list_mode))     :: res_mass_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: un_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: momentum_m1_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode))     :: res_ns_gauss
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode))                           :: res_ns
    INTEGER,          POINTER, DIMENSION(:)                     :: vv_3_ifrom
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8)    :: norm_vel_L2, norm_mom_L2
    !REAL(KIND=8)    :: norm_res_ns_L2
    !INTEGER         :: rank
    INTEGER         :: i, k, nu_mat, mode, m, l, TYPE, index, n
    INTEGER         :: nb_procs, code, bloc_size, m_max_pad
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Vec,                        SAVE :: vb_3_145, vb_3_236, vx_3, vx_3_ghost
    Mat, DIMENSION(:), POINTER, SAVE :: LES_mat
    KSP, DIMENSION(:), POINTER, SAVE :: LES_ksp

    IF (once) THEN

       once = .FALSE.

       !===CREATE PETSC VECTORS AND GHOSTS
       CALL create_my_ghost(vv_mesh,vv_3_LA,vv_3_ifrom)
       n = 3*vv_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(vv_3_ifrom), vv_3_ifrom, vx_3, ierr)
       CALL VecGhostGetLocalForm(vx_3, vx_3_ghost, ierr)
       CALL VecDuplicate(vx_3, vb_3_145, ierr)
       CALL VecDuplicate(vx_3, vb_3_236, ierr)

       !===Compute Volume
       CALL twoD_volume(comm_one_d(1),vv_mesh,Volume_3D)
       Volume_3D = Volume_3D*2*ACOS(-1.d0)

       !===ASSEMBLING NS RESIDUAL MATRICES
       ALLOCATE(LES_mat(2),LES_ksp(2))
       DO k = 1, 2
          nu_mat = k
          CALL create_local_petsc_matrix(comm_one_d(1), vv_3_LA,LES_mat(nu_mat), clean=.FALSE.)
          CALL qs_mass_vect_3x3(vv_3_LA, vv_mesh, 1.d0, LES_mat(nu_mat))
          CALL init_solver(inputs%my_par_vv,LES_ksp(nu_mat),LES_mat(nu_mat),comm_one_d(1),&
               solver=inputs%my_par_vv%solver,precond=inputs%my_par_vv%precond)
       END DO

       iteration = 0
    END IF !end of once

    iteration = iteration + 1

    !===Compute Strain rate tensor
    CALL smb_explicit_strain_rate_tensor(vv_mesh, list_mode, un_m1, visc_grad_vel_m1)
    visc_grad_vel_m1 = 1/inputs%Re* visc_grad_vel_m1
    CALL smb_explicit_strain_rate_tensor_bdy(vv_mesh, list_mode, un_m1, visc_grad_vel_m1_scal_n_bdy)
    visc_grad_vel_m1_scal_n_bdy = -1/inputs%Re*visc_grad_vel_m1_scal_n_bdy
    !===End Compute Strain rate tensor

    !===Compute tensor scalar normal on bdy
    CALL smb_explicit_tensor_bdy(vv_mesh, list_mode, tensor_m1, tensor_m1_scal_n_bdy)

    !===Computation of rhs at Gauss points for every mode without tensors
    CALL rhs_residual_ns_gauss_3x3(vv_mesh, pp_mesh, comm_one_d(2), list_mode, time-inputs%dt, &
         (momentum-momentum_m2)/(2*inputs%dt), pn_m1, -rotb_b_m1, rhs_gauss, density, tempn, concn)
    !===End Computation of rhs

    !===Computation of tensor_m1(=m x u) on gauss points
    DO i = 1, SIZE(list_mode)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO TYPE = 1, 6
                DO k = 1, 3
                   tensor_m1_gauss(k,index,TYPE,i)=SUM(tensor_m1(k,j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                END DO
             END DO
          END DO
       END DO
    END DO
    !===End computation of tensor_m1(=m x u) on gauss points

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)

       !===Assemble vb_3_145, vb_3_236 using rhs_gauss
       CALL rhs_3x3(vv_mesh, vv_3_LA, mode, rhs_gauss(:,:,i), vb_3_145, vb_3_236,&
            opt_tensor=-tensor_m1_gauss(:,:,:,i)+visc_grad_vel_m1(:,:,:,i), &
            opt_tensor_scaln_bdy=visc_grad_vel_m1_scal_n_bdy(:,:,i)+tensor_m1_scal_n_bdy(:,:,i))
       !===End Assemble vb_3_145, vb_3_236 using rhs_gauss


       !===Solve linear system for momentum equation
       !Solve system 1, ur_c, ut_s, uz_c
       nu_mat  = 1
       CALL solver(LES_ksp(nu_mat),vb_3_145,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,1,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,4,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,5,i))

       !Solve system 2, ur_s, ut_c, uz_s
       nu_mat = 2
       CALL solver(LES_ksp(nu_mat),vb_3_236,vx_3,reinit=.FALSE.,verbose=inputs%my_par_vv%verbose)
       CALL VecGhostUpdateBegin(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(vx_3,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(vx_3_ghost,1,1,vv_3_LA,res_ns(:,2,i))
       CALL extract(vx_3_ghost,2,2,vv_3_LA,res_ns(:,3,i))
       CALL extract(vx_3_ghost,3,3,vv_3_LA,res_ns(:,6,i))
       !===End Solve linear system for momentum equation
    END DO

    !===Compute un_m1 and res_ns on gauss points
    DO i = 1, SIZE(list_mode)
       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)
          DO l = 1, vv_mesh%gauss%l_G
             index = index + 1
             DO TYPE = 1, 6
                un_m1_gauss(index,TYPE,i)  = SUM(un_m1(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                momentum_m1_gauss(index,TYPE,i)  = SUM(momentum_m1(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
                res_ns_gauss(index,TYPE,i) = SUM(res_ns(j_loc,TYPE,i)*vv_mesh%gauss%ww(:,l))
             END DO
          END DO
       END DO
    END DO
    !===End compute un_m1 and res_ns on gauss points

    !===Compute entropy viscosity for momentum and level set equations on gauss points
    norm_vel_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_m1)
    norm_mom_L2 = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, momentum_m1)

    CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
    bloc_size = SIZE(un_m1_gauss,1)/nb_procs+1
    bloc_size = vv_mesh%gauss%l_G*(bloc_size/vv_mesh%gauss%l_G)+vv_mesh%gauss%l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===Entropy viscosity for momentum
    res_mass_gauss=0.d0
    CALL FFT_COMPUTE_ENTROPY_VISC_MOM(comm_one_d(2), comm_one_d(1), un_m1_gauss, momentum_m1_gauss, res_ns_gauss, &
         res_mass_gauss, vv_mesh%hloc_gauss, visc_entro_real, nb_procs, bloc_size, m_max_pad,&
         vv_mesh%gauss%l_G)

  END SUBROUTINE compute_entropy_viscosity_mom_no_level_set

  SUBROUTINE smb_explicit_grad_vel_LES(mesh, list_mode, vel, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))                :: grad1_vel, grad2_vel, grad3_vel
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8)                                :: ray, hh, hm
    INTEGER                                     :: m, l , i, mode, index, k

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vel_loc(:,k) = vel(j_loc,k,i)
          END DO
          DO l = 1, mesh%gauss%l_G
             index = index + 1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !===Compute local mesh sizes
             hh=mesh%hloc_gauss(index)
             hm=MIN(mesh%hm(i),hh) !(JLG April 7 2017)
             !hm=0.5d0/inputs%m_max

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

             !-----------------Grad u_r on Gauss points------------------------------------
             grad1_vel(index,1,i) = SUM(vel_loc(:,1)*dw_loc(1,:))*hh
             grad1_vel(index,2,i) = SUM(vel_loc(:,2)*dw_loc(1,:))*hh
             grad1_vel(index,3,i) =  (mode*SUM(vel_loc(:,2)*mesh%gauss%ww(:,l)) - &
                  SUM(vel_loc(:,3)*mesh%gauss%ww(:,l)))/ray*hm
             grad1_vel(index,4,i) = (-mode*SUM(vel_loc(:,1)*mesh%gauss%ww(:,l)) - &
                  SUM(vel_loc(:,4)*mesh%gauss%ww(:,l)))/ray*hm
             grad1_vel(index,5,i) =  SUM(vel_loc(:,1)*dw_loc(2,:))*hh
             grad1_vel(index,6,i) =  SUM(vel_loc(:,2)*dw_loc(2,:))*hh

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(index,1,i) = SUM(vel_loc(:,3)*dw_loc(1,:))*hh
             grad2_vel(index,2,i) = SUM(vel_loc(:,4)*dw_loc(1,:))*hh
             grad2_vel(index,3,i) =  (mode*SUM(vel_loc(:,4)*mesh%gauss%ww(:,l)) + &
                  SUM(vel_loc(:,1)*mesh%gauss%ww(:,l)))/ray*hm
             grad2_vel(index,4,i) = (-mode*SUM(vel_loc(:,3)*mesh%gauss%ww(:,l)) + &
                  SUM(vel_loc(:,2)*mesh%gauss%ww(:,l)))/ray*hm
             grad2_vel(index,5,i) = SUM(vel_loc(:,3)*dw_loc(2,:))*hh
             grad2_vel(index,6,i) = SUM(vel_loc(:,4)*dw_loc(2,:))*hh

             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(index,1,i) = SUM(vel_loc(:,5)*dw_loc(1,:))*hh
             grad3_vel(index,2,i) = SUM(vel_loc(:,6)*dw_loc(1,:))*hh
             grad3_vel(index,3,i) = mode*SUM(vel_loc(:,6)*mesh%gauss%ww(:,l))/ray*hm
             grad3_vel(index,4,i) = -mode*SUM(vel_loc(:,5)*mesh%gauss%ww(:,l))/ray*hm
             grad3_vel(index,5,i) = SUM(vel_loc(:,5)*dw_loc(2,:))*hh
             grad3_vel(index,6,i) = SUM(vel_loc(:,6)*dw_loc(2,:))*hh
          ENDDO
       ENDDO
    END DO

    V_out(1,:,:,:) = grad1_vel
    V_out(2,:,:,:) = grad2_vel
    V_out(3,:,:,:) = grad3_vel

  END SUBROUTINE smb_explicit_grad_vel_LES

  SUBROUTINE smb_explicit_strain_rate_tensor(mesh, list_mode, vel, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(3,mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))                :: grad1_vel, grad2_vel, grad3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode))                :: gradT1_vel, gradT2_vel, gradT3_vel
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, l , i, mode, index, k

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vel_loc(:,k) = vel(j_loc,k,i)
          END DO
          DO l = 1, mesh%gauss%l_G
             index = index + 1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

             !-----------------Grad u_r on Gauss points------------------------------------
             grad1_vel(index,1,i) = SUM(vel_loc(:,1)*dw_loc(1,:))
             grad1_vel(index,2,i) = SUM(vel_loc(:,2)*dw_loc(1,:))
             grad1_vel(index,3,i) =  (mode*SUM(vel_loc(:,2)*mesh%gauss%ww(:,l)) - &
                  SUM(vel_loc(:,3)*mesh%gauss%ww(:,l)))/ray
             grad1_vel(index,4,i) = (-mode*SUM(vel_loc(:,1)*mesh%gauss%ww(:,l)) - &
                  SUM(vel_loc(:,4)*mesh%gauss%ww(:,l)))/ray
             grad1_vel(index,5,i) =  SUM(vel_loc(:,1)*dw_loc(2,:))
             grad1_vel(index,6,i) =  SUM(vel_loc(:,2)*dw_loc(2,:))

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(index,1,i) = SUM(vel_loc(:,3)*dw_loc(1,:))
             grad2_vel(index,2,i) = SUM(vel_loc(:,4)*dw_loc(1,:))
             grad2_vel(index,3,i) =  (mode*SUM(vel_loc(:,4)*mesh%gauss%ww(:,l)) + &
                  SUM(vel_loc(:,1)*mesh%gauss%ww(:,l)))/ray
             grad2_vel(index,4,i) = (-mode*SUM(vel_loc(:,3)*mesh%gauss%ww(:,l)) + &
                  SUM(vel_loc(:,2)*mesh%gauss%ww(:,l)))/ray
             grad2_vel(index,5,i) = SUM(vel_loc(:,3)*dw_loc(2,:))
             grad2_vel(index,6,i) = SUM(vel_loc(:,4)*dw_loc(2,:))

             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(index,1,i) = SUM(vel_loc(:,5)*dw_loc(1,:))
             grad3_vel(index,2,i) = SUM(vel_loc(:,6)*dw_loc(1,:))
             grad3_vel(index,3,i) = mode*SUM(vel_loc(:,6)*mesh%gauss%ww(:,l))/ray
             grad3_vel(index,4,i) = -mode*SUM(vel_loc(:,5)*mesh%gauss%ww(:,l))/ray
             grad3_vel(index,5,i) = SUM(vel_loc(:,5)*dw_loc(2,:))
             grad3_vel(index,6,i) = SUM(vel_loc(:,6)*dw_loc(2,:))
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
    V_out(1,:,:,:) = grad1_vel + gradT1_vel
    V_out(2,:,:,:) = grad2_vel + gradT2_vel
    V_out(3,:,:,:) = grad3_vel + gradT3_vel

  END SUBROUTINE smb_explicit_strain_rate_tensor

  SUBROUTINE smb_explicit_strain_rate_tensor_bdy(mesh, list_mode, vel, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode))              :: grad1_vel, grad2_vel, grad3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode))              :: gradT1_vel, gradT2_vel, gradT3_vel
    INTEGER,      DIMENSION(mesh%gauss%n_ws)                 :: js_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dws_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6)               :: vels_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, ms, ls , i, mode, indexs, k

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       indexs = 0
       DO ms = 1, mesh%dom_mes
          js_loc = mesh%jjs(:,ms)
          m = mesh%neighs(ms)
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vels_loc(:,k) = vel(js_loc,k,i)
             vel_loc(:,k)  = vel(j_loc,k,i)
          END DO
          DO ls = 1, mesh%gauss%l_Gs
             indexs = indexs + 1
             dws_loc = mesh%gauss%dw_s(:,:,ls,ms)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,js_loc)*mesh%gauss%wws(:,ls))

             !===Dont compute on r = 0
             IF (ray.LT.1.d-10) THEN
                V_out(indexs,:,i) = 0.d0
                CYCLE
             END IF

             !-----------------Grad u_r on Gauss points------------------------------------
             grad1_vel(indexs,1,i) = SUM(vel_loc(:,1)*dws_loc(1,:))
             grad1_vel(indexs,2,i) = SUM(vel_loc(:,2)*dws_loc(1,:))
             grad1_vel(indexs,3,i) =  (mode*SUM(vels_loc(:,2)*mesh%gauss%wws(:,ls)) - &
                  SUM(vels_loc(:,3)*mesh%gauss%wws(:,ls)))/ray
             grad1_vel(indexs,4,i) = (-mode*SUM(vels_loc(:,1)*mesh%gauss%wws(:,ls)) - &
                  SUM(vels_loc(:,4)*mesh%gauss%wws(:,ls)))/ray
             grad1_vel(indexs,5,i) =  SUM(vel_loc(:,1)*dws_loc(2,:))
             grad1_vel(indexs,6,i) =  SUM(vel_loc(:,2)*dws_loc(2,:))

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(indexs,1,i) = SUM(vel_loc(:,3)*dws_loc(1,:))
             grad2_vel(indexs,2,i) = SUM(vel_loc(:,4)*dws_loc(1,:))
             grad2_vel(indexs,3,i) =  (mode*SUM(vels_loc(:,4)*mesh%gauss%wws(:,ls)) + &
                  SUM(vels_loc(:,1)*mesh%gauss%wws(:,ls)))/ray
             grad2_vel(indexs,4,i) = (-mode*SUM(vels_loc(:,3)*mesh%gauss%wws(:,ls)) + &
                  SUM(vels_loc(:,2)*mesh%gauss%wws(:,ls)))/ray
             grad2_vel(indexs,5,i) = SUM(vel_loc(:,3)*dws_loc(2,:))
             grad2_vel(indexs,6,i) = SUM(vel_loc(:,4)*dws_loc(2,:))

             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(indexs,1,i) = SUM(vel_loc(:,5)*dws_loc(1,:))
             grad3_vel(indexs,2,i) = SUM(vel_loc(:,6)*dws_loc(1,:))
             grad3_vel(indexs,3,i) =  mode*SUM(vels_loc(:,6)*mesh%gauss%wws(:,ls))/ray
             grad3_vel(indexs,4,i) = -mode*SUM(vels_loc(:,5)*mesh%gauss%wws(:,ls))/ray
             grad3_vel(indexs,5,i) = SUM(vel_loc(:,5)*dws_loc(2,:))
             grad3_vel(indexs,6,i) = SUM(vel_loc(:,6)*dws_loc(2,:))

             !-----------------GradT u_r on Gauss points------------------------------------
             gradT1_vel(indexs,1,i) = grad1_vel(indexs,1,i)
             gradT1_vel(indexs,2,i) = grad1_vel(indexs,2,i)
             gradT1_vel(indexs,3,i) = grad2_vel(indexs,1,i)
             gradT1_vel(indexs,4,i) = grad2_vel(indexs,2,i)
             gradT1_vel(indexs,5,i) = grad3_vel(indexs,1,i)
             gradT1_vel(indexs,6,i) = grad3_vel(indexs,2,i)
             !-----------------GradT u_th on Gauss points-----------------------------------
             gradT2_vel(indexs,1,i) = grad1_vel(indexs,3,i)
             gradT2_vel(indexs,2,i) = grad1_vel(indexs,4,i)
             gradT2_vel(indexs,3,i) = grad2_vel(indexs,3,i)
             gradT2_vel(indexs,4,i) = grad2_vel(indexs,4,i)
             gradT2_vel(indexs,5,i) = grad3_vel(indexs,3,i)
             gradT2_vel(indexs,6,i) = grad3_vel(indexs,4,i)
             !-----------------GradT u_z on Gauss points------------------------------------
             gradT3_vel(indexs,1,i) = grad1_vel(indexs,5,i)
             gradT3_vel(indexs,2,i) = grad1_vel(indexs,6,i)
             gradT3_vel(indexs,3,i) = grad2_vel(indexs,5,i)
             gradT3_vel(indexs,4,i) = grad2_vel(indexs,6,i)
             gradT3_vel(indexs,5,i) = grad3_vel(indexs,5,i)
             gradT3_vel(indexs,6,i) = grad3_vel(indexs,6,i)

             !-----------------Grad_sym_u scalar normal-------------------------------------
             V_out(indexs,1,i) = (grad1_vel(indexs,1,i)+gradT1_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad1_vel(indexs,5,i)+gradT1_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,2,i) = (grad1_vel(indexs,2,i)+gradT1_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad1_vel(indexs,6,i)+gradT1_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,3,i) = (grad2_vel(indexs,1,i)+gradT2_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad2_vel(indexs,5,i)+gradT2_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,4,i) = (grad2_vel(indexs,2,i)+gradT2_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad2_vel(indexs,6,i)+gradT2_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,5,i) = (grad3_vel(indexs,1,i)+gradT3_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad3_vel(indexs,5,i)+gradT3_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,6,i) = (grad3_vel(indexs,2,i)+gradT3_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad3_vel(indexs,6,i)+gradT3_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
          ENDDO !l_Gs
       ENDDO !mes
    END DO !i

  END SUBROUTINE smb_explicit_strain_rate_tensor_bdy

  SUBROUTINE smb_explicit_strain_rate_tensor_bdy_mom(communicator, mesh, list_mode, visc_dyn, vel, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: visc_dyn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode))    :: grad1_vel, grad2_vel, grad3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode))    :: gradT1_vel, gradT2_vel, gradT3_vel
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode))    :: V_bdy
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,2,SIZE(list_mode))    :: visc_dyn_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_ws)                 :: js_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dws_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: vel_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6)               :: vels_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2)               :: visc_dyns_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, ms, ls , i, mode, indexs, k
    INTEGER                                     :: code, bloc_size, m_max_pad, nb_procs
    MPI_Comm          :: communicator

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       indexs = 0
       DO ms = 1, mesh%dom_mes
          js_loc = mesh%jjs(:,ms)
          m = mesh%neighs(ms)
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             vels_loc(:,k) = vel(js_loc,k,i)
             vel_loc(:,k)  = vel(j_loc,k,i)
          END DO
          DO k = 1, 2
             visc_dyns_loc(:,k) = visc_dyn(js_loc,k,i)
          END DO
          DO ls = 1, mesh%gauss%l_Gs
             indexs = indexs + 1
             dws_loc = mesh%gauss%dw_s(:,:,ls,ms)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,js_loc)*mesh%gauss%wws(:,ls))

             !===Dont compute on r = 0
             IF (ray.LT.1.d-10) THEN
                visc_dyn_gauss(indexs,:,i)=0.d0
                V_bdy(indexs,:,i) =0.d0
                CYCLE
             END IF

             !-----------------Visc_dyn on Gauss points------------------------------------
             visc_dyn_gauss(indexs,1,i) = SUM(visc_dyns_loc(:,1)*mesh%gauss%wws(:,ls))
             visc_dyn_gauss(indexs,2,i) = SUM(visc_dyns_loc(:,2)*mesh%gauss%wws(:,ls))

             !-----------------Grad u_r on Gauss points------------------------------------
             grad1_vel(indexs,1,i) = SUM(vel_loc(:,1)*dws_loc(1,:))
             grad1_vel(indexs,2,i) = SUM(vel_loc(:,2)*dws_loc(1,:))
             grad1_vel(indexs,3,i) =  (mode*SUM(vels_loc(:,2)*mesh%gauss%wws(:,ls)) - &
                  SUM(vels_loc(:,3)*mesh%gauss%wws(:,ls)))/ray
             grad1_vel(indexs,4,i) = (-mode*SUM(vels_loc(:,1)*mesh%gauss%wws(:,ls)) - &
                  SUM(vels_loc(:,4)*mesh%gauss%wws(:,ls)))/ray
             grad1_vel(indexs,5,i) =  SUM(vel_loc(:,1)*dws_loc(2,:))
             grad1_vel(indexs,6,i) =  SUM(vel_loc(:,2)*dws_loc(2,:))

             !-----------------Grad u_th on Gauss points-----------------------------------
             grad2_vel(indexs,1,i) = SUM(vel_loc(:,3)*dws_loc(1,:))
             grad2_vel(indexs,2,i) = SUM(vel_loc(:,4)*dws_loc(1,:))
             grad2_vel(indexs,3,i) =  (mode*SUM(vels_loc(:,4)*mesh%gauss%wws(:,ls)) + &
                  SUM(vels_loc(:,1)*mesh%gauss%wws(:,ls)))/ray
             grad2_vel(indexs,4,i) = (-mode*SUM(vels_loc(:,3)*mesh%gauss%wws(:,ls)) + &
                  SUM(vels_loc(:,2)*mesh%gauss%wws(:,ls)))/ray
             grad2_vel(indexs,5,i) = SUM(vel_loc(:,3)*dws_loc(2,:))
             grad2_vel(indexs,6,i) = SUM(vel_loc(:,4)*dws_loc(2,:))

             !-----------------Grad u_z on Gauss points------------------------------------
             grad3_vel(indexs,1,i) = SUM(vel_loc(:,5)*dws_loc(1,:))
             grad3_vel(indexs,2,i) = SUM(vel_loc(:,6)*dws_loc(1,:))
             grad3_vel(indexs,3,i) =  mode*SUM(vels_loc(:,6)*mesh%gauss%wws(:,ls))/ray
             grad3_vel(indexs,4,i) = -mode*SUM(vels_loc(:,5)*mesh%gauss%wws(:,ls))/ray
             grad3_vel(indexs,5,i) = SUM(vel_loc(:,5)*dws_loc(2,:))
             grad3_vel(indexs,6,i) = SUM(vel_loc(:,6)*dws_loc(2,:))

             !-----------------GradT u_r on Gauss points------------------------------------
             gradT1_vel(indexs,1,i) = grad1_vel(indexs,1,i)
             gradT1_vel(indexs,2,i) = grad1_vel(indexs,2,i)
             gradT1_vel(indexs,3,i) = grad2_vel(indexs,1,i)
             gradT1_vel(indexs,4,i) = grad2_vel(indexs,2,i)
             gradT1_vel(indexs,5,i) = grad3_vel(indexs,1,i)
             gradT1_vel(indexs,6,i) = grad3_vel(indexs,2,i)
             !-----------------GradT u_th on Gauss points-----------------------------------
             gradT2_vel(indexs,1,i) = grad1_vel(indexs,3,i)
             gradT2_vel(indexs,2,i) = grad1_vel(indexs,4,i)
             gradT2_vel(indexs,3,i) = grad2_vel(indexs,3,i)
             gradT2_vel(indexs,4,i) = grad2_vel(indexs,4,i)
             gradT2_vel(indexs,5,i) = grad3_vel(indexs,3,i)
             gradT2_vel(indexs,6,i) = grad3_vel(indexs,4,i)
             !-----------------GradT u_z on Gauss points------------------------------------
             gradT3_vel(indexs,1,i) = grad1_vel(indexs,5,i)
             gradT3_vel(indexs,2,i) = grad1_vel(indexs,6,i)
             gradT3_vel(indexs,3,i) = grad2_vel(indexs,5,i)
             gradT3_vel(indexs,4,i) = grad2_vel(indexs,6,i)
             gradT3_vel(indexs,5,i) = grad3_vel(indexs,5,i)
             gradT3_vel(indexs,6,i) = grad3_vel(indexs,6,i)

             !-----------------Grad_sym_u scalar normal-------------------------------------
             V_bdy(indexs,1,i) = (grad1_vel(indexs,1,i)+gradT1_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad1_vel(indexs,5,i)+gradT1_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_bdy(indexs,2,i) = (grad1_vel(indexs,2,i)+gradT1_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad1_vel(indexs,6,i)+gradT1_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
             V_bdy(indexs,3,i) = (grad2_vel(indexs,1,i)+gradT2_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad2_vel(indexs,5,i)+gradT2_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_bdy(indexs,4,i) = (grad2_vel(indexs,2,i)+gradT2_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad2_vel(indexs,6,i)+gradT2_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
             V_bdy(indexs,5,i) = (grad3_vel(indexs,1,i)+gradT3_vel(indexs,1,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad3_vel(indexs,5,i)+gradT3_vel(indexs,5,i))*mesh%gauss%rnorms(2,ls,ms)
             V_bdy(indexs,6,i) = (grad3_vel(indexs,2,i)+gradT3_vel(indexs,2,i))*mesh%gauss%rnorms(1,ls,ms) &
                  + (grad3_vel(indexs,6,i)+gradT3_vel(indexs,6,i))*mesh%gauss%rnorms(2,ls,ms)
          ENDDO !l_Gs
       ENDDO !mes
    END DO !i

    !===Compute visc_dyn*(Grad_sym_u scalar normal)
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    bloc_size = SIZE(V_bdy,1)/nb_procs+1
    CALL FFT_SCALAR_VECT_DCL(communicator, V_bdy, visc_dyn_gauss, V_out, 1, nb_procs, bloc_size, m_max_pad)
    !===End Compute visc_dyn*(Grad_sym_u scalar normal)

  END SUBROUTINE smb_explicit_strain_rate_tensor_bdy_mom

  SUBROUTINE smb_explicit_tensor_bdy(mesh, list_mode, tensor, V_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                       :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)    :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(IN)  :: tensor
    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*mesh%dom_mes,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    INTEGER,      DIMENSION(mesh%gauss%n_ws)                 :: js_loc
    REAL(KIND=8), DIMENSION(3,mesh%gauss%n_ws,6)             :: tensors_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, ms, ls , i, mode, indexs, k, n

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       indexs = 0
       DO ms = 1, mesh%dom_mes
          js_loc = mesh%jjs(:,ms)
          m = mesh%neighs(ms)
          DO k = 1, 6
             DO n = 1, 3
                tensors_loc(n,:,k) = tensor(n,js_loc,k,i)
             END DO
          END DO
          DO ls = 1, mesh%gauss%l_Gs
             indexs = indexs + 1

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,js_loc)*mesh%gauss%wws(:,ls))

             !===Dont compute on r = 0
             IF (ray.LT.1.d-10) THEN
                V_out(indexs,:,i) = 0.d0
                CYCLE
             END IF
             !-----------------Tensor scalar normal on gauss points------------------------
             V_out(indexs,1,i) = SUM(tensors_loc(1,:,1)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(1,:,5)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,2,i) = SUM(tensors_loc(1,:,2)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(1,:,6)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,3,i) = SUM(tensors_loc(2,:,1)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(2,:,5)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,4,i) = SUM(tensors_loc(2,:,2)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(2,:,6)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,5,i) = SUM(tensors_loc(3,:,1)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(3,:,5)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
             V_out(indexs,6,i) = SUM(tensors_loc(3,:,2)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(1,ls,ms) &
                  + SUM(tensors_loc(3,:,6)*mesh%gauss%wws(:,ls))*mesh%gauss%rnorms(2,ls,ms)
          ENDDO !l_Gs
       ENDDO !mes
    END DO !i

  END SUBROUTINE smb_explicit_tensor_bdy

  SUBROUTINE compute_res_mass_gauss(mesh, list_mode, density_m2, density, momentum_m1, c_out)
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: density_m2, density
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: momentum_m1
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: Div
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: dens_m2_gauss, dens_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                :: mom_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                :: dens_m2_loc, dens_loc
    REAL(KIND=8)                                :: ray
    INTEGER                                     :: m, l , i, mode, index, k

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          DO k = 1, 6
             mom_loc(:,k) = momentum_m1(j_loc,k,i)
          END DO
          DO k = 1, 2
             dens_m2_loc(:,k) = density_m2(j_loc,k,i)
             dens_loc(:,k)    = density(j_loc,k,i)
          END DO
          DO l = 1, mesh%gauss%l_G
             index = index + 1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

             !===Compute density_m2 and density on gauss point
             dens_m2_gauss(index,1,i) = SUM(dens_m2_loc(:,1)*mesh%gauss%ww(:,l))
             dens_m2_gauss(index,2,i) = SUM(dens_m2_loc(:,2)*mesh%gauss%ww(:,l))

             dens_gauss(index,1,i) = SUM(dens_loc(:,1)*mesh%gauss%ww(:,l))
             dens_gauss(index,2,i) = SUM(dens_loc(:,2)*mesh%gauss%ww(:,l))

             !===Compute divergence of momentum on gauss point
             Div(index,1,i) = SUM(mom_loc(:,1)*dw_loc(1,:)) + SUM(mom_loc(:,1)*mesh%gauss%ww(:,l))/ray &
                  + (mode/ray)*SUM(mom_loc(:,4)*mesh%gauss%ww(:,l)) +  SUM(mom_loc(:,5)*dw_loc(2,:))
             Div(index,2,i) = SUM(mom_loc(:,2)*dw_loc(1,:)) + SUM(mom_loc(:,2)*mesh%gauss%ww(:,l))/ray &
                  - (mode/ray)*SUM(mom_loc(:,3)*mesh%gauss%ww(:,l)) +  SUM(mom_loc(:,6)*dw_loc(2,:))
          END DO
       END DO
    END DO

    c_out = (dens_gauss-dens_m2_gauss)/(2*inputs%dt) + Div

  END SUBROUTINE compute_res_mass_gauss

  SUBROUTINE twoD_volume(communicator,mesh,RESLT)
    !===========================
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
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
  END SUBROUTINE twoD_volume

END MODULE entropy_viscosity
