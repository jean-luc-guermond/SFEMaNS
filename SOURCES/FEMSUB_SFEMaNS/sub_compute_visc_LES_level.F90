!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_compute_visc_LES_level
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC :: compute_visc_LES_level
  PRIVATE

CONTAINS

  SUBROUTINE compute_visc_LES_level(comm_one_d,time, cc_1_LA, list_mode, cc_mesh, cn_m1, cn, &
       my_par_cc, cc_list_dirichlet_sides, cc_per, nb_inter, visc_entro_level, ff_entro)
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
    USE sft_parallele
    USE input_data
    USE subroutine_level_set
    IMPLICIT NONE
    REAL(KIND=8)                                        :: time
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: list_mode
    INTEGER,                        INTENT(IN)          :: nb_inter
    TYPE(mesh_type),                INTENT(IN)          :: cc_mesh
    type(petsc_csr_LA)                                  :: cc_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: cc_per
    TYPE(solver_param),             INTENT(IN)          :: my_par_cc
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: cn_m1, cn
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: cc_list_dirichlet_sides
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)          :: visc_entro_level
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)         :: ff_entro
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: cc_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: cc_mode_global_js_D
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: m_max_c
    INTEGER,     DIMENSION(:),   POINTER,          SAVE :: cc_js_D ! Dirichlet nodes
    INTEGER,                                       SAVE :: my_petscworld_rank, nb_procs
    REAL(KIND=8),                                  SAVE :: LES_coeff1_in_level
    !----------FIN SAVE--------------------------------------------------------------------

    !----------Declaration sans save-------------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: cc_1_ifrom
    INTEGER                                  :: i, n
    INTEGER                                  :: code, mode
    INTEGER                                  :: bloc_size, m_max_pad
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(cc_mesh%np)                      :: ff
    REAL(KIND=8), DIMENSION(cc_mesh%np,2,SIZE(list_mode))    :: cext
    REAL(KIND=8), DIMENSION(cc_mesh%np,2,SIZE(list_mode))    :: cext_reg
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 2, SIZE(list_mode)) :: ff_comp
!    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 6, SIZE(list_mode)) :: ff_entro
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 2, SIZE(list_mode)) :: ff_phi_1mphi
    REAL(KIND=8) :: tps, tps_tot, tps_cumul
    !Communicators for Petsc, in space and Fourier------------------------------
    !#include "petsc/finclude/petsc.h"
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: cc_mat
    Vec,                        SAVE :: cb_1, cb_2, cx_1, cx_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: cc_ksp
    Vec,                        SAVE :: cb_reg_1, cb_reg_2 !vectors for level set regularization
    Mat, DIMENSION(:), POINTER, SAVE :: cc_reg_mat
    KSP, DIMENSION(:), POINTER, SAVE :: cc_reg_ksp
    !------------------------------END OF DECLARATION--------------------------------------

    IF (once) THEN

       once = .FALSE.

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)

       !-----CREATE PETSC VECTORS AND GHOSTS-----------------------------------------
       CALL create_my_ghost(cc_mesh,cc_1_LA,cc_1_ifrom)
       n = cc_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(cc_1_ifrom), cc_1_ifrom, cx_1, ierr)
       CALL VecGhostGetLocalForm(cx_1, cx_1_ghost, ierr)
       CALL VecDuplicate(cx_1, cb_1, ierr)
       CALL VecDuplicate(cx_1, cb_2, ierr)
       CALL VecDuplicate(cx_1, cb_reg_1, ierr)
       CALL VecDuplicate(cx_1, cb_reg_2, ierr)
       !------------------------------------------------------------------------------

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------

       !---------PREPARE cc_js_D ARRAY FOR PHASE--------------------------------------
       CALL dirichlet_nodes_parallel(cc_mesh, cc_list_dirichlet_sides, cc_js_D)
       !===JLG June 9 2017, replaced scalar_glob_js_D by scalar_with_bc_glob_js_D
       !CALL scalar_glob_js_D(cc_mesh, list_mode, cc_1_LA, cc_mode_global_js_D)
       CALL scalar_with_bc_glob_js_D(cc_mesh, list_mode, cc_1_LA, cc_js_D, cc_mode_global_js_D)
       ALLOCATE(cc_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(cc_global_D(i)%DRL(SIZE(cc_mode_global_js_D(i)%DIL)))
       END DO
       !------------------------------------------------------------------------------

       !-------------ASSEMBLE PHASE MATRICES------------------------------------------
       ALLOCATE(cc_mat(m_max_c),cc_ksp(m_max_c))

       !===Assembling regularization left hand side
       IF (inputs%if_compression_mthd_JLG) THEN
          ALLOCATE(cc_reg_mat(m_max_c),cc_reg_ksp(m_max_c))
          DO i = 1, m_max_c
             mode = list_mode(i)

             !---PHASE MATRIX for level set regularization
             CALL create_local_petsc_matrix(comm_one_d(1), cc_1_LA, cc_reg_mat(i), clean=.FALSE.)
             CALL qs_regul_M (cc_mesh, cc_1_LA, 3.d0, i, mode, cc_reg_mat(i))
             IF (cc_per%n_bord/=0) THEN
                CALL periodic_matrix_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cc_reg_mat(i), cc_1_LA)
             END IF
             CALL Dirichlet_M_parallel(cc_reg_mat(i),cc_mode_global_js_D(i)%DIL)
             CALL init_solver(my_par_cc,cc_reg_ksp(i),cc_reg_mat(i),comm_one_d(1),&
                  solver=my_par_cc%solver,precond=my_par_cc%precond)
          ENDDO
       END IF
       !===End Assembling regularization left hand side

       IF (inputs%if_level_set_P2) THEN
          LES_coeff1_in_level=inputs%LES_coeff1
       ELSE
          LES_coeff1_in_level=4.d0*inputs%LES_coeff1
       END IF
    ENDIF

    tps_tot = user_time()
    tps_cumul = 0

    IF (inputs%if_level_bdf2) THEN
       cext = 2.d0*cn-cn_m1
    ELSE
       cext = cn
    END IF

    IF (inputs%if_compression_mthd_JLG) THEN
       !---------------REGULARIZATION OF LEVEL SET FOR COMPRESSION---------------------------
       DO i = 1, m_max_c
          !===Compute rhs for level set regularization
          CALL qs_00 (cc_mesh,cc_1_LA, cext(:,1,i), cb_reg_1)
          CALL qs_00 (cc_mesh,cc_1_LA, cext(:,2,i), cb_reg_2)

          !===RHS periodicity
          IF (cc_per%n_bord/=0) THEN
             CALL periodic_rhs_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cb_reg_1, cc_1_LA)
             CALL periodic_rhs_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cb_reg_2, cc_1_LA)
          END IF

          !===RHS Dirichlet
          n = SIZE(cc_js_D)
          cc_global_D(i)%DRL(1:n) = level_set_exact(nb_inter,1,cc_mesh%rr(:,cc_js_D), mode, time)
          cc_global_D(i)%DRL(n+1:) = 0.d0
          CALL dirichlet_rhs(cc_mode_global_js_D(i)%DIL-1,cc_global_D(i)%DRL,cb_reg_1)
          cc_global_D(i)%DRL(1:n) = level_set_exact(nb_inter,2,cc_mesh%rr(:,cc_js_D), mode, time)
          cc_global_D(i)%DRL(n+1:) = 0.d0
          CALL dirichlet_rhs(cc_mode_global_js_D(i)%DIL-1,cc_global_D(i)%DRL,cb_reg_2)

          !===Solve level set regularization equation
          tps = user_time()
          !Solve system cc_c
          CALL solver(cc_reg_ksp(i),cb_reg_1,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
          CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL extract(cx_1_ghost,1,1,cc_1_LA,cext_reg(:,1,i))

          !Solve system cc_s
          CALL solver(cc_reg_ksp(i),cb_reg_2,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
          CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
          CALL extract(cx_1_ghost,1,1,cc_1_LA,cext_reg(:,2,i))
          tps = user_time() - tps; tps_cumul=tps_cumul+tps
          !WRITE(*,*) ' Tps solution des pb de vitesse', tps, 'for mode ', mode

       END DO
       !--------------- END REGULARIZATION OF LEVEL SET FOR COMPRESSION----------------------
    END IF

    !===Compute visc_LES_level (includes c1_LES and compression terms)
    IF (inputs%if_compression_mthd_JLG) THEN
       CALL smb_compr_visc_entro_gauss_fft_par(comm_one_d, cc_mesh, &
            list_mode, cext, cext_reg, visc_entro_level, &
            LES_coeff1_in_level, nb_procs, ff_entro, ff_phi_1mphi)
    ELSE
       CALL smb_visc_entro_gauss_fft_par(comm_one_d, cc_mesh, list_mode, cext, visc_entro_level, &
            LES_coeff1_in_level, nb_procs, ff_entro, ff_phi_1mphi)
    END IF

  END SUBROUTINE compute_visc_LES_level
  !============================================

END MODULE subroutine_compute_visc_LES_level
