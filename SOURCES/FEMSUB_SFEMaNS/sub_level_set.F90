!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_level_set
  USE my_util
#include "petsc/finclude/petsc.h"
  USE petsc
!  PUBLIC :: three_level_level_set
!TEST LC LES_SUITE 2024/06
  PUBLIC :: three_level_level_set, smb_compr_visc_entro_gauss_fft_par, smb_visc_entro_gauss_fft_par, qs_regul_M
!TEST LC LES_SUITE 2024/06
  PRIVATE
  REAL(KIND=8)  :: max_velocity_at_tn
  !===JLG Sept 27, 2016
  !  LOGICAL :: compression_mthd_JLG=.TRUE.
  !===JLG Sept 27, 2016

CONTAINS

  SUBROUTINE three_level_level_set(comm_one_d,time, cc_1_LA, dt, list_mode, cc_mesh, cn_m1, cn, &
       chmp_vit, max_vel, my_par_cc, cc_list_dirichlet_sides, cc_per, nb_inter, &
       visc_entro_level, cext_reg, visc_LES_level)
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
    IMPLICIT NONE
    REAL(KIND=8)                                        :: time, dt
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: list_mode
    INTEGER,                        INTENT(IN)          :: nb_inter
    TYPE(mesh_type),                INTENT(IN)          :: cc_mesh
    type(petsc_csr_LA)                                  :: cc_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: cc_per
    TYPE(solver_param),             INTENT(IN)          :: my_par_cc
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)       :: cn_m1, cn
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: cc_list_dirichlet_sides
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: chmp_vit
    REAL(KIND=8),                   INTENT(INOUT)       :: max_vel
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)          :: visc_entro_level
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)         :: cext_reg
!TEST LC LES_SUITE 2024/06
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: visc_LES_level
!TEST LC LES_SUITE 2024/06
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: cc_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: cc_mode_global_js_D
    LOGICAL,                                       SAVE :: once = .TRUE., once_vel=.TRUE.
!TEST LC LES_SUITE 2024/06
    LOGICAL,                                       SAVE :: once_LES = .TRUE.
!TEST LC LES_SUITE 2024/06
    LOGICAL,                                       SAVE :: re_init=.FALSE.
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
    REAL(KIND=8), DIMENSION(cc_mesh%np, 2)                   :: cn_p1
    REAL(KIND=8), DIMENSION(cc_mesh%np,2,SIZE(list_mode))    :: cext
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 2, SIZE(list_mode)) :: ff_conv
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 2, SIZE(list_mode)) :: ff_comp
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 6, SIZE(list_mode)) :: ff_entro
    REAL(KIND=8), DIMENSION(cc_mesh%gauss%l_G*cc_mesh%dom_me, 2, SIZE(list_mode)) :: ff_phi_1mphi
    REAL(KIND=8) :: int_mass_correct
    REAL(KIND=8) :: tps, tps_tot, tps_cumul
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/
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
       max_velocity_at_tn = 0.d0

       IF (inputs%if_level_set_P2) THEN
          LES_coeff1_in_level=inputs%LES_coeff1
       ELSE
          LES_coeff1_in_level=4.d0*inputs%LES_coeff1
       END IF
    ENDIF

    !===Assembling left hand side
    IF (once_vel) THEN
       re_init=.FALSE.
       once_vel = .FALSE.

       ! Take care of once_vel
       max_vel = MAX(1.1d0*max_velocity_at_tn,max_vel)
       IF (my_petscworld_rank==0) WRITE(*,*) ' Recomputing matrix for level set function'
       IF (my_petscworld_rank==0) WRITE(*,*) ' NEW MAX VEL test 1', time, max_vel

       DO i = 1, m_max_c
          mode = list_mode(i)

          !---PHASE MATRIX
          CALL create_local_petsc_matrix(comm_one_d(1), cc_1_LA, cc_mat(i), clean=.FALSE.)
          IF (inputs%if_level_bdf2) THEN
             CALL qs_diff_mass_scal_M_level(cc_mesh, cc_1_LA, 0.d0, 1.5d0/dt, &
                  LES_coeff1_in_level, i, mode, cc_mat(i))  !===Coefficient 4 because P1 approximation
          ELSE
             CALL qs_diff_mass_scal_M_level(cc_mesh, cc_1_LA, 0.d0, 1.d0/dt, &
                  LES_coeff1_in_level, i, mode, cc_mat(i)) !===Coefficient 4 because P1 approximation
          END IF
          IF (cc_per%n_bord/=0) THEN
             CALL periodic_matrix_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cc_mat(i), cc_1_LA)
          END IF
          CALL Dirichlet_M_parallel(cc_mat(i),cc_mode_global_js_D(i)%DIL)
          CALL init_solver(my_par_cc,cc_ksp(i),cc_mat(i),comm_one_d(1),&
               solver=my_par_cc%solver,precond=my_par_cc%precond, opt_re_init=re_init)
       END DO
    END IF
    !===End assembling left hand side

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

    !===Compute rhs by FFT at Gauss points
    CALL smb_ugradc_gauss_fft_par(comm_one_d(2), cc_mesh, list_mode, chmp_vit, cext, nb_procs, ff_conv)

    IF (inputs%if_compression_mthd_JLG) THEN
       CALL smb_compr_visc_entro_gauss_fft_par(comm_one_d, cc_mesh, &
            list_mode, cext, cext_reg, visc_entro_level,&
            LES_coeff1_in_level, nb_procs, ff_entro, ff_phi_1mphi)
    ELSE
       CALL smb_compression_gauss_fft_par(comm_one_d, cc_mesh, list_mode, &
            chmp_vit, cext, nb_procs, ff_comp)
       ff_conv=ff_conv-ff_comp

       CALL smb_visc_entro_gauss_fft_par(comm_one_d, cc_mesh, list_mode, cext, visc_entro_level,&
            LES_coeff1_in_level, nb_procs, ff_entro, ff_phi_1mphi)
    END IF

!TEST LC LES_SUITE 2024/06
    !===Use restart to initialize ff_entro (if restart_LES=.t.)
    IF (once_LES) THEN
       once_LES = .FALSE.
       IF (inputs%irestart_LES) THEN
          ff_entro = visc_LES_level
       END IF
    END IF
    !===End Use restart to initialize ff_entro (if restart_LES=.t.)
!TEST LC LES_SUITE 2024/06

    IF (inputs%if_mass_correction) THEN
       CALL compute_int_mass_correct(comm_one_d, cc_mesh, list_mode, ff_conv, ff_phi_1mphi, int_mass_correct)
    ELSE
       int_mass_correct = 0.d0
    END IF

    !------------DEBUT BOUCLE SUR LES MODES----------------
    DO i = 1, m_max_c
       mode = list_mode(i)

       !===RHS phase
       IF (inputs%if_level_bdf2) THEN
          ff = (2/dt)*cn(:,1,i) - (1/(2*dt))*cn_m1(:,1,i) &
               +source_in_level_set(nb_inter,1, cc_mesh%rr, mode, time)
          CALL qs_00_level_set_gauss (cc_mesh, cc_1_LA, ff, -ff_conv(:,1,i), &
               mode, 1, cb_1, cext(:,1,i),  &
               ff_entro(:,:,i), -ff_phi_1mphi(:,1,i), int_mass_correct)

          ff = (2/dt)*cn(:,2,i) - (1/(2*dt))*cn_m1(:,2,i) &
               +source_in_level_set(nb_inter,2, cc_mesh%rr, mode, time)
          CALL qs_00_level_set_gauss (cc_mesh, cc_1_LA, ff, -ff_conv(:,2,i),&
               mode, 2, cb_2, cext(:,2,i), &
               ff_entro(:,:,i), -ff_phi_1mphi(:,2,i), int_mass_correct)
       ELSE !BDF1
          ff = (1/dt)*cn(:,1,i) + source_in_level_set(nb_inter,1, cc_mesh%rr, mode, time)
          CALL qs_00_level_set_gauss (cc_mesh, cc_1_LA, ff, -ff_conv(:,1,i), &
               mode, 1, cb_1, cext(:,1,i),  &
               ff_entro(:,:,i), -ff_phi_1mphi(:,1,i), int_mass_correct)

          ff = (1/dt)*cn(:,2,i) + source_in_level_set(nb_inter,2, cc_mesh%rr, mode, time)
          CALL qs_00_level_set_gauss (cc_mesh, cc_1_LA, ff, -ff_conv(:,2,i),&
               mode, 2, cb_2, cext(:,2,i), &
               ff_entro(:,:,i), -ff_phi_1mphi(:,2,i), int_mass_correct)
       END IF

       !===RHS periodicity
       IF (cc_per%n_bord/=0) THEN
          CALL periodic_rhs_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cb_1, cc_1_LA)
          CALL periodic_rhs_petsc(cc_per%n_bord, cc_per%list, cc_per%perlist, cb_2, cc_1_LA)
       END IF
       !----------------------------------------------------------

       n = SIZE(cc_js_D)
       cc_global_D(i)%DRL(1:n) = level_set_exact(nb_inter,1,cc_mesh%rr(:,cc_js_D), mode, time)
       cc_global_D(i)%DRL(n+1:) = 0.d0
       CALL dirichlet_rhs(cc_mode_global_js_D(i)%DIL-1,cc_global_D(i)%DRL,cb_1)
       cc_global_D(i)%DRL(1:n) = level_set_exact(nb_inter,2,cc_mesh%rr(:,cc_js_D), mode, time)
       cc_global_D(i)%DRL(n+1:) = 0.d0
       CALL dirichlet_rhs(cc_mode_global_js_D(i)%DIL-1,cc_global_D(i)%DRL,cb_2)

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps second membre vitesse', tps
       !-------------------------------------------------------------------------------------

       !--------------------INVERSION DES OPERATEURS--------------
       tps = user_time()
       !Solve system cc_c
       CALL solver(cc_ksp(i),cb_1,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,cc_1_LA,cn_p1(:,1))

       !Solve system cc_s
       CALL solver(cc_ksp(i),cb_2,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,cc_1_LA,cn_p1(:,2))
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de vitesse', tps, 'for mode ', mode
       !-------------------------------------------------------------------------------------

       !---------------UPDATES-----------------------
       tps = user_time()

       IF (mode==0) THEN
          cn_p1 (:,2) = 0.d0
       END IF

       cn_m1(:,:,i) = cn(:,:,i)
       cn   (:,:,i) = cn_p1

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps  des updates', tps
       !-------------------------------------------------------------------------------------
    ENDDO

    bloc_size = SIZE(cn,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$    IF (inputs%if_kill_overshoot) THEN
!!$       IF (nb_procs==1.AND.SIZE(list_mode)==1.AND.list_mode(1)==0) THEN !case axisym
!!$          cn = MIN(1.d0, cn)
!!$          cn = MAX(0.d0, cn)
!!$       ELSE !level set depends of theta
!!$          CALL FFT_NO_OVERSHOOT_LEVEL_SET(comm_one_d(2), cn, nb_procs, bloc_size, m_max_pad)
!!$       END IF
!!$    END IF

    tps_tot = user_time() - tps_tot
    !WRITE(*,'(A,2(f13.3,2x))') '  Tps boucle en temps Navier_stokes', tps_tot, tps_cumul
    !WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE three_level_level_set
  !============================================

  SUBROUTINE smb_ugradc_gauss_fft_par(communicator,mesh,list_mode,V_in,c_in, nb_procs, c_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in, c_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: Gradc, W
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: Div, Cgauss
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: cint
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)   :: cs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray, tps
    REAL(KIND=8), DIMENSION(3)                  :: temps
    INTEGER                                     :: m_max_pad, bloc_size
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    tps = user_time()
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
          END DO
          DO k = 1, 2
             cs(:,k) = c_in(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Vs(:,1)*ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*ww(:,l))

             Div(index,1,i) = SUM(Vs(:,1)*dw_loc(1,:)) + SUM(Vs(:,1)*ww(:,l))/ray &
                  + (mode/ray)*SUM(Vs(:,4)*ww(:,l)) +  SUM(Vs(:,5)*dw_loc(2,:))
             Div(index,2,i) = SUM(Vs(:,2)*dw_loc(1,:)) + SUM(Vs(:,2)*ww(:,l))/ray &
                  - (mode/ray)*SUM(Vs(:,3)*ww(:,l)) +  SUM(Vs(:,6)*dw_loc(2,:))

             !-----------------gradient de c sur les points de Gauss---------------------------
             !coeff sur les cosinus et sinus
             Gradc(index,1,i) = SUM(cs(:,1)*dw_loc(1,:))
             Gradc(index,2,i) = SUM(cs(:,2)*dw_loc(1,:))
             Gradc(index,3,i) =  mode/ray*SUM(cs(:,2)*ww(:,l))
             Gradc(index,4,i) = -mode/ray*SUM(cs(:,1)*ww(:,l))
             Gradc(index,5,i) = SUM(cs(:,1)*dw_loc(2,:))
             Gradc(index,6,i) = SUM(cs(:,2)*dw_loc(2,:))

             Cgauss(index,1,i) = SUM(cs(:,1)*ww(:,l))
             Cgauss(index,2,i) = SUM(cs(:,2)*ww(:,l))

          ENDDO
       ENDDO
    END DO

    !tps = user_time() - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    !tps = user_time()
    temps = 0

    bloc_size = SIZE(Gradc,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    CALL FFT_PAR_DOT_PROD_DCL(communicator, Gradc, W, c_out, nb_procs, bloc_size, m_max_pad, temps)
!!$    !===Conservative term (phi div(u))
!!$     bloc_size = SIZE(Div,1)/nb_procs+1
!!$     CALL FFT_PAR_PROD_DCL(communicator, Div, Cgauss, cint, nb_procs, bloc_size, m_max_pad, temps)
!!$     c_out = c_out + cint
!!$    !===End Conservative term (phi div(u))
    tps = user_time() - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_ugradc_gauss_fft_par

  SUBROUTINE smb_compression_gauss_fft_par(communicator,mesh,list_mode,V_in,c_in, nb_procs, c_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE tn_axi
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in, c_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: Gradc, W
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: Cgauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)   :: cs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray, tps, norm_vel_L2, Volume_3D
    INTEGER                                     :: m_max_pad, bloc_size
    !#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(2)       :: communicator

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    tps = user_time()
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = jj(:,m)
          DO k = 1, 6
             Vs(:,k) = V_in(j_loc,k,i)
          END DO
          DO k = 1, 2
             cs(:,k) = c_in(j_loc,k,i)
          END DO
          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------vitesse sur les points de Gauss---------------------------
             W(index,1,i) = SUM(Vs(:,1)*ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*ww(:,l))

             !-----------------gradient de c sur les points de Gauss---------------------------
             !coeff sur les cosinus et sinus
             Gradc(index,1,i) = SUM(cs(:,1)*dw_loc(1,:))
             Gradc(index,2,i) = SUM(cs(:,2)*dw_loc(1,:))
             Gradc(index,3,i) =  mode/ray*SUM(cs(:,2)*ww(:,l))
             Gradc(index,4,i) = -mode/ray*SUM(cs(:,1)*ww(:,l))
             Gradc(index,5,i) = SUM(cs(:,1)*dw_loc(2,:))
             Gradc(index,6,i) = SUM(cs(:,2)*dw_loc(2,:))

             Cgauss(index,1,i) = SUM(cs(:,1)*ww(:,l))
             Cgauss(index,2,i) = SUM(cs(:,2)*ww(:,l))
          ENDDO
       ENDDO
    END DO

    bloc_size = mesh%gauss%l_G*mesh%me/nb_procs+1
    bloc_size = mesh%gauss%l_G*(bloc_size/mesh%gauss%l_G)+mesh%gauss%l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !TEST JLG LC vel L2
    CALL twoD_volume(communicator(1),mesh,Volume_3D)
    Volume_3D = Volume_3D*2*ACOS(-1.d0)
    norm_vel_L2 = norm_SF(communicator, 'L2', mesh, list_mode, V_in)/Volume_3D
    !TEST JLG LC vel L2

    CALL FFT_COMPRESSION_LEVEL_SET_DCL(communicator(2), communicator(1),Gradc, W, Cgauss, c_out, &
         mesh%hloc_gauss, mesh%gauss%l_G, nb_procs, bloc_size, m_max_pad)

    tps = user_time() - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_compression_gauss_fft_par

  SUBROUTINE twoD_volume(communicator,mesh,RESLT)
    !===========================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    REAL(KIND=8),                INTENT(OUT) :: RESLT
    REAL(KIND=8)                             :: vol_loc, vol_out
    INTEGER,      DIMENSION(mesh%gauss%n_w)  :: j_loc
    INTEGER                                  :: m, l , i , ni, code
    REAL(KIND=8)                             :: ray
    !#include "petsc/finclude/petsc.h"
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

  SUBROUTINE qs_regul_M (mesh, LA, stab, i_mode, mode, matrix)
    !=================================================
    USE def_type_mesh
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    type(petsc_csr_LA)                                     :: LA
    REAL(KIND=8),                 INTENT(IN)               :: stab
    INTEGER,                      INTENT(IN)               :: mode, i_mode
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc
    REAL(KIND=8)                                           :: ray
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, n_w
    REAL(KIND=8) :: viscolm, xij, viscomode, hm, hh
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    n_w = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       a_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          hh=mesh%hloc(m)
          hm=MIN(mesh%hm(i_mode),hh)!WRONG choice
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)

          viscolm  = (stab*hh)**2*mesh%gauss%rj(l,m)
          viscomode = (stab*hm)**2*mesh%gauss%rj(l,m)
          DO nj = 1, n_w
             DO ni = 1, n_w
                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij  &
                     + ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix,n_w,idxm,n_w,idxn,a_loc,ADD_VALUES,ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_regul_M

  SUBROUTINE smb_compr_visc_entro_gauss_fft_par(communicator, mesh, list_mode, c_in, c_reg_in, visc_entro_real,&
       coeff1_in_level, nb_procs, V_out, c_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: c_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: c_reg_in
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: visc_entro_real
    REAL(KIND=8),                   INTENT(IN)  :: coeff1_in_level
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)) :: Gradc_in, Gradc_reg_in
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)) :: c_in_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                :: c_in_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                :: c_reg_in_loc
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)                                :: ray, hh, hm
    INTEGER                                     :: m, l , i, mode, index, k
    INTEGER                                     :: m_max_pad, bloc_size
    !#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(2)       :: communicator

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = jj(:,m)
          DO k = 1, 2
             c_in_loc(:,k) = c_in(j_loc,k,i)
             c_reg_in_loc(:,k) = c_reg_in(j_loc,k,i)
          END DO

          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute local mesh sizes
             hh=mesh%hloc_gauss(index)
             hm=MIN(mesh%hm(i),hh)!WRONG choice
             !===Hm must have a dimension (JLG, April 7 2017)
             !hm=0.5d0/inputs%m_max
             !hm=mesh%hm(i) !(JLG April 7 2017)
             !===Hm did not have a dimension (JLG, April 7 2017)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------c_in on Gauss points------------------
             c_in_gauss(index,1,i) = SUM(c_in_loc(:,1)*ww(:,l))
             c_in_gauss(index,2,i) = SUM(c_in_loc(:,2)*ww(:,l))

             !----------Gradient of c_in on Gauss points---------
             Gradc_in(index,1,i) = SUM(c_in_loc(:,1)*dw_loc(1,:))*hh
             Gradc_in(index,2,i) = SUM(c_in_loc(:,2)*dw_loc(1,:))*hh
             Gradc_in(index,3,i) = mode/ray*SUM(c_in_loc(:,2)*ww(:,l))*hm
             Gradc_in(index,4,i) = -mode/ray*SUM(c_in_loc(:,1)*ww(:,l))*hm
             Gradc_in(index,5,i) = SUM(c_in_loc(:,1)*dw_loc(2,:))*hh
             Gradc_in(index,6,i) = SUM(c_in_loc(:,2)*dw_loc(2,:))*hh

             !----------Gradient of c_reg_in on Gauss points---------
             Gradc_reg_in(index,1,i) = SUM(c_reg_in_loc(:,1)*dw_loc(1,:))*hh
             Gradc_reg_in(index,2,i) = SUM(c_reg_in_loc(:,2)*dw_loc(1,:))*hh
             Gradc_reg_in(index,3,i) = mode/ray*SUM(c_reg_in_loc(:,2)*ww(:,l))*hm
             Gradc_reg_in(index,4,i) = -mode/ray*SUM(c_reg_in_loc(:,1)*ww(:,l))*hm
             Gradc_reg_in(index,5,i) = SUM(c_reg_in_loc(:,1)*dw_loc(2,:))*hh
             Gradc_reg_in(index,6,i) = SUM(c_reg_in_loc(:,2)*dw_loc(2,:))*hh
          END DO
       END DO
    END DO
    bloc_size = SIZE(c_in_gauss,1)/nb_procs+1
    bloc_size = l_G*(bloc_size/l_G)+l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===Compute c_out = max(0,phi*(1-phi))
    !===Compute V_out = (coeff1_in_level-visc_entro)*Grad(phi) + comp*visc_entro/(h_loc*|Grad(phi_reg)(x,t)|)*c_out*Grad(phi_reg)
    CALL FFT_PAR_COMPR_ENTRO_VISC_DCL(communicator(2), Gradc_in, Gradc_reg_in, c_in_gauss, visc_entro_real, &
         mesh%hloc_gauss, coeff1_in_level, V_out, c_out, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_compr_visc_entro_gauss_fft_par

  SUBROUTINE smb_visc_entro_gauss_fft_par(communicator, mesh, list_mode, c_in, visc_entro_real,&
       coeff1_in_level, nb_procs, V_out, c_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: nb_procs
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: c_in
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: visc_entro_real
    REAL(KIND=8),                   INTENT(IN)  :: coeff1_in_level
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: V_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,6,SIZE(list_mode)) :: Gradc_in
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%dom_me,2,SIZE(list_mode)) :: c_in_gauss
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                :: c_in_loc
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)                                :: ray, hh, hm
    INTEGER                                     :: m, l , i, mode, index, k
    INTEGER                                     :: m_max_pad, bloc_size
    !#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(2)       :: communicator

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          j_loc = jj(:,m)
          DO k = 1, 2
             c_in_loc(:,k) = c_in(j_loc,k,i)
          END DO

          DO l = 1, l_G
             index = index + 1
             dw_loc = dw(:,:,l,m)

             !===Compute local mesh sizes
             hh=mesh%hloc_gauss(index)
             hm=MIN(mesh%hm(i),hh)!WRONG choice
             !hm=0.5d0/inputs%m_max
             !hm=mesh%hm(i) !(JLG April 7 2017)

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,j_loc)*ww(:,l))

             !-----------------c_in on Gauss points------------------
             c_in_gauss(index,1,i) = SUM(c_in_loc(:,1)*ww(:,l))
             c_in_gauss(index,2,i) = SUM(c_in_loc(:,2)*ww(:,l))

             !----------Gradient of c_in on Gauss points---------
             Gradc_in(index,1,i) = SUM(c_in_loc(:,1)*dw_loc(1,:))*hh
             Gradc_in(index,2,i) = SUM(c_in_loc(:,2)*dw_loc(1,:))*hh
             Gradc_in(index,3,i) = mode/ray*SUM(c_in_loc(:,2)*ww(:,l))*hm
             Gradc_in(index,4,i) = -mode/ray*SUM(c_in_loc(:,1)*ww(:,l))*hm
             Gradc_in(index,5,i) = SUM(c_in_loc(:,1)*dw_loc(2,:))*hh
             Gradc_in(index,6,i) = SUM(c_in_loc(:,2)*dw_loc(2,:))*hh
          END DO
       END DO
    END DO
    !CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(c_in_gauss,1)/nb_procs+1
    bloc_size = l_G*(bloc_size/l_G)+l_G
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===Compute c_out = max(0,phi*(1-phi))
    !===Compute V_out = (coeff1_in_level-visc_entro)*Grad(phi)
    CALL FFT_PAR_ENTRO_VISC_DCL(communicator(2), Gradc_in, c_in_gauss, visc_entro_real, &
         mesh%hloc_gauss, coeff1_in_level, V_out, c_out, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_visc_entro_gauss_fft_par

  SUBROUTINE compute_int_mass_correct(communicator, mesh, list_mode, c1_in, c2_in, int_out)
    !=================================
    USE Gauss_points
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: c1_in  !=ff_conv
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: c2_in  !=ff_phi_1mphi
    REAL(KIND=8),                   INTENT(OUT) :: int_out
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: j_loc
    REAL(KIND=8)                                :: ray
    REAL(KIND=8)                                :: int_c1_in, int_c1_in_loc, int_c1_in_F
    REAL(KIND=8)                                :: int_c2_in, int_c2_in_loc, int_c2_in_F
    INTEGER                                     :: m, l , i, index, code
    !#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(2)       :: communicator

    int_c1_in_loc = 0.d0
    int_c2_in_loc = 0.d0

    DO i = 1, SIZE(list_mode)
       index = 0
       IF (list_mode(i)==0) THEN
          DO m = 1, mesh%dom_me
             j_loc = mesh%jj(:,m)
             DO l = 1, mesh%gauss%l_G
                index = index + 1
                !===Compute radius of Gauss point
                ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

                !===Compute integral of c1_in and c2_in
                int_c1_in_loc = int_c1_in_loc + c1_in(index,1,i)*ray*mesh%gauss%rj(l,m)
                int_c2_in_loc = int_c2_in_loc + c2_in(index,1,i)*ray*mesh%gauss%rj(l,m)
             END DO
          END DO
       END IF
    END DO

    int_c1_in_loc = int_c1_in_loc*2*ACOS(-1.d0)
    CALL  MPI_ALLREDUCE(int_c1_in_loc, int_c1_in_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         communicator(2), code)
    CALL MPI_ALLREDUCE(int_c1_in_F, int_c1_in, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         communicator(1), code)

    int_c2_in_loc = int_c2_in_loc*2*ACOS(-1.d0)
    CALL  MPI_ALLREDUCE(int_c2_in_loc, int_c2_in_F, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         communicator(2), code)
    CALL MPI_ALLREDUCE(int_c2_in_F, int_c2_in, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         communicator(1), code)

    IF (int_c2_in.LE.1.d-14) THEN
       int_out=0.d0
    ELSE
       int_out=-int_c1_in/int_c2_in
    END IF

  END SUBROUTINE compute_int_mass_correct

  SUBROUTINE qs_00_level_set_gauss (mesh, LA, ff,  ff_gauss, mode, type, vect, level_set_ext, &
       fcompr, ff_phi_1mphi, stab_mass)
    !=================================
    USE def_type_mesh
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    TYPE(petsc_csr_LA)                          :: LA
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff, ff_gauss
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: level_set_ext
    INTEGER     ,                 INTENT(IN)    :: mode
    INTEGER     ,                 INTENT(IN)    :: type ! 1 = cosine, 2 = sine
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: fcompr
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff_phi_1mphi
    REAL(KIND=8),                 INTENT(IN)    :: stab_mass
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: ff_loc, level_set_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: jj_loc
    REAL(KIND=8), DIMENSION(3)                  :: fcomprl
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: v_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: idxm
    INTEGER                                     :: i, m, l, ni, iglob, index
    REAL(KIND=8)                                :: fl, ray
    !#include "petsc/finclude/petsc.h"
    Vec                                         :: vect
    PetscErrorCode                              :: ierr

    CALL VecSet(vect, 0.d0, ierr)

    index = 0
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       ff_loc = ff(jj_loc)
       level_set_loc = level_set_ext(jj_loc)

       DO ni = 1, mesh%gauss%n_w
          i = mesh%jj(ni,m)
          iglob = LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          index = index + 1
          ray = 0
          DO ni = 1, mesh%gauss%n_w
             i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          ! Compute ff on gauss points + ff_gauss
          fl = (SUM(ff_loc(:)*mesh%gauss%ww(:,l)) + ff_gauss(index)+stab_mass*ff_phi_1mphi(index))*mesh%gauss%rj(l,m)*ray

          ! Compute compressive term on gauss points
          fcomprl=0.d0
          IF (type==1) THEN
             fcomprl(1) = fcompr(index,1)*mesh%gauss%rj(l,m)*ray
             fcomprl(2) = -mode*fcompr(index,4)*mesh%gauss%rj(l,m)
             fcomprl(3) = fcompr(index,5)*mesh%gauss%rj(l,m)*ray

          ELSE IF (type==2) THEN
             fcomprl(1) = fcompr(index,2)*mesh%gauss%rj(l,m)*ray
             fcomprl(2) = mode*fcompr(index,3)*mesh%gauss%rj(l,m)
             fcomprl(3) = fcompr(index,6)*mesh%gauss%rj(l,m)*ray
          ELSE
             CALL error_petsc('error in type while calling qs_00_level_set_gauss')
          END IF

          DO ni = 1,  mesh%gauss%n_w
             ! Add time derivative, advection and forcing term
             v_loc(ni) = v_loc(ni) + mesh%gauss%ww(ni,l) *fl
             ! Add fcompr=(c1_LES-nu_entro)*Grad(phi) + compression
             v_loc(ni) = v_loc(ni) + (fcomprl(1)*mesh%gauss%dw(1,ni,l,m) + fcomprl(2)*mesh%gauss%ww(ni,l) &
                  + fcomprl(3)*mesh%gauss%dw(2,ni,l,m))
          END DO
       END DO

       CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(vect,ierr)
    CALL VecAssemblyEnd(vect,ierr)
  END SUBROUTINE qs_00_level_set_gauss

END MODULE subroutine_level_set
