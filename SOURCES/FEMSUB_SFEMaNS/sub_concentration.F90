!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_concentration
  USE my_util
  USE boundary

  PUBLIC :: three_level_concentration
  PRIVATE
CONTAINS

  SUBROUTINE three_level_concentration(comm_one_d,time, conc_1_LA, dt, list_mode, &
       conc_mesh, concn_m1, concn, vel,& !chmp_mag,  &
       conc_diffusivity, my_par_cc, conc_list_dirichlet_sides, &
       conc_list_robin_sides, convection_coeff_conc_lhs, convection_coeff_conc_rhs, exterior_concentration, conc_per,&
       j_H_to_conc)
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
    USE chaine_caractere
    USE sub_plot
    USE st_matrix
    USE sft_parallele
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    REAL(KIND=8)                                        :: time, dt
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: list_mode
    TYPE(mesh_type),                INTENT(IN)          :: conc_mesh
    type(petsc_csr_LA)                                  :: conc_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: conc_per
    TYPE(solver_param),             INTENT(IN)          :: my_par_cc
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)       :: concn_m1, concn
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: conc_list_dirichlet_sides
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: conc_list_robin_sides
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)          :: conc_diffusivity
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)          :: convection_coeff_conc_lhs, convection_coeff_conc_rhs
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)          :: exterior_concentration
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)          :: j_H_to_conc
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN) :: vel
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: m_max_c
    INTEGER,     DIMENSION(:),   POINTER,          SAVE :: conc_js_D ! Dirichlet nodes
    INTEGER,                                       SAVE :: my_petscworld_rank
    REAL(KIND=8),                                  SAVE :: mass0, hmoy
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: conc_global_D ! MODIFICATION: axis BC
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: conc_mode_global_js_D ! MODIFICATION: axis BC
    !----------FIN SAVE--------------------------------------------------------------------

    !----------Declaration sans save-------------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: conc_1_ifrom
    INTEGER                                  :: i, m, n, l, index
    INTEGER                                  :: code, mode
    !allocations des variables locales
    REAL(KIND=8), DIMENSION(conc_mesh%np)                      :: ff
    REAL(KIND=8), DIMENSION(conc_mesh%np, 2)                   :: concn_p1
    REAL(KIND=8), DIMENSION(conc_mesh%gauss%l_G*conc_mesh%me,2, SIZE(list_mode)) :: ff_conv
    REAL(KIND=8)   ::tps, tps_tot, tps_cumul
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/
    REAL(KIND=8), DIMENSION(2,conc_mesh%gauss%l_G*conc_mesh%me)                :: rr_gauss
    INTEGER,      DIMENSION(conc_mesh%gauss%n_w)                               :: j_loc
    !Communicators for Petsc, in space and Fourier------------------------------
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: conc_mat
    Vec,                        SAVE :: cb_1, cb_2, cx_1, cx_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: conc_ksp
    !------------------------------END OF DECLARATION--------------------------------------

    IF (once) THEN

       once = .FALSE.

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)

       !-----CREATE PETSC VECTORS AND GHOSTS-----------------------------------------
       CALL create_my_ghost(conc_mesh,conc_1_LA,conc_1_ifrom)
       n = conc_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(conc_1_ifrom), conc_1_ifrom, cx_1, ierr)
       CALL VecGhostGetLocalForm(cx_1, cx_1_ghost, ierr)
       CALL VecDuplicate(cx_1, cb_1, ierr)
       CALL VecDuplicate(cx_1, cb_2, ierr)
       !------------------------------------------------------------------------------

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------

       !---------PREPARE pp_js_D ARRAY FOR TEMPERATURE--------------------------------------
       CALL dirichlet_nodes_parallel(conc_mesh, conc_list_dirichlet_sides, conc_js_D)
       CALL scalar_with_bc_glob_js_D(conc_mesh, list_mode, conc_1_LA, conc_js_D, conc_mode_global_js_D) ! MODIFICATION: axis BC
       ALLOCATE(conc_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(conc_global_D(i)%DRL(SIZE(conc_mode_global_js_D(i)%DIL)))
       END DO
       !------------------------------------------------------------------------------

       !--------------------------------------------------------------------------------
       hmoy = 0
       DO m = 1, conc_mesh%dom_me
          hmoy = hmoy + SQRT(SUM(conc_mesh%gauss%rj(:,m)))/2
       END DO
       hmoy =  hmoy/conc_mesh%dom_me
       mass0 = 0.d0
       DO i = 1, m_max_c
          mode = list_mode(i)
          IF (mode == 0)  THEN
             CALL mass_tot(comm_one_d(1),conc_mesh, concn(:,1,i), mass0)
          ENDIF
       END DO
       !--------------------------------------------------------------------------------

       !-------------ASSEMBLE CONCENTRATION MATRICES--------------------------------------------
       ALLOCATE(conc_mat(m_max_c),conc_ksp(m_max_c))

       DO i = 1, m_max_c
          mode = list_mode(i)

          !---CONCENTRATION MATRIX
          CALL create_local_petsc_matrix(comm_one_d(1), conc_1_LA, conc_mat(i), clean=.FALSE.)
          CALL qs_diff_mass_scal_M_conc(conc_mesh, conc_1_LA, conc_diffusivity, &
               1.5d0/dt, conc_list_robin_sides, convection_coeff_conc_lhs, zero, mode, conc_mat(i))
          IF (conc_per%n_bord/=0) THEN
             CALL periodic_matrix_petsc(conc_per%n_bord, conc_per%list, conc_per%perlist, conc_mat(i), conc_1_LA)
          END IF
          CALL Dirichlet_M_parallel(conc_mat(i),conc_mode_global_js_D(i)%DIL) ! MODIFICATION: axis BC
          CALL init_solver(my_par_cc,conc_ksp(i),conc_mat(i),comm_one_d(1),&
               solver=my_par_cc%solver,precond=my_par_cc%precond)
       ENDDO

    ENDIF
    tps_tot = user_time()
    tps_cumul = 0

    !===Compute rhs by FFT at Gauss points
    tps = user_time()
    CALL smb_ugradc_gauss_fft_par_conc(comm_one_d(2),conc_mesh,list_mode,vel,2*concn-concn_m1,ff_conv)
    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Tps fft vitesse', tps
    !------------CONSTRUCTION OF rr_gauss------------------
    index = 0
    DO m = 1, conc_mesh%me
       j_loc = conc_mesh%jj(:,m)
       DO l = 1, conc_mesh%gauss%l_G
          index = index + 1
          rr_gauss(1,index) = SUM(conc_mesh%rr(1,j_loc)*conc_mesh%gauss%ww(:,l))
          rr_gauss(2,index) = SUM(conc_mesh%rr(2,j_loc)*conc_mesh%gauss%ww(:,l))
       END DO
    END DO
    !------------DEBUT BOUCLE SUR LES MODES----------------
    DO i = 1, m_max_c
       mode = list_mode(i)

       !===RHS concentration
       ff = (2d0/dt)*concn(:,1,i) - (1d0/(2*dt))*concn_m1(:,1,i)
       CALL qs_00_gauss_conc (conc_mesh, conc_1_LA,ff,&
            -ff_conv(:,1,i), cb_1)

       ff = (2d0/dt)*concn(:,2,i) - (1d0/(2*dt))*concn_m1(:,2,i)
       CALL qs_00_gauss_conc (conc_mesh, conc_1_LA,ff,&
            -ff_conv(:,2,i), cb_2)

       IF (SIZE(j_H_to_conc,1).GT.1) THEN
          !===RHS Neumann BCs from current
          CALL qs_00_gauss_H_conc(conc_mesh, j_H_to_conc(:,:,i), conc_1_LA, cb_1, cb_2)
       ELSE
          !===RHS Robins BCs
          IF (mode == 0) THEN ! exterior temperature = constant
             !MODIFICATION: implementation of the term int_(partial Omega) h*Text*v
             CALL qs_00_gauss_surface_conc(conc_mesh, conc_1_LA, conc_list_robin_sides, &
                  convection_coeff_conc_rhs, exterior_concentration, cb_1)
          END IF
       END IF

       !===RHS periodicity
       IF (conc_per%n_bord/=0) THEN
          CALL periodic_rhs_petsc(conc_per%n_bord, conc_per%list, conc_per%perlist, cb_1, conc_1_LA)
          CALL periodic_rhs_petsc(conc_per%n_bord, conc_per%list, conc_per%perlist, cb_2, conc_1_LA)
       END IF
       !----------------------------------------------------------

       n = SIZE(conc_js_D) ! MODIFICATION: axis BC
       conc_global_D(i)%DRL(n+1:) = 0.d0
       conc_global_D(i)%DRL(1:n) = concentration_exact(1,conc_mesh%rr(:,conc_js_D), mode, time)
       CALL dirichlet_rhs(conc_mode_global_js_D(i)%DIL-1,conc_global_D(i)%DRL,cb_1)
       conc_global_D(i)%DRL(1:n) = concentration_exact(2,conc_mesh%rr(:,conc_js_D), mode, time)
       CALL dirichlet_rhs(conc_mode_global_js_D(i)%DIL-1,conc_global_D(i)%DRL,cb_2)

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps second membre vitesse', tps
       !-------------------------------------------------------------------------------------

       !--------------------INVERSION DES OPERATEURS--------------
       tps = user_time()
       !Solve system conc_c
       CALL solver(conc_ksp(i),cb_1,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,conc_1_LA,concn_p1(:,1))

       !Solve system conc_s
       CALL solver(conc_ksp(i),cb_2,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,conc_1_LA,concn_p1(:,2))
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solution des pb de vitesse', tps, 'for mode ', mode
       !-------------------------------------------------------------------------------------

       !---------------UPDATES-----------------------
       tps = user_time()

       IF (mode==0) THEN
          concn_p1 (:,2) = 0.d0
       END IF

       concn_m1(:,:,i) = concn(:,:,i)
       concn   (:,:,i) = concn_p1

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps  des updates', tps
       !-------------------------------------------------------------------------------------
    ENDDO

    tps_tot = user_time() - tps_tot
    !WRITE(*,'(A,2(f13.3,2x))') '  Tps boucle en temps Navier_stokes', tps_tot, tps_cumul
    !WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE three_level_concentration
  !============================================

  SUBROUTINE smb_ugradc_gauss_fft_par_conc(communicator,mesh,list_mode,V_in,c_in,c_out)
    !=================================
    USE Gauss_points
    USE sft_parallele
    USE chaine_caractere
    USE boundary
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: V_in, c_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: c_out
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: Gradc, W
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: Div, Cgauss, cint
    INTEGER,      DIMENSION(mesh%gauss%n_w)                  :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)   :: dw_loc
    INTEGER                                                  ::  m, l , i, mode, index, k
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)   :: Vs
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)   :: cs
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8)   :: ray, tps
    REAL(KIND=8), DIMENSION(3)                  :: temps
    INTEGER                                     :: code, m_max_pad, bloc_size, nb_procs
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

             !-----------------Gradient of c on Gauss points---------------------------
             !coeff sur les cosinus et sinus
             Gradc(index,1,i) = SUM(cs(:,1)*dw_loc(1,:))
             Gradc(index,2,i) = SUM(cs(:,2)*dw_loc(1,:))
             Gradc(index,3,i) =  mode/ray*SUM(cs(:,2)*ww(:,l))
             Gradc(index,4,i) = -mode/ray*SUM(cs(:,1)*ww(:,l))
             Gradc(index,5,i) = SUM(cs(:,1)*dw_loc(2,:))
             Gradc(index,6,i) = SUM(cs(:,2)*dw_loc(2,:))

             !-----------------c on Gauss points---------------------------------------
             Cgauss(index,1,i) = SUM(cs(:,1)*ww(:,l))
             Cgauss(index,2,i) = SUM(cs(:,2)*ww(:,l))
          ENDDO
       ENDDO
    END DO

    !tps = user_time() - tps
    !WRITE(*,*) ' Tps dans la grande boucle', tps
    !tps = user_time()
    temps = 0

    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(Gradc,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    CALL FFT_PAR_DOT_PROD_DCL(communicator, Gradc, W, c_out, nb_procs, bloc_size, m_max_pad, temps)
    bloc_size = SIZE(Div,1)/nb_procs+1
    CALL FFT_PAR_PROD_DCL(communicator, Div, Cgauss, cint, nb_procs, bloc_size, m_max_pad, temps)
    c_out = c_out + cint
    tps = user_time() - tps
    !WRITE(*,*) ' Tps dans FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Temps de Comm   ', temps(1)
    !write(*,*) ' Temps de Calc   ', temps(2)
    !write(*,*) ' Temps de Chan   ', temps(3)

  END SUBROUTINE smb_ugradc_gauss_fft_par_conc

  SUBROUTINE mass_tot(communicator,mesh,concn,RESLT)
    !===========================
    !moyenne
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(:)  ,   INTENT(IN)  :: concn
    REAL(KIND=8)                ,   INTENT(OUT) :: RESLT
    REAL(KIND=8)                                :: r_loc, r_out
    INTEGER ::  m, l , i , ni, code
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: j_loc
    REAL(KIND=8)   :: ray
    MPI_Comm                                    :: communicator
    r_loc = 0.d0

    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)
       DO l = 1, mesh%gauss%l_G
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = j_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          r_loc = r_loc +  SUM(concn(j_loc(:))*mesh%gauss%ww(:,l))*ray*mesh%gauss%rj(l,m)
       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(r_loc,r_out,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    RESLT = r_out

  END SUBROUTINE mass_tot


END MODULE subroutine_concentration
