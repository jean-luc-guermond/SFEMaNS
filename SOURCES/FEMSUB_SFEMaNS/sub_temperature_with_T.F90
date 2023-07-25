!
!Authors Jean-Luc Guermond, Raphael Laguerre, Caroline Nore, Copyrights 2005
!Revised for PETSC, Jean-Luc Guermond, Franky Luddens, January 2011
!
MODULE subroutine_temperature_with_T
  USE my_util
  USE boundary

  PUBLIC :: three_level_temperature_with_T
  PRIVATE
CONTAINS

  SUBROUTINE three_level_temperature_with_T(comm_one_d,time, temp_1_LA, dt, list_mode, &
       temp_mesh, tempn_m1, tempn, vel_field, mag_field, pdt_H_field, vol_heat_capacity, &
       temp_diffusivity, my_par_cc, temp_list_dirichlet_sides, &
       temp_list_robin_sides, convection_coeff, exterior_temperature, temp_per)
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
    TYPE(mesh_type),                INTENT(IN)          :: temp_mesh
    type(petsc_csr_LA)                                  :: temp_1_LA
    TYPE(periodic_type),            INTENT(IN)          :: temp_per
    TYPE(solver_param),             INTENT(IN)          :: my_par_cc
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)       :: tempn_m1, tempn
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: temp_list_dirichlet_sides
    INTEGER,      DIMENSION(:),     INTENT(IN)          :: temp_list_robin_sides
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)          :: vol_heat_capacity, temp_diffusivity
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)          :: convection_coeff, exterior_temperature
    REAL(KIND=8), DIMENSION(:,:,:),          INTENT(IN) :: vel_field, mag_field, pdt_H_field
    LOGICAL,                                       SAVE :: once = .TRUE.
    INTEGER,                                       SAVE :: m_max_c
    INTEGER,     DIMENSION(:),   POINTER,          SAVE :: temp_js_D ! Dirichlet nodes
    INTEGER,                                       SAVE :: my_petscworld_rank
    REAL(KIND=8),                                  SAVE :: mass0, hmoy
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE :: temp_global_D ! axis BC
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE :: temp_mode_global_js_D ! axis BC
    !----------END SAVE--------------------------------------------------------------------

    !----------Declaration without save----------------------------------------------------
    INTEGER,          POINTER, DIMENSION(:)  :: temp_1_ifrom
    INTEGER                                  :: i, m, n, l, index
    INTEGER                                  :: code, mode
    !Allocation of local variables
    REAL(KIND=8), DIMENSION(temp_mesh%np)                      :: ff
    REAL(KIND=8), DIMENSION(temp_mesh%np, 2)                   :: tempn_p1
    REAL(KIND=8), DIMENSION(temp_mesh%gauss%l_G*temp_mesh%me,2, SIZE(list_mode)) :: ff_conv, pyromag_term
    REAL(KIND=8)   ::tps, tps_tot, tps_cumul
    REAL(KIND=8) :: one, zero, three
    DATA zero, one, three/0.d0,1.d0,3.d0/
    REAL(KIND=8), DIMENSION(2,temp_mesh%gauss%l_G*temp_mesh%me)                  :: rr_gauss
    INTEGER,      DIMENSION(temp_mesh%gauss%n_w)                               :: j_loc
    !Communicators for Petsc, in space and Fourier----------------------------------------
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: temp_mat
    Vec,                        SAVE :: cb_1, cb_2, cx_1, cx_1_ghost
    KSP, DIMENSION(:), POINTER, SAVE :: temp_ksp
    !------------------------------END OF DECLARATION--------------------------------------

    IF (once) THEN

       once = .FALSE.

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)

       !-----CREATE PETSC VECTORS AND GHOSTS-------------------------------------------
       CALL create_my_ghost(temp_mesh,temp_1_LA,temp_1_ifrom)
       n = temp_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, &
            PETSC_DETERMINE, SIZE(temp_1_ifrom), temp_1_ifrom, cx_1, ierr)
       CALL VecGhostGetLocalForm(cx_1, cx_1_ghost, ierr)
       CALL VecDuplicate(cx_1, cb_1, ierr)
       CALL VecDuplicate(cx_1, cb_2, ierr)
       !--------------------------------------------------------------------------------

       !-------------DIMENSIONS---------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !--------------------------------------------------------------------------------

       !---------PREPARE pp_js_D ARRAY FOR TEMPERATURE----------------------------------
       CALL dirichlet_nodes_parallel(temp_mesh, temp_list_dirichlet_sides, temp_js_D)
       CALL scalar_with_bc_glob_js_D(temp_mesh, list_mode, temp_1_LA, temp_js_D, temp_mode_global_js_D) ! MODIFICATION: axis BC
       ALLOCATE(temp_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(temp_global_D(i)%DRL(SIZE(temp_mode_global_js_D(i)%DIL)))
       END DO
       !--------------------------------------------------------------------------------

       !------------------------------------------------------------------------------
       hmoy = 0
       DO m = 1, temp_mesh%dom_me
          hmoy = hmoy + SQRT(SUM(temp_mesh%gauss%rj(:,m)))/2
       END DO
       hmoy =  hmoy/temp_mesh%dom_me
       mass0 = 0.d0
       DO i = 1, m_max_c
          mode = list_mode(i)
          IF (mode == 0)  THEN
             CALL mass_tot(comm_one_d(1),temp_mesh, tempn(:,1,i), mass0)
          ENDIF
       END DO
       !--------------------------------------------------------------------------------

       !-------------ASSEMBLE TEMPERATURE MATRICES--------------------------------------
       ALLOCATE(temp_mat(m_max_c),temp_ksp(m_max_c))

       DO i = 1, m_max_c
          mode = list_mode(i)

          !---TEMPERATURE MATRIX
          CALL create_local_petsc_matrix(comm_one_d(1), temp_1_LA, temp_mat(i), clean=.FALSE.)
          CALL qs_diff_mass_scal_M_variant(temp_mesh, temp_1_LA, vol_heat_capacity, temp_diffusivity, &
               1.5d0/dt, temp_list_robin_sides, convection_coeff, zero, mode, temp_mat(i))

          IF (temp_per%n_bord/=0) THEN
             CALL periodic_matrix_petsc(temp_per%n_bord, temp_per%list, temp_per%perlist, temp_mat(i), temp_1_LA)
          END IF

          CALL Dirichlet_M_parallel(temp_mat(i),temp_mode_global_js_D(i)%DIL)

          CALL init_solver(my_par_cc,temp_ksp(i),temp_mat(i),comm_one_d(1),&
               solver=my_par_cc%solver,precond=my_par_cc%precond)
       ENDDO

    ENDIF
    tps_tot = user_time()
    tps_cumul = 0

    !===Compute convection term at Gauss points
    tps = user_time()
    CALL smb_ugradc_gauss_fft_par(comm_one_d(2),temp_mesh,list_mode,vol_heat_capacity,vel_field,2*tempn-tempn_m1,ff_conv)

    !===Compute pyromagnetic term at Gauss points if fhd
    IF (inputs%type_pb=='fhd') THEN
       CALL smb_pyromag_gauss_fft_par(comm_one_d(2),temp_mesh,list_mode,2*tempn-tempn_m1,vel_field,mag_field, &
            pdt_H_field,pyromag_term)
       ff_conv = ff_conv + pyromag_term
    END IF
    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Time FFT in temperature equation', tps

    !===Compute radius at Gauss points
    index = 0
    DO m = 1, temp_mesh%me
       j_loc = temp_mesh%jj(:,m)
       DO l = 1, temp_mesh%gauss%l_G
          index = index + 1
          rr_gauss(1,index) = SUM(temp_mesh%rr(1,j_loc)*temp_mesh%gauss%ww(:,l))
          rr_gauss(2,index) = SUM(temp_mesh%rr(2,j_loc)*temp_mesh%gauss%ww(:,l))
       END DO
    END DO

    !------------BEGIN LOOP ON FOURIER MODES---------------
    DO i = 1, m_max_c
       mode = list_mode(i)

       !===RHS temperature
       ff = (2.d0/dt)*tempn(:,1,i) - (1.d0/(2*dt))*tempn_m1(:,1,i)
       CALL qs_00_gauss (temp_mesh, temp_1_LA, vol_heat_capacity, ff, &
            -ff_conv(:,1,i) + source_in_temperature(1, rr_gauss, mode, time), cb_1)

       ff = (2.d0/dt)*tempn(:,2,i) - (1.d0/(2*dt))*tempn_m1(:,2,i)
       CALL qs_00_gauss (temp_mesh, temp_1_LA, vol_heat_capacity, ff, &
            -ff_conv(:,2,i) + source_in_temperature(2, rr_gauss, mode, time), cb_2)

       !===RHS Robins BCs
       IF (mode == 0) THEN ! exterior temperature = constant
          CALL qs_00_gauss_surface(temp_mesh, temp_1_LA, temp_list_robin_sides, convection_coeff, exterior_temperature, cb_1)
       END IF

       !===RHS periodicity
       IF (temp_per%n_bord/=0) THEN
          CALL periodic_rhs_petsc(temp_per%n_bord, temp_per%list, temp_per%perlist, cb_1, temp_1_LA)
          CALL periodic_rhs_petsc(temp_per%n_bord, temp_per%list, temp_per%perlist, cb_2, temp_1_LA)
       END IF
       !------------------------------------------------------

       n = SIZE(temp_js_D)
       temp_global_D(i)%DRL(n+1:) = 0.d0
       temp_global_D(i)%DRL(1:n) = temperature_exact(1,temp_mesh%rr(:,temp_js_D), mode, time)
       CALL dirichlet_rhs(temp_mode_global_js_D(i)%DIL-1,temp_global_D(i)%DRL,cb_1)
       temp_global_D(i)%DRL(1:n) = temperature_exact(2,temp_mesh%rr(:,temp_js_D), mode, time)
       CALL dirichlet_rhs(temp_mode_global_js_D(i)%DIL-1,temp_global_D(i)%DRL,cb_2)

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Time computing RHS temperature problem', tps
       !------------------------------------------------------

       !------------INVERTING OPERATORS-----------------------
       tps = user_time()
       !Solve system temp_cosine
       CALL solver(temp_ksp(i),cb_1,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,temp_1_LA,tempn_p1(:,1))

       !Solve system temp_sine
       CALL solver(temp_ksp(i),cb_2,cx_1,reinit=.FALSE.,verbose=my_par_cc%verbose)
       CALL VecGhostUpdateBegin(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL VecGhostUpdateEnd(cx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       CALL extract(cx_1_ghost,1,1,temp_1_LA,tempn_p1(:,2))
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Time inverting temperature problem', tps, 'for mode ', mode
       !------------------------------------------------------

       !------------UPDATES-----------------------------------
       tps = user_time()

       !Force Sine Fourier coefficient to zero for mode=0
       IF (mode==0) THEN
          tempn_p1 (:,2) = 0.d0
       END IF

       tempn_m1(:,:,i) = tempn(:,:,i)
       tempn   (:,:,i) = tempn_p1

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps  des updates', tps
       !------------------------------------------------------
    ENDDO

    tps_tot = user_time() - tps_tot
    !WRITE(*,'(A,2(f13.3,2x))') ' Time for loop in Temperature', tps_tot, tps_cumul
    !WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE three_level_temperature_with_T
  !============================================

  SUBROUTINE smb_ugradc_gauss_fft_par(communicator,mesh,list_mode,heat_capa_in,V_in,c_in,c_out)
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
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: heat_capa_in
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
    REAL(KIND=8)                                :: ray, tps
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

             !------------Velocity at Gauss points------------------
             W(index,1,i) = SUM(Vs(:,1)*ww(:,l))
             W(index,3,i) = SUM(Vs(:,3)*ww(:,l))
             W(index,5,i) = SUM(Vs(:,5)*ww(:,l))

             W(index,2,i) = SUM(Vs(:,2)*ww(:,l))
             W(index,4,i) = SUM(Vs(:,4)*ww(:,l))
             W(index,6,i) = SUM(Vs(:,6)*ww(:,l))

             !------------Divergence Vecocity at Gauss points-------
             Div(index,1,i) = SUM(Vs(:,1)*dw_loc(1,:)) + SUM(Vs(:,1)*ww(:,l))/ray &
                  + (mode/ray)*SUM(Vs(:,4)*ww(:,l)) +  SUM(Vs(:,5)*dw_loc(2,:))
             Div(index,2,i) = SUM(Vs(:,2)*dw_loc(1,:)) + SUM(Vs(:,2)*ww(:,l))/ray &
                  - (mode/ray)*SUM(Vs(:,3)*ww(:,l)) +  SUM(Vs(:,6)*dw_loc(2,:))

             !------------Temperature gradient at Gauss points------
             Gradc(index,1,i) = SUM(cs(:,1)*dw_loc(1,:))
             Gradc(index,2,i) = SUM(cs(:,2)*dw_loc(1,:))
             Gradc(index,3,i) =  mode/ray*SUM(cs(:,2)*ww(:,l))
             Gradc(index,4,i) = -mode/ray*SUM(cs(:,1)*ww(:,l))
             Gradc(index,5,i) = SUM(cs(:,1)*dw_loc(2,:))
             Gradc(index,6,i) = SUM(cs(:,2)*dw_loc(2,:))

             Gradc(index,:,i) = heat_capa_in(m) * Gradc(index,:,i)

             !------------Temperature at Gauss points---------------
             Cgauss(index,1,i) = SUM(cs(:,1)*ww(:,l))
             Cgauss(index,2,i) = SUM(cs(:,2)*ww(:,l))

             Cgauss(index,:,i) = heat_capa_in(m) * Cgauss(index,:,i)

          ENDDO
       ENDDO
    END DO

    !tps = user_time() - tps
    !WRITE(*,*) ' Time in big loop', tps
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
    !WRITE(*,*) ' Time in FFT_PAR_PROD_VECT', tps
    !write(*,*) ' Communication time   ', temps(1)
    !write(*,*) ' Computation time   ', temps(2)
    !write(*,*) ' Change Format Time   ', temps(3)
  END SUBROUTINE smb_ugradc_gauss_fft_par

  SUBROUTINE mass_tot(communicator,mesh,tempn,RESLT)
    !===========================
    !moyenne
    USE def_type_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                             :: mesh
    REAL(KIND=8), DIMENSION(:)  ,   INTENT(IN)  :: tempn
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

          r_loc = r_loc +  SUM(tempn(j_loc(:))*mesh%gauss%ww(:,l))*ray*mesh%gauss%rj(l,m)

       ENDDO
    ENDDO

    CALL MPI_ALLREDUCE(r_loc,r_out,1,MPI_DOUBLE_PRECISION, MPI_SUM, communicator, code)
    RESLT = r_out

  END SUBROUTINE mass_tot

  SUBROUTINE smb_pyromag_gauss_fft_par(communicator,mesh,list_mode,scal_in,vect_1_in,vect_2_in,vect_3_in,scal_out)
    USE Gauss_points
    USE sft_parallele
    USE def_type_mesh
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN) :: mesh
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: scal_in, vect_1_in, vect_2_in, vect_3_in
    REAL(KIND=8), DIMENSION(:,:,:)             :: scal_out
    REAL(KIND=8), DIMENSION(mesh%np,2,SIZE(list_mode)) :: T_dchi_dT_coeff, vect_2_in_square, v2_dot_v3
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: T_dchi_dT_coeff_gauss, pdtH2, ugradH2, DH2_Dt
    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: vect_1_in_gauss, grad_vect_2_in_square
    INTEGER                                                           :: i, mode, index, m, k, l
    INTEGER,      DIMENSION(:,:), POINTER                             :: jj
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc
    REAL(KIND=8)                                                      :: rad
    !===FOR FFT
    INTEGER                                     :: nb_procs, bloc_size, m_max_pad, code
    MPI_Comm                                    :: communicator

    ! We compute the pyromagnetic coefficient term: T_dchi_dT_coeff(T) * 1/2 * D(H**2)/Dt

    !===nb_procs and m_max_pad calculus for FFT
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2

    !===T_dchi_dT_coeff(T) on Gauss nodes computation
    T_dchi_dT_coeff = scal_in
    bloc_size = SIZE(T_dchi_dT_coeff,1)/nb_procs+1
    CALL FFT_PAR_SCAL_FUNCT(communicator, T_dchi_dT_coeff, T_dchi_dT_coeff_law, nb_procs, bloc_size, m_max_pad)
    CALL gauss(mesh)
    jj => mesh%jj
    DO i = 1, SIZE(list_mode)
       index = 0 ! global index of Gauss node
       DO m = 1, mesh%dom_me
          IF (MINVAL(ABS(mesh%i_d(m) - inputs%list_dom_ns)) == 0) THEN ! dchi/dT non null inside the (ferro)fluid domain only
             DO l=1, l_G
                index = index + 1
                DO k = 1, 2
                   T_dchi_dT_coeff_gauss(index,k,i) = SUM(T_dchi_dT_coeff(jj(:,m),k,i)*ww(:,l))
                END DO
             END DO
          ELSE
             T_dchi_dT_coeff_gauss(index+1:index+l_G,:,i) = 0.d0
             index = index + l_G
          END IF
       END DO
    END DO

    !===pdt(H**2) on Gauss nodes computation
    IF (inputs%if_steady_current_fhd) THEN
       pdtH2 = 0.d0
    ELSE
       bloc_size = SIZE(vect_2_in,1)/nb_procs+1
       CALL FFT_PAR_DOT_PROD_DCL(communicator, vect_2_in, vect_3_in, v2_dot_v3, nb_procs, bloc_size, m_max_pad)
       DO i = 1, SIZE(list_mode)
          index = 0 ! global index of Gauss node
          DO m = 1, mesh%dom_me
             IF (MINVAL(ABS(mesh%i_d(m) - inputs%list_dom_ns)) == 0) THEN ! dchi/dT non null inside the (ferro)fluid domain only
                DO l=1, l_G
                   index = index + 1
                   DO k = 1, 2
                      pdtH2(index,k,i) = 2*SUM(v2_dot_v3(jj(:,m),k,i)*ww(:,l))
                   END DO
                END DO
             ELSE
                pdtH2(index+1:index+l_G,:,i) = 0.d0
                index = index + l_G
             END IF
          END DO
       END DO
    END IF

    !===u on Gauss nodes computation
    DO i = 1, SIZE(list_mode)
       index = 0 ! global index of Gauss node
       DO m = 1, mesh%dom_me
          IF (MINVAL(ABS(mesh%i_d(m) - inputs%list_dom_ns)) == 0) THEN ! dchi/dT non null inside the (ferro)fluid domain only
             DO l=1, l_G
                index = index + 1
                DO k = 1, 6
                   vect_1_in_gauss(index,k,i) = SUM(vect_1_in(jj(:,m),k,i)*ww(:,l))
                END DO
             END DO
          ELSE
             vect_1_in_gauss(index+1:index+l_G,:,i) = 0.d0
             index = index + l_G
          END IF
       END DO
    END DO

    !===grad(H**2) on Gauss nodes computation
    bloc_size = SIZE(vect_2_in,1)/nb_procs+1
    CALL FFT_PAR_DOT_PROD_DCL(communicator, vect_2_in, vect_2_in, vect_2_in_square, nb_procs, bloc_size, m_max_pad)
    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0
       DO m = 1, mesh%dom_me
          IF (MINVAL(ABS(mesh%i_d(m) - inputs%list_dom_ns)) == 0) THEN ! dchi/dT non null inside the (ferro)fluid domain only
             DO l=1, l_G
                index = index + 1
                dw_loc = dw(:,:,l,m)
                rad = SUM(mesh%rr(1,jj(:,m))*ww(:,l)) ! radius of Gauss node
                grad_vect_2_in_square(index,1,i) = SUM(vect_2_in_square(jj(:,m),1,i)*dw_loc(1,:))
                grad_vect_2_in_square(index,2,i) = SUM(vect_2_in_square(jj(:,m),2,i)*dw_loc(1,:))
                grad_vect_2_in_square(index,3,i) = mode/rad * SUM(vect_2_in_square(jj(:,m),2,i)*ww(:,l))
                grad_vect_2_in_square(index,4,i) = - mode/rad * SUM(vect_2_in_square(jj(:,m),1,i)*ww(:,l))
                grad_vect_2_in_square(index,5,i) = SUM(vect_2_in_square(jj(:,m),1,i)*dw_loc(2,:))
                grad_vect_2_in_square(index,6,i) = SUM(vect_2_in_square(jj(:,m),2,i)*dw_loc(2,:))
             END DO
          ELSE
             grad_vect_2_in_square(index+1:index+l_G,:,i) = 0.d0
             index = index + l_G
          END IF
       END DO
    END DO

    !===u.grad(H**2) on Gauss nodes computation
    bloc_size = SIZE(vect_1_in_gauss,1)/nb_procs+1
    CALL FFT_PAR_DOT_PROD_DCL(communicator, vect_1_in_gauss, grad_vect_2_in_square, ugradH2, nb_procs, bloc_size, m_max_pad)

    !===D(H**2)/Dt on Gauss nodes computation
    DH2_Dt = pdtH2 + ugradH2

    !===T_dchi_dT_coeff(T) * 1/2 * D(H**2)/Dt
    bloc_size = SIZE(T_dchi_dT_coeff_gauss,1)/nb_procs+1
    CALL FFT_PAR_PROD_DCL(communicator, T_dchi_dT_coeff_gauss, 0.5*DH2_Dt, scal_out, nb_procs, bloc_size, m_max_pad)

  END SUBROUTINE smb_pyromag_gauss_fft_par

END MODULE subroutine_temperature_with_T
