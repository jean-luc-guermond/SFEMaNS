!$
!Authors Jean-Luc Guermond, Raphael Laguerre, Copyrights 2005
!Revised June 2008, Jean-Luc Guermond
!Revised Jan/Feb 2009, Caroline Nore, Jean-Luc Guermond, Franky Luddens
!
MODULE update_maxwell_mxs_with_H

  PUBLIC:: maxwell_mxs_with_H
  PRIVATE
  REAL(KIND=8), PARAMETER, PRIVATE :: alpha=0.6d0
  INTEGER,  DIMENSION(:),  ALLOCATABLE :: Neumann_bdy_H_sides
  INTEGER,  DIMENSION(:),  ALLOCATABLE :: Neumann_bdy_pmag_sides
  INTEGER,  DIMENSION(:),  ALLOCATABLE :: Neumann_bdy_phi_sides
  !INTEGER,  DIMENSION(:),  ALLOCATABLE :: Dirichlet_bdy_H_sides
CONTAINS

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

  SUBROUTINE maxwell_mxs_with_H(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, stab_jump_h, sigma_in, &
       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt_in, Rem, list_mode, &
       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns_in, jj_v_to_H, conc_to_H)
    USE def_type_mesh
    USE chaine_caractere
    USE solve_petsc
    USE boundary
    USE tn_axi
    USE prep_maill
    USE Dir_nodes_petsc
    USE st_matrix
    USE Dir_nodes
    USE my_util
    USE sft_parallele
    USE sub_plot
    USE periodic
    USE input_data
    USE verbose
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscvec.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: H_mesh, phi_mesh, pmag_mesh
    TYPE(interface_type),           INTENT(IN)     :: interface_H_phi, interface_H_mu
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: vel  
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Hn, Hn1
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Bn, Bn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: phin, phin1
    REAL(KIND=8), DIMENSION(3),     INTENT(IN)     :: stab_in 
    REAL(KIND=8),                   INTENT(IN)     :: stab_jump_h
    REAL(KIND=8),                   INTENT(IN)     :: R_fourier
    INTEGER,                        INTENT(IN)     :: index_fourier
    REAL(KIND=8),                   INTENT(IN)     :: mu_phi, time, dt_in, Rem
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: sigma_in, mu_H_field
    TYPE(periodic_type),            INTENT(IN)     :: H_phi_per
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: one_over_sigma_ns_in
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: conc_to_H
    INTEGER,      DIMENSION(:)    , INTENT(IN)     :: jj_v_to_H
    TYPE(petsc_csr_LA)                             :: LA_H, LA_pmag, LA_phi, LA_mhd
    REAL(KIND=8),                                  SAVE  :: dt
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE,       SAVE  :: sigma_ns_bar
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE,       SAVE  :: sigma_np
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE,       SAVE  :: sigma
    REAL(KIND=8),                                  SAVE  :: sigma_min
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE  :: H_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE  :: H_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE  :: pmag_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE  :: pmag_global_D
    TYPE(dyn_int_line), DIMENSION(:), POINTER,     SAVE  :: phi_mode_global_js_D
    TYPE(dyn_real_line),DIMENSION(:), ALLOCATABLE, SAVE  :: phi_global_D
    INTEGER,     DIMENSION(:),   POINTER,          SAVE  :: pmag_js_D, phi_js_D
    INTEGER,     DIMENSION(:),  ALLOCATABLE,       SAVE  :: Dirichlet_bdy_H_sides
    LOGICAL,                                       SAVE  :: once=.TRUE.
    INTEGER,                                       SAVE  :: m_max_c
    REAL(KIND=8), DIMENSION(3),                    SAVE  :: stab 
    INTEGER,                                       SAVE  :: my_petscworld_rank
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: sigma_curl_gauss_bdy
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: J_over_sigma_gauss_bdy
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: sigma_curl_gauss_inter_mu
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: J_over_sigma_gauss_inter_mu
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: sigma_tot_gauss_Neumann
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:,:),  SAVE  :: minus_grad_times_der_pot_rhoLi
    REAL(KIND=8), ALLOCATABLE , DIMENSION(:,:),    SAVE  :: sigma_nj_m
    REAL(KIND=8), DIMENSION(H_mesh%gauss%l_G*H_mesh%me,6,SIZE(list_mode))  :: sigma_curl_gauss
    REAL(KIND=8), DIMENSION(H_mesh%gauss%l_G*H_mesh%me,6,SIZE(list_mode))  :: J_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(SIZE(Hn,1),6,SIZE(Hn,3))                       :: H_ns
!!$    REAL(KIND=8), DIMENSION(SIZE(Hn,1),2,SIZE(Hn,3))                       :: one_over_sigma_tot
    REAL(KIND=8), DIMENSION(H_mesh%np, 6, SIZE(list_mode))  :: ff
    REAL(KIND=8), DIMENSION(H_mesh%np, 2, SIZE(list_mode))  :: rhoLi_node, der_pot_rhoLi_node 
    LOGICAL, ALLOCATABLE, DIMENSION(:)                   :: Dir_pmag 
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                 :: rhs_H
    REAL(KIND=8), DIMENSION(phi_mesh%np,2)               :: rhs_phi
    REAL(KIND=8), DIMENSION(SIZE(Hn,1),6,SIZE(Hn,3))     :: NL, H_ext, B_ext
    REAL(KIND=8), DIMENSION(3)                           :: temps_par
    INTEGER,          POINTER, DIMENSION(:)              :: H_ifrom, pmag_ifrom, phi_ifrom, H_p_phi_ifrom
    REAL(KIND=8), DIMENSION(phi_mesh%np, 2)              :: phin_p1
    REAL(KIND=8), DIMENSION(H_mesh%np, 6)                :: Hn_p1
    LOGICAL,      DIMENSION(H_mesh%mes)                  :: virgin1
    LOGICAL,      DIMENSION(phi_mesh%mes)                :: virgin2
    INTEGER        :: mode, k, i, n, m, ms, code, nj, j, count
    INTEGER                                              :: nb_procs, bloc_size, m_max_pad
    REAL(KIND=8)   ::  tps,  nr_vel, tps_tot, tps_cumul, norm
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one_and_half
    DATA one_and_half/1.5d0/
    !April 17th, 2008, JLG
    PetscErrorCode                   :: ierr
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d
    Mat, DIMENSION(:), POINTER, SAVE :: H_p_phi_mat1, H_p_phi_mat2
    Mat                              :: tampon1, tampon2, precond1, precond2
    KSP, DIMENSION(:), POINTER, SAVE :: H_p_phi_ksp1, H_p_phi_ksp2
    Vec,                        SAVE :: vx_1, vb_1, vx_1_ghost, vx_2, vb_2, vx_2_ghost
    !------------------------------END OF DECLARATION--------------------------------------

    IF (once) THEN

       once = .FALSE.

!!$       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
!!$          per = .TRUE.
!!$       ELSE
!!$          per = .FALSE.
!!$       END IF

       CALL MPI_COMM_RANK(PETSC_COMM_WORLD,my_petscworld_rank,code)

       CALL create_my_ghost(H_mesh,LA_H,H_ifrom)
       CALL create_my_ghost(pmag_mesh,LA_pmag,pmag_ifrom)
       CALL create_my_ghost(phi_mesh,LA_phi,phi_ifrom)

!!$       !===Test if fhd ! MODIFICATION: fhd => CROSS H = js
!!$       IF (inputs%type_pb=='fhd') THEN
!!$          dt = 1.d20*dt_in
!!$       ELSE
!!$          dt = dt_in
!!$       END IF
       !===Modify dt for mxs/mhs problem
       dt = 1.d20*dt_in

       n  = SIZE(H_ifrom)+SIZE(pmag_ifrom)+SIZE(phi_ifrom)
       ALLOCATE(H_p_phi_ifrom(n))
       IF (SIZE(H_ifrom)/=0) THEN
          H_p_phi_ifrom(1:SIZE(H_ifrom)) = H_ifrom
       END IF
       IF (SIZE(pmag_ifrom)/=0) THEN
          H_p_phi_ifrom(SIZE(H_ifrom)+1:SIZE(H_ifrom)+SIZE(pmag_ifrom)) = pmag_ifrom
       END IF
       IF (SIZE(phi_ifrom)/=0) THEN
          H_p_phi_ifrom(SIZE(H_ifrom)+SIZE(pmag_ifrom)+1:)=phi_ifrom
       END IF

       n = 3*H_mesh%dom_np + pmag_mesh%dom_np + phi_mesh%dom_np
       CALL VecCreateGhost(comm_one_d(1), n, & 
            PETSC_DETERMINE, SIZE(H_p_phi_ifrom), H_p_phi_ifrom, vx_1, ierr)
       CALL VecGhostGetLocalForm(vx_1, vx_1_ghost, ierr)
       CALL VecDuplicate(vx_1, vb_1, ierr)
       CALL VecCreateGhost(comm_one_d(1), n, & 
            PETSC_DETERMINE, SIZE(H_p_phi_ifrom), H_p_phi_ifrom, vx_2, ierr)
       CALL VecGhostGetLocalForm(vx_2, vx_2_ghost, ierr)
       CALL VecDuplicate(vx_2, vb_2, ierr)
       !------------------------------------------------------------------------------

       !-------------RESCALING DE SIGMA-----------------------------------------------
       ALLOCATE(sigma(SIZE(sigma_in)))
       sigma = sigma_in * Rem

       CALL MPI_ALLREDUCE(MINVAL(sigma),sigma_min,1,MPI_DOUBLE_PRECISION, MPI_MIN,comm_one_d(1), ierr)
       !------------------------------------------------------------------------------

       !-------------RESCALING DE STAB------------------------------------------------

       stab = stab_in / Rem ! MODIFICATION: stab_in = data coefficients, normalization by Rm

!!$       !MARCH, 2010
!!$       IF (inputs%type_pb=='mhd') THEN
!!$          ! FL, 31/03/11
!!$          !stab = stab_in*(1/MINVAL(sigma)+1.d0) 
!!$          stab = stab_in*(1/sigma_min+1.d0) 
!!$          ! FL, 31/03/11
!!$	  ! Velocity assume to be used as reference scale
!!$!LC 2016/02/29
!!$          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
!!$             stab = stab_in*(1/(MINVAL(inputs%sigma_fluid)*Rem)+1.d0)
!!$          END IF
!!$!LC 2016/02/29
!!$       ELSE
!!$          nr_vel = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, vel) 
!!$
!!$          IF (nr_vel .LE. 1.d-10) THEN
!!$             ! FL, 31/03/11
!!$             !stab = stab_in*(1/MINVAL(sigma))
!!$             stab = stab_in*(1/sigma_min)
!!$             ! FL, 31/03/11
!!$             !WRITE(*,*) 'case 1, stab = ',stab
!!$          ELSE
!!$             ! FL, 31/03/11
!!$             !stab = stab_in*(1/MINVAL(sigma)+1.d0)
!!$             stab = stab_in*(1/sigma_min+1.d0)
!!$             ! FL, 31/03/11
!!$             !WRITE(*,*) 'case 2, stab = ',stab
!!$          ENDIF
!!$          ! Velocity could be zero in case of Ohmic decay
!!$       END IF
!!$       WRITE(*,*) 'stab = ',stab
!!$       !MARCH, 2010
       !------------------------------------------------------------------------------

       !-------------DIMENSIONS-------------------------------------------------------
       m_max_c = SIZE(list_mode)
       !------------------------------------------------------------------------------

       !------------SIGMA IF LEVEL SET------------------------------------------------
       ALLOCATE(sigma_nj_m(H_mesh%gauss%n_w,H_mesh%me))
       IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
          ALLOCATE(sigma_ns_bar(SIZE(Hn,1)))
          sigma_ns_bar = sigma_bar_in_fourier_space(H_mesh)*Rem

          !===check if j=H_mesh%jj(nj,m) is in ns domain or not and define sigma in consequence
          DO m = 1, H_mesh%me
             DO nj = 1, H_mesh%gauss%n_w
                j = H_mesh%jj(nj,m)
                IF (jj_v_to_H(j) == -1) THEN
                   sigma_nj_m(nj,m) = sigma(m)
                ELSE
                   sigma_nj_m(nj,m) = sigma_ns_bar(j)
                END IF
             END DO
          END DO
       ELSE
          DO m = 1, H_mesh%me
             sigma_nj_m(:,m) = sigma(m)
          END DO
       END IF

       ALLOCATE(sigma_np(SIZE(Hn,1)))
       sigma_np = 0.d0
       DO m = 1, H_mesh%me
          DO nj = 1, H_mesh%gauss%n_w
             sigma_np(H_mesh%jj(nj,m)) = sigma_nj_m(nj,m)
          END DO
       END DO
       !------------------------------------------------------------------------------

       !---------------BOUNDARY CONDITIONS FOR pmag-----------------------------------
       ! Creation of Dirichlet boundary conditions for the magnetic pressure
       ! Only on the boundary that is not periodic...
       ALLOCATE (Dir_pmag(MAXVAL(pmag_mesh%sides)))
       Dir_pmag = .FALSE.
       DO ms = 1, SIZE(Dir_pmag)
          IF (MINVAL(ABS(inputs%list_dirichlet_sides_H-ms)) == 0) THEN
             Dir_pmag(ms) = .TRUE.
          END IF
          IF (MINVAL(ABS(inputs%list_inter_H_phi-ms)) == 0) THEN
             Dir_pmag(ms) = .TRUE.
          END IF
       END DO

       CALL dirichlet_nodes(pmag_mesh%jjs, pmag_mesh%sides, Dir_pmag, pmag_js_D)
       DEALLOCATE(Dir_pmag)
       !ALLOCATE(pmag_bv_D(SIZE(pmag_js_D)))
       !pmag_bv_D = 0.d0
       CALL scalar_with_bc_glob_js_D(pmag_mesh, list_mode, LA_pmag, pmag_js_D, pmag_mode_global_js_D)
       ALLOCATE(pmag_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(pmag_global_D(i)%DRL(SIZE(pmag_mode_global_js_D(i)%DIL)))
          pmag_global_D(i)%DRL = 0.d0
       END DO
       ! End creation of Dirichlet boundary conditions for the magnetic pressure

       !===JLG+CN July 2017
       !===Neuman BC for H
       virgin1=.TRUE.
       virgin2=.TRUE.
       IF (interface_H_phi%mes/=0) THEN    
          virgin1(interface_H_phi%mesh1) = .FALSE.
          virgin2(interface_H_phi%mesh2) = .FALSE.
       END IF
       IF (interface_H_mu%mes/=0) THEN
          virgin1(interface_H_mu%mesh1) = .FALSE.
          virgin1(interface_H_mu%mesh2) = .FALSE.
       END IF
       !===Create Neumann_bdy_H_sides
       count = 0
       DO ms = 1, H_mesh%mes
          IF (MAXVAL(ABS(H_mesh%rr(1,H_mesh%jjs(:,ms)))).LT.1d-12*H_mesh%global_diameter) CYCLE ! No Neumann BC on the z-axis
          IF (.NOT.virgin1(ms)) CYCLE ! No Neumann BC on H_mu interface
          IF(MINVAL(ABS(H_mesh%sides(ms)-inputs%list_dirichlet_sides_H))==0) CYCLE ! Dirichlet boundary
          !===JLG Jan 22 2018
          IF (inputs%my_periodic%nb_periodic_pairs /=0) THEN
             IF (MINVAL(ABS(inputs%my_periodic%list_periodic-H_mesh%sides(ms))) == 0) CYCLE !Periodic Boundary !JLG Jan 20 2018
          END IF
          !===JLG Jan 22 2018
          count =  count + 1
       END DO
       ALLOCATE(Neumann_bdy_H_sides(count))
       count = 0
       DO ms = 1, H_mesh%mes
          IF (MAXVAL(ABS(H_mesh%rr(1,H_mesh%jjs(:,ms)))).LT.1d-12*H_mesh%global_diameter) CYCLE 
          IF (.NOT.virgin1(ms)) CYCLE
          IF(MINVAL(ABS(H_mesh%sides(ms)-inputs%list_dirichlet_sides_H))==0) CYCLE
          !===JLG Jan 22 2018
          IF (inputs%my_periodic%nb_periodic_pairs /=0) THEN
             IF (MINVAL(ABS(inputs%my_periodic%list_periodic-H_mesh%sides(ms))) == 0) CYCLE !Periodic Boundary !JLG Jan 20 2018
          END IF
          !===JLG Jan 22 2018
          count =  count + 1
          Neumann_bdy_H_sides(count) = ms
       END DO
       !===Create Neumann_bdy_pmag_sides
       count = 0
       DO ms = 1, pmag_mesh%mes
          IF (MAXVAL(ABS(pmag_mesh%rr(1,pmag_mesh%jjs(:,ms)))).LT.1d-12*pmag_mesh%global_diameter) CYCLE ! No Neumann BC on the z-axis
          IF(MINVAL(ABS(pmag_mesh%sides(ms)-inputs%list_dirichlet_sides_H))==0) CYCLE ! Dirichlet boundary
          IF(MINVAL(ABS(pmag_mesh%sides(ms)-inputs%list_inter_H_phi))==0) CYCLE ! No Neumann BC on H-phi interface
          !===JLG Jan 22 2018
          IF (inputs%my_periodic%nb_periodic_pairs /=0) THEN
             IF (MINVAL(ABS(inputs%my_periodic%list_periodic-pmag_mesh%sides(ms))) == 0) CYCLE !Periodic Boundary !JLG Jan 20 2018
          END IF
          !===JLG Jan 22 2018
          count =  count + 1
       END DO
       ALLOCATE(Neumann_bdy_pmag_sides(count))
       count = 0
       DO ms = 1, pmag_mesh%mes
          IF (MAXVAL(ABS(pmag_mesh%rr(1,pmag_mesh%jjs(:,ms)))).LT.1d-12*pmag_mesh%global_diameter) CYCLE 
          IF(MINVAL(ABS(pmag_mesh%sides(ms)-inputs%list_dirichlet_sides_H))==0) CYCLE
          IF(MINVAL(ABS(pmag_mesh%sides(ms)-inputs%list_inter_H_phi))==0) CYCLE
          !===JLG Jan 22 2018
          IF (inputs%my_periodic%nb_periodic_pairs /=0) THEN
             IF (MINVAL(ABS(inputs%my_periodic%list_periodic-pmag_mesh%sides(ms))) == 0) CYCLE !Periodic Boundary !JLG Jan 20 2018
          END IF
          !===JLG Jan 22 2018
          count =  count + 1
          Neumann_bdy_pmag_sides(count) = ms
       END DO
       !===Create Neumann_bdy_phi_sides
       count = 0
       DO ms = 1, phi_mesh%mes
          !IF (PRESENT(index_fourier)) THEN 
          IF (phi_mesh%sides(ms)==index_fourier) CYCLE ! No Neumann BC on Fourier boundary
          !END IF
          IF (.NOT.virgin2(ms)) CYCLE ! No Neumann BC on H_phi interface
          IF (MAXVAL(ABS(phi_mesh%rr(1,phi_mesh%jjs(:,ms)))).LT.1d-12*phi_mesh%global_diameter) CYCLE ! No Neumann BC on the z-axis
          IF (MINVAL(ABS(phi_mesh%sides(ms)-inputs%phi_list_dirichlet_sides))==0) CYCLE ! Dirichlet boundary
          count =  count + 1
       END DO
       ALLOCATE(Neumann_bdy_phi_sides(count))
       count = 0
       DO ms = 1, phi_mesh%mes
          !IF (PRESENT(index_fourier)) THEN
          IF (phi_mesh%sides(ms)==index_fourier) CYCLE
          !END IF
          IF (.NOT.virgin2(ms)) CYCLE
          IF (MAXVAL(ABS(phi_mesh%rr(1,phi_mesh%jjs(:,ms)))).LT.1d-12*phi_mesh%global_diameter) CYCLE 
          IF (MINVAL(ABS(phi_mesh%sides(ms)-inputs%phi_list_dirichlet_sides))==0) CYCLE ! Dirichlet boundary
          count =  count + 1
          Neumann_bdy_phi_sides(count) = ms
       END DO
       !===End Neuman BC for H

       !---------------BOUNDARY CONDITIONS FOR Hxn------------------------------------        
       !===Compute sides that are on Dirichlet boundary (H-H_D)xn=0 
       n = 0
       DO ms = 1, H_mesh%mes
          IF (MINVAL(ABS(H_mesh%sides(ms)-inputs%list_dirichlet_sides_H))/=0) CYCLE
          IF (MAXVAL(ABS(H_mesh%rr(1,H_mesh%jjs(:,ms)))) .LT.1d-12*H_mesh%global_diameter) CYCLE 
          n = n + 1
       END DO
       ALLOCATE(Dirichlet_bdy_H_sides(n))
       n = 0
       DO ms = 1, H_mesh%mes
          IF (MINVAL(ABS(H_mesh%sides(ms)-inputs%list_dirichlet_sides_H))/=0) CYCLE
          IF (MAXVAL(ABS(H_mesh%rr(1,H_mesh%jjs(:,ms)))) .LT.1d-12*H_mesh%global_diameter) CYCLE 
          n = n + 1
          Dirichlet_bdy_H_sides(n) = ms
       END DO
       !===BCs on axis for magnetic field
       CALL vector_without_bc_glob_js_D(H_mesh, list_mode, LA_H, H_mode_global_js_D)
       ALLOCATE(H_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(H_global_D(i)%DRL(SIZE(H_mode_global_js_D(i)%DIL)))
       END DO

       !---------PREPARE phi_js_D ARRAY FOR POTENTIAL---------------------------------
       CALL dirichlet_nodes_parallel(phi_mesh, inputs%phi_list_dirichlet_sides, phi_js_D)
       !SB-CN_LC 2022/01/25
!!$       CALL dirichlet_cavities(comm_one_d(1), interface_H_phi, phi_mesh, phi_js_D)
!!$       ALLOCATE(phi_bv1_D(SIZE(phi_js_D)), phi_bv2_D(SIZE(phi_js_D)))
       !===Account for BCs on axis
       CALL scalar_with_bc_glob_js_D(phi_mesh, list_mode, LA_phi, phi_js_D, phi_mode_global_js_D)
       ALLOCATE(phi_global_D(m_max_c))
       DO i = 1, m_max_c
          ALLOCATE(phi_global_D(i)%DRL(SIZE(phi_mode_global_js_D(i)%DIL)))
          phi_global_D(i)%DRL = 0.d0
       END DO
       !------------------------------------------------------------------------------

       !-------------MATRIX ALLOCATION------------------------------------------------
       ALLOCATE(H_p_phi_mat1(m_max_c),H_p_phi_ksp1(m_max_c))
       ALLOCATE(H_p_phi_mat2(m_max_c),H_p_phi_ksp2(m_max_c))

       IF (SIZE(Dirichlet_bdy_H_sides).GE.1) THEN
          ALLOCATE(sigma_curl_gauss_bdy(H_mesh%gauss%l_Gs*SIZE(Dirichlet_bdy_H_sides),6,SIZE(list_mode)))
          ALLOCATE(J_over_sigma_gauss_bdy(H_mesh%gauss%l_Gs*SIZE(Dirichlet_bdy_H_sides),6,SIZE(list_mode)))
       ELSE
          ALLOCATE(sigma_curl_gauss_bdy(0,0,0))
          ALLOCATE(J_over_sigma_gauss_bdy(0,0,0))
          sigma_curl_gauss_bdy = 0.d0
          J_over_sigma_gauss_bdy = 0.d0
       END IF

       IF (interface_H_mu%mes.GE.1) THEN
          ALLOCATE(sigma_curl_gauss_inter_mu(2*H_mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode)))
          ALLOCATE(J_over_sigma_gauss_inter_mu(2*H_mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode)))
          ALLOCATE(minus_grad_times_der_pot_rhoLi(H_mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode)))
       ELSE
          ALLOCATE(sigma_curl_gauss_inter_mu(0,0,0))
          ALLOCATE(J_over_sigma_gauss_inter_mu(0,0,0))
          ALLOCATE(minus_grad_times_der_pot_rhoLi(0,0,0))
          sigma_curl_gauss_inter_mu = 0.d0
          J_over_sigma_gauss_inter_mu = 0.d0
          minus_grad_times_der_pot_rhoLi=0.d0
       END IF

       IF (SIZE(Neumann_bdy_H_sides).GE.1) THEN
          IF(my_petscworld_rank==0) THEN
             WRITE(*,*) "WARNING, Neumann BC: either sigma or CURL(H)_Neumann is axisymmetric."
          END IF
          ALLOCATE(sigma_tot_gauss_Neumann(H_mesh%gauss%l_Gs*SIZE(Neumann_bdy_H_sides),2,SIZE(list_mode)))
       ELSE
          ALLOCATE(sigma_tot_gauss_Neumann(0,0,0))
          sigma_tot_gauss_Neumann = 0.d0
       END IF
       !------------------------------------------------------------------------------

       DO i = 1, m_max_c !Boucle sur les modes
          mode = list_mode(i)

          tps = user_time()
          CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, H_p_phi_mat1(i), clean=.FALSE.)
          CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, H_p_phi_mat2(i), clean=.FALSE.)
          IF (i == 1) THEN
             CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, tampon1, clean=.FALSE.)
             CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, tampon2, clean=.FALSE.)
             CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, precond1, clean=.FALSE.)
             CALL create_local_petsc_matrix(comm_one_d(1), LA_mhd, precond2, clean=.FALSE.)
          END IF
          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps create_local_petsc_matrix', tps

          tps = user_time()
          CALL mat_H_p_phi_maxwell(H_mesh,pmag_mesh,phi_mesh,interface_H_phi, &
               mode,mu_H_field, mu_phi, one_and_half/dt, stab, R_fourier, index_fourier, &
               LA_H, LA_pmag, LA_phi, H_p_phi_mat1(i), H_p_phi_mat2(i), sigma_nj_m, sigma)
          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps mat_H_p_phi_maxwell', tps

          !Take care of discontinuous mu
          tps = user_time()
          CALL mat_maxwell_mu(H_mesh, jj_v_to_H, interface_H_mu, mode, stab, stab_jump_h,&
               mu_H_field, sigma, LA_H, H_p_phi_mat1(i), H_p_phi_mat2(i), sigma_np)
          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps mat_maxwell_mu', tps

          !Take care of discontinuous rot H
          tps = user_time()
          CALL mat_dirichlet_maxwell(H_mesh, jj_v_to_H, Dirichlet_bdy_H_sides, &
               mode, stab, LA_H, H_p_phi_mat1(i), H_p_phi_mat2(i), sigma_np, sigma)
!!$          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps mat_dirichlet_maxwell', tps

!!$          IF (per) THEN
          IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
             CALL periodic_matrix_petsc(H_phi_per%n_bord, H_phi_per%list, &
                  H_phi_per%perlist, H_p_phi_mat1(i), LA_mhd)
             CALL periodic_matrix_petsc(H_phi_per%n_bord, H_phi_per%list, &
                  H_phi_per%perlist, H_p_phi_mat2(i), LA_mhd)
          END IF
          tps = user_time()
!!$          CALL Dirichlet_M_parallel(H_p_phi_mat1(i),LA_pmag%loc_to_glob(1,pmag_js_D))
!!$          CALL Dirichlet_M_parallel(H_p_phi_mat1(i),LA_phi%loc_to_glob(1,phi_js_D))
          CALL Dirichlet_M_parallel(H_p_phi_mat1(i),pmag_mode_global_js_D(i)%DIL)
          CALL Dirichlet_M_parallel(H_p_phi_mat1(i),phi_mode_global_js_D(i)%DIL)
          CALL Dirichlet_M_parallel(H_p_phi_mat1(i),H_mode_global_js_D(i)%DIL)
!!$          CALL Dirichlet_M_parallel(H_p_phi_mat2(i),LA_pmag%loc_to_glob(1,pmag_js_D))
!!$          CALL Dirichlet_M_parallel(H_p_phi_mat2(i),LA_phi%loc_to_glob(1,phi_js_D))
          CALL Dirichlet_M_parallel(H_p_phi_mat2(i),pmag_mode_global_js_D(i)%DIL)
          CALL Dirichlet_M_parallel(H_p_phi_mat2(i),phi_mode_global_js_D(i)%DIL)
          CALL Dirichlet_M_parallel(H_p_phi_mat2(i),H_mode_global_js_D(i)%DIL)
          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps Dirichlet_M_parallel', tps

          tps = user_time()
          CALL init_solver(inputs%my_par_H_p_phi,H_p_phi_ksp1(i),H_p_phi_mat1(i),comm_one_d(1),&
               solver=inputs%my_par_H_p_phi%solver,precond=inputs%my_par_H_p_phi%precond)
          CALL init_solver(inputs%my_par_H_p_phi,H_p_phi_ksp2(i),H_p_phi_mat2(i),comm_one_d(1),&
               solver=inputs%my_par_H_p_phi%solver,precond=inputs%my_par_H_p_phi%precond)
          tps = user_time() - tps
!!$          WRITE(*,*) ' Tps init_solver', tps

!!$          !==================TEST===================
          CALL MatDestroy(H_p_phi_mat1(i),ierr)
          CALL MatDestroy(H_p_phi_mat2(i),ierr)
       ENDDO
       !------------------------------------------------------------------------------
    ENDIF ! End of once
    tps_tot = user_time()
    tps_cumul = 0
    CALL MPI_COMM_RANK(PETSC_COMM_WORLD, my_petscworld_rank, code)

    IF (inputs%if_coupling_H_x) THEN
       !===Initialize rhoLi_node using input tempn/rho
       rhoLi_node = conc_to_H
       !===Compute derivative of potentiel of rhoLi via FFT
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       bloc_size = SIZE(rhoLi_node,1)/nb_procs+1
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       der_pot_rhoLi_node = rhoLi_node
       CALL FFT_PAR_SCAL_FUNCT(comm_one_d(2), der_pot_rhoLi_node,&
            Derivative_of_potential_from_rhoLi,nb_procs, bloc_size, m_max_pad)
       CALL compute_minus_grad_times_der_pot_rhoLi(comm_one_d(2),H_mesh,interface_H_mu,&
            list_mode,rhoLi_node,der_pot_rhoLi_node,minus_grad_times_der_pot_rhoLi)
    END IF

    !-------------TRANSPORT TERM---------------------------------------------------
    tps = user_time()
    nr_vel = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, vel) 
    H_ext = 2*Hn - Hn1
    B_ext = 2*Bn - Bn1
    IF (nr_vel .LE. 1.d-10) THEN
       NL = 0.d0
    ELSE IF (inputs%type_pb=="fhd" .OR. inputs%type_pb=="mhs") THEN
       NL=0.d0
    ELSE
       CALL MPI_COMM_SIZE(comm_one_d(2), nb_procs, code)
       bloc_size = SIZE(vel,1)/nb_procs+1
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       CALL FFT_PAR_CROSS_PROD_DCL(comm_one_d(2), vel, H_ext, NL, nb_procs, bloc_size, m_max_pad, temps_par)
    ENDIF

    !SB-CN_LC 2022/01/25
!!$    IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
!!$       H_ns = 0.d0
!!$       one_over_sigma_tot = 0.d0
!!$       DO m = 1, H_mesh%me
!!$          DO nj = 1, H_mesh%gauss%n_w
!!$             j = H_mesh%jj(nj,m)
!!$             !Check if node is in Navier-Stokes domain(s)
!!$             IF (jj_v_to_H(j) /= -1) THEN 
!!$                H_ns(j,:,:)      = 2*Hn(j,:,:)- Hn1(j,:,:)
!!$                one_over_sigma_tot(j,:,:) = one_over_sigma_ns_in(jj_v_to_H(j),:,:)/Rem
!!$             ELSE
!!$                DO i = 1, SIZE(list_mode)
!!$                   mode = list_mode(i)
!!$                   IF (mode==0) THEN
!!$                      one_over_sigma_tot(j,1,i) = 1.d0/sigma(m)
!!$                   END IF
!!$                END DO
!!$             END IF
!!$          END DO
!!$       END DO
!!$
!!$       !===Compute (1/sigma_ns_bar - 1/sigma)*CURL(H_ns) in fluid domain and 0 elsewhere
!!$       CALL smb_sigma_prod_curl(comm_one_d(2), H_mesh, jj_v_to_H, list_mode, H_ns, &
!!$            one_over_sigma_tot, sigma_nj_m, sigma, sigma_curl_gauss)
!!$       IF (SIZE(Dirichlet_bdy_H_sides).GE.1) THEN
!!$          CALL smb_sigma_prod_curl_bdy(comm_one_d(2), H_mesh, jj_v_to_H, Dirichlet_bdy_H_sides, list_mode, H_ns, &
!!$               one_over_sigma_tot, sigma_np, sigma, sigma_curl_gauss_bdy)          
!!$       ELSE
!!$          sigma_curl_gauss_bdy = 0.d0
!!$       END IF
!!$       IF (interface_H_mu%mes.GE.1) THEN
!!$          CALL smb_sigma_prod_curl_inter_mu(comm_one_d(2), H_mesh, jj_v_to_H, interface_H_mu, list_mode, H_ns, &
!!$               one_over_sigma_tot, sigma_np, sigma, sigma_curl_gauss_inter_mu)
!!$       ELSE
!!$          sigma_curl_gauss_inter_mu = 0.d0
!!$       END IF
!!$
!!$       !===Compute J/sigma
!!$       CALL smb_current_over_sigma(comm_one_d(2), H_mesh, jj_v_to_H, list_mode, B_ext,&
!!$            mu_H_field, mu_phi, one_over_sigma_tot, time, sigma, J_over_sigma_gauss)
!!$       IF (SIZE(Dirichlet_bdy_H_sides).GE.1) THEN
!!$          CALL smb_current_over_sigma_bdy(comm_one_d(2), H_mesh, jj_v_to_H, Dirichlet_bdy_H_sides,&
!!$               list_mode, B_ext, mu_H_field, mu_phi, one_over_sigma_tot, time, sigma,&
!!$               J_over_sigma_gauss_bdy) 
!!$       ELSE
!!$          J_over_sigma_gauss_bdy = 0.d0
!!$       END IF
!!$       IF (interface_H_mu%mes.GE.1) THEN
!!$          CALL smb_current_over_sigma_inter_mu(comm_one_d(2), H_mesh, jj_v_to_H, interface_H_mu,&
!!$               list_mode, B_ext, mu_H_field, mu_phi, one_over_sigma_tot, time, sigma,&
!!$               J_over_sigma_gauss_inter_mu)  
!!$       ELSE
!!$          J_over_sigma_gauss_inter_mu=0.d0
!!$       END IF
!!$
!!$       !===Compute sigma at the gauss points on Neumann bdy
!!$       IF (SIZE(Neumann_bdy_H_sides).GE.1) THEN
!!$          CALL smb_sigma_Neumann(comm_one_d(2), H_mesh, Neumann_bdy_H_sides,&
!!$               list_mode, one_over_sigma_tot, sigma_tot_gauss_Neumann) 
!!$       ELSE
!!$          sigma_tot_gauss_Neumann = 0.d0
!!$       END IF
!!$    ELSE
!!$       sigma_curl_gauss           = 0.d0
!!$       sigma_curl_gauss_bdy       = 0.d0
!!$       sigma_curl_gauss_inter_mu  = 0.d0
!!$       J_over_sigma_gauss         = 0.d0
!!$       J_over_sigma_gauss_bdy     = 0.d0
!!$       J_over_sigma_gauss_inter_mu= 0.d0
!!$       sigma_tot_gauss_Neumann    = 0.d0
!!$    END IF
    sigma_curl_gauss           = 0.d0
    sigma_curl_gauss_bdy       = 0.d0
    sigma_curl_gauss_inter_mu  = 0.d0
    J_over_sigma_gauss         = 0.d0
    J_over_sigma_gauss_bdy     = 0.d0
    J_over_sigma_gauss_inter_mu= 0.d0
    sigma_tot_gauss_Neumann    = 0.d0
    !SB-CN-LC 2022/01/25

    tps = user_time() - tps; tps_cumul=tps_cumul+tps
    !WRITE(*,*) ' Tps NLS_SFT Maxwell', tps
    !------------------------------------------------------------------------------

    !-------------SOLUTION OF LINEAR SYSTEMS---------------------------------------
    DO i = 1, m_max_c

       ! LC-SB-CN 08/12/2021
       CALL VecZeroEntries(vb_1, ierr)
       CALL VecZeroEntries(vb_2, ierr)
       ! LC-SB-CN 08/12/2021

       mode = list_mode(i)

       !-------------SOURCES TERMS----------------------------------------------------
       tps = user_time()
       !SB-CN-LC 2022/01/25
!!$       DO k = 1, 6
!!$          !rhs_H (:,k) = mu_H_field*(4*Hn(:,k,i)-Hn1(:,k,i))/(2*dt)
!!$          rhs_H (:,k) = (4*Bn(:,k,i)-Bn1(:,k,i))/(2*dt)
!!$       END DO
!!$       DO k = 1, 2
!!$          rhs_phi (:,k) = mu_phi*(4*phin(:,k,i)-phin1(:,k,i))/(2*dt)
!!$       END DO
       !SB-CN-LC 2022/01/25
       DO k = 1, 6
          rhs_H (:,k) = (4*Bn(:,k,i)-Bn1(:,k,i))/(2*dt)
       END DO
       rhs_phi = 0.d0

       !SB-CN-LC 2022/01/25       
!!$       !-------------Integration by parts of the scalar potential------------------
!!$       CALL courant_int_by_parts(H_mesh,phi_mesh,interface_H_phi,sigma,mu_phi,mu_H_field,time,mode, &
!!$            rhs_H, NL(:,:,i), LA_H, LA_phi, vb_1, vb_2, B_ext(:,:,i),& 
!!$            sigma_curl_gauss(:,:,i), J_over_sigma_gauss(:,:,i))
!!$       !! Feb 2010, JLG + FL
!!$       !CALL courant_int_by_parts(H_mesh,phi_mesh,interface_H_phi,sigma,mu_phi,mu_H_field,time,mode, &
!!$       !     rhs_H, NL(:,:,i), LA_H, LA_phi, vb_1, vb_2, H_ext(:,:,i))
       !SB-CN-LC 2022/01/25

       CALL courant_int_by_parts(H_mesh,phi_mesh,interface_H_phi,sigma,mu_phi,mu_H_field,time,mode, &
            rhs_H, NL(:,:,i), LA_H, LA_phi, vb_1, vb_2, B_ext(:,:,i),& 
            sigma_curl_gauss(:,:,i), J_over_sigma_gauss(:,:,i))
       !-------------Integration by parts of the scalar potential------------------

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps courant', tps
       !Takes care of discontinuous mu
       tps = user_time()
       CALL  courant_mu(H_mesh, interface_H_mu, sigma, mu_H_field, time, mode, NL(:,:,i), &
            LA_H, vb_1, vb_2, B_ext(:,:,i), J_over_sigma_gauss_inter_mu(:,:,i), &
            sigma_curl_gauss_inter_mu(:,:,i))

       !SB-CN_LC 2022/01/25
       IF (inputs%if_coupling_H_x) THEN
          CALL jump_rot_H_consistant_rhoLi(H_mesh, jj_v_to_H, interface_H_mu, stab_jump_h, sigma, &
               minus_grad_times_der_pot_rhoLi(:,:,i), LA_H, vb_1, vb_2,sigma_np) 
       ELSE IF (inputs%if_coupling_analytical) THEN
          ff = rot_H_jump_interface(H_mesh,H_mesh%rr,list_mode) 
          CALL jump_rot_H_consistant(H_mesh, jj_v_to_H, interface_H_mu, stab_jump_h, sigma, &
               ff(:,:,i), LA_H, vb_1, vb_2,sigma_np) 
       ELSE
          ff = 0.d0 
          CALL jump_rot_H_consistant(H_mesh, jj_v_to_H, interface_H_mu, stab_jump_h, sigma, &
               ff(:,:,i), LA_H, vb_1, vb_2,sigma_np) 
       END IF
       !SB-CN_LC 2022/01/25

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
!!$       !WRITE(*,*) ' Tps courant_mu', tps

       !JLG, FL, Feb 10, 2011
       !Take care of Dirichlet conditions on H (H x n = Hd x n)
       CALL rhs_dirichlet(H_mesh, Dirichlet_bdy_H_sides, &
            sigma, mu_H_field, time, mode, NL(:,:,i), stab, LA_H, vb_1,vb_2, B_ext(:,:,i), &
            J_over_sigma_gauss_bdy(:,:,i), sigma_curl_gauss_bdy(:,:,i))
       !------------------------------------------------------------------------------

       !-------------INTERFACE INTEGRAL-----------------------------------------------
       tps = user_time()
       CALL surf_int(H_mesh, phi_mesh, pmag_mesh, interface_H_phi, interface_H_mu, inputs%list_dirichlet_sides_H, &
            sigma, mu_phi, mu_H_field, time, mode, LA_H, LA_phi, LA_pmag, vb_1, vb_2, &
            sigma_tot_gauss_Neumann(:,:,i), R_fourier, index_fourier)
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
!!$       WRITE(*,*) ' Tps surf_int', tps
       !------------------------------------------------------------------------------

       !---------------------PERIODIC-------------------
!!$       IF (per) THEN
       IF (inputs%my_periodic%nb_periodic_pairs/=0) THEN
          CALL periodic_rhs_petsc(H_phi_per%n_bord, H_phi_per%list, H_phi_per%perlist, vb_1, LA_mhd)
          CALL periodic_rhs_petsc(H_phi_per%n_bord, H_phi_per%list, H_phi_per%perlist, vb_2, LA_mhd)
       END IF

       !-------------DIRICHLET BOUNDARY CONDITIONS-------------------------------------
       tps = user_time()
!!$       CALL dirichlet_rhs(LA_pmag%loc_to_glob(1,pmag_js_D)-1, pmag_bv_D,vb_1)
!!$       CALL dirichlet_rhs(LA_pmag%loc_to_glob(1,pmag_js_D)-1, pmag_bv_D,vb_2)
       pmag_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(pmag_mode_global_js_D(i)%DIL-1,pmag_global_D(i)%DRL,vb_1)
       CALL dirichlet_rhs(pmag_mode_global_js_D(i)%DIL-1,pmag_global_D(i)%DRL,vb_2)

       IF (SIZE(phi_js_D)>0) THEN
!!$          phi_bv1_D = Phiexact(1,phi_mesh%rr(1:2,phi_js_D), mode, mu_phi, time)
!!$          phi_bv2_D = Phiexact(2,phi_mesh%rr(:,phi_js_D), mode, mu_phi, time)
!!$          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_bv1_D, vb_1)
!!$          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_bv2_D, vb_2)
          !===Recall that axis nodes are at the end of the array
          n = SIZE(phi_js_D)
          phi_global_D(i)%DRL(1:n) = Phiexact(1,phi_mesh%rr(1:2,phi_js_D), mode, mu_phi, time)
          IF (SIZE(phi_global_D(i)%DRL)>n) phi_global_D(i)%DRL(n+1:)=0.d0
          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_global_D(i)%DRL, vb_1)
          phi_global_D(i)%DRL(1:n) = Phiexact(2,phi_mesh%rr(1:2,phi_js_D), mode, mu_phi, time)
          IF (SIZE(phi_global_D(i)%DRL)>n) phi_global_D(i)%DRL(n+1:)=0.d0
          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_global_D(i)%DRL, vb_2)
       ELSE
!!$          phi_bv1_D = 0.d0
!!$          phi_bv2_D = 0.d0
!!$          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_bv1_D, vb_1)
!!$          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_bv2_D, vb_2)
          phi_global_D(i)%DRL=0.d0
          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_global_D(i)%DRL, vb_1)
          CALL dirichlet_rhs(LA_phi%loc_to_glob(1,phi_js_D)-1, phi_global_D(i)%DRL, vb_2)
       END IF

       !===Axis boundary conditions on magnetic field
       H_global_D(i)%DRL = 0.D0
       CALL dirichlet_rhs(H_mode_global_js_D(i)%DIL-1,H_global_D(i)%DRL,vb_1)
       CALL dirichlet_rhs(H_mode_global_js_D(i)%DIL-1,H_global_D(i)%DRL,vb_2)

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
!!$       WRITE(*,*) ' Tps bcs', tps
       !-------------------------------------------------------------------------------

       !-------------SOLVING THE LINEAR SYSTEMS----------------------------------------
       IF (inputs%my_par_H_p_phi%verbose .AND. (i==1)) WRITE(*,*) 'start solving'
       tps = user_time()

       CALL solver(H_p_phi_ksp1(i),vb_1,vx_1,reinit=.FALSE.,verbose=inputs%my_par_H_p_phi%verbose)

       CALL VecGhostUpdateBegin(vx_1,INSERT_VALUES,SCATTER_FORWARD,ierr) 
       CALL VecGhostUpdateEnd(vx_1,INSERT_VALUES,SCATTER_FORWARD,ierr)
       IF (H_mesh%me/=0) THEN
          CALL extract(vx_1_ghost,1,1,LA_mhd,Hn_p1(:,1))
          CALL extract(vx_1_ghost,2,2,LA_mhd,Hn_p1(:,4))
          CALL extract(vx_1_ghost,3,3,LA_mhd,Hn_p1(:,5))
       END IF
       IF (phi_mesh%me/=0) THEN
          CALL extract(vx_1_ghost,5,5,LA_mhd,phin_p1(:,1))
       END IF

       CALL solver(H_p_phi_ksp2(i),vb_2,vx_2,reinit=.FALSE.,verbose=inputs%my_par_H_p_phi%verbose)
       CALL VecGhostUpdateBegin(vx_2,INSERT_VALUES,SCATTER_FORWARD,ierr) 
       CALL VecGhostUpdateEnd(vx_2,INSERT_VALUES,SCATTER_FORWARD,ierr)
       IF (H_mesh%me/=0) THEN
          CALL extract(vx_2_ghost,1,1,LA_mhd,Hn_p1(:,2))
          CALL extract(vx_2_ghost,2,2,LA_mhd,Hn_p1(:,3))
          CALL extract(vx_2_ghost,3,3,LA_mhd,Hn_p1(:,6))
       END IF
       IF (phi_mesh%me/=0) THEN
          CALL extract(vx_2_ghost,5,5,LA_mhd,phin_p1(:,2))
       END IF

       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps solve Maxwell', tps
       !-------------------------------------------------------------------------------


       !-------------UPDATE------------------------------------------------------------
       !JLG AR, Dec 18 2008
       IF (mode==0) THEN
          IF (H_mesh%me /=0) THEN
             Hn_p1 (:,2) = 0.d0
             Hn_p1 (:,4) = 0.d0
             Hn_p1 (:,6) = 0.d0
          END IF
          IF (phi_mesh%me /=0 ) THEN
             phin_p1 (:,2) = 0.d0
          END IF
       END IF
       !JLG AR, Dec 18 2008 

       tps = user_time()
       IF (H_mesh%me /=0) THEN
          Hn1(:,:,i) = Hn(:,:,i)

          Hn(:,1,i)   = Hn_p1(:,1)  
          Hn(:,4,i)   = Hn_p1(:,4)
          Hn(:,5,i)   = Hn_p1(:,5)

          Hn(:,2,i)   = Hn_p1(:,2)
          Hn(:,3,i)   = Hn_p1(:,3)
          Hn(:,6,i)   = Hn_p1(:,6)

          DO k = 1, 6
             Bn1(:,k,i) = Bn(:,k,i) 
             Bn(:,k,i)  = mu_H_field*Hn(:,k,i)
          END DO

       END IF

       IF (phi_mesh%me /= 0) THEN
          phin1(:,:,i) = phin(:,:,i)

          phin(:,1,i) = phin_p1(:,1)

          phin(:,2,i) = phin_p1(:,2)
       END IF
       tps = user_time() - tps; tps_cumul=tps_cumul+tps
       !WRITE(*,*) ' Tps update', tps
       !------------------------------------------------------------------------------ 

    ENDDO

    !===Verbose divergence of velocity
    IF (inputs%verbose_divergence) THEN
       norm = norm_SF(comm_one_d, 'L2',  H_mesh, list_mode, Bn)
       talk_to_me%div_B_L2  = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)/norm
       talk_to_me%time=time
    END IF

    tps_tot = user_time() - tps_tot
!!$    WRITE(*,'(A,2(f13.3,2x),10(I3,x))') ' Tps boucle en temps Maxwell', tps_tot, tps_cumul, list_mode
!!$    WRITE(*,*) ' TIME = ', time, '========================================'

  END SUBROUTINE maxwell_mxs_with_H

  SUBROUTINE  mat_H_p_phi_maxwell(H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       mode, mu_H_field, mu_phi, c_mass, stab, R_fourier, index_fourier, &
       LA_H, LA_pmag, LA_phi, H_p_phi_mat1, H_p_phi_mat2, sigma_nj_m, sigma)
    USE def_type_mesh
    USE Dir_nodes
    USE gauss_points
    USE boundary
    USE input_data ! MODIFICATION: to call sigma_min and mu_min
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE    
    TYPE(mesh_type),              INTENT(IN)    :: H_mesh
    TYPE(mesh_type),              INTENT(IN)    :: pmag_mesh
    TYPE(mesh_type),              INTENT(IN)    :: phi_mesh
    TYPE(interface_type),         INTENT(IN)    :: interface_H_phi 
    INTEGER,                      INTENT(IN)    :: mode   
    REAL(KIND=8),                 INTENT(IN)    :: mu_phi, c_mass
    REAL(KIND=8), DIMENSION(3),   INTENT(IN)    :: stab
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: mu_H_field
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_w,H_mesh%me), INTENT(IN) :: sigma_nj_m
    REAL(KIND=8), DIMENSION(H_mesh%me),INTENT(IN):: sigma
    REAL(KIND=8), OPTIONAL :: R_fourier
    INTEGER,      OPTIONAL :: index_fourier
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_ws,phi_mesh%gauss%l_Gs)   :: w_cs
    REAL(KIND=8), DIMENSION(2,  H_mesh%gauss%n_w, phi_mesh%gauss%l_Gs, H_mesh%mes) :: dw_cs
    INTEGER :: m, l, ms, ls, ni, nj, k, i, j, &
         n_ws1, n_ws2, n_w2, n_w1, m1, m2, ki, kj,ib,jb, ms1, ms2
    REAL(KIND=8) :: x, y, hm1, stab_div, stab_colle_H_phi
    REAL(KIND=8) :: ray, error
!!$    LOGICAL :: mark=.FALSE.
    REAL(KIND=8), DIMENSION(3,H_mesh%gauss%n_w,pmag_mesh%gauss%n_w) :: THpmag   
    REAL(KIND=8), DIMENSION(pmag_mesh%gauss%n_w,pmag_mesh%gauss%n_w) :: Tpmag   
    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w)  :: TH
!!$    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_w,phi_mesh%gauss%n_w):: TPhi

    !MATRICES POUR LES TERMES DE VOLUMES c_mass*mu_H*H + Rot((1/sigma)Rot(H)) - Grad(Div(H))
    !                                    -c_mass*mu_phi*Lap(Phi)
    !========================================================================
    !Le probleme est decouple en deux sous groupes de variables : 
    !H1, H4, H5 et Phi1 d'une part et H2, H3, H6 et Phi2 d'autre part.
    !Les matrices (symetriques sans terme de bord) s'ecrivent : 

    !MATRICE 1 :: 
    ! (------------------------------) 
    ! ( TH1 | TH2 | TH3 |       |    )   H1
    ! (     | TH4 | TH5 |       |    )   H4
    ! (           | TH6 |       |    )   H5
    ! (                 | Tpmag |    )   P1
    ! (                         |TPhi)   Phi1
    ! (------------------------------)

    !MATRICE 2 (TH2 => TH8 et TH5 => TH9:: 
    ! (------------------------) 
    ! ( TH1 | TH8 | TH3 |      )   H2
    ! (     | TH4 | TH9 |      )   H3
    ! (           | TH6 |      )   H6
    ! (                 | TPhi )   Phi2
    ! (------------------------)
    !=========================================================================


    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w)      :: Hsij
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_w,phi_mesh%gauss%n_w)    :: Phisij
    REAL(KIND=8), DIMENSION(6,phi_mesh%gauss%n_w,phi_mesh%gauss%n_w)  :: Sij

    ! MATRICES POUR LES TERMES DE BORDS Hsij et Phisij
    !=================================================
    ! (--------------------------------------------------------------------)
    ! ( Hsij(1)        |        Hsij(2) | Hsij(4)        || Sij(1)         )
    ! (        Hsij(1) | Hsij(3)        |        Hsij(4) ||        Sij(2)  )
    ! (--------------------------------------------------------------------)
    ! (                | Hsij(5)        |                ||        Sij(3)  )
    ! (                |        Hsij(5) |                || Sij(4)         )
    ! (--------------------------------------------------------------------)
    ! ( Hsij(7)        |        Hsij(9) | Hsij(6)        || Sij(5)         )             
    ! (        Hsij(7) | Hsij(8)        |        Hsij(6) ||        Sij(6)  ) 
    ! (====================================================================)
    ! ( Sij'(1)        |        Sij'(3) | Sij'(5)        || Phisij         )
    ! (        Sij'(2) | Sij'(4)        |        Sij'(6) ||        Phisij  )
    ! (------------------------------------------------------------------- )
    !
    ! L'autre partie des termes croises est la symetrique de la premiere
    ! juste apres le calcul du terme de bord dissymetrique    

    !fonctions de forme propres a H_mesh
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: ww_H
    !derivees des fonctions de forme propres a H_mesh
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER :: dw_H
    !jacobien pour H
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: rj_H
    !fonctions de forme propres a phi_mesh
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: ww_phi  
    !derivees des fonctions de forme propres a phi_mesh
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER :: dw_phi 
    !REAL(KIND=8), DIMENSION(2,H_mesh%gauss%n_w,H_mesh%gauss%l_G) :: dwp !JLG Jan 22 2018
    !REAL(KIND=8), DIMENSION(H_mesh%gauss%n_w,H_mesh%gauss%l_G)   :: wwp !JLG Jan 22 2018
    REAL(KIND=8), DIMENSION(2,pmag_mesh%gauss%n_w,H_mesh%gauss%l_G) :: dwp
    REAL(KIND=8), DIMENSION(pmag_mesh%gauss%n_w,H_mesh%gauss%l_G)   :: wwp
    !jacobian for phi
    REAL(KIND=8), DIMENSION(:,:),     POINTER :: rj_phi   

!!$    REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%l_Gs) :: gauss1, gauss2
!!$    INTEGER  :: ls1, ls2
    !REAL(KIND=8) :: ref, diff, mu_H, c_mu_H, c_mu_phi, muhl, &
    REAL(KIND=8) :: mu_H, c_mu_H, c_mu_phi, muhl, &
         dzmuhl, drmuhl, c_div, hloc, viscolm, xij, eps
    !June 8 2008
    REAL(KIND=8) :: c_sym=.0d0 ! Symmetrization of the bilinear form
    !June 8 2008
    !June 2009, JLG, CN, Normalization
    REAL(KIND=8) :: c_lap
    !June 2009, JLG, CN
!!$ FL + CN 22/03/2013
!!$    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: mat_loc1, mat_loc2
!!$    INTEGER     , DIMENSION(:),   ALLOCATABLE   :: idxn, jdxn
    REAL(KIND=8), DIMENSION(3*H_mesh%gauss%n_w+pmag_mesh%gauss%n_w+ &
         phi_mesh%gauss%n_w , 3*H_mesh%gauss%n_w+pmag_mesh%gauss%n_w+ &
         phi_mesh%gauss%n_w)                    :: mat_loc1, mat_loc2
    INTEGER     , DIMENSION(3*H_mesh%gauss%n_w+pmag_mesh%gauss%n_w+ &
         phi_mesh%gauss%n_w)                    :: idxn, jdxn
!!$ FL + CN 22/03/2013
    TYPE(petsc_csr_LA)                          :: LA_H, LA_pmag, LA_phi
    INTEGER                                     :: n_wpmag, n_wH, n_wphi, ix, jx
    !LC 2016/03/25
    REAL(KIND=8) :: sigma_np_gauss
    !LC 2016/03/25
    Mat                                         :: H_p_phi_mat1, H_p_phi_mat2
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (H_p_phi_mat1, ierr)
    CALL MatZeroEntries (H_p_phi_mat2, ierr)
    CALL MatSetOption (H_p_phi_mat1, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
    CALL MatSetOption (H_p_phi_mat2, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    !June 2009, JLG, CN, Normalization
    c_lap = .1d0 
    stab_colle_H_phi = stab(2) 
    stab_div = stab(1)
    !Jan 2010, JLG, CN, Normalization,

    c_mu_phi = c_mass*mu_phi 

    ww_H   => H_mesh%gauss%ww
    dw_H   => H_mesh%gauss%dw
    rj_H   => H_mesh%gauss%rj 
    ww_phi => phi_mesh%gauss%ww
    dw_phi => phi_mesh%gauss%dw
    rj_phi => phi_mesh%gauss%rj

    n_wH = H_mesh%gauss%n_w
    n_wpmag = pmag_mesh%gauss%n_w
    n_wphi = phi_mesh%gauss%n_w

    !==Block on H
!!$    ALLOCATE(mat_loc1(3*n_wH+n_wpmag+n_wphi,3*n_wH+n_wpmag+n_wphi))
!!$    ALLOCATE(mat_loc2(3*n_wH+n_wpmag+n_wphi,3*n_wH+n_wpmag+n_wphi))
!!$    ALLOCATE(jdxn(3*n_wH+n_wpmag+n_wphi),idxn(3*n_wH+n_wpmag+n_wphi))
    DO m = 1, H_mesh%me

       TH = 0.d0

       DO l = 1, H_mesh%gauss%l_G
          !hloc = SQRT(SUM(H_mesh%gauss%rj(:,m)))**(2*alpha)
          hloc = (SQRT(SUM(H_mesh%gauss%rj(:,m)))/H_mesh%global_diameter)**(2*alpha) ! MODIFICATION: normalization for stabilization term for divergence
          !===Compute radius of Gauss point
          !Feb 8 2007, muhl
          muhl = SUM(mu_H_field(H_mesh%jj(:,m))*ww_H(:,l))
          drmuhl = SUM(mu_H_field(H_mesh%jj(:,m))*dw_H(1,:,l,m))
          dzmuhl = SUM(mu_H_field(H_mesh%jj(:,m))*dw_H(2,:,l,m))
          c_mu_H   = c_mass*muhl
          !Feb 8 2007, muhl
          sigma_np_gauss = SUM(sigma_nj_m(:,m)*ww_H(:,l))

          !June 7 2008, Normalization, JLG, FL, May, 28, 2009
          !c_div = stab_div*hloc
          !c_div = stab_div*hloc/muhl
          !c_div = stab_div*hloc/muhl**2
          !c_div    = stab_div
          !c_div    = stab_div/muhl 
          !c_div    = stab_div/muhl**2
          !June 7 2008, Normalization
          c_div = stab_div*hloc/(inputs%mu_min**2*inputs%sigma_min) ! MODIFICATION: normalization for penalization term for divergence

          ray = 0
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             ray = ray + H_mesh%rr(1,i)*ww_H(ni,l)
          END DO

          DO ni = 1, H_mesh%gauss%n_w     
             DO nj = 1, H_mesh%gauss%n_w
                j = H_mesh%jj(nj,m)

                ! TEST
                ! mu_H * <bi,bj> + <Div bi,Div bj> + <(1/sigma) Rot bi,Rot bj>
                TH(1,ni,nj) = TH(1,ni,nj) +  rj_H(l,m) * ray* ( &
                     !DCQ + JLG (Nov 13 2013). Mass integration done same way on LHS and RHS.
                     !                     c_mu_H*ww_H(ni,l)*ww_H(nj,l) &
                     c_mass*mu_H_field(j)*ww_H(nj,l)*ww_H(ni,l) & 
                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) + mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl) &
                     *(muhl*(ww_H(nj,l)/ray+dw_H(1,nj,l,m)) + ww_H(nj,l)*drmuhl))
                !+ stab_div*(ww_H(ni,l)*ww_H(nj,l)/ray**2+dw_H(1,ni,l,m)*dw_H(1,nj,l,m) &
                !+ 1/ray*(ww_H(ni,l)*dw_H(1,nj,l,m)+ww_H(nj,l)*dw_H(1,ni,l,m))))
                !                       

                TH(2,ni,nj) = TH(2,ni,nj)+ rj_H(l,m) * ray* (  &
                     mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma_np_gauss &
                                !DIVERGENCE , June 8 2008
                     + c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) &
                     + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
                !+ stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
                !             

                TH(8,ni,nj) = TH(8,ni,nj)+ rj_H(l,m) * ray* (  &
                     - mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     - c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) &
                     + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
                !-stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
                !           

                TH(3,ni,nj) = TH(3,ni,nj)+ rj_H(l,m) * ray* ( &
                     - dw_H(2,ni,l,m)*dw_H(1,nj,l,m)/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl)*&
                     (muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+ stab_div*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*dw_H(2,nj,l,m))
                !        

                TH(4,ni,nj) = TH(4,ni,nj) + rj_H(l,m) * ray* ( &
                     !                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  & 
                     !DCQ + JLG (Nov 13 2013). Mass integration done same way on LHS and RHS.
                     c_mass*mu_H_field(j)*ww_H(nj,l)*ww_H(ni,l) & 
                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) &
                     + 1/ray**2 *(ww_H(ni,l)+ray*dw_H(1,ni,l,m))*(ww_H(nj,l)&
                     +ray*dw_H(1,nj,l,m)))/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     +c_div*muhl**2*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))
                !+stab_div*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l)) 
                !

                TH(5,ni,nj) = TH(5,ni,nj)  + rj_H(l,m) * ray* (&
                     + mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     +c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))
                !    

                TH(9,ni,nj) = TH(9,ni,nj)  + rj_H(l,m) * ray* (&
                     - mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     - c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !- stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))             
                !

                TH(6,ni,nj) = TH(6,ni,nj) + rj_H(l,m) * ray* ( &
                     !                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  &
                     !DCQ + JLG (Nov 13 2013). Mass integration done same way on LHS and RHS.
                     c_mass*mu_H_field(j)*ww_H(nj,l)*ww_H(ni,l) & 
                     + (mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l) + dw_H(1,ni,l,m)*dw_H(1,nj,l,m))/sigma_np_gauss &
                                !DIVERGENCE, June 8 2008
                     + c_div*(muhl*dw_H(2,ni,l,m) + ww_H(ni,l)*dzmuhl) &
                     *(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
                !+ stab_div*dw_H(2,ni,l,m)*dw_H(2,nj,l,m))
                !
                ! TEST

!!$                ! mu_H * <bi,bj> + <Div bi,Div bj> + <(1/sigma) Rot bi,Rot bj>
!!$                TH(1,ni,nj) = TH(1,ni,nj) +  rj_H(l,m) * ray* ( &
!!$                     c_mu_H*ww_H(ni,l)*ww_H(nj,l) & 
!!$                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) + mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl) &
!!$                     *(muhl*(ww_H(nj,l)/ray+dw_H(1,nj,l,m)) + ww_H(nj,l)*drmuhl))
!!$                !+ stab_div*(ww_H(ni,l)*ww_H(nj,l)/ray**2+dw_H(1,ni,l,m)*dw_H(1,nj,l,m) &
!!$                !+ 1/ray*(ww_H(ni,l)*dw_H(1,nj,l,m)+ww_H(nj,l)*dw_H(1,ni,l,m))))
!!$                !                       
!!$
!!$                TH(2,ni,nj) = TH(2,ni,nj)+ rj_H(l,m) * ray* (  &
!!$                     mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma(m) &
!!$                                !DIVERGENCE , June 8 2008
!!$                     + c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) &
!!$                     + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
!!$                !+ stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
!!$                !             
!!$
!!$                TH(8,ni,nj) = TH(8,ni,nj)+ rj_H(l,m) * ray* (  &
!!$                     - mode/ray**2 * ww_H(ni,l)*(ww_H(nj,l)+ray*dw_H(1,nj,l,m))/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     - c_div*mode/ray*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) &
!!$                     + ww_H(ni,l)*drmuhl)*muhl*ww_H(nj,l))
!!$                !-stab_div*mode/ray*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*ww_H(nj,l))
!!$                !           
!!$
!!$                TH(3,ni,nj) = TH(3,ni,nj)+ rj_H(l,m) * ray* ( &
!!$                     - dw_H(2,ni,l,m)*dw_H(1,nj,l,m)/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     + c_div*(muhl*(ww_H(ni,l)/ray+dw_H(1,ni,l,m)) + ww_H(ni,l)*drmuhl)*&
!!$                     (muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
!!$                !+ stab_div*(ww_H(ni,l)/ray+dw_H(1,ni,l,m))*dw_H(2,nj,l,m))
!!$                !        
!!$
!!$                TH(4,ni,nj) = TH(4,ni,nj) + rj_H(l,m) * ray* ( &
!!$                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  & 
!!$                     + (dw_H(2,ni,l,m)*dw_H(2,nj,l,m) &
!!$                     + 1/ray**2 *(ww_H(ni,l)+ray*dw_H(1,ni,l,m))*(ww_H(nj,l)&
!!$                     +ray*dw_H(1,nj,l,m)))/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     +c_div*muhl**2*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l))
!!$                !+stab_div*mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l)) 
!!$                !
!!$
!!$                TH(5,ni,nj) = TH(5,ni,nj)  + rj_H(l,m) * ray* (&
!!$                     + mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     +c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
!!$                !+stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))
!!$                !    
!!$
!!$                TH(9,ni,nj) = TH(9,ni,nj)  + rj_H(l,m) * ray* (&
!!$                     - mode/ray*dw_H(2,ni,l,m)*ww_H(nj,l)/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     - c_div*mode/ray*muhl*ww_H(ni,l)*(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
!!$                !- stab_div*mode/ray*ww_H(ni,l)*dw_H(2,nj,l,m))             
!!$                !
!!$
!!$                TH(6,ni,nj) = TH(6,ni,nj) + rj_H(l,m) * ray* ( &
!!$                     c_mu_H*ww_H(ni,l)*ww_H(nj,l)  &
!!$                     + (mode**2/ray**2*ww_H(ni,l)*ww_H(nj,l) + dw_H(1,ni,l,m)*dw_H(1,nj,l,m))/sigma(m) &
!!$                                !DIVERGENCE, June 8 2008
!!$                     + c_div*(muhl*dw_H(2,ni,l,m) + ww_H(ni,l)*dzmuhl) &
!!$                     *(muhl*dw_H(2,nj,l,m) + ww_H(nj,l)*dzmuhl))
!!$                !+ stab_div*dw_H(2,ni,l,m)*dw_H(2,nj,l,m))
!!$                !               
             ENDDO
          END DO

       END DO

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3  
          DO ni = 1, n_wH
             i = H_mesh%jj(ni, m)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_wH+ni
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_wH
                   j = H_mesh%jj(nj, m)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_wH+nj
                   jdxn(jx) = jb - 1

                   IF   ((ki == 1) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = TH(1,ni,nj)
                      mat_loc2(ix,jx) = TH(1,ni,nj)
                   ELSEIF   ((ki == 1) .AND. (kj == 2)) THEN
                      mat_loc1(ix,jx) = TH(2,ni,nj)
                      mat_loc2(ix,jx) = TH(8,ni,nj)
                   ELSEIF   ((ki == 2) .AND. (kj == 1)) THEN  
                      mat_loc1(ix,jx) = TH(2,nj,ni)
                      mat_loc2(ix,jx) = TH(8,nj,ni)
                   ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = TH(3,ni,nj)
                      mat_loc2(ix,jx) = TH(3,ni,nj)                          
                   ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = TH(3,nj,ni)  
                      mat_loc2(ix,jx) = TH(3,nj,ni)  
                   ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                      mat_loc1(ix,jx) = TH(4,ni,nj)
                      mat_loc2(ix,jx) = TH(4,ni,nj)
                   ELSEIF   ((ki == 2) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = TH(5,ni,nj)
                      mat_loc2(ix,jx) = TH(9,ni,nj)
                   ELSEIF   ((ki == 3) .AND. (kj == 2)) THEN
                      mat_loc1(ix,jx) = TH(5,nj,ni)   
                      mat_loc2(ix,jx) = TH(9,nj,ni)
                   ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                      mat_loc1(ix,jx) = TH(6,ni,nj)   
                      mat_loc2(ix,jx) = TH(6,ni,nj) 
                   ENDIF

                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_wH, idxn(1:3*n_wH), 3*n_wH, jdxn(1:3*n_wH), &
            mat_loc1(1:3*n_wH,1:3*n_wH), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_wH, idxn(1:3*n_wH), 3*n_wH, jdxn(1:3*n_wH), &
            mat_loc2(1:3*n_wH,1:3*n_wH), ADD_VALUES, ierr)
    END DO

    ! Block on Pmag
    DO m = 1, pmag_mesh%me
       !hloc = stab_div*SQRT(SUM(pmag_mesh%gauss%rj(:,m)))**(2*(1-alpha))
       hloc = stab_div*(SQRT(SUM(pmag_mesh%gauss%rj(:,m)))/pmag_mesh%global_diameter)**(2*(1-alpha)) ! MODIFICATION: normalization for magnetic pressure term
       Tpmag = 0.d0
       DO l = 1, pmag_mesh%gauss%l_G
          !Normalization
          muhl = 1 ! SUM(mu_H_field(H_mesh%jj(1:3,m))*pmag_mesh%gauss%ww(1:3,l))
          !ATTENTION ATTENTION: above line should be replaced by the next one
          !JLG DCQ, July 17, 2013 (this should be the proper normalization)
          !muhl = SUM(mu_H_field(H_mesh%jj(1:3,m))*pmag_mesh%gauss%ww(1:3,l))
          !Normalization
          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, pmag_mesh%gauss%n_w
             i = pmag_mesh%jj(ni,m)
             ray = ray + pmag_mesh%rr(1,i)*pmag_mesh%gauss%ww(ni,l)
          END DO
          !viscolm  = hloc*muhl*pmag_mesh%gauss%rj(l,m)
          viscolm  = (pmag_mesh%global_diameter)**2*inputs%mu_min**2*inputs%sigma_min*hloc*pmag_mesh%gauss%rj(l,m) ! MODIFICATION: normalization for magnetic pressure term
          DO nj = 1, pmag_mesh%gauss%n_w
             j = pmag_mesh%jj(nj, m)
             DO ni = 1, pmag_mesh%gauss%n_w
                i = pmag_mesh%jj(ni, m)
                !grad(u).grad(v) en r et z
                xij = 0.d0
                DO k = 1, 2
                   xij =  xij + pmag_mesh%gauss%dw(k,nj,l,m) * pmag_mesh%gauss%dw(k,ni,l,m)
                END DO
                !blocs diagonaux
                Tpmag(ni,nj) =  Tpmag(ni,nj) + ray * viscolm* xij    &
                     + viscolm*mode**2*pmag_mesh%gauss%ww(ni,l)*pmag_mesh%gauss%ww(nj,l)/ray
             ENDDO
          ENDDO
       ENDDO

       DO ni = 1, pmag_mesh%gauss%n_w 
          i = pmag_mesh%jj(ni, m)
          ib = LA_pmag%loc_to_glob(1,i)
          idxn(ni) = ib - 1 
          DO nj = 1, pmag_mesh%gauss%n_w
             j = pmag_mesh%jj(nj, m)
             jb =  LA_pmag%loc_to_glob(1,j)
             jdxn(nj) = jb - 1 
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_wpmag, idxn(1:n_wpmag), n_wpmag, jdxn(1:n_wpmag), &
            Tpmag(1:n_wpmag,1:n_wpmag), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_wpmag, idxn(1:n_wpmag), n_wpmag, jdxn(1:n_wpmag), &
            Tpmag(1:n_wpmag,1:n_wpmag), ADD_VALUES, ierr)
    ENDDO
    ! End Block on PmagxPmag

    ! Block on PmagxH and HxPmag
    DO m = 1, pmag_mesh%me
       IF (H_mesh%gauss%n_w==3) THEN
          dwp=H_mesh%gauss%dw(:,:,:,m)
          wwp=H_mesh%gauss%ww
       ELSE
          dwp(:,1,:) = H_mesh%gauss%dw(:,1,:,m) + 0.5d0*(H_mesh%gauss%dw(:,5,:,m)+H_mesh%gauss%dw(:,6,:,m))
          dwp(:,2,:) = H_mesh%gauss%dw(:,2,:,m) + 0.5d0*(H_mesh%gauss%dw(:,6,:,m)+H_mesh%gauss%dw(:,4,:,m))
          dwp(:,3,:) = H_mesh%gauss%dw(:,3,:,m) + 0.5d0*(H_mesh%gauss%dw(:,4,:,m)+H_mesh%gauss%dw(:,5,:,m))
          wwp(1,:)   = H_mesh%gauss%ww(1,:) + 0.5d0*(H_mesh%gauss%ww(5,:)+H_mesh%gauss%ww(6,:))
          wwp(2,:)   = H_mesh%gauss%ww(2,:) + 0.5d0*(H_mesh%gauss%ww(6,:)+H_mesh%gauss%ww(4,:))
          wwp(3,:)   = H_mesh%gauss%ww(3,:) + 0.5d0*(H_mesh%gauss%ww(4,:)+H_mesh%gauss%ww(5,:))
       END IF

       THpmag = 0.d0
       DO l = 1, H_mesh%gauss%l_G
          ray = 0.d0
          DO ni = 1, H_mesh%gauss%n_w
             i = H_mesh%jj(ni,m)
             ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
          END DO
          muhl = stab_div*ray*H_mesh%gauss%rj(l,m)*SUM(mu_H_field(H_mesh%jj(:,m))*H_mesh%gauss%ww(:,l))
          DO nj = 1, pmag_mesh%gauss%n_w
             j = pmag_mesh%jj(nj, m)
             DO ni = 1, H_mesh%gauss%n_w
                i = H_mesh%jj(ni, m)
                THpmag(1,ni,nj) = THpmag(1,ni,nj) + muhl*dwp(1,nj,l)*H_mesh%gauss%ww(ni,l)
                THpmag(2,ni,nj) = THpmag(2,ni,nj) - muhl*mode*wwp(nj,l)*H_mesh%gauss%ww(ni,l)/ray
                THpmag(3,ni,nj) = THpmag(3,ni,nj) + muhl*dwp(2,nj,l)*H_mesh%gauss%ww(ni,l)
             END DO
          END DO
       END DO

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       idxn = 0
       jdxn = 0
       DO ni = 1, n_wH
          i = H_mesh%jj(ni, m)
          DO k = 1, 3
             IF (k==2) THEN
                eps=-1
             ELSE
                eps=1
             END IF
             ib = LA_H%loc_to_glob(k,i)
             ix = (k-1)*n_wH + ni
             idxn(ix) = ib - 1
             DO nj = 1, n_wpmag
                j = pmag_mesh%jj(nj, m)
                jb = LA_pmag%loc_to_glob(1,j)
                jx = nj
                jdxn(jx) = jb - 1
                mat_loc1(ix,jx) = THpmag(k,ni,nj)
                mat_loc2(ix,jx) = eps*THpmag(k,ni,nj)
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_wH, idxn(1:3*n_wH), n_wpmag, jdxn(1:n_wpmag), &
            mat_loc1(1:3*n_wH,1:n_wpmag), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_wH, idxn(1:3*n_wH), n_wpmag, jdxn(1:n_wpmag), &
            mat_loc2(1:3*n_wH,1:n_wpmag), ADD_VALUES, ierr)

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ni = 1, n_wpmag
          i = pmag_mesh%jj(ni, m)
          ib =  LA_pmag%loc_to_glob(1,i)
          ix = ni!+3*n_wH
          idxn(ix) = ib - 1
          DO k = 1, 3
             IF (k==2) THEN
                eps=-1
             ELSE
                eps=1
             END IF
             DO nj = 1, n_wH
                j = H_mesh%jj(nj, m)
                jb = LA_H%loc_to_glob(k,j)
                jx = (k-1)*n_wH + nj
                jdxn(jx) = jb - 1
                mat_loc1(ix,jx) = - THpmag(k,nj,ni)
                mat_loc2(ix,jx) = - eps*THpmag(k,nj,ni)
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_wpmag, idxn(1:n_wpmag), 3*n_wH, jdxn(1:3*n_wH), &
            mat_loc1(1:n_wpmag,1:3*n_wH), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_wpmag, idxn(1:n_wpmag), 3*n_wH, jdxn(1:3*n_wH), &
            mat_loc2(1:n_wpmag,1:3*n_wH), ADD_VALUES, ierr)
    END DO
    ! End Block on PmagxH and HxPmag

!!$    !==Block on phi 
!!$    DO m = 1,phi_mesh%me
!!$
!!$       TPhi = 0.d0
!!$
!!$       DO l = 1, phi_mesh%gauss%l_G
!!$
!!$          !===Compute radius of Gauss point
!!$          ray = 0
!!$          DO ni = 1, phi_mesh%gauss%n_w;  i = phi_mesh%jj(ni,m)
!!$             ray = ray + phi_mesh%rr(1,i)*ww_phi(ni,l)
!!$          END DO
!!$
!!$          DO ni = 1, phi_mesh%gauss%n_w
!!$             DO nj = 1, phi_mesh%gauss%n_w
!!$
!!$                !mu_phi * <Grad bi, Grad bj>
!!$                !JLG, FL May 28, 2009
!!$                !On ajoute le laplacien de phi.
!!$                !TPhi(ni,nj) = TPhi(ni,nj) + rj_phi(l,m) * ray* (c_mu_phi) & 
!!$                !     *(dw_phi(1,ni,l,m)*dw_phi(1,nj,l,m)+dw_phi(2,ni,l,m)*dw_phi(2,nj,l,m) &
!!$                !     +mode**2/ray**2*ww_phi(ni,l)*ww_phi(nj,l))
!!$                TPhi(ni,nj) = TPhi(ni,nj) + rj_phi(l,m) * ray* (c_mass+c_lap)*mu_phi & 
!!$                     *(dw_phi(1,ni,l,m)*dw_phi(1,nj,l,m)+dw_phi(2,ni,l,m)*dw_phi(2,nj,l,m) &
!!$                     +mode**2/ray**2*ww_phi(ni,l)*ww_phi(nj,l))
!!$                !JLG, FL May 28, 2009
!!$             ENDDO
!!$          END DO
!!$       END DO
!!$
!!$       !TEST      
!!$       !TPhi = 0.d0
!!$       !TEST
!!$
!!$       DO ni = 1, phi_mesh%gauss%n_w 
!!$          i = phi_mesh%jj(ni, m)
!!$          ib = LA_phi%loc_to_glob(1,i)
!!$          idxn(ni) = ib - 1
!!$          DO nj = 1, phi_mesh%gauss%n_w
!!$             j = phi_mesh%jj(nj, m)
!!$             jb = LA_phi%loc_to_glob(1,j)
!!$             jdxn(nj) = jb - 1
!!$          END DO
!!$       END DO
!!$       CALL MatSetValues(H_p_phi_mat1, n_wphi, idxn(1:n_wphi), n_wphi, jdxn(1:n_wphi), &
!!$            TPhi(1:n_wphi,1:n_wphi), ADD_VALUES, ierr)
!!$       CALL MatSetValues(H_p_phi_mat2, n_wphi, idxn(1:n_wphi), n_wphi, jdxn(1:n_wphi), &
!!$            TPhi(1:n_wphi,1:n_wphi), ADD_VALUES, ierr)
!!$    END DO

!!$    !*********************************************************************************
!!$    !--------------------TERMS on interface_H_phi SIGMA-------------------------------
!!$    !**********************************************************************************
!!$
!!$    !WRITE(*,*) 'Assembling interface_H_phi '
!!$    CALL gauss(phi_mesh)
!!$    n_ws1 = H_mesh%gauss%n_ws
!!$    n_ws2 = phi_mesh%gauss%n_ws
!!$    n_w1  = H_mesh%gauss%n_w
!!$    n_w2  = phi_mesh%gauss%n_w
!!$
!!$    IF (H_mesh%gauss%n_ws == n_ws) THEN
!!$
!!$       DO ms = 1, interface_H_phi%mes
!!$
!!$          ms2 = interface_H_phi%mesh2(ms)
!!$          m2 = phi_mesh%neighs(ms2)
!!$          ms1 = interface_H_phi%mesh1(ms)
!!$          m1 = H_mesh%neighs(ms1)
!!$
!!$          ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
!!$          diff = SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - phi_mesh%rr(:,phi_mesh%jjs(1,ms2)))**2)
!!$          IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
!!$             w_cs = wws
!!$          ELSE                ! 1 = 2
!!$             DO ls = 1, l_Gs
!!$                w_cs(1,ls)= wws(2,ls)
!!$                w_cs(2,ls)= wws(1,ls)
!!$                w_cs(3,ls)= wws(3,ls) 
!!$                WRITE(*,*) ' Ouaps! oder of shape functions changed?'
!!$             END DO
!!$          END IF
!!$
!!$          DO ls = 1, l_Gs
!!$             gauss2(1,ls) = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))*phi_mesh%gauss%wws(:,ls))
!!$             gauss2(2,ls) = SUM(phi_mesh%rr(2,phi_mesh%jjs(:,ms2))*phi_mesh%gauss%wws(:,ls))
!!$             gauss1(1,ls) = SUM(  H_mesh%rr(1,  H_mesh%jjs(:,ms1))*  H_mesh%gauss%wws(:,ls))
!!$             gauss1(2,ls) = SUM(  H_mesh%rr(2,  H_mesh%jjs(:,ms1))*  H_mesh%gauss%wws(:,ls))
!!$          END DO
!!$
!!$          DO ls2 = 1, l_Gs
!!$             ref = SQRT(1.d-8+SUM(gauss2(:,ls2)**2))
!!$             mark = .FALSE.
!!$             DO ls1 = 1, l_Gs
!!$                diff = SQRT(SUM((gauss2(:,ls2)-gauss1(:,ls1))**2))
!!$                IF (diff .LT. 1.d-10) THEN
!!$                   dw_cs(:,:,ls2,ms1) =  H_mesh%gauss%dw_s(:,:,ls1,ms1)
!!$                   mark = .TRUE.
!!$                   EXIT
!!$                END IF
!!$             END DO
!!$             IF (.NOT.mark) WRITE(*,*) ' BUG '
!!$          END DO
!!$
!!$       END DO
!!$
!!$    ELSE
!!$       DO ms = 1, interface_H_phi%mes
!!$
!!$          ms2 = interface_H_phi%mesh2(ms)
!!$          m2 = phi_mesh%neighs(ms2)
!!$          ms1 = interface_H_phi%mesh1(ms)
!!$          m1 = H_mesh%neighs(ms1)
!!$
!!$          ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
!!$          diff = SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - phi_mesh%rr(:,phi_mesh%jjs(1,ms2)))**2)
!!$          IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
!!$             DO ls = 1, l_Gs
!!$                w_cs(1,ls)= wws(1,ls)+0.5*wws(3,ls)
!!$                w_cs(2,ls)= wws(2,ls)+0.5*wws(3,ls)
!!$                w_cs(3,ls)= 0  
!!$             END DO
!!$          ELSE                ! 1 = 2
!!$             DO ls = 1, l_Gs
!!$                w_cs(1,ls)= wws(2,ls)+0.5*wws(3,ls)
!!$                w_cs(2,ls)= wws(1,ls)+0.5*wws(3,ls)
!!$                w_cs(3,ls)= 0 
!!$                WRITE(*,*) ' Ouaps! oder of shape functions changed?'
!!$             END DO
!!$          END IF
!!$
!!$          DO ls = 1, l_Gs
!!$             dw_cs(1,:,ls,ms1) = H_mesh%gauss%dw(1,:,1,m1)
!!$             dw_cs(2,:,ls,ms1) = H_mesh%gauss%dw(2,:,1,m1)
!!$          END DO
!!$
!!$       END DO
!!$    END IF

    error = 0
    DO ms = 1, interface_H_phi%mes

       ms2 = interface_H_phi%mesh2(ms)
       ms1 = interface_H_phi%mesh1(ms)
       m2 = phi_mesh%neighs(ms2)
       m1 =   H_mesh%neighs(ms1)
       mu_H = SUM(mu_H_field(H_mesh%jj(:,m1)))/H_mesh%gauss%n_w
       !JLG, FL, May, 28, 2009
       !hm1 = stab_colle_H_phi/SUM(rjs(:,ms2)) 
       !hm1 = stab_colle_H_phi*(((mu_phi+mu_H)/mu_H)/SUM(rjs(:,ms2)))       
       !JLG, FL, May, 28, 2009
       hm1 = stab_colle_H_phi/(SUM(rjs(:,ms2))*inputs%sigma_min) ! MODIFICATION: normalization for interface H/phi term  

       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC H----------------------------
       !====================================================================================

       !-------------------------------hm1 (bi x ni) . (bj x nj)----------------------------
       !====================================================================================

       Hsij = 0.d0  
       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          DO ni = 1, n_ws1
             DO nj = 1, n_ws1
                y = x * w_cs(ni,ls)*w_cs(nj,ls)
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(rnorms(2,ls,ms2)**2) 
                Hsij(4,ni,nj) = Hsij(4,ni,nj) - y*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)                        
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y                                                
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(rnorms(1,ls,ms2)**2)
             ENDDO
          ENDDO

       ENDDO


       !TEST
       !Hsij = 0.d0
       !Hsij = Hsij / hm1
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3 
          DO ni = 1, n_ws1 
             i = interface_H_phi%jjs1(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1+ni
             idxn(ix) = ib - 1 
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface_H_phi%jjs1(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_ws1+nj
                   jdxn(jx) = jb - 1 
                   IF  ((ki == 1) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(1,ni,nj)
                      mat_loc2(ix,jx) = Hsij(1,ni,nj)
                   ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,ni,nj)
                      mat_loc2(ix,jx) = Hsij(4,ni,nj)
                   ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(4,nj,ni)
                      mat_loc2(ix,jx) = Hsij(4,nj,ni)                                             
                   ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)                                            
                   ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(6,ni,nj)
                      mat_loc2(ix,jx) = Hsij(6,ni,nj)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)

       !====================================================================================
       !------------------------(1/sigma) (Rot bj) . (bi x ni)------------------------------
       !====================================================================================

       Hsij = 0.d0
       DO ls = 1, phi_mesh%gauss%l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme sans derivees
          DO ni = 1,n_ws1
             DO nj = 1, n_ws1
                y = x*w_cs(ni,ls)*w_cs(nj,ls)
                Hsij(2,ni,nj) = Hsij(2,ni,nj) + y * (-mode/ray)*(-rnorms(1,ls,ms2))
                Hsij(3,ni,nj) = Hsij(3,ni,nj) + y *   mode/ray *(-rnorms(1,ls,ms2))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y * (-1/ray)   *(-rnorms(1,ls,ms2))
                Hsij(8,ni,nj) = Hsij(8,ni,nj) + y * (-mode/ray)*(-rnorms(2,ls,ms2))
                Hsij(9,ni,nj) = Hsij(9,ni,nj) + y *   mode/ray *(-rnorms(2,ls,ms2))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1 + ni
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface_H_phi%jjs1(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_ws1 + nj
                   jdxn(jx) = jb - 1
                   IF  ( (ki == 2) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(2,ni,nj)
                      mat_loc2(ix,jx) = Hsij(3,ni,nj)
                   ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)
                   ELSEIF  ( (ki == 2) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(8,ni,nj)
                      mat_loc2(ix,jx) = Hsij(9,ni,nj) 
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)

       !Feb 2 2007
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       Hsij=c_sym*Hsij !SYM
       DO ki= 1, 3
          DO ni = 1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1 + ni
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = interface_H_phi%jjs1(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_ws1 + nj
                   jdxn(jx) = jb - 1
                   IF  ( (kj == 2) .AND. (ki == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(2,nj,ni)
                      mat_loc2(ix,jx) = Hsij(3,nj,ni)
                   ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                      mat_loc1(ix,jx) = Hsij(5,nj,ni)
                      mat_loc2(ix,jx) = Hsij(5,nj,ni)
                   ELSEIF  ( (kj == 2) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(8,nj,ni)
                      mat_loc2(ix,jx) = Hsij(9,nj,ni)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       !feb 2 2007
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)

       Hsij = 0.d0
       DO ls = 1, phi_mesh%gauss%l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray /sigma(m1)

          !termes avec derivees
          DO ni = 1,n_ws1
             y = x*w_cs(ni,ls)
             DO nj = 1, n_w1
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(-dw_cs(2,nj,ls,ms1))*(-rnorms(2,ls,ms2))
                Hsij(4,ni,nj) = Hsij(4,ni,nj) + y*  dw_cs(1,nj,ls,ms1) *(-rnorms(2,ls,ms2))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + &
                     y*(-dw_cs(2,nj,ls,ms1)*(-rnorms(2,ls,ms2))-dw_cs(1,nj,ls,ms1)*(-rnorms(1,ls,ms2)))
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(-dw_cs(1,nj,ls,ms1))*(-rnorms(1,ls,ms2))
                Hsij(7,ni,nj) = Hsij(7,ni,nj) + y*  dw_cs(2,nj,ls,ms1) *(-rnorms(1,ls,ms2))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             ib =  LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1 + ni
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_w1
                   j = H_mesh%jj(nj,m1)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w1 + nj
                   jdxn(jx) = jb - 1
                   IF  ((ki == 1) .AND. (kj == 1))  THEN
                      mat_loc1(ix,jx) = Hsij(1,ni,nj)
                      mat_loc2(ix,jx) = Hsij(1,ni,nj)
                   ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,ni,nj)
                      mat_loc2(ix,jx) = Hsij(4,ni,nj)
                   ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)
                   ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                      mat_loc1(ix,jx) = Hsij(6,ni,nj)   
                      mat_loc2(ix,jx) = Hsij(6,ni,nj) 
                   ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(7,ni,nj)
                      mat_loc2(ix,jx) = Hsij(7,ni,nj)
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc1(1:3*n_ws1,1:3*n_w1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc2(1:3*n_ws1,1:3*n_w1), ADD_VALUES, ierr)

       !Feb 2 2007
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       Hsij=c_sym*Hsij !SYM
       DO ki = 1, 3
          DO ni = 1, n_w1
             i = H_mesh%jj(ni,m1)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_w1 + ni
             idxn(ix) = ib - 1
             DO kj= 1, 3
                DO nj = 1, n_ws1
                   j = interface_H_phi%jjs1(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = (kj-1)*n_ws1 + nj
                   jdxn(jx) = jb - 1
                   IF  ((kj == 1) .AND. (ki == 1))  THEN
                      mat_loc1(ix,jx) = Hsij(1,nj,ni)
                      mat_loc2(ix,jx) = Hsij(1,nj,ni)
                   ELSEIF  ((kj == 1) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,nj,ni)
                      mat_loc2(ix,jx) = Hsij(4,nj,ni)
                   ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                      mat_loc1(ix,jx) = Hsij(5,nj,ni)
                      mat_loc2(ix,jx) = Hsij(5,nj,ni)
                   ELSEIF  ((kj == 3) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(6,nj,ni)
                      mat_loc2(ix,jx) = Hsij(6,nj,ni)
                   ELSEIF  ((kj == 3) .AND. (ki == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(7,nj,ni)
                      mat_loc2(ix,jx) = Hsij(7,nj,ni)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_w1, idxn(1:3*n_w1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_w1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_w1, idxn(1:3*n_w1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_w1,1:3*n_ws1), ADD_VALUES, ierr)
       !Feb 2 2007


       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC PHI--------------------------
       !====================================================================================

       !------------------------hm1 (Grad(phi_i) x ni).(Grad(phi_j) x nj)-------------------
       !====================================================================================

       Phisij = 0.d0

       DO ls = 1, phi_mesh%gauss%l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme sans derivee
          DO ni=1, n_ws2
             DO nj=1, n_ws2
                Phisij(ni,nj) = Phisij(ni,nj) + x*mode**2/ray**2*wws(ni,ls)*wws(nj,ls) 
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Phisij = 0.d0
       !Phisij = Phisij/hm1
       !TEST
       DO ni = 1, n_ws2
          i =  interface_H_phi%jjs2(ni,ms)
          ib = LA_phi%loc_to_glob(1,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_ws2
             j = interface_H_phi%jjs2(nj,ms)
             jb = LA_phi%loc_to_glob(1,j)
             jdxn(nj) = jb - 1
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), n_ws2, jdxn(1:n_ws2), &
            Phisij(1:n_ws2,1:n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), n_ws2, jdxn(1:n_ws2), &
            Phisij(1:n_ws2,1:n_ws2), ADD_VALUES, ierr)

       Phisij = 0.d0
       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme avec derivee
          DO ni = 1, n_w2
             DO nj = 1, n_w2
                Phisij(ni,nj) = Phisij(ni,nj) + x*( &
                     (dw_s(2,ni,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,ni,ls,ms2)*rnorms(2,ls,ms2))* &
                     (dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,nj,ls,ms2)*rnorms(2,ls,ms2)))
             ENDDO
          ENDDO
       ENDDO

       !Phisij = 0.d0
       !Phisij = Phisij/hm1
       !TEST

       DO ni = 1, n_w2
          i = phi_mesh%jj(ni, m2)
          ib = LA_phi%loc_to_glob(1,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_w2
             j = phi_mesh%jj(nj, m2)
             jb = LA_phi%loc_to_glob(1,j)
             jdxn(nj) = jb - 1
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_w2, idxn(1:n_w2), n_w2, jdxn(1:n_w2), &
            Phisij(1:n_w2,1:n_w2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_w2, idxn(1:n_w2), n_w2, jdxn(1:n_w2), &
            Phisij(1:n_w2,1:n_w2), ADD_VALUES, ierr)
       !====================================================================================   
       !------------------------------------TERMES CROISES----------------------------------
       !====================================================================================

       !====================================================================================
       !------------------------hm1 (bi x ni) . (Grad(phi_j) x nj)--------------------------
       !------------------      + hm1(Grad(phi_i) x ni).(bj x nj)---------------------------
       !====================================================================================

       Sij = 0.d0
       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme sans derivee
          DO ni = 1, n_ws1
             DO nj = 1, n_ws2
                Sij(3,ni,nj) = Sij(3,ni,nj) + x*(mode/ray)*w_cs(ni,ls)*wws(nj,ls)
             ENDDO
          ENDDO
       ENDDO
       Sij(4,:,:) = -Sij(3,:,:)

       !TEST
       !Sij = 0.d0
       !Sij = Sij /hm1
       !TEST

       ki = 2
       DO ni = 1, n_ws1 
          i = interface_H_phi%jjs1(ni,ms)
          ib = LA_H%loc_to_glob(ki,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_ws2
             j = interface_H_phi%jjs2(nj,ms)
             jb = LA_phi%loc_to_glob(1,j)
             jdxn(nj) = jb - 1
          END DO
       ENDDO
       CALL MatSetValues(H_p_phi_mat1, n_ws1, idxn(1:n_ws1), n_ws2, jdxn(1:n_ws2), &
            Sij(3,1:n_ws1,1:n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws1, idxn(1:n_ws1), n_ws2, jdxn(1:n_ws2), &
            Sij(4,1:n_ws1,1:n_ws2), ADD_VALUES, ierr)

       !TEST SYM
       !Feb 2 2003
       !Sij = 0.d0
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       kj = 2
       DO ni = 1, n_ws2
          i = interface_H_phi%jjs2(ni,ms)
          ib = LA_phi%loc_to_glob(1,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_ws1
             j = interface_H_phi%jjs1(nj,ms)
             jb = LA_H%loc_to_glob(kj,j)
             jdxn(nj) = jb - 1
             mat_loc1(ni,nj) = Sij(3,nj,ni)
             mat_loc2(ni,nj) = Sij(4,nj,ni)
          END DO
       ENDDO
       CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), n_ws1, jdxn(1:n_ws1), &
            mat_loc1(1:n_ws2,1:n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), n_ws1, jdxn(1:n_ws1), &
            mat_loc2(1:n_ws2,1:n_ws1), ADD_VALUES, ierr)

       !Feb 2 2003
       !TEST SYM
       Sij = 0.d0
       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray 

          !terme avec derivee
          DO ni = 1, n_ws1
             y = x * w_cs(ni,ls)
             DO nj = 1, n_w2
                Sij(1,ni,nj) = Sij(1,ni,nj) + &
                     y*(-dw_s(1,nj,ls,ms2)*rnorms(2,ls,ms2)**2 + dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2)*rnorms(2,ls,ms2))
                Sij(5,ni,nj) = Sij(5,ni,nj) + & 
                     y*(-dw_s(2,nj,ls,ms2)*rnorms(1,ls,ms2)**2 + dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)) 
             ENDDO
          ENDDO
       ENDDO

       !TEST
       !Sij = 0.d0
       !Sij = Sij /hm1
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3 
          DO ni = 1, n_ws1 
             i = interface_H_phi%jjs1(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1 + ni
             idxn(ix) = ib - 1
             DO nj = 1, n_w2
                j = phi_mesh%jj(nj,m2)
                jb = LA_phi%loc_to_glob(1,j)
                jx = nj
                jdxn(jx) = jb - 1
                IF (ki == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,ni,nj)
                   mat_loc2(ix,jx) = Sij(1,ni,nj)
                ELSEIF  (ki == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,ni,nj)
                   mat_loc2(ix,jx) = Sij(5,ni,nj)
                END IF
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), n_w2, jdxn(1:n_w2), &
            mat_loc1(1:3*n_ws1,1:n_w2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), n_w2, jdxn(1:n_w2), &
            mat_loc2(1:3*n_ws1,1:n_w2), ADD_VALUES, ierr)

       !TEST SYM
       !Feb 2 2003
       !Sij = 0.d0
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ni = 1, n_w2
          i = phi_mesh%jj(ni,m2)
          ib = LA_phi%loc_to_glob(1,i)
          ix = ni 
          idxn(ix) = ib - 1
          DO kj=1,3
             DO nj = 1, n_ws1
                j = interface_H_phi%jjs1(nj,ms)
                jb = LA_H%loc_to_glob(kj,j)
                jx = (kj-1)*n_ws1 + nj
                jdxn(jx) = jb - 1
                IF (kj == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,nj,ni)
                   mat_loc2(ix,jx) = Sij(1,nj,ni)
                ELSEIF  (kj == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,nj,ni)
                   mat_loc2(ix,jx) = Sij(5,nj,ni)
                ENDIF
             END DO
          END DO
       ENDDO
       CALL MatSetValues(H_p_phi_mat1, n_w2, idxn(1:n_w2), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:n_w2,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_w2, idxn(1:n_w2), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:n_w2,1:3*n_ws1), ADD_VALUES, ierr)

       !TEST SYM
       !Feb 2 2003

       !====================================================================================
       !----------------------(1/sigma) (Rot bj).(Grad(phi_i) x ni)-------------------------
       !====================================================================================
       !        GOTO 200


       Sij = 0.d0
       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme sans derivee
          DO ni = 1, n_ws2
             DO nj = 1, n_ws1
                y = x * wws(ni,ls)*w_cs(nj,ls)
                Sij(1,ni,nj) = Sij(1,ni,nj) + y*( mode/ray)**2*rnorms(1,ls,ms2)
                Sij(3,ni,nj) = Sij(3,ni,nj) + y*( mode/ray**2)*rnorms(1,ls,ms2)
                Sij(4,ni,nj) = Sij(4,ni,nj) + y*(-mode/ray**2)*rnorms(1,ls,ms2)
                Sij(5,ni,nj) = Sij(5,ni,nj) + y*( mode/ray)**2*rnorms(2,ls,ms2)
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Sij = 0.d0
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ni = 1, n_ws2 
          i = interface_H_phi%jjs2(ni,ms)
          ib = LA_phi%loc_to_glob(1,i)
          ix = ni 
          idxn(ix) = ib - 1
          DO kj =1,3
             DO nj = 1, n_ws1
                j = interface_H_phi%jjs1(nj,ms)
                jb = LA_H%loc_to_glob(kj,j)
                jx = (kj-1)*n_ws1 + nj
                jdxn(jx) = jb - 1
                IF (kj == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,ni,nj)
                   mat_loc2(ix,jx) = Sij(1,ni,nj)                          
                ELSEIF (kj == 2)  THEN
                   mat_loc1(ix,jx) = Sij(3,ni,nj)
                   mat_loc2(ix,jx) = Sij(4,ni,nj)  
                ELSEIF  (kj == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,ni,nj)
                   mat_loc2(ix,jx) = Sij(5,ni,nj)
                ENDIF
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:n_ws2,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:n_ws2,1:3*n_ws1), ADD_VALUES, ierr)

       !Feb 2 2007
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       Sij = c_sym*Sij !SYM
       DO ki =1,3
          DO ni = 1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_ws1 + ni
             idxn(ix) = ib - 1
             DO nj = 1, n_ws2
                j = interface_H_phi%jjs2(nj,ms)
                jb = LA_phi%loc_to_glob(1,j)
                jx = nj 
                jdxn(jx) = jb - 1
                IF (ki == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,nj,ni)
                   mat_loc2(ix,jx) = Sij(1,nj,ni)
                ELSEIF (ki == 2)  THEN
                   mat_loc1(ix,jx) = Sij(3,nj,ni)
                   mat_loc2(ix,jx) = Sij(4,nj,ni)
                ELSEIF  (ki == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,nj,ni)
                   mat_loc2(ix,jx) = Sij(5,nj,ni)
                ENDIF
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1, idxn(1:3*n_ws1), n_ws2, jdxn(1:n_ws2), &
            mat_loc1(1:3*n_ws1,1:n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1, idxn(1:3*n_ws1), n_ws2, jdxn(1:n_ws2), &
            mat_loc2(1:3*n_ws1,1:n_ws2), ADD_VALUES, ierr)
       !Feb 2 2007

       Sij = 0.d0

       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme avec derivee de bi seulement
          DO ni = 1, n_ws2
             y =  x*wws(ni,ls)*mode/ray
             DO nj = 1, n_w1 
                Sij(3,ni,nj) = Sij(3,ni,nj) + &
                     y*(dw_cs(2,nj,ls,ms1)*rnorms(2,ls,ms2) + dw_cs(1,nj,ls,ms1)*rnorms(1,ls,ms2))
             ENDDO
          ENDDO
       ENDDO
       Sij(4,:,:) = -Sij(3,:,:)
       !TEST
       !Sij = 0.d0
       !TEST
       kj=2
       DO ni = 1, n_ws2 
          i = interface_H_phi%jjs2(ni,ms)
          ib = LA_phi%loc_to_glob(1,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_w1
             j = H_mesh%jj(nj,m1)
             jb = LA_H%loc_to_glob(kj,j)
             jdxn(nj) = jb - 1
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), n_w1, jdxn(1:n_w1), &
            Sij(3,1:n_ws2,1:n_w1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), n_w1, jdxn(1:n_w1), &
            Sij(4,1:n_ws2,1:n_w1), ADD_VALUES, ierr)

       !Feb 2 2007
       Sij = c_sym*Sij !SYM

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       ki=2
       DO ni = 1, n_w1
          i = H_mesh%jj(ni,m1)
          ib =  LA_H%loc_to_glob(ki,i)
          idxn(ni) = ib - 1
          DO nj = 1, n_ws2
             j = interface_H_phi%jjs2(nj,ms)
             jb = LA_phi%loc_to_glob(1,j)
             jdxn(nj) = jb - 1
             mat_loc1(ix,jx) = Sij(3,nj,ni)
             mat_loc2(ix,jx) = Sij(4,nj,ni)
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_w1, idxn(1:n_w1), n_ws2, jdxn(1:n_ws2), &
            mat_loc1(1:n_w1,1:n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_w1, idxn(1:n_w1), n_ws2, jdxn(1:n_ws2), &
            mat_loc2(1:n_w1,1:n_ws2), ADD_VALUES, ierr)
       !Feb 2 2007

       Sij = 0.d0
       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray/sigma(m1)

          !terme avec derivee de phi et derivee de bi
          DO ni = 1, n_w2
             y =  x*(dw_s(2,ni,ls,ms2)*rnorms(1,ls,ms2) - dw_s(1,ni,ls,ms2)*rnorms(2,ls,ms2))
             DO nj = 1, n_w1
                Sij(1,ni,nj) = Sij(1,ni,nj) +   y *dw_cs(2,nj,ls,ms1) 
                Sij(5,ni,nj) = Sij(5,ni,nj) + (-y)*dw_cs(1,nj,ls,ms1)
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Sij = 0.d0
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ni = 1, n_w2 
          i = phi_mesh%jj(ni,m2)
          ib =  LA_phi%loc_to_glob(1,i)
          ix = ni
          idxn(ix) = ib - 1
          DO nj = 1, n_w1
             j = H_mesh%jj(nj,m1)
             DO kj=1,3
                jb = LA_H%loc_to_glob(kj,j)
                jx = (kj-1)*n_w1 + nj
                jdxn(jx) = jb - 1           
                IF (kj == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,ni,nj)
                   mat_loc2(ix,jx) = Sij(1,ni,nj)
                ELSEIF  (kj == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,ni,nj)
                   mat_loc2(ix,jx) = Sij(5,ni,nj)
                ENDIF
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_w2, idxn(1:n_w2), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc1(1:n_w2,1:3*n_w1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_w2, idxn(1:n_w2), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc2(1:n_w2,1:3*n_w1), ADD_VALUES, ierr)

       !Feb 2 2007
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       Sij=c_sym*Sij !SYM
       DO ki=1,3
          DO ni = 1, n_w1
             i = H_mesh%jj(ni,m1)
             ib = LA_H%loc_to_glob(ki,i)
             ix = (ki-1)*n_w1 + ni
             idxn(ix) = ib - 1 
             DO nj = 1, n_w2
                j = phi_mesh%jj(nj,m2)
                jb =  LA_phi%loc_to_glob(1,j)
                jx = nj
                jdxn(jx) = jb - 1
                IF (ki == 1)  THEN
                   mat_loc1(ix,jx) = Sij(1,nj,ni)
                   mat_loc2(ix,jx) = Sij(1,nj,ni)
                ELSEIF  (ki == 3)  THEN
                   mat_loc1(ix,jx) = Sij(5,nj,ni)
                   mat_loc2(ix,jx) = Sij(5,nj,ni)
                ENDIF
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_w1, idxn(1:3*n_w1), n_w2, jdxn(1:n_w2), &
            mat_loc1(1:3*n_w1,1:n_w2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_w1, idxn(1:3*n_w1), n_w2, jdxn(1:n_w2), &
            mat_loc2(1:3*n_w1,1:n_w2), ADD_VALUES, ierr)

       !JLG, FL, May, 28, 2009
       !Ajout du laplacien de phi
       Sij = 0.d0
       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* phi_mesh%gauss%wws(:,ls))
          muhl = SUM(mu_H_field(interface_H_phi%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
          x = c_lap*muhl*rjs(ls,ms2)*ray
          DO ni = 1, n_ws2 
             DO nj = 1, n_ws1
                Sij(1,ni,nj) = Sij(1,ni,nj) - x*w_cs(nj,ls)*wws(ni,ls)*rnorms(1,ls,ms2)
                Sij(5,ni,nj) = Sij(5,ni,nj) - x*w_cs(nj,ls)*wws(ni,ls)*rnorms(2,ls,ms2)
             ENDDO
          END DO
       END DO

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ni = 1, n_ws2
          i = interface_H_phi%jjs2(ni,ms)
          ib =  LA_phi%loc_to_glob(1,i)
          ix = ni
          idxn(ix) = ib - 1
          DO nj = 1, n_ws1
             j = interface_H_phi%jjs1(nj,ms)
             jb = LA_H%loc_to_glob(1,j)
             jx = nj         !(1-1)*n_ws1 + nj
             jdxn(jx) = jb - 1 
             mat_loc1(ix,jx) = Sij(1,ni,nj)
             mat_loc2(ix,jx) = Sij(1,ni,nj)

             jb = LA_H%loc_to_glob(3,j)
             jx = n_ws1 + nj !(3-1)*n_ws1 + nj
             jdxn(jx) = jb - 1 
             mat_loc1(ix,jx) = Sij(5,ni,nj)
             mat_loc2(ix,jx) = Sij(5,ni,nj)
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), 2*n_ws1, jdxn(1:2*n_ws1), &
            mat_loc1(1:n_ws2,1:2*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), 2*n_ws1, jdxn(1:2*n_ws1), &
            mat_loc2(1:n_ws2,1:2*n_ws1), ADD_VALUES, ierr)
       !JLG, FL, May, 28, 2009

       !Feb 2 2007
       !==================

       !(use .true. for convergence tests)
       !June 6 2008, I put back (.true.) always.
       !Works much better when mu is discontinuous.
       !Mars 22 2007

       IF (stab(2) > 1.d-12) THEN
          !IF (.FALSE.) THEN
          !Mars 22 2007 
          !Enforcing weak continuity on the normal components
          Hsij   = 0.d0 
          Sij    = 0.d0
          Phisij = 0.d0


          ms2 = interface_H_phi%mesh2(ms)
          !hm1 = SUM(rjs(:,ms2))**(2*alpha-1) 
          hm1 =(SUM(rjs(:,ms2))/H_mesh%global_diameter)**(2*alpha-1)/(inputs%sigma_min*inputs%mu_min**2*H_mesh%global_diameter) ! MODIFICATION: normalization for divergence stabilization term

          DO ls = 1, l_Gs

             !Feb 8 2007, muhl
             muhl = SUM(mu_H_field(interface_H_phi%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
             !Feb 8 2007, muhl
             ray = 0.d0
             DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,ms2)
                ray = ray + phi_mesh%rr(1,i)* phi_mesh%gauss%wws(ni,ls)
             END DO


             !ray = ray*hm1*rjs(ls,ms2)
             !June 8, 2008, Normalization, JLG, FL, May, 28, 2009
             ray = stab_div*ray*hm1*rjs(ls,ms2)
             !ray = stab_div*ray*hm1*rjs(ls,ms2)/muhl 
             !ray = stab_div*ray*hm1*rjs(ls,ms2)/muhl**2
             !June 8, 2008, Normalization, JLG, FL, May, 28, 2009
             DO ni = 1, n_ws1
                DO nj = 1, n_ws1
                   x = muhl**2*w_cs(ni,ls)*w_cs(nj,ls)*ray
                   Hsij(1,ni,nj) = Hsij(1,ni,nj) + x*rnorms(1,ls,ms2)**2
                   Hsij(4,ni,nj) = Hsij(4,ni,nj) + x*rnorms(1,ls,ms2)*rnorms(2,ls,ms2)
                   Hsij(6,ni,nj) = Hsij(6,ni,nj) + x*rnorms(2,ls,ms2)**2
                END DO

                DO nj = 1, n_w2
                   x = muhl*mu_phi*w_cs(ni,ls)*(dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2) +&
                        dw_s(2,nj,ls,ms2)*rnorms(2,ls,ms2))*ray
                   Sij(1,ni,nj) = Sij(1,ni,nj) - x*rnorms(1,ls,ms2)
                   Sij(5,ni,nj) = Sij(5,ni,nj) - x*rnorms(2,ls,ms2)
                ENDDO
             ENDDO

             DO ni = 1, n_w2
                DO nj = 1, n_w2
                   x = mu_phi**2*(dw_s(1,ni,ls,ms2)*rnorms(1,ls,ms2) + dw_s(2,ni,ls,ms2)*rnorms(2,ls,ms2))* &
                        (dw_s(1,nj,ls,ms2)*rnorms(1,ls,ms2) + dw_s(2,nj,ls,ms2)*rnorms(2,ls,ms2))*ray
                   Phisij(ni,nj) = Phisij(ni,nj) + x
                ENDDO
             ENDDO

          END DO
          Sij(2,:,:) = Sij(1,:,:)
          Sij(6,:,:) = Sij(5,:,:)


          mat_loc1 = 0.d0
          mat_loc2 = 0.d0
          DO ni = 1, n_ws1
             i = H_mesh%jjs(ni,ms1)
             DO ki= 1, 3, 2
                ib = LA_H%loc_to_glob(ki,i) 
                ix = (ki/2)*n_ws1 + ni
                idxn(ix) = ib - 1 
                DO nj = 1, n_ws1
                   j = H_mesh%jjs(nj,ms1)
                   DO kj = 1, 3, 2
                      jb = LA_H%loc_to_glob(kj,j) 
                      jx = (kj/2)*n_ws1 + nj
                      jdxn(jx) = jb - 1 
                      IF (ki*kj==1) THEN
                         mat_loc1(ix,jx) = Hsij(1,ni,nj)
                         mat_loc2(ix,jx) = Hsij(1,ni,nj)
                      ELSE IF (ki*kj==9) THEN
                         mat_loc1(ix,jx) = Hsij(6,ni,nj)
                         mat_loc2(ix,jx) = Hsij(6,ni,nj)
                      ELSE IF (ki*kj==3) THEN
                         mat_loc1(ix,jx) = Hsij(4,ni,nj)
                         mat_loc2(ix,jx) = Hsij(4,ni,nj)
                      END IF
                   END DO
                END DO

                DO nj = 1, n_w2
                   j = phi_mesh%jj(nj,m2)
                   jb = LA_phi%loc_to_glob(1,j)
                   jx = 2*n_ws1 + nj
                   jdxn(jx) = jb - 1 
                   mat_loc1(ix,jx) = Sij(2*ki-1,ni,nj)
                   mat_loc2(ix,jx) = Sij(2*ki-1,ni,nj)
                END DO
             ENDDO
          ENDDO
          CALL MatSetValues(H_p_phi_mat1, 2*n_ws1, idxn(1:2*n_ws1), 2*n_ws1+n_w2, jdxn(1:2*n_ws1+n_w2), &
               mat_loc1(1:2*n_ws1,1:2*n_ws1+n_w2), ADD_VALUES, ierr)
          CALL MatSetValues(H_p_phi_mat2, 2*n_ws1, idxn(1:2*n_ws1), 2*n_ws1+n_w2, jdxn(1:2*n_ws1+n_w2), &
               mat_loc2(1:2*n_ws1,1:2*n_ws1+n_w2), ADD_VALUES, ierr)

          mat_loc1 = 0.d0
          mat_loc2 = 0.d0
          DO ni = 1, n_w2
             i = phi_mesh%jj(ni,m2)
             ib = LA_phi%loc_to_glob(1,i)
             ix = ni
             idxn(ix) = ib -1
             DO nj = 1, n_ws1
                j = H_mesh%jjs(nj,ms1)
                DO kj = 1, 3, 2
                   jb =  LA_H%loc_to_glob(kj,j)
                   jx = (kj/2)*n_ws1 + nj
                   jdxn(jx) = jb - 1
                   mat_loc1(ix,jx) = Sij(2*kj-1,nj,ni)
                   mat_loc2(ix,jx) = Sij(2*kj-1,nj,ni)
                END DO
             END DO

             DO nj = 1, n_w2
                j = phi_mesh%jj(nj,m2)
                jb =  LA_phi%loc_to_glob(1,j)
                jx = 2*n_ws1 + nj
                jdxn(jx) = jb - 1 
                mat_loc1(ix,jx) = Phisij(ni,nj)  
                mat_loc2(ix,jx) = Phisij(ni,nj)
             END DO
          END DO
          CALL MatSetValues(H_p_phi_mat1, n_w2, idxn(1:n_w2), 2*n_ws1+n_w2, jdxn(1:2*n_ws1+n_w2), &
               mat_loc1(1:n_w2,1:2*n_ws1+n_w2), ADD_VALUES, ierr)
          CALL MatSetValues(H_p_phi_mat2, n_w2, idxn(1:n_w2), 2*n_ws1+n_w2, jdxn(1:2*n_ws1+n_w2), &
               mat_loc2(1:n_w2,1:2*n_ws1+n_w2), ADD_VALUES, ierr)
       END IF
       !FIN TEST

    ENDDO


    !=========================================================
    !--- Artificial boundary condition: d(phi)/dR + (1/R)*phi = 0
    !=========================================================

    IF (.NOT.PRESENT(index_fourier) .OR. .NOT.PRESENT(R_fourier)) RETURN
    IF (R_fourier.GT.0.d0) THEN 
       !WRITE(*,*) ' Assembling the Fourier condition'
       DO ms = 1, phi_mesh%mes
          IF (phi_mesh%sides(ms) /= index_fourier) CYCLE ! Not on the artificial boundary

          Phisij = 0.d0

          DO ls = 1, phi_mesh%gauss%l_Gs

             !===Compute radius of Gauss point
             ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms))* phi_mesh%gauss%wws(:,ls))

             x = c_mu_phi*rjs(ls,ms)*ray/R_fourier

             DO ni=1,  phi_mesh%gauss%n_ws
                DO nj=1,  phi_mesh%gauss%n_ws
                   Phisij(ni,nj) = Phisij(ni,nj) + x*wws(ni,ls)*wws(nj,ls) 
                ENDDO
             ENDDO

          ENDDO


          DO ni = 1, phi_mesh%gauss%n_ws
             i =  phi_mesh%jjs(ni,ms)
             ib = LA_phi%loc_to_glob(1,i)
             idxn(ni) = ib - 1
             DO nj = 1, phi_mesh%gauss%n_ws
                j = phi_mesh%jjs(nj,ms)
                jb = LA_phi%loc_to_glob(1,j)
                jdxn(nj) = jb - 1
             END DO
          END DO
          CALL MatSetValues(H_p_phi_mat1, n_ws2, idxn(1:n_ws2), n_ws2, jdxn(1:n_ws2), &
               Phisij(1:n_ws2,1:n_ws2), ADD_VALUES, ierr)
          CALL MatSetValues(H_p_phi_mat2, n_ws2, idxn(1:n_ws2), n_ws2, jdxn(1:n_ws2), &
               Phisij(1:n_ws2,1:n_ws2), ADD_VALUES, ierr)
       END DO
    END IF

    CALL MatAssemblyBegin(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyBegin(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)

!!$    DEALLOCATE(mat_loc1, mat_loc2, idxn, jdxn)

  END SUBROUTINE mat_H_p_phi_maxwell

  SUBROUTINE mat_dirichlet_maxwell(H_mesh, jj_v_to_H, Dirichlet_bdy_H_sides, &
       mode, stab, LA_H, H_p_phi_mat1, H_p_phi_mat2, sigma_np, sigma)
    USE def_type_mesh
    USE Dir_nodes
    USE gauss_points
    USE boundary
    USE my_util
    USE input_data ! MODIFICATION: to call sigma_min and mu_min
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE    
    TYPE(mesh_type),            INTENT(IN)    :: H_mesh
    INTEGER,      DIMENSION(:), INTENT(IN)    :: jj_v_to_H
    INTEGER,      DIMENSION(:), INTENT(IN)    :: Dirichlet_bdy_H_sides
    INTEGER,                    INTENT(IN)    :: mode   
    REAL(KIND=8), DIMENSION(3), INTENT(IN)    :: stab
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: sigma_np
    REAL(KIND=8), DIMENSION(H_mesh%me), INTENT(IN) :: sigma

    INTEGER :: ms, ls, ni, nj, i, j, &
         n_ws1, n_w1, m1, ki, kj, ib, jb
    REAL(KIND=8) :: x, y, hm1
    REAL(KIND=8) :: ray, error, stab_colle_h_mu
    REAL(KIND=8), DIMENSION(9,H_mesh%gauss%n_w,H_mesh%gauss%n_w)      :: Hsij
    ! MATRICES POUR LES TERMES DE BORDS Hsij et Phisij
    !=================================================
    ! (--------------------------------------------------------------------)
    ! ( Hsij(1)        |        Hsij(2) | Hsij(4)        || Sij(1)         )
    ! (        Hsij(1) | Hsij(3)        |        Hsij(4) ||        Sij(2)  )
    ! (--------------------------------------------------------------------)
    ! (                | Hsij(5)        |                ||        Sij(3)  )
    ! (                |        Hsij(5) |                || Sij(4)         )
    ! (--------------------------------------------------------------------)
    ! ( Hsij(7)        |        Hsij(9) | Hsij(6)        || Sij(5)         )             
    ! (        Hsij(7) | Hsij(8)        |        Hsij(6) ||        Sij(6)  ) 
    ! (====================================================================)
    ! ( Sij'(1)        |        Sij'(3) | Sij'(5)        || Phisij         )
    ! (        Sij'(2) | Sij'(4)        |        Sij'(6) ||        Phisij  )
    ! (------------------------------------------------------------------- )
    !
    ! L'autre partie des termes croises est la symetrique de la premiere
    ! juste apres le calcsrhs_maul du terme de bord dissymetrique    
    !June 8 2008
    REAL(KIND=8) :: c_sym=.0d0 ! Symmetrization of the bilinear form
    !June 8 2008
!!$ FL + CN 22/03 2013
!!$    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: mat_loc1, mat_loc2
!!$    INTEGER     , DIMENSION(:),   ALLOCATABLE   :: idxn, jdxn
    REAL(KIND=8), DIMENSION(3*H_mesh%gauss%n_w,3*H_mesh%gauss%n_w) :: mat_loc1, mat_loc2
    INTEGER     , DIMENSION(3*H_mesh%gauss%n_w) :: idxn, jdxn
!!$ FL + CN 22/03 2013
    TYPE(petsc_csr_LA)                          :: LA_H
    INTEGER                                     :: ix, jx
    INTEGER                                     :: count
    PetscErrorCode                              :: ierr
    Mat                                         :: H_p_phi_mat1, H_p_phi_mat2

    !June 2009, JLG, CN, Normalization
    stab_colle_H_mu = stab(3)
    !Jan 2010, JLG, CN, Normalization,

    !*********************************************************************************
    !--------------------TERMS ON DIRICHLET BOUNDARY-----------------------------
    !**********************************************************************************
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w

!!$    ALLOCATE(mat_loc1(3*n_w1,3*n_w1))
!!$    ALLOCATE(mat_loc2(3*n_w1,3*n_w1))
!!$    ALLOCATE(idxn(3*n_w1))
!!$    ALLOCATE(jdxn(3*n_w1))

    error = 0
    DO count = 1, SIZE(Dirichlet_bdy_H_sides)
       ms = Dirichlet_bdy_H_sides(count)
       !hm1 = stab_colle_H_mu/SUM(H_mesh%gauss%rjs(:,ms))  
       hm1 = stab_colle_H_mu/(SUM(H_mesh%gauss%rjs(:,ms))*inputs%sigma_min) ! MODIFICATION: normalization for dirichlet term LHS
       m1 = H_mesh%neighs(ms)
       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC H----------------------------
       !====================================================================================

       !-------------------------------hm1 (bi x ni) . (bj x nj)----------------------------
       !====================================================================================

       Hsij = 0.d0  
       DO ls = 1, H_mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          x = hm1*H_mesh%gauss%rjs(ls,ms)*ray 

          DO ni = 1, H_mesh%gauss%n_ws
             DO nj = 1, H_mesh%gauss%n_ws
                y = x * H_mesh%gauss%wws(ni,ls)*H_mesh%gauss%wws(nj,ls)
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(H_mesh%gauss%rnorms(2,ls,ms)**2) 
                Hsij(4,ni,nj) = Hsij(4,ni,nj) - y*H_mesh%gauss%rnorms(1,ls,ms)*H_mesh%gauss%rnorms(2,ls,ms)
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y                                                
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(H_mesh%gauss%rnorms(1,ls,ms)**2)
             ENDDO
          ENDDO

       ENDDO


       !TEST
       !Hsij = 0.d0
       !Hsij = Hsij / hm1
       !TEST
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3 
          DO ni = 1, n_ws1 
             i = H_mesh%jjs(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = ni + (ki-1)*n_ws1
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = H_mesh%jjs(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = nj + (kj-1)*n_ws1
                   jdxn(jx) = jb - 1
                   IF  ((ki == 1) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(1,ni,nj)
                      mat_loc2(ix,jx) = Hsij(1,ni,nj)
                   ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,ni,nj)
                      mat_loc2(ix,jx) = Hsij(4,ni,nj)
                   ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(4,nj,ni)
                      mat_loc2(ix,jx) = Hsij(4,nj,ni)
                   ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)
                   ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(6,ni,nj)
                      mat_loc2(ix,jx) = Hsij(6,ni,nj)
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)

       !====================================================================================
       !------------------------(1/sigma) (Rot bj) . (bi x ni)------------------------------
       !====================================================================================

       !JLG+FL: Jan 18 2013
       !There was a bug on the sign of the normal
       !The sign before rnorms has been changed everywhere in this loop.
       Hsij = 0.d0

       DO ls = 1, H_mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))* H_mesh%gauss%wws(:,ls))

          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             x = H_mesh%gauss%rjs(ls,ms)*ray/sigma(m1)
          ELSE
             x = H_mesh%gauss%rjs(ls,ms)*ray/SUM(sigma_np(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          END IF

          !terme sans derivees
          DO ni = 1,n_ws1
             DO nj = 1, n_ws1
                y = x*H_mesh%gauss%wws(ni,ls)*H_mesh%gauss%wws(nj,ls)
                Hsij(2,ni,nj) = Hsij(2,ni,nj) + y * (-mode/ray)*(rnorms(1,ls,ms))
                Hsij(3,ni,nj) = Hsij(3,ni,nj) + y *   mode/ray *(rnorms(1,ls,ms))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + y * (-1/ray)   *(rnorms(1,ls,ms))
                Hsij(8,ni,nj) = Hsij(8,ni,nj) + y * (-mode/ray)*(rnorms(2,ls,ms))
                Hsij(9,ni,nj) = Hsij(9,ni,nj) + y *   mode/ray *(rnorms(2,ls,ms))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = H_mesh%jjs(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = ni + (ki-1)*n_ws1
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = H_mesh%jjs(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = nj + (kj-1)*n_ws1
                   jdxn(jx) = jb - 1
                   IF  ( (ki == 2) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(2,ni,nj)
                      mat_loc2(ix,jx) = Hsij(3,ni,nj)
                   ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)
                   ELSEIF  ( (ki == 2) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(8,ni,nj)
                      mat_loc2(ix,jx) = Hsij(9,ni,nj) 
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)

       !Feb 2 2007
       Hsij=c_sym*Hsij !SYM
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_ws1
             i = H_mesh%jjs(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = ni + (ki-1)*n_ws1
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_ws1
                   j = H_mesh%jjs(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = nj + (kj-1)*n_ws1
                   jdxn(jx) = jb - 1
                   IF  ( (kj == 2) .AND. (ki == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(2,nj,ni)
                      mat_loc2(ix,jx) = Hsij(3,nj,ni)
                   ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                      mat_loc1(ix,jx) = Hsij(5,nj,ni)
                      mat_loc2(ix,jx) = Hsij(5,nj,ni)
                   ELSEIF  ( (kj == 2) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(8,nj,ni)
                      mat_loc2(ix,jx) = Hsij(9,nj,ni)
                   ENDIF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_ws1,1:3*n_ws1), ADD_VALUES, ierr)
       !feb 2 2007


       Hsij = 0.d0

       DO ls = 1, H_mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))* H_mesh%gauss%wws(:,ls))

          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             x = H_mesh%gauss%rjs(ls,ms)*ray/(sigma(m1))
          ELSE
             x = H_mesh%gauss%rjs(ls,ms)*ray/(SUM(sigma_np(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls)))
          END IF

          !termes avec derivees
          DO ni = 1,n_ws1
             y = x*H_mesh%gauss%wws(ni,ls)
             DO nj = 1, n_w1
                Hsij(1,ni,nj) = Hsij(1,ni,nj) + y*(-H_mesh%gauss%dw_s(2,nj,ls,ms))*(rnorms(2,ls,ms))
                Hsij(4,ni,nj) = Hsij(4,ni,nj) + y*  H_mesh%gauss%dw_s(1,nj,ls,ms) *(rnorms(2,ls,ms))
                Hsij(5,ni,nj) = Hsij(5,ni,nj) + &
                     y*(-H_mesh%gauss%dw_s(2,nj,ls,ms)*(rnorms(2,ls,ms))-H_mesh%gauss%dw_s(1,nj,ls,ms)*(rnorms(1,ls,ms)))
                Hsij(6,ni,nj) = Hsij(6,ni,nj) + y*(-H_mesh%gauss%dw_s(1,nj,ls,ms))*(rnorms(1,ls,ms))
                Hsij(7,ni,nj) = Hsij(7,ni,nj) + y*  H_mesh%gauss%dw_s(2,nj,ls,ms) *(rnorms(1,ls,ms))
             ENDDO
          ENDDO

       ENDDO

       !TEST
       !Hsij = 0.d0
       !TEST


       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki= 1, 3  
          DO ni = 1, n_ws1
             i = H_mesh%jjs(ni,ms)
             ib = LA_H%loc_to_glob(ki,i)
             ix = ni + (ki-1)*n_ws1
             idxn(ix) = ib - 1
             DO kj = 1, 3
                DO nj = 1, n_w1
                   j = H_mesh%jj(nj,m1)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = nj + (kj-1)*n_w1
                   jdxn(jx) = jb - 1
                   IF  ((ki == 1) .AND. (kj == 1))  THEN
                      mat_loc1(ix,jx) = Hsij(1,ni,nj)
                      mat_loc2(ix,jx) = Hsij(1,ni,nj)
                   ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,ni,nj)
                      mat_loc2(ix,jx) = Hsij(4,ni,nj)
                   ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                      mat_loc1(ix,jx) = Hsij(5,ni,nj)
                      mat_loc2(ix,jx) = Hsij(5,ni,nj)
                   ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                      mat_loc1(ix,jx) = Hsij(6,ni,nj)   
                      mat_loc2(ix,jx) = Hsij(6,ni,nj) 
                   ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(7,ni,nj)
                      mat_loc2(ix,jx) = Hsij(7,ni,nj)
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc1(1:3*n_ws1,1:3*n_w1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_ws1 , idxn(1:3*n_ws1), 3*n_w1, jdxn(1:3*n_w1), &
            mat_loc2(1:3*n_ws1,1:3*n_w1), ADD_VALUES, ierr)

       !Feb 2 2007
       Hsij=c_sym*Hsij !SYM
       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ki = 1, 3
          DO ni = 1, n_w1
             i = H_mesh%jj(ni,m1)
             ib = LA_H%loc_to_glob(ki,i)
             ix = ni + (ki-1)*n_w1
             idxn(ix) = ib - 1
             DO kj= 1, 3
                DO nj = 1, n_ws1
                   j = H_mesh%jjs(nj,ms)
                   jb = LA_H%loc_to_glob(kj,j)
                   jx = nj + (kj-1)*n_ws1
                   jdxn(jx) = jb - 1
                   IF  ((kj == 1) .AND. (ki == 1))  THEN
                      mat_loc1(ix,jx) = Hsij(1,nj,ni)
                      mat_loc2(ix,jx) = Hsij(1,nj,ni)
                   ELSEIF  ((kj == 1) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(4,nj,ni)
                      mat_loc2(ix,jx) = Hsij(4,nj,ni)
                   ELSEIF  ((kj == 2) .AND. (ki == 2))  THEN
                      mat_loc1(ix,jx) = Hsij(5,nj,ni)
                      mat_loc2(ix,jx) = Hsij(5,nj,ni)
                   ELSEIF  ((kj == 3) .AND. (ki == 3)) THEN
                      mat_loc1(ix,jx) = Hsij(6,nj,ni)
                      mat_loc2(ix,jx) = Hsij(6,nj,ni)
                   ELSEIF  ((kj == 3) .AND. (ki == 1)) THEN
                      mat_loc1(ix,jx) = Hsij(7,nj,ni)
                      mat_loc2(ix,jx) = Hsij(7,nj,ni)
                   ENDIF
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 3*n_w1 , idxn(1:3*n_w1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc1(1:3*n_w1,1:3*n_ws1), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 3*n_w1 , idxn(1:3*n_w1), 3*n_ws1, jdxn(1:3*n_ws1), &
            mat_loc2(1:3*n_w1,1:3*n_ws1), ADD_VALUES, ierr)

    ENDDO

    CALL MatAssemblyBegin(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyBegin(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)

!!$    IF (ALLOCATED(mat_loc1)) DEALLOCATE(mat_loc1)
!!$    IF (ALLOCATED(mat_loc2)) DEALLOCATE(mat_loc2)
!!$    IF (ALLOCATED(idxn)) DEALLOCATE(idxn)
!!$    IF (ALLOCATED(jdxn)) DEALLOCATE(jdxn)


  END SUBROUTINE mat_dirichlet_maxwell

  SUBROUTINE courant_int_by_parts(H_mesh,phi_mesh,interface_H_phi,sigma,mu_phi,mu_H_field,time,mode,&
       rhs_H,nl, LA_H, LA_phi, vb_1, vb_2, B_ext, sigma_curl_gauss, J_over_sigma_gauss)
    !forcage faisant intervenir J, volumique et interface_H_phi
    !pour le probleme en entier

    USE def_type_mesh
    USE gauss_points
    USE boundary
    USE my_util
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh, phi_mesh
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_phi
    REAL(KIND=8),                          INTENT(IN)   :: mu_phi, time
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: mu_H_field
    INTEGER,                               INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: nl
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: B_ext
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: rhs_H
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: J_over_sigma_gauss !Used only if sigma variable in fluid
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: sigma_curl_gauss !Used only if sigma variable in fluid
    INTEGER                                             :: index
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                :: src_H
    REAL(KIND=8), DIMENSION(phi_mesh%np,2)              :: src_phi
    ! CN possible faute POINTER
    !REAL(KIND=8), DIMENSION(:,:), POINTER               :: src_H, src_phi
    !REAL(KIND=8), DIMENSION(:,:), POINTER               :: nl
    !REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE           :: src_H, src_phi
!!$    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE           :: nl
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_ws,phi_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: gaussp
    REAL(KIND=8) :: ray
    INTEGER :: m, l, i, ni, k, ms, ls, n_ws1, n_ws2, ms1, ms2, H_bloc_size, n_w2, m1
    INTEGER :: mesh_id1
    REAL(KIND=8), DIMENSION(6)             :: JsolH_anal, rhs_Hl
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%n_w) :: dwH
    !REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%n_w) :: dwphi
    !REAL(KIND=8), DIMENSION(2,phi_mesh%gauss%n_w) :: src_phil
    REAL(KIND=8) :: ray_rjl, muhl
    !REAL(KIND=8) :: moderay2
    !REAL(KIND=8) :: tps, dummy 
!!$ FL + CN, 22/03/2013
!!$    INTEGER, DIMENSION(:), ALLOCATABLE          :: idxn
    INTEGER, DIMENSION(H_mesh%np)               :: idxn_H
    INTEGER, DIMENSION(phi_mesh%np)             :: idxn_phi
!!$ FL + CN, 22/03/2013
    TYPE(petsc_csr_LA)                          :: LA_H, LA_phi
    REAL(KIND=8), DIMENSION(6)                  :: B_ext_l
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2

    !ALLOCATE(src_H(H_mesh%np,6), src_phi(phi_mesh%np,2))

    !CALL VecZeroEntries(vb_1, ierr)
    !CALL VecZeroEntries(vb_2, ierr)

    !forcage volumique  
    !attention on comprime le calcul sur les points de Gauss et integration !!
    !j/sigma *(Rot(b))

    !tps = user_time(dummy)
    src_H    = 0.d0
    src_phi  = 0.d0
    index    = 0

    DO m = 1, H_mesh%me
       mesh_id1 = H_mesh%i_d(m) 
       DO l = 1, H_mesh%gauss%l_G
          index = index + 1
          !Feb 8 2007, muhl
          muhl=SUM(mu_H_field(H_mesh%jj(:,m))*H_mesh%gauss%ww(:,l))
          !Feb 8 2007, muhl
          dwH = H_mesh%gauss%dw(:,:,l,m)
          !===Compute radius of Gauss point
          DO k=1, 6
             B_ext_l(k) = SUM(B_ext(H_mesh%jj(:,m),k)*H_mesh%gauss%ww(:,l))
          END DO

          JsolH_anal = 0.d0
          rhs_Hl = 0.d0
          gaussp = 0.d0
          DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
             gaussp = gaussp + H_mesh%rr(:,i)*H_mesh%gauss%ww(ni,l)
             JsolH_anal(:) = JsolH_anal(:) + muhl*NL(i,:)*H_mesh%gauss%ww(ni,l)
             rhs_Hl(:) = rhs_Hl(:) + rhs_H(i,:)*H_mesh%gauss%ww(ni,l)
          ENDDO
          ray = gaussp(1)
          ray_rjl = H_mesh%gauss%rj(l,m)*ray

          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
             DO k = 1, 6
                JsolH_anal(k) =  JsolH_anal(k) + J_over_sigma_gauss(index,k) + sigma_curl_gauss(index,k)
             END DO
          ELSE
             DO k = 1, 6
                !JsolH_anal(k) = muhl*JsolH_anal(k) + & ! BUG Jan 2010, JLG, CN, FL
                JsolH_anal(k) = JsolH_anal(k) + &
                     Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m), muhl, time, mesh_id1, B_ext_l)/sigma(m)
             END DO
          END IF

          DO ni = 1,H_mesh%gauss%n_w

             i = H_mesh%jj(ni,m)

             !--------Composante r------
             src_H(i,1) = src_H(i,1)+ ray_rjl &
                  *(JsolH_anal(3)*dwH(2,ni) &
                  + mode/ray*JsolH_anal(6)*H_mesh%gauss%ww(ni,l) &
                  + rhs_Hl(1)*H_mesh%gauss%ww(ni,l))

             src_H(i,2) = src_H(i,2)+ ray_rjl  &
                  *(JsolH_anal(4)*dwH(2,ni) &
                  - mode/ray*JsolH_anal(5)*H_mesh%gauss%ww(ni,l) &
                  + rhs_Hl(2)*H_mesh%gauss%ww(ni,l))   

             !--------Composante theta------
             src_H(i,3) = src_H(i,3)+ ray_rjl  &
                  * (-JsolH_anal(1)*dwH(2,ni)  &
                  + 1/ray*JsolH_anal(5)*(ray*dwH(1,ni) + H_mesh%gauss%ww(ni,l)) &
                  + rhs_Hl(3)*H_mesh%gauss%ww(ni,l)) 

             src_H(i,4) = src_H(i,4)+ ray_rjl &
                  * (-JsolH_anal(2)*dwH(2,ni) &
                  + 1/ray*JsolH_anal(6)*(ray*dwH(1,ni) + H_mesh%gauss%ww(ni,l)) &
                  + rhs_Hl(4)*H_mesh%gauss%ww(ni,l))

             !--------Composante z------
             src_H(i,5) = src_H(i,5)+ ray_rjl* &
                  (-mode/ray*JsolH_anal(2)*H_mesh%gauss%ww(ni,l) &
                  - JsolH_anal(3)*dwH(1,ni) &
                  + rhs_Hl(5)*H_mesh%gauss%ww(ni,l))

             src_H(i,6) = src_H(i,6)+ ray_rjl* &
                  (mode/ray*JsolH_anal(1)*H_mesh%gauss%ww(ni,l) &
                  - JsolH_anal(4)*dwH(1,ni) &
                  + rhs_Hl(6)*H_mesh%gauss%ww(ni,l)) 
          ENDDO

       END DO
    END DO
    !tps = user_time(dummy)- tps
    !WRITE(*,*) ' Temps in courant boucle me H', tps
    !tps = user_time(dummy)

    ! We integrate by parts this term
    ! JLG + FL, FEB 10, 2010
    !DO m = 1, phi_mesh%me
    !   src_phil=0
    !   DO l = 1, phi_mesh%gauss%l_G
    !      dwphi = phi_mesh%gauss%dw(:,:,l,m)
    !      !===Compute radius of Gauss point
    !      rhs_dphil=0
    !      rhs_phil=0
    !      ray = 0
    !      DO ni = 1, phi_mesh%gauss%n_w;  i = phi_mesh%jj(ni,m)
    !         ray = ray + phi_mesh%rr(1,i)*phi_mesh%gauss%ww(ni,l)
    !         rhs_phil(:) = rhs_phil(:) + rhs_phi(i,:)*phi_mesh%gauss%ww(ni,l)
    !         DO k =1 ,2
    !            rhs_dphil(:,k) = rhs_dphil(:,k) + rhs_phi(i,:)*dwphi(k,ni)
    !         END DO
    !      END DO
    !      ray_rjl = phi_mesh%gauss%rj(l,m)*ray
    !      moderay2 = (mode/ray)**2

    !      DO ni = 1, phi_mesh%gauss%n_w

    !         src_phil(1,ni) =  src_phil(1,ni) + ray_rjl* &
    !              (rhs_dphil(1,1)*dwphi(1,ni) + &
    !              moderay2*rhs_phil(1)*phi_mesh%gauss%ww(ni,l) + &
    !              rhs_dphil(1,2)*dwphi(2,ni))

    !         src_phil(2,ni) =  src_phil(2,ni) + ray_rjl* &
    !              (rhs_dphil(2,1)*dwphi(1,ni) + &
    !              moderay2*rhs_phil(2)*phi_mesh%gauss%ww(ni,l) + &
    !              rhs_dphil(2,2)*dwphi(2,ni))
    !      END DO

    !   END DO
    !   DO ni = 1, phi_mesh%gauss%n_w
    !      i = phi_mesh%jj(ni,m)
    !      src_phi(i,:) = src_phi(i,:) + src_phil(:,ni) 
    !   END DO
    !END DO
    ! End integration by parts
    ! JLG + FL, FEB 10, 2010


    !==interface_H_phi
    !forcage sur l'interface_H_phi
    !attention on comprime le calcul sur les points de Gauss et integration !!
    !j/sigma*(b x nc + grad(phi) x nv)

    CALL gauss(phi_mesh)

    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = phi_mesh%gauss%n_ws
    n_w2 = phi_mesh%gauss%n_w

    H_bloc_size = H_mesh%np

    IF (interface_H_phi%mes /=0) THEN  ! Ajout du test pour les grands nb de domaines
       IF (H_mesh%gauss%n_ws == n_ws) THEN
          w_cs = wws
       ELSE    
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(1,ls)+0.5*wws(3,ls)
             w_cs(2,ls)= wws(2,ls)+0.5*wws(3,ls)
             w_cs(3,ls)=0
          ENDDO
       END IF
    END IF


    !WRITE(*,*) ' Courant: init gauss'
    DO ms = 1, interface_H_phi%mes

       ms2 = interface_H_phi%mesh2(ms)
       ms1 = interface_H_phi%mesh1(ms)
       m = phi_mesh%neighs(ms2)           
       m1 = H_mesh%neighs(ms1)
       mesh_id1 = H_mesh%i_d(m1)
       DO ls = 1,l_Gs
          !Feb 9 2007, muhl
          muhl=SUM(mu_H_field(interface_H_phi%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
          !Feb 9 2007, muhl
          DO k=1, 6
             B_ext_l(k) = SUM(B_ext(interface_H_phi%jjs1(1:n_ws1,ms),k)*w_cs(1:n_ws1,ls))
          END DO

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,interface_H_phi%mesh2(ms))
             ray = ray + phi_mesh%rr(1,i)* wws(ni,ls)
          END DO

          gaussp = 0.d0    
          DO ni=1, n_ws2
             i=phi_mesh%jjs(ni,ms2)
             gaussp = gaussp + phi_mesh%rr(:,i)*phi_mesh%gauss%wws(ni,ls)
          ENDDO

          DO k=1, 6
             JsolH_anal(k) = Jexact_gauss(k, gaussp, mode, mu_phi ,sigma(m1), &
                  muhl, time, mesh_id1, B_ext_l)/sigma(m1) &
                  + muhl * SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
          ENDDO
!!$! TO DO : to do before using H with phi 
!!$          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
!!$             DO k = 1, 6
!!$                JsolH_anal(k) = J_over_sigma_gauss(k) + sigma_curl(index,k) &
!!$                  + SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
!!$             END DO
!!$          ELSE
!!$             DO k = 1, 6
!!$                JsolH_anal(k) = Jexact_gauss(k, gaussp, mode, mu_phi ,sigma(m1), &
!!$                  muhl, time, mesh_id1, B_ext_l)/sigma(m1) &
!!$                  + SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
!!$             END DO
!!$          END IF
!!$! TO DO

          !---------forcage pour H            

          DO ni=1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             src_H(i,1) = src_H(i,1)+rjs(ls,ms2)*ray*( &
                  -JsolH_anal(3)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)))

             src_H(i,2) = src_H(i,2)+rjs(ls,ms2)*ray*( &
                  -JsolH_anal(4)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)))

             src_H(i,3) = src_H(i,3)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(1)*w_cs(ni,ls)*(-rnorms(2,ls,ms2)) &
                  -JsolH_anal(5)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,4) = src_H(i,4)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(2)*w_cs(ni,ls)*(-rnorms(2,ls,ms2))  &
                  -JsolH_anal(6)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,5) = src_H(i,5)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(3)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))

             src_H(i,6) = src_H(i,6)+rjs(ls,ms2)*ray*( &
                  JsolH_anal(4)*w_cs(ni,ls)*(-rnorms(1,ls,ms2)))
          ENDDO

          !---------forcage pour phi            
          !terme sans derivee de phi
          DO ni=1,n_ws2           
             i = interface_H_phi%jjs2(ni,ms)
             !attention si on force sur l'axe, il faut retirer les 1/ray
             !There was a BUG here. There was w_cs instead of wws
             src_phi(i,1) = src_phi(i,1)+rjs(ls,ms2)*( &
                  - mode*JsolH_anal(2)*wws(ni,ls) * rnorms(2,ls,ms2) &
                  + mode*JsolH_anal(6)*wws(ni,ls) * rnorms(1,ls,ms2))

             src_phi(i,2) = src_phi(i,2)+rjs(ls,ms2)*( &
                  + mode*JsolH_anal(1)*wws(ni,ls) * rnorms(2,ls,ms2) &
                  - mode*JsolH_anal(5)*wws(ni,ls) * rnorms(1,ls,ms2))

          ENDDO

          !terme avec derivee de phi
          DO ni=1,n_w2           
             i = phi_mesh%jj(ni,m)
             src_phi(i,1) = src_phi(i,1)+rjs(ls,ms2)*ray*( &
                  + JsolH_anal(3) *(dw_s(2,ni,ls,ms2) * rnorms(1,ls,ms2)&
                  -dw_s(1,ni,ls,ms2) * rnorms(2,ls,ms2)))  

             src_phi(i,2) = src_phi(i,2)+rjs(ls,ms2)*ray*( &
                  + JsolH_anal(4)*(dw_s(2,ni,ls,ms2) * rnorms(1,ls,ms2)&
                  -dw_s(1,ni,ls,ms2) * rnorms(2,ls,ms2)))

          ENDDO

          ! Integration by parts of int(GRAD rhs_phi Grad psi)
          rhs_Hl = 0.d0
          DO ni=1, n_ws1
             i = interface_H_phi%jjs1(ni,ms)
             rhs_Hl(:) = rhs_Hl(:) + rhs_H(i,:)*w_cs(ni,ls)
          ENDDO

          DO ni=1, n_ws2
             i = interface_H_phi%jjs2(ni,ms)
             src_phi(i,1) = src_phi(i,1)+rjs(ls,ms2)*ray*wws(ni,ls)*( &
                  rhs_Hl(1)*rnorms(1,ls,ms2) + rhs_Hl(5)*rnorms(2,ls,ms2))
             src_phi(i,2) = src_phi(i,2)+rjs(ls,ms2)*ray*wws(ni,ls)*( &
                  rhs_Hl(2)*rnorms(1,ls,ms2) + rhs_Hl(6)*rnorms(2,ls,ms2))
          END DO
          ! End integration by parts of int(GRAD rhs_phi Grad psi)
       END DO
    END DO
    !tps = user_time(dummy)- tps
    !WRITE(*,*) ' Courant: init interface_H_phi'

    IF (H_mesh%np /= 0) THEN
!!$       ALLOCATE(idxn(H_mesh%np))
       idxn_H = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,2), ADD_VALUES, ierr)
       idxn_H = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,3), ADD_VALUES, ierr)
       idxn_H = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,6), ADD_VALUES, ierr)
!!$       DEALLOCATE(idxn)
    END IF
    IF (phi_mesh%np /=0) THEN
!!$       ALLOCATE(idxn(phi_mesh%np))
       idxn_phi = LA_phi%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, phi_mesh%np, idxn_phi, src_phi(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, phi_mesh%np, idxn_phi, src_phi(:,2), ADD_VALUES, ierr)       
!!$       DEALLOCATE(idxn)
    END IF

    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)

!!$    IF (H_mesh%me /=0) THEN
!!$       DEALLOCATE(nl)
!!$    END IF
    !DEALLOCATE(src_H, src_phi)

  END SUBROUTINE courant_int_by_parts

  !===JLG Jan 22 2018
  !SUBROUTINE surf_int(H_mesh,phi_mesh, interface_H_phi, interface_H_mu, &
  SUBROUTINE surf_int(H_mesh, phi_mesh, pmag_mesh, interface_H_phi, interface_H_mu, &
       list_dirichlet_sides_H, sigma,mu_phi, mu_H_field, time, mode, &
       LA_H, LA_phi, LA_pmag, vb_1, vb_2, sigma_tot_gauss, R_fourier, index_fourier)
    !calcul du forcage a la frontiere exterieure
    USE my_util
    USE def_type_mesh
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)    :: H_mesh, phi_mesh, pmag_mesh
    TYPE(interface_type),         INTENT(IN)    :: interface_H_phi, interface_H_mu
    INTEGER,     DIMENSION(:),    INTENT(IN)    :: list_dirichlet_sides_H
    REAL(KIND=8),                 INTENT(IN)    :: mu_phi, time
    REAL(KIND=8),DIMENSION(H_mesh%me),INTENT(IN):: sigma
    REAL(KIND=8),DIMENSION(:),    INTENT(IN)    :: mu_H_field 
    INTEGER,                      INTENT(IN)    :: mode
    REAL(KIND=8), DIMENSION(:,:)                :: sigma_tot_gauss
    REAL(KIND=8),                 OPTIONAL      :: R_fourier
    INTEGER,                      OPTIONAL      :: index_fourier
    REAL(KIND=8), DIMENSION(H_mesh%np, 6)       :: src_H 
    REAL(KIND=8), DIMENSION(phi_mesh%np, 2)     :: src_phi
    REAL(KIND=8), DIMENSION(pmag_mesh%np, 2)    :: src_pmag !===JLG Jan 22 2018
    !REAL(KIND=8), DIMENSION(pmag_mesh%gauss%n_ws,H_mesh%gauss%l_Gs) :: wwps !===JLG Jan 22 2018
    REAL(KIND=8), DIMENSION(4,H_mesh%gauss%l_Gs):: Banal
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%l_Gs):: rloc
    REAL(KIND=8), DIMENSION(H_mesh%gauss%l_Gs)  :: B_dot_n_cos, B_dot_n_sin, muloc
    REAL(KIND=8), DIMENSION(4,pmag_mesh%gauss%l_Gs):: Banal_pmesh
    REAL(KIND=8), DIMENSION(2,pmag_mesh%gauss%l_Gs):: rloc_pmesh
    REAL(KIND=8), DIMENSION(pmag_mesh%gauss%l_Gs)  :: B_dot_n_cos_pmesh, B_dot_n_sin_pmesh, muloc_pmesh
    REAL(KIND=8)                                :: ray, muhl, norm, y
    INTEGER                                     :: ms, ls, ns, i, k, m, n, ni, count, index
    REAL(KIND=8), DIMENSION(2)                  :: gaussp
    REAL(KIND=8), DIMENSION(6)                  :: EsolH_anal, Esolphi_anal
    INTEGER, DIMENSION(H_mesh%np)               :: idxn_H
    INTEGER, DIMENSION(phi_mesh%np)             :: idxn_phi
    INTEGER, DIMENSION(pmag_mesh%np)            :: idxn_pmag !===JLG Jan 22 2018
    TYPE(petsc_csr_LA)                          :: LA_H, LA_phi, LA_pmag !===JLG Jan 22 2018
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2

    src_H = 0.d0
    src_phi = 0.d0
    src_pmag = 0.d0

!!$    !LC 2019/04/29: error dimension pmag_mesh%jjs and H_mesh%jjs
!!$    !===JLG Jan 22 2018
!!$    !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
!!$    IF (H_mesh%gauss%n_ws==2) THEN
!!$       wwps=H_mesh%gauss%wws
!!$    ELSE
!!$       wwps(1,:)   = H_mesh%gauss%wws(1,:) + 0.5d0*H_mesh%gauss%wws(3,:)
!!$       wwps(2,:)   = H_mesh%gauss%wws(2,:) + 0.5d0*H_mesh%gauss%wws(3,:)
!!$    END IF
!!$    !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
!!$    !===JLG Jan 22 2018
!!$    !LC 2019/04/29: error dimension pmag_mesh%jjs and H_mesh%jjs

    index = 0
    DO count = 1, SIZE(Neumann_bdy_H_sides)
       ms = Neumann_bdy_H_sides(count)
       m = H_mesh%neighs(ms)
       !===JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)
       DO ls = 1, H_mesh%gauss%l_Gs
          muhl = SUM(mu_H_field(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          muloc(ls) = muhl
          rloc(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls)) 
          rloc(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls)) 
       END DO
       Banal(1,:) = muhl*Hexact(H_mesh, 1, rloc, mode, muloc, time)
       Banal(2,:) = muhl*Hexact(H_mesh, 2, rloc, mode, muloc, time)
       Banal(3,:) = muhl*Hexact(H_mesh, 5, rloc, mode, muloc, time)
       Banal(4,:) = muhl*Hexact(H_mesh, 6, rloc, mode, muloc, time)
       B_dot_n_cos = Banal(1,:)*H_mesh%gauss%rnorms(1,:,ms) + Banal(3,:)*H_mesh%gauss%rnorms(2,:,ms)
       B_dot_n_sin = Banal(2,:)*H_mesh%gauss%rnorms(1,:,ms) + Banal(4,:)*H_mesh%gauss%rnorms(2,:,ms)
       !===JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)

       DO ls = 1, H_mesh%gauss%l_Gs
          index = index + 1
          !Feb 8 2007, mmuhl
          muhl = SUM(mu_H_field(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          !Feb 8 2007, mmuhl

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, H_mesh%gauss%n_ws;  i = H_mesh%jjs(ni,ms)
             ray = ray + H_mesh%rr(1,i)* H_mesh%gauss%wws(ni,ls)
          END DO

          IF (ray.LT.1.d-12*H_mesh%global_diameter) CYCLE !ATTENTION Axe

          gaussp = 0.d0    
          DO ns=1, H_mesh%gauss%n_ws
             i=H_mesh%jjs(ns,ms)
             gaussp = gaussp + H_mesh%rr(:,i)*H_mesh%gauss%wws(ns,ls)
          ENDDO

!!$          DO k=1, 6
!!$             EsolH_anal(k) = Eexact_gauss(k,gaussp,mode,mu_phi,sigma(m),muhl, time)
!!$          ENDDO
          !LC-JLG-CN 2018/04       
          DO k=1, 6
             EsolH_anal(k) = Eexact_gauss(k,gaussp,mode,mu_phi,sigma_tot_gauss(index,mod(k+1,2)+1),muhl, time)
          ENDDO
          !LC-JLG-CN 2018/04

!!$          !LC 2019/04/29: error dimension pmag_mesh%jjs and H_mesh%jjs
!!$          !===JLG Jan 22 2018
!!$          !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
!!$          DO ns=1, pmag_mesh%gauss%n_ws
!!$             i = pmag_mesh%jjs(ns,ms)
!!$             src_pmag(i,1) = src_pmag(i,1) - inputs%stab(1)*B_dot_n_cos(ls)*wwps(ns,ls)*H_mesh%gauss%rjs(ls,ms)*ray
!!$             src_pmag(i,2) = src_pmag(i,2) - inputs%stab(1)*B_dot_n_sin(ls)*wwps(ns,ls)*H_mesh%gauss%rjs(ls,ms)*ray
!!$          END DO
!!$          !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
!!$          !===JLG Jan 22 2018
!!$          !LC 2019/04/29: error dimension pmag_mesh%jjs and H_mesh%jjs

          !===Forcing at the boundary
          !=== - E.(b x nc)
          DO ns=1, H_mesh%gauss%n_ws
             i = H_mesh%jjs(ns,ms)
             src_H(i,1) = src_H(i,1)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  -EsolH_anal(3)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)))

             src_H(i,2) = src_H(i,2)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  -EsolH_anal(4)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)))

             src_H(i,3) = src_H(i,3)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(1)*H_mesh%gauss%wws(ns,ls)* &
                  (H_mesh%gauss%rnorms(2,ls,ms)) - &
                  EsolH_anal(5)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,4) = src_H(i,4)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(2)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(2,ls,ms)) - &
                  EsolH_anal(6)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,5) = src_H(i,5)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(3)*H_mesh%gauss%wws(ns,ls)*  &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

             src_H(i,6) = src_H(i,6)-H_mesh%gauss%rjs(ls,ms)*ray*( &
                  EsolH_anal(4)*H_mesh%gauss%wws(ns,ls) * &
                  (H_mesh%gauss%rnorms(1,ls,ms)))

          ENDDO

          !=== JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)
          norm = (inputs%stab(1)/inputs%Rem)*(SUM(H_mesh%gauss%rjs(:,ms))/H_mesh%global_diameter)**(2*alpha-1)&
               /(H_mesh%global_diameter*inputs%sigma_min*inputs%mu_min**2)
          !TEST JLG RZ, Jan 22 2018  Removed stabilization
          !norm =0.d0
          !TEST JLG RZ, Jan 22 2018  Removed stabilization
          DO ns = 1, H_mesh%gauss%n_ws
             i = H_mesh%jjs(ns,ms)
             y = norm*muhl*H_mesh%gauss%rjs(ls,ms)*ray
             src_H(i,1) = src_H(i,1) + y*B_dot_n_cos(ls)*H_mesh%gauss%wws(ns,ls)*H_mesh%gauss%rnorms(1,ls,ms)
             src_H(i,2) = src_H(i,2) + y*B_dot_n_sin(ls)*H_mesh%gauss%wws(ns,ls)*H_mesh%gauss%rnorms(1,ls,ms)
             src_H(i,5) = src_H(i,5) + y*B_dot_n_cos(ls)*H_mesh%gauss%wws(ns,ls)*H_mesh%gauss%rnorms(2,ls,ms)
             src_H(i,6) = src_H(i,6) + y*B_dot_n_sin(ls)*H_mesh%gauss%wws(ns,ls)*H_mesh%gauss%rnorms(2,ls,ms)
          END DO
          !=== JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)

       ENDDO
    ENDDO

    !===Neumann boundary pmag_mesh
    DO count = 1, SIZE(Neumann_bdy_pmag_sides)
       ms = Neumann_bdy_pmag_sides(count)
       m = pmag_mesh%neighs(ms)
       !===JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)
       DO ls = 1, pmag_mesh%gauss%l_Gs
          muhl = SUM(mu_H_field(pmag_mesh%jjs(:,ms))*pmag_mesh%gauss%wws(:,ls))
          muloc_pmesh(ls) = muhl
          rloc_pmesh(1,ls) = SUM(pmag_mesh%rr(1,pmag_mesh%jjs(:,ms))*pmag_mesh%gauss%wws(:,ls)) 
          rloc_pmesh(2,ls) = SUM(pmag_mesh%rr(2,pmag_mesh%jjs(:,ms))*pmag_mesh%gauss%wws(:,ls)) 
       END DO
       Banal_pmesh(1,:) = muhl*Hexact(pmag_mesh, 1, rloc_pmesh, mode, muloc_pmesh, time)
       Banal_pmesh(2,:) = muhl*Hexact(pmag_mesh, 2, rloc_pmesh, mode, muloc_pmesh, time)
       Banal_pmesh(3,:) = muhl*Hexact(pmag_mesh, 5, rloc_pmesh, mode, muloc_pmesh, time)
       Banal_pmesh(4,:) = muhl*Hexact(pmag_mesh, 6, rloc_pmesh, mode, muloc_pmesh, time)
       B_dot_n_cos_pmesh = Banal_pmesh(1,:)*pmag_mesh%gauss%rnorms(1,:,ms) &
            + Banal_pmesh(3,:)*pmag_mesh%gauss%rnorms(2,:,ms)
       B_dot_n_sin_pmesh = Banal_pmesh(2,:)*pmag_mesh%gauss%rnorms(1,:,ms) &
            + Banal_pmesh(4,:)*pmag_mesh%gauss%rnorms(2,:,ms)
       !===JLG+CN July 20 2017. Non-homogeneous B.n (due to divergence stabilization)

       DO ls = 1, pmag_mesh%gauss%l_Gs
          !Feb 8 2007, mmuhl
          muhl = SUM(mu_H_field(pmag_mesh%jjs(:,ms))*pmag_mesh%gauss%wws(:,ls))
          !Feb 8 2007, mmuhl

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, pmag_mesh%gauss%n_ws;  i = pmag_mesh%jjs(ni,ms)
             ray = ray + pmag_mesh%rr(1,i)* pmag_mesh%gauss%wws(ni,ls)
          END DO

          IF (ray.LT.1.d-12*H_mesh%global_diameter) CYCLE !ATTENTION Axe

          !===JLG Jan 22 2018
          !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
          DO ns=1, pmag_mesh%gauss%n_ws
             i = pmag_mesh%jjs(ns,ms)
             src_pmag(i,1) = src_pmag(i,1) - inputs%stab(1)*B_dot_n_cos_pmesh(ls)* &
                  pmag_mesh%gauss%wws(ns,ls)*pmag_mesh%gauss%rjs(ls,ms)*ray
             src_pmag(i,2) = src_pmag(i,2) - inputs%stab(1)*B_dot_n_sin_pmesh(ls)* &
                  pmag_mesh%gauss%wws(ns,ls)*pmag_mesh%gauss%rjs(ls,ms)*ray
          END DO
          !===Non-homogeneous Neumann BC. -\int_{\Gamma_N} (B_anal.n) q ds
          !===JLG Jan 22 2018
       ENDDO
    ENDDO

    !===Neumann boundary phi_mesh
    DO count = 1, SIZE(Neumann_bdy_phi_sides)
       ms = Neumann_bdy_phi_sides(count)
       m = phi_mesh%neighs(ms)
       DO ls = 1, phi_mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, phi_mesh%gauss%n_ws;  i = phi_mesh%jjs(ni,ms)
             ray = ray + phi_mesh%rr(1,i)* phi_mesh%gauss%wws(ni,ls)
          END DO

          gaussp = 0.d0    
          DO ns=1, phi_mesh%gauss%n_ws
             i=phi_mesh%jjs(ns,ms)
             gaussp = gaussp + phi_mesh%rr(:,i)*phi_mesh%gauss%wws(ns,ls)
          ENDDO

          DO k=1, 6
             !===Here sigma and mu_H_field should not intervene
             !===I put boggus values for sigma and mu_H_field, Feb 8 2007, Jean-Luc Guermond       
             Esolphi_anal(k) = Eexact_gauss(k,gaussp,mode,mu_phi,sigma(1),mu_H_field(1),time)
          ENDDO
          !TO DEBUG

          !TO DEBUG

          !===Nemnann forcing for phi in rhs:  - E.(grad(phi) x nv)
          DO ns=1, phi_mesh%gauss%n_ws           
             i = phi_mesh%jjs(ns,ms) 
             DO n = 1, phi_mesh%gauss%n_w
                IF (phi_mesh%jj(n,m) == i) EXIT
             END DO
             !===There should not be any Neumann forcing on z-axis (1/ray would be infinite) 
             src_phi(i,1) = src_phi(i,1)-phi_mesh%gauss%rjs(ls,ms)*ray*( &
                  +Esolphi_anal(3)*(phi_mesh%gauss%dw_s(2,n,ls,ms)*phi_mesh%gauss%rnorms(1,ls,ms) &
                  -phi_mesh%gauss%dw_s(1,n,ls,ms)*phi_mesh%gauss%rnorms(2,ls,ms))) &
                  -phi_mesh%gauss%rjs(ls,ms)*(&
                  -mode*Esolphi_anal(2)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(2,ls,ms) &
                  +mode*Esolphi_anal(6)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(1,ls,ms))

             src_phi(i,2) = src_phi(i,2)-phi_mesh%gauss%rjs(ls,ms)*ray*( &
                  +Esolphi_anal(4)*(phi_mesh%gauss%dw_s(2,n,ls,ms)*phi_mesh%gauss%rnorms(1,ls,ms) &
                  -phi_mesh%gauss%dw_s(1,n,ls,ms)*phi_mesh%gauss%rnorms(2,ls,ms))) & 
                  -phi_mesh%gauss%rjs(ls,ms)*(&
                  mode*Esolphi_anal(1)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(2,ls,ms) &
                  -mode*Esolphi_anal(5)*phi_mesh%gauss%wws(ns,ls)*phi_mesh%gauss%rnorms(1,ls,ms))
          ENDDO
       ENDDO
    ENDDO
    !===ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION
    !JLG, FL, FEB, 10, 2010
    !We assume that integral int_Gammav n.GRAD (4phi_n-phi_(n-1))/(2dt) psi dsigma = 0
    !JLG, FL, FEB, 10, 2010
    !JLG, FL, May, 28, 2009
    !We assume that integral int_Gammav (Hinfty . n) psi ds = 0!
    !JLG, FL, May, 28, 2009
    !===ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION ATTENTION

    !=========================================================
    !--- Artificial boundary condition: d(phi_t)/dR + (1/R)*phi_t = 0
    !=========================================================

    IF (PRESENT(R_fourier)) THEN
       IF (R_fourier.GE.0.d0) CALL error_Petsc('maxwell_update_time_with_H: R_fourier should be -1')
    END IF
    !IF (.NOT.present(index_fourier) .OR. .NOT.present(R_fourier)) RETURN
    !IF (R_fourier.le.0.d0) RETURN
    !DO ms = 1, phi_mesh%mes
    !   IF (phi_mesh%sides(ms) /= index_fourier) CYCLE ! Not on the artificial boundary

    !   DO ls = 1, phi_mesh%gauss%l_Gs

    !===Compute radius of Gauss point
    !      ray = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms))* phi_mesh%gauss%wws(:,ls))
    !      x = phi_mesh%gauss%rjs(ls,ms)*ray/R_fourier
    !      y1 = x* SUM(rhs(phi_mesh%jjs(:,ms),1)* phi_mesh%gauss%wws(:,ls))
    !      y2 = x* SUM(rhs(phi_mesh%jjs(:,ms),2)* phi_mesh%gauss%wws(:,ls))
    !      DO ns =1, phi_mesh%gauss%n_ws    
    !         src_phi(1,phi_mesh%jjs(ns,ms)) = src_phi(1,phi_mesh%jjs(ns,ms)) + &
    !               y1*phi_mesh%gauss%wws(ns,ls)
    !         src_phi(2,phi_mesh%jjs(ns,ms)) = src_phi(2,phi_mesh%jjs(ns,ms)) + &
    !               y2*phi_mesh%gauss%wws(ns,ls)
    !      ENDDO
    !      
    !   ENDDO
    !END DO
    IF (H_mesh%mes/=0) THEN
!!$       ALLOCATE(idxn(H_mesh%np))
       idxn_H = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,2), ADD_VALUES, ierr)
       idxn_H = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,3), ADD_VALUES, ierr)
       idxn_H = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn_H, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn_H, src_H(:,6), ADD_VALUES, ierr)

       !===JLG Jan 22 2018
       idxn_pmag = LA_pmag%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, pmag_mesh%np, idxn_pmag, src_pmag(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, pmag_mesh%np, idxn_pmag, src_pmag(:,2), ADD_VALUES, ierr)
       !===JLG Jan 22 2018

!!$       DEALLOCATE(idxn)
    END IF
    IF (phi_mesh%mes/=0) THEN
!!$       ALLOCATE(idxn(phi_mesh%np))
       idxn_phi = LA_phi%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, phi_mesh%np, idxn_phi, src_phi(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, phi_mesh%np, idxn_phi, src_phi(:,2), ADD_VALUES, ierr)       
!!$       DEALLOCATE(idxn)
    END IF

    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)

    !DEALLOCATE(src_H,src_phi)

    !===Dummies variables to avoid warning
    count=index_fourier; count=interface_H_mu%mes; count=interface_H_phi%mes
    count=SIZE(list_dirichlet_sides_H)
    !===Dummies variables to avoid warning
  END SUBROUTINE surf_int

  SUBROUTINE  mat_maxwell_mu(H_mesh, jj_v_to_H, interface_H_mu, mode, stab, stab_jump,&
       mu_H_field, sigma, LA_H, H_p_phi_mat1, H_p_phi_mat2, sigma_np)
    USE def_type_mesh
    USE gauss_points
    USE my_util
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc  
    IMPLICIT NONE
    TYPE(mesh_type),            INTENT(IN)    :: H_mesh
    INTEGER,      DIMENSION(:), INTENT(IN)    :: jj_v_to_H
    TYPE(interface_type),       INTENT(IN)    :: interface_H_mu
    INTEGER,                    INTENT(IN)    :: mode
    REAL(KIND=8), DIMENSION(3), INTENT(IN)    :: stab
    REAL(KIND=8),               INTENT(IN)    :: stab_jump
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: sigma, mu_H_field
    REAL(KIND=8), DIMENSION(:), INTENT(IN)    :: sigma_np
    INTEGER :: ms, ls, ni, nj, i, j, &
         n_ws1, n_ws2, n_w2, n_w1, m1, m2, ki, kj,ib,jb, ms1, ms2
    REAL(KIND=8) :: x, y, z, norm, hm1
    REAL(KIND=8) :: ray, stab_colle_H_mu
    LOGICAL :: mark=.FALSE.
    REAL(KIND=8), DIMENSION(9,SIZE(H_mesh%jj,1),SIZE(H_mesh%jj,1),2,2)      :: Hsij, Gsij

    ! MATRICES POUR LES TERMES DE BORDS Hsij et Gsij 
    !=================================================
    ! (--------------------------------------------------)
    ! ( Hsij(1) +G     | GSij(2)        | Hsij(4) +G     )
    ! ( Hsij(1) +G     | GSij(3)        | Hsij(4) +G     )
    ! (--------------------------------------------------)
    ! ( Hsij(2)        | Hsij(5) +G     | Hsij(8)        )
    ! ( Hsij(3)        | Hsij(5) +G     | Hsij(9)        )
    ! (--------------------------------------------------)
    ! ( Hsij(7) +G     |  GSij(8)       | Hsij(6) +G     )
    ! ( Hsij(7) +G     |  GSij(9)       | Hsij(6) +G     )
    ! (==================================================)
!!$ FL+CN 22/03/2013
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs
    REAL(KIND=8), DIMENSION(2, H_mesh%gauss%n_w, H_mesh%gauss%l_Gs, H_mesh%mes) :: dw_cs
    REAL(KIND=8), DIMENSION(2, H_mesh%gauss%n_w) :: dwsi,dwsj
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%l_Gs)  :: gauss1, gauss2
!!$    REAL(KIND=8), DIMENSION(:,:)    , ALLOCATABLE :: w_cs
!!$    REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dw_cs
!!$    REAL(KIND=8), DIMENSION(:,:)    , ALLOCATABLE :: dwsi,dwsj
!!$    REAL(KIND=8), DIMENSION(:,:)    , ALLOCATABLE :: gauss1, gauss2
!!$ FL+CN 22/03/2013
    REAL(KIND=8), DIMENSION(2)                    :: normi, normj
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%jjs,1))   :: wwsi, wwsj 
    INTEGER                                       :: n_wsi, n_wsj, ci, cj, n_wi, n_wj

    INTEGER      :: ls1, ls2
    REAL(KIND=8) :: ref, diff,  mu_H, muhl1, muhl2, muhi, muhj, sigmai, sigmaj
    REAL(KIND=8) :: thetai, thetaj, indice_i, indice_j
    REAL(KIND=8) :: sigmal1, sigmal2!, lambda
    ! June 14 2008
    REAL(KIND=8) :: c_sym =.0d0 ! (c_sym=1.d0 symmetrizes the bilinear form)
    REAL(KIND=8) :: wwiwwj, normt, stab_div
    ! June 14 2008
!!$ FL +CN 22/03/2013
!!$    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE   :: mat_loc1, mat_loc2
!!$    INTEGER     , DIMENSION(:),   ALLOCATABLE   :: idxn, jdxn
    REAL(KIND=8), DIMENSION(6*H_mesh%gauss%n_w,6*H_mesh%gauss%n_w)   :: mat_loc1, mat_loc2
    INTEGER     , DIMENSION(6*H_mesh%gauss%n_w) :: idxn, jdxn
!!$ FL +CN 22/03/2013
    TYPE(petsc_csr_LA)                          :: LA_H
    INTEGER                                     :: ix, jx
    Mat                                         :: H_p_phi_mat1, H_p_phi_mat2
    PetscErrorCode                              :: ierr

    ! June 2009, JLG, CN, normalization
    stab_colle_H_mu = stab(3)
    stab_div = stab(1)
    ! June 2009, JLG, CN

    !**********************************************************************************
    !--------------------TERMS ON SIGMA_MU-------------------------------
    !**********************************************************************************

    !WRITE(*,*) 'Assembling interface_H_mu'
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

!!$    ALLOCATE(w_cs(n_ws1,l_Gs))
!!$    ALLOCATE(dw_cs(2, n_w1, l_Gs, H_mesh%mes))
!!$    ALLOCATE(dwsi(2, n_w1),dwsj(2, n_w2))
!!$    ALLOCATE(gauss1(2,l_Gs),gauss2(2,l_Gs))

!!$    ALLOCATE(mat_loc1(6*n_w1,6*n_w2))
!!$    ALLOCATE(mat_loc2(6*n_w1,6*n_w2))
!!$    ALLOCATE(idxn(6*n_w1))
!!$    ALLOCATE(jdxn(6*n_w2))

    DO ms = 1, interface_H_mu%mes

       ms2 = interface_H_mu%mesh2(ms)
       m2 = H_mesh%neighs(ms2)
       ms1 = interface_H_mu%mesh1(ms)
       m1 = H_mesh%neighs(ms1)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF

       DO ls = 1, l_Gs
          gauss2(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))*H_mesh%gauss%wws(:,ls))
          gauss2(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms2))*H_mesh%gauss%wws(:,ls))
          gauss1(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms1))*H_mesh%gauss%wws(:,ls))
          gauss1(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms1))*H_mesh%gauss%wws(:,ls))
       END DO

       DO ls2 = 1, l_Gs
          ref = SQRT(1.d-8+SUM(gauss2(:,ls2)**2))
          mark = .FALSE.
          DO ls1 = 1, l_Gs
             diff = SQRT(SUM((gauss2(:,ls2)-gauss1(:,ls1))**2))
             IF (diff .LT. 1.d-10) THEN
                dw_cs(:,:,ls2,ms1) =  H_mesh%gauss%dw_s(:,:,ls1,ms1)
                mark = .TRUE.
                EXIT
             END IF
          END DO
          IF (.NOT.mark) WRITE(*,*) ' BUG '
       END DO

    END DO

    DO ms = 1, interface_H_mu%mes

       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       mu_H = SUM(mu_H_field(H_mesh%jj(:,m1)))/H_mesh%gauss%n_w
       !JLG, FL, May, 28, 2009
       !hm1 = stab_colle_H_mu*(((mu_H+mu_H)/mu_H)/SUM(rjs(:,ms2)))
       hm1 = 1/SUM(rjs(:,ms2))
       !JLG, FL, May, 28, 2009
       !====================================================================================
       !------------------------------------TERMES SUR LE BLOC H----------------------------
       !====================================================================================

       !-------------------------------hm1 (bi x ni) . (bj x nj)----------------------------
       !---------------------------------+ (mui bi.ni) (muj bj.nj)--------------------------
       !====================================================================================
       Hsij = 0.d0  
       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = hm1*rjs(ls,ms2)*ray

          !June 14 2008, muhl
          muhl1 = SUM(mu_H_field(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          muhl2 = SUM(mu_H_field(H_mesh%jjs(:,ms2))* wws(:,ls))
          !JLG, FL, May, 28, 2009, Normalization
          normt =stab_colle_H_mu
          normt = stab_colle_H_mu/inputs%sigma_min ! MODIFICATION: normalization for interface H/H term with H x n
          !norm = stab_div*SUM(rjs(:,ms2))**(2*alpha)
          !norm = stab_div*SUM(rjs(:,ms2))**(2*alpha)/MAX(muhl1,muhl2)
          !norm = stab_div*SUM(rjs(:,ms2))**(2*alpha)/MAX(muhl1,muhl2)**2
          !norm = 1.d0/MAX(muhl1,muhl2)
          !norm = 1.d0/MAX(muhl1,muhl2)**2
          !JLG, FL, May, 28, 2009, Normalization
          !June 14 2008, muhl
          norm = stab_div*(SUM(rjs(:,ms2))/H_mesh%global_diameter)**(2*alpha)/(inputs%sigma_min*inputs%mu_min**2) ! MODIFICATION: normalization for divergence stabilization term

          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                muhi = muhl1
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                muhi = muhl2
             END IF
             DO cj = 1, 2 
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   n_wsj = n_ws1
                   muhj = muhl1
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   n_wsj = n_ws2
                   muhj = muhl2
                END IF
                DO ni = 1, n_wsi
                   DO nj = 1, n_wsj
                      wwiwwj = x * wwsi(ni)*wwsj(nj)
                      y = normt * wwiwwj
                      ! June 14 2008, added z
                      z = norm * muhi * muhj * wwiwwj
                      ! June 14 2008, added z
                      Hsij(1,ni,nj,ci,cj) = Hsij(1,ni,nj,ci,cj) + y*normi(2)*normj(2) &
                           + z*normi(1)*normj(1)
                      Hsij(4,ni,nj,ci,cj) = Hsij(4,ni,nj,ci,cj) - y*normj(1)*normi(2) &
                           + z*normi(1)*normj(2)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y*(normi(1)*normj(1) + normi(2)*normj(2)) 
                      Hsij(6,ni,nj,ci,cj) = Hsij(6,ni,nj,ci,cj) + y*normi(1)*normj(1) &
                           + z*normi(2)*normj(2)
                   END DO
                END DO
             END DO
          END DO
       END DO

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ci = 1, 2
          DO ki = 1, 3 
             DO ni = 1, n_ws1 
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = LA_H%loc_to_glob(ki,i)
                ix = ni + n_ws1*((ki-1) + 3*(ci-1))
                idxn(ix) = ib-1

                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_ws2
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = LA_H%loc_to_glob(kj,j)
                         jx = nj + n_ws2*((kj-1) + 3*(cj-1))
                         jdxn(jx) = jb-1
                         IF  ((ki == 1) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Hsij(1,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(1,ni,nj,ci,cj)
                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Hsij(4,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(4,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Hsij(4,nj,ni,cj,ci)
                            mat_loc2(ix,jx) = Hsij(4,nj,ni,cj,ci)
                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                            mat_loc1(ix,jx) = Hsij(5,ni,nj,ci,cj) 
                            mat_loc2(ix,jx) = Hsij(5,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Hsij(6,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(6,ni,nj,ci,cj)
                         ENDIF
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc1(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc2(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)

       !====================================================================================
       !------------------------(1/sigma) (Rot bj) . (bi x ni)------------------------------
       !====================================================================================

       !terme sans derivee
       Hsij = 0.d0
       Gsij = 0.d0
       DO ls = 1, H_mesh%gauss%l_Gs

          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray
          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             sigmal1 = sigma(m1)
          ELSE
             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          END IF
          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
             sigmal2 = sigma(m2)
          ELSE
             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
          END IF
          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                sigmai = sigmal1
                thetai = sigmal2/(sigmal1+sigmal2)
                indice_i = 1.0
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                sigmai = sigmal2
                thetai = sigmal1/(sigmal1+sigmal2)
                indice_i = -1.0
             END IF
             DO cj = 1, 2
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   n_wsj = n_ws1
                   sigmaj = sigmal1
                   thetaj = sigmal1/(sigmal1+sigmal2)
                   indice_j = 1.d0
                   !WRITE(*,*) 'sigma 1', sigmaj
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   n_wsj = n_ws2
                   sigmaj = sigmal2
                   thetaj = sigmal2/(sigmal1+sigmal2)
                   indice_j = -1.d0
                   !WRITE(*,*) 'sigma 2', sigmaj
                END IF

                DO ni = 1,n_wsi  !
                   DO nj = 1, n_wsj!
                      y = x*wwsi(ni)*wwsj(nj)*(thetaj/(sigmaj))!+ thetai*indice_i*indice_j/(sigmaj))
                      Hsij(2,ni,nj,ci,cj) = Hsij(2,ni,nj,ci,cj) + y * (-mode/ray)*normi(1)
                      Hsij(3,ni,nj,ci,cj) = Hsij(3,ni,nj,ci,cj) + y *   mode/ray *normi(1)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y * (-1/ray)   *normi(1)
                      Hsij(8,ni,nj,ci,cj) = Hsij(8,ni,nj,ci,cj) + y * (-mode/ray)*normi(2)
                      Hsij(9,ni,nj,ci,cj) = Hsij(9,ni,nj,ci,cj) + y *   mode/ray *normi(2)
                      !                      y = x*wwsi(ni)*wwsj(nj)/(2*sigmai)
                      !                      Gsij(2,ni,nj,ci,cj) = Gsij(2,ni,nj,ci,cj) + y * (-mode/ray)*normj(1)
                      !                      Gsij(3,ni,nj,ci,cj) = Gsij(3,ni,nj,ci,cj) + y * ( mode/ray)*normj(1)
                      !                      Gsij(5,ni,nj,ci,cj) = Gsij(5,ni,nj,ci,cj) + y * (-1/ray)   *normj(1)
                      !                      Gsij(8,ni,nj,ci,cj) = Gsij(8,ni,nj,ci,cj) + y * (-mode/ray)*normj(2)
                      !                      Gsij(9,ni,nj,ci,cj) = Gsij(9,ni,nj,ci,cj) + y *   mode/ray *normj(2)
                   ENDDO
                ENDDO
             ENDDO
          END DO
       END DO

       !June 14 2008
       Gsij = c_sym*Gsij
       !June 14 2008

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0

       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = LA_H%loc_to_glob(ki,i)
                ix = ni + n_wsi*((ki-1) + 3*(ci-1))
                idxn(ix) = ib - 1
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wsj
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = LA_H%loc_to_glob(kj,j)
                         jx = nj + n_wsj*((kj-1) + 3*(cj-1))
                         jdxn(jx) = jb - 1
                         IF  ((ki == 2) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Hsij(2,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(3,ni,nj,ci,cj)
                         ELSE IF((ki == 1) .AND. (kj == 2)) THEN
                            mat_loc1(ix,jx) = Gsij(2,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(3,ni,nj,ci,cj)
                         ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
                            mat_loc1(ix,jx) = Hsij(5,ni,nj,ci,cj)+Gsij(5,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(5,ni,nj,ci,cj)+Gsij(5,ni,nj,ci,cj)
                         ELSEIF ((ki == 2) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Hsij(8,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(9,ni,nj,ci,cj)
                         ELSEIF ((ki == 3) .AND. (kj == 2)) THEN
                            mat_loc1(ix,jx) = Gsij(8,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(9,ni,nj,ci,cj)
                         ENDIF
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc1(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc2(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)

       !terme avec derivees
       Hsij = 0.d0
       Gsij = 0.d0
       DO ls = 1, H_mesh%gauss%l_Gs

          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
          x = rjs(ls,ms2)*ray
          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             sigmal1 = sigma(m1)
          ELSE
             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          END IF
          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
             sigmal2 = sigma(m2)
          ELSE
             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
          END IF

          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                dwsi = dw_cs(:,:,ls,ms1)
                n_wsi = n_ws1
                n_wi = n_w1
                sigmai = sigmal1
                thetai = sigmal2/(sigmal1+sigmal2)
                indice_i = 1.d0
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                dwsi = dw_s(:,:,ls,ms2)
                n_wsi = n_ws2
                n_wi = n_w2
                sigmai = sigmal2
                thetai = sigmal1/(sigmal1+sigmal2)
                indice_i = -1.d0
             END IF
             DO cj = 1, 2
                IF (cj==1) THEN
                   normj = rnorms(:,ls,ms1)
                   wwsj = w_cs(:,ls)
                   dwsj = dw_cs(:,:,ls,ms1)
                   n_wsj = n_ws1
                   n_wj = n_w1
                   sigmaj = sigmal1                   
                   thetaj = sigmal1/(sigmal1+sigmal2)
                   indice_j = 1.d0
                ELSE
                   normj = rnorms(:,ls,ms2)
                   wwsj = wws(:,ls)
                   dwsj = dw_s(:,:,ls,ms2) 
                   n_wsj = n_ws2
                   n_wj = n_w2
                   sigmaj = sigmal2
                   thetaj = sigmal2/(sigmal1+sigmal2) 
                   indice_j = -1.d0
                END IF

                !termes avec derivees
                DO ni = 1,n_wsi
                   DO nj = 1, n_wj
                      y = x*wwsi(ni)*(thetaj/(sigmaj))
                      Hsij(1,ni,nj,ci,cj) = Hsij(1,ni,nj,ci,cj) + y*(-dwsj(2,nj))*normi(2)
                      Hsij(4,ni,nj,ci,cj) = Hsij(4,ni,nj,ci,cj) + y*  dwsj(1,nj) *normi(2)
                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y*(-dwsj(2,nj) *normi(2)-dwsj(1,nj)*normi(1))
                      Hsij(6,ni,nj,ci,cj) = Hsij(6,ni,nj,ci,cj) + y*(-dwsj(1,nj))*normi(1)
                      Hsij(7,ni,nj,ci,cj) = Hsij(7,ni,nj,ci,cj) + y*  dwsj(2,nj) *normi(1)
                   ENDDO
                END DO
                !                DO ni = 1,n_wi
                !                   DO nj = 1, n_wsj
                !                      y = x*wwsj(nj)/(2*sigmai)
                !                      Gsij(1,ni,nj,ci,cj) = Gsij(1,ni,nj,ci,cj) + y*(-dwsi(2,ni))*normj(2)
                !                      Gsij(4,ni,nj,ci,cj) = Gsij(4,ni,nj,ci,cj) + y*  dwsi(2,ni) *normj(1)
                !                      Gsij(5,ni,nj,ci,cj) = Gsij(5,ni,nj,ci,cj) + y*(-dwsi(2,ni) *normj(2)-dwsi(1,ni)*normj(1))
                !                      Gsij(6,ni,nj,ci,cj) = Gsij(6,ni,nj,ci,cj) + y*(-dwsi(1,ni))*normj(1)
                !                      Gsij(7,ni,nj,ci,cj) = Gsij(7,ni,nj,ci,cj) + y*  dwsi(1,ni) *normj(2)
                !                   ENDDO
                !                END DO

             ENDDO
          ENDDO
       ENDDO

       !June 14 2008
       Gsij = c_sym*Gsij
       !June 14 2008

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                ib = LA_H%loc_to_glob(ki,i)
                ix = ni + n_wsi*((ki-1) + 3*(ci-1))
                idxn(ix) = ib - 1
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wj
                         IF (cj==1) THEN
                            j = H_mesh%jj(nj,m1)
                         ELSE
                            j = H_mesh%jj(nj,m2)
                         END IF
                         jb = LA_H%loc_to_glob(kj,j)
                         jx = nj + n_wj*((kj-1) + 3*(cj-1))
                         jdxn(jx) = jb - 1
                         IF      ((ki == 1) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Hsij(1,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(1,ni,nj,ci,cj)
                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Hsij(4,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(4,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Hsij(7,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(7,ni,nj,ci,cj)
                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN  
                            mat_loc1(ix,jx) = Hsij(5,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Hsij(5,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
                            mat_loc1(ix,jx) = Hsij(6,ni,nj,ci,cj)   
                            mat_loc2(ix,jx) = Hsij(6,ni,nj,ci,cj)
                         ENDIF

                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 6*n_ws1, idxn(1:6*n_ws1), 6*n_w2, jdxn(1:6*n_w2), &
            mat_loc1(1:6*n_ws1,1:6*n_w2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 6*n_ws1, idxn(1:6*n_ws1), 6*n_w2, jdxn(1:6*n_w2), &
            mat_loc2(1:6*n_ws1,1:6*n_w2), ADD_VALUES, ierr)

       mat_loc1 = 0.d0
       mat_loc2 = 0.d0
       DO ci = 1, 2
          DO ki = 1, 3
             DO ni = 1, n_wi
                IF (ci==1) THEN
                   i = H_mesh%jj(ni,m1)
                ELSE
                   i = H_mesh%jj(ni,m2)
                END IF
                ib = LA_H%loc_to_glob(ki,i)
                ix = ni + n_wi*((ki-1) + 3*(ci-1))
                idxn(ix) = ib - 1
                DO cj = 1, 2
                   DO kj = 1, 3
                      DO nj = 1, n_wsj
                         IF (cj==1) THEN
                            j = interface_H_mu%jjs1(nj,ms)
                         ELSE
                            j = interface_H_mu%jjs2(nj,ms)
                         END IF
                         jb = LA_H%loc_to_glob(kj,j)
                         jx = nj + n_wsj*((kj-1) + 3*(cj-1))
                         jdxn(jx) = jb-1 
                         IF      ((ki == 1) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Gsij(1,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(1,ni,nj,ci,cj)
                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Gsij(4,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(4,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
                            mat_loc1(ix,jx) = Gsij(7,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(7,ni,nj,ci,cj)
                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN
                            mat_loc1(ix,jx) = Gsij(5,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(5,ni,nj,ci,cj)
                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN
                            mat_loc1(ix,jx) = Gsij(6,ni,nj,ci,cj)
                            mat_loc2(ix,jx) = Gsij(6,ni,nj,ci,cj)
                         ENDIF
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO

       CALL MatSetValues(H_p_phi_mat1, 6*n_w1, idxn(1:6*n_w1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc1(1:6*n_w1,1:6*n_ws2), ADD_VALUES, ierr)
       CALL MatSetValues(H_p_phi_mat2, 6*n_w1, idxn(1:6*n_w1), 6*n_ws2, jdxn(1:6*n_ws2), &
            mat_loc2(1:6*n_w1,1:6*n_ws2), ADD_VALUES, ierr)

       !    END DO ! end loop on ms

       !    DO ms = 1, interface_H_mu%mes

       !       ms2 = interface_H_mu%mesh2(ms)
       !       ms1 = interface_H_mu%mesh1(ms)
       !       m2 = H_mesh%neighs(ms2)
       !       m1 = H_mesh%neighs(ms1)
       !       hm1 = 1/SUM(rjs(:,ms2))
       !       IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(ms1))) > 0) CYCLE
       !WRITE(*,*) 'H_mesh sides', H_mesh%sides(ms)
       !
       !       !====================================================================================
       !       !---------------------(lambda/h) (sign) (Rot bj) . (bi x ni)-------------------------
       !       !---- (lambda/h) (Curl(H1)/sigma1 - curl(H2)/sigma2) . (b1 x n1 + b2xn2)
       !       !---- where lambda = 2.d0 / (sigma1 + sigma2)
       !       !====================================================================================
       !
       !       !terms without spatial derivatives
       !       Hsij = 0.d0
       !       DO ls = 1, H_mesh%gauss%l_Gs
       !
       !          !===Compute radius of Gauss point
       !          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
       !          x = rjs(ls,ms2)*ray
       !
       !          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
       !             sigmal1 = sigma(m1)
       !          ELSE
       !             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
       !          END IF
       !          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
       !             sigmal2 = sigma(m2)
       !          ELSE
       !             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
       !          END IF   
       !
       !          lambda = 2.d0/(sigmal1+sigmal2)
       !
       !          DO ci = 1, 2
       !             IF (ci==1) THEN
       !                normi = rnorms(:,ls,ms1)
       !                wwsi = w_cs(:,ls)
       !                n_wsi = n_ws1
       !                sigmai = sigmal1
       !                indice_i = 1.d0
       !             ELSE
       !                normi = rnorms(:,ls,ms2)
       !                wwsi = wws(:,ls)
       !                n_wsi = n_ws2
       !                sigmai = sigmal2
       !                indice_i =-1.d0
       !             END IF
       !             DO cj = 1, 2
       !                IF (cj==1) THEN
       !                   normj = rnorms(:,ls,ms1)
       !                   wwsj = w_cs(:,ls)
       !                   n_wsj = n_ws1
       !                   sigmaj = sigmal1
       !                   indice_j = 1.d0
       !                ELSE
       !                   normj = rnorms(:,ls,ms2)
       !                   wwsj = wws(:,ls)
       !                   n_wsj = n_ws2
       !                   sigmaj = sigmal2
       !                   indice_j = -1.d0
       !                END IF
       !
       !                DO ni = 1,n_wsi  !
       !                   DO nj = 1, n_wsj!
       !                      y = stab_jump*lambda*hm1*indice_j*x*wwsi(ni)*wwsj(nj)/(sigmaj)
       !                      Hsij(2,ni,nj,ci,cj) = Hsij(2,ni,nj,ci,cj) + y * (-mode/ray)*normi(1)
       !                      Hsij(3,ni,nj,ci,cj) = Hsij(3,ni,nj,ci,cj) + y *   mode/ray *normi(1)
       !                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y * (-1/ray)   *normi(1)
       !                      Hsij(8,ni,nj,ci,cj) = Hsij(8,ni,nj,ci,cj) + y * (-mode/ray)*normi(2)
       !                      Hsij(9,ni,nj,ci,cj) = Hsij(9,ni,nj,ci,cj) + y *   mode/ray *normi(2)
       !                   ENDDO
       !                ENDDO
       !             ENDDO
       !          END DO
       !       END DO
       !       mat_loc1 = 0.d0
       !       mat_loc2 = 0.d0
       !       DO ci = 1, 2
       !          DO ki = 1, 3
       !             DO ni = 1, n_wsi
       !                IF (ci==1) THEN
       !                   i = interface_H_mu%jjs1(ni,ms)
       !                ELSE
       !                   i = interface_H_mu%jjs2(ni,ms)
       !                END IF
       !                ib = LA_H%loc_to_glob(ki,i)
       !                ix = ni + n_wsi*((ki-1) + 3*(ci-1))
       !                idxn(ix) = ib - 1
       !                DO cj = 1, 2
       !                   DO kj = 1, 3
       !                      DO nj = 1, n_wsj
       !                         IF (cj==1) THEN
       !                            j = interface_H_mu%jjs1(nj,ms)
       !                         ELSE
       !                            j = interface_H_mu%jjs2(nj,ms)
       !                         END IF
       !                         jb = LA_H%loc_to_glob(kj,j)
       !                         jx = nj + n_wsj*((kj-1) + 3*(cj-1))
       !                         jdxn(jx) = jb - 1
       !                         IF  ((ki == 2) .AND. (kj == 1)) THEN
       !                            mat_loc1(ix,jx) = Hsij(2,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(3,ni,nj,ci,cj)
       !                         ELSEIF  ((ki == 2) .AND. (kj == 2))  THEN  
       !                            mat_loc1(ix,jx) = Hsij(5,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(5,ni,nj,ci,cj)
       !                         ELSEIF ((ki == 2) .AND. (kj == 3)) THEN
       !                            mat_loc1(ix,jx) = Hsij(8,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(9,ni,nj,ci,cj)
       !                         ENDIF
       !                      END DO
       !                   END DO
       !                END DO
       !             END DO
       !          END DO
       !       END DO
       !
       !       CALL MatSetValues(H_p_phi_mat1, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
       !            mat_loc1(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)
       !       CALL MatSetValues(H_p_phi_mat2, 6*n_ws1, idxn(1:6*n_ws1), 6*n_ws2, jdxn(1:6*n_ws2), &
       !            mat_loc2(1:6*n_ws1,1:6*n_ws2), ADD_VALUES, ierr)
       !
       !       !terms with spatial derivatives
       !       Hsij = 0.d0
       !       DO ls = 1, H_mesh%gauss%l_Gs
       !
       !          !===Compute radius of Gauss point
       !          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))
       !          x = rjs(ls,ms2)*ray
       !
       !          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
       !             sigmal1 = sigma(m1)
       !          ELSE
       !             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
       !          END IF
       !          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
       !             sigmal2 = sigma(m2)
       !          ELSE
       !             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
       !          END IF
       !
       !          lambda = 2.d0/(sigmal1+sigmal2)
       !
       !          DO ci = 1, 2
       !             IF (ci==1) THEN
       !                normi = rnorms(:,ls,ms1)
       !                wwsi = w_cs(:,ls)
       !                dwsi = dw_cs(:,:,ls,ms1)
       !                n_wsi = n_ws1
       !                n_wi = n_w1
       !                sigmai = sigmal1
       !                indice_i = 1.d0
       !             ELSE
       !                normi = rnorms(:,ls,ms2)
       !                wwsi = wws(:,ls)
       !                dwsi = dw_s(:,:,ls,ms2)
       !                n_wsi = n_ws2
       !                n_wi = n_w2
       !                sigmai = sigmal2
       !                indice_i =-1.d0
       !             END IF
       !             DO cj = 1, 2
       !                IF (cj==1) THEN
       !                   normj = rnorms(:,ls,ms1)
       !                   wwsj = w_cs(:,ls)
       !                   dwsj = dw_cs(:,:,ls,ms1)
       !                   n_wsj = n_ws1
       !                   n_wj = n_w1
       !                   sigmaj = sigmal1
       !                   indice_j = 1.d0
       !                ELSE
       !                   normj = rnorms(:,ls,ms2)
       !                   wwsj = wws(:,ls)
       !                   dwsj = dw_s(:,:,ls,ms2) 
       !                   n_wsj = n_ws2
       !                   n_wj = n_w2
       !                   sigmaj = sigmal2
       !                   indice_j = -1.d0
       !                END IF
       !
       !                !termes avec derivees
       !                DO ni = 1,n_wsi
       !                   DO nj = 1, n_wj
       !                      y = stab_jump*lambda*hm1*indice_j*x*wwsi(ni)/(sigmaj)
       !                      Hsij(1,ni,nj,ci,cj) = Hsij(1,ni,nj,ci,cj) + y*(-dwsj(2,nj))*normi(2)
       !                      Hsij(4,ni,nj,ci,cj) = Hsij(4,ni,nj,ci,cj) + y*  dwsj(1,nj) *normi(2)
       !                      Hsij(5,ni,nj,ci,cj) = Hsij(5,ni,nj,ci,cj) + y*(-dwsj(2,nj) *normi(2)-dwsj(1,nj)*normi(1))
       !                      Hsij(6,ni,nj,ci,cj) = Hsij(6,ni,nj,ci,cj) + y*(-dwsj(1,nj))*normi(1)
       !                      Hsij(7,ni,nj,ci,cj) = Hsij(7,ni,nj,ci,cj) + y*  dwsj(2,nj) *normi(1)
       !                   ENDDO
       !                END DO
       !
       !             ENDDO
       !          ENDDO
       !       ENDDO
       !
       !       mat_loc1 = 0.d0
       !       mat_loc2 = 0.d0
       !       DO ci = 1, 2
       !          DO ki = 1, 3
       !             DO ni = 1, n_wsi
       !                IF (ci==1) THEN
       !                   i = interface_H_mu%jjs1(ni,ms)
       !                ELSE
       !                   i = interface_H_mu%jjs2(ni,ms)
       !                END IF
       !                ib = LA_H%loc_to_glob(ki,i)
       !                ix = ni + n_wsi*((ki-1) + 3*(ci-1))
       !                idxn(ix) = ib - 1
       !                DO cj = 1, 2
       !                   DO kj = 1, 3
       !                      DO nj = 1, n_wj
       !                         IF (cj==1) THEN
       !                            j = H_mesh%jj(nj,m1)
       !                         ELSE
       !                            j = H_mesh%jj(nj,m2)
       !                         END IF
       !                         jb = LA_H%loc_to_glob(kj,j)
       !                         jx = nj + n_wj*((kj-1) + 3*(cj-1))
       !                         jdxn(jx) = jb - 1
       !                         IF      ((ki == 1) .AND. (kj == 1)) THEN
       !                            mat_loc1(ix,jx) = Hsij(1,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(1,ni,nj,ci,cj)
       !                         ELSEIF  ((ki == 1) .AND. (kj == 3)) THEN
       !                            mat_loc1(ix,jx) = Hsij(4,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(4,ni,nj,ci,cj)
       !                         ELSEIF  ((ki == 3) .AND. (kj == 1)) THEN
       !                            mat_loc1(ix,jx) = Hsij(7,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(7,ni,nj,ci,cj)
       !                         ELSEIF  ((ki == 2) .AND. (kj == 2)) THEN  
       !                            mat_loc1(ix,jx) = Hsij(5,ni,nj,ci,cj)
       !                            mat_loc2(ix,jx) = Hsij(5,ni,nj,ci,cj)
       !                         ELSEIF  ((ki == 3) .AND. (kj == 3)) THEN  
       !                            mat_loc1(ix,jx) = Hsij(6,ni,nj,ci,cj)   
       !                            mat_loc2(ix,jx) = Hsij(6,ni,nj,ci,cj)
       !                         ENDIF
       !
       !                      END DO
       !                   END DO
       !                END DO
       !             END DO
       !          END DO
       !       END DO
       !
       !       CALL MatSetValues(H_p_phi_mat1, 6*n_ws1, idxn(1:6*n_ws1), 6*n_w2, jdxn(1:6*n_w2), &
       !            mat_loc1(1:6*n_ws1,1:6*n_w2), ADD_VALUES, ierr)
       !       CALL MatSetValues(H_p_phi_mat2, 6*n_ws1, idxn(1:6*n_ws1), 6*n_w2, jdxn(1:6*n_w2), &
       !            mat_loc2(1:6*n_ws1,1:6*n_w2), ADD_VALUES, ierr)
    END DO

    CALL MatAssemblyBegin(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat1,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyBegin(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(H_p_phi_mat2,MAT_FINAL_ASSEMBLY,ierr)

!!$    DEALLOCATE(mat_loc1, mat_loc2, idxn, jdxn)
!!$    DEALLOCATE(w_cs, dw_cs, gauss1, gauss2, dwsi, dwsj)

  END SUBROUTINE  mat_maxwell_mu

  SUBROUTINE courant_mu(H_mesh,interface_H_mu,sigma,mu_H_field,time,mode,nl, &
       LA_H, vb_1, vb_2, B_ext, J_over_sigma_gauss, sigma_curl_gauss)
    !forcage faisant intervenir J, volumique et interface pour H et phi
    !pour le probleme en entier
    USE def_type_mesh
    USE gauss_points
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_mu
    REAL(KIND=8),                          INTENT(IN)   :: time
    REAL(KIND=8), DIMENSION(H_mesh%me),    INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: mu_H_field
    INTEGER,                               INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: nl
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: B_ext
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: J_over_sigma_gauss
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: sigma_curl_gauss
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                :: src_H

    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: normi, gaussp1, gaussp2
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws) :: wwsi 
    REAL(KIND=8) ::  x, ray
    INTEGER :: i, ni, ms, k, ls, n_ws1, n_ws2, ms1, ms2, n_w1, n_w2, m1, m2, ci, n_wsi
    INTEGER :: mesh_id1, mesh_id2, index
    REAL(KIND=8), DIMENSION(6)             :: JsolH_anal, test, B_ext_l, J_over_sigma_l
    REAL(KIND=8) :: muhl1, muhl2, ref, diff
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one
    DATA one/1.d0/
    !April 17th, 2008, JLG
    !$$ FL+CN 22/03/2013
!!$    INTEGER, DIMENSION(:), ALLOCATABLE          :: idxn
    INTEGER, DIMENSION(H_mesh%np)               :: idxn
    !$$ FL+CN 22/03/2013
    TYPE(petsc_csr_LA)                          :: LA_H
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2


    !**********************************************************************************
    !--------------------TERMS ON SIGMA_MU-------------------------------
    !**********************************************************************************
    src_H = 0.d0
    !WRITE(*,*) 'Assembling rhs interface_H_mu'
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

    DO ms = 1, interface_H_mu%mes
       ms1 = interface_H_mu%mesh1(ms)
       ms2 = interface_H_mu%mesh2(ms)
       m1 = H_mesh%neighs(ms1)
       m2 = H_mesh%neighs(ms2)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF
    END DO

    index=0
    DO ms = 1, interface_H_mu%mes
       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       mesh_id1 = H_mesh%i_d(m1)
       mesh_id2 = H_mesh%i_d(m2)
       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))


          ! Side 1
          index=index+1
          DO k=1, 6
             B_ext_l(k) = SUM(B_ext(H_mesh%jjs(:,ms1),k)*H_mesh%gauss%wws(:,ls))
             J_over_sigma_l(k) = J_over_sigma_gauss(index,k) + sigma_curl_gauss(index,k)
          END DO
          muhl1=SUM(mu_H_field(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          gaussp1(1) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms1))*w_cs(:,ls))    
          gaussp1(2) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms1))*w_cs(:,ls))

          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
             DO k=1, 6
                JsolH_anal(k) = J_over_sigma_l(k) &
                     + muhl1 *SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
             ENDDO
          ELSE
             DO k=1, 6
                JsolH_anal(k) = Jexact_gauss(k, gaussp1, mode, one ,sigma(m1), &
                     muhl1, time, mesh_id1, B_ext_l)/sigma(m1) &
                     + muhl1 * SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
             ENDDO
          END IF

          ! Side 2
          index=index+1
          DO k=1, 6
             B_ext_l(k) = SUM(B_ext(H_mesh%jjs(:,ms2),k)*H_mesh%gauss%wws(:,ls))
             J_over_sigma_l(k) = J_over_sigma_gauss(index,k) + sigma_curl_gauss(index,k)
          END DO
          muhl2=SUM(mu_H_field(H_mesh%jjs(:,ms2))*wws(:,ls))
          gaussp2(1) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))*wws(:,ls))    
          gaussp2(2) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms2))*wws(:,ls))    
          IF (MAXVAL(ABS(gaussp1-gaussp2)) > 1.d-11) THEN
             WRITE(*,*) ' BUG courant_mu '
             STOP
          END IF
          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
             DO k=1, 6
                test(k) = J_over_sigma_l(k) &
                     + muhl2 * SUM(NL(H_mesh%jjs(1:n_ws2,ms2),k)*wws(1:n_ws2,ls))
                JsolH_anal(k) = JsolH_anal(k) + test(k)
             ENDDO
          ELSE
             DO k=1, 6
                test(k) = Jexact_gauss(k, gaussp2, mode, one ,sigma(m2), &
                     muhl2, time, mesh_id2, B_ext_l)/sigma(m2) &
                     + muhl2 * SUM(NL(H_mesh%jjs(1:n_ws2,ms2),k)*wws(1:n_ws2,ls))
                JsolH_anal(k) = JsolH_anal(k) + test(k)
             ENDDO
          END IF
          ! Division by 2 to get the mean is in definition of x below.

          !---------forcage pour H            
          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
             END IF
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                x = rjs(ls,ms2)*ray*wwsi(ni)/2
                src_H(i,1) = src_H(i,1)+x*(-JsolH_anal(3)*normi(2))
                src_H(i,2) = src_H(i,2)+x*(-JsolH_anal(4)*normi(2))
                src_H(i,3) = src_H(i,3)+x*(JsolH_anal(1)*normi(2)-JsolH_anal(5)*normi(1))
                src_H(i,4) = src_H(i,4)+x*(JsolH_anal(2)*normi(2)-JsolH_anal(6)*normi(1))
                src_H(i,5) = src_H(i,5)+x*(JsolH_anal(3)*normi(1))
                src_H(i,6) = src_H(i,6)+x*(JsolH_anal(4)*normi(1))
             END DO
          ENDDO
       END DO
    END DO

    IF (H_mesh%np /= 0) THEN
!!$       ALLOCATE(idxn(H_mesh%np))
       idxn = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,2), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,3), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,6), ADD_VALUES, ierr)
!!$       DEALLOCATE(idxn)
    END IF
    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)

  END SUBROUTINE courant_mu

  SUBROUTINE jump_rot_H_consistant(H_mesh, jj_v_to_H, interface_H_mu, stab_jump, sigma, NL, &
       LA_H, vb_1, vb_2, sigma_np)
    USE def_type_mesh
    USE gauss_points
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    INTEGER,      DIMENSION(:), INTENT(IN)              :: jj_v_to_H
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_mu
    REAL(KIND=8),                          INTENT(IN)   :: stab_jump
    REAL(KIND=8), DIMENSION(:),    INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: NL
    REAL(KIND=8), DIMENSION(:), INTENT(IN)              :: sigma_np
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                :: src_H

    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: normi
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws) :: wwsi 
    REAL(KIND=8) ::  x, ray, hm1,y
    INTEGER :: i, ni, ms, k, ls, n_ws1, n_ws2, ms1, ms2, n_w1, n_w2, m1, m2, ci, n_wsi, indice_i
    INTEGER :: mesh_id1, mesh_id2
    REAL(KIND=8), DIMENSION(6)             :: jump_sol 
    REAL(KIND=8) :: ref, diff
    REAL(KIND=8) :: sigmal1, sigmal2, thetai
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one
    DATA one/1.d0/
    !April 17th, 2008, JLG
    INTEGER, DIMENSION(H_mesh%np)               :: idxn
    TYPE(petsc_csr_LA)                          :: LA_H
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2
    src_H = 0.d0
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

    DO ms = 1, interface_H_mu%mes
       ms1 = interface_H_mu%mesh1(ms)
       ms2 = interface_H_mu%mesh2(ms)
       m1 = H_mesh%neighs(ms1)
       m2 = H_mesh%neighs(ms2)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF
    END DO

    DO ms = 1, interface_H_mu%mes
       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       mesh_id1 = H_mesh%i_d(m1)
       mesh_id2 = H_mesh%i_d(m2)
       hm1 = 1/SUM(rjs(:,ms2))


       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))

          x = rjs(ls,ms2)*ray

          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             sigmal1 = sigma(m1)
          ELSE
             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          END IF
          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
             sigmal2 = sigma(m2)
          ELSE
             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
          END IF

          DO k=1, 6
             jump_sol(k) =+SUM(NL(H_mesh%jjs(1:n_ws1,ms1),k)*w_cs(1:n_ws1,ls))
          END DO

          !---------forcage pour H            
          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                indice_i=1.d0
                thetai=sigmal2/(sigmal1+sigmal2)
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                indice_i=-1.d0
                thetai=sigmal1/(sigmal1+sigmal2)
             END IF
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                y = -thetai*indice_i*x*wwsi(ni)
                src_H(i,1) = src_H(i,1)+y*(-jump_sol(3)*normi(2))
                src_H(i,2) = src_H(i,2)+y*(-jump_sol(4)*normi(2))
                src_H(i,3) = src_H(i,3)+y*(jump_sol(1)*normi(2)-jump_sol(5)*normi(1))
                src_H(i,4) = src_H(i,4)+y*(jump_sol(2)*normi(2)-jump_sol(6)*normi(1))
                src_H(i,5) = src_H(i,5)+y*(jump_sol(3)*normi(1))
                src_H(i,6) = src_H(i,6)+y*(jump_sol(4)*normi(1))
             END DO
          ENDDO
       END DO
    END DO
    IF (H_mesh%np /= 0) THEN
       idxn = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,2), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,3), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,6), ADD_VALUES, ierr)
    END IF
    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)

  END SUBROUTINE jump_rot_H_consistant

  SUBROUTINE jump_rot_H_consistant_rhoLi(H_mesh, jj_v_to_H, interface_H_mu, stab_jump, sigma, &
       function_of_rhoLi, LA_H, vb_1, vb_2, sigma_np)
    USE def_type_mesh
    USE gauss_points
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    INTEGER,      DIMENSION(:),            INTENT(IN)   :: jj_v_to_H
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_mu
    REAL(KIND=8),                          INTENT(IN)   :: stab_jump
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: function_of_rhoLi 
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: sigma_np
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                :: src_H
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_cs   
    REAL(KIND=8), DIMENSION(2) :: normi
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws) :: wwsi 
    REAL(KIND=8) ::  x, ray, hm1,y
    INTEGER :: i, ni, ms, k, ls, n_ws1, n_ws2, ms1, ms2, n_w1, n_w2, m1, m2, ci, n_wsi, indice_i, index
    INTEGER :: mesh_id1, mesh_id2
    REAL(KIND=8), DIMENSION(6)             :: jump_sol 
    REAL(KIND=8) :: ref, diff
    REAL(KIND=8) :: sigmal1, sigmal2, thetai
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one
    DATA one/1.d0/
    !April 17th, 2008, JLG
    INTEGER, DIMENSION(H_mesh%np)               :: idxn
    TYPE(petsc_csr_LA)                          :: LA_H
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2

    src_H = 0.d0
    CALL gauss(H_mesh)
    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws
    n_w1  = H_mesh%gauss%n_w
    n_w2  = H_mesh%gauss%n_w

    DO ms = 1, interface_H_mu%mes
       !     IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(ms))) > 0) CYCLE
       ms1 = interface_H_mu%mesh1(ms)
       ms2 = interface_H_mu%mesh2(ms)
       m1 = H_mesh%neighs(ms1)
       m2 = H_mesh%neighs(ms2)

       IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(ms1))) > 0) CYCLE
       !WRITE(*,*) 'sides ms', H_mesh%sides(ms)
       !WRITE(*,*) 'sides ms1', H_mesh%sides(ms1)

       ref = 1.d-8+SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(2,ms1)))**2)
       diff =      SUM((H_mesh%rr(:,H_mesh%jjs(1,ms1)) - H_mesh%rr(:,H_mesh%jjs(1,ms2)))**2)
       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
          w_cs = wws
       ELSE                ! 1 = 2
          DO ls = 1, l_Gs
             w_cs(1,ls)= wws(2,ls)
             w_cs(2,ls)= wws(1,ls)
             IF (n_ws1==3) w_cs(n_ws1,ls) = wws(n_ws1,ls) 
             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
          END DO
       END IF
    END DO

    index = 0
    DO ms = 1, interface_H_mu%mes
       !     IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(ms))) > 0) CYCLE
       ms2 = interface_H_mu%mesh2(ms)
       ms1 = interface_H_mu%mesh1(ms)
       m2 = H_mesh%neighs(ms2)
       m1 = H_mesh%neighs(ms1)
       mesh_id1 = H_mesh%i_d(m1)
       mesh_id2 = H_mesh%i_d(m2)
       hm1 = 1/SUM(rjs(:,ms2))
       IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(ms1))) > 0) CYCLE

       DO ls = 1, l_Gs
          index = index + 1
          !===Compute radius of Gauss point
          ray = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* H_mesh%gauss%wws(:,ls))

          x = rjs(ls,ms2)*ray

          IF (jj_v_to_H(H_mesh%jj(1,m1)) == -1) THEN
             sigmal1 = sigma(m1)
          ELSE
             sigmal1 = SUM(sigma_np(H_mesh%jjs(:,ms1))*w_cs(:,ls))
          END IF
          IF (jj_v_to_H(H_mesh%jj(1,m2)) == -1) THEN
             sigmal2 = sigma(m2)
          ELSE
             sigmal2 = SUM(sigma_np(H_mesh%jjs(:,ms2))* wws(:,ls))
          END IF

          DO k = 1, 6
             jump_sol(k) = function_of_rhoLi(index,k)
          END DO

          !---------forcage pour H            
          DO ci = 1, 2
             IF (ci==1) THEN
                normi = rnorms(:,ls,ms1)
                wwsi = w_cs(:,ls)
                n_wsi = n_ws1
                indice_i=1.d0
                thetai=sigmal2/(sigmal1+sigmal2)
             ELSE
                normi = rnorms(:,ls,ms2)
                wwsi = wws(:,ls)
                n_wsi = n_ws2
                indice_i=-1.d0
                thetai=sigmal1/(sigmal1+sigmal2)
             END IF
             DO ni = 1, n_wsi
                IF (ci==1) THEN
                   i = interface_H_mu%jjs1(ni,ms)
                ELSE
                   i = interface_H_mu%jjs2(ni,ms)
                END IF
                y = -thetai*indice_i*x*wwsi(ni)
! check LC-CN: with 1,2,5,6, 3D runs OK; only 3,4 KO
                src_H(i,1) = src_H(i,1)+y*(-jump_sol(3)*normi(2))
                src_H(i,2) = src_H(i,2)+y*(-jump_sol(4)*normi(2))
                src_H(i,3) = src_H(i,3)+y*(jump_sol(1)*normi(2)-jump_sol(5)*normi(1))
                src_H(i,4) = src_H(i,4)+y*(jump_sol(2)*normi(2)-jump_sol(6)*normi(1))
                src_H(i,5) = src_H(i,5)+y*(jump_sol(3)*normi(1))
                src_H(i,6) = src_H(i,6)+y*(jump_sol(4)*normi(1))
             END DO
          ENDDO
       END DO
    END DO

    IF (H_mesh%np /= 0) THEN
       idxn = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,2), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,3), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,6), ADD_VALUES, ierr)
    END IF
    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)
  END SUBROUTINE jump_rot_H_consistant_rhoLi

  SUBROUTINE compute_minus_grad_times_der_pot_rhoLi(communicator,H_mesh,interface_H_mu,list_mode,rhoLi,&
       der_pot_rhoLi_node,minus_grad_times_der_pot_rhoLi)
    USE def_type_mesh
    USE gauss_points
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    USE sft_parallele
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    TYPE(interface_type),                  INTENT(IN)   :: interface_H_mu
    INTEGER,      DIMENSION(:),            INTENT(IN)   :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),        INTENT(IN)   :: rhoLi, der_pot_rhoLi_node 
    REAL(KIND=8), DIMENSION(:,:,:),        INTENT(OUT)  :: minus_grad_times_der_pot_rhoLi

    INTEGER,      DIMENSION(H_mesh%gauss%n_w)           :: j_loc
    INTEGER,      DIMENSION(H_mesh%gauss%n_ws)          :: js_loc
    !REAL(KIND=8), DIMENSION(2, H_mesh%gauss%n_ws)       :: dws_loc
    REAL(KIND=8), DIMENSION(2, H_mesh%gauss%n_w)       :: dws_loc
    REAL(KIND=8) ::  ray
    INTEGER :: i, ms, k, ls, index, mode, m, ms1, ms2, m1, m2, mms, mesh_id1, mesh_id2
    REAL(KIND=8), DIMENSION(H_mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode))  :: minus_grad_rhoLi 
    REAL(KIND=8), DIMENSION(H_mesh%gauss%l_Gs*interface_H_mu%mes,2,SIZE(list_mode))  :: der_pot_rhoLi_gauss 
    REAL(KIND=8), DIMENSION(3)                  :: temps
    INTEGER                                     :: code, m_max_pad, bloc_size, nb_procs
    MPI_Comm       :: communicator

    CALL gauss(H_mesh)

    minus_grad_rhoLi=0.d0
    der_pot_rhoLi_gauss=0.d0
    minus_grad_times_der_pot_rhoLi=0.d0

    DO i = 1, SIZE(list_mode)
       mode = list_mode(i)
       index = 0

       DO ms = 1, interface_H_mu%mes
          ms2 = interface_H_mu%mesh2(ms)
          ms1 = interface_H_mu%mesh1(ms)
          m2 = H_mesh%neighs(ms2)
          m1 = H_mesh%neighs(ms1)
          mesh_id1 = H_mesh%i_d(m1)
          mesh_id2 = H_mesh%i_d(m2)

          IF (MINVAL(ABS(mesh_id1 - inputs%list_dom_conc)) == 0) THEN
             m = m1
             mms = ms1
          ELSE IF (MINVAL(ABS(mesh_id2 - inputs%list_dom_conc)) == 0) THEN
             m = m2
             mms = ms2
          ELSE
             CYCLE
          END IF

          IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - H_mesh%sides(mms))) > 0) CYCLE

          j_loc = H_mesh%jj(:,m)
          js_loc = H_mesh%jjs(:,m)
          DO ls = 1, l_Gs
             index = index + 1
             !===Compute radius of Gauss point
             ray = SUM(H_mesh%rr(1,js_loc)* H_mesh%gauss%wws(:,ls))
             dws_loc = dw_s(:,:,ls,mms)
             !--------- Compute grad(rhoLi) in fourier space
             minus_grad_rhoLi(index,1,i) = -SUM(rhoLi(j_loc,1,i)*dws_loc(1,:))
             minus_grad_rhoLi(index,2,i) = -SUM(rhoLi(j_loc,2,i)*dws_loc(1,:))
             minus_grad_rhoLi(index,3,i) = -SUM(rhoLi(js_loc,2,i)*wws(:,ls))*mode/ray
             minus_grad_rhoLi(index,4,i) =  SUM(rhoLi(js_loc,1,i)*wws(:,ls))*mode/ray
             minus_grad_rhoLi(index,5,i) = -SUM(rhoLi(j_loc,1,i)*dws_loc(2,:))
             minus_grad_rhoLi(index,6,i) = -SUM(rhoLi(j_loc,2,i)*dws_loc(2,:))
             !--------- Compute der_pot(rhoLi) on gauss points
             DO k = 1,2
                der_pot_rhoLi_gauss(index,k,i) = SUM(der_pot_rhoLi_node(js_loc,k,i)*wws(:,ls))
             ENDDO
          ENDDO
       END DO
    END DO

    temps = 0
    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
    bloc_size = SIZE(minus_grad_rhoLi,1)/nb_procs+1
    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
    CALL FFT_SCALAR_VECT_DCL(communicator, minus_grad_rhoLi, der_pot_rhoLi_gauss,minus_grad_times_der_pot_rhoLi,&
         1, nb_procs,bloc_size, m_max_pad, temps) 

  END SUBROUTINE compute_minus_grad_times_der_pot_rhoLi

  SUBROUTINE rhs_dirichlet(H_mesh, Dirichlet_bdy_H_sides, sigma,&
       mu_H_field, time, mode, nl, stab, LA_H, vb_1, vb_2, B_ext, J_over_sigma_bdy, sigma_curl_bdy)
    !forcage faisant intervenir J, volumique et surfacique
    !pour le probleme en entier
    USE def_type_mesh
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: H_mesh
    INTEGER,      DIMENSION(:),            INTENT(IN)   :: Dirichlet_bdy_H_sides
    REAL(KIND=8),                          INTENT(IN)   :: time
    REAL(KIND=8), DIMENSION(H_mesh%me),    INTENT(IN)   :: sigma
    REAL(KIND=8), DIMENSION(:),            INTENT(IN)   :: mu_H_field
    INTEGER,                               INTENT(IN)   :: mode
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: nl
    REAL(KIND=8), DIMENSION(3),            INTENT(IN)   :: stab
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: B_ext
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: J_over_sigma_bdy !Used only if sigma variable in fluid
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: sigma_curl_bdy !Used only if sigma variable in fluid
    REAL(KIND=8), DIMENSION(H_mesh%np,6)                :: src_H
    REAL(KIND=8), DIMENSION(2) :: gaussp1
    REAL(KIND=8) ::   x, ray, stab_colle_H_mu
    INTEGER :: i, ni, ms, k, ls, m1, count
    INTEGER :: mesh_id1
    REAL(KIND=8), DIMENSION(6)                   :: JsolH_anal, B_ext_l
    REAL(KIND=8) :: muhl1, hm1
    REAL(KIND=8), DIMENSION(6,H_mesh%gauss%l_Gs) :: Hloc, Hlocxn 
    REAL(KIND=8), DIMENSION(2,H_mesh%gauss%l_Gs) :: rloc
    REAL(KIND=8), DIMENSION(1)                   :: muloc
    INTEGER                                      :: index
    !April 17th, 2008, JLG
    REAL(KIND=8) :: one
    DATA one/1.d0/
    !April 17th, 2008, JLG
!!$ FL+CN 22/03/2013
!!$    INTEGER, DIMENSION(:), ALLOCATABLE           :: idxn
    INTEGER, DIMENSION(H_mesh%np)                :: idxn
!!$ FL+CN 22/03/2013
    TYPE(petsc_csr_LA)                           :: LA_H
    PetscErrorCode          :: ierr
    Vec                     :: vb_1, vb_2

    !IF (SIZE(Dirichlet_bdy_H_sides)==0) THEN
    !   RETURN
    !END IF

    src_H = 0.d0

    !IF (SIZE(Dirichlet_bdy_H_sides)==0) THEN
    !   IF (ASSOCIATED(nl)) DEALLOCATE(nl)
    !   RETURN
    !END IF

    stab_colle_H_mu = stab(3)
    index = 0

    DO count = 1, SIZE(Dirichlet_bdy_H_sides)
       ms = Dirichlet_bdy_H_sides(count)
       !hm1 = stab_colle_H_mu/SUM(H_mesh%gauss%rjs(:,ms))     
       hm1 = stab_colle_H_mu/(SUM(H_mesh%gauss%rjs(:,ms))*inputs%sigma_min) ! MODIFICATION: normalization for dirichlet term RHS
       m1 = H_mesh%neighs(ms)
       mesh_id1 = H_mesh%i_d(m1)
       muloc(1) = mu_H_field(H_mesh%jj(1,m1))

       DO ls = 1, H_mesh%gauss%l_Gs
          rloc(1,ls) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls)) 
          rloc(2,ls) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls)) 
       END DO

       DO k = 1, 6
          Hloc(k,:) = Hexact(H_mesh, k, rloc, mode, muloc, time)
       END DO

       Hlocxn(1,:) = Hloc(3,:)*H_mesh%gauss%rnorms(2,:,ms)
       Hlocxn(2,:) = Hloc(4,:)*H_mesh%gauss%rnorms(2,:,ms)
       Hlocxn(3,:) = Hloc(5,:)*H_mesh%gauss%rnorms(1,:,ms)-Hloc(1,:)*H_mesh%gauss%rnorms(2,:,ms)
       Hlocxn(4,:) = Hloc(6,:)*H_mesh%gauss%rnorms(1,:,ms)-Hloc(2,:)*H_mesh%gauss%rnorms(2,:,ms)
       Hlocxn(5,:) = -Hloc(3,:)*H_mesh%gauss%rnorms(1,:,ms)
       Hlocxn(6,:) = -Hloc(4,:)*H_mesh%gauss%rnorms(1,:,ms)

       DO ls = 1, H_mesh%gauss%l_Gs
          index = index + 1
          !===Compute radius of Gauss point
          ray = rloc(1,ls) !SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))* H_mesh%gauss%wws(:,ls))
          DO k = 1, 6
             B_ext_l(k) = SUM(B_ext(H_mesh%jjs(:,ms),k)*H_mesh%gauss%wws(:,ls))
          END DO

          ! Side 1
          muhl1=SUM(mu_H_field(H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))
          gaussp1(1) = rloc(1,ls) !SUM(H_mesh%rr(1,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))    
          gaussp1(2) = rloc(2,ls) !SUM(H_mesh%rr(2,H_mesh%jjs(:,ms))*H_mesh%gauss%wws(:,ls))    
!!$          DO k=1, 6
!!$             JsolH_anal(k) = Jexact_gauss(k, gaussp1, mode, one ,sigma(m1), muhl1, &
!!$                  time, mesh_id1, B_ext_l)/sigma(m1) &
!!$                  + muhl1 * SUM(NL(H_mesh%jjs(:,ms),k)*H_mesh%gauss%wws(:,ls)) &
!!$                  + hm1*Hlocxn(k,ls)
!!$          ENDDO
          IF (inputs%if_level_set.AND.inputs%variation_sigma_fluid) THEN
             DO k = 1, 6
                !SB-CN-LC 2022/01/25
!!$                JsolH_anal(k) = J_over_sigma_bdy(index,k) &
!!$                     + sigma_curl_bdy(index,k) &
!!$                     + muhl1 * SUM(NL(H_mesh%jjs(:,ms),k)*H_mesh%gauss%wws(:,ls)) &
!!$                     + hm1*Hlocxn(k,ls)
                JsolH_anal(k) = hm1*Hlocxn(k,ls)
                !SB-CN-LC 2022/01/25
             END DO
          ELSE
             DO k = 1, 6
                !SB-CN-LC 2022/01/25
!!$                JsolH_anal(k) = Jexact_gauss(k, gaussp1, mode, one ,sigma(m1), muhl1, &
!!$                     time,  mesh_id1, B_ext_l)/sigma(m1) &
!!$                     + muhl1 * SUM(NL(H_mesh%jjs(:,ms),k)*H_mesh%gauss%wws(:,ls)) &
!!$                     + hm1*Hlocxn(k,ls)
                JsolH_anal(k) = Jexact_gauss(k, gaussp1, mode, one ,sigma(m1), muhl1, &
                     time,  mesh_id1, B_ext_l)/sigma(m1) &
                     + hm1*Hlocxn(k,ls)
                !SB-CN-LC 2022/01/25
             END DO
          END IF

          DO ni = 1, H_mesh%gauss%n_ws
             i = H_mesh%jjs(ni,ms)
             x = H_mesh%gauss%rjs(ls,ms)*ray*H_mesh%gauss%wws(ni,ls)

             src_H(i,1) = src_H(i,1)+x*(-JsolH_anal(3)*H_mesh%gauss%rnorms(2,ls,ms))
             src_H(i,2) = src_H(i,2)+x*(-JsolH_anal(4)*H_mesh%gauss%rnorms(2,ls,ms))
             src_H(i,3) = src_H(i,3)+x*(JsolH_anal(1)*H_mesh%gauss%rnorms(2,ls,ms)&
                  -JsolH_anal(5)*H_mesh%gauss%rnorms(1,ls,ms))
             src_H(i,4) = src_H(i,4)+x*(JsolH_anal(2)*H_mesh%gauss%rnorms(2,ls,ms)&
                  -JsolH_anal(6)*H_mesh%gauss%rnorms(1,ls,ms))
             src_H(i,5) = src_H(i,5)+x*(JsolH_anal(3)*H_mesh%gauss%rnorms(1,ls,ms))
             src_H(i,6) = src_H(i,6)+x*(JsolH_anal(4)*H_mesh%gauss%rnorms(1,ls,ms))

          END DO
       ENDDO

    END DO

    IF (H_mesh%np /= 0) THEN
!!$       ALLOCATE(idxn(H_mesh%np))
       idxn = LA_H%loc_to_glob(1,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,2), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(2,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,4), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,3), ADD_VALUES, ierr)
       idxn = LA_H%loc_to_glob(3,:)-1
       CALL VecSetValues(vb_1, H_mesh%np, idxn, src_H(:,5), ADD_VALUES, ierr)
       CALL VecSetValues(vb_2, H_mesh%np, idxn, src_H(:,6), ADD_VALUES, ierr)
!!$       DEALLOCATE(idxn)
    END IF
    CALL VecAssemblyBegin(vb_1,ierr)
    CALL VecAssemblyEnd(vb_1,ierr)
    CALL VecAssemblyBegin(vb_2,ierr)
    CALL VecAssemblyEnd(vb_2,ierr)

  END SUBROUTINE rhs_dirichlet

!!$  SUBROUTINE dirichlet_cavities(communicator, interface_H_phi, mesh, js_D)
!!$    USE def_type_mesh
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE my_util
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc  
!!$    IMPLICIT NONE
!!$    TYPE(interface_type),           INTENT(IN)    :: interface_H_phi
!!$    TYPE(mesh_type),                INTENT(IN)    :: mesh
!!$    INTEGER, POINTER, DIMENSION(:)                :: js_D
!!$    INTEGER, ALLOCATABLE, DIMENSION(:)            :: on_proc_loc, on_proc, not_cav_loc, not_cav
!!$    INTEGER, ALLOCATABLE, DIMENSION(:)            :: is_ok, j_tmp
!!$    INTEGER, DIMENSION(1)                         :: loc
!!$    INTEGER                                       :: m, ms, i, nb_dom, idx, nb_cav, ni
!!$    LOGICAL                                       :: okay
!!$    MPI_Comm,                       INTENT(IN)    :: communicator
!!$    PetscInt                                      :: rank
!!$    PetscErrorCode                                :: ierr
!!$
!!$    IF (inputs%nb_dom_phi==0) RETURN
!!$
!!$    CALL MPI_COMM_RANK(communicator, rank, ierr)
!!$
!!$    nb_dom = inputs%nb_dom_phi 
!!$    ALLOCATE(on_proc_loc(nb_dom), on_proc(nb_dom))
!!$    ALLOCATE(not_cav_loc(nb_dom), not_cav(nb_dom))
!!$    on_proc_loc = -1
!!$    on_proc = -1
!!$    not_cav_loc = -1
!!$    not_cav = -1
!!$
!!$    DO m = 1, mesh%me
!!$       i = mesh%i_d(m)
!!$       IF (MINVAL(ABS(inputs%list_dom_phi-i)) /= 0) THEN
!!$          WRITE(*,*) 'error in dirichlet cavities'
!!$       END IF
!!$       loc = MINLOC(ABS(inputs%list_dom_phi-i))
!!$       on_proc_loc(loc(1)) = rank
!!$    END DO
!!$    IF (mesh%mes /= 0) THEN
!!$       ALLOCATE(is_ok(mesh%mes))
!!$       is_ok = mesh%i_d(mesh%neighs)
!!$       IF (interface_H_phi%mes /=0) THEN
!!$          is_ok(interface_H_phi%mesh2) = 0
!!$       END IF
!!$       DO ms = 1, mesh%mes
!!$          IF (SUM(ABS(mesh%rr(1,mesh%jjs(:,ms)))) .LT. 1.d-12*mesh%global_diameter) THEN
!!$             is_ok(ms) = 0
!!$          END IF
!!$          IF (inputs%my_periodic%nb_periodic_pairs /=0) THEN
!!$             IF (MINVAL(ABS(inputs%my_periodic%list_periodic-mesh%sides(ms))) == 0) THEN
!!$                is_ok(ms) = 0
!!$             END IF
!!$          END IF
!!$       END DO
!!$
!!$       DO ms = 1, mesh%mes
!!$          IF (is_ok(ms) == 0) CYCLE
!!$          i = is_ok(ms)
!!$          IF (MINVAL(ABS(inputs%list_dom_phi-i)) /= 0) THEN
!!$             WRITE(*,*) 'error in dirichlet cavities'
!!$          END IF
!!$          loc = MINLOC(ABS(inputs%list_dom_phi-i))
!!$          not_cav_loc(loc(1)) = rank
!!$       END DO
!!$    END IF
!!$    CALL MPI_ALLREDUCE(on_proc_loc, on_proc, nb_dom, MPI_INTEGER, MPI_MAX, communicator, ierr)
!!$    CALL MPI_ALLREDUCE(not_cav_loc, not_cav, nb_dom, MPI_INTEGER, MPI_MAX, communicator, ierr)
!!$
!!$    ALLOCATE(j_tmp(SIZE(js_D)+nb_dom))
!!$    j_tmp(1:SIZE(js_D)) = js_D
!!$    idx = SIZE(js_D)
!!$    DO i = 1, nb_dom
!!$       IF ( (not_cav(i)==-1) .AND. (on_proc(i)==rank) ) THEN
!!$          idx = idx + 1
!!$          okay = .FALSE.
!!$          DO m = 1, mesh%me
!!$             IF (mesh%i_d(m) == inputs%list_dom_phi(i)) THEN
!!$                DO ni = 1, mesh%gauss%n_w
!!$                   IF (MINVAL(ABS(mesh%jjs-mesh%jj(ni,m))) /=0) THEN
!!$                      j_tmp(idx) = mesh%jj(ni,m)
!!$                      okay = .TRUE.
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$                IF (okay) THEN
!!$                   WRITE(*,*) 'add ', j_tmp(idx), 'in dom ', inputs%list_dom_phi(i), ' : proc ', rank
!!$                   WRITE(*,*) 'add ', mesh%rr(:,j_tmp(idx)), mesh%i_d(m)
!!$                   EXIT
!!$                END IF
!!$             END IF
!!$          END DO
!!$       END IF
!!$    END DO
!!$
!!$    nb_cav = idx - SIZE(js_D)
!!$    IF (nb_cav /= 0) THEN
!!$       DEALLOCATE(js_D)
!!$       ALLOCATE(js_D(idx))
!!$       js_D = j_tmp(1:idx)
!!$    END IF
!!$
!!$    DEALLOCATE(on_proc_loc, on_proc, j_tmp)
!!$    DEALLOCATE(not_cav_loc, not_cav)
!!$    IF (ALLOCATED(is_ok)) DEALLOCATE(is_ok)
!!$
!!$    WRITE(*,'(a,x,i2,x,a,x,i2)') 'I have detected', nb_cav, ' cavity(ies) on proc', rank
!!$
!!$  END SUBROUTINE dirichlet_cavities

!!$  SUBROUTINE smb_sigma_prod_curl(communicator, mesh, jj_v_to_H, list_mode, H_in, one_over_sigma_in, sigma_nj_m,&
!!$       sigma, V_out)
!!$    !=================================
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE user_data
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode 
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: H_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_in
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%me), INTENT(IN)  :: sigma_nj_m
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:)              :: V_out
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: H_gauss, RotH
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,2,SIZE(list_mode)) :: one_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_G*mesh%me,6,SIZE(list_mode)) :: RotH_bar
!!$    INTEGER,      DIMENSION(mesh%gauss%n_w)                           :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc     
!!$    INTEGER                                                           :: m, l , i, mode, index, k
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6)                         :: H_in_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2)                         :: one_over_sigma_in_loc
!!$    REAL(KIND=8)   :: ray, sigma_np_gauss
!!$    INTEGER       :: nb_procs, bloc_size, m_max_pad, code
!!$    MPI_Comm       :: communicator
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO m = 1, mesh%me
!!$          j_loc = mesh%jj(:,m)
!!$          DO k = 1, 6
!!$             H_in_loc(:,k) = H_in(j_loc,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc(:,k) = one_over_sigma_in(j_loc,k,i)
!!$          END DO
!!$
!!$          DO l = 1, mesh%gauss%l_G
!!$             index = index + 1
!!$             dw_loc = mesh%gauss%dw(:,:,l,m)
!!$
!!$             !===Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
!!$
!!$             !-----------------magnetic field on gauss points---------------------------
!!$             H_gauss(index,1,i) = SUM(H_in_loc(:,1)*mesh%gauss%ww(:,l))
!!$             H_gauss(index,3,i) = SUM(H_in_loc(:,3)*mesh%gauss%ww(:,l))
!!$             H_gauss(index,5,i) = SUM(H_in_loc(:,5)*mesh%gauss%ww(:,l))
!!$
!!$             H_gauss(index,2,i) = SUM(H_in_loc(:,2)*mesh%gauss%ww(:,l))
!!$             H_gauss(index,4,i) = SUM(H_in_loc(:,4)*mesh%gauss%ww(:,l))
!!$             H_gauss(index,6,i) = SUM(H_in_loc(:,6)*mesh%gauss%ww(:,l))
!!$             !-----------------Curl of H on gauss points--------------------------------
!!$             !coeff sur les cosinus 
!!$             RotH(index,1,i) = mode/ray*H_gauss(index,6,i) &
!!$                  -SUM(H_in_loc(:,3)*dw_loc(2,:))
!!$             RotH(index,4,i) = SUM(H_in_loc(:,2)*dw_loc(2,:)) &
!!$                  -SUM(H_in_loc(:,6)*dw_loc(1,:))
!!$             RotH(index,5,i) = 1/ray*H_gauss(index,3,i) &
!!$                  +SUM(H_in_loc(:,3)*dw_loc(1,:)) &
!!$                  -mode/ray*H_gauss(index,2,i)
!!$             !coeff sur les sinus       
!!$             RotH(index,2,i) =-mode/ray*H_gauss(index,5,i) &
!!$                  -SUM(H_in_loc(:,4)*dw_loc(2,:))
!!$             RotH(index,3,i) = SUM(H_in_loc(:,1)*dw_loc(2,:)) &
!!$                  -SUM(H_in_loc(:,5)*dw_loc(1,:))
!!$             RotH(index,6,i) = 1/ray*H_gauss(index,4,i) &
!!$                  +SUM(H_in_loc(:,4)*dw_loc(1,:))&
!!$                  +mode/ray*H_gauss(index,1,i)
!!$             !-----------------one over sigma on gauss points---------------------------
!!$             IF (jj_v_to_H(mesh%jj(1,m)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m)
!!$                END IF
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc(:,1)*mesh%gauss%ww(:,l))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc(:,2)*mesh%gauss%ww(:,l))
!!$             END IF
!!$             !-----------------RotHbar on gauss points----------------------------------
!!$             sigma_np_gauss = SUM(sigma_nj_m(:,m)*mesh%gauss%ww(:,l))
!!$             DO k = 1, 6
!!$                RotH_bar(index,k,i) = RotH(index,k,i)/sigma_np_gauss             
!!$             END DO
!!$          ENDDO
!!$       ENDDO
!!$    END DO
!!$
!!$    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$    bloc_size = SIZE(one_over_sigma_gauss,1)/nb_procs+1
!!$    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$    CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$         RotH, one_over_sigma_gauss, V_out, 1, nb_procs, bloc_size, m_max_pad)
!!$
!!$    V_out = RotH_bar - V_out
!!$
!!$  END SUBROUTINE smb_sigma_prod_curl

!!$  SUBROUTINE smb_sigma_prod_curl_bdy(communicator, mesh, jj_v_to_H, Dirichlet_bdy_H_sides, list_mode, &
!!$       H_in, one_over_sigma_in, sigma_np, sigma, V_out)
!!$    !=================================
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE user_data
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: Dirichlet_bdy_H_sides
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode 
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: H_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_in
!!$    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: sigma_np
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:)              :: V_out
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*SIZE(Dirichlet_bdy_H_sides),6,SIZE(list_mode)) :: H_gauss, RotH
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*SIZE(Dirichlet_bdy_H_sides),6,SIZE(list_mode)) :: RotH_bar
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%l_Gs*SIZE(Dirichlet_bdy_H_sides),2,SIZE(list_mode)) :: one_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc     
!!$    INTEGER                                                           :: ms, ls , i, mode, index, k
!!$    INTEGER,      DIMENSION(mesh%gauss%n_ws)                          :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6)                        :: H_in_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2)                        :: one_over_sigma_in_loc
!!$    REAL(KIND=8)   :: ray
!!$    INTEGER       :: nb_procs, bloc_size, m_max_pad, code, count, m1
!!$    MPI_Comm       :: communicator
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO count = 1, SIZE(Dirichlet_bdy_H_sides)
!!$          ms = Dirichlet_bdy_H_sides(count)
!!$          m1 = mesh%neighs(ms)
!!$
!!$          j_loc = mesh%jjs(:,ms)
!!$          DO k = 1, 6
!!$             H_in_loc(:,k) = H_in(j_loc,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc(:,k) = one_over_sigma_in(j_loc,k,i)
!!$          END DO
!!$
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             index = index + 1
!!$             dw_loc = mesh%gauss%dw_s(:,:,ls,ms)
!!$
!!$             !===Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
!!$
!!$             !-----------------magnetic field on bdy gauss points---------------------------
!!$             H_gauss(index,1,i) = SUM(H_in_loc(:,1)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,3,i) = SUM(H_in_loc(:,3)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,5,i) = SUM(H_in_loc(:,5)*mesh%gauss%wws(:,ls))
!!$
!!$             H_gauss(index,2,i) = SUM(H_in_loc(:,2)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,4,i) = SUM(H_in_loc(:,4)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,6,i) = SUM(H_in_loc(:,6)*mesh%gauss%wws(:,ls))
!!$             !-----------------Curl of H on bdy gauss points--------------------------------
!!$             !coeff sur les cosinus 
!!$             RotH(index,1,i) = mode/ray*H_gauss(index,6,i) &
!!$                  -SUM(H_in(mesh%jj(:,m1),3,i)*dw_loc(2,:))
!!$
!!$             RotH(index,4,i) = SUM(H_in(mesh%jj(:,m1),2,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m1),6,i)*dw_loc(1,:))
!!$
!!$             RotH(index,5,i) = 1/ray*H_gauss(index,3,i) &
!!$                  +SUM(H_in(mesh%jj(:,m1),3,i)*dw_loc(1,:)) &
!!$                  -mode/ray*H_gauss(index,2,i)
!!$
!!$             !coeff sur les sinus       
!!$             RotH(index,2,i) =-mode/ray*H_gauss(index,5,i) &
!!$                  -SUM(H_in(mesh%jj(:,m1),4,i)*dw_loc(2,:))
!!$
!!$             RotH(index,3,i) = SUM(H_in(mesh%jj(:,m1),1,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m1),5,i)*dw_loc(1,:))
!!$
!!$             RotH(index,6,i) = 1/ray*H_gauss(index,4,i) &
!!$                  +SUM(H_in(mesh%jj(:,m1),4,i)*dw_loc(1,:))&
!!$                  +mode/ray*H_gauss(index,1,i)
!!$             !-----------------one over sigma and RotH_bar on bdy gauss points--------------
!!$             IF (jj_v_to_H(mesh%jj(1,m1)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m1)
!!$                END IF
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/sigma(m1)
!!$                END DO
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc(:,1)*mesh%gauss%wws(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc(:,2)*mesh%gauss%wws(:,ls))
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/SUM(sigma_np(mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
!!$                END DO
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF ( SIZE(Dirichlet_bdy_H_sides).GE.1) THEN
!!$       CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$       bloc_size = SIZE(one_over_sigma_gauss,1)/nb_procs+1
!!$       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$       CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$            RotH, one_over_sigma_gauss, V_out, 1, nb_procs, bloc_size, m_max_pad)
!!$
!!$       V_out =  RotH_bar - V_out
!!$    END IF
!!$
!!$  END SUBROUTINE smb_sigma_prod_curl_bdy

!!$  SUBROUTINE smb_sigma_prod_curl_inter_mu(communicator, mesh, jj_v_to_H, interface_H_mu, list_mode, &
!!$       H_in, one_over_sigma_in, sigma_np, sigma, V_out)
!!$    !=================================
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE user_data
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    TYPE(interface_type),           INTENT(IN)  :: interface_H_mu
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode 
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: H_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_in
!!$    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: sigma_np
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:)              :: V_out
!!$    REAL(KIND=8), DIMENSION(2*mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode)) :: H_gauss, RotH
!!$    REAL(KIND=8), DIMENSION(2*mesh%gauss%l_Gs*interface_H_mu%mes,6,SIZE(list_mode)) :: RotH_bar
!!$    REAL(KIND=8), DIMENSION(2*mesh%gauss%l_Gs*interface_H_mu%mes,2,SIZE(list_mode)) :: one_over_sigma_gauss  
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)            :: dw_loc     
!!$    INTEGER                                                           :: ms, ls , i, mode, index, k
!!$    INTEGER,      DIMENSION(mesh%gauss%n_ws)                          :: j_loc1, j_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6)                        :: H_in_loc1, H_in_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2)                        :: one_over_sigma_in_loc1, one_over_sigma_in_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,mesh%gauss%l_Gs)      :: w_cs
!!$    REAL(KIND=8)  :: ray, diff, ref
!!$    INTEGER       :: nb_procs, bloc_size, m_max_pad, code
!!$    INTEGER       :: ms1, ms2, m1, m2, mesh_id1, mesh_id2
!!$    MPI_Comm       :: communicator
!!$
!!$    DO ms = 1, interface_H_mu%mes
!!$       ms1 = interface_H_mu%mesh1(ms)
!!$       ms2 = interface_H_mu%mesh2(ms)
!!$       m1 = mesh%neighs(ms1)
!!$       m2 = mesh%neighs(ms2)
!!$       ref = 1.d-8+SUM((mesh%rr(:,mesh%jjs(1,ms1)) - mesh%rr(:,mesh%jjs(2,ms1)))**2)
!!$       diff =      SUM((mesh%rr(:,mesh%jjs(1,ms1)) - mesh%rr(:,mesh%jjs(1,ms2)))**2)
!!$       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
!!$          w_cs = mesh%gauss%wws
!!$       ELSE                ! 1 = 2
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             w_cs(1,ls)= mesh%gauss%wws(2,ls)
!!$             w_cs(2,ls)= mesh%gauss%wws(1,ls)
!!$             IF (mesh%gauss%n_ws==3) w_cs(mesh%gauss%n_ws,ls) = mesh%gauss%wws(mesh%gauss%n_ws,ls) 
!!$             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
!!$          END DO
!!$       END IF
!!$    END DO
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO ms = 1, interface_H_mu%mes
!!$          ms2 = interface_H_mu%mesh2(ms)
!!$          ms1 = interface_H_mu%mesh1(ms)
!!$          m2 = mesh%neighs(ms2)
!!$          m1 = mesh%neighs(ms1)
!!$          mesh_id1 = mesh%i_d(m1)
!!$          mesh_id2 = mesh%i_d(m2)
!!$          j_loc1 = mesh%jjs(:,ms1)
!!$          j_loc2 = mesh%jjs(:,ms2)          
!!$          DO k = 1, 6
!!$             H_in_loc1(:,k) = H_in(j_loc1,k,i)
!!$             H_in_loc2(:,k) = H_in(j_loc2,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc1(:,k) = one_over_sigma_in(j_loc1,k,i)
!!$             one_over_sigma_in_loc2(:,k) = one_over_sigma_in(j_loc2,k,i)
!!$          END DO
!!$
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             !===Side 1
!!$             index = index + 1
!!$             dw_loc = mesh%gauss%dw_s(:,:,ls,ms1)
!!$             !Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms1))*w_cs(:,ls))
!!$             !-----------------magnetic field on bdy gauss points---------------------------
!!$             H_gauss(index,1,i) = SUM(H_in_loc1(:,1)*w_cs(:,ls))
!!$             H_gauss(index,3,i) = SUM(H_in_loc1(:,3)*w_cs(:,ls))
!!$             H_gauss(index,5,i) = SUM(H_in_loc1(:,5)*w_cs(:,ls))
!!$             H_gauss(index,2,i) = SUM(H_in_loc1(:,2)*w_cs(:,ls))
!!$             H_gauss(index,4,i) = SUM(H_in_loc1(:,4)*w_cs(:,ls))
!!$             H_gauss(index,6,i) = SUM(H_in_loc1(:,6)*w_cs(:,ls))
!!$             !-----------------Curl of H on bdy gauss points--------------------------------
!!$             !coeff sur les cosinus 
!!$             RotH(index,1,i) = mode/ray*H_gauss(index,6,i) &
!!$                  -SUM(H_in(mesh%jj(:,m1),3,i)*dw_loc(2,:))
!!$             RotH(index,4,i) = SUM(H_in(mesh%jj(:,m1),2,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m1),6,i)*dw_loc(1,:))
!!$             RotH(index,5,i) = 1/ray*H_gauss(index,3,i) &
!!$                  +SUM(H_in(mesh%jj(:,m1),3,i)*dw_loc(1,:)) &
!!$                  -mode/ray*H_gauss(index,2,i)
!!$             !coeff sur les sinus       
!!$             RotH(index,2,i) =-mode/ray*H_gauss(index,5,i) &
!!$                  -SUM(H_in(mesh%jj(:,m1),4,i)*dw_loc(2,:))
!!$             RotH(index,3,i) = SUM(H_in(mesh%jj(:,m1),1,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m1),5,i)*dw_loc(1,:))
!!$             RotH(index,6,i) = 1/ray*H_gauss(index,4,i) &
!!$                  +SUM(H_in(mesh%jj(:,m1),4,i)*dw_loc(1,:))&
!!$                  +mode/ray*H_gauss(index,1,i)
!!$             !-----------------one over sigma and RotH_bar on bdy gauss points--------------
!!$             IF (jj_v_to_H(mesh%jj(1,m1)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m1)
!!$                END IF
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/sigma(m1)
!!$                END DO
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc1(:,1)*w_cs(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc1(:,2)*w_cs(:,ls))
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/SUM(sigma_np(mesh%jjs(:,ms1))*w_cs(:,ls))
!!$                END DO
!!$             END IF
!!$
!!$             !===Side 2
!!$             index = index + 1
!!$             dw_loc = mesh%gauss%dw_s(:,:,ls,ms2)
!!$             !Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms2))*mesh%gauss%wws(:,ls))
!!$             !-----------------magnetic field on bdy gauss points---------------------------
!!$             H_gauss(index,1,i) = SUM(H_in_loc2(:,1)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,3,i) = SUM(H_in_loc2(:,3)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,5,i) = SUM(H_in_loc2(:,5)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,2,i) = SUM(H_in_loc2(:,2)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,4,i) = SUM(H_in_loc2(:,4)*mesh%gauss%wws(:,ls))
!!$             H_gauss(index,6,i) = SUM(H_in_loc2(:,6)*mesh%gauss%wws(:,ls))
!!$             !-----------------Curl of H on bdy gauss points--------------------------------
!!$             !coeff sur les cosinus 
!!$             RotH(index,1,i) = mode/ray*H_gauss(index,6,i) &
!!$                  -SUM(H_in(mesh%jj(:,m2),3,i)*dw_loc(2,:))
!!$             RotH(index,4,i) = SUM(H_in(mesh%jj(:,m2),2,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m2),6,i)*dw_loc(1,:))
!!$             RotH(index,5,i) = 1/ray*H_gauss(index,3,i) &
!!$                  +SUM(H_in(mesh%jj(:,m2),3,i)*dw_loc(1,:)) &
!!$                  -mode/ray*H_gauss(index,2,i)
!!$             !coeff sur les sinus       
!!$             RotH(index,2,i) =-mode/ray*H_gauss(index,5,i) &
!!$                  -SUM(H_in(mesh%jj(:,m2),4,i)*dw_loc(2,:))
!!$             RotH(index,3,i) = SUM(H_in(mesh%jj(:,m2),1,i)*dw_loc(2,:)) &
!!$                  -SUM(H_in(mesh%jj(:,m2),5,i)*dw_loc(1,:))
!!$             RotH(index,6,i) = 1/ray*H_gauss(index,4,i) &
!!$                  +SUM(H_in(mesh%jj(:,m2),4,i)*dw_loc(1,:))&
!!$                  +mode/ray*H_gauss(index,1,i)
!!$             !-----------------one over sigma and RotH_bar on bdy gauss points--------------
!!$             IF (jj_v_to_H(mesh%jj(1,m2)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m2)
!!$                END IF
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/sigma(m2)
!!$                END DO
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc2(:,1)*mesh%gauss%wws(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc2(:,2)*mesh%gauss%wws(:,ls))
!!$                DO k = 1, 6
!!$                   RotH_bar(index,k,i) = RotH(index,k,i)/SUM(sigma_np(mesh%jjs(:,ms2))*mesh%gauss%wws(:,ls))
!!$                END DO
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF (interface_H_mu%mes.GE.1) THEN
!!$       CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$       bloc_size = SIZE(one_over_sigma_gauss,1)/nb_procs+1
!!$       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$       CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$            RotH, one_over_sigma_gauss, V_out, 1, nb_procs, bloc_size, m_max_pad)
!!$
!!$       V_out =  RotH_bar - V_out
!!$    END IF
!!$
!!$  END SUBROUTINE smb_sigma_prod_curl_inter_mu

!!$  SUBROUTINE smb_current_over_sigma(communicator, mesh, jj_v_to_H, list_mode, B_in, &
!!$       mu_H_field, mu_phi, one_over_sigma_tot, time, sigma, J_over_sigma_gauss)
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE boundary
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: B_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_tot
!!$    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: mu_H_field
!!$    REAL(KIND=8),                   INTENT(IN)  :: mu_phi, time
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: J_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(mesh%me*mesh%gauss%l_G,6,SIZE(list_mode)) :: J_exact_gauss
!!$    REAL(KIND=8), DIMENSION(mesh%me*mesh%gauss%l_G,2,SIZE(list_mode)) :: one_over_sigma_gauss
!!$    INTEGER,      DIMENSION(mesh%gauss%n_w)   :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,6) :: B_in_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2) :: one_over_sigma_in_loc
!!$    REAL(KIND=8), DIMENSION(6)  :: B_ext_l
!!$    REAL(KIND=8)                :: muhl, ray
!!$    REAL(KIND=8), DIMENSION(2) :: gaussp
!!$    INTEGER               :: mode, mesh_id1, k, i, index, l, m, ni
!!$    INTEGER               :: nb_procs, bloc_size, m_max_pad, code
!!$    MPI_Comm       :: communicator
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO m = 1, mesh%me
!!$          mesh_id1 = mesh%i_d(m) 
!!$          j_loc = mesh%jj(:,m)
!!$          DO k = 1, 6
!!$             B_in_loc(:,k) = B_in(j_loc,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc(:,k) = one_over_sigma_tot(j_loc,k,i)
!!$          END DO
!!$          DO l = 1, mesh%gauss%l_G
!!$             index = index + 1
!!$
!!$             !===Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
!!$
!!$             !-----------------Variable for Jexact on gauss points----------------------
!!$             B_ext_l(1) = SUM(B_in_loc(:,1)*mesh%gauss%ww(:,l))
!!$             B_ext_l(3) = SUM(B_in_loc(:,3)*mesh%gauss%ww(:,l))
!!$             B_ext_l(5) = SUM(B_in_loc(:,5)*mesh%gauss%ww(:,l))
!!$             B_ext_l(2) = SUM(B_in_loc(:,2)*mesh%gauss%ww(:,l))
!!$             B_ext_l(4) = SUM(B_in_loc(:,4)*mesh%gauss%ww(:,l))
!!$             B_ext_l(6) = SUM(B_in_loc(:,6)*mesh%gauss%ww(:,l))
!!$             gaussp = 0.d0
!!$             DO ni = 1, mesh%gauss%n_w
!!$                gaussp = gaussp + mesh%rr(:,mesh%jj(ni,m))*mesh%gauss%ww(ni,l)         
!!$             ENDDO
!!$             muhl=SUM(mu_H_field(mesh%jj(:,m))*mesh%gauss%ww(:,l))
!!$             !-----------------J_exact on gauss points----------------------------------
!!$             DO k = 1, 6
!!$                J_exact_gauss(index,k,i)=Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m), &
!!$                     muhl, time, mesh_id1, B_ext_l)
!!$             END DO
!!$             !-----------------one over sigma on gauss points---------------------------
!!$             IF (jj_v_to_H(mesh%jj(1,m)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m)
!!$                END IF
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc(:,1)*mesh%gauss%ww(:,l))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc(:,2)*mesh%gauss%ww(:,l))
!!$             END IF
!!$
!!$          ENDDO
!!$       ENDDO
!!$    END DO
!!$
!!$    CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$    bloc_size = SIZE(J_exact_gauss,1)/nb_procs+1
!!$    m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$    CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$         J_exact_gauss, one_over_sigma_gauss, J_over_sigma_gauss,&
!!$         1, nb_procs, bloc_size, m_max_pad)
!!$
!!$  END SUBROUTINE smb_current_over_sigma

!!$  SUBROUTINE smb_current_over_sigma_bdy(communicator, mesh, jj_v_to_H, Dirichlet_bdy_H_sides, list_mode, B_in, &
!!$       mu_H_field, mu_phi, one_over_sigma_tot, time, sigma, J_over_sigma_gauss)
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE boundary
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: Dirichlet_bdy_H_sides
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: B_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_tot
!!$    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: mu_H_field
!!$    REAL(KIND=8),                   INTENT(IN)  :: mu_phi, time
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: J_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(SIZE(Dirichlet_bdy_H_sides)*mesh%gauss%l_Gs,6,SIZE(list_mode)) :: J_exact_gauss
!!$    REAL(KIND=8), DIMENSION(SIZE(Dirichlet_bdy_H_sides)*mesh%gauss%l_Gs,2,SIZE(list_mode)) :: one_over_sigma_gauss
!!$    INTEGER,      DIMENSION(mesh%gauss%n_ws)   :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6) :: B_in_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2) :: one_over_sigma_in_loc
!!$    REAL(KIND=8), DIMENSION(6) :: B_ext_l
!!$    REAL(KIND=8)               :: muhl, ray
!!$    REAL(KIND=8), DIMENSION(2) :: gaussp
!!$    INTEGER               :: mode, mesh_id1, k, i, count, index, ls, ms, ni
!!$    INTEGER               :: nb_procs, bloc_size, m_max_pad, code, m1
!!$    MPI_Comm       :: communicator
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO count = 1, SIZE(Dirichlet_bdy_H_sides)
!!$          ms = Dirichlet_bdy_H_sides(count)
!!$          m1 = mesh%neighs(ms)
!!$          mesh_id1 = mesh%i_d(m1) 
!!$          j_loc = mesh%jjs(:,ms)
!!$          DO k = 1, 6
!!$             B_in_loc(:,k) = B_in(j_loc,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc(:,k) = one_over_sigma_tot(j_loc,k,i)
!!$          END DO
!!$
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             index = index + 1
!!$
!!$             !===Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
!!$
!!$             !-----------------Variable for Jexact on gauss points----------------------
!!$             B_ext_l(1) = SUM(B_in_loc(:,1)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(3) = SUM(B_in_loc(:,3)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(5) = SUM(B_in_loc(:,5)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(2) = SUM(B_in_loc(:,2)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(4) = SUM(B_in_loc(:,4)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(6) = SUM(B_in_loc(:,6)*mesh%gauss%wws(:,ls))
!!$             gaussp = 0.d0
!!$             DO ni = 1, mesh%gauss%n_ws
!!$                gaussp = gaussp + mesh%rr(:,mesh%jjs(ni,ms))*mesh%gauss%wws(ni,ls)         
!!$             ENDDO
!!$             muhl=SUM(mu_H_field(mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
!!$             !-----------------J_exact on gauss points----------------------------------
!!$             DO k = 1, 6
!!$                J_exact_gauss(index,k,i)=Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m1),&
!!$                     muhl, time, mesh_id1, B_ext_l)
!!$             END DO
!!$             !-----------------one over sigma on gauss points---------------------------
!!$             IF (jj_v_to_H(mesh%jj(1,m1)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m1)
!!$                END IF
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc(:,1)*mesh%gauss%wws(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc(:,2)*mesh%gauss%wws(:,ls))
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF ( SIZE(Dirichlet_bdy_H_sides).GE.1) THEN
!!$       CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$       bloc_size = SIZE(J_exact_gauss,1)/nb_procs+1
!!$       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$       CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$            J_exact_gauss, one_over_sigma_gauss, J_over_sigma_gauss,&
!!$            1, nb_procs, bloc_size, m_max_pad)
!!$    END IF
!!$
!!$  END SUBROUTINE smb_current_over_sigma_bdy

!!$  SUBROUTINE smb_current_over_sigma_inter_mu(communicator, mesh, jj_v_to_H, interface_H_mu, list_mode, B_in, &
!!$       mu_H_field, mu_phi, one_over_sigma_tot, time, sigma, J_over_sigma_gauss)
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE boundary
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:)    , INTENT(IN)  :: jj_v_to_H
!!$    TYPE(interface_type),           INTENT(IN)  :: interface_H_mu    
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: B_in
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_tot
!!$    REAL(KIND=8), DIMENSION(:),     INTENT(IN)  :: mu_H_field
!!$    REAL(KIND=8),                   INTENT(IN)  :: mu_phi, time
!!$    REAL(KIND=8), DIMENSION(mesh%me),INTENT(IN) :: sigma
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: J_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(2*interface_H_mu%mes*mesh%gauss%l_Gs,6,SIZE(list_mode)) :: J_exact_gauss
!!$    REAL(KIND=8), DIMENSION(2*interface_H_mu%mes*mesh%gauss%l_Gs,2,SIZE(list_mode)) :: one_over_sigma_gauss
!!$    REAL(KIND=8), DIMENSION(6) :: B_ext_l
!!$    REAL(KIND=8)               :: muhl, diff, ref, ray
!!$    REAL(KIND=8), DIMENSION(2) :: gaussp
!!$    INTEGER               :: mode, k, i, mesh_id1, mesh_id2, ni
!!$    INTEGER               :: nb_procs, bloc_size, m_max_pad, code     
!!$    INTEGER                                         :: ms, ms1, ms2, m1, m2, ls, index
!!$    INTEGER,      DIMENSION(mesh%gauss%n_ws)        :: j_loc1, j_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,6)      :: B_in_loc1, B_in_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2)      :: one_over_sigma_in_loc1,one_over_sigma_in_loc2
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,mesh%gauss%l_Gs) :: w_cs
!!$    MPI_Comm       :: communicator
!!$
!!$
!!$    DO ms = 1, interface_H_mu%mes
!!$       ms1 = interface_H_mu%mesh1(ms)
!!$       ms2 = interface_H_mu%mesh2(ms)
!!$       m1 = mesh%neighs(ms1)
!!$       m2 = mesh%neighs(ms2)
!!$       ref = 1.d-8+SUM((mesh%rr(:,mesh%jjs(1,ms1)) - mesh%rr(:,mesh%jjs(2,ms1)))**2)
!!$       diff =      SUM((mesh%rr(:,mesh%jjs(1,ms1)) - mesh%rr(:,mesh%jjs(1,ms2)))**2)
!!$       IF (diff/ref .LT. 1.d-10) THEN ! 1 = 1
!!$          w_cs = mesh%gauss%wws
!!$       ELSE                ! 1 = 2
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             w_cs(1,ls)= mesh%gauss%wws(2,ls)
!!$             w_cs(2,ls)= mesh%gauss%wws(1,ls)
!!$             IF (mesh%gauss%n_ws==3) w_cs(mesh%gauss%n_ws,ls) = mesh%gauss%wws(mesh%gauss%n_ws,ls) 
!!$             WRITE(*,*) ' Ouaps! oder of shape functions changed?'
!!$          END DO
!!$       END IF
!!$    END DO
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO ms = 1, interface_H_mu%mes
!!$          ms2 = interface_H_mu%mesh2(ms)
!!$          ms1 = interface_H_mu%mesh1(ms)
!!$          m2 = mesh%neighs(ms2)
!!$          m1 = mesh%neighs(ms1)
!!$          mesh_id1 = mesh%i_d(m1)
!!$          mesh_id2 = mesh%i_d(m2)
!!$          j_loc1 = mesh%jjs(:,ms1)
!!$          j_loc2 = mesh%jjs(:,ms2)          
!!$          DO k = 1, 6
!!$             B_in_loc1(:,k) = B_in(j_loc1,k,i)
!!$             B_in_loc2(:,k) = B_in(j_loc2,k,i)
!!$          END DO
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc1(:,k) = one_over_sigma_tot(j_loc1,k,i)
!!$             one_over_sigma_in_loc2(:,k) = one_over_sigma_tot(j_loc2,k,i)
!!$          END DO
!!$
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             !===Side 1
!!$             index = index + 1
!!$             !Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms1))*w_cs(:,ls))
!!$             !-----------------Variable for Jexact on gauss points----------------------
!!$             B_ext_l(1) = SUM(B_in_loc1(:,1)*w_cs(:,ls))
!!$             B_ext_l(3) = SUM(B_in_loc1(:,3)*w_cs(:,ls))
!!$             B_ext_l(5) = SUM(B_in_loc1(:,5)*w_cs(:,ls))
!!$             B_ext_l(2) = SUM(B_in_loc1(:,2)*w_cs(:,ls))
!!$             B_ext_l(4) = SUM(B_in_loc1(:,4)*w_cs(:,ls))
!!$             B_ext_l(6) = SUM(B_in_loc1(:,6)*w_cs(:,ls))
!!$             gaussp = 0.d0
!!$             DO ni = 1, mesh%gauss%n_ws
!!$                gaussp = gaussp + mesh%rr(:,mesh%jjs(ni,ms1))*w_cs(ni,ls)
!!$             ENDDO
!!$             muhl=SUM(mu_H_field(mesh%jjs(:,ms1))*w_cs(:,ls))
!!$             !-----------------J_exact on gauss points----------------------------------
!!$             DO k = 1, 6
!!$                J_exact_gauss(index,k,i)=Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m1),&
!!$                     muhl, time, mesh_id1, B_ext_l)
!!$             END DO
!!$             !-----------------one over sigma on gauss points---------------------------
!!$             IF (jj_v_to_H(mesh%jj(1,m1)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m1)
!!$                END IF
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc1(:,1)*w_cs(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc1(:,2)*w_cs(:,ls))
!!$             END IF
!!$
!!$             !===Side 2
!!$             index = index + 1
!!$             !Compute radius of Gauss point
!!$             ray = SUM(mesh%rr(1,mesh%jjs(:,ms2))*mesh%gauss%wws(:,ls))
!!$             !-----------------Variable for Jexact on gauss points----------------------
!!$             B_ext_l(1) = SUM(B_in_loc2(:,1)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(3) = SUM(B_in_loc2(:,3)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(5) = SUM(B_in_loc2(:,5)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(2) = SUM(B_in_loc2(:,2)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(4) = SUM(B_in_loc2(:,4)*mesh%gauss%wws(:,ls))
!!$             B_ext_l(6) = SUM(B_in_loc2(:,6)*mesh%gauss%wws(:,ls))
!!$             gaussp = 0.d0
!!$             DO ni = 1, mesh%gauss%n_ws
!!$                gaussp = gaussp + mesh%rr(:,mesh%jjs(ni,ms2))*mesh%gauss%wws(ni,ls)
!!$             ENDDO
!!$             muhl=SUM(mu_H_field(mesh%jjs(:,ms2))*mesh%gauss%wws(:,ls))
!!$             !-----------------J_exact on gauss points----------------------------------
!!$             DO k = 1, 6
!!$                J_exact_gauss(index,k,i)=Jexact_gauss(k, gaussp, mode, mu_phi, sigma(m2),&
!!$                     muhl, time, mesh_id2, B_ext_l)
!!$             END DO
!!$             !-----------------one over sigma on gauss points---------------------------
!!$             IF (jj_v_to_H(mesh%jj(1,m2)) == -1) THEN
!!$                one_over_sigma_gauss(index,1,i) = 0.d0
!!$                one_over_sigma_gauss(index,2,i) = 0.d0
!!$                IF (mode==0) THEN
!!$                   one_over_sigma_gauss(index,1,i) = 1.d0/sigma(m2)
!!$                END IF
!!$             ELSE
!!$                one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc2(:,1)*mesh%gauss%wws(:,ls))
!!$                one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc2(:,2)*mesh%gauss%wws(:,ls))
!!$             END IF
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF (interface_H_mu%mes.GE.1) THEN
!!$       CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$       bloc_size = SIZE(J_exact_gauss,1)/nb_procs+1
!!$       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$       CALL FFT_SCALAR_VECT_NO_OVERSHOOT(communicator, 1.d0/(inputs%sigma_fluid*inputs%Rem),  &
!!$            J_exact_gauss, one_over_sigma_gauss, J_over_sigma_gauss,&
!!$            1, nb_procs, bloc_size, m_max_pad)
!!$    END IF
!!$
!!$  END SUBROUTINE smb_current_over_sigma_inter_mu

!!$  SUBROUTINE smb_sigma_Neumann(communicator, mesh, Neumann_bdy_H_sides, list_mode, &
!!$       one_over_sigma_tot, sigma_tot_gauss_Neumann)
!!$    USE sft_parallele
!!$    USE chaine_caractere
!!$    USE input_data
!!$    USE def_type_mesh
!!$    USE boundary
!!$#include "petsc/finclude/petsc.h"
!!$    USE petsc
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),                INTENT(IN)  :: mesh
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: Neumann_bdy_H_sides
!!$    INTEGER,      DIMENSION(:),     INTENT(IN)  :: list_mode
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)  :: one_over_sigma_tot
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: sigma_tot_gauss_Neumann
!!$    REAL(KIND=8), DIMENSION(SIZE(Neumann_bdy_H_sides)*mesh%gauss%l_Gs,2,SIZE(list_mode)):: one_over_sigma_gauss
!!$    INTEGER,      DIMENSION(mesh%gauss%n_ws)   :: j_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2) :: one_over_sigma_in_loc
!!$    INTEGER               :: mode, k, i, count, index, ls, ms
!!$    INTEGER               :: nb_procs, bloc_size, m_max_pad, code
!!$    MPI_Comm       :: communicator
!!$
!!$    DO i = 1, SIZE(list_mode)
!!$       mode = list_mode(i)
!!$       index = 0
!!$       DO count = 1, SIZE(Neumann_bdy_H_sides)
!!$          ms = Neumann_bdy_H_sides(count)
!!$          j_loc = mesh%jjs(:,ms)
!!$          DO k = 1, 2
!!$             one_over_sigma_in_loc(:,k) = one_over_sigma_tot(j_loc,k,i)
!!$          END DO
!!$
!!$          DO ls = 1, mesh%gauss%l_Gs
!!$             index = index + 1
!!$             one_over_sigma_gauss(index,1,i) = SUM(one_over_sigma_in_loc(:,1)*mesh%gauss%wws(:,ls))
!!$             one_over_sigma_gauss(index,2,i) = SUM(one_over_sigma_in_loc(:,2)*mesh%gauss%wws(:,ls))
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    sigma_tot_gauss_Neumann=one_over_sigma_gauss
!!$
!!$    IF ( SIZE(Neumann_bdy_H_sides).GE.1) THEN
!!$       CALL MPI_COMM_SIZE(communicator, nb_procs, code)
!!$       bloc_size = SIZE(one_over_sigma_gauss,1)/nb_procs+1
!!$       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
!!$       CALL FFT_PAR_SCAL_FUNCT(communicator, sigma_tot_gauss_Neumann, one_over_x, nb_procs, bloc_size, m_max_pad)
!!$    END IF
!!$
!!$  CONTAINS
!!$    FUNCTION one_over_x(x) RESULT(vv)
!!$      IMPLICIT NONE
!!$      REAL(KIND=8) :: x
!!$      REAL(KIND=8) :: vv
!!$
!!$      vv = x/ MAX(x*x, 1.d-20)
!!$
!!$    END FUNCTION one_over_x
!!$
!!$  END SUBROUTINE smb_sigma_Neumann

END MODULE update_maxwell_mxs_with_H
