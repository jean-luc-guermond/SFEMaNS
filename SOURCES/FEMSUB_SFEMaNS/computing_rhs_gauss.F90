MODULE rhs_gauss_computing
  PUBLIC:: rhs_ns_gauss_3x3, rhs_residual_ns_gauss_3x3, rhs_residual_ns_gauss_3x3_mom,&
       rhs_ns_gauss_3x3_art_comp_mom
  PRIVATE
CONTAINS

  SUBROUTINE rhs_ns_gauss_3x3(vv_mesh, pp_mesh, communicator, list_mode, time, V1m, pn, pn_inc, rotv_v, &
       rhs_gauss, density, tempn, concn)
    !=================================
    !RHS for Navier-Stokes
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
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: pn, pn_inc
    INTEGER,      DIMENSION(:),                 INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)                                   :: fs, ft, fp_inc
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8), DIMENSION(vv_mesh%np,6)                        :: ff, imposed_vel
    REAL(KIND=8), DIMENSION(vv_mesh%np,2)                        :: P, P_inc
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

    DO i = 1, SIZE(list_mode)
       DO k= 1, 6
          ff(:,k) = source_in_NS_momentum(k, vv_mesh%rr, list_mode(i), i, time,  inputs%Re, 'ns', &
               density, tempn, concn)
       END DO
       DO k = 1, 2
          !===JLG+HF Dec 10 2019
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn_inc(:,k,i), P_inc(:,k))
          !CALL inject(pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
          !CALL inject(pp_mesh%jj, vv_mesh%jj, pn_inc(:,k,i), P_inc(:,k))
       ENDDO

       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)

          DO l = 1, vv_mesh%gauss%l_G
             index  = index +1
             dw_loc = vv_mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
             rr_gauss(1,index) = ray
             rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))

             !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
             fs(1) = SUM(ff(j_loc,1) * vv_mesh%gauss%ww(:,l))
             ft(1) = SUM(V1m(j_loc,1,i) * vv_mesh%gauss%ww(:,l))
             fp(index,1,i) = -SUM(P(j_loc,1)*dw_loc(1,:))
             fp_inc(1) = -SUM(P_inc(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
             fs(2) = SUM(ff(j_loc,2) * vv_mesh%gauss%ww(:,l))
             ft(2) = SUM(V1m(j_loc,2,i) * vv_mesh%gauss%ww(:,l))
             fp(index,2,i) = -SUM(P(j_loc,2)*dw_loc(1,:))
             fp_inc(2) = -SUM(P_inc(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
             fs(3) = SUM(ff(j_loc,3) * vv_mesh%gauss%ww(:,l))
             ft(3) = SUM(V1m(j_loc,3,i) * vv_mesh%gauss%ww(:,l))
             fp(index,3,i) = -SUM(P(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             fp_inc(3) = -SUM(P_inc(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
             fs(4) = SUM(ff(j_loc,4) * vv_mesh%gauss%ww(:,l))
             ft(4) = SUM(V1m(j_loc,4,i) * vv_mesh%gauss%ww(:,l))
             fp(index,4,i) = SUM(P(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             fp_inc(4) = SUM(P_inc(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
             fs(5) = SUM(ff(j_loc,5) * vv_mesh%gauss%ww(:,l))
             ft(5) = SUM(V1m(j_loc,5,i) * vv_mesh%gauss%ww(:,l))
             fp(index,5,i) = -SUM(P(j_loc,1)*dw_loc(2,:))
             fp_inc(5) = -SUM(P_inc(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
             fs(6) = SUM(ff(j_loc,6) * vv_mesh%gauss%ww(:,l))
             ft(6) = SUM(V1m(j_loc,6,i) * vv_mesh%gauss%ww(:,l))
             fp(index,6,i) = -SUM(P(j_loc,2)*dw_loc(2,:))
             fp_inc(6) = -SUM(P_inc(j_loc,2)*dw_loc(2,:))

             rhs_gauss(index,:,i) =  (ft+fp_inc+fs-rotv_v(index,:,i))

          ENDDO
       ENDDO
       IF (inputs%if_ns_penalty) THEN
          IF(inputs%if_impose_vel_in_solids) THEN
             IF (list_mode(i)==0) THEN
                !Velocity imposed by penalty in solids (using BDF2).
                imposed_vel(:,:) = 3.d0*imposed_velocity_by_penalty(vv_mesh%rr,time)/(2*inputs%dt)
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
             ENDIF
          END IF
       END IF
    END DO

    IF (inputs%if_ns_penalty) THEN
       CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       bloc_size = SIZE(rhs_gauss,1)/nb_procs+1
       CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator, penal_in_real_space, vv_mesh, &
            rhs_gauss, rhs_gauss_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)

       rhs_gauss = rhs_gauss_penal

       IF(inputs%if_impose_vel_in_solids) THEN
          DO i = 1, SIZE(list_mode)
             IF (list_mode(i)==0) THEN
                rhs_gauss(:,:,i) = rhs_gauss_penal(:,:,i) + imposed_vel_gauss(:,:)
             END IF
          END DO
       END IF
    END IF

    rhs_gauss = rhs_gauss + fp

  END SUBROUTINE rhs_ns_gauss_3x3

  SUBROUTINE rhs_residual_ns_gauss_3x3(vv_mesh, pp_mesh, communicator, list_mode, time, du_dt,&
       pn, rotv_v, rhs_gauss, density, tempn, concn)
    !=================================
    !RHS for Residual of Navier-Stokes
    USE def_type_mesh
    USE my_util
    USE input_data
    USE fem_rhs_axi
    USE sft_parallele
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type)                                        :: vv_mesh, pp_mesh
    INTEGER,      DIMENSION(:),                 INTENT(IN) :: list_mode
    REAL(KIND=8),                               INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: du_dt
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: pn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: rotv_v
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)                                   :: fs
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8), DIMENSION(vv_mesh%np,6)                        :: ff
    REAL(KIND=8), DIMENSION(vv_mesh%np,2)                        :: P
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(vv_mesh%dom_me*vv_mesh%gauss%l_G,6,SIZE(list_mode)) :: rhs_gauss_penal, fp, ft
    REAL(KIND=8), DIMENSION(2,vv_mesh%gauss%l_G*vv_mesh%dom_me)                 :: rr_gauss
    REAL(KIND=8) :: ray
    INTEGER      :: m, l , i, k, index
    INTEGER      :: nb_procs, m_max_pad, bloc_size
#include "petsc/finclude/petsc.h"
    PetscErrorCode                   :: ierr
    MPI_Comm                         :: communicator

    DO i = 1, SIZE(list_mode)
       DO k = 1, 6
          ff(:,k) = source_in_NS_momentum(k, vv_mesh%rr, list_mode(i), i, time,  inputs%Re, 'ns', &
               density, tempn, concn)
       END DO
       DO k = 1, 2
          !===JLG+HF Dec 10 2019
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
          !CALL inject(pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
       ENDDO

       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)

          DO l = 1, vv_mesh%gauss%l_G
             index  = index +1
             dw_loc = vv_mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
             rr_gauss(1,index) = ray
             rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))

             !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
             fs(1) = SUM(ff(j_loc,1) * vv_mesh%gauss%ww(:,l))
             ft(index,1,i) = SUM(du_dt(j_loc,1,i) * vv_mesh%gauss%ww(:,l))
             fp(index,1,i)    = SUM(P(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
             fs(2) = SUM(ff(j_loc,2) * vv_mesh%gauss%ww(:,l))
             ft(index,2,i) = SUM(du_dt(j_loc,2,i) * vv_mesh%gauss%ww(:,l))
             fp(index,2,i)    = SUM(P(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
             fs(3) = SUM(ff(j_loc,3) * vv_mesh%gauss%ww(:,l))
             ft(index,3,i) = SUM(du_dt(j_loc,3,i) * vv_mesh%gauss%ww(:,l))
             fp(index,3,i)    = SUM(P(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
             fs(4) = SUM(ff(j_loc,4) * vv_mesh%gauss%ww(:,l))
             ft(index,4,i) = SUM(du_dt(j_loc,4,i) * vv_mesh%gauss%ww(:,l))
             fp(index,4,i)    = -SUM(P(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
             fs(5) = SUM(ff(j_loc,5) * vv_mesh%gauss%ww(:,l))
             ft(index,5,i) = SUM(du_dt(j_loc,5,i) * vv_mesh%gauss%ww(:,l))
             fp(index,5,i)    = SUM(P(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
             fs(6) = SUM(ff(j_loc,6) * vv_mesh%gauss%ww(:,l))
             ft(index,6,i) = SUM(du_dt(j_loc,6,i) * vv_mesh%gauss%ww(:,l))
             fp(index,6,i)    = SUM(P(j_loc,2)*dw_loc(2,:))

             rhs_gauss(index,:,i) =  -fs+rotv_v(index,:,i)
          ENDDO
       ENDDO
    END DO

    IF (inputs%if_ns_penalty) THEN
       CALL MPI_COMM_SIZE(communicator, nb_procs, ierr)
       m_max_pad = 3*SIZE(list_mode)*nb_procs/2
       bloc_size = SIZE(rhs_gauss,1)/nb_procs+1
       CALL FFT_PAR_VAR_ETA_PROD_GAUSS_DCL(communicator, penal_in_real_space, vv_mesh, &
            rhs_gauss, rhs_gauss_penal, nb_procs, bloc_size, m_max_pad, rr_gauss, time)
       rhs_gauss = rhs_gauss_penal
    END IF

    !===Add term time derivative and grad pressure
    rhs_gauss = rhs_gauss + ft + fp

  END SUBROUTINE rhs_residual_ns_gauss_3x3

  SUBROUTINE rhs_residual_ns_gauss_3x3_mom(vv_mesh, pp_mesh, list_mode, time, du_dt, pn, &
       density, rotb_b, rhs_gauss, tempn, concn)
    !=================================
    !RHS for Navier-Stokes
    USE def_type_mesh
    USE my_util
    USE input_data
    USE fem_rhs_axi
    USE sft_parallele
    USE boundary
    IMPLICIT NONE
    TYPE(mesh_type)                                        :: vv_mesh, pp_mesh
    INTEGER,      DIMENSION(:),                 INTENT(IN) :: list_mode
    REAL(KIND=8),                               INTENT(IN) :: time
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: du_dt
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: pn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: rotb_b
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)                                   :: fs, ft, fp
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8), DIMENSION(vv_mesh%np,6)                        :: ff
    REAL(KIND=8), DIMENSION(vv_mesh%np,2)                        :: P
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(2,vv_mesh%gauss%l_G*vv_mesh%dom_me)  :: rr_gauss
    REAL(KIND=8) :: ray
    INTEGER :: m, l , i, k, index

    DO i = 1, SIZE(list_mode)
       DO k = 1, 6
          ff(:,k) = source_in_NS_momentum(k, vv_mesh%rr, list_mode(i), i, time, inputs%Re, 'ns', &
               density, tempn, concn)
       END DO
       DO k = 1, 2
          !===JLG+HF Dec 10 2019
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
          !CALL inject(pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
       ENDDO

       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)

          DO l = 1, vv_mesh%gauss%l_G
             index  = index +1
             dw_loc = vv_mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
             rr_gauss(1,index) = ray
             rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))

             !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
             fs(1) = SUM(ff(j_loc,1) * vv_mesh%gauss%ww(:,l))
             ft(1) = SUM(du_dt(j_loc,1,i) * vv_mesh%gauss%ww(:,l))
             fp(1) = SUM(P(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
             fs(2) = SUM(ff(j_loc,2) * vv_mesh%gauss%ww(:,l))
             ft(2) = SUM(du_dt(j_loc,2,i) * vv_mesh%gauss%ww(:,l))
             fp(2) = SUM(P(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
             fs(3) = SUM(ff(j_loc,3) * vv_mesh%gauss%ww(:,l))
             ft(3) = SUM(du_dt(j_loc,3,i) * vv_mesh%gauss%ww(:,l))
             fp(3) = SUM(P(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
             fs(4) = SUM(ff(j_loc,4) * vv_mesh%gauss%ww(:,l))
             ft(4) = SUM(du_dt(j_loc,4,i) * vv_mesh%gauss%ww(:,l))
             fp(4) = -SUM(P(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
             fs(5) = SUM(ff(j_loc,5) * vv_mesh%gauss%ww(:,l))
             ft(5) = SUM(du_dt(j_loc,5,i) * vv_mesh%gauss%ww(:,l))
             fp(5) = SUM(P(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
             fs(6) = SUM(ff(j_loc,6) * vv_mesh%gauss%ww(:,l))
             ft(6) = SUM(du_dt(j_loc,6,i) * vv_mesh%gauss%ww(:,l))
             fp(6) = SUM(P(j_loc,2)*dw_loc(2,:))

             rhs_gauss(index,:,i) =  ft+fp-fs-rotb_b(index,:,i)
          ENDDO
       ENDDO
    END DO

  END SUBROUTINE rhs_residual_ns_gauss_3x3_mom

  SUBROUTINE rhs_ns_gauss_3x3_art_comp_mom(vv_mesh, pp_mesh, communicator, list_mode, time, V1m, pn, rotv_v, &
       rhs_gauss, tempn, concn, density, buoyancy)
    !=================================
    !RHS for Navier-Stokes
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
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: concn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: density
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: buoyancy
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%l_G*vv_mesh%dom_me,6,SIZE(list_mode)), INTENT(OUT) :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)                                   :: fs, ft
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: j_loc
    REAL(KIND=8), DIMENSION(vv_mesh%np,6)                        :: ff
    REAL(KIND=8), DIMENSION(vv_mesh%np,2)                        :: P
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(vv_mesh%dom_me*vv_mesh%gauss%l_G,6,SIZE(list_mode)) :: fp
    REAL(KIND=8), DIMENSION(2,vv_mesh%gauss%l_G*vv_mesh%dom_me)  :: rr_gauss
    REAL(KIND=8) :: ray
    INTEGER :: m, l , i, k, index
#include "petsc/finclude/petsc.h"
    PetscErrorCode                   :: ierr
    PetscMPIInt                      :: rank
    MPI_Comm                         :: communicator
    CALL MPI_Comm_rank(communicator,rank,ierr)
    DO i = 1, SIZE(list_mode)
       DO k = 1, 6
          ff(:,k) = source_in_NS_momentum(k, vv_mesh%rr, list_mode(i), i, time, inputs%Re, 'ns', &
               density, tempn, concn) + buoyancy(:,k,i)
       END DO
       DO k = 1, 2
          CALL inject_generic(inputs%type_fe_velocity, pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
          !CALL inject(pp_mesh%jj, vv_mesh%jj, pn(:,k,i), P(:,k))
       ENDDO

       index = 0
       DO m = 1, vv_mesh%dom_me
          j_loc = vv_mesh%jj(:,m)

          DO l = 1, vv_mesh%gauss%l_G
             index  = index +1
             dw_loc = vv_mesh%gauss%dw(:,:,l,m)

             !===Compute radius of Gauss point
             ray = SUM(vv_mesh%rr(1,j_loc)*vv_mesh%gauss%ww(:,l))
             rr_gauss(1,index) = ray
             rr_gauss(2,index) = SUM(vv_mesh%rr(2,j_loc)*vv_mesh%gauss%ww(:,l))

             !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
             fs(1) = SUM(ff(j_loc,1) * vv_mesh%gauss%ww(:,l))
             ft(1) = SUM(V1m(j_loc,1,i) * vv_mesh%gauss%ww(:,l))
             fp(index,1,i) = -SUM(P(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
             fs(2) = SUM(ff(j_loc,2) * vv_mesh%gauss%ww(:,l))
             ft(2) = SUM(V1m(j_loc,2,i) * vv_mesh%gauss%ww(:,l))
             fp(index,2,i) = -SUM(P(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
             fs(3) = SUM(ff(j_loc,3) * vv_mesh%gauss%ww(:,l))
             ft(3) = SUM(V1m(j_loc,3,i) * vv_mesh%gauss%ww(:,l))
             fp(index,3,i) = -SUM(P(j_loc,2)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
             fs(4) = SUM(ff(j_loc,4) * vv_mesh%gauss%ww(:,l))
             ft(4) = SUM(V1m(j_loc,4,i) * vv_mesh%gauss%ww(:,l))
             fp(index,4,i) = SUM(P(j_loc,1)*vv_mesh%gauss%ww(:,l))/ray*list_mode(i)
             !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
             fs(5) = SUM(ff(j_loc,5) * vv_mesh%gauss%ww(:,l))
             ft(5) = SUM(V1m(j_loc,5,i) * vv_mesh%gauss%ww(:,l))
             fp(index,5,i) = -SUM(P(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
             fs(6) = SUM(ff(j_loc,6) * vv_mesh%gauss%ww(:,l))
             ft(6) = SUM(V1m(j_loc,6,i) * vv_mesh%gauss%ww(:,l))
             fp(index,6,i) = -SUM(P(j_loc,2)*dw_loc(2,:))

             rhs_gauss(index,:,i) =  (ft+fs-rotv_v(index,:,i))

          ENDDO
       ENDDO
    END DO

    rhs_gauss = rhs_gauss + fp

  END SUBROUTINE rhs_ns_gauss_3x3_art_comp_mom

  !===Local subroutine
!!$  SUBROUTINE inject(jj_c, jj_f, pp_c, pp_f)
!!$    IMPLICIT NONE
!!$    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
!!$    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
!!$    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
!!$    REAL(KIND=8) :: half = 0.5
!!$    INTEGER:: m
!!$    IF (SIZE(jj_c,1)==3) THEN
!!$       DO m = 1, SIZE(jj_f,2)
!!$          pp_f(jj_f(1:3,m)) =  pp_c(jj_c(:,m))
!!$          pp_f(jj_f(4,m)) = (pp_c(jj_c(2,m)) + pp_c(jj_c(3,m)))*half
!!$          pp_f(jj_f(5,m)) = (pp_c(jj_c(3,m)) + pp_c(jj_c(1,m)))*half
!!$          pp_f(jj_f(6,m)) = (pp_c(jj_c(1,m)) + pp_c(jj_c(2,m)))*half
!!$       END DO
!!$
!!$    ELSE
!!$       DO m = 1, SIZE(jj_f,2)
!!$          pp_f(jj_f(1:4,m)) =  pp_c(jj_c(:,m))
!!$       END DO
!!$       pp_f(jj_f(5,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(4,:)))*half
!!$       pp_f(jj_f(6,:)) = (pp_c(jj_c(4,:)) + pp_c(jj_c(2,:)))*half
!!$       pp_f(jj_f(7,:)) = (pp_c(jj_c(2,:)) + pp_c(jj_c(3,:)))*half
!!$       pp_f(jj_f(8,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(4,:)))*half
!!$       pp_f(jj_f(9,:)) = (pp_c(jj_c(3,:)) + pp_c(jj_c(1,:)))*half
!!$       pp_f(jj_f(10,:)) = (pp_c(jj_c(1,:)) + pp_c(jj_c(2,:)))*half
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE inject
END MODULE rhs_gauss_computing
