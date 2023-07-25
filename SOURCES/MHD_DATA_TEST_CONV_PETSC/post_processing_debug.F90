MODULE post_processing_debug

  PUBLIC :: compute_error, post_proc_test, regression

  PRIVATE
CONTAINS
  !---------------------------------------------------------------------------

  SUBROUTINE regression(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
       un, pn, Hn, Bn, phin, tempn, level_setn, concn, mu_H_field, &
       time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)
    USE boundary
    USE def_type_mesh
    USE input_data
    USE my_util
    USE tn_axi
    USE subroutine_ns_with_u
    USE sft_parallele
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                    :: conc_mesh
    TYPE(mesh_type), POINTER                    :: pp_mesh, vv_mesh
    TYPE(mesh_type), POINTER                    :: temp_mesh
    TYPE(mesh_type), POINTER                    :: H_mesh, phi_mesh
    INTEGER,      POINTER,  DIMENSION(:)        :: list_mode
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)    :: concn, un, pn, Hn, Bn, phin, tempn
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:)  :: level_setn
    REAL(KIND=8), POINTER,  DIMENSION(:)        :: mu_H_field
    REAL(KIND=8)                                :: time
    INTEGER                                     :: m_max_c, ierr, nd
    REAL(KIND=8) :: error, error_ref, error_cumul, norm
    INTEGER :: error_out
    MPI_Comm, DIMENSION(:), POINTER         :: comm_one_d, comm_one_d_ns
    MPI_Comm, DIMENSION(:), POINTER         :: comm_one_d_temp, comm_one_d_conc

    !inputs%numero_du_test_debug = 42
    !CALL post_proc_test(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
    !   un, pn, Hn, Bn, phin, tempn, level_setn, mu_H_field, &
    !   time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)

    error_cumul=0.d0
    norm=0.d0

    OPEN(UNIT = 21, FILE =  'regression_reference', FORM = 'formatted', STATUS = 'unknown')
    OPEN(UNIT = 22, FILE =  'current_regression_reference', FORM = 'formatted', STATUS = 'unknown')

    IF (conc_mesh%np/=0) THEN
       error=norm_SF(comm_one_d_conc, 'L2', conc_mesh, list_mode, concn)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)
    END IF

    IF (vv_mesh%np/=0) THEN
       error=norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)

       error=norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)

       IF (inputs%if_level_set) THEN
          IF (inputs%if_level_set_P2) THEN
             error=norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode,level_setn(1,:,:,:))
          ELSE
             error=norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode,level_setn(1,:,:,:))
          END IF
          WRITE(22,*) error
          READ(21,*) error_ref
          error_cumul = error_cumul + ABS(error-error_ref)
          norm = norm + ABS(error_ref)
       END IF
    END IF

    IF (temp_mesh%np/=0) THEN
       error=norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)
    END IF

    IF (H_mesh%np/=0) THEN
       error=norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)
    END IF

    IF (phi_mesh%np/=0) THEN
       error=norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phin)
       WRITE(22,*) error
       READ(21,*) error_ref
       error_cumul = error_cumul + ABS(error-error_ref)
       norm = norm + ABS(error_ref)
    END IF


    IF (isnan(error_cumul)) THEN
       error_out = 1 !===Test Failed
    ELSE IF (error_cumul/norm .GT. 1.d-7) THEN
       error_out = 2 !===Test Failed
    ELSE
       error_out= 123 !===Test Passed
    END IF

    CLOSE(21)
    CLOSE(22)
    CALL PetscFinalize(ierr)
    CALL EXIT(error_out)

    !===Dummy variables to avoid warning when post_proc_test call commented
    nd=SIZE(Bn,1); nd=m_max_c; nd=SIZE(mu_H_field); nd=FLOOR(time)
    !===Dummy variables to avoid warning when post_proc_test call commented
  END SUBROUTINE regression

  SUBROUTINE post_proc_test(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
       un, pn, Hn, Bn, phin, tempn, level_setn, concn, mu_H_field, &
       time, m_max_c, comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc)
    USE boundary
    USE def_type_mesh
    USE input_data
    USE my_util
    USE tn_axi
    USE subroutine_ns_with_u
    USE sft_parallele
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                    :: conc_mesh
    TYPE(mesh_type), POINTER                    :: pp_mesh, vv_mesh
    TYPE(mesh_type), POINTER                    :: temp_mesh
    TYPE(mesh_type), POINTER                    :: H_mesh, phi_mesh
    INTEGER,      POINTER,  DIMENSION(:)        :: list_mode
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)    :: un, pn, Hn, Bn, phin, tempn, concn
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:)  :: level_setn
    REAL(KIND=8), POINTER,  DIMENSION(:)        :: mu_H_field
    REAL(KIND=8)                                :: time
    INTEGER                                     :: m_max_c
    REAL(KIND=8), DIMENSION(SIZE(un,1),  SIZE(un,2),  SIZE(un,3))   :: un_m1, un_ex, un_error
    REAL(KIND=8), DIMENSION(SIZE(pn,1),  SIZE(pn,2),  SIZE(pn,3))   :: pn_m1, pn_ex, pn_error
    REAL(KIND=8), DIMENSION(SIZE(Hn,1),  SIZE(Hn,2),  SIZE(Hn,3))   :: Hn1, Hn_ex, Hn_error
    REAL(KIND=8), DIMENSION(SIZE(phin,1),SIZE(phin,2),SIZE(phin,3)) :: phin1
    REAL(KIND=8), DIMENSION(SIZE(tempn,1),  SIZE(tempn,2),  SIZE(tempn,3))   :: tempn_m1, tempn_ex, tempn_error
    REAL(KIND=8), DIMENSION(SIZE(concn,1),  SIZE(concn,2),  SIZE(concn,3))   :: concn_ex, concn_error
    REAL(KIND=8), DIMENSION(SIZE(level_setn,1),SIZE(level_setn,2),SIZE(level_setn,3),SIZE(level_setn,4)):: level_setn_m1

    INTEGER :: i, k, code, rank, int_nb, n, TYPE
    REAL(KIND=8), DIMENSION(4) :: norm_err_loc, norm_err
    REAL(KIND=8) :: err, norm
    REAL(KIND=8) :: moyenne
    MPI_Comm, DIMENSION(:), POINTER         :: comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)

    SELECT CASE(inputs%numero_du_test_debug)
    CASE(1,2,8,9,19,20,24,25,39,40,41)
       DO i = 1, m_max_c
          DO k= 1, 6
             un_m1(:,k,i) = un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
          END DO
          DO k= 1, 2
             pn_m1(:,k,i) = pn(:,k,i) - pp_exact(k,pp_mesh%rr,list_mode(i),time)
          END DO
          IF (list_mode(i) == 0)  THEN
             CALL Moy(comm_one_d(1),pp_mesh, pn_m1(:,1,i),moyenne)
             pn_m1(:,1,i) = pn_m1(:,1,i) - moyenne
          ENDIF
       END DO

       !norm_err(1) = norm_SF(comm_one_d_NS, 'L2', vv_mesh, list_mode, un_m1)
       norm_err(1) = SQRT(dot_product_SF(comm_one_d_NS,vv_mesh, list_mode, un_m1, un_m1))
       norm_err(2) = norm_SF(comm_one_d_NS, 'sH1', vv_mesh, list_mode, un_m1)
       norm_err(3) = norm_SF(comm_one_d_NS, 'div', vv_mesh, list_mode, un)
       norm_err(4) = norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn_m1)
       IF (rank==0) THEN
          WRITE(10,*) 'Velocity field   #####################'
          WRITE(10,*) 'L2 error on velocity  = ', norm_err(1)
          WRITE(10,*) 'H1 error on velocity  = ', norm_err(2)
          WRITE(10,*) 'L2 norm of divergence = ', norm_err(3)
          WRITE(10,*) 'Pressure field   #####################'
          WRITE(10,*) 'L2 error on pressure  = ', norm_err(4)
       END IF

       IF (inputs%if_temperature) THEN
          DO i = 1, m_max_c
             DO k= 1, 2
                tempn_m1(:,k,i) = tempn(:,k,i) - temperature_exact(k,temp_mesh%rr, list_mode(i), time)
             END DO
          END DO
          norm_err(1) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_m1)
          norm_err(2) = norm_SF(comm_one_d_temp, 'sH1', temp_mesh, list_mode, tempn_m1)
          IF (rank==0) THEN
             WRITE(10,*) 'Temperature field#####################'
             WRITE(10,*) 'L2 error on temperature = ', norm_err(1)
             WRITE(10,*) 'H1 error on temperature = ', norm_err(2)
             WRITE(10,*) 'Temperature field#####################'
          END IF
       END IF

       IF (inputs%if_level_set) THEN
          IF (inputs%if_level_set_P2) THEN
             DO i = 1, m_max_c
                DO k = 1, 2
                   DO int_nb = 1, inputs%nb_fluid - 1
                      level_setn_m1(int_nb,:,k,i) = level_setn(int_nb,:,k,i) &
                           -level_set_exact(int_nb,k,vv_mesh%rr,list_mode(i),time)
                   END DO
                END DO
             END DO
             IF (inputs%numero_du_test_debug==39) THEN
                norm_err(3) = norm_S_L1_zero_mode(comm_one_d_ns, vv_mesh, list_mode,level_setn_m1(1,:,:,:))
                IF (rank==0) THEN
                   WRITE(10,*) 'Level set field#####################'
                   WRITE(10,*) 'L1 error on level set', norm_err(3)
                   WRITE(10,*) 'Level set field#####################'
                END IF
             ELSE
                norm_err(3) = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode,level_setn_m1(1,:,:,:))
                IF (rank==0) THEN
                   WRITE(10,*) 'Level set field#####################'
                   WRITE(10,*) 'L2 error on level set', norm_err(3)
                   WRITE(10,*) 'Level set field#####################'
                END IF
             END IF
          ELSE
             DO i = 1, m_max_c
                DO k = 1, 2
                   DO int_nb = 1, inputs%nb_fluid - 1
                      level_setn_m1(int_nb,:,k,i) = level_setn(int_nb,:,k,i) &
                           -level_set_exact(int_nb,k,pp_mesh%rr,list_mode(i),time)
                   END DO
                END DO
             END DO
             IF (inputs%numero_du_test_debug==39) THEN
                norm_err(3) = norm_S_L1_zero_mode(comm_one_d_ns, pp_mesh, list_mode,level_setn_m1(1,:,:,:))
                IF (rank==0) THEN
                   WRITE(10,*) 'Level set field#####################'
                   WRITE(10,*) 'L1 error on level set', norm_err(3)
                   WRITE(10,*) 'Level set field#####################'
                END IF
             ELSE
                norm_err(3) = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode,level_setn_m1(1,:,:,:))
                IF (rank==0) THEN
                   WRITE(10,*) 'Level set field#####################'
                   WRITE(10,*) 'L2 error on level set', norm_err(3)
                   WRITE(10,*) 'Level set field#####################'
                END IF
             END IF
          END IF
       END IF

       IF (inputs%numero_du_test_debug==24) THEN
          DO i = 1, m_max_c
             DO k= 1, 2
                DO n = 1, SIZE(pp_mesh%rr,2)
                   IF (pp_mesh%rr(1,n).GT..5d0) THEN
                      pn_m1(n,k,i) = pn_m1(n,k,i)
                   ELSE
                      pn_m1(n,k,i) = 0.d0
                   END IF
                END DO
             END DO
          END DO

          norm_err(4) = norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn_m1)
          IF (rank==0) THEN
             WRITE(10,*) 'L2 error on pressure outter obstacle  = ', norm_err(4)
          END IF
       END IF


    CASE(3:7,10,12)
       DO i = 1, m_max_c
          DO k =1, 6
             Hn1(:,k,i) = Hn(:,k,i) - Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
          END DO
       END DO

       IF (inputs%nb_dom_phi/=0) THEN
          DO i = 1, m_max_c
             DO k =1, 2
                phin1(:,k,i) = phin(:,k,i) - Phiexact(k, phi_mesh%rr, list_mode(i), inputs%mu_phi, time)
             END DO
          END DO
       END IF

       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn1)
       norm_err(1) = err/norm

       norm = norm_SF(comm_one_d, 'sH1', H_mesh, list_mode, Hn)
       err  = norm_SF(comm_one_d, 'curl', H_mesh, list_mode, Hn1)
       norm_err(2) = err/norm

       norm = norm_SF(comm_one_d, 'sH1', H_mesh, list_mode, Bn)
       err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)

       norm_err(3) = err/norm

       IF (inputs%nb_dom_phi/=0) THEN
          norm = norm_SF(comm_one_d, 'sH1', phi_mesh, list_mode, phin)
          err  = norm_SF(comm_one_d, 'sH1', phi_mesh, list_mode, phin1)
          norm_err(4) = err/norm
       ELSE
          norm_err(4) = inputs%norm_ref(4)
       END IF

       IF (rank==0) THEN
          WRITE(rank+10,*) 'Magnetic field #####################'
          WRITE(rank+10,*) 'L2 norm of error on Hn          = ', norm_err(1)
          WRITE(rank+10,*) 'L2 norm of error on Curl(Hn)    = ', norm_err(2)
          WRITE(rank+10,*) 'L2 norm of Div(mu Hn)           = ', norm_err(3)
          IF (inputs%nb_dom_phi/=0) THEN
             WRITE(10,*) 'Scal potential #####################'
             WRITE(10,*) 'H1 norm of error on phin        = ',  norm_err(4)
          END IF
       END IF

    CASE(11)
       norm_err(1) = norm_SF(comm_one_d_NS, 'sH1', vv_mesh, list_mode, un)
       norm_err(2) = norm_SF(comm_one_d,    'L2', H_mesh, list_mode, Hn)
       norm_err(3) = norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn)
       norm_err(4) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn) ! MODIFICATION: comm_one_d_temp instead of ns
       IF (rank==0) THEN
          WRITE(10,*) '######################################'
          WRITE(10,*) 'H1 norm on velocity     = ', norm_err(1)
          WRITE(10,*) 'L2 norm of Hn           = ', norm_err(2)
          WRITE(10,*) 'L2 norm of pressure     = ', norm_err(3)
          WRITE(10,*) 'L2 norm of temperature  = ', norm_err(4)
          WRITE(10,*) '######################################'
       END IF

    CASE(13)
       norm_err(1) = norm_SF(comm_one_d_NS, 'sH1', vv_mesh, list_mode, un)
       norm_err(2) = norm_SF(comm_one_d,    'div', H_mesh,  list_mode, Hn)
       norm_err(3) = norm_SF(comm_one_d,    'L2',  H_mesh,  list_mode, Hn)
       norm_err(4) = norm_SF(comm_one_d_NS, 'L2',  pp_mesh, list_mode, pn)
       IF (rank==0) THEN
          WRITE(10,*) '######################################'
          WRITE(10,*) 'H1 norm on velocity     = ', norm_err(1)
          WRITE(10,*) 'L2 norm of div(Hn)      = ', norm_err(2)
          WRITE(10,*) 'L2 norm of Hn           = ', norm_err(3)
          WRITE(10,*) 'L2 norm of pressure     = ', norm_err(4)
          WRITE(10,*) '######################################'
       END IF

    CASE(14)
       IF (rank==0) THEN
          READ(10,*) norm_err(1)
          READ(10,*) norm_err(2)
          READ(11,*) norm_err(3)
          READ(11,*) norm_err(4)
       END IF

    CASE(15)

       norm_err_loc = 0.d0
       DO i = 1, m_max_c
          norm_err_loc(rank+1)= 0.5*(norme_L2_champ_par(comm_one_d_ns(1), &
               vv_mesh, list_mode(i:i), un(:,:,i:i)))**2
       END DO
       norm_err_loc(4)= norm_SF(comm_one_d_NS, 'div', vv_mesh, list_mode, un)/&
            norm_SF(comm_one_d_NS, 'sH1', vv_mesh, list_mode, un)
       CALL MPI_ALLREDUCE(norm_err_loc,norm_err,4,MPI_DOUBLE_PRECISION, &
            MPI_MAX, comm_one_d(2), code)

       IF (rank==0) THEN
          WRITE(10,*) '######################################'
          WRITE(10,*) 'L2 norm of mode 0       = ', norm_err(1)
          WRITE(10,*) 'L2 norm of mode 1       = ', norm_err(2)
          WRITE(10,*) 'L2 norm of mode 2       = ', norm_err(3)
          WRITE(10,*) 'Relative L2 norm of div.= ', norm_err(4)
          WRITE(10,*) '######################################'
       END IF

    CASE(16)
       norm_err(1) = 0.5*norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)**2
       CALL angular_momentum(vv_mesh, list_mode, un, norm_err(2:4))
       IF (rank==0) THEN
          WRITE(10,*) '######################################'
          WRITE(10,*) 'Total kinetic energy = ', norm_err(1)
          WRITE(10,*) 'Mx                   = ', norm_err(2)
          WRITE(10,*) 'My                   = ', norm_err(3)
          WRITE(10,*) 'Mz                   = ', norm_err(4)
          WRITE(10,*) '######################################'
       END IF

    CASE(17,22,27)
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
       err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
       norm_err(1) = err/norm
       DO i = 1, m_max_c
          DO k = 1, 6
             Hn1(:,k,i) = Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
             Hn1(:,k,i) = Hn1(:,k,i) - Hn(:,k,i)
          ENDDO
       ENDDO
       err = norm_SF(comm_one_d, 'curl', H_mesh, list_mode, Hn1)
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       norm_err(2) = err/norm
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn1)
       norm_err(3) = err/norm
       DO i = 1, m_max_c
          DO k = 1, 2
             phin1(:,k,i) = phin(:,k,i) - Phiexact(k, phi_mesh%rr, list_mode(i), inputs%mu_phi, time)
          END DO
       END DO
       err = norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phin1)
       norm = norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phin)
       norm_err(4) = err/norm
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm on div of B err/norm           =', norm_err(1)
          WRITE(10,*) ' L2-norm on curl of (H-Hexact) err/norm =', norm_err(2)
          WRITE(10,*) ' L2-norm of (H-Hexact) err/norm         =', norm_err(3)
          WRITE(10,*) ' L2-norm on phi err/norm                =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(18,23,26)
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
       err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
       norm_err(1) = err/norm
       DO i = 1, m_max_c
          DO k = 1, 6
             Hn1(:,k,i) = Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
             Hn1(:,k,i) = Hn1(:,k,i) - Hn(:,k,i)
          ENDDO
       ENDDO
       err = norm_SF(comm_one_d, 'curl', H_mesh, list_mode, Hn1)
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       norm_err(2) = err/norm
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn1)
       norm_err(3) = err/norm
       norm_err(4) = norm
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm on div of B err/norm           =', norm_err(1)
          WRITE(10,*) ' L2-norm on curl of (H-Hexact) err/norm =', norm_err(2)
          WRITE(10,*) ' L2-norm of (H-Hexact) err/norm         =', norm_err(3)
          WRITE(10,*) ' L2-norm of H                           =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(21)
       DO i = 1, m_max_c
          DO k= 1, 6
             un_m1(:,k,i) = un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
             Hn1(:,k,i) = Hn(:,k,i) - Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
          END DO
          DO k= 1, 2
             pn_m1(:,k,i) = pn(:,k,i) - pp_exact(k,pp_mesh%rr,list_mode(i),time)
             DO int_nb = 1, inputs%nb_fluid - 1
                level_setn_m1(int_nb,:,k,i) = level_setn(int_nb,:,k,i) &
                     - level_set_exact(int_nb,k,pp_mesh%rr,list_mode(i),time)
             END DO
          END DO
          IF (list_mode(i) == 0)  THEN
             CALL Moy(comm_one_d(1),pp_mesh, pn_m1(:,1,i),moyenne)
             pn_m1(:,1,i) = pn_m1(:,1,i) - moyenne
          ENDIF
       END DO
       norm_err(1) = norm_SF(comm_one_d_NS, 'L2', vv_mesh, list_mode, un_m1)
       norm_err(2) = norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn_m1)
       norm_err(3) = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, level_setn_m1(1,:,:,:))
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn1)
       norm_err(4) = err
       IF (rank==0) THEN
          WRITE(10,*) 'Velocity field   #####################'
          WRITE(10,*) 'L2 error on velocity  = ', norm_err(1)
          WRITE(10,*) 'Pressure field   #####################'
          WRITE(10,*) 'L2 error on pressure  = ', norm_err(2)
          WRITE(10,*) 'Level Set        #####################'
          WRITE(10,*) 'L2 error on level set = ', norm_err(3)
          WRITE(10,*) 'Magnetic field   #####################'
          WRITE(10,*) 'L2 error on Hn        = ', norm_err(4)
       END IF


    CASE(28)
       err = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)
       norm = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un)
       norm_err(1) = err
       norm_err(2) = err/norm
       norm = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un)
       norm_err(3) = norm
       norm = norm_SF(comm_one_d, 'H1', pp_mesh, list_mode, pn)
       norm_err(4) = norm

       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm on div of  u                    =', norm_err(1)
          WRITE(10,*) ' L2-norm of div of  u err/norm           =', norm_err(2)
          WRITE(10,*) ' L2-norm  of u                           =', norm_err(3)
          WRITE(10,*) ' H1-norm  of p                           =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(29)

       err = norm_SF(comm_one_d, 'div', vv_mesh, list_mode, un)
       norm = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un)
       norm_err(1) = err/norm
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Bn)
       err = norm_SF(comm_one_d, 'div', H_mesh, list_mode, Bn)
       norm_err(2) = err/norm
       norm_err(3) = norm
       err=  dot_product_SF(comm_one_d, H_mesh, list_mode, Hn, Bn)
       norm_err(4) = 0.5*err
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of div of  u err/norm          =', norm_err(1)
          WRITE(10,*) ' L2-norm on div of B err/norm           =', norm_err(2)
          WRITE(10,*) ' L2-norm of of  B                       =', norm_err(3)
          WRITE(10,*) ' integral of 0.5*B.H                    =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(30)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(2) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(3) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex)
       norm_err(4) = norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_ex)
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) ' L2-norm of error on p                      =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(3)
          WRITE(10,*) ' H1-norm of error on T / H1-norm of T exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(31,32)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(2) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error) / norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(3) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex) ! MODIFICATION: comm_one_d_temp instead of ns
       norm_err(4) = norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_ex) ! MODIFICATION: comm_one_d_temp instead of ns
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) ' L2-norm of error on p / L2-norm of p exact =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(3)
          WRITE(10,*) ' H1-norm of error on T / H1-norm of T exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(33)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(2) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(3) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex) ! MODIFICATION: comm_one_d_temp instead of ns
       DO TYPE=1,6
          DO i=1,size(list_mode)
             Hn_ex(:,TYPE,i) = Hexact(H_mesh,TYPE,H_mesh%rr,list_mode(i),mu_H_field,time)
          END DO
       END DO
       Hn_error = Hn - Hn_ex
       norm_err(4) = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_error) / norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_ex)
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) ' L2-norm of error on p                      =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(3)
          WRITE(10,*) ' L2-norm of error on H / L2-norm of H exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(34,37,38)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(2) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error) / norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(3) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex) ! MODIFICATION: comm_one_d_temp instead of ns
       DO TYPE=1,6
          DO i=1,size(list_mode)
             Hn_ex(:,TYPE,i) = Hexact(H_mesh,TYPE,H_mesh%rr,list_mode(i),mu_H_field,time)
          END DO
       END DO
       Hn_error = Hn - Hn_ex
       norm_err(4) = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_error) / norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_ex)
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) ' L2-norm of error on p / L2-norm of p exact =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(3)
          WRITE(10,*) ' L2-norm of error on H / L2-norm of H exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(35)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(2) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(3) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex)
       norm_err(4) = norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'H1', temp_mesh, list_mode, tempn_ex)
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u                      =', norm_err(1)
          WRITE(10,*) ' L2-norm of error on p                      =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(3)
          WRITE(10,*) ' H1-norm of error on T / H1-norm of T exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(36)

       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(1) = norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'L2', vv_mesh, list_mode, un_ex)
       DO TYPE=1,6
          DO i=1,size(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm_err(2) = norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un_error) / norm_SF(comm_one_d, 'H1', vv_mesh, list_mode, un_ex)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       pn_error = pn - pn_ex
       norm_err(3) = norm_SF(comm_one_d, 'L2', pp_mesh, list_mode, pn_error)
       DO TYPE=1,2
          DO i=1,size(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm_err(4) = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error) / &
            norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex) ! MODIFICATION: comm_one_d_temp instead of ns
       IF (rank==0) THEN
          WRITE(10,*) '########################################################'
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) ' H1-norm of error on u / H1-norm of u exact =', norm_err(2)
          WRITE(10,*) ' L2-norm of error on p                      =', norm_err(3)
          WRITE(10,*) ' L2-norm of error on T / L2-norm of T exact =', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE(43)
       DO i = 1, m_max_c
          DO k= 1, 6
             un_m1(:,k,i) = un(:,k,i) - vv_exact(k,vv_mesh%rr,list_mode(i),time)
          END DO
          DO k= 1, 2
             pn_m1(:,k,i) = pn(:,k,i) - pp_exact(k,pp_mesh%rr,list_mode(i),time)
          END DO
          IF (list_mode(i) == 0)  THEN
             CALL Moy(comm_one_d(1),pp_mesh, pn_m1(:,1,i),moyenne)
             pn_m1(:,1,i) = pn_m1(:,1,i) - moyenne
          ENDIF
       END DO
       norm_err(1) = SQRT(dot_product_SF(comm_one_d_NS,vv_mesh, list_mode, un_m1, un_m1))
       norm_err(2) = norm_SF(comm_one_d_NS, 'L2', pp_mesh, list_mode, pn_m1)

       DO i = 1, m_max_c
          DO k =1, 6
             Hn1(:,k,i) = Hn(:,k,i) - Hexact(H_mesh, k, H_mesh%rr, list_mode(i), mu_H_field, time)
          END DO
       END DO
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn)
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn1)
       norm_err(3) = err/norm

       DO TYPE=1,2
          DO i=1,size(list_mode)
             concn_ex(:,TYPE,i) = concentration_exact(TYPE,conc_mesh%rr,list_mode(i),time)
          END DO
       END DO
       concn_error = concn - concn_ex
       norm_err(4) = norm_SF(comm_one_d_conc, 'L2', conc_mesh, list_mode, concn_error) / &
            norm_SF(comm_one_d_conc, 'L2', conc_mesh, list_mode, concn_ex)

       IF (rank==0) THEN
          WRITE(10,*) ' L2-norm of error on u / L2-norm of u exact =', norm_err(1)
          WRITE(10,*) '########################################################'
          WRITE(10,*) 'L2 error on velocity  = ', norm_err(1)
          WRITE(10,*) 'L2 error on pressure  = ', norm_err(2)
          WRITE(10,*) 'L2 error on magnetic field = ', norm_err(3)
          WRITE(10,*) 'L2 error on concentration  = ', norm_err(4)
          WRITE(10,*) '########################################################'
       END IF

    CASE DEFAULT
       CALL error_petsc(' BUG in post_proc_test: We should not be here')

    END SELECT

    IF (rank==0) THEN
       ! Compare with references
!!$       IF (MAXVAL(ABS(inputs%norm_ref-norm_err)/inputs%norm_ref)<error_on_tests) THEN
       IF (MAXVAL(ABS(inputs%norm_ref-norm_err))<1.d-8) THEN
          WRITE(57,'(A,I2,A)') 'test #',inputs%numero_du_test_debug,' OK'
       ELSE
!!$          WRITE(57,'(A,I2,2x,e15.7)') 'Problem with test #', inputs%numero_du_test_debug, &
!!$               MAXVAL(ABS(inputs%norm_ref-norm_err)/inputs%norm_ref)
          WRITE(57,'(A,I2,2x,e15.7)') 'Problem with test #', inputs%numero_du_test_debug, &
               MAXVAL(ABS(inputs%norm_ref-norm_err))
       END IF
       ! End compare
    END IF

  END SUBROUTINE post_proc_test

  SUBROUTINE compute_error(conc_mesh, vv_mesh, pp_mesh, temp_mesh, H_mesh, phi_mesh, list_mode, &
       concn, un, pn, tempn, Hn, phin, level_setn, time, mu_H_field, comm_one_d_conc, &
       comm_one_d_ns, comm_one_d_temp, comm_one_d)
    USE boundary
    USE def_type_mesh
    USE input_data
    USE my_util
    USE tn_axi
    USE subroutine_ns_with_u
    USE sft_parallele
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type), POINTER                    :: conc_mesh
    TYPE(mesh_type), POINTER                    :: pp_mesh, vv_mesh
    TYPE(mesh_type), POINTER                    :: temp_mesh
    TYPE(mesh_type), POINTER                    :: H_mesh, phi_mesh
    INTEGER,      POINTER,  DIMENSION(:)        :: list_mode
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:)    :: concn, un, pn, tempn, Hn, phin
    REAL(KIND=8), POINTER,  DIMENSION(:,:,:,:)  :: level_setn
    REAL(KIND=8)                                :: time
    REAL(KIND=8), POINTER,  DIMENSION(:)        :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(un,1),  SIZE(un,2),  SIZE(un,3))   :: un_ex, un_error
    REAL(KIND=8), DIMENSION(SIZE(pn,1),  SIZE(pn,2),  SIZE(pn,3))   :: pn_ex, pn_error
    REAL(KIND=8), DIMENSION(SIZE(Hn,1),  SIZE(Hn,2),  SIZE(Hn,3))   :: Hn_ex, Hn_error
    REAL(KIND=8), DIMENSION(SIZE(phin,1),SIZE(phin,2),SIZE(phin,3)) :: phin_ex, phin_error
    REAL(KIND=8), DIMENSION(SIZE(tempn,1),  SIZE(tempn,2),  SIZE(tempn,3))   :: tempn_ex, tempn_error
    REAL(KIND=8), DIMENSION(SIZE(concn,1),  SIZE(concn,2),  SIZE(concn,3))   :: concn_ex, concn_error
    REAL(KIND=8), DIMENSION(SIZE(level_setn,2),SIZE(level_setn,3),SIZE(level_setn,4)):: level_setn_ex
    REAL(KIND=8), DIMENSION(SIZE(level_setn,2),SIZE(level_setn,3),SIZE(level_setn,4)):: level_setn_error
    INTEGER :: i, k, code, rank, int_nb, TYPE
    REAL(KIND=8) :: err, norm, err_H1, norm_H1
    REAL(KIND=8) :: moyenne
    MPI_Comm, DIMENSION(:), POINTER         :: comm_one_d, comm_one_d_ns, comm_one_d_temp, comm_one_d_conc

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,code)
    IF (rank==0) THEN
       WRITE(11,*) '########################################################'
       WRITE(11,*) 'Time = ', time
    END IF

    IF (inputs%if_concentration) THEN
       DO TYPE = 1, 2
          DO i = 1, SIZE(list_mode)
             concn_ex(:,TYPE,i) = concentration_exact(TYPE,conc_mesh%rr,list_mode(i),time)
          END DO
       END DO
       concn_error = concn - concn_ex
       norm = norm_SF(comm_one_d_conc, 'L2', temp_mesh, list_mode, concn_ex)
       err = norm_SF(comm_one_d_conc, 'L2', temp_mesh, list_mode, concn_error)       
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on concentration  = ', err, err/norm
       END IF
    END IF

    IF (vv_mesh%np>0) THEN
       !Velocity
       DO TYPE = 1, 6
          DO i = 1, SIZE(list_mode)
             un_ex(:,TYPE,i) = vv_exact(TYPE,vv_mesh%rr,list_mode(i),time)
          END DO
       END DO
       un_error = un - un_ex
       norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un_ex)
       err = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, un_error)
       norm_H1 = norm_SF(comm_one_d_ns, 'H1', vv_mesh, list_mode, un_ex)
       err_H1 = norm_SF(comm_one_d_ns, 'H1', vv_mesh, list_mode, un_error)
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on velocity  = ', err, err/norm
          WRITE(11,*) 'H1 error/relative error on velocity  = ', err_H1, err_H1/norm_H1
       END IF
       !Pressure
       DO TYPE = 1, 2
          DO i = 1, SIZE(list_mode)
             pn_ex(:,TYPE,i) = pp_exact(TYPE,pp_mesh%rr,list_mode(i),time)
             IF (list_mode(i) == 0)  THEN
                CALL Moy(comm_one_d(1),pp_mesh, pn_ex(:,1,i),moyenne)
                pn_ex(:,1,i) = pn_ex(:,1,i) - moyenne
             ENDIF
          END DO
       END DO
       pn_error = pn - pn_ex
       norm = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, pn_ex)
       err = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, pn_error)
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on pressure  = ', err, err/norm
       END IF
       !Level set
       IF (inputs%if_level_set) THEN
          IF (inputs%if_level_set_P2) THEN
             DO int_nb = 1, inputs%nb_fluid - 1
                DO TYPE = 1, 2
                   DO i = 1, SIZE(list_mode)
                      level_setn_ex(:,k,i) = &
                           level_set_exact(int_nb,k,vv_mesh%rr,list_mode(i),time)
                   END DO
                END DO
                level_setn_error=level_setn(int_nb,:,:,:) - level_setn_ex
                norm = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, level_setn_ex)
                err = norm_SF(comm_one_d_ns, 'L2', vv_mesh, list_mode, level_setn_error)
                IF (rank==0) THEN
                   WRITE(11,*) "Level set interface int_nb = ", int_nb
                   WRITE(11,*) 'L2 error/relative error on level set  = ', err, err/norm
                END IF
             END DO
          ELSE
             DO int_nb = 1, inputs%nb_fluid - 1
                DO TYPE = 1, 2
                   DO i = 1, SIZE(list_mode)
                      level_setn_ex(:,k,i) = &
                           level_set_exact(int_nb,k,pp_mesh%rr,list_mode(i),time)
                   END DO
                END DO
                level_setn_error=level_setn(int_nb,:,:,:) - level_setn_ex
                norm = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, level_setn_ex)
                err = norm_SF(comm_one_d_ns, 'L2', pp_mesh, list_mode, level_setn_error)
                IF (rank==0) THEN
                   WRITE(11,*) "Level set interface int_nb = ", int_nb
                   WRITE(11,*) 'L2 error/relative error on level set  = ', err, err/norm
                END IF
             END DO
          END IF
       END IF
    END IF

    IF (inputs%if_temperature) THEN
        DO TYPE = 1, 2
          DO i = 1, SIZE(list_mode)
             tempn_ex(:,TYPE,i) = temperature_exact(TYPE,temp_mesh%rr,list_mode(i),time)
          END DO
       END DO
       tempn_error = tempn - tempn_ex
       norm = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_ex)
       err = norm_SF(comm_one_d_temp, 'L2', temp_mesh, list_mode, tempn_error)
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on temperature  = ', err, err/norm
       END IF
    END IF

    IF (H_mesh%np>0) THEN
        DO TYPE = 1, 6
          DO i = 1, SIZE(list_mode)
             Hn_ex(:,TYPE,i) = Hexact(H_mesh,TYPE,H_mesh%rr,list_mode(i),mu_H_field,time)
          END DO
       END DO
       Hn_error = Hn - Hn_ex
       norm = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_ex)
       err = norm_SF(comm_one_d, 'L2', H_mesh, list_mode, Hn_error)
       norm_H1 = norm_SF(comm_one_d, 'H1', H_mesh, list_mode, Hn_ex)
       err_H1 = norm_SF(comm_one_d, 'H1', H_mesh, list_mode, Hn_error)
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on magnetic field = ', err, norm/err 
          WRITE(11,*) 'H1 error/relativeerror on magnetic field = ', err_H1, norm_H1/err_H1
       END IF
    END IF

    IF (phi_mesh%np>0) THEN
       DO TYPE = 1, 2
          DO i = 1, SIZE(list_mode)
             phin_ex(:,k,i) = Phiexact(TYPE,phi_mesh%rr,list_mode(i),inputs%mu_phi,time)
          END DO
       END DO
       phin_error = phin - phin_ex
       norm = norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phin_ex)
       err = norm_SF(comm_one_d, 'L2', phi_mesh, list_mode, phin_error)
       IF (rank==0) THEN
          WRITE(11,*) 'L2 error/relative error on magnetic field = ', err, norm/err 
       END IF
    END IF

    IF (rank==0) THEN
       WRITE(11,*) '########################################################'
    END IF

  END SUBROUTINE compute_error
  !---------------------------------------------------------------------------

END MODULE post_processing_debug
