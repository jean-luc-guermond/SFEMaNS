MODULE test_18
  IMPLICIT NONE
  !TEST 18
  REAL (KIND=8), PARAMETER, PUBLIC  :: lambda_mu_T18=10.0d0
CONTAINS
  FUNCTION mu_bar_in_fourier_space_anal_T18(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    USE def_type_mesh
    USE input_data
    USE my_util

    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER                                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL  :: pts
    INTEGER,DIMENSION(ne-nb+1),OPTIONAL         :: pts_ids
    INTEGER, DIMENSION(H_mesh%np)                :: global_ids
    INTEGER, DIMENSION(ne-nb+1)                  :: local_ids
    INTEGER                                      :: n,m
    REAL(KIND=8),DIMENSION(ne-nb+1)              :: r,z

    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
       local_ids=pts_ids
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
       DO m = 1, H_mesh%me
          global_ids(H_mesh%jj(:,m)) = H_mesh%i_d(m)
       END DO
       local_ids=global_ids(nb:ne)
    END IF

    DO n  = 1, ne - nb + 1
       IF(local_ids(n)==1) THEN
          vv(n) = 1.d0 + r(n) !mu1 , DCQ: If you change mu1_bar, you have to change
          !Jexact_gauss() as well
       ELSE
          vv(n) = 1.d0 + r(n) + 2*lambda_mu_T18*(1+r(n))/(z(n)**2*(3*r(n)+2))   !mu2
       END IF
    END DO
    RETURN
  END FUNCTION mu_bar_in_fourier_space_anal_T18

  !===Analytical grad_bar_mu_in_fourier_space (if needed)
  FUNCTION grad_mu_bar_in_fourier_space_anal_T18(pt,pt_id) RESULT(vv)
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2)           :: pt,vv
    INTEGER,DIMENSION(1)                :: pt_id
    REAL(KIND=8)                        :: r,z

    r=pt(1)
    z=pt(2)

    IF(pt_id(1)==1) THEN !grad_mu_1
       vv(1)= 1.d0
       vv(2)= 0.d0
    ELSE                 !grad_mu_2
       vv(1)=1.d0 +  ( (3*r+2)*(2*lambda_mu_T18)  -  ( (2*lambda_mu_T18*(1+r)))*(3)  )  /( z*(3*r+2) )**2
       vv(2)=  (2*lambda_mu_T18*(1+r))/(3*r+2)*(-2/z**3)
    END IF
    RETURN
  END FUNCTION grad_mu_bar_in_fourier_space_anal_T18

END MODULE test_18
