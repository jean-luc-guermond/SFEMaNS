MODULE test_23
  IMPLICIT NONE
  !TEST 23
  REAL (kind=8),PARAMETER, PUBLIC :: ratio_mu_T23 = 50.0d0 ! the variation of mu1
  REAL (kind=8),           PUBLIC :: b_factor_T23 = (1.d0/0.00016) * (ratio_mu_T23 - 1.d0)/(ratio_mu_T23 +1.d0)
  REAL (kind=8),           PUBLIC :: lambda_mu_T23 = 1.d0
  INTEGER,                 PUBLIC :: mode_mu_T23 = 3

CONTAINS

  FUNCTION s_test_T23(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)      :: r, z
    REAL(KIND=8)      :: vv
    vv = b_factor_T23*(  r*(r-1.d0)*(r-2.d0)*(z-0.25)*(z-1.d0)  )**3
    RETURN
  END FUNCTION s_test_T23

  FUNCTION Ds_test_T23(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)      :: r, z
    REAL(KIND=8),DIMENSION(2)      :: vv
    vv(1) = b_factor_T23*((z-0.25)*(z-1.d0))**3 * &
         (  3*( r*(r-1.d0)*(r-2.d0) )**2*( r*(  (r-1)+(r-2) ) +  (r-1)*(r-2)    ) )

    vv (2) = b_factor_T23*(  r*(r-1.d0)*(r-2.d0))**3 * &
         ( 3* (z-0.25)*(z-1.d0)  )**2 *( (z-1.d0) + (z-0.25) )
    RETURN
  END FUNCTION Ds_test_T23

  FUNCTION mu_bar_in_fourier_space_anal_T23(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    USE def_type_mesh
    USE input_data
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER                                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,DIMENSION(ne-nb+1),OPTIONAL        :: pts_ids
    INTEGER, DIMENSION(H_mesh%np)              :: global_ids
    INTEGER, DIMENSION(ne-nb+1)                :: local_ids
    INTEGER                                    :: n,m
    REAL(KIND=8),DIMENSION(ne-nb+1)            :: r,z
    REAL(KIND=8)                               :: s

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
       s   = s_test_T23(r(n),z(n))
       IF(local_ids(n)==1) THEN
          vv(n) = 1.d0/(1.d0 + ABS(s))  !mu1_bar, DCQ: If you change mu1_bar, you have to change
          !Jexact_gauss() as well
       ELSE
          vv(n) = 1.d0/(  (1.d0 + ABS(s)))*(1.d0 + lambda_mu_T23/z(n))   !mu2_bar
       END IF
    END DO
    RETURN
  END FUNCTION mu_bar_in_fourier_space_anal_T23

  !===Analytical grad_mu_in_fourier_space (if needed)
  FUNCTION grad_mu_bar_in_fourier_space_anal_T23(pt,pt_id) RESULT(vv)
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2)           :: vv
    REAL(KIND=8),DIMENSION(2)           :: pt
    INTEGER,DIMENSION(1)                :: pt_id
    REAL(KIND=8)                        :: r,z,sign,s
    REAL(KIND=8),DIMENSION(2)           :: tmp

    r=pt(1)
    z=pt(2)

    s=s_test_T23(r,z)
    IF (s .GE. 0.d0 ) THEN
       sign =1.0
    ELSE
       sign =-1.0
    END IF

    tmp=Ds_test_T23(r,z)!derivative

    IF(pt_id(1)==1) THEN !grad_mu_1
       vv(1)=-sign*tmp(1)/(1.d0 +abs(s))**2
       vv(2)=-sign*tmp(2)/(1.d0 +abs(s))**2
    ELSE                 !grad_mu_2
       vv(1)=-sign*tmp(1)/(1.d0 +abs(s))**2*(1+lambda_mu_T23/z)
       vv(2)=-sign*tmp(2)/(1.d0 +abs(s))**2*(1+lambda_mu_T23/z) + 1.d0/(1.d0+abs(s))*(-lambda_mu_T23/z**2)
    END IF
    RETURN
  END FUNCTION grad_mu_bar_in_fourier_space_anal_T23

  FUNCTION mu_in_real_space_anal_T23(H_mesh,angles,nb_angles,nb,ne) RESULT(vv)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(:)                 :: angles
    INTEGER                                    :: nb_angles
    INTEGER                                    :: nb, ne
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
    INTEGER                                    :: ang, n_loc, m, n
    REAL(KIND=8)                               :: tmp
    INTEGER, DIMENSION(H_mesh%np)              :: id
    DO m = 1, H_mesh%me
       id(H_mesh%jj(:,m)) = H_mesh%i_d(m)  !DCQ: Speed Efficient but requires more memory
    END DO
    DO n = nb, ne
       n_loc = n - nb + 1
       tmp   = s_test_T23(H_mesh%rr(1,n),H_mesh%rr(2,n))
       DO ang = 1, nb_angles

          IF (id(n)==1) THEN
             vv(ang,n_loc)  = 1.d0/(1.d0 + tmp*COS(mode_mu_T23*angles(ang)) )!mu_1
          ELSE
             vv(ang,n_loc)  = 1.d0/( 1.d0 + tmp*COS(mode_mu_T23*angles(ang)) ) &
                  *(  1.d0 + lambda_mu_T23/(H_mesh%rr(2,n))  )  !mu_2
          ENDIF
       END DO
    END DO
  END FUNCTION mu_in_real_space_anal_T23

END MODULE test_23
