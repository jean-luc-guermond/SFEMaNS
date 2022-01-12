MODULE test_22
  IMPLICIT NONE
  !TEST 22
  REAL (KIND=8), PARAMETER, PUBLIC :: ratio_mu_T22 = 50.d0 ! the variation of mu
  REAL (KIND=8), PUBLIC            :: b_factor_T22 = (2**6) * (ratio_mu_T22-1.d0)/(ratio_mu_T22+1.d0)
  INTEGER,       PUBLIC            :: mode_mu_T22 = 4

CONTAINS

  FUNCTION f_test_T22(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: r, z
    REAL(KIND=8), DIMENSION(SIZE(r))       :: vv
    vv = b_factor_T22*(r*(1-r)*(z**2-1))**3
    RETURN
  END FUNCTION f_test_T22

  FUNCTION dfdr_test_T22(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: r, z
    REAL(KIND=8)            :: vv
    vv = 3 * b_factor_T22 * (z**2-1)**3 * (r*(1-r))**2 * (1-2*r)
    RETURN
  END FUNCTION dfdr_test_T22

  FUNCTION dfdz_test_T22(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: r, z
    REAL(KIND=8)            :: vv
    vv = 3*b_factor_T22*(r*(1-r))**3*(z**2-1)**2*(2*z)
    RETURN
  END FUNCTION dfdz_test_T22

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION mu_bar_in_fourier_space_anal_T22(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    USE def_type_mesh
    USE input_data
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER                                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1),  OPTIONAL :: pts_ids
    REAL(KIND=8),DIMENSION(ne-nb+1)            :: r,z

    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
    END IF

    vv=1.d0/(1.d0+abs(f_test_T22(r,z)))
    RETURN
  END FUNCTION mu_bar_in_fourier_space_anal_T22

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION grad_mu_bar_in_fourier_space_anal_T22(pt,pt_id) RESULT(vv)
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2)           :: pt,vv
    INTEGER,DIMENSION(1)                :: pt_id
    REAL(KIND=8),DIMENSION(1)           :: tmp,r,z
    REAL(KIND=8)                        :: sign
    INTEGER      :: n

    r(1)=pt(1)
    z(1)=pt(2)
    tmp=f_test_T22(r,z)
    IF (tmp(1) .GE. 0.d0 ) THEN
       sign =1.0
    ELSE
       sign =-1.0
    END IF

    vv(1)=-sign*dfdr_test_T22(r(1),z(1))/(1.d0 +abs(tmp(1)))**2
    vv(2)=-sign*dfdz_test_T22(r(1),z(1))/(1.d0 +abs(tmp(1)))**2
    RETURN

    !===Dummies variables to avoid warning
    n=pt_id(1)
    !===Dummies variables to avoid warning
  END FUNCTION grad_mu_bar_in_fourier_space_anal_T22

  FUNCTION mu_in_real_space_anal_T22(H_mesh,angles,nb_angles,nb,ne) RESULT(vv)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(:)                 :: angles
    INTEGER                                    :: nb_angles
    INTEGER                                    :: nb, ne
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
    INTEGER                                    :: ang

    DO ang = 1, nb_angles
       vv(ang,:) = 1/(1+f_test_T22(H_mesh%rr(1,nb:ne),H_mesh%rr(2,nb:ne))*COS(mode_mu_T22*angles(ang)))
    END DO
    RETURN
  END FUNCTION mu_in_real_space_anal_T22

END MODULE test_22
