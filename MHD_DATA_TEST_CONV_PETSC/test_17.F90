MODULE test_17
  IMPLICIT NONE
  !test 17
  REAL (KIND=8), PARAMETER, PRIVATE :: ratio_mu_T17 = 50.d0 ! the variation of mu
  REAL (KIND=8), PRIVATE            :: b_factor_T17_anal = (2**6) * (1.d0 - 1/ratio_mu_T17)

CONTAINS

  !====================================================
  !=======================Extra subroutines for test 17
  FUNCTION f_test_T17(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: r, z
    REAL(KIND=8), DIMENSION(SIZE(r))       :: vv
    vv = b_factor_t17_anal*(r*(1-r)*(z**2-1))**3
    RETURN
  END FUNCTION f_test_T17

  FUNCTION dfdr_test_T17(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: r, z
    REAL(KIND=8)            :: vv
    vv = 3 * b_factor_t17_anal * (z**2-1)**3 * (r*(1-r))**2 * (1-2*r)
    RETURN
  END FUNCTION dfdr_test_T17

  FUNCTION dfdz_test_T17(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN):: r, z
    REAL(KIND=8)            :: vv
    vv = 3*b_factor_t17_anal*(r*(1-r))**3*(z**2-1)**2*(2*z)
    RETURN
  END FUNCTION dfdz_test_T17

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION mu_bar_in_fourier_space_anal_T17(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
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
    vv=1.d0/(1.d0+f_test_T17(r,z))
    RETURN
  END FUNCTION mu_bar_in_fourier_space_anal_T17

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION grad_mu_bar_in_fourier_space_anal_T17(pt,pt_id) RESULT(vv)
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2)           :: vv
    REAL(KIND=8),DIMENSION(2)           :: pt
    INTEGER,DIMENSION(1)                :: pt_id
    REAL(KIND=8),DIMENSION(1)           :: r,z, tmp
    INTEGER      :: n

    r=pt(1)
    z=pt(2)
    tmp=(1.d0 +f_test_T17(r,z))**2
    vv(1)=-dfdr_test_T17(r(1),z(1))/tmp(1)
    vv(2)=-dfdz_test_T17(r(1),z(1))/tmp(1)
    RETURN

    !===Dummies variables to avoid warning
    n=pt_id(1)
    !===Dummies variables to avoid warning
  END FUNCTION grad_mu_bar_in_fourier_space_anal_T17

END MODULE test_17
