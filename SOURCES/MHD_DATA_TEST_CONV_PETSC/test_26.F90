MODULE test_26
  IMPLICIT NONE
  !test 26
  REAL(KIND=8),  PARAMETER:: mu_disk_T26 = 5.d0
  REAL(KIND=8),  PARAMETER:: alpha_T26 = 1.d0
  INTEGER,  PARAMETER     :: test_mode_T26 = 4;
  REAL(KIND=8),  PARAMETER:: wjump_T26 = 0.032d0*(1.0d0)
  REAL(KIND=8),  PARAMETER:: zm_T26 = 0.d0, z1_T26 = zm_T26-wjump_T26
CONTAINS

  !====================================================
  !=======================Extra subroutines for test 26

  FUNCTION smooth_jump_down_T26(x,x0,x1) RESULT(vv)
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: a0,a1,a2,a3
    !This function is 1 at x0,
    !This function is 0 at x1,
    !Its derivative is 0  at (x1+x0)/2,
    !Cubic (Factorized)

    a0 = x1**2*(3*x0-x1)/(x0-x1)**3;
    a1 = -6.0*x0*x1/(x0-x1)**3;
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;

    vv = a0+a1*x+a2*x*x + a3*x*x*x
    RETURN
  END FUNCTION smooth_jump_down_T26

  !derivative with respect to x
  FUNCTION Dsmooth_jump_down_T26(x,x0,x1) RESULT(vv)
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    REAL(KIND=8)                               :: a0,a1,a2,a3
    !This function is 1 at x0,
    !This function is 0 at x1,
    !Its derivative is 0  at (x1+x0)/2,
    !Cubic Factorized

    a0 = x1**2*(3*x0-x1)/(x0-x1)**3;
    a1 = -6.0*x0*x1/(x0-x1)**3;
    a2 = (3.0*(x0+x1))/(x0-x1)**3;
    a3 = -2.0/(x0-x1)**3;

    vv = a1+2.d0*a2*x + 3.d0*a3*x*x
    RETURN
  END FUNCTION Dsmooth_jump_down_T26

  FUNCTION smooth_jump_up_T26(x,x0,x1) RESULT(vv)
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    !This function is 0 at x0,
    !This function is 1 at x1,
    !Its derivative is 0  at (x1+x0)/2,

    vv = 1.d0 - smooth_jump_down_T26( x,x0,x1 );
    RETURN
  END FUNCTION smooth_jump_up_T26

  !derivative with respect to x
  FUNCTION Dsmooth_jump_up_T26(x,x0,x1) RESULT(vv)
    REAL(KIND=8)                               :: x,x0,x1
    REAL(KIND=8)                               :: vv
    !This function is 0 at x0,
    !This function is 1 at x1,
    !Its derivative is 0  at (x1+x0)/2,

    vv =  - Dsmooth_jump_down_T26( x,x0,x1 );
    RETURN
  END FUNCTION Dsmooth_jump_up_T26

  FUNCTION mu_func_T26(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)       :: r,z
    REAL(KIND=8)       :: vv
    REAL(KIND=8)       :: Fz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! A smooth jump in z

    Fz=smooth_jump_up_T26(z,zm_T26,z1_T26)
    IF ( z.GE.zm_T26 ) THEN
       vv  = 1.d0
    ELSE  IF ( z.LE. z1_T26 ) THEN
       vv  = mu_disk_T26
    ELSE
       vv = Fz*(mu_disk_T26-1.d0) + 1.d0
    END IF
    RETURN

    !===Dummies variables to avoid warning
    Fz=r
    !===Dummies variables to avoid warning
  END FUNCTION mu_func_T26

  FUNCTION Dmu_func_T26(r,z) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8)       :: r,z
    REAL(KIND=8),DIMENSION(2)       :: vv
    REAL(KIND=8)       :: DFz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! A smooth jump in z

    DFz=Dsmooth_jump_up_T26(z,zm_T26,z1_T26)
    IF ( z.GE.zm_T26 ) THEN
       vv  = 0.d0
    ELSE  IF ( z.LE. z1_T26 ) THEN
       vv  = 0.d0
    ELSE
       vv(1)  = 0
       vv(2)  = DFz*(mu_disk_T26-1.d0)
    END IF
    RETURN

    !===Dummies variables to avoid warning
    DFz=r
    !===Dummies variables to avoid warning
  END FUNCTION Dmu_func_T26

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION mu_bar_in_fourier_space_anal_T26(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    USE def_type_mesh
    USE input_data
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER                                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,DIMENSION(ne-nb+1),OPTIONAL        :: pts_ids
    REAL(KIND=8),DIMENSION(ne-nb+1)            :: r,z
    INTEGER                                    :: n

    IF( PRESENT(pts) .AND. PRESENT(pts_ids) ) THEN !Computing mu at  pts
       r=pts(1,nb:ne)
       z=pts(2,nb:ne)
    ELSE
       r=H_mesh%rr(1,nb:ne) !Computing mu  at nodes
       z=H_mesh%rr(2,nb:ne)
    END IF

    DO n = 1, ne - nb + 1
       vv(n) = mu_func_T26(r(n),z(n))
    END DO
    RETURN
  END FUNCTION mu_bar_in_fourier_space_anal_T26

  !===Analytical mu_in_fourier_space (if needed)
  FUNCTION grad_mu_bar_in_fourier_space_anal_T26(pt,pt_id) RESULT(vv)
    USE input_data
    USE my_util
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2)           :: vv
    REAL(KIND=8),DIMENSION(2)           :: pt
    INTEGER,DIMENSION(1)                :: pt_id
    REAL(KIND=8)                        :: r,z
    REAL(KIND=8),DIMENSION(2)           :: tmp
    INTEGER      :: n

    r=pt(1)
    z=pt(2)
    tmp=Dmu_func_T26(r,z)
    vv(1)=tmp(1)
    vv(2)=tmp(2)
    RETURN

    !===Dummies variables to avoid warning
    n=pt_id(1)
    !===Dummies variables to avoid warning
  END FUNCTION grad_mu_bar_in_fourier_space_anal_T26

END MODULE test_26
