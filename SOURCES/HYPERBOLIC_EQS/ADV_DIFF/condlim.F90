MODULE condlim
  PUBLIC :: uexact, vexact, source, block_source
  PRIVATE
  REAL(KIND=8), PUBLIC :: mu1=1.d0, mu2=1.d0
CONTAINS

  FUNCTION uexact(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: uu
    REAL(KIND=8) :: pi
    pi = ACOS(-1.d0)
    uu = COS(mu1*2*pi*rr(1,:))*COS(mu2*2*pi*rr(2,:))
  END FUNCTION uexact

  FUNCTION vexact(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: uu
    REAL(KIND=8) :: pi
    pi = ACOS(-1.d0)
    uu = COS(mu2*2*pi*rr(1,:))*COS(mu1*2*pi*rr(2,:))
  END FUNCTION vexact

  FUNCTION source(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rr
    REAL(KIND=8), DIMENSION(SIZE(rr,2)) :: uu
    REAL(KIND=8) :: pi
    pi = ACOS(-1.d0)
    uu = uexact(rr) + ((mu1*2*pi)**2+(mu2*2*pi)**2)*uexact(rr)
  END FUNCTION source

  FUNCTION block_source(rr) RESULT(uu)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:) :: rr
    REAL(KIND=8), DIMENSION(2,SIZE(rr,2)) :: uu
    REAL(KIND=8) :: pi
    pi = ACOS(-1.d0)
    uu(1,:) = ((mu1*2*pi)**2+(mu2*2*pi)**2)*uexact(rr) + uexact(rr) + vexact(rr)
    uu(2,:) = ((mu1*2*pi)**2+(mu2*2*pi)**2)*vexact(rr) + vexact(rr) + uexact(rr)
  END FUNCTION block_source

END MODULE condlim
