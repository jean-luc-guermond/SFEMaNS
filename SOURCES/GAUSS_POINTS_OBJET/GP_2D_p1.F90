!===
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!Author:  Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE GP_2d_p1
  PRIVATE
  PUBLIC element_2d_p1, element_2d_p1_boundary, element_1d_p1, element_1d_p1_at_nodes
CONTAINS
  SUBROUTINE element_2d_p1 (w, d, p, n_w, l_G)
    !===Triangular element with linear interpolation
    !===and three Gauss integration points
    !===w(n_w, l_G)    : values of shape functions at Gauss points
    !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
    !===p(l_G)         : weight for Gaussian quadrature at Gauss points
    IMPLICIT NONE
    INTEGER,                              INTENT(IN)  :: n_w, l_G
    REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p
    REAL(KIND=8), DIMENSION(l_G) :: xx, yy
    INTEGER :: j
    REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  six = 6
    REAL(KIND=8) :: f1, f2, f3, x, y
    f1(x, y) = one - x - y
    f2(x, y) = x
    f3(x, y) = y
    xx(1) = one/six;  xx(2) = two/three;  xx(3) = one/six
    yy(1) = one/six;  yy(2) = one/six;    yy(3) = two/three
    DO j = 1, l_G
       w(1, j) = f1(xx(j), yy(j))
       d(1, 1, j) = - one
       d(2, 1, j) = - one
       w(2, j) = f2(xx(j), yy(j))
       d(1, 2, j) = one
       d(2, 2, j) = zero
       w(3, j) = f3(xx(j), yy(j))
       d(1, 3, j) = zero
       d(2, 3, j) = one
       p(j) = one/six
    ENDDO
  END SUBROUTINE element_2d_p1

  SUBROUTINE element_2d_p1_boundary (face, d, w, n_w, l_Gs)
    !===Triangular element with linear interpolation
    !===and three Gauss integration points
    !===w(n_w, l_Gs)    : values of shape functions at Gauss points
    !===d(2, n_w, l_Gs) : derivatives values of shape functions at Gauss points
    !===p(l_Gs)         : weight for Gaussian quadrature at Gauss points
    ! 3
    ! 1 2  with orientation 1->2, 1->3, 2->3 (lowest to highest index convention)
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)  :: n_w, l_Gs
    INTEGER,                                INTENT(IN)  :: face
    REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(   n_w, l_Gs),  INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(l_Gs) :: xx, yy
    INTEGER :: j
    REAL(KIND=8) :: zero = 0,  one = 1, three=3, half = 0.5d0
    REAL(KIND=8) :: f1, f2, f3, x, y
    f1(x, y) = 1-x-y
    f2(x, y) = x
    f3(x, y) = y
    IF (face==1) THEN
       xx(1) = half + half/SQRT(three)
       xx(2) = half - half/SQRT(three)
       yy(1) = half - half/SQRT(three)
       yy(2) = half + half/SQRT(three)
    ELSE IF (face==2) THEN
       xx(1) = zero
       xx(2) = zero
       yy(1) = half - half/SQRT(three)
       yy(2) = half + half/SQRT(three)
    ELSE
       xx(1) = half - half/SQRT(three)
       xx(2) = half + half/SQRT(three)
       yy(1) = zero
       yy(2) = zero
    END IF
    DO j = 1, l_Gs
       w(1,j) = f1(xx(j), yy(j))
       d(1, 1, j) = -one
       d(2, 1, j) = -one
       w(2,j) = f2(xx(j), yy(j))
       d(1, 2, j) = one
       d(2, 2, j) = zero
       w(3,j) = f3(xx(j), yy(j))
       d(1, 3, j) = zero
       d(2, 3, j) = one
    ENDDO
  END SUBROUTINE element_2d_p1_boundary

  SUBROUTINE element_1d_p1 (w, d, p, n_ws, l_Gs)
    !===One-dimensional element with linear interpolation
    !===and two Gauss integration points
    !===w(n_w, l_G)    : values of shape functions at Gauss points
    !===d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
    !===p(l_G)         : weight for Gaussian quadrature at Gauss points
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)  :: n_ws, l_Gs
    REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p
    REAL(KIND=8), DIMENSION(l_Gs) :: xx
    INTEGER :: j
    REAL(KIND=8) ::  one = 1.d0,  two = 2.d0,  three = 3.d0
    REAL(KIND=8) :: f1, f2, x
    f1(x) = (one - x)/two
    f2(x) = (x + one)/two
    xx(1) = - one/SQRT(three)
    xx(2) = + one/SQRT(three)
    DO j = 1, l_Gs
       w(1, j) = f1(xx(j))
       d(1, 1, j) = - one/two
       w(2, j) = f2(xx(j))
       d(1, 2, j) = + one/two
       p(j) = one
    ENDDO
  END SUBROUTINE element_1d_p1

  SUBROUTINE element_1d_p1_at_nodes (d, n_ws)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)  :: n_ws
    REAL(KIND=8), DIMENSION(n_ws, n_ws), INTENT(OUT) :: d
    INTEGER :: j
    REAL(KIND=8) ::  one = 1.d0,  two = 2.d0
    DO j = 1, n_ws
       d(1, j) = - one/two
       d(2, j) = + one/two
    ENDDO
  END SUBROUTINE element_1d_p1_at_nodes

END MODULE GP_2d_p1
