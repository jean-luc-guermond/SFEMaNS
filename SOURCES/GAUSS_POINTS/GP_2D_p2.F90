!===
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!Author:  Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE GP_2d_p2
  PRIVATE
  PUBLIC element_2d_p2, element_2d_p2_boundary, element_1d_p2
CONTAINS
  SUBROUTINE element_2d_p2 (w, d, p, n_w, l_G)
    !===Triangular element with quadratic  interpolation
    !===and seven Gauss integration points
    !===w(n_w, l_G)    : values of shape functions at Gauss points
    !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
    !===p(l_G)         : weight for Gaussian quadrature at Gauss points
    ! 3
    ! 5 4     with orientation 1->2, 1->3, 2->3 (lowest to highest index convention)
    ! 1 6 2
    IMPLICIT NONE
    INTEGER,                              INTENT(IN)  :: n_w, l_G
    REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p
    REAL(KIND=8), DIMENSION(l_G) :: xx, yy
    INTEGER :: j
    REAL(KIND=8) :: zero = 0,  half  = 0.5,  one  = 1,  &
         two  = 2,  three = 3,    four = 4,  &
         five = 5,   nine = 9
    REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
         df1x, df2x, df3x, df4x, df5x, df6x, &
         df1y, df2y, df3y, df4y, df5y, df6y, &
         x, y, r, a,  s1, s2, t1, t2, b1, b2,  area, sq

    f1(x, y) = (half - x - y) * (one - x - y) * two
    f2(x, y) = x * (x - half) * two
    f3(x, y) = y * (y - half) * two
    f4(x, y) = x * y * four
    f5(x, y) = y * (one - x - y) * four
    f6(x, y) = x * (one - x - y) * four

    df1x(x, y) = -three + four * (x + y)
    df2x(x, y) = (two*x - half) * two
    df3x(x, y) = zero
    df4x(x, y) =  y * four
    df5x(x, y) = -y * four
    df6x(x, y) = (one - two*x - y) * four

    df1y(x, y) = -three + four * (x + y)
    df2y(x, y) = zero
    df3y(x, y) = (two*y - half) * two
    df4y(x, y) =  x * four
    df5y(x, y) = (one - x - two*y) * four
    df6y(x, y) = -x * four
    !===Degree 5; 7 Points;  Stroud: p. 314, Approximate calculation of
    !===Multiple integrals (Prentice--Hall), 1971.
    area = one/two
    sq = SQRT(three*five)
    r  = one/three;                          a = area * nine/40
    s1 = (6 - sq)/21;  t1 = (9 + 2*sq)/21;  b1 = area * (155 - sq)/1200
    s2 = (6 + sq)/21;  t2 = (9 - 2*sq)/21;  b2 = area * (155 + sq)/1200
    xx(1) = r;    yy(1) = r;    p(1) = a
    xx(2) = s1;   yy(2) = s1;   p(2) = b1
    xx(3) = s1;   yy(3) = t1;   p(3) = b1
    xx(4) = t1;   yy(4) = s1;   p(4) = b1
    xx(5) = s2;   yy(5) = s2;   p(5) = b2
    xx(6) = s2;   yy(6) = t2;   p(6) = b2
    xx(7) = t2;   yy(7) = s2;   p(7) = b2

    DO j = 1, l_G
       w(1, j) =  f1 (xx(j), yy(j))
       d(1, 1, j) = df1x(xx(j), yy(j))
       d(2, 1, j) = df1y(xx(j), yy(j))
       w(2, j) =  f2 (xx(j), yy(j))
       d(1, 2, j) = df2x(xx(j), yy(j))
       d(2, 2, j) = df2y(xx(j), yy(j))
       w(3, j) =  f3 (xx(j), yy(j))
       d(1, 3, j) = df3x(xx(j), yy(j))
       d(2, 3, j) = df3y(xx(j), yy(j))
       w(4, j) =  f4 (xx(j), yy(j))
       d(1, 4, j) = df4x(xx(j), yy(j))
       d(2, 4, j) = df4y(xx(j), yy(j))
       w(5, j) =  f5 (xx(j), yy(j))
       d(1, 5, j) = df5x(xx(j), yy(j))
       d(2, 5, j) = df5y(xx(j), yy(j))
       w(6, j) =  f6 (xx(j), yy(j))
       d(1, 6, j) = df6x(xx(j), yy(j))
       d(2, 6, j) = df6y(xx(j), yy(j))
    ENDDO
  END SUBROUTINE element_2d_p2

  SUBROUTINE element_2d_p2_boundary (face, d, w, n_w, l_Gs)
    !===Triangular element with quadratic interpolation
    !===and seven Gauss integration points
    !===w(n_w, l_Gs)    : values of shape functions at Gauss points
    !===d(2, n_w, l_Gs) : derivatives values of shape functions at Gauss points
    !===p(l_Gs)         : weight for Gaussian quadrature at Gauss points
    ! 3
    ! 5 4     with orientation 1->2, 1->3, 2->3 (lowest to highest index convention)
    ! 1 6 2
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)  :: n_w, l_Gs
    INTEGER,                                INTENT(IN)  :: face
    REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(   n_w, l_Gs),  INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(l_Gs) :: xx, yy
    INTEGER :: j
    REAL(KIND=8) :: zero = 0,  half  = 0.5d0,  one  = 1,  &
         two  = 2,  three = 3,    four = 4,  &
         five = 5
    REAL(KIND=8) ::  f1,   f2,   f3,   f4,   f5,   f6,  &
         df1x, df2x, df3x, df4x, df5x, df6x, &
         df1y, df2y, df3y, df4y, df5y, df6y, &
         x, y
    f1(x, y) = (half - x - y) * (one - x - y) * two
    f2(x, y) = x * (x - half) * two
    f3(x, y) = y * (y - half) * two
    f4(x, y) = x * y * four
    f5(x, y) = y * (one - x - y) * four
    f6(x, y) = x * (one - x - y) * four
    df1x(x, y) = -three + four * (x + y)
    df2x(x, y) = (two*x - half) * two
    df3x(x, y) = zero
    df4x(x, y) =  y * four
    df5x(x, y) = -y * four
    df6x(x, y) = (one - two*x - y) * four
    df1y(x, y) = -three + four * (x + y)
    df2y(x, y) = zero
    df3y(x, y) = (two*y - half) * two
    df4y(x, y) =  x * four
    df5y(x, y) = (one - x - two*y) * four
    df6y(x, y) = -x * four
    IF (face==1) THEN
       xx(1) = half + half*SQRT(three/five)
       xx(2) = half
       xx(3) = half - half*SQRT(three/five)
       yy(1) = half - half*SQRT(three/five)
       yy(2) = half
       yy(3) = half + half*SQRT(three/five)
    ELSE IF (face==2) THEN
       xx(1) = zero
       xx(2) = zero
       xx(3) = zero
       yy(1) = half - half*SQRT(three/five)
       yy(2) = half
       yy(3) = half + half*SQRT(three/five)
    ELSE
       xx(1) = half - half*SQRT(three/five)
       xx(2) = half
       xx(3) = half + half*SQRT(three/five)
       yy(1) = zero
       yy(2) = zero
       yy(3) = zero
    END IF

    DO j = 1, l_Gs
       w(1, j) =  f1 (xx(j), yy(j))
       d(1, 1, j) = df1x(xx(j), yy(j))
       d(2, 1, j) = df1y(xx(j), yy(j))
       w(2, j) =  f2 (xx(j), yy(j))
       d(1, 2, j) = df2x(xx(j), yy(j))
       d(2, 2, j) = df2y(xx(j), yy(j))
       w(3, j) =  f3 (xx(j), yy(j))
       d(1, 3, j) = df3x(xx(j), yy(j))
       d(2, 3, j) = df3y(xx(j), yy(j))
       w(4, j) =  f4 (xx(j), yy(j))
       d(1, 4, j) = df4x(xx(j), yy(j))
       d(2, 4, j) = df4y(xx(j), yy(j))
       w(5, j) =  f5 (xx(j), yy(j))
       d(1, 5, j) = df5x(xx(j), yy(j))
       d(2, 5, j) = df5y(xx(j), yy(j))
       w(6, j) =  f6 (xx(j), yy(j))
       d(1, 6, j) = df6x(xx(j), yy(j))
       d(2, 6, j) = df6y(xx(j), yy(j))
    ENDDO
  END SUBROUTINE element_2d_p2_boundary

  SUBROUTINE element_1d_p2(w, d, p, n_ws, l_Gs)
    !===One-dimensional element with quadratic interpolation
    !===and three Gauss integration points
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
    REAL(KIND=8) :: zero = 0,  one = 1,  two = 2,  three = 3,  &
         five = 5, eight= 8,   nine = 9
    REAL(KIND=8) :: f1, f2, f3, df1, df2, df3, x
    f1(x) = (x - one)*x/two
    f2(x) = (x + one)*x/two
    f3(x) = (x + one)*(one - x)
    df1(x) = (two*x - one)/two
    df2(x) = (two*x + one)/two
    df3(x) = -two*x
    xx(1) = -SQRT(three/five)
    xx(2) =  zero
    xx(3) =  SQRT(three/five)
    p(1)  =  five/nine
    p(2)  =  eight/nine
    p(3)  =  five/nine
    DO j = 1, l_Gs
       w(1, j) =  f1(xx(j))
       d(1, 1, j) = df1(xx(j))
       w(2, j) =  f2(xx(j))
       d(1, 2, j) = df2(xx(j))
       w(3, j) =  f3(xx(j))
       d(1, 3, j) = df3(xx(j))
    ENDDO
  END SUBROUTINE element_1d_p2
END MODULE GP_2d_p2
