!===
!Author: Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE GP_2d_p3
  PRIVATE
  PUBLIC element_2d_p3, element_2d_p3_boundary, element_1d_p3
CONTAINS
  SUBROUTINE element_2d_p3 (w, d, p, n_w, l_G)
    !===triangular element with cubic interpolation
    !===and seven Gauss integration points
    !===w(n_w, l_G) : values of shape functions at Gauss points
    !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
    !===p(l_G) : weight for Gaussian quadrature at Gauss points
    !=== 3
    !=== 7  5
    !=== 6  10 4
    !=== 1  8  9  2
    IMPLICIT NONE
    INTEGER,                              INTENT(IN)  :: n_w, l_G
    REAL(KIND=8), DIMENSION(   n_w, l_G), INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(2, n_w, l_G), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(l_G),         INTENT(OUT) :: p
    REAL(KIND=8), DIMENSION(l_G)                      :: xx, yy
    INTEGER :: j
    REAL(KIND=8) :: one=1.d0, two=2.d0, three=3.d0, five=5.d0, nine=9.d0
    REAL(KIND=8) ::  f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,   &
         df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x, df9x, df10x, &
         df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y, df9y, df10y, &
         x, y, r, a, s1, s2, t1, t2, b1, b2, area, sq

    f1(x,y) = -0.9d1/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y - 0.27d2/0.2d1*x*y**2 &
         - 0.9d1/0.2d1*y**3 + 0.9d1*x**2 + 0.18d2*x*y + 0.9d1*y**2 &
         - 0.11d2/0.2d1*x - 0.11d2/0.2d1*y + 0.1d1
    f2(x,y) = 0.9d1/0.2d1*x**3 - 0.9d1/0.2d1*x**2 + x
    f3(x,y) = 0.9d1/0.2d1*y**3 - 0.9d1/0.2d1*y**2 + y
    f4(x,y) = 0.27d2/0.2d1*x**2*y - 0.9d1/0.2d1*x*y
    f5(x,y) = 0.27d2/0.2d1*x*y**2 - 0.9d1/0.2d1*x*y
    f6(x,y) = 0.27d2/0.2d1*x**2*y + 0.27d2*x*y**2 + 0.27d2 / 0.2d1*y**3 &
         - 0.45d2/0.2d1*x*y - 0.45d2/0.2d1*y**2 + 0.9d1*y
    f7(x,y) = -0.27d2/0.2d1*x*y**2 - 0.27d2/0.2d1*y**3 + 0.9d1/0.2d1*x*y &
         + 0.18d2*y**2 - 0.9d1/0.2d1*y
    f8(x,y) = 0.27d2/0.2d1*x**3 + 0.27d2*x**2*y + 0.27d2/0.2d1*x*y**2 &
         - 0.45d2/0.2d1*x**2 - 0.45d2/0.2d1*x*y + 0.9d1*x
    f9(x,y) = -0.27d2/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y + 0.18d2*x**2 &
         + 0.9d1/0.2d1*x*y - 0.9d1/0.2d1*x
    f10(x,y) = -27*x**2*y - 27*x*y**2 + 27*x*y

    df1x(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
         + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
    df2x(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1*x + 0.1d1
    df3x(x,y) = 0
    df4x(x,y) = 27*x*y - 0.9d1/0.2d1 *y
    df5x(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1/0.2d1*y
    df6x(x,y) = 27*x*y + 27*y**2 - 0.45d2/0.2d1*y
    df7x(x,y) = -0.27d2/0.2d1*y**2 + 0.9d1/0.2d1*y
    df8x(x,y) = 0.81d2/0.2d1*x**2 + 0.54d2*x*y - 0.45d2*x &
         + 0.27d2/0.2d1*y**2 - 0.45d2/0.2d1*y + 0.9d1
    df9x(x,y) = -0.81d2/0.2d1*x**2 + 0.36d2*x - 0.27d2*x*y &
         + 0.9d1/0.2d1*y - 0.9d1/0.2d1
    df10x(x,y) = -54*x*y - 27*y**2 + 27*y

    df1y(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1* y**2 &
         + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
    df2y(x,y) = 0.d0
    df3y(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1*y + 0.1d1
    df4y(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1/0.2d1*x
    df5y(x,y) = 27*x*y - 0.9d1/0.2d1*x
    df6y(x,y) = 54*x*y + 0.81d2/0.2d1*y**2 &
         - 45*y + 0.27d2/0.2d1*x**2 - 0.45d2/0.2d1*x + 0.9d1
    df7y(x,y) = -0.81d2/0.2d1*y**2 + 0.36d2*y - 0.27d2*x*y + 0.9d1/0.2d1*x - 0.9d1/0.2d1
    df8y(x,y) = 27*x**2 + 27*x*y - 0.45d2/0.2d1*x
    df9y(x,y) = -0.27d2/0.2d1*x**2 + 0.9d1/0.2d1*x
    df10y(x,y) = -27*x**2 - 54*x*y + 27*x

    !===Degree 5; 7 Points;  Stroud: p. 314, Approximate calculation of
    !===Multiple integrals (Prentice--Hall), 1971.
    area = one/two
    sq = SQRT(three*five)
    r  = one/three;                          a = area * nine/40
    s1 = (6 - sq)/21;  t1 = (9 + 2*sq)/21;  b1 = area * (155 - sq)/1200
    s2 = (6 + sq)/21;  t2 = (9 - 2*sq)/21;  b2 = area * (155 + sq)/1200
    IF (l_G==7) THEN !===5th order accurate
       xx(1) = r;    yy(1) = r;    p(1) = a
       xx(2) = s1;   yy(2) = s1;   p(2) = b1
       xx(3) = s1;   yy(3) = t1;   p(3) = b1
       xx(4) = t1;   yy(4) = s1;   p(4) = b1
       xx(5) = s2;   yy(5) = s2;   p(5) = b2
       xx(6) = s2;   yy(6) = t2;   p(6) = b2
       xx(7) = t2;   yy(7) = s2;   p(7) = b2
    ELSE IF (l_G==12) THEN !===6th orer accurate
       xx(1) = 0.501426509658179d0; yy(1) = 0.249286745170910d0; p(1) = 0.116786275726379d0
       xx(2) = 0.249286745170910d0; yy(2) = 0.249286745170910d0; p(2) = 0.116786275726379d0
       xx(3) = 0.249286745170910d0; yy(3) = 0.501426509658179d0; p(3) = 0.116786275726379d0
       xx(4) = 0.873821971016996d0; yy(4) = 0.063089014491502d0; p(4) = 0.050844906370207d0
       xx(5) = 0.063089014491502d0; yy(5) = 0.063089014491502d0; p(5) = 0.050844906370207d0
       xx(6) = 0.063089014491502d0; yy(6) = 0.873821971016996d0; p(6) = 0.050844906370207d0
       xx(7) = 0.053145049844817d0; yy(7) = 0.310352451033784d0; p(7) = 0.082851075618374d0
       xx(8) = 0.310352451033784d0; yy(8) = 0.636502499121399d0; p(8) = 0.082851075618374d0
       xx(9) = 0.636502499121399d0; yy(9) = 0.053145049844817d0; p(9) = 0.082851075618374d0
       xx(10)= 0.053145049844817d0; yy(10)= 0.636502499121399d0; p(10)= 0.082851075618374d0
       xx(11)= 0.636502499121399d0; yy(11)= 0.310352451033784d0; p(11)= 0.082851075618374d0
       xx(12)= 0.310352451033784d0; yy(12)= 0.053145049844817d0; p(12) = 0.082851075618374d0
       p = p*area
    END IF

    DO j = 1, l_G
       w(1, j)    =   f1(xx(j), yy(j))
       d(1, 1, j) = df1x(xx(j), yy(j))
       d(2, 1, j) = df1y(xx(j), yy(j))
       w(2, j)    =   f2(xx(j), yy(j))
       d(1, 2, j) = df2x(xx(j), yy(j))
       d(2, 2, j) = df2y(xx(j), yy(j))
       w(3, j)    =   f3(xx(j), yy(j))
       d(1, 3, j) = df3x(xx(j), yy(j))
       d(2, 3, j) = df3y(xx(j), yy(j))
       w(4, j)    =   f4(xx(j), yy(j))
       d(1, 4, j) = df4x(xx(j), yy(j))
       d(2, 4, j) = df4y(xx(j), yy(j))
       w(5, j)    =   f5(xx(j), yy(j))
       d(1, 5, j) = df5x(xx(j), yy(j))
       d(2, 5, j) = df5y(xx(j), yy(j))
       w(6, j)    =   f6(xx(j), yy(j))
       d(1, 6, j) = df6x(xx(j), yy(j))
       d(2, 6, j) = df6y(xx(j), yy(j))
       w(7, j)    =   f7(xx(j), yy(j))
       d(1, 7, j) = df7x(xx(j), yy(j))
       d(2, 7, j) = df7y(xx(j), yy(j))
       w(8, j)    =   f8(xx(j), yy(j))
       d(1, 8, j) = df8x(xx(j), yy(j))
       d(2, 8, j) = df8y(xx(j), yy(j))
       w(9, j)    =   f9(xx(j), yy(j))
       d(1, 9, j) = df9x(xx(j), yy(j))
       d(2, 9, j) = df9y(xx(j), yy(j))
       w(10, j)   =  f10(xx(j), yy(j))
       d(1, 10, j)=df10x(xx(j), yy(j))
       d(2, 10, j)=df10y(xx(j), yy(j))
    ENDDO
  END SUBROUTINE element_2d_p3

  SUBROUTINE element_2d_p3_boundary (face, d, w, n_w, l_Gs)
    !===triangular element with cubic interpolation
    !===and seven Gauss integration points
    !===d(2, n_w, l_G) : derivatives values of shape functions at Gauss points
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)  :: n_w, l_Gs
    INTEGER,                                INTENT(IN)  :: face
    REAL(KIND=8), DIMENSION(2, n_w,  l_Gs), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(   n_w, l_Gs),  INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(l_Gs)                       :: xx, yy, gp
    INTEGER :: j, l
    REAL(KIND=8) :: half=0.5d0
    REAL(KIND=8) ::  f1, f2, f3, f4, f5, f6, f7, f8, f9, f10,   &
         df1x, df2x, df3x, df4x, df5x, df6x, df7x, df8x, df9x, df10x, &
         df1y, df2y, df3y, df4y, df5y, df6y, df7y, df8y, df9y, df10y, &
         x, y

    f1(x,y) = -0.9d1/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y - 0.27d2/0.2d1*x*y**2 &
         - 0.9d1/0.2d1*y**3 + 0.9d1*x**2 + 0.18d2*x*y + 0.9d1*y**2 &
         - 0.11d2/0.2d1*x - 0.11d2/0.2d1*y + 0.1d1
    f2(x,y) = 0.9d1/0.2d1*x**3 - 0.9d1/0.2d1*x**2 + x
    f3(x,y) = 0.9d1/0.2d1*y**3 - 0.9d1/0.2d1*y**2 + y
    f4(x,y) = 0.27d2/0.2d1*x**2*y - 0.9d1/0.2d1*x*y
    f5(x,y) = 0.27d2/0.2d1*x*y**2 - 0.9d1/0.2d1*x*y
    f6(x,y) = 0.27d2/0.2d1*x**2*y + 0.27d2*x*y**2 + 0.27d2/0.2d1*y**3 &
         - 0.45d2/0.2d1*x*y - 0.45d2/0.2d1*y**2 + 0.9d1*y
    f7(x,y) = -0.27d2/0.2d1*x*y**2 - 0.27d2/0.2d1*y**3 + 0.9d1/0.2d1*x*y &
         + 0.18d2*y**2 - 0.9d1/0.2d1*y
    f8(x,y) = 0.27d2/0.2d1*x**3 + 0.27d2*x**2*y + 0.27d2/0.2d1*x*y**2 &
         - 0.45d2/0.2d1*x**2 - 0.45d2/0.2d1*x*y + 0.9d1*x
    f9(x,y) = -0.27d2/0.2d1*x**3 - 0.27d2/0.2d1*x**2*y + 0.18d2*x**2 &
         + 0.9d1/0.2d1*x*y - 0.9d1/0.2d1*x
    f10(x,y) = -27*x**2*y - 27*x*y**2 + 27*x*y

    df1x(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
         + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
    df2x(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1*x + 0.1d1
    df3x(x,y) = 0
    df4x(x,y) = 27*x*y - 0.9d1/0.2d1 *y
    df5x(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1/0.2d1*y
    df6x(x,y) = 27*x*y + 27*y**2 - 0.45d2/0.2d1*y
    df7x(x,y) = -0.27d2/0.2d1*y**2 + 0.9d1/0.2d1*y
    df8x(x,y) = 0.81d2/0.2d1*x**2 + 0.54d2*x*y - 0.45d2*x &
         + 0.27d2/0.2d1*y**2 - 0.45d2/0.2d1*y + 0.9d1
    df9x(x,y) = -0.81d2/0.2d1*x**2 + 0.36d2*x - 0.27d2*x*y &
         + 0.9d1/0.2d1*y - 0.9d1/0.2d1
    df10x(x,y) = -54*x*y - 27*y**2 + 27*y

    df1y(x,y) = -0.27d2/0.2d1*x**2 - 0.27d2*x*y - 0.27d2/0.2d1*y**2 &
         + 0.18d2*x + 0.18d2*y - 0.11d2/0.2d1
    df2y(x,y) = 0.d0
    df3y(x,y) = 0.27d2/0.2d1*y**2 - 0.9d1*y + 0.1d1
    df4y(x,y) = 0.27d2/0.2d1*x**2 - 0.9d1/0.2d1*x
    df5y(x,y) = 27*x*y - 0.9d1/0.2d1*x
    df6y(x,y) = 54*x*y + 0.81d2/0.2d1*y**2 &
         - 45*y + 0.27d2/0.2d1*x**2 - 0.45d2/0.2d1*x + 0.9d1
    df7y(x,y) = -0.81d2/0.2d1*y**2 + 0.36d2*y - 0.27d2*x*y + 0.9d1/0.2d1*x - 0.9d1/0.2d1
    df8y(x,y) = 27*x**2 + 27*x*y - 0.45d2/0.2d1*x
    df9y(x,y) = -0.27d2/0.2d1*x**2 + 0.9d1/0.2d1*x
    df10y(x,y) = -27*x**2 - 54*x*y + 27*x

    gp(1) = -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
    gp(2) = -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    gp(3) =  SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    gp(4) =  SQRT((15.d0+2*SQRT(30.d0))/35.d0)

    IF (face==1) THEN
       DO l = 1, l_Gs
          xx(l) = 1.d0-(gp(l)+1.d0)*half
          yy(l) = 0.d0+(gp(l)+1.d0)*half
       END DO
    ELSE IF (face==2) THEN
       DO l = 1, l_Gs
          xx(l) = 0.d0
          yy(l) = (gp(l)+1.d0)*half
       END DO
    ELSE
       DO l = 1, l_Gs
          xx(l) = (gp(l)+1.d0)*half
          yy(l) = 0.d0
       END DO
    END IF

    DO j = 1, l_Gs
       w(1, j)    =   f1(xx(j), yy(j))
       d(1, 1, j) = df1x(xx(j), yy(j))
       d(2, 1, j) = df1y(xx(j), yy(j))
       w(2, j)    =   f2(xx(j), yy(j))
       d(1, 2, j) = df2x(xx(j), yy(j))
       d(2, 2, j) = df2y(xx(j), yy(j))
       w(3, j)    =   f3(xx(j), yy(j))
       d(1, 3, j) = df3x(xx(j), yy(j))
       d(2, 3, j) = df3y(xx(j), yy(j))
       w(4, j)    =   f4(xx(j), yy(j))
       d(1, 4, j) = df4x(xx(j), yy(j))
       d(2, 4, j) = df4y(xx(j), yy(j))
       w(5, j)    =   f5(xx(j), yy(j))
       d(1, 5, j) = df5x(xx(j), yy(j))
       d(2, 5, j) = df5y(xx(j), yy(j))
       w(6, j)    =   f6(xx(j), yy(j))
       d(1, 6, j) = df6x(xx(j), yy(j))
       d(2, 6, j) = df6y(xx(j), yy(j))
       w(7, j)    =   f7(xx(j), yy(j))
       d(1, 7, j) = df7x(xx(j), yy(j))
       d(2, 7, j) = df7y(xx(j), yy(j))
       w(8, j)    =   f8(xx(j), yy(j))
       d(1, 8, j) = df8x(xx(j), yy(j))
       d(2, 8, j) = df8y(xx(j), yy(j))
       w(9, j)    =   f9(xx(j), yy(j))
       d(1, 9, j) = df9x(xx(j), yy(j))
       d(2, 9, j) = df9y(xx(j), yy(j))
       w(10, j)   =  f10(xx(j), yy(j))
       d(1, 10, j)=df10x(xx(j), yy(j))
       d(2, 10, j)=df10y(xx(j), yy(j))
    ENDDO

  END SUBROUTINE element_2d_p3_boundary

  SUBROUTINE element_1d_p3(w, d, p, n_ws, l_Gs)
    !===one-dimensional element with cubic interpolation
    !===and 4 Gauss integration points
    !===w(n_w, l_G)    : values of shape functions at Gauss points
    !===d(1, n_w, l_G) : derivatives values of shape functions at Gauss points
    !===p(l_G)         : weight for Gaussian quadrature at Gauss points
    !===Enumeration: 1  3  4  2
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)  :: n_ws, l_Gs
    REAL(KIND=8), DIMENSION(   n_ws, l_Gs), INTENT(OUT) :: w
    REAL(KIND=8), DIMENSION(1, n_ws, l_Gs), INTENT(OUT) :: d
    REAL(KIND=8), DIMENSION(l_Gs),          INTENT(OUT) :: p
    REAL(KIND=8), DIMENSION(l_Gs) :: xx
    INTEGER :: j
    REAL(KIND=8) :: f1, f2, f3, f4, df1, df2, df3, df4, x

    f1(x) = -0.1D1/0.16D2 + x/0.16D2 + 0.9D1/0.16D2*x**2 - 0.9D1/0.16D2*x**3
    f2(x) = -x/0.16D2 + 0.9D1/0.16D2*x**2 + 0.9D1/0.16D2*x**3 - 0.1D1/0.16D2
    f3(x) = -0.9D1/0.16D2*x**2 + 0.27D2/0.16D2*x**3 + 0.9D1/0.16D2 - 0.27D2/0.16D2*x
    f4(x) = -0.9D1/0.16D2*x**2 - 0.27D2/0.16D2*x**3 + 0.9D1/0.16D2 + 0.27D2/0.16D2*x

    df1(x) = 0.1D1/0.16D2 - 0.27D2/0.16D2*x**2 + 0.9D1/0.8D1*x
    df2(x)  = -0.1D1/0.16D2 + 0.27D2/0.16D2*x**2 + 0.9D1/0.8D1*x
    df3(x) = -0.27D2/0.16D2 - 0.9D1/0.8D1*x + 0.81D2/0.16D2*x**2
    df4(x) = -0.9D1/0.8D1*x - 0.81D2/0.16D2*x**2 + 0.27D2/0.16D2

    xx(1) = -SQRT((15.d0+2*SQRT(30.d0))/35.d0)
    xx(2) = -SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    xx(3) =  SQRT((15.d0-2*SQRT(30.d0))/35.d0)
    xx(4) =  SQRT((15.d0+2*SQRT(30.d0))/35.d0)
    p(1) = (18.d0-SQRT(30.d0))/36.d0
    p(2) = (18.d0+SQRT(30.d0))/36.d0
    p(3) = (18.d0+SQRT(30.d0))/36.d0
    p(4) = (18.d0-SQRT(30.d0))/36.d0

    DO j = 1, l_Gs
       w(1, j) =  f1(xx(j))
       d(1, 1, j) = df1(xx(j))
       w(2, j) =  f2(xx(j))
       d(1, 2, j) = df2(xx(j))
       w(3, j) =  f3(xx(j))
       d(1, 3, j) = df3(xx(j))
       w(4, j) =  f4(xx(j))
       d(1, 4, j) = df4(xx(j))
    ENDDO
  END SUBROUTINE element_1d_p3
END MODULE GP_2d_p3
