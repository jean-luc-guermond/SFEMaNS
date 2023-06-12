MODULE basis_change
  PUBLIC :: p1_p2, p2_p3
  PRIVATE
  !===Compute change of basis matrices for p1 to p2 and p2 to p3
  !===P(k-1) basis: phi_i; P_k basis: psi_j
  !===phi_i = sum_j (a_ij psi_j)
  !===a_ij = phi_i(x_j), where x_j is the jth Lagrange node such that psi_k(x_j)=delta_kj
CONTAINS
  subroutine p1_p2(aij)
    IMPLICIT NONE
    INTEGER, PARAMETER :: type_fe_p2 = 1, nw_p2=(type_fe_p2+2)*(type_fe_p2+1)/2, &
         type_fe_p3 = 2, nw_p3=(type_fe_p3+2)*(type_fe_p3+1)/2
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
    INTEGER, DIMENSION(nw_p3) :: Cart_FE
    REAL(KIND=8) ::  f1, f2, f3, x, y, delta, one = 1.d0
    REAL(KIND=8), DIMENSION(nw_p3) :: xx, yy, xxp, yyp
    INTEGER :: j,  h, k

    f1(x, y) = one - x - y
    f2(x, y) = x
    f3(x, y) = y

    ALLOCATE(aij(nw_p2,nw_p3))
    delta = 1.d0/type_fe_p3
    j = 0
    DO h = 1, type_fe_p3+1
       DO k = 1, type_fe_p3+2-h
          j = j+1
          xxp(j) = (h-1)*delta
          yyp(j) = (k-1)*delta
       END DO
    END DO

    Cart_FE(1) = 1
    Cart_FE(2) = 5
    Cart_FE(3) = 3
    Cart_FE(4) = 6
    Cart_FE(5) = 4
    Cart_FE(6) = 2
    xx(Cart_FE) = xxp
    yy(Cart_FE) = yyp

    DO j = 1, nw_p3
       aij(1,j) = f1(xx(j),yy(j))
       aij(2,j) = f2(xx(j),yy(j))
       aij(3,j) = f3(xx(j),yy(j))
    END DO
  END subroutine p1_p2

  subroutine p2_p3(aij)
    IMPLICIT NONE
    INTEGER, PARAMETER :: type_fe_p2 = 2, nw_p2=(type_fe_p2+2)*(type_fe_p2+1)/2, &
         type_fe_p3 = 3, nw_p3=(type_fe_p3+2)*(type_fe_p3+1)/2
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
    INTEGER, DIMENSION(nw_p3) :: Cart_FE
    REAL(KIND=8) ::  f1, f2, f3, f4, f5, f6, x, y, delta, half  = 0.5d0, one = 1.d0, two=2.d0, four=4.d0
    REAL(KIND=8), DIMENSION(nw_p3) :: xx, yy, xxp, yyp
    INTEGER :: j, h, k

    f1(x, y) = (half - x - y) * (one - x - y) * two
    f2(x, y) = x * (x - half) * two
    f3(x, y) = y * (y - half) * two
    f4(x, y) = x * y * four
    f5(x, y) = y * (one - x - y) * four
    f6(x, y) = x * (one - x - y) * four

    ALLOCATE(aij(nw_p2,nw_p3))
    delta = 1.d0/type_fe_p3
    j = 0
    DO h = 1, type_fe_p3+1
       DO k = 1, type_fe_p3+2-h
          j = j+1
          xxp(j) = (h-1)*delta
          yyp(j) = (k-1)*delta
       END DO
    END DO

    Cart_FE(1) = 1
    Cart_FE(2) = 6
    Cart_FE(3) = 7
    Cart_FE(4) = 3
    Cart_FE(5) = 8
    Cart_FE(6) = 10
    Cart_FE(7) = 5
    Cart_FE(8) = 9
    Cart_FE(9) = 4
    Cart_FE(10) = 2
    xx(Cart_FE) = xxp
    yy(Cart_FE) = yyp

    DO j = 1, nw_p3
       aij(1,j) = f1(xx(j),yy(j))
       aij(2,j) = f2(xx(j),yy(j))
       aij(3,j) = f3(xx(j),yy(j))
       aij(4,j) = f4(xx(j),yy(j))
       aij(5,j) = f5(xx(j),yy(j))
       aij(6,j) = f6(xx(j),yy(j))
    END DO
  END subroutine p2_p3

END MODULE basis_change

!!$program t
!!$  USE basis_change
!!$  REAL(KIND=8), DIMENSION(:,:), POINTER :: aij_p2, aij_p3
!!$  INTEGER :: i
!!$  CALL p1_p2(aij_p2)
!!$  DO i = 1, 3
!!$     write(*,'(7(f8.3,3x))') aij_p2(i,:)
!!$  END DO
!!$  DO i = 1, 6
!!$     write(*,'(f8.3)') sum(aij_p2(:,i))
!!$  END DO
!!$  CALL p2_p3(aij_p3)
!!$  write(*,*) aij_p2(:,4)
!!$  !DO i = 1, 6
!!$  !   write(*,'(11(f8.3,3x))') aij_p3(i,:)
!!$  !END DO
!!$  !DO i = 1, 10
!!$  !   write(*,'(f8.3)') sum(aij_p3(:,i))
!!$  !END DO
!!$end program t
