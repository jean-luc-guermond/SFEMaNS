!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE Dir_nodes

CONTAINS

  SUBROUTINE precond_mat(ia, aa, vect)
    IMPLICIT NONE
    INTEGER,      DIMENSION(:), INTENT(IN)    :: ia
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: aa
    REAL(KIND=8), DIMENSION(:), INTENT(OUT)   :: vect
    INTEGER                                   :: ii,j
    REAL(KIND=8)                              :: S


    !    write(*,*) 'SIZE(ia)=',SIZE(ia)
    !    write(*,*) 'SIZE(aa)=',SIZE(aa)
    !    write(*,*) 'SIZE(vect)=',SIZE(vect)

    DO ii=1, SIZE(ia)-1     !pour toutes les lignes
       S = 0.d0
       DO j=ia(ii), ia(ii+1)-1
          S = S + ABS(aa(j))
          !S = MAX(S,ABS(aa(j)))
          !MODIF 16 MARS 2007
          !S = 1.d0
          !MODIF 16 MARS 2007
       ENDDO
       !WRITE(*,*) 'S= ', S
       vect(ii) = 1.d0/S
       DO j=ia(ii), ia(ii+1)-1
          aa(j) = aa(j)/S
       ENDDO
    ENDDO


  END SUBROUTINE precond_mat



  SUBROUTINE Dirichlet_M (js, ia, ja,  a0)
    !=================================================

    !  Enforces Dirichlet boundary conditions

    IMPLICIT NONE
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: js
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    INTEGER :: n, i, p

    DO n = 1, SIZE(js); i = js(n)

       DO p = ia(i), ia(i+1) - 1
          IF (ja(p) == i) THEN
             a0(p) = 1.d0

          ELSE
             a0(p)=0.d0
          END IF
       END DO

    END DO

  END SUBROUTINE Dirichlet_M

  SUBROUTINE Dirichlet_rescale_M (js, ia, ja,  alpha, a0)
    !=================================================

    !  Enforces Dirichlet boundary conditions

    IMPLICIT NONE
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: js
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8),   INTENT(IN) :: alpha

    INTEGER :: n, i, p

    DO n = 1, SIZE(js); i = js(n)
       DO p = ia(i), ia(i+1) - 1
          IF (ja(p) == i) THEN
             a0(p) = alpha
          ELSE
             a0(p)=0.d0
          END IF
       END DO
    END DO

  END SUBROUTINE Dirichlet_rescale_M

  SUBROUTINE Dirichlet_bloc_M (js_D, ia, ja,  a0)
    !=================================================

    !  Enforces Dirichlet boundary conditions

    USE dyn_line;
    IMPLICIT NONE
    TYPE(dyn_int_line),  DIMENSION(:)           :: js_D   ! Dirichlet nodes
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: n, i, k, p, np

    np = (SIZE(ia) - 1)/SIZE(js_D)
    DO k = 1, SIZE(js_D)
       DO n = 1, SIZE(js_D(k)%DIL); i = js_D(k)%DIL(n) + (k-1)*np
          DO p = ia(i), ia(i+1) - 1
             IF (ja(p) == i) THEN
                a0(p) = 1.d0
             ELSE
                a0(p)=0.d0
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE Dirichlet_bloc_M

  SUBROUTINE Dir_nodes_size(js, sides, Dir, j_D_size)

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: js
    INTEGER, DIMENSION(:),   INTENT(IN)  :: sides
    LOGICAL, DIMENSION(:),   INTENT(IN)  :: Dir
    INTEGER,                 INTENT(OUT) :: j_D_size

    INTEGER, DIMENSION(:), ALLOCATABLE :: j_all, j_dif
    INTEGER :: nj, ms, ns

    !  count all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    !  If the same node belongs to different sides of the boundary
    !  the algorithm takes it into account only once.

    nj = 0
    DO ms = 1, SIZE(js, 2)
       IF ( Dir(sides(ms)) )  nj = nj + SIZE(js, 1)
    ENDDO

    ALLOCATE (j_all(nj), j_dif(nj))

    !  collect all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    nj = 0
    DO ms = 1, SIZE(js, 2)

       IF ( Dir(sides(ms)) ) THEN

          DO ns = 1, SIZE(js, 1)
             nj = nj + 1
             j_all(nj) = js(ns, ms)
          ENDDO

          !TEST***************! problem with version.1 of F90 on CRAY90
          !   PRINT*, ' ms, j_all(ms), js', ms,j_all(ms),js(1,ms),js(2,ms)
          !TEST***************
       ENDIF

    ENDDO


    CALL sort_diff (j_all,  j_dif, j_D_size)

    DEALLOCATE (j_all, j_dif)

  CONTAINS

    SUBROUTINE sort_diff (a,  a_d, n_a_d)

      !  sort in ascending order of the integer array  a  and generation
      !  of the integer array  a_d  whose first  n_a_d  leading entries
      !  contain different values in ascending order, while all the
      !  remaining entries are set to zero

      !  sorting by Shell's method.

      IMPLICIT NONE

      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
      INTEGER,               INTENT(OUT)   :: n_a_d

      INTEGER :: n, na, inc, i, j, k, ia

      na = SIZE(a)

      !  sort phase

      IF (na == 0) THEN
         n_a_d = 0
         RETURN
      ENDIF

      inc = 1
      DO WHILE (inc <= na)
         inc = inc * 3
         inc = inc + 1
      ENDDO

      DO WHILE (inc > 1)
         inc = inc/3
         DO i = inc + 1, na
            ia = a(i)
            j = i
            DO WHILE (a(j-inc) > ia)
               a(j) = a(j-inc)
               j = j - inc
               IF (j <= inc) EXIT
            ENDDO
            a(j) = ia
         ENDDO
      ENDDO

      !  compression phase

      n = 1
      a_d(n) = a(1)
      DO k = 2, na
         IF (a(k) > a(k-1)) THEN
            n = n + 1
            a_d(n) = a(k)
         ENDIF
      ENDDO

      n_a_d = n

      a_d(n_a_d + 1 : na) = 0

    END SUBROUTINE sort_diff

  END SUBROUTINE Dir_nodes_size



  SUBROUTINE Dir_nodes_gen(js, sides, Dir, j_D)

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: js
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides
    LOGICAL, DIMENSION(:),   INTENT(IN) :: Dir
    INTEGER, DIMENSION(:),  INTENT(OUT) :: j_D

    INTEGER, DIMENSION(:), ALLOCATABLE :: j_all, j_dif
    INTEGER :: nj, ms, ns, nj_dif

    !  count all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    !  If the same node belongs to different sides of the boundary
    !  the algorithm takes it into account only once.

    nj = 0
    DO ms = 1, SIZE(js, 2)
       IF ( Dir(sides(ms)) )  nj = nj + SIZE(js, 1)
    ENDDO

    ALLOCATE (j_all(nj), j_dif(nj))

    !  collect all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    nj = 0
    DO ms = 1, SIZE(js, 2)

       IF ( Dir(sides(ms)) ) THEN

          DO ns = 1, SIZE(js, 1)
             nj = nj + 1
             j_all(nj) = js(ns, ms)
          ENDDO

       ENDIF

    ENDDO

    CALL sort_diff (j_all,  j_dif, nj_dif)

    j_D = j_dif(1:nj_dif)

    DEALLOCATE (j_all, j_dif)


  CONTAINS

    SUBROUTINE sort_diff (a,  a_d, n_a_d)

      !  sort in ascending order of the integer array  a  and generation
      !  of the integer array  a_d  whose first  n_a_d  leading entries
      !  contain different values in ascending order, while all the
      !  remaining entries are set to zero

      !  sorting by Shell's method.

      IMPLICIT NONE

      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
      INTEGER,               INTENT(OUT)   :: n_a_d

      INTEGER :: n, na, inc, i, j, k, ia

      na = SIZE(a)

      !  sort phase

      IF (na == 0) THEN
         n_a_d = 0
         RETURN
      ENDIF

      inc = 1
      DO WHILE (inc <= na)
         inc = inc * 3
         inc = inc + 1
      ENDDO

      DO WHILE (inc > 1)
         inc = inc/3
         DO i = inc + 1, na
            ia = a(i)
            j = i
            DO WHILE (a(j-inc) > ia)
               a(j) = a(j-inc)
               j = j - inc
               IF (j <= inc) EXIT
            ENDDO
            a(j) = ia
         ENDDO
      ENDDO

      !  compression phase

      n = 1
      a_d(n) = a(1)
      DO k = 2, na
         IF (a(k) > a(k-1)) THEN
            n = n + 1
            a_d(n) = a(k)
         ENDIF
      ENDDO

      n_a_d = n

      a_d(n_a_d + 1 : na) = 0

    END SUBROUTINE sort_diff

  END SUBROUTINE Dir_nodes_gen


  SUBROUTINE Dir_nodes_sides_gen(js, sides, Dir, j_D, s_D)

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: js
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides
    LOGICAL, DIMENSION(:),   INTENT(IN) :: Dir
    INTEGER, DIMENSION(:),  INTENT(OUT) :: j_D, s_D

    INTEGER, DIMENSION(:), ALLOCATABLE :: j_all, j_dif, s_all
    INTEGER :: nj, ms, ns, nj_dif

    !  count all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    !  If the same node belongs to different sides of the boundary
    !  the algorithm takes it into account only once.

    nj = 0
    DO ms = 1, SIZE(js, 2)
       IF ( Dir(sides(ms)) )  nj = nj + SIZE(js, 1)
    ENDDO

    ALLOCATE (j_all(nj), j_dif(nj))

    !  collect all nodes of the boundary elements which belong to
    !  the sides where a Dirichlet boundary condition is prescribed

    nj = 0
    DO ms = 1, SIZE(js, 2)

       IF ( Dir(sides(ms)) ) THEN

          DO ns = 1, SIZE(js, 1)
             nj = nj + 1
             j_all(nj) = js(ns, ms)
          ENDDO

       ENDIF

    ENDDO

    CALL sort_diff (j_all,  j_dif, nj_dif)

    j_D = j_dif(1:nj_dif)

    nj = MAXVAL(j_D)
    ALLOCATE (s_all(nj))
    nj = 0
    DO ms = 1, SIZE(js, 2)

       IF ( Dir(sides(ms)) ) THEN

          DO ns = 1, SIZE(js, 1)
             nj = nj + 1
             s_all(js(ns, ms)) = sides(ms)
          ENDDO

       ENDIF

    ENDDO

    s_D = s_all(j_D(1:nj_dif))

    DEALLOCATE (j_all, j_dif, s_all)

  CONTAINS

    SUBROUTINE sort_diff (a,  a_d, n_a_d)

      !  sort in ascending order of the integer array  a  and generation
      !  of the integer array  a_d  whose first  n_a_d  leading entries
      !  contain different values in ascending order, while all the
      !  remaining entries are set to zero

      !  sorting by Shell's method.

      IMPLICIT NONE

      INTEGER, DIMENSION(:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
      INTEGER,               INTENT(OUT)   :: n_a_d

      INTEGER :: n, na, inc, i, j, k, ia

      na = SIZE(a)

      !  sort phase

      IF (na == 0) THEN
         n_a_d = 0
         RETURN
      ENDIF

      inc = 1
      DO WHILE (inc <= na)
         inc = inc * 3
         inc = inc + 1
      ENDDO

      DO WHILE (inc > 1)
         inc = inc/3
         DO i = inc + 1, na
            ia = a(i)
            j = i
            DO WHILE (a(j-inc) > ia)
               a(j) = a(j-inc)
               j = j - inc
               IF (j <= inc) EXIT
            ENDDO
            a(j) = ia
         ENDDO
      ENDDO

      !  compression phase

      n = 1
      a_d(n) = a(1)
      DO k = 2, na
         IF (a(k) > a(k-1)) THEN
            n = n + 1
            a_d(n) = a(k)
         ENDIF
      ENDDO

      n_a_d = n

      a_d(n_a_d + 1 : na) = 0

    END SUBROUTINE sort_diff

  END SUBROUTINE Dir_nodes_sides_gen

  !-------------------------------------------------------------------------------

  SUBROUTINE dirichlet_nodes(jjs_in, sides_in, dir_in, js_d)

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs_in
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides_in
    LOGICAL, DIMENSION(:),   INTENT(IN) :: dir_in
    INTEGER, DIMENSION(:),   POINTER    :: js_d
    INTEGER:: j_size

    CALL Dir_nodes_size(jjs_in, sides_in, dir_in,  j_size)
    !      IF (ASSOCIATED(js_d)) THEN
    !         DEALLOCATE(js_d)
    !      END IF
    ALLOCATE (js_d(j_size))
    CALL Dir_nodes_gen(jjs_in, sides_in, dir_in, js_d)

  END SUBROUTINE dirichlet_nodes


  SUBROUTINE dirichlet_nodes_sides(jjs_in, sides_in, dir_in, js_d, sides_D)

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs_in
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides_in
    LOGICAL, DIMENSION(:),   INTENT(IN) :: dir_in
    INTEGER, DIMENSION(:),   POINTER    :: js_d, sides_D
    INTEGER:: j_size

    CALL Dir_nodes_size(jjs_in, sides_in, dir_in,  j_size)
    !      IF (ASSOCIATED(js_d)) THEN
    !         DEALLOCATE(js_d)
    !      END IF
    ALLOCATE (js_d(j_size), sides_D(j_size))
    CALL Dir_nodes_sides_gen(jjs_in, sides_in, dir_in, js_d, sides_D)

  END SUBROUTINE dirichlet_nodes_sides

  SUBROUTINE dirichlet_nodes_bloc(jjs_in, sides_in, dir_in, &
       nb_bloc, bloc_size, js_d)

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jjs_in
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides_in
    LOGICAL, DIMENSION(:),   INTENT(IN) :: dir_in
    INTEGER,                 INTENT(IN) :: nb_bloc, bloc_size
    INTEGER, DIMENSION(:),   POINTER    :: js_d_s, js_d
    INTEGER:: j_size, n0, n_D, k

    CALL Dir_nodes_size(jjs_in, sides_in, dir_in,  j_size)
    ALLOCATE (js_d_s(j_size))
    CALL Dir_nodes_gen(jjs_in, sides_in, dir_in, js_d_s)
    n_D = SIZE(js_d_s)
    ALLOCATE(js_D(nb_bloc*n_D))
    DO k = 1, nb_bloc
       n0 = (k-1)*n_D
       js_D(n0+1:n0+n_D) = js_D_s + (k-1)*bloc_size
    END DO
    DEALLOCATE (js_d_s)

  END SUBROUTINE dirichlet_nodes_bloc


END MODULE Dir_nodes
