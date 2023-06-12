!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE sub_plot

CONTAINS

  SUBROUTINE plot_vit_2d(jj, rr, uu)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation

    INTEGER :: m, n1, n2, n3

    OPEN(UNIT = 20, FILE ='vit_u', FORM ='formatted', STATUS ='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20,*) '% meshplot  = False'

    DO m = 1, SIZE(jj, 2)
       !      n1 = jj(1, m)
       !      n2 = jj(2, m)
       !      n3 = jj(3, m)
       !      WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       !      WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       !      WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       !      WRITE (20, *)

       n1 = jj(1, m)
       n2 = jj(6, m)
       n3 = jj(5, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)


       n1 = jj(2, m)
       n2 = jj(4, m)
       n3 = jj(6, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)

       n1 = jj(4, m)
       n2 = jj(3, m)
       n3 = jj(5, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)

       n1 = jj(5, m)
       n2 = jj(6, m)
       n3 = jj(4, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO

  END SUBROUTINE plot_vit_2d

  !----------------------------------------------------------------

  SUBROUTINE plot_arrow(jj, rr, vv)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr, vv

    INTEGER :: n
    REAL(KIND=8) :: z_d, w_d

    OPEN(UNIT = 20,FILE ='velocity',FORM ='formatted',STATUS ='unknown')

    WRITE (20, 9000) SIZE(rr,2), SIZE(jj,1), SIZE(jj,2)
9000 FORMAT('#',3i10)
    WRITE (20, *) '$ DATA = VECTOR'
    WRITE (20, *) '% axisscale = FALSE'
    WRITE (20, *) '% vscale = 1.'

    !  WRITE (20, *) '% contstyle = 2'
    !  WRITE (20, *) '% meshplot  = True'

    SELECT CASE(SIZE(rr,1))

    CASE(2)
       z_d = 0
       w_d = 0
       DO n = 1, SIZE(rr, 2)
          WRITE (20,*) rr(1,n), rr(2,n), z_d, vv(1,n), vv(2,n), w_d
       ENDDO

    CASE(3)
       DO n = 1, SIZE(rr, 2)
          WRITE (20,*) rr(:,n), vv(:,n)
       ENDDO

    END SELECT

    CLOSE(20)

  END SUBROUTINE plot_arrow

  !------------------------------------------------------------------------------

  SUBROUTINE plot_pressure_2d(jj, rr, uu)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation

    INTEGER :: na, nb, n, m, n1, n2, n3
    REAL(KIND=8) :: xa, xb, xm, r, p_midpoint

    OPEN (UNIT = 20, FILE = 'pressure_2D.plt', FORM = 'formatted', STATUS = 'unknown')

    !  set pressure field so as to have a zero pressure at
    !  the mid-point of the bottom unit side of the cavity

    !  find the mid-point

    xa = 0;  xb = 1;  xm = 0.5
    na = 1;  nb = 1
    DO n = 1, SIZE(uu)

       IF (rr(2,n) < 1.0d-4) THEN

          r = rr(1,n)
          IF (r < xm  .AND.  r > xa) THEN;  xa = r;  na = n
          ENDIF
          IF (r > xm  .AND.  r < xb) THEN;  xb = r;  nb = n
          ENDIF

       ENDIF

    ENDDO

    p_midpoint = uu(na) + (uu(nb) - uu(na)) * (xm-xa)/(xb-xa)

    !  reset pressure field

    uu = uu - p_midpoint

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 1'
    WRITE (20, *) '% contours = "( -.12  -.06   -.03 )"'
    WRITE (20, *) '% contours = "( -.015 -.0075 -.00375 )"'
    WRITE (20, *) '% contours = "( 0.0 )"'
    WRITE (20, *) '% contours = "(  .12   .06    .03 )"'
    WRITE (20, *) '% contours = "(  .015  .0075  .00375 )"'

    WRITE (20,*) '% meshplot  = False'

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO
    CLOSE(20)

  END SUBROUTINE plot_pressure_2d

  !------------------------------------------------------------------------------

  SUBROUTINE plot_loc_rel_var(jj, rr, vo, vv)
    !---FORMAT PLTMTV

    !  plot local relative variation of speed


    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr, vo, vv

    REAL(KIND=8), DIMENSION(SIZE(vv,2)) :: uu
    REAL(KIND=8) :: vo_n_mod
    INTEGER :: n, m, n1, n2, n3

    OPEN (UNIT=20, FILE ='speed_var',FORM ='formatted',STATUS ='unknown')


    DO n = 1, SIZE(vv,2)

       vo_n_mod = SQRT(SUM(vo(:,n)**2))

       IF (vo_n_mod < 1.0d-15)  THEN

          uu(n) = 0

       ELSE

          uu(n) = SQRT(SUM( (vo(:,n) - vv(:,n))**2 ))/vo_n_mod

       ENDIF

    ENDDO


    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% meshplot  = True'

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO

  END SUBROUTINE plot_loc_rel_var

  !------------------------------------------------------------------------------

  SUBROUTINE plot_pressure(jj, rr, uu)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation
    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE='pressure', FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% cstep = 20'
    WRITE (20, *) '% meshplot  = False'

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, *) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, *) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, *) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO

  END SUBROUTINE plot_pressure

  SUBROUTINE plot_const_p1_label(jj, rr, uu, file_name )
    !---FORMAT PLTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: uu
    INTEGER ::  m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = true'
    WRITE (20, *)

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(m), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(m), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(m), n3
       WRITE (20, *)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)
    CLOSE(20)

  END SUBROUTINE plot_const_p1_label

  SUBROUTINE plot_pressure_p1_label(jj, rr, uu, file_name )
    !---FORMAT PLTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation
    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = false'
    WRITE (20, *)

    DO m = 1, SIZE(jj, 2)

       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)
    CLOSE(20)

  END SUBROUTINE plot_pressure_p1_label

  SUBROUTINE plot_p1_cont_label(jj, jjs, sides, list, rr, uu, file_name )
    !---FORMAT PLOTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj, jjs
    INTEGER, DIMENSION(:),   INTENT(IN) :: sides, list
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation

    INTEGER :: m, n1, n2, n3, ms

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = f'
    WRITE (20, *)

    DO ms = 1, SIZE(sides)
       IF(MINVAL(ABS(list-sides(ms))) /= 0 ) CYCLE
       n1 = jjs(1,ms)
       n2 = jjs(2,ms)
       WRITE(20,110) '@line x1=',rr(1,n1), 'y1=',rr(2,n1), 'z1=',0.,  &
            'x2=',rr(1,n2), 'y2=',rr(2,n2), 'z2=',0.
    END DO
110 FORMAT(6(A,x,e12.5,3x))

    WRITE (20, *)
    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)
    CLOSE(20)

  END SUBROUTINE plot_p1_cont_label

  SUBROUTINE plot_p1_matiere_label(jj, neigh, i_d, rr, uu, file_name )
    !---FORMAT PLOTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj, neigh
    INTEGER, DIMENSION(:),   INTENT(IN) :: i_d
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation

    INTEGER :: n, m, n1, n2, n3, nghm

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = f'
    WRITE (20, *)

    DO m = 1, SIZE(neigh,2)
       DO n = 1, 3
          nghm = neigh(n,m)

          IF (nghm == 0) THEN
             CONTINUE
          ELSE IF (i_d(nghm) /= i_d(m)) THEN
             CONTINUE
          ELSE
             CYCLE
          END IF

          n1 = jj(MODULO(n,3) + 1,m)
          n2 = jj(MODULO(n+1,3) + 1,m)
          WRITE(20,110) '@line x1=',rr(1,n1), 'y1=',rr(2,n1), 'z1=',0.,  &
               'x2=',rr(1,n2), 'y2=',rr(2,n2), 'z2=',0.
       END DO
    END DO
110 FORMAT(6(A,x,e12.5,3x))

    WRITE (20, *)
    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)
    CLOSE(20)

  END SUBROUTINE plot_p1_matiere_label

  SUBROUTINE plot_pressure_label(jj, rr, uu, file_name )
    !---FORMAT PLTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER, DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation

    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = False'
    WRITE (20, *)

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, *)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)

    CLOSE(20)

  END SUBROUTINE plot_pressure_label

  SUBROUTINE plot_arrow_label(jj, rr, vv, file_name)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    CHARACTER(*) :: file_name
    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr, vv

    INTEGER :: n
    REAL(KIND=8) :: z_d, w_d

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, 9000) SIZE(rr,2), SIZE(jj,1), SIZE(jj,2)
9000 FORMAT('#',3i10)
    WRITE (20, *) '$ DATA = VECTOR'
    WRITE (20, *) '% axisscale = FALSE'
    WRITE (20, *) '% vscale = 1.'

    !  WRITE (20, *) '% contstyle = 2'
    !  WRITE (20, *) '% meshplot  = True'

    SELECT CASE(SIZE(rr,1))

    CASE(2)
       z_d = 0
       w_d = 0
       IF (SIZE(vv,1)/=2) THEN
          DO n = 1, SIZE(rr, 2)
             WRITE (20,100) rr(1,n), rr(2,n), z_d, vv(n,1), vv(n,2), w_d
          ENDDO
       ELSE
          DO n = 1, SIZE(rr, 2)
             WRITE (20,100) rr(1,n), rr(2,n), z_d, vv(1,n), vv(2,n), w_d
          ENDDO
       END IF
    CASE(3)
       IF (SIZE(vv,1)/=3) THEN
          DO n = 1, SIZE(rr, 2)
             WRITE (20,100) rr(:,n), vv(:,n)
          ENDDO
       ELSE
          DO n = 1, SIZE(rr, 2)
             WRITE (20,100) rr(:,n), vv(:,n)
          ENDDO
       END IF
    END SELECT

    CLOSE(20)
100 FORMAT(6(e12.5,3x))

  END SUBROUTINE plot_arrow_label

  SUBROUTINE plot_pressure_P2_label(jj, rr, uu, file_name)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: uu  ! rigid translation
    CHARACTER(*) :: file_name

    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20,*) '% meshplot  = false'

    DO m = 1, SIZE(jj, 2)

       n1 = jj(1, m)
       n2 = jj(6, m)
       n3 = jj(5, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, 100)


       n1 = jj(2, m)
       n2 = jj(4, m)
       n3 = jj(6, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, 100)

       n1 = jj(4, m)
       n2 = jj(3, m)
       n3 = jj(5, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, 100)

       n1 = jj(5, m)
       n2 = jj(6, m)
       n3 = jj(4, m)
       WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
       WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
       WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
       WRITE (20, 100)
    ENDDO
100 FORMAT(3(e12.5,3x),i5)

    CLOSE(20)

  END SUBROUTINE plot_pressure_p2_label

  !----------------------------------------------------------------

  SUBROUTINE plot_ENSIGHT_vecteur(u8, vit)
    !--FORMAT ENSIGHT

    IMPLICIT NONE

    CHARACTER(LEN=*),             INTENT(IN) :: vit
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: u8
    CHARACTER(LEN=80)                        :: sketuve
    INTEGER                                  :: i, j

    OPEN(UNIT=33,FILE=vit,STATUS='unknown',FORM='UNFORMATTED')

    sketuve = 'rien'
    WRITE(33) sketuve

    WRITE(33) ((REAL(u8(i,j),KIND=4),i=1,SIZE(u8,1)),j=1,SIZE(u8,2))

    CLOSE(33)

  END SUBROUTINE plot_ENSIGHT_vecteur

  !----------------------------------------------------------------

  SUBROUTINE plot_ENSIGHT_scalaire(p8, pres)
    !--FORMAT ENSIGHT

    IMPLICIT NONE

    CHARACTER(LEN=*),             INTENT(IN) :: pres
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: p8
    CHARACTER(LEN=80)                        :: sketuve
    INTEGER                                  :: i

    OPEN(UNIT=44,FILE=pres,STATUS='unknown',FORM='UNFORMATTED')

    sketuve = 'rien'
    WRITE(44) sketuve

    WRITE(44) (REAL(p8(i),KIND =4),i=1,SIZE(p8))

    CLOSE(44)

  END SUBROUTINE plot_ENSIGHT_scalaire

  SUBROUTINE  plot_vorticity (jj, rr, zz, i)

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: zz
    INTEGER,                      INTENT(IN) :: i

    INTEGER :: m, n1, n2, n3

    !   OPEN (UNIT = 20, FILE = 'vorticity', FORM = 'formatted', STATUS = 'unknown')

    WRITE (i, *) '$ DATA = CONTCURVE'
    WRITE (i, *) '% contstyle = 1'
    WRITE (i, *) '% nsteps = 40'

    !   WRITE (i, *) '% contours = "(-6.0 -5.0 -4.0 -3.0 -2.0 -1.0)"'
    !   WRITE (i, *) '% contours = "(-0.5  0.0  0.5)"'
    !   WRITE (i, *) '% contours = "( 1.0  2.0  3.0  4.0  5.0  6.0)"'

    WRITE (i, *) '% meshplot  = False'

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (i, *) rr(1,n1), rr(2,n1), zz(n1), n1
       WRITE (i, *) rr(1,n2), rr(2,n2), zz(n2), n2
       WRITE (i, *) rr(1,n3), rr(2,n3), zz(n3), n3
       WRITE (i, *)
    ENDDO

    CLOSE (UNIT = i)

  END SUBROUTINE  plot_vorticity

  SUBROUTINE  plot_stream (jj, rr, pp, i)

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: pp
    INTEGER,                      INTENT(IN) :: i

    INTEGER :: m, n1, n2, n3

    !   OPEN (UNIT = 20, FILE = 'stream', FORM = 'formatted', STATUS = 'unknown')

    WRITE (i, *) '$ DATA = CONTCURVE'
    WRITE (i, *) '% contstyle = 1'
    WRITE (i, *) '% nsteps = 30'

    !   WRITE (i, *) '% contours = "(-1.e-10  -1.e-7  -1.e-5  -1.e-4 )"'
    !   WRITE (i, *) '% contours = "(-0.01  -0.03  -0.05  -0.07  -0.09 )"'
    !   WRITE (i, *) '% contours = "(-0.100  -0.110  -0.115  -0.1175 )"'
    !   WRITE (i, *) '% contours = "(+1.0e-8 +1.0e-7 +1.0e-6 +1.0e-5 )"'
    !   WRITE (i, *) '% contours = "(+5.0e-5 +1.0e-4 +2.5e-4 +5.0e-4 )"'
    !   WRITE (i, *) '% contours = "(+1.0e-3 +1.5e-3 +3.0e-3 )"'
    WRITE (i, *) '% meshplot  = False'

    DO m = 1, SIZE(jj, 2)
       n1 = jj(1, m)
       n2 = jj(2, m)
       n3 = jj(3, m)
       WRITE (i, *) rr(1,n1), rr(2,n1), pp(n1), n1
       WRITE (i, *) rr(1,n2), rr(2,n2), pp(n2), n2
       WRITE (i, *) rr(1,n3), rr(2,n3), pp(n3), n3
       WRITE (i, *)
    ENDDO

    CLOSE (UNIT = i)

  END SUBROUTINE  plot_stream

  SUBROUTINE plot_scalar_field(jj, rr, uu, file_name)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    REAL(KIND=8), DIMENSION(:),   INTENT(IN) :: uu  ! rigid translation
    CHARACTER(*) :: file_name

    INTEGER :: m, n1, n2, n3, unit_w=47

100 FORMAT(3(e15.8,3x),i5)
    OPEN (UNIT=unit_w, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (unit_w, '(A)') '$ DATA = CONTCURVE'
    WRITE (unit_w, '(A)') '% contstyle = 2'
    WRITE (unit_w, '(A)') '% nsteps = 50'
    WRITE (unit_w, '(A)') '% meshplot  = false'
    WRITE (unit_w, '(A)')

    IF (SIZE(jj,1)==3) THEN
       DO m = 1, SIZE(jj, 2)
          n1 = jj(1, m)
          n2 = jj(2, m)
          n3 = jj(3, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)
       ENDDO
    ELSE IF (SIZE(jj,1)==6) THEN
       DO m = 1, SIZE(jj, 2)
          n1 = jj(1, m)
          n2 = jj(6, m)
          n3 = jj(5, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(2, m)
          n2 = jj(4, m)
          n3 = jj(6, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(4, m)
          n2 = jj(3, m)
          n3 = jj(5, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(5, m)
          n2 = jj(6, m)
          n3 = jj(4, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)
       ENDDO

    ELSE IF (SIZE(jj,1)==10) THEN

       DO m = 1, SIZE(jj, 2)
          n1 = jj(1, m)
          n2 = jj(8, m)
          n3 = jj(6, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(8, m)
          n2 = jj(10, m)
          n3 = jj(6, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(8, m)
          n2 = jj(9, m)
          n3 = jj(10, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(9, m)
          n2 = jj(4, m)
          n3 = jj(10, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(9, m)
          n2 = jj(2, m)
          n3 = jj(4, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(6, m)
          n2 = jj(10, m)
          n3 = jj(7, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(10, m)
          n2 = jj(5, m)
          n3 = jj(7, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(10, m)
          n2 = jj(4, m)
          n3 = jj(5, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)

          n1 = jj(7, m)
          n2 = jj(5, m)
          n3 = jj(3, m)
          WRITE (unit_w, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (unit_w, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (unit_w, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (unit_w, 100)
       ENDDO
    ELSE
       WRITE(*,*) ' Problem in plot_scalar_field ', SIZE(jj,1)
       STOP
    END IF

    CLOSE(UNIT=unit_w)

  END SUBROUTINE plot_scalar_field

  SUBROUTINE plot_two_scalar_field(jj, rr, uu, jj2, rr2, uu2, file_name)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj, jj2
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr, rr2
    REAL(KIND=8), DIMENSION(:)               :: uu, uu2
    CHARACTER(*) :: file_name

    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = false'
    WRITE (20, *)

    IF (SIZE(jj,1)==3) THEN
       DO m = 1, SIZE(jj, 2)

          n1 = jj(1, m)
          n2 = jj(2, m)
          n3 = jj(3, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (20, *)
       ENDDO
       DO m = 1, SIZE(jj2, 2)

          n1 = jj2(1, m)
          n2 = jj2(2, m)
          n3 = jj2(3, m)
          WRITE (20, 100) rr2(1,n1), rr2(2,n1), uu2(n1), n1
          WRITE (20, 100) rr2(1,n2), rr2(2,n2), uu2(n2), n2
          WRITE (20, 100) rr2(1,n3), rr2(2,n3), uu2(n3), n3
          WRITE (20, *)
       ENDDO


    ELSE IF (SIZE(jj,1)==6) THEN

       DO m = 1, SIZE(jj, 2)

          n1 = jj(1, m)
          n2 = jj(6, m)
          n3 = jj(5, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (20, 100)


          n1 = jj(2, m)
          n2 = jj(4, m)
          n3 = jj(6, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (20, 100)

          n1 = jj(4, m)
          n2 = jj(3, m)
          n3 = jj(5, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (20, 100)

          n1 = jj(5, m)
          n2 = jj(6, m)
          n3 = jj(4, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), uu(n1), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), uu(n2), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), uu(n3), n3
          WRITE (20, 100)
       ENDDO
    ELSE IF (SIZE(jj,1)==10) THEN
       DO m = 1, SIZE(jj2, 2)

          n1 = jj2(1, m)
          n2 = jj2(6, m)
          n3 = jj2(5, m)
          WRITE (20, 100) rr2(1,n1), rr2(2,n1), uu2(n1), n1
          WRITE (20, 100) rr2(1,n2), rr2(2,n2), uu2(n2), n2
          WRITE (20, 100) rr2(1,n3), rr2(2,n3), uu2(n3), n3
          WRITE (20, 100)


          n1 = jj2(2, m)
          n2 = jj2(4, m)
          n3 = jj2(6, m)
          WRITE (20, 100) rr2(1,n1), rr2(2,n1), uu2(n1), n1
          WRITE (20, 100) rr2(1,n2), rr2(2,n2), uu2(n2), n2
          WRITE (20, 100) rr2(1,n3), rr2(2,n3), uu2(n3), n3
          WRITE (20, 100)

          n1 = jj2(4, m)
          n2 = jj2(3, m)
          n3 = jj2(5, m)
          WRITE (20, 100) rr2(1,n1), rr2(2,n1), uu2(n1), n1
          WRITE (20, 100) rr2(1,n2), rr2(2,n2), uu2(n2), n2
          WRITE (20, 100) rr2(1,n3), rr2(2,n3), uu2(n3), n3
          WRITE (20, 100)

          n1 = jj2(5, m)
          n2 = jj2(6, m)
          n3 = jj2(4, m)
          WRITE (20, 100) rr2(1,n1), rr2(2,n1), uu2(n1), n1
          WRITE (20, 100) rr2(1,n2), rr2(2,n2), uu2(n2), n2
          WRITE (20, 100) rr2(1,n3), rr2(2,n3), uu2(n3), n3
          WRITE (20, 100)
       ENDDO
    ELSE
       WRITE(*,*) ' Problem in plot_two_scalar_field ', SIZE(jj,1)
       STOP
    END IF

100 FORMAT(3(e12.5,3x),i5)

    CLOSE(20)

  END SUBROUTINE plot_two_scalar_field

  SUBROUTINE plot_scalar_field_domain(jj, rr, id, index_dom, file_name)
    !---FORMAT PLTMTV

    IMPLICIT NONE

    INTEGER,      DIMENSION(:,:), INTENT(IN) :: jj
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: rr
    INTEGER,      DIMENSION(:),   INTENT(IN) :: id
    INTEGER,                      INTENT(IN) :: index_dom
    CHARACTER(*) :: file_name

    INTEGER :: m, n1, n2, n3

    OPEN (UNIT=20, FILE=file_name, FORM='formatted', STATUS='unknown')

    WRITE (20, *) '$ DATA = CONTCURVE'
    WRITE (20, *) '% contstyle = 2'
    WRITE (20, *) '% nsteps = 50'
    WRITE (20, *) '% meshplot  = true'
    WRITE (20, *)

    IF (SIZE(jj,1)==3) THEN
       DO m = 1, SIZE(jj, 2)
          IF (id(m)/=index_dom) CYCLE
          n1 = jj(1, m)
          n2 = jj(2, m)
          n3 = jj(3, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), id(m), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), id(m), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), id(m), n3
          WRITE (20, *)
       ENDDO

    ELSE IF (SIZE(jj,1)==6) THEN

       DO m = 1, SIZE(jj, 2)
          IF (id(m)/=index_dom) CYCLE
          n1 = jj(1, m)
          n2 = jj(6, m)
          n3 = jj(5, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), id(m), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), id(m), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), id(m), n3
          WRITE (20, 100)

          n1 = jj(2, m)
          n2 = jj(4, m)
          n3 = jj(6, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), id(m), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), id(m), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), id(m), n3
          WRITE (20, 100)

          n1 = jj(4, m)
          n2 = jj(3, m)
          n3 = jj(5, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), id(m), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), id(m), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), id(m), n3
          WRITE (20, 100)

          n1 = jj(5, m)
          n2 = jj(6, m)
          n3 = jj(4, m)
          WRITE (20, 100) rr(1,n1), rr(2,n1), id(m), n1
          WRITE (20, 100) rr(1,n2), rr(2,n2), id(m), n2
          WRITE (20, 100) rr(1,n3), rr(2,n3), id(m), n3
          WRITE (20, 100)
       ENDDO

    ELSE
       WRITE(*,*) ' Problem in plot_scalar_field_domain ', SIZE(jj,1)
       STOP
    END IF

100 FORMAT(2(e12.5,3x),i5,3x,i5)

    CLOSE(20)

  END SUBROUTINE plot_scalar_field_domain

  SUBROUTINE vtk_p1_2d(mesh, champ, unit_file, file_name)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: champ
    INTEGER,                                    INTENT(IN) :: unit_file
    CHARACTER(*),                               INTENT(IN) :: file_name
    INTEGER :: i

    write(*,*) 'ecriture header'
    OPEN (UNIT=unit_file,FILE=file_name,FORM = 'formatted',&
         STATUS = 'unknown')

    write(unit_file,'(A)') '# vtk DataFile Version 3.0'
    write(unit_file,'(A)') 'vtk '//file_name//''
    write(unit_file,'(A)')'ASCII'
    write(unit_file,'(A)')'DATASET UNSTRUCTURED_GRID'
    write(unit_file,'(A,I7,A)')'POINTS ', mesh%np, ' float'
    write(*,*) 'points ...'
    DO i=1, mesh%np
       write(unit_file,'(2(e14.7,2x),A)') mesh%rr(1,i), &
            mesh%rr(2,i), ' 0.0 '
    ENDDO
    write(*,*) 'cells ...'
    write(unit_file,'(A,I7,I8)') 'CELLS ', mesh%me, mesh%me*4
    DO i=1, mesh%me
       write(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(1,i)-1,  &
            mesh%jj(2,i)-1, mesh%jj(3,i)-1
    ENDDO
    write(unit_file,'(A,I7)') 'CELL_TYPES ', mesh%me
    DO i=1, mesh%me
       write(unit_file,'(A)') '5'
    ENDDO

    write(*,*) 'data ...'
    write(unit_file,'(A,I7)') 'POINT_DATA ',mesh%np
    write(unit_file,'(A)') 'SCALARS scalars float 1'
    write(unit_file,'(A)') 'LOOKUP_TABLE default'
    DO i=1, mesh%np
       write(unit_file,'(e14.7,2x)') champ(i)
    ENDDO
    CLOSE(unit_file)
  END SUBROUTINE vtk_p1_2d

  SUBROUTINE vtk_2d(mesh, champ, unit_file, file_name)

    USE def_type_mesh

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: champ
    INTEGER,                                    INTENT(IN) :: unit_file
    CHARACTER(*),                               INTENT(IN) :: file_name
    INTEGER :: i, vtk_cell, n, m


    write(*,*) 'ecriture header'
    OPEN (UNIT=unit_file,FILE=file_name,FORM = 'formatted',&
         STATUS = 'unknown')

    write(unit_file,'(A)') '# vtk DataFile Version 3.0'
    write(unit_file,'(A)') 'vtk '//file_name//''
    write(unit_file,'(A)')'ASCII'
    write(unit_file,'(A)')'DATASET UNSTRUCTURED_GRID'
    write(unit_file,'(A,I7,A)')'POINTS ', mesh%np, ' float'
    write(*,*) 'points ...'
    DO i=1, mesh%np
       write(unit_file,'(2(e14.7,2x),A)') mesh%rr(1,i), &
            mesh%rr(2,i), ' 0.0 '
    ENDDO
    write(*,*) 'cells ...'

    write(unit_file,'(A,I7,I8)') 'CELLS ', mesh%me, mesh%me*(mesh%gauss%n_w+1)
    IF (mesh%gauss%n_w==3) THEN
       vtK_cell = 5
       DO m=1, mesh%me
          write(unit_file,'(I2,6(I8,1x))') mesh%gauss%n_w,  (mesh%jj(n,m)-1, n=1,mesh%gauss%n_w)
       ENDDO
    ELSE IF (mesh%gauss%n_w==6) THEN
       vtK_cell = 22
       DO m=1, mesh%me
          write(unit_file,'(I2,6(I8,1x))') mesh%gauss%n_w,  &
               mesh%jj(1,m)-1,  mesh%jj(2,m)-1, mesh%jj(3,m)-1,  mesh%jj(6,m)-1, mesh%jj(4,m)-1,  mesh%jj(5,m)-1
       ENDDO
    ELSE
       WRITE(*,*) ' BUG  in  vtk_2d '
    END IF
    write(unit_file,'(A,I7)') 'CELL_TYPES ', mesh%me
    DO i=1, mesh%me
       write(unit_file,'(i2)')  vtK_cell
    ENDDO

    write(*,*) 'data ...'
    write(unit_file,'(A,I7)') 'POINT_DATA ',mesh%np
    write(unit_file,'(A)') 'SCALARS scalars float 1'
    write(unit_file,'(A)') 'LOOKUP_TABLE default'
    DO i=1, mesh%np
       write(unit_file,'(e14.7,2x)') champ(i)
    ENDDO
    CLOSE(unit_file)


    write(*,*) 'ecriture header'
    OPEN (UNIT=unit_file,FILE=file_name,FORM = 'formatted',&
         STATUS = 'unknown')

!!$    write(unit_file,'(A)') '# vtk DataFile Version 3.0'
!!$    write(unit_file,'(A)') 'vtk '//file_name//''
!!$    write(unit_file,'(A)')'ASCII'
!!$    write(unit_file,'(A)')'DATASET UNSTRUCTURED_GRID'
!!$    write(unit_file,'(A,I7,A)')'POINTS ', mesh%np, ' float'
!!$    write(*,*) 'points ...'
!!$    DO i = 1, mesh%np
!!$       write(unit_file,'(2(e14.7,2x),A)') mesh%rr(1,i), &
!!$            mesh%rr(2,i), ' 0.0 '
!!$    ENDDO
!!$    write(*,*) 'cells ...'
!!$    write(unit_file,'(A,I7,I8)') 'CELLS ', mesh%me, mesh%me*(mesh%gauss%n_w+1)
!!$    DO m = 1, mesh%me
!!$       write(unit_file,'(A,3(I8,1x))') mesh%gauss%n_w,  (mesh%jj(n,m)-1, n=1,mesh%gauss%n_w)
!!$    ENDDO
!!$    IF (mesh%gauss%n_w==3) THEN
!!$       vtK_cell = 5
!!$    ELSE IF (mesh%gauss%n_w==6) THEN
!!$       vtK_cell = 22
!!$    ELSE
!!$       WRITE(*,*) ' BUG  in  vtk_2d '
!!$    END IF
!!$    write(unit_file,'(A,I7)') 'CELL_TYPES ', mesh%me
!!$    DO m = 1, mesh%me
!!$       write(unit_file,'(I2)') vtk_cell
!!$    ENDDO
!!$
!!$    write(*,*) 'data ...'
!!$    write(unit_file,'(A,I7)') 'POINT_DATA ',mesh%np
!!$    write(unit_file,'(A)') 'SCALARS scalars float 1'
!!$    write(unit_file,'(A)') 'LOOKUP_TABLE default'
!!$    DO i=1, mesh%np
!!$       write(unit_file,'(e14.7,2x)') champ(i)
!!$    ENDDO
!!$    CLOSE(unit_file)
  END SUBROUTINE vtk_2d

  SUBROUTINE trace_profile(mesh, v, it, freq_plot, list_mode, nom_champ, num_dom)

    USE chaine_caractere
    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type)                                 :: mesh
    REAL(KIND=8)   , DIMENSION(:,:,:), INTENT(INOUT):: v
    INTEGER,                           INTENT(IN)   :: it, freq_plot
    INTEGER, DIMENSION(:),             INTENT(IN)   :: list_mode
    CHARACTER(len=3),                  INTENT(IN)   :: nom_champ
    INTEGER, OPTIONAL,                 INTENT(IN)   :: num_dom

    INTEGER                           :: l, lblank, rang_S
    CHARACTER(len=3)                  :: tit, tmode, tit_S
    CHARACTER(len=1)                  :: tcomp
    INTEGER                           :: i, k

    IF (PRESENT(num_dom)) THEN
       rang_S = num_dom
    ELSE
       rang_S = 0
    END IF

    WRITE(tit_S,'(i3)') rang_S
    lblank = eval_blank(3,tit_S)
    DO l = 1, lblank - 1
       tit_S(l:l) = '0'
    END DO

    WRITE(tit,'(i3)') it/freq_plot
    lblank = eval_blank(3,tit)
    DO l = 1, lblank - 1
       tit(l:l) = '0'
    END DO

    DO i= 1, SIZE(v,3)

       WRITE(tmode,'(i3)') list_mode(i)
       lblank = eval_blank(3,tmode)
       DO l = 1, lblank - 1
          tmode(l:l) = '0'
       END DO

       DO k= 1, size(v,2)
          WRITE(tcomp,'(i1)') k
          CALL plot_scalar_field(mesh%jj, mesh%rr, &
               v(:,k,i) , nom_champ//tcomp//'_m='//tmode//'_'//tit_S//'_'//tit//'.plt')
       ENDDO

    END DO

  END SUBROUTINE trace_profile

END MODULE sub_plot
