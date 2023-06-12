!
!Authors: Jean-Luc Guermond, Raphael Laguerre, Luigi Quartapelle, Copyrights 1996, 2000, 2004
!
MODULE  fem_tn_axi
CONTAINS


  SUBROUTINE dot_product (mesh, ff, gg,  t)
    !===============================

    !  sqrt(< f^2 >)   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff, gg
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l ,i ,ni

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0

    DO m = 1, me
       DO l = 1, l_G

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          t = t + SUM(ff(jj(:,m)) * ww(:,l))*SUM(gg(jj(:,m)) * ww(:,l))*rj(l,m)*ray

       ENDDO
    ENDDO

  END SUBROUTINE dot_product


  SUBROUTINE ns_0 (mesh, ff,  t)
    !===============================

    !  sqrt(< f^2 >)   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l ,i ,ni

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0

    DO m = 1, me
       DO l = 1, l_G

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          t = t + SUM(ff(jj(:,m)) * ww(:,l))**2 * rj(l,m)*ray

       ENDDO
    ENDDO

    t = SQRT(t)

  END SUBROUTINE ns_0

  SUBROUTINE ns_l1 (mesh, ff,  t)
    !===============================

    !  < f >   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l ,i ,ni

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0

    DO m = 1, me
       DO l = 1, l_G

          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          t = t + SUM(ff(jj(:,m)) * ww(:,l))* rj(l,m)*ray

       ENDDO
    ENDDO

  END SUBROUTINE ns_l1

  SUBROUTINE average (mesh, ff,  t)
    !===============================

    !  < f >/vol   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l ,i ,ni

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray, vol


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0
    vol = 0.d0
    IF (SIZE(ff)==mesh%np) THEN
       DO m = 1, me
          DO l = 1, l_G

             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + rr(1,i)*ww(ni,l)
             END DO

             t = t + SUM(ff(jj(:,m)) * ww(:,l))* rj(l,m)*ray
             vol = vol + rj(l,m)*ray
          ENDDO
       ENDDO
    ELSE IF (SIZE(ff)==mesh%me) THEN
       DO m = 1, me
          DO l = 1, l_G

             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + rr(1,i)*ww(ni,l)
             END DO

             t = t + ff(m)* rj(l,m)*ray
             vol = vol + rj(l,m)*ray
          ENDDO
       ENDDO
    ELSE
       WRITE(*,*) ' BUG in average '
       STOP
    END IF
    t =  t / vol

  END SUBROUTINE average

  SUBROUTINE ns_anal_0 (mesh, ff, ff_anal, t)
    !===============================

    !  sqrt(< f^2 >)   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff_anal
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l ,i , ni, index

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray, fl

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr
    t = 0

    DO m = 1, me
       index = (m-1)*l_G
       DO l = 1, l_G; index = index + 1

          !===Compute radius of Gauss point
          ray = 0
          fl = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
             fl = fl + ff(i) * ww(ni,l)
          END DO

          t = t + (fl - ff_anal(index))**2 * ray* rj(l,m)

       ENDDO
    ENDDO

    t = SQRT(t)

  END SUBROUTINE ns_anal_0

  SUBROUTINE ns_1 (mesh, ff,  t)
    !===============================

    !  sqrt(<< (Df).(Df) >>)   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    REAL(KIND=8),                 INTENT(OUT) :: t

    INTEGER ::  m, l, k,i,ni
    REAL(KIND=8) :: s

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr

    t = 0

    DO m = 1, me
       DO l = 1, l_G
          s = 0
          DO k = 1, k_d

             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + rr(1,i)*ww(ni,l)
             END DO

             s = s + SUM(ff(jj(:,m)) * dw(k,:,l,m))**2*ray

          ENDDO

          t = t + s * rj(l,m)

       ENDDO
    ENDDO

    t = SQRT(t)

  END SUBROUTINE ns_1

  SUBROUTINE nv_1(mesh, mod_max, ff, t)
    !semi-norme vectoriel H1 pour une representation de Fourier
    !(v_rc, v_rs, v_thc, v_ths, v_zc, v_zs)
    !sqrt(<< (Dv).(Dv) >>)   ===>   t

    USE Gauss_points

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,                        INTENT(IN)  :: mod_max
    REAL(KIND=8), DIMENSION(6,mesh%np,0:mod_max), INTENT(IN)  :: ff
    REAL(KIND=8),                   INTENT(OUT) :: t

    INTEGER ::  m, l,i,ni ,j

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray
    REAL(KIND=8), DIMENSION(2)                  :: div

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr

    t = 0.d0
    div = 0.d0

    DO j = 0, mod_max
       DO m = 1, me
          DO l = 1, l_G
             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, n_w;  i = jj(ni,m)
                ray = ray + rr(1,i)*ww(ni,l)
             END DO

             div(1) = (SUM(ff(1,jj(:,m),j) * dw(1,:,l,m)) + &
                  SUM(ff(1,jj(:,m),j) * ww(:,l))/ray + &
                  j/ray * SUM(ff(4,jj(:,m),j) * ww(:,l)) + &
                  SUM(ff(5,jj(:,m),j) * dw(2,:,l,m)))* ray*rj(l,m)

             IF (j > 0) THEN
                div(2) = (SUM(ff(2,jj(:,m),j) * dw(1,:,l,m)) + &
                     SUM(ff(2,jj(:,m),j) * ww(:,l))/ray - &
                     j/ray * SUM(ff(3,jj(:,m),j) * ww(:,l)) + &
                     SUM(ff(6,jj(:,m),j) * dw(2,:,l,m)))* ray*rj(l,m)
             ENDIF

             t = t +(div(1)**2+div(2)**2)

          ENDDO
       ENDDO

    ENDDO

    t = SQRT(t)


  END SUBROUTINE nv_1

  SUBROUTINE nv_0_cn(mesh, ff, p, t)
    !semi-norme vectoriel H1 pour une representation de Fourier
    !(v_rc, v_rs, v_zc, v_zs)

    USE Gauss_points

    IMPLICIT NONE

    TYPE(mesh_type), TARGET                     :: mesh
    REAL(KIND=8), DIMENSION(mesh%np,6), INTENT(IN)  :: ff
    REAL(KIND=8),                   INTENT(OUT) :: t, p

    INTEGER ::  m, l, i, ni
    REAL(KIND=8) :: s, rp, rt

    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me
    REAL(KIND=8) ,DIMENSION(:,:), POINTER       :: rr
    REAL(KIND=8)                                :: ray

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me
    rr=> mesh%rr

    rp = 0.d0
    rt = 0.d0
    s = 0.d0
    DO m = 1, me
       DO l = 1, l_G
          !===Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj(ni,m)
             ray = ray + rr(1,i)*ww(ni,l)
          END DO

          rp = rp + sqrt(SUM(ff(jj(:,m),1)* ww(:,l))**2 + SUM(ff(jj(:,m),5)* ww(:,l))**2)*ray*rj(l,m)
          rt = rt + ABS(SUM(ff(jj(:,m),3)* ww(:,l)))*ray*rj(l,m)
          s = s + ray*rj(l,m)

       ENDDO
    ENDDO

    p = rp/s
    t = rt/s

  END SUBROUTINE nv_0_cn

END MODULE fem_tn_axi
