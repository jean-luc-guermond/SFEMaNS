MODULE fem_s_M

CONTAINS


  SUBROUTINE qs_00_M(mesh, alpha, ia, ja,  a0)
    !=================================================

    !  alpha < w, _ >   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me
       DO l = 1, l_G

          al = alpha * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                !               IF (j >= i) THEN
                x = ww(nj,l) * al * ww(ni,l)
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO
                !               ENDIF

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_00_M

  SUBROUTINE qs_100_M(mesh, vv, ia, ja,  a0)
    !=================================================

    !  < D w, vv _ >   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: vv
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: m, l, ni, nj, i, j, p, k
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: vl

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me
       DO l = 1, l_G

          vl = 0
          DO ni = 1, n_w
             vl = vl + vv(:,jj(ni,m)) * ww(ni,l)
          END DO


          vl = vl * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                x = 0
                DO k = 1, k_d
                   x = x + vl(k) * dw(k, ni, l, m)
                END DO
                x = x * ww(nj,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_100_M

  !------------------------------------------------------------------------------
  SUBROUTINE qs_000_M(mesh, ff, ia, ja,  a0)
    !=================================================

    !  alpha < w, f _ >   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: m, l, ni, nj, n, i, j, p
    REAL(KIND=8) :: al, x ,  ffl

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          ffl =0.d0
          DO n = 1, n_w
             ffl =  ffl + ff(jj(n,m)) * ww(n,l)
          ENDDO

          al = ffl * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                !               IF (j >= i) THEN
                x = ww(nj,l) * al * ww(ni,l)
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO
                !               ENDIF

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_000_M


  SUBROUTINE qs_11_M (mesh, alpha, ia, ja,  a0)
    !=================================================

    !  alpha << (Dw), (D_) >>   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER :: m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          al = alpha * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                !               IF (j >= i) THEN
                x = al * SUM(dw(:,nj,l,m) * dw(:,ni,l,m))
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO
                !               ENDIF

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_11_M

  SUBROUTINE qs_11_bb_p1_M (mesh, alpha, bcd)
    !=================================================

    !  alpha << (Dw_h^H), (D(_)_h^H) >>   ===>   bcd    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: bcd

    INTEGER :: m, m_mother, l, ni, nj, i, j
    REAL(KIND=8) :: al, x, mest

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       m_mother = (m-1)/n_w + 1

       mest = 0
       DO l = 1, l_G
          mest = mest + rj(l,m)
       END DO
       mest = alpha*(n_w*mest)**(1.d0/k_d)

       DO l = 1, l_G

          al = mest * rj(l,m)

          nj = n_w;  j = jj(nj, m)
          ni = n_w;  i = jj(ni, m)

          x = al * SUM(dw(:,nj,l,m) * dw(:,ni,l,m))
          bcd(m_mother,n_w+1) = bcd(m_mother,n_w+1) + x

       ENDDO
    ENDDO

  END SUBROUTINE qs_11_bb_p1_M


  !=====================================================
  SUBROUTINE qs_v_grad_v_grad_w (mesh, gg, ia, ja, a0)
    !=====================================================

    !  << (g.D)w, (g.D)_ >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Declaration des variables globales
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    !--------------------------------------------------------------------------
    ! Declaration des variables locales
    INTEGER                        :: l, k, ni, nj, p, m, i, j
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)   :: gl
    REAL(KIND=8)                   :: x
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)   :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    boucle_mm :  DO m = 1, me

       boucle_l : DO l = 1, l_G

          gl = 0
          boucle_k : DO k = 1, k_d
             boucle_ni : DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO boucle_ni
          ENDDO  boucle_k

          y = 0
          boucle_ni_2 : DO ni = 1, n_w
             boucle_k_2 : DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO boucle_k_2
          END DO  boucle_ni_2

          boucle_ni_3 : DO ni = 1, n_w
             i = jj(ni, m)
             boucle_nj : DO nj = 1, n_w
                j = jj(nj, m)
                x =  rj(l,m) * y(nj) * y(ni)
                boucle_p : DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN
                      a0(p) = a0(p) + x
                      EXIT
                   ENDIF
                ENDDO  boucle_p
             ENDDO  boucle_nj
          ENDDO  boucle_ni_3

       ENDDO  boucle_l

    ENDDO boucle_mm

  END SUBROUTINE qs_v_grad_v_grad_w


  !=====================================================
  SUBROUTINE qs_LS_mass_adv (mesh, alpha, gg, ia, ja, a0)
    !=====================================================

    !  <<alpha w + (g.D)w, alpha_ +(g.D)_ >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Declaration des variables globales
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    !--------------------------------------------------------------------------
    ! Declaration des variables locales
    INTEGER                        :: l, k, ni, nj, p, m, i, j
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)   :: gl
    REAL(KIND=8)                   :: x
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)   :: y

    !--------------------------------------------------------------------------

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    boucle_mm :  DO m = 1, me

       boucle_l : DO l = 1, l_G

          gl = 0
          boucle_k : DO k = 1, k_d
             boucle_ni : DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO boucle_ni
          ENDDO  boucle_k

          y = alpha*ww(:,l)
          boucle_ni_2 : DO ni = 1, n_w
             boucle_k_2 : DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO boucle_k_2
          END DO  boucle_ni_2

          boucle_ni_3 : DO ni = 1, n_w
             i = jj(ni, m)
             boucle_nj : DO nj = 1, n_w
                j = jj(nj, m)
                x =  rj(l,m) * y(nj) * y(ni)
                boucle_p : DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN
                      a0(p) = a0(p) + x
                      EXIT
                   ENDIF
                ENDDO  boucle_p
             ENDDO  boucle_nj
          ENDDO  boucle_ni_3

       ENDDO  boucle_l

    ENDDO boucle_mm

  END SUBROUTINE qs_LS_mass_adv


  !=====================================================
  SUBROUTINE qs_GALS_mass_adv_M (mesh, param, alpha, gg, ia, ja, a0)
    !=====================================================

    !  <<w + h(alpha w + (g.D)w), alpha_ +(g.D)_ >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Declaration des variables globales
    REAL(KIND=8),                 INTENT(IN)    :: alpha, param
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    !--------------------------------------------------------------------------
    ! Declaration des variables locales
    INTEGER                        :: l, k, ni, nj, p, m, i, j
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)   :: gl
    REAL(KIND=8)                   :: x, mest, hloc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)   :: y

    !--------------------------------------------------------------------------

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       mest = 0
       DO l = 1, l_G
          mest = mest + rj(l,m)
       END DO
       hloc = param*mest**(1.d0/k_d)

       DO l = 1, l_G

          gl = 0
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = alpha*ww(:,l)
          DO ni = 1, n_w
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO

          DO ni = 1, n_w
             i = jj(ni, m)
             DO nj = 1, n_w
                j = jj(nj, m)
                x =  rj(l,m) * y(nj) * (ww(ni,l) + hloc*y(ni))
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN
                      a0(p) = a0(p) + x
                      EXIT
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ENDDO

    ENDDO

  END SUBROUTINE qs_GALS_mass_adv_M

  !=====================================================
  SUBROUTINE qs_h_v_grad_v_grad_w (mesh, gg, ia, ja, a0)
    !=====================================================

    !  << h (g.D)w, (g.D)_ >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Declaration des variables globales
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    !--------------------------------------------------------------------------
    ! Declaration des variables locales
    INTEGER                        :: l, k, ni, nj, p, m, i, j
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d)   :: gl
    REAL(KIND=8)                   :: x, mest, hloc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)   :: y

    !--------------------------------------------------------------------------

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    boucle_mm :  DO m = 1, me

       mest = 0
       DO l = 1, l_G
          mest = mest + rj(l,m)
       END DO
       hloc = mest**(1.d0/k_d)

       boucle_l : DO l = 1, l_G

          gl = 0
          boucle_k : DO k = 1, k_d
             boucle_ni : DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO boucle_ni
          ENDDO  boucle_k

          y = 0
          boucle_ni_2 : DO ni = 1, n_w
             boucle_k_2 : DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO boucle_k_2
          END DO  boucle_ni_2

          boucle_ni_3 : DO ni = 1, n_w
             i = jj(ni, m)
             boucle_nj : DO nj = 1, n_w
                j = jj(nj, m)
                x =  rj(l,m) * y(nj) * y(ni)
                boucle_p : DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN
                      a0(p) = a0(p) + x*hloc
                      EXIT
                   ENDIF
                ENDDO  boucle_p
             ENDDO  boucle_nj
          ENDDO  boucle_ni_3

       ENDDO  boucle_l

    ENDDO boucle_mm

  END SUBROUTINE qs_h_v_grad_v_grad_w


  SUBROUTINE qs_dif_mass_adv_M(mesh, visco, alpha, gg, ia, ja, a0)
    !==========================================================

    !  << visco (Dw), D_) >>
    !  +  alpha * < w, _ >
    !  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2
    !  ===>  a0

    USE Gauss_points

    REAL(KIND=8),                 INTENT(IN)    :: visco, alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: dg, xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       aij = 0
       DO l = 1, l_G

          dg = 0.
          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) * dw(k,ni,l,m)
             END DO
          ENDDO
          viscolm = visco*rj(l,m)
          masslm  = (alpha+0.5*dg)*rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + viscolm*xij + masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l)

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT;
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_dif_mass_adv_M

  SUBROUTINE qs_diff_mass_adv_M(mesh, visco, alpha, gg, ia, ja, a0)
    !==========================================================

    !  << visco (Dw), D_) >>
    !  +  alpha * < w, _ >
    !  +  < w, (g.D)_ >
    !  ===>  a0

    USE gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)

    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       aij = 0
       DO l = 1, l_G

          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO
          viscolm = visco*rj(l,m)
          masslm  = alpha*rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                aij(ni,nj) =  aij(ni,nj) + viscolm*xij + masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l)

             ENDDO
          ENDDO

       ENDDO

       DO ni = 1, n_w; i = jj(ni, m)
          DO nj = 1, n_w;  j = jj(nj, m)
             DO p = ia(i),  ia(i+1) - 1
                IF (ja(p) == j) THEN;  a0(p) = a0(p) + aij(ni,nj);  EXIT;
                ENDIF
             ENDDO
          END DO
       END DO

    ENDDO

  END SUBROUTINE qs_diff_mass_adv_M

  SUBROUTINE qs_dif_mass_adv_skew_M(mesh, visco, alpha, gg, ia, ja, a0)
    !==========================================================

    !  << visco (Dw), D_) >>
    !  +  alpha * < w, _ >
    !  +  < w, (g.D)_ >/2  -  < _, (g.D) w > / 2
    !  ===>  a0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: visco, alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: x, xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO l = 1, l_G

          viscolm = visco*rj(l,m)
          masslm  = alpha*rj(l,m)

          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = 0.5*rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = viscolm*xij + masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l) - ww(nj,l)*y(ni)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_dif_mass_adv_skew_M

  SUBROUTINE qs_dif_mass_M(mesh, visco, alpha, ia, ja, a0)
    !==========================================================

    !  << visco (Dw), D_) >>
    !  +  alpha * < w, _ >
    !  ===>  a0

    USE Gauss_points

    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: visco, alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8) :: x, xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO l = 1, l_G

          viscolm = visco*rj(l,m)
          masslm  = alpha*rj(l,m)

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = viscolm*xij + masslm*wwprod(ni,nj,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_dif_mass_M

  SUBROUTINE qs_adv_M_bis(mesh, gg, ia, ja, a0)
    !==========================================================

    !  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2
    !  ===>  a0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: dg, x,  masslm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO l = 1, l_G

          dg = 0.
          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) * dw(k,ni,l,m)
             END DO
          ENDDO
          masslm  = 0.5*dg*rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x = masslm*wwprod(ni,nj,l) + y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO


             ENDDO
          ENDDO

       ENDDO
    ENDDO
  END SUBROUTINE qs_adv_M_bis

  SUBROUTINE qs_adv_M(mesh, gg, ia, ja, a0)
    !==========================================================

    !  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2
    !  ===>  a0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER :: k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: dg, x, masslm, rjlm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    INTEGER,      DIMENSION(mesh%gauss%n_w) :: j_loc

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = 0.5d0*ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO ni = 1, n_w
          j_loc(ni) = jj(ni,m)
       END DO

       DO l = 1, l_G

          rjlm = rj(l,m)
          dg = 0.
          DO k = 1, k_d
             gl(k) = 0.
             DO ni =1 ,n_w
                dw_loc(k,ni) = dw(k,ni,l,m)
                gl(k) = gl(k) + gg(k, j_loc(ni)) * ww(ni,l)
                dg = dg + gg(k, j_loc(ni)) * dw_loc(k,ni)
             END DO
             gl(k) = gl(k)*rjlm
          ENDDO
          masslm  = dg*rjlm

          DO ni = 1, n_w
             y(ni) = 0.
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw_loc(k,ni)
             END DO
          END DO

          DO ni = 1, n_w; i = j_loc(ni)
             DO nj = 1, n_w;  j = j_loc(nj)

                x = masslm*wwprod(ni,nj,l) + y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO


             ENDDO
          ENDDO

       ENDDO
    ENDDO
  END SUBROUTINE qs_adv_M

  SUBROUTINE qs_001_M(mesh, gg, ia, ja, a0)
    !==========================================================

    !  < w, (g.D)_ >   ===>  a0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x = y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO


             ENDDO
          ENDDO

       ENDDO
    ENDDO
  END SUBROUTINE qs_001_M

  SUBROUTINE qs_h_100_M(mesh, alpha, gg, ia, ja, a0)
    !==========================================================

    !  < H_loc alpha _, (g.D)w >   ===>  a0

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: x, mest, hloc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       mest = 0
       DO l = 1, l_G
          mest = mest + rj(l,m)
       END DO
       hloc = alpha*mest**(1.d0/k_d)

       DO l = 1, l_G

          gl = 0
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
             END DO
          ENDDO

          y = 0
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x = y(ni)*ww(nj,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + hloc*x;  EXIT;
                   ENDIF
                ENDDO


             ENDDO
          ENDDO

       ENDDO
    ENDDO
  END SUBROUTINE qs_h_100_M

  SUBROUTINE qs_1_sgs_1_M (mesh, ff, ia, ja, a0)
    !=======================================

    !  << (Dw), h* f (D_) >>   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: ja
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: a0

    INTEGER :: m, l, ni, nj, i, j, k, p
    REAL(KIND=8) :: fl, x, xij, h, exp

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    m = SIZE(jj,1)
    SELECT CASE(m)
    CASE(3); exp = 1.d0/2
    CASE(6); exp = 1.d0/2
    CASE(4); exp = 1.d0/3
    CASE(10); exp = 1.d0/3
    END SELECT

    DO m = 1, me

       h = 0
       DO l = 1, l_G
          h = h + rj(l,m)
       END DO
       h = h**exp

       DO l = 1, l_G

          fl = 0.
          DO ni =1 ,n_w
             fl = fl + ff(jj(ni,m)) * ww(ni,l)
          END DO
          fl  = h*fl*rj(l,m)

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = fl*xij
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_1_SGS_1_M

  SUBROUTINE qs_101_M (mesh, ff, ia, ja, a0)
    !=======================================

    !  << (Dw), f (D_) >>   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)  :: ja
    REAL(KIND=8), DIMENSION(:), INTENT(INOUT) :: a0

    INTEGER ::  m, l, ni, nj, i, j, k, p
    REAL(KIND=8) :: fl, x, xij

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          fl = 0.
          DO ni =1 ,n_w
             fl = fl + ff(jj(ni,m)) * ww(ni,l)
          END DO
          fl  = fl*rj(l,m)

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = fl*xij
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_101_M

  SUBROUTINE qs_v_plap_M (mesh, p, vstar, ff, ia, ja,  a0)
    !=======================================

    !  << (Dw), |vstar.D(ff)|^(p-2) (D_) >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: p
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: vstar
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    REAL(KIND=8) :: x, xrjlm, xij
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: vl, gl
    INTEGER ::  ni, nj, i, j, k, l, m, q

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          vl = 0.
          gl = 0.
          DO k= 1, k_d
             DO ni =1 ,n_w
                vl(k) = vl(k) + vstar(k,jj(ni,m)) * ww(ni,l)
                gl(k) = gl(k) + ff(jj(ni,m)) * dw(k,ni,l,m)
             END DO
          END DO

          x = 0.
          DO k= 1, k_d
             x = x + vl(k)*gl(k)
          END DO
          IF (ABS(x) .LT. 1.d-12) THEN
             xrjlm = 0.
          ELSE
             xrjlm = (ABS(x)**(p-2.d0))* rj(l,m)
          ENDIF

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w; j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = xrjlm * xij
                DO q = ia(i),  ia(i+1) - 1
                   IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_v_plap_M

  SUBROUTINE qs_plap_M (mesh, p, ff, ia, ja,  a0)
    !=======================================

    !  << (Dw), |D(ff)|^(p-2) (D_) >> ===> A0 incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: p
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    REAL(KIND=8) :: x, xrjlm, xij
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    INTEGER ::  ni, nj, i, j, k, l, m, q

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          gl = 0.
          DO k= 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + ff(jj(ni,m)) * dw(k,ni,l,m)
             END DO
          END DO

          x = 0.
          DO k= 1, k_d
             x = x + gl(k)*gl(k)
          END DO
          IF (ABS(x) .LT. 1.d-20) THEN
             xrjlm = 0.
          ELSE
             xrjlm = (SQRT(x)**(p-2.d0))* rj(l,m)
          ENDIF

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w; j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = xrjlm * xij
                DO q = ia(i),  ia(i+1) - 1
                   IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_plap_M

  SUBROUTINE qs_adv_stab_M(mesh, gg, stab, istab, ia, ja, a0)
    !==========================================================

    !  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2
    !if (istab .eq. 1)   +  c*h^k << (Dw), |D(gg)|^p (D_) >>
    !if (istab .eq. 2)   +  c*h^k << (gg.Dw), |gg.D(gg)|^p (gg.D_) >>
    !  ===>  a0
    !

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: stab
    INTEGER,                      INTENT(IN)    :: istab
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, kp, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: vl
    REAL(KIND=8) :: dg, x, xij, masslm, rjlm
    REAL(KIND=8) :: h,hmu,xrjlm, exponent
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%k_d) :: gradv
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: vdw
    INTEGER,      DIMENSION(mesh%gauss%n_w) :: j_loc

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    exponent = 1./k_d

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = 0.5*ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me

       DO ni = 1, n_w
          j_loc(ni) = jj(ni,m)
       END DO

       h = 0.
       DO l = 1, l_G
          h = h + rj(l,m)
       END DO
       !      h = h**exponent
       hmu = stab(1)*h**(stab(2)*exponent)

       DO l = 1, l_G

          rjlm = rj(l,m)
          dg = 0.
          gradv = 0.

          dw_loc = dw(:,:,l,m)

          DO k = 1, k_d
             vl(k) = 0.
             DO ni =1 ,n_w
                vl(k) = vl(k) + gg(k, j_loc(ni)) * ww(ni,l)
                DO kp = 1,  k_d
                   gradv(k,kp) = gradv(k,kp) + gg(k, j_loc(ni)) * dw_loc(kp,ni)
                END DO
             END DO
             dg = dg + gradv(k,k)
          ENDDO

          masslm  = dg*rjlm

          x = 0.
          IF (istab .EQ. 1) THEN
             DO k= 1, k_d
                DO kp = 1, k_d
                   x = x + (gradv(k,kp) + gradv(kp,k))**2
                END DO
             END DO
          ELSE IF (istab .EQ. 2) THEN
             DO k= 1, k_d
                DO kp = 1, k_d
                   x = x +(vl(kp)*gradv(k,kp))**2
                END DO
             END DO
          ELSE
             WRITE(*,*) 'qs_adv_stab_M: istab > 2'
             STOP
          END IF

          IF (ABS(x) .LT. 1.d-12) THEN
             xrjlm = 0.
          ELSE
             xrjlm = hmu*(SQRT(x)**(stab(3)))* rjlm
          ENDIF

          DO ni = 1, n_w
             vdw(ni) = 0.
             DO k = 1, k_d
                vdw(ni) = vdw(ni) + vl(k)*dw_loc(k,ni)
             END DO
          END DO

          DO ni = 1, n_w; i = j_loc(ni)
             DO nj = 1, n_w;  j = j_loc(nj)


                xij = 0.
                IF (istab .EQ. 1) THEN
                   DO k = 1, k_d
                      xij =  xij + dw_loc(k,nj) * dw_loc(k,ni)
                   END DO
                ELSE IF (istab .EQ. 2) THEN
                   xij =  vdw(ni)*vdw(nj)
                END IF

                x = masslm*wwprod(ni,nj,l) + vdw(nj)*ww(ni,l)*rjlm + xrjlm * xij
                !               x =  x + xrjlm * xij

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO


             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_adv_stab_M

  SUBROUTINE bs_101_M (mesh, ff,  ia, ja, a0)
    !==============================================

    !  - < J(w,_), f >   ===>   a0    incremental accumulation
    !  < (Dw) x k, ff (D_) >   ===>   a0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    REAL(KIND=8) :: f, x, s
    INTEGER ::  m, l, n, n1, i, j, p

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me

       DO l = 1, l_G

          f = SUM(ff(jj(:,m)) * ww(:,l)) * rj(l,m)

          DO n = 1, n_w;  i = jj(n, m)
             DO n1 = 1, n_w;  j = jj(n1, m)
                s = dw(1,n,l,m) * dw(2,n1,l,m) - dw(2,n,l,m) * dw(1,n1,l,m)
                x = -s * f

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                END DO

             ENDDO
          ENDDO

       ENDDO
    ENDDO


  END SUBROUTINE bs_101_M


  !==============================================

  SUBROUTINE qs_00_s_M (mesh, fs, ia, ja, a0)
    !==========================================

    !  < ws, fs _ >_s   ===>   A0    incremental accumulation of boundary terms

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: fs
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  ms, ls, ns, ns1, i, j, q
    REAL(KIND=8) :: fls, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: js, is
    INTEGER,                      POINTER       :: mes

    CALL gauss(mesh)
    js => mesh%jjs
    is => mesh%iis
    mes => mesh%mes

    DO ms = 1, mes

       DO ls = 1, l_Gs
          fls = 0
          DO ns = 1,  n_ws
             fls = fls + fs(is(ns,ms)) * wws(ns,ls) * rjs(ls,ms)
          END DO
          fls = fls * rjs(ls,ms)

          DO ns = 1, n_ws;  i = js(ns, ms)
             DO ns1 = 1, n_ws;  j = js(ns1, ms)
                x = wws(ns,ls) * fls * wws(ns1,ls)

                DO q = ia(i),  ia(i+1) - 1
                   IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE qs_00_s_M

  SUBROUTINE cv_11_M (mesh, alpha, ia, ja, a0)
    !=======================================

    !  alpha << (D x w), (D x _) >>  ===>   a0   incremental acc. in 3D

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  m, l, n, q, n0, k, k1, k2, h, h1, h2, i, j, i_b, j_b, &
         bloc_size
    REAL(KIND=8) :: fla, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    bloc_size = (SIZE(ia) - 1)/3

    DO m = 1, me
       DO l = 1, l_G

          fla = alpha * rj(l,m)

          DO n = 1, n_w;  i_b = jj(n, m)
             DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
                i = i_b + (k-1)*bloc_size

                DO n0 = 1, n_w;  j_b = jj(n0,m)
                   DO h = 1, k_d;  h1 = MODULO(h,k_d) + 1;  h2 = MODULO(h+1,k_d) + 1
                      j = j_b + (h-1)*bloc_size

                      IF (h.EQ.k) THEN
                         x = dw(k1,n,l,m)*dw(h1,n0,l,m) + &
                              dw(k2,n,l,m)*dw(h2,n0,l,m)
                      ELSE
                         x = - dw(h,n,l,m)*dw(k,n0,l,m)
                      END IF

                      x = x * fla

                      DO q = ia(i),  ia(i+1) - 1
                         IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE cv_11_M

  SUBROUTINE cc_101_M (mesh, gg, ia, ja, a0)
    !=======================================

    !  << (D x w), gg x (D x _) >>  ===>   a0   incremental acc. in 3D

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  m, l, n, q, n0, k, k1, k2, h, h1, h2, i, j, i_b, j_b, &
         bloc_size, s , o
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(3) :: gl

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    bloc_size = (SIZE(ia) - 1)/3

    DO m = 1, me
       DO l = 1, l_G

          DO k = 1, k_d
             gl(k) = SUM(gg(k, jj(:,m)) * ww(:,l)) * rj(l,m)
          ENDDO

          DO n = 1, n_w;  i_b = jj(n, m)
             DO k = 1, k_d;  k1 = MODULO(k,k_d) + 1;  k2 = MODULO(k+1,k_d) + 1
                i = i_b + (k-1)*bloc_size

                DO n0 = 1, n_w;  j_b = jj(n0,m)
                   DO h = 1, k_d;  h1 = MODULO(h,k_d) + 1;  h2 = MODULO(h+1,k_d) + 1
                      j = j_b + (h-1)*bloc_size

                      IF (h == k) THEN
                         x = gl(h) * ( dw(k2,n,l,m) * dw(h1,n0,l,m)  &
                              - dw(k1,n,l,m) * dw(h2,n0,l,m) )
                      ELSE
                         o = 6 - k - h
                         s = 0; IF (k > h) s = 1
                         x = (-1)**(k+h+s) * &
                              ( dw(o,n,l,m) * (gl(h1)*dw(h1,n0,l,m) + gl(h2)*dw(h2,n0,l,m)) &
                              + dw(h,n,l,m) * gl(h) * dw(o,n0,l,m) )
                      END IF

                      x = x * rj(l,m)

                      DO q = ia(i),  ia(i+1) - 1
                         IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE cc_101_M


  SUBROUTINE cv_00_s_M (mesh, fs, ia, ja, a0)
    !==========================================

    !  < ws x ns, fs (_ x ns) >_s   ===>   A0    incremental accumulation of boundary terms

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: fs
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  ms, ls, ns, ns1, i, j, q, i_b, j_b, bloc_size, h, k
    REAL(KIND=8) :: fls, x, wwprod
    REAL(KIND=8), DIMENSION(3,3) :: tab

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: js, is
    INTEGER,                      POINTER       :: mes

    CALL gauss(mesh)
    js => mesh%jjs
    is => mesh%iis
    mes => mesh%mes

    bloc_size = (SIZE(ia) - 1)/3

    DO ms = 1, mes
       DO ls = 1, l_Gs

          fls = SUM(fs(is(:,ms)) * wws(:,ls)) * rjs(ls,ms)

          tab(1,1) = rnorms(2,ls,ms)**2 + rnorms(3,ls,ms)**2
          tab(2,2) = rnorms(1,ls,ms)**2 + rnorms(3,ls,ms)**2
          tab(3,3) = rnorms(1,ls,ms)**2 + rnorms(2,ls,ms)**2

          tab(1,2) = -rnorms(2,ls,ms) * rnorms(1,ls,ms)
          tab(1,3) = -rnorms(3,ls,ms) * rnorms(1,ls,ms)
          tab(2,1) = -rnorms(1,ls,ms) * rnorms(2,ls,ms)
          tab(2,3) = -rnorms(3,ls,ms) * rnorms(2,ls,ms)
          tab(3,1) = -rnorms(1,ls,ms) * rnorms(3,ls,ms)
          tab(3,2) = -rnorms(2,ls,ms) * rnorms(3,ls,ms)


          DO ns = 1, n_ws;  i_b = js(ns, ms)
             DO ns1 = 1, n_ws;  j_b = js(ns1, ms)

                wwprod = wws(ns, ls) * wws(ns1, ls) * fls

                DO k = 1, k_d
                   DO h = 1, k_d

                      i = i_b + (k-1)*bloc_size;  j = j_b + (h-1)*bloc_size

                      x = wwprod * tab(k,h)

                      DO q = ia(i),  ia(i+1) - 1
                         IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE cv_00_s_M

  SUBROUTINE cc_00_s_M (mesh, alpha, ia, ja, a0)
    !==========================================

    !  < ws x ns, (_ x ns) >_s   ===>   A0    incremental accumulation of boundary terms

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  ms, ls, ns, ns1, i, j, q, i_b, j_b, bloc_size, h, k
    REAL(KIND=8) :: x, wwprod
    REAL(KIND=8), DIMENSION(3,3) :: tab

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: js
    INTEGER,                      POINTER       :: mes

    CALL gauss(mesh)
    js => mesh%jjs
    mes => mesh%mes

    bloc_size = (SIZE(ia) - 1)/3

    DO ms = 1, mes
       DO ls = 1, l_Gs

          tab(1,1) = rnorms(2,ls,ms)**2 + rnorms(3,ls,ms)**2
          tab(2,2) = rnorms(1,ls,ms)**2 + rnorms(3,ls,ms)**2
          tab(3,3) = rnorms(1,ls,ms)**2 + rnorms(2,ls,ms)**2

          tab(1,2) = -rnorms(2,ls,ms) * rnorms(1,ls,ms)
          tab(1,3) = -rnorms(3,ls,ms) * rnorms(1,ls,ms)
          tab(2,1) = -rnorms(1,ls,ms) * rnorms(2,ls,ms)
          tab(2,3) = -rnorms(3,ls,ms) * rnorms(2,ls,ms)
          tab(3,1) = -rnorms(1,ls,ms) * rnorms(3,ls,ms)
          tab(3,2) = -rnorms(2,ls,ms) * rnorms(3,ls,ms)

          tab = tab * rjs(ls,ms) * alpha

          DO ns = 1, n_ws;  i_b = js(ns, ms)
             DO ns1 = 1, n_ws;  j_b = js(ns1, ms)

                wwprod = wws(ns, ls) * wws(ns1, ls)

                DO k = 1, k_d
                   DO h = 1, k_d

                      i = i_b + (k-1)*bloc_size;  j = j_b + (h-1)*bloc_size

                      x = wwprod * tab(k,h)

                      DO q = ia(i),  ia(i+1) - 1
                         IF (ja(q) == j) THEN;  a0(q) = a0(q) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE cc_00_s_M

  SUBROUTINE qs_dif_massvar_M (mesh, visco, ff, ia, ja,  a0)
    !===========================================================

    !  << visco (Dw), D_) >>
    !  +  < w, ff _ >
    !  ===>  a0   ! incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: visco
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8) :: x, xij, masslm, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO l = 1, l_G

          viscolm = visco*rj(l,m)

          masslm = 0
          DO ni = 1, n_w
             masslm  = masslm + ff(jj(ni,m))*ww(ni,l)
          END DO
          masslm = masslm * rj(l,m)

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                xij = 0.
                DO k = 1, k_d
                   xij =  xij + dw(k,nj,l,m) * dw(k,ni,l,m)
                END DO

                x = viscolm*xij + masslm*wwprod(ni,nj,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_dif_massvar_M

  SUBROUTINE qs_massvar_adv_M (mesh, ff, gg, ia, ja,  a0, a_skew)
    !===================================================================

    !  +  < w, ff _ >
    !  +  < w, (g.D)_ > ! +  skew * < w, (D.g) _ >
    !  ===>  a0   ! incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8), OPTIONAL,       INTENT(IN)    :: a_skew

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: fl, dg, x, masslm, skew=1
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    !   REAL(KIND=8), DIMENSION(k_d,n_w) :: dw_loc

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    IF (PRESENT(a_skew)) THEN
       skew =  a_skew
    ELSE
       skew = 0.5d0
    END IF

    DO m = 1, me
       DO l = 1, l_G

          !        dw_loc = dw(:,:,l,m)

          dg = 0.
          gl = 0.
          fl = 0
          DO ni =1 ,n_w
             fl = fl + ff(jj(ni,m)) * ww(ni,l)
             DO k = 1, k_d
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) *  dw(k,ni,l,m)   ! divergence
             END DO
          ENDDO
          masslm  = (fl + skew*dg)*rj(l,m)   ! skew form
          !         masslm  = fl * rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y


          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x =  masslm*wwprod(ni,nj,l) + y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_massvar_adv_M

  SUBROUTINE qs_varden_adv_M (mesh,  mass, ff, gg, ia, ja,  a0, a_skew)
    !===================================================================

    !  +  < w, ff _ >
    !  +  < w, ( ff*g.D)_ > ! +  skew * < w, ff*(D.g) _ >
    !  ===>  a0   ! incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: mass
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)    :: ff
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8), OPTIONAL,       INTENT(IN)    :: a_skew

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: fl, dg, x, masslm, skew=1
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y


    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    IF (PRESENT(a_skew)) THEN
       skew =  a_skew
    ELSE
       skew = 0.25d0
    END IF

    DO m = 1, me
       DO l = 1, l_G

          dg = 0.
          gl = 0.
          fl = 0
          DO ni =1 ,n_w
             fl = fl + ff(jj(ni,m)) * ww(ni,l)
             DO k = 1, k_d
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) *  dw(k,ni,l,m)   ! divergence
             END DO
          ENDDO
          masslm  = (mass*fl + fl*skew*dg)*rj(l,m)   ! skew form

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*fl*y


          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x =  masslm*wwprod(ni,nj,l) + y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_varden_adv_M

  SUBROUTINE qs_mass_adv_M (mesh, alpha, gg, ia, ja,  a0, a_skew)
    !===================================================================

    !  +  alpha * < w, _ >
    !  +  < w, (g.D)_ >  +  < w, (D.g) _ > / 2
    !  ===>  a0   ! incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8), OPTIONAL,       INTENT(IN)    :: a_skew

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: dg, x, masslm, skew
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    IF (PRESENT(a_skew)) THEN
       skew =  a_skew
    ELSE
       skew = 0.5d0
    END IF


    DO m = 1, me
       DO l = 1, l_G


          dg = 0.
          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) * dw(k,ni,l,m)
             END DO
          ENDDO
          masslm  = (alpha+skew*dg)*rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) * dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x =  masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_mass_adv_M


  SUBROUTINE qs_mass_div_M (mesh, alpha, gg, ia, ja,  a0)
    !===================================================================

    !  +  alpha * < w, _ >
    !  +  < w, (g.D)_ >  +  < w, (D.g) _ >
    !  ===>  a0   ! incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)    :: gg
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia, ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    INTEGER ::  k, l, m, ni, nj, i, j, p
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d) :: gl
    REAL(KIND=8) :: dg, x, masslm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w) :: y
    !   REAL(KIND=8), DIMENSION(k_d,n_w) :: dw_loc

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO l = 1, l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = ww(ni,l)*ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, me
       DO l = 1, l_G

          !         dw_loc = dw(:,:,l,m)

          dg = 0.
          gl = 0.
          DO k = 1, k_d
             DO ni =1 ,n_w
                gl(k) = gl(k) + gg(k, jj(ni,m)) * ww(ni,l)
                dg = dg + gg(k, jj(ni,m)) * dw(k,ni,l,m)
             END DO
          ENDDO
          masslm  = (alpha+dg)*rj(l,m)

          y = 0.
          DO ni = 1, n_w;
             DO k = 1, k_d
                y(ni) =  y(ni) + gl(k) *  dw(k,ni,l,m)
             END DO
          END DO
          y = rj(l,m)*y

          DO ni = 1, n_w; i = jj(ni, m)
             DO nj = 1, n_w;  j = jj(nj, m)

                x =  masslm*wwprod(ni,nj,l) + &
                     y(nj)*ww(ni,l)

                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_mass_div_M

  SUBROUTINE elast_M (mesh, alpha, ia, ja, a0)
    !=======================================

    !   << (Dw), (D_) >>
    !  << (Dw), (D_)^t >>
    !  alpha << (D.w), (D._) >>
    !  ===>   a0   incremental acc. in k_d Dimensions

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  m, l, p, ni, nj, k, k1, h, i, j, i_b, j_b, &
         bloc_size
    REAL(KIND=8) :: x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    bloc_size = (SIZE(ia) - 1)/k_d

    DO m = 1, me
       DO l = 1, l_G

          DO ni = 1, n_w;  i_b = jj(ni, m)
             DO k = 1, k_d;
                i = i_b + (k-1)*bloc_size

                DO nj = 1, n_w;  j_b = jj(nj,m)
                   DO h = 1, k_d;
                      j = j_b + (h-1)*bloc_size

                      x =       dw(h,ni,l,m)*dw(k,nj,l,m) &
                           + alpha*dw(k,ni,l,m)*dw(h,nj,l,m)
                      IF (h.EQ.k) THEN
                         DO k1 = 1, k_d
                            x = x + dw(k1,ni,l,m)*dw(k1,nj,l,m)
                         END DO
                      END IF

                      x = x * rj(l,m)

                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE elast_M

  SUBROUTINE curl_div_2D_M (mesh, alpha, stabh, expdiv, ia, ja, a0)
    !=======================================
    !          << w , _ >>
    !  +       << (Dxw), (Dx_) >>
    !  + alpha << (D.w), (D._) >>
    !  ===>   a0   incremental acc. in k_d Dimensions

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: stabh, expdiv
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER ::  m, l, p, ni, nj, k, k1, h, i, j, i_b, j_b, &
         bloc_size
    REAL(KIND=8) :: x, alphah

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    bloc_size = (SIZE(ia) - 1)/2

    DO m = 1, me
       alphah = stabh*SUM(rj(:,m))**(expdiv/2)
       DO l = 1, l_G

          DO ni = 1, n_w;  i_b = jj(ni, m)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                k1 = MODULO(k,2) + 1

                DO nj = 1, n_w;  j_b = jj(nj,m)
                   DO h = 1, 2
                      j = j_b + (h-1)*bloc_size

                      x = alphah*dw(k,ni,l,m)*dw(h,nj,l,m)
                      IF (h==k) THEN
                         x = x + alpha*ww(ni,l)*ww(nj,l) + dw(k1,ni,l,m)*dw(k1,nj,l,m)
                      ELSE
                         x = x - dw(h,ni,l,m)*dw(k,nj,l,m)
                      END IF

                      x = x * rj(l,m)

                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE curl_div_2D_M

  SUBROUTINE curl_grad_2D_M (mesh, alpha, stab, exp_sta, ia, ja, a0)
    !=======================================
    ! alpha(phij, phii) + (Curl phij, Curl phii) + (Grad psij, phii) + sta*h^(2*exp_sta)*(div phij, Div phii)
    ! - (Grad phij, psii) + h^(2(1-exp_sta))*(Grad psij, Grad psii)
    !  ===>   a0   incremental acc. in k_d Dimensions

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: stab, exp_sta
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER ::  m, l, p, ni, nj, k, k1, h, i, j, i_b, j_b, &
         bloc_size, type_fe
    REAL(KIND=8) :: x, alphah, betah, hloc

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    bloc_size = (SIZE(ia) - 1)/3
    IF (mesh%gauss%n_w == 3) THEN
       type_fe = 1
    ELSE IF (mesh%gauss%n_w == 6) THEN
       type_fe = 2
    ELSE
       WRITE(*,*) ' BUG, bad FE'
    END IF

    DO m = 1, me
       hloc = SQRT(SUM(rj(:,m)))/type_fe
       alphah = stab*hloc**(2*exp_sta)
       betah = hloc**(2*(1-exp_sta))*(1.d0/stab)
       DO l = 1, l_G

          DO ni = 1, n_w;  i_b = jj(ni, m)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                k1 = MODULO(k,2) + 1

                DO nj = 1, n_w;  j_b = jj(nj,m)
                   DO h = 1, 2
                      j = j_b + (h-1)*bloc_size

                      x = alphah*dw(k,ni,l,m)*dw(h,nj,l,m)
                      IF (h==k) THEN
                         x = x + alpha*ww(ni,l)*ww(nj,l) + dw(k1,ni,l,m)*dw(k1,nj,l,m)
                      ELSE
                         x = x - dw(h,ni,l,m)*dw(k,nj,l,m)
                      END IF
                      x = x * rj(l,m)
                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   END DO
                ENDDO

             ENDDO
          ENDDO

          DO ni = 1, n_w;  i_b = jj(ni, m)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                DO nj = 1, n_w;  j_b = jj(nj,m)
                   h = 3
                   j = j_b + (h-1)*bloc_size
                   x = dw(k,nj,l,m)*ww(ni,l)* rj(l,m)
                   DO p = ia(i),  ia(i+1) - 1
                      IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;
                         EXIT
                      ENDIF
                   ENDDO

                ENDDO
             ENDDO
          ENDDO

          DO ni = 1, n_w;  i_b = jj(ni, m)
             k = 3
             i = i_b + (k-1)*bloc_size
             DO nj = 1, n_w;  j_b = jj(nj,m)
                DO h = 1, 3
                   j = j_b + (h-1)*bloc_size
                   IF (h<3) THEN
                      x = -ww(nj,l) * dw(h,ni,l,m)* rj(l,m)
                   ELSE
                      x = betah*SUM(dw(:,ni,l,m) * dw(:,nj,l,m))* rj(l,m)
                   END IF
                   DO p = ia(i),  ia(i+1) - 1
                      IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;
                         EXIT
                      ENDIF
                   ENDDO
                END DO
             END DO

          ENDDO
       ENDDO
    END DO
  END SUBROUTINE curl_grad_2D_M

  SUBROUTINE curl_surf_2D_M (mesh, stab_b, exp_b, consist, adj, ia, ja, a0)
    !=======================================
    !         alpha*h << wxn , _xn >>
    !  -       << (Dxw)xn, _ >>
    !  ===>   a0   incremental acc. in k_d Dimensions

    USE Gauss_points

    IMPLICIT NONE
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0
    REAL(KIND=8),                 INTENT(IN)    :: stab_b, consist, adj, exp_b
    INTEGER ::  m, ms, ls, p, ni, nj, k, k1, h, h1, i, j, i_b, j_b, &
         bloc_size
    REAL(KIND=8) :: x, alphah

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jjs
    INTEGER,                      POINTER       :: mes

    CALL gauss(mesh)
    jjs => mesh%jjs
    mes => mesh%mes

    bloc_size = mesh%np

    DO ms = 1, mes
       alphah = stab_b*SUM(rjs(:,ms))**(exp_b)
       m  =mesh%neighs(ms)

       DO ls = 1, l_Gs

          DO ni = 1, n_ws;  i_b = jjs(ni, ms)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                k1 = MODULO(k,2) + 1

                DO nj = 1, n_ws;  j_b = jjs(nj,ms)
                   DO h = 1, 2
                      j = j_b + (h-1)*bloc_size
                      h1 = MODULO(h,2) + 1

                      x = alphah*(-1)**(k1+h1)*wws(ni,ls)*rnorms(k1,ls,ms)*wws(nj,ls)*rnorms(h1,ls,ms)
                      x = x * rjs(ls,ms)

                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

          DO ni = 1, n_ws;  i_b = jjs(ni, ms)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                k1 = MODULO(k,2) + 1

                DO nj = 1, n_w;  j_b = mesh%jj(nj,m)
                   DO h = 1, 2
                      j = j_b + (h-1)*bloc_size
                      h1 = MODULO(h,2) + 1

                      x =consist*(-1)**(k1+h) *wws(ni,ls)*rnorms(k1,ls,ms)*dw_s(h1,nj,ls,ms)
                      x = x * rjs(ls,ms)

                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

          DO ni = 1, n_w;  i_b = mesh%jj(ni,m)
             DO k = 1, 2
                i = i_b + (k-1)*bloc_size
                k1 = MODULO(k,2) + 1

                DO nj = 1, n_ws;  j_b = jjs(nj,ms)
                   DO h = 1, 2
                      j = j_b + (h-1)*bloc_size
                      h1 = MODULO(h,2) + 1

                      x = adj*(-1)**(k+h1) *wws(nj,ls)*rnorms(h1,ls,ms)*dw_s(k1,ni,ls,ms)
                      x = x * rjs(ls,ms)

                      DO p = ia(i),  ia(i+1) - 1
                         IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                         ENDIF
                      ENDDO

                   ENDDO
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    END DO

  END SUBROUTINE curl_surf_2D_M

  SUBROUTINE qs_1x1x_M (mesh, alpha, ia, ja,  a0)
    !=================================================

    !  alpha << (Dw), (D_) >>   ===>   A0    incremental accumulation

    USE Gauss_points

    IMPLICIT NONE
    REAL(KIND=8),                 INTENT(IN)    :: alpha
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: a0

    INTEGER ::  m, l, ni, nj, i, j, p
    REAL(KIND=8) :: al, x

    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER,      DIMENSION(:,:), POINTER       :: jj
    INTEGER,                      POINTER       :: me

    CALL gauss(mesh)
    jj => mesh%jj
    me => mesh%me

    DO m = 1, me
       DO l = 1, l_G

          al = alpha * rj(l,m)

          DO nj = 1, n_w;  j = jj(nj, m)
             DO ni = 1, n_w;  i = jj(ni, m)

                x = al * dw(1,nj,l,m) * dw(1,ni,l,m)
                DO p = ia(i),  ia(i+1) - 1
                   IF (ja(p) == j) THEN;  a0(p) = a0(p) + x;  EXIT;
                   ENDIF
                ENDDO

             ENDDO
          ENDDO

       ENDDO
    ENDDO

  END SUBROUTINE qs_1x1x_M

  SUBROUTINE edge_stab_M(mesh, coeff_visc, ia, ja, aa)
    USE def_type_mesh
    !USE solve_sp
    IMPLICIT NONE
    TYPE(mesh_type)    :: mesh
    !TYPE(csr_matrix)   :: mat
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ia
    INTEGER,      DIMENSION(:),   INTENT(IN)    :: ja
    REAL(KIND=8), DIMENSION(:),   INTENT(INOUT) :: aa
    REAL(KIND=8)       :: coeff_visc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%l_Gs, 2) :: dwni_loc
    INTEGER, DIMENSION(mesh%gauss%n_w,2) :: jji_loc
    INTEGER :: ls, ms, cotei, cotej, p, ni, nj, i, j, type_fe
    REAL(KIND=8) :: x, h2, coeff !=  0.02! 0.003!  .003
    IF (mesh%gauss%n_w==3) THEN
       type_fe = 1
    ELSE
       type_fe = 2
    END IF
    coeff = coeff_visc/type_fe**2
    DO ms = 1, mesh%mi
       dwni_loc = mesh%gauss%dwni(:,:,:,ms)
       jji_loc = mesh%jji(:,:,ms)
       h2 = coeff*SUM(mesh%gauss%rji(:,ms))**2
       DO cotei = 1, 2
          DO ni = 1, mesh%gauss%n_w
             i = jji_loc(ni,cotei)
             DO cotej = 1, 2
                DO nj = 1, mesh%gauss%n_w
                   j = jji_loc(nj,cotej)
                   x = 0.d0
                   DO ls = 1, mesh%gauss%l_Gs
                      x =  x + dwni_loc(ni,ls,cotei)*dwni_loc(nj,ls,cotej)*mesh%gauss%rji(ls,ms)
                   END DO
                   DO p = ia(i),  ia(i+1) - 1
                      IF (ja(p) == j) THEN;
                         aa(p) = aa(p) + x*h2;
                         EXIT
                      ENDIF
                   ENDDO
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE edge_stab_M

END MODULE fem_s_M
