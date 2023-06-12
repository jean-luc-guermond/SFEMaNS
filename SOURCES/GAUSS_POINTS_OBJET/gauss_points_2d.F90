!===
!Author: Jean-Luc Guermond, Copyright December 24th 2018
!===
MODULE mod_gauss_points_2d
  PRIVATE
  PUBLIC gauss_points_2d
  INTEGER, PARAMETER :: k_d=2
CONTAINS
  SUBROUTINE gauss_points_2d(mesh,type_fe)
    USE def_type_mesh
    USE GP_2d_p1
    USE GP_2d_p2
    USE GP_2d_p3
    USE sub_plot
    IMPLICIT NONE
    TYPE(mesh_type),                  TARGET    :: mesh
    INTEGER                                     :: type_fe
    INTEGER,                          POINTER   :: me, mes
    INTEGER,      DIMENSION(:,:),     POINTER   :: jj
    INTEGER,      DIMENSION(:,:),     POINTER   :: js
    REAL(KIND=8), DIMENSION(:,:),     POINTER   :: rr
    REAL(KIND=8), DIMENSION(:,:),     POINTER   :: ww
    REAL(KIND=8), DIMENSION(:,:),     POINTER   :: wws
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER   :: dw
    REAL(KIND=8), DIMENSION(:,:,:),   POINTER   :: rnorms
    !REAL(KIND=8), DIMENSION(:,    :, :),POINTER :: rnorms_v !normal at vertices
    REAL(KIND=8), DIMENSION(:,:),     POINTER   :: rj
    REAL(KIND=8), DIMENSION(:,:),     POINTER   :: rjs
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER   :: dw_s
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER   :: dwps
    REAL(KIND=8), DIMENSION(:,:,:,:), POINTER   :: dws
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dd
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: dds
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dds_v
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: pp
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: pps
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: r
    REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: rs
    REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: ww_s
    INTEGER,      DIMENSION(:),     ALLOCATABLE :: ngauss
    REAL(KIND=8), DIMENSION(k_d, k_d)           :: dr
    REAL(KIND=8), DIMENSION(1,   k_d)           :: drs
    REAL(KIND=8) :: rjac, rjacs, x
    INTEGER      :: m, l, k, k1, n, n1, n2, ms, ns, ls, face, cote, orient
    INTEGER      :: n_w, n_ws, l_G, l_Gs
    REAL(KIND=8), DIMENSION(k_d) :: rnor, rsd

    SELECT CASE(type_fe)
    CASE(1)
       n_w=3;  n_ws=2; l_G=3; l_Gs=2
    CASE(2)
       n_w=6;  n_ws=3; l_G=7; l_Gs=3
    CASE(3)
       ! n_w=10; n_ws=4; l_G=7; l_Gs=4
       n_w=10; n_ws=4; l_G=12; l_Gs=4
    CASE DEFAULT
       WRITE(*,*) ' FE not programmed yet', type_fe
       STOP
    END SELECT
    ALLOCATE(dd(k_d,n_w,l_G),dds(1,n_ws,l_Gs),pp(l_G),pps(l_Gs),&
         r(n_w),rs(n_ws),ww_s(n_w,l_Gs),dds_v(n_ws,n_ws),ngauss(l_Gs))
    me => mesh%me
    mes => mesh%mes
    jj => mesh%jj
    js => mesh%jjs
    rr => mesh%rr
    IF (mesh%edge_stab) THEN
       ALLOCATE(mesh%gauss%dwni(n_w, l_Gs, 2, mesh%mi))
       ALLOCATE(mesh%gauss%rji(l_Gs, mesh%mi))
       ALLOCATE(mesh%gauss%rnormsi(k_d, l_Gs, 2, mesh%mi))
       ALLOCATE(mesh%gauss%wwsi(n_ws, l_gs))
    END IF
    ALLOCATE(mesh%gauss%ww(n_w, l_g))
    ALLOCATE(mesh%gauss%wws(n_ws, l_gs))
    ALLOCATE(mesh%gauss%dw(k_d, n_w, l_G,  me ))
    ALLOCATE(mesh%gauss%rnorms(k_d,  l_Gs, mes))
    ALLOCATE(mesh%gauss%rnorms_v(k_d,n_ws, mes))
    ALLOCATE(mesh%gauss%rj(l_G, me ))
    ALLOCATE(mesh%gauss%rjs(l_Gs, mes))
    ALLOCATE(mesh%gauss%dw_s(k_d, n_w,  l_Gs,  mesh%mes))
    ALLOCATE(mesh%gauss%dwps(1, n_ws, l_Gs, mes))
    ALLOCATE(mesh%gauss%dws(1, n_ws, l_Gs, mes))
    ww => mesh%gauss%ww
    wws => mesh%gauss%wws
    dw => mesh%gauss%dw
    rnorms => mesh%gauss%rnorms
    rj => mesh%gauss%rj
    rjs => mesh%gauss%rjs
    dw_s => mesh%gauss%dw_s
    dwps => mesh%gauss%dwps
    dws => mesh%gauss%dws
    mesh%gauss%k_d  = k_d
    mesh%gauss%n_w  = n_w
    mesh%gauss%l_G  = l_G
    mesh%gauss%n_ws = n_ws
    mesh%gauss%l_Gs = l_Gs
    !===Evaluate and store the values of derivatives and of the
    !===Jacobian determinant at Gauss points of all volume elements
    SELECT CASE(type_fe)
    CASE(1)
       CALL element_2d_p1(ww, dd, pp, n_w, l_G)
    CASE(2)
       CALL element_2d_p2(ww, dd, pp, n_w, l_G)
    CASE(3)
       CALL element_2d_p3(ww, dd, pp, n_w, l_G)
    END SELECT
    DO m = 1, me
       DO l = 1, l_G
          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,l))
             ENDDO
          ENDDO
          rjac = dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1)
          DO n = 1, n_w
             dw(1, n, l, m)   &
                  = (+ dd(1,n,l)*dr(2,2) - dd(2,n,l)*dr(2,1))/rjac
             dw(2, n, l, m)   &
                  = (- dd(1,n,l)*dr(1,2) + dd(2,n,l)*dr(1,1))/rjac
          ENDDO
          rj(l, m) = ABS(rjac)*pp(l) !===Sign of Jacobian determinant unknown
       ENDDO
    ENDDO

    SELECT CASE(type_fe)
    CASE(1)
       CALL element_1d_p1(wws, dds, pps, n_ws, l_Gs)
    CASE(2)
       CALL element_1d_p2(wws, dds, pps, n_ws, l_Gs)
    CASE(3)
       CALL element_1d_p3(wws, dds, pps, n_ws, l_Gs)
    END SELECT
    !===surface elements
    DO ms = 1, mes
       m = mesh%neighs(ms)
       !===Determine which face of the reference element is associated with ms
       DO n = 1, 3
          IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
          face = n
       END DO
       SELECT CASE(type_fe)
       CASE(1)
          CALL element_2d_p1_boundary(face, dd, ww_s, n_w, l_Gs)
       CASE(2)
          CALL element_2d_p2_boundary(face, dd, ww_s, n_w, l_Gs)
       CASE(3)
          CALL element_2d_p3_boundary(face, dd, ww_s, n_w, l_Gs)
       END SELECT
       !===Fix for meshes not oriented correctly
       n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
       IF (js(1,ms)==jj(n1,m) .AND. js(2,ms)==jj(n2,m)) THEN
          orient = 1
       ELSE
          orient = -1
       END IF
       IF (face==2) THEN
          orient = -orient
       END IF
       !=========TEST
       !orient = 1
       !=========TEST
       DO ls = 1, l_Gs
          ngauss(ls) = orient*ls+(l_Gs+1)*(1-orient)/2 !===rearrangement of gauss points
       END DO
       !===End fix
       DO ls = 1, l_Gs
          DO k = 1, k_d
             r = rr(k, jj(:,m))
             DO k1 = 1, k_d
                dr(k, k1) = SUM(r * dd(k1,:,ls))
             ENDDO
          ENDDO
          rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))
          DO n = 1, n_w
             dw_s(1,n,ngauss(ls),ms) = (+dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac
             dw_s(2,n,ngauss(ls),ms) = (-dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
          ENDDO
       ENDDO
    ENDDO
    DO ms = 1, mes
       DO ls = 1, l_Gs
          DO k = 1, k_d
             rs = rr(k, js(:,ms))
             drs(1, k) = SUM(rs*dds(1,:,ls))
          ENDDO
          rjacs = SQRT(drs(1,1)**2 + drs(1,2)**2)
          rnorms(1, ls, ms) = - drs(1,2)/rjacs
          rnorms(2, ls, ms) = + drs(1,1)/rjacs
          rjs(ls, ms) = rjacs * pps(ls)
          m = mesh%neighs(ms)
          DO n = 1, n_w
             IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
             face = n
          END DO
          rsd(1:k_d) = rr(:,jj(face,m)) - (rr(:,js(1,ms))+rr(:,js(2,ms)))/2
          x = SUM(rnorms(:,ls,ms)*rsd(1:k_d))
          IF (x>0) THEN
             rnorms(:,ls,ms) = -rnorms(:,ls,ms)
          END IF
       ENDDO
    ENDDO
    DO ms = 1, mes  !===Evaluation of tangent gradient. Obsolete
       DO ls = 1, l_Gs
          dwps(1,:,ls,ms) = dds(1,:,ls)*pps(ls)
          dws(1,:,ls,ms)  = dds(1,:,ls)
       ENDDO
    ENDDO

    !===Normal vector at Lagrange nodes
    SELECT CASE(type_fe)
    CASE(1)
       CALL element_1d_p1_at_nodes (dds_v, n_ws)
       DO ms = 1, mes
          DO ns = 1, n_ws
             DO k = 1, k_d
                rs = rr(k, js(:,ms))
                drs(1, k) = SUM(rs * dds_v(:,ns))
             ENDDO
             rjacs = SQRT( drs(1,1)**2 + drs(1,2)**2 )
             mesh%gauss%rnorms_v(1, ns, ms) = -drs(1,2)/rjacs
             mesh%gauss%rnorms_v(2, ns, ms) = drs(1,1)/rjacs
             m = mesh%neighs(ms)
             !===Find correct face and orient normal outwards
             DO n = 1, n_w
                IF (MINVAL(ABS(js(:,ms)-jj(n,m)))==0) CYCLE
                face = n
             END DO
             rs = rr(:,jj(face,m)) - (rr(:,js(1,ms))+rr(:,js(2,ms)))/2
             x = SUM(mesh%gauss%rnorms_v(:,ns,ms)*rs)
             IF (x>0) THEN
                mesh%gauss%rnorms_v(:,ns,ms) = - mesh%gauss%rnorms_v(:,ns,ms)
             END IF
          END DO
       ENDDO
    CASE DEFAULT
       !Not programmed yet
       mesh%gauss%rnorms_v = 0.d0
       !WRITE(*,*) 'BUG gauss_points_2d, rnorms_v not programmed yet'
    END SELECT


    !===Cell interface (JLG, April 2009)
    IF (mesh%edge_stab) THEN
       SELECT CASE(type_fe)
       CASE(1)
          CALL element_1d_p1(mesh%gauss%wwsi, dds, pps, n_ws, l_Gs)
       CASE(2)
          CALL element_1d_p2(mesh%gauss%wwsi, dds, pps, n_ws, l_Gs)
       CASE(3)
          CALL element_1d_p3(mesh%gauss%wwsi, dds, pps, n_ws, l_Gs)
       END SELECT
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             DO n = 1, 3
                IF (MINVAL(ABS(mesh%jjsi(1:2,ms)-jj(n,m)))==0) CYCLE
                face = n
             END DO
             SELECT CASE(type_fe)
             CASE(1)
                CALL element_2d_p1_boundary(face, dd, ww_s, n_w, l_Gs)
             CASE(2)
                CALL element_2d_p2_boundary(face, dd, ww_s, n_w, l_Gs)
             CASE(3)
                CALL element_2d_p3_boundary(face, dd, ww_s, n_w, l_Gs)
             END SELECT
             !===Fix for meshes not oriented correctly
             n1 = MODULO(face,3)+1; n2 = MODULO(face+1,3)+1
             IF (mesh%jjsi(1,ms)==jj(n1,m) .AND. mesh%jjsi(2,ms)==jj(n2,m)) THEN
                orient = 1
             ELSE
                orient = -1
             END IF
             IF (face==2) THEN
                orient = -orient
             END IF
             DO ls = 1, l_Gs
                ngauss(ls) = orient*ls+(l_Gs+1)*(1-orient)/2 !===rearrangement of gauss points
             END DO
             !===End fix
             DO ls = 1, l_Gs
                DO k = 1, k_d
                   rs = rr(k,mesh%jjsi(:,ms))
                   drs(1, k) = SUM(rs*dds(1,:,ls))
                ENDDO
                rjacs = SQRT(drs(1,1)**2 + drs(1,2)**2)
                mesh%gauss%rji(ls, ms) = rjacs*pps(ls)
                rnor(1) =  drs(1,2)/rjacs
                rnor(2) = - drs(1,1)/rjacs
                mesh%gauss%rnormsi(1, ls, cote, ms) = rnor(1)
                mesh%gauss%rnormsi(2, ls, cote, ms) = rnor(2)
                rsd = rr(:,jj(face,m)) - (rr(:,mesh%jjsi(1,ms))+rr(:,mesh%jjsi(2,ms)))/2
                x = SUM(rnor(:)*rsd)
                IF (x>0) THEN !===Verify the normal is outward
                   rnor = -rnor
                END IF
                DO k = 1, k_d
                   r = rr(k, jj(:,m))
                   DO k1 = 1, k_d
                      dr(k, k1) = SUM(r * dd(k1,:,ls))
                   ENDDO
                ENDDO
                rjac = (dr(1,1)*dr(2,2) - dr(1,2)*dr(2,1))
                DO n = 1, n_w !===Fix with ngauss
                   mesh%gauss%dwni(n, ngauss(ls), cote, ms)   &
                        = rnor(1)*(+ dd(1,n,ls)*dr(2,2) - dd(2,n,ls)*dr(2,1))/rjac &
                        + rnor(2)*(- dd(1,n,ls)*dr(1,2) + dd(2,n,ls)*dr(1,1))/rjac
                ENDDO
             ENDDO
          END DO
       ENDDO
       !===Some verifications
       DO ms = 1, mesh%mi
          DO cote = 1, 2
             m = mesh%neighi(cote,ms)
             DO ls = 1, l_Gs
                x = ABS(SQRT(SUM(mesh%gauss%rnormsi(:, ls, cote, ms)**2))-1.d0)
                IF (x.GE.1.d-13) THEN
                   write(*,*) 'Bug1 in normal '
                   STOP
                END IF
                rjac = (rr(1, mesh%jjsi(2,ms))-rr(1, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(1, ls, cote, ms) &
                     + (rr(2, mesh%jjsi(2,ms))-rr(2, mesh%jjsi(1,ms)))*mesh%gauss%rnormsi(2, ls, cote, ms)
                x = SQRT((rr(1, mesh%jjsi(2,ms))-rr(1, mesh%jjsi(1,ms)))**2 +&
                     (rr(2, mesh%jjsi(2,ms))-rr(2, mesh%jjsi(1,ms)))**2)
                IF (ABS(rjac/x) .GE. 1.d-10)  THEN
                   write(*,*) 'Bug2 in normal ', ms, cote, ls, ABS(rjac/x)
                   STOP
                END IF
             END DO
          END DO
       END DO
    END IF

    DEALLOCATE(dd,dds,pp,pps,r,rs,ww_s,dds_v,ngauss)
  END SUBROUTINE gauss_points_2d
END MODULE mod_gauss_points_2d
