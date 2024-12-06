!
!Authors: Jean-Luc Guermond Copyrights 2005
!
MODULE fem_tn_NS_MHD
  IMPLICIT NONE

CONTAINS

  !DCQ, compute L1_norm using only list_mode(mode_idx) mode
  FUNCTION norme_L1_one_mode(mesh, mode_idx, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_axi
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v
    INTEGER                                     :: mode_idx
    REAL(KIND=8)                        :: err1, s1, s2,norm
    INTEGER                             :: nn

    err1 = 0.d0
    s2 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO nn = 1,SIZE(v,1)
          CALL ns_l1(mesh , abs(v(nn,:,mode_idx)), s1)
          s2 = s2 + s1
       END DO
       err1 = err1 + s2
       ! CN-AR Tue Jan 13 2009
    ELSE
       DO nn = 1,SIZE(v,2)
          CALL ns_l1(mesh , abs(v(:,nn,mode_idx)), s1)
          s2 = s2 + s1
       END DO
       ! CN-AR Tue Jan 13 2009
       ! JLG/CN correction du bug CN-AR April 7, 2010
       err1 = err1 + s2
       ! CN-AR Tue Jan 13 2009
    END IF
    norm = err1
  END FUNCTION norme_L1_one_mode

  FUNCTION norme_L2_champ(mesh, list_mode, v) RESULT(norm)
    USE def_type_mesh
    USE fem_tn_axi
    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v
    REAL(KIND=8)                        :: err1, s1, s2, norm
    INTEGER                             :: k, nn

    err1 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO k = 1, SIZE(list_mode)
          s2 = 0.d0
          DO nn = 1,SIZE(v,1)
             s1 = 0.d0
             CALL ns_0(mesh , v(nn,:,k), s1)
             s2 = s2 + s1**2
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2 = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009

       ENDDO
    ELSE
       DO k = 1, SIZE(list_mode)
          s2 = 0.d0
          DO nn = 1,SIZE(v,2)
             s1 = 0.d0
             CALL ns_0(mesh , v(:,nn,k), s1)
             s2 = s2 + s1**2
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2 = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009

       ENDDO
    END IF
    norm = SQRT(err1)

  END FUNCTION norme_L2_champ

  FUNCTION dot_product_champ(mesh, list_mode, v, w) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_axi

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v, w

    REAL(KIND=8)                        :: err1, s1, s2, norm
    INTEGER                             :: k, nn

    err1 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO k = 1, SIZE(list_mode)
          s2 = 0.d0
          DO nn = 1,SIZE(v,1)
             s1 = 0.d0
             CALL dot_product(mesh , v(nn,:,k),  w(nn,:,k), s1)
             s2 = s2 + s1
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2 = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009

       ENDDO
    ELSE
       DO k = 1, SIZE(list_mode)
          s2 = 0.d0
          DO nn = 1,SIZE(v,2)
             s1 = 0.d0
             CALL dot_product(mesh , v(:,nn,k), w(:,nn,k), s1)
             s2 = s2 + s1
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2 = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009

       ENDDO
    END IF
    norm = err1

  END FUNCTION dot_product_champ

  FUNCTION norme_H1_champ(mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_axi

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) ::  v

    REAL(KIND=8)                        :: err1, s1, s2, s0, norm
    INTEGER                             :: k,nn

    err1 = 0.d0
    IF (SIZE(v,2)==mesh%np) THEN
       DO k = 1, SIZE(list_mode)
          s2 = 0.d0
          DO nn = 1,SIZE(v,1)
             CALL ns_1(mesh , v(nn,:,k), s1)
             CALL ns_0(mesh , v(nn,:,k), s0)
             s2 = s2+ s1**2 + list_mode(k)**2*s0**2
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2  = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009
       ENDDO
    ELSE
       DO k = 1, SIZE(list_mode)
          s2=0.d0
          DO nn = 1,SIZE(v,2)
             CALL ns_1(mesh , v(:,nn,k), s1)
             CALL ns_0(mesh , v(:,nn,k), s0)
             s2 = s2 + s1**2 + list_mode(k)**2*s0**2
          END DO
          ! CN-AR Tue Jan 13 2009
          ! JLG/CN correction du bug CN-AR April 7, 2010
          IF (list_mode(k) /= 0) THEN
             s2  = s2/2.d0
          END IF
          err1 = err1 + s2
          ! CN-AR Tue Jan 13 2009
       ENDDO
    END IF

    norm = SQRT(err1)

  END FUNCTION norme_H1_champ

  FUNCTION norme_div(H_mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_axi

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN), TARGET :: v


    INTEGER              :: m_max_c, mode, k, m, l, ni, i
    REAL(KIND=8)         :: norm, err, div1, div2, jr, ray

    err = 0.d0

    m_max_c = SIZE(list_mode)

    IF (SIZE(v,2)==H_mesh%np) THEN

       DO m = 1, H_mesh%me
          DO l = 1, H_mesh%gauss%l_G

             !===Compute radius of Gauss point
             ray = 0.d0
             DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
                ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             END DO
             jr = ray * H_mesh%gauss%rj(l,m)


             DO k=1, m_max_c
                mode = list_mode(k)
                div1 = 0.d0
                div2 = 0.d0

                DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)

                   div1 = div1 + v(1,i,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + v(1,i,k)*H_mesh%gauss%dw(1,ni,l,m)  &
                        + mode/ray*v(4,i,k)*H_mesh%gauss%ww(ni,l) &
                        + v(5,i,k)*H_mesh%gauss%dw(2,ni,l,m)

                   div2 = div2 + v(2,i,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + v(2,i,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        - mode/ray*v(3,i,k)*H_mesh%gauss%ww(ni,l) &
                        + v(6,i,k)*H_mesh%gauss%dw(2,ni,l,m)

                ENDDO

                err = err + (div1**2+div2**2)*jr
             ENDDO

          END DO
       END DO
    ELSE
       DO m = 1, H_mesh%me
          DO l = 1, H_mesh%gauss%l_G

             !===Compute radius of Gauss point
             ray = 0.d0
             DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
                ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             END DO
             jr = ray * H_mesh%gauss%rj(l,m)


             DO k=1, m_max_c
                mode = list_mode(k)
                div1 = 0.d0
                div2 = 0.d0

                DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)

                   div1 = div1   + v(i,1,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + v(i,1,k)*H_mesh%gauss%dw(1,ni,l,m)  &
                        + mode/ray*v(i,4,k)*H_mesh%gauss%ww(ni,l) &
                        + v(i,5,k)*H_mesh%gauss%dw(2,ni,l,m)

                   div2 = div2   + v(i,2,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + v(i,2,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        - mode/ray*v(i,3,k)*H_mesh%gauss%ww(ni,l) &
                        + v(i,6,k)*H_mesh%gauss%dw(2,ni,l,m)

                ENDDO

                err = err + (div1**2+div2**2)*jr
             ENDDO

          END DO
       END DO

    END IF


    norm = SQRT(err)

  END FUNCTION norme_div

  FUNCTION norme_curl(H_mesh, list_mode, v) RESULT(norm)

    USE def_type_mesh
    USE fem_tn_axi

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: H_mesh !type de maillage
    INTEGER, DIMENSION(:),           INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:),  INTENT(IN) :: v
    REAL(KIND=8)         :: err, norm, jr, ray
    INTEGER              :: m_max_c, mode, k, m, l, ni, i
    REAL(KIND=8), DIMENSION(6) :: c

    m_max_c = SIZE(list_mode)
    err = 0.d0

    IF (SIZE(v,2)==H_mesh%np) THEN

       DO m = 1, H_mesh%me
          DO l = 1, H_mesh%gauss%l_G

             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
                ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             END DO
             jr = ray * H_mesh%gauss%rj(l,m)

             DO k=1, m_max_c
                mode = list_mode(k)
                c = 0
                DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)
                   !--------Composante r------
                   c(1) = c(1) + ( mode/ray*v(6,i,k)*H_mesh%gauss%ww(ni,l) &
                        - v(3,i,k)*H_mesh%gauss%dw(2,ni,l,m))
                   c(2) = c(2) + (-mode/ray*v(5,i,k)*H_mesh%gauss%ww(ni,l) &
                        - v(4,i,k)*H_mesh%gauss%dw(2,ni,l,m))
                   !--------Composante theta------
                   c(3) = c(3) + (v(1,i,k)*H_mesh%gauss%dw(2,ni,l,m) &
                        - v(5,i,k)*H_mesh%gauss%dw(1,ni,l,m))
                   c(4) = c(4) + (v(2,i,k)*H_mesh%gauss%dw(2,ni,l,m) &
                        - v(6,i,k)*H_mesh%gauss%dw(1,ni,l,m))
                   !--------Composante z------
                   c(5) = c(5) + (v(3,i,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        + v(3,i,k)*H_mesh%gauss%ww(ni,l)/ray &
                        - mode/ray*v(2,i,k)*H_mesh%gauss%ww(ni,l))
                   c(6) = c(6) + (v(4,i,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        + v(4,i,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + mode/ray*v(1,i,k)*H_mesh%gauss%ww(ni,l))
                ENDDO
                err = err + SUM(c**2)*jr
             END DO
          END DO
       END DO

    ELSE
       DO m = 1, H_mesh%me
          DO l = 1, H_mesh%gauss%l_G

             !===Compute radius of Gauss point
             ray = 0
             DO ni = 1, H_mesh%gauss%n_w;  i = H_mesh%jj(ni,m)
                ray = ray + H_mesh%rr(1,i)*H_mesh%gauss%ww(ni,l)
             END DO
             jr = ray * H_mesh%gauss%rj(l,m)

             DO k=1, m_max_c
                mode = list_mode(k)
                c = 0
                DO ni = 1,H_mesh%gauss%n_w; i = H_mesh%jj(ni,m)
                   !--------Composante r------
                   c(1) = c(1) + ( mode/ray*v(i,6,k)*H_mesh%gauss%ww(ni,l) &
                        - v(i,3,k)*H_mesh%gauss%dw(2,ni,l,m))
                   c(2) = c(2) + (-mode/ray*v(i,5,k)*H_mesh%gauss%ww(ni,l) &
                        - v(i,4,k)*H_mesh%gauss%dw(2,ni,l,m))
                   !--------Composante theta------
                   c(3) = c(3) + (v(i,1,k)*H_mesh%gauss%dw(2,ni,l,m) &
                        - v(i,5,k)*H_mesh%gauss%dw(1,ni,l,m))
                   c(4) = c(4) + (v(i,2,k)*H_mesh%gauss%dw(2,ni,l,m) &
                        - v(i,6,k)*H_mesh%gauss%dw(1,ni,l,m))
                   !--------Composante z------
                   c(5) = c(5) + (v(i,3,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        + v(i,3,k)*H_mesh%gauss%ww(ni,l)/ray &
                        - mode/ray*v(i,2,k)*H_mesh%gauss%ww(ni,l))
                   c(6) = c(6) + (v(i,4,k)*H_mesh%gauss%dw(1,ni,l,m) &
                        + v(i,4,k)*H_mesh%gauss%ww(ni,l)/ray &
                        + mode/ray*v(i,1,k)*H_mesh%gauss%ww(ni,l))
                ENDDO
                err = err + SUM(c**2)*jr
             END DO
          END DO
       END DO
    END IF

    norm = SQRT(err)

  END FUNCTION norme_curl

  SUBROUTINE  norme_interface(H_mesh,phi_mesh,INTERFACE,mu_H_field,mu_phi,H,phi,mode,x)
    USE def_type_mesh
    USE Dir_nodes
    USE gauss_points

    IMPLICIT NONE

    TYPE(mesh_type),              INTENT(IN)  :: H_mesh, phi_mesh
    TYPE(interface_type),         INTENT(IN)  :: INTERFACE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: H
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: phi
    INTEGER,                      INTENT(IN)  :: mode
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: mu_H_field
    REAL(KIND=8),                 INTENT(IN)  :: mu_phi
    REAL(KIND=8),                 INTENT(OUT) :: x
    REAL(KIND=8), DIMENSION(2) :: rgauss
    REAL(KIND=8), DIMENSION(phi_mesh%gauss%n_ws,phi_mesh%gauss%l_Gs)   :: w_cs


    INTEGER :: ms, ls, ms2, n_ws1, n_ws2, m, i, ni
    REAL(KIND=8), DIMENSION(6) :: b, mub, grd
    REAL(KIND=8) :: z, zmu, ray, err, muhl
    CALL gauss(phi_mesh)

    IF (H_mesh%gauss%n_ws == n_ws) THEN
       w_cs = wws
    ELSE
       DO ls = 1, l_Gs
          w_cs(1,ls)= wws(1,ls)+0.5*wws(3,ls)
          w_cs(2,ls)= wws(2,ls)+0.5*wws(3,ls)
          w_cs(3,ls)=0
       END DO
    END IF

    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = phi_mesh%gauss%n_ws

    err = 0
    x = 0
    z = 0
    zmu = 0

    IF (SIZE(H,2)==H_mesh%np) THEN
       DO ms = 1, interface%mes
          ms2 = interface%mesh2(ms)
          m = phi_mesh%neighs(ms2)

          DO ls = 1, l_Gs

             !===Compute radius of Gauss point
             ray = 0.d0
             DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,ms2)
                ray = ray + phi_mesh%rr(1,i)* wws(ni,ls)
             END DO

             rgauss(1) = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* wws(:,ls))
             rgauss(2) = SUM(phi_mesh%rr(2,phi_mesh%jjs(:,ms2))* wws(:,ls))

             DO i = 1, 6
                b(i) = SUM(H(i,interface%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
                mub(i) = SUM(mu_H_field(interface%jjs1(1:n_ws1,ms))*H(i,interface%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
             ENDDO

             grd(1) = SUM(phi(1,phi_mesh%jj(:,m))*dw_s(1,:,ls,ms2))
             grd(2) = SUM(phi(2,phi_mesh%jj(:,m))*dw_s(1,:,ls,ms2))
             grd(3) = mode/ray * SUM(phi(2,interface%jjs2(:,ms))*wws(:,ls))
             grd(4) = -mode/ray * SUM(phi(1,interface%jjs2(:,ms))*wws(:,ls))
             grd(5) = SUM(phi(1,phi_mesh%jj(:,m))*dw_s(2,:,ls,ms2))
             grd(6) = SUM(phi(2,phi_mesh%jj(:,m))*dw_s(2,:,ls,ms2))

             z = z + SUM(b(:)**2)*rjs(ls,ms2)*ray
             zmu = zmu + SUM(mub(:)**2)*rjs(ls,ms2)*ray

             !Error on tangential component (magnetic induction)
             x  = x + ray* rjs(ls,ms2)*( &
                  ((b(5)-grd(5))*rnorms(1,ls,ms2) - &
                  (b(1)-grd(1))*rnorms(2,ls,ms2))**2 + &
                  ((b(6)-grd(6))*rnorms(1,ls,ms2) - &
                  (b(2)-grd(2))*rnorms(2,ls,ms2))**2  + &
                  (b(3)-grd(3))**2 + (b(4)-grd(4))**2)

             !Error on normal component (magnetic field)
             err = err + ray* rjs(ls,ms2)*( &
                  ((mub(1)-mu_phi*grd(1))*rnorms(1,ls,ms2) +&
                  (mub(5)-mu_phi*grd(5))*rnorms(2,ls,ms2))**2 + &
                  ((mub(2)-mu_phi*grd(2))*rnorms(1,ls,ms2) +&
                  (mub(6)-mu_phi*grd(6))*rnorms(2,ls,ms2))**2)

          END DO
       END DO
    ELSE
       DO ms = 1, interface%mes
          ms2 = interface%mesh2(ms)
          m = phi_mesh%neighs(ms2)

          DO ls = 1, l_Gs
             !June 6 2008, muhl
             muhl = SUM(mu_H_field(interface%jjs1(1:n_ws1,ms))*w_cs(1:n_ws1,ls))
             !June 6 2008, muhl
             !===Compute radius of Gauss point
             ray = 0.d0
             DO ni = 1, n_ws2;  i = phi_mesh%jjs(ni,ms2)
                ray = ray + phi_mesh%rr(1,i)* wws(ni,ls)
             END DO

             rgauss(1) = SUM(phi_mesh%rr(1,phi_mesh%jjs(:,ms2))* wws(:,ls))
             rgauss(2) = SUM(phi_mesh%rr(2,phi_mesh%jjs(:,ms2))* wws(:,ls))

             DO i = 1, 6
                b(i) = SUM(H(interface%jjs1(1:n_ws1,ms),i)*w_cs(1:n_ws1,ls))
                mub(i) = b(i)*muhl
             ENDDO

             grd(1) = SUM(phi(phi_mesh%jj(:,m),1)*dw_s(1,:,ls,ms2))
             grd(2) = SUM(phi(phi_mesh%jj(:,m),2)*dw_s(1,:,ls,ms2))
             grd(3) = mode/ray * SUM(phi(interface%jjs2(:,ms),2)*wws(:,ls))
             grd(4) =-mode/ray * SUM(phi(interface%jjs2(:,ms),1)*wws(:,ls))
             grd(5) = SUM(phi(phi_mesh%jj(:,m),1)*dw_s(2,:,ls,ms2))
             grd(6) = SUM(phi(phi_mesh%jj(:,m),2)*dw_s(2,:,ls,ms2))

             z = z + SUM(b(:)**2)*rjs(ls,ms2)*ray
             zmu = zmu + SUM(mub(:)**2)*rjs(ls,ms2)*ray

             !Error on tangential component (magnetic induction)
             x  = x + ray* rjs(ls,ms2)*( &
                  ((b(5)-grd(5))*rnorms(1,ls,ms2) - &
                  (b(1)-grd(1))*rnorms(2,ls,ms2))**2 + &
                  ((b(6)-grd(6))*rnorms(1,ls,ms2) - &
                  (b(2)-grd(2))*rnorms(2,ls,ms2))**2   + &
                  (b(3)-grd(3))**2 + (b(4)-grd(4))**2)

             !Error on normal component (magnetic field)
             err = err + ray* rjs(ls,ms2)*( &
                  ((mub(1)-mu_phi*grd(1))*rnorms(1,ls,ms2) +&
                  (mub(5)-mu_phi*grd(5))*rnorms(2,ls,ms2))**2 + &
                  ((mub(2)-mu_phi*grd(2))*rnorms(1,ls,ms2) +&
                  (mub(6)-mu_phi*grd(6))*rnorms(2,ls,ms2))**2)

          END DO
       END DO
    END IF
    WRITE(*,'(A,e12.5)') 'Collage tangent Sigma_phi ',  SQRT(x)/(SQRT(z)+1.d-16)
    WRITE(*,'(A,e12.5)') 'Collage normal Sigma_phi  ',  SQRT(err)/(SQRT(zmu)+1.d-16)
    x = SQRT(x)/(SQRT(z)+1.d-16) + SQRT(err)/(SQRT(zmu)+1.d-16)

  END SUBROUTINE  norme_interface

  SUBROUTINE  norme_interface_H_mu(H_mesh,INTERFACE,mu_H_field,H,x)
    USE def_type_mesh
    USE Dir_nodes
    USE gauss_points

    IMPLICIT NONE

    TYPE(mesh_type),              INTENT(IN)  :: H_mesh
    TYPE(interface_type),         INTENT(IN)  :: INTERFACE
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)  :: H
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: mu_H_field
    REAL(KIND=8),                 INTENT(OUT) :: x
    REAL(KIND=8), DIMENSION(2) :: rgauss
    REAL(KIND=8), DIMENSION(H_mesh%gauss%n_ws,H_mesh%gauss%l_Gs)   :: w_1s, w_2s

    INTEGER :: ms, ls, ms2, n_ws1, n_ws2, m, i, ni
    REAL(KIND=8), DIMENSION(6) :: b1, mub1, b2, mub2
    REAL(KIND=8) :: z, zmu, ray, err
    CALL gauss(H_mesh)

    w_1s = wws
    w_2s = wws

    n_ws1 = H_mesh%gauss%n_ws
    n_ws2 = H_mesh%gauss%n_ws

    err = 0
    x = 0
    z = 0
    zmu = 0

    DO ms = 1, interface%mes
       ms2 = interface%mesh2(ms)
       m = H_mesh%neighs(ms2)

       DO ls = 1, l_Gs

          !===Compute radius of Gauss point
          ray = 0.d0
          DO ni = 1, n_ws2;  i = H_mesh%jjs(ni,ms2)
             ray = ray + H_mesh%rr(1,i)* wws(ni,ls)
          END DO

          rgauss(1) = SUM(H_mesh%rr(1,H_mesh%jjs(:,ms2))* wws(:,ls))
          rgauss(2) = SUM(H_mesh%rr(2,H_mesh%jjs(:,ms2))* wws(:,ls))

          DO i = 1, 6
             b1(i) = SUM(H(interface%jjs1(1:n_ws1,ms),i)*w_1s(1:n_ws1,ls))
             mub1(i) = SUM(mu_H_field(interface%jjs1(1:n_ws1,ms))*H(interface%jjs1(1:n_ws1,ms),i)*w_1s(1:n_ws1,ls))
             b2(i) = SUM(H(interface%jjs2(1:n_ws2,ms),i)*w_2s(1:n_ws2,ls))
             mub2(i) = SUM(mu_H_field(interface%jjs2(1:n_ws2,ms))*H(interface%jjs2(1:n_ws2,ms),i)*w_2s(1:n_ws2,ls))
          ENDDO

          z = z +  (SUM(b2(:)**2)+SUM(b1(:)**2))*rjs(ls,ms2)*ray
          zmu = zmu + (SUM(mub1(:)**2)+SUM(mub2(:)**2))*rjs(ls,ms2)*ray

          !Error on tangential component (magnetic field)
          x  = x + ray* rjs(ls,ms2)*( &
               ((b1(5)-b2(5))*rnorms(1,ls,ms2) - &
               (b1(1)-b2(1))*rnorms(2,ls,ms2))**2 + &
               ((b1(6)-b2(6))*rnorms(1,ls,ms2) - &
               (b1(2)-b2(2))*rnorms(2,ls,ms2))**2  + &
               (b1(3)-b2(3))**2 + (b1(4)-b2(4))**2)
          !Error on normal component (magnetic induction)
          err = err + ray* rjs(ls,ms2)*( &
               ((mub1(1)-mub2(1))*rnorms(1,ls,ms2) +&
               (mub1(5)-mub2(5))*rnorms(2,ls,ms2))**2 + &
               ((mub1(2)-mub2(2))*rnorms(1,ls,ms2) +&
               (mub1(6)-mub2(6))*rnorms(2,ls,ms2))**2)

       END DO
    END DO
    WRITE(*,'(A,e12.5)') 'Collage tangent Sigma_mu ',  SQRT(x)/(SQRT(z)+1.d-16)
    WRITE(*,'(A,e12.5)') 'Collage normal Sigma_mu  ',  SQRT(err)/(SQRT(zmu)+1.d-16)
    x = SQRT(x)/(SQRT(z)+1.d-16) + SQRT(err)/(SQRT(zmu)+1.d-16)

  END SUBROUTINE  norme_interface_H_mu

END MODULE fem_tn_NS_MHD
