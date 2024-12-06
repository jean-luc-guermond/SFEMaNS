!
!Authors Jean-Luc Guermond, Copyrights 2000
!
MODULE prep_interface

   IMPLICIT NONE

   PUBLIC :: load_interface
   PRIVATE

CONTAINS

   !------------------------------------------------------------------------------


   SUBROUTINE load_interface(mesh_master, mesh_slave, list_inter, mesh_INTERFACE, disjoint)

      USE def_type_mesh

      TYPE(mesh_type), INTENT(IN) :: mesh_master, mesh_slave
      INTEGER, DIMENSION(:), INTENT(IN) :: list_inter
      TYPE(interface_type), INTENT(OUT) :: mesh_INTERFACE
      LOGICAL, INTENT(IN) :: disjoint

      INTEGER :: dim, ms, ms1, ms2, ns, k, nws_master, nws_slave, m1, m2, k1, k2, cell_g, n, nw
      INTEGER, DIMENSION(:), ALLOCATABLE :: list, interface_mesh1, interface_mesh2
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: interface_jjs1, interface_jjs2
      REAL(KIND = 8) :: eps_ref = 1.d-7, r_norm, epsilon
      INTEGER, DIMENSION(2) :: n1_ks, n2_ks
      LOGICAL :: okay
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin_elem

      ! Removing gauss, FL, Mar. 23
      !!$    nws_master = mesh_master%gauss%n_ws
      !!$    nws_slave  = mesh_slave%gauss%n_ws

      nws_master = SIZE(mesh_master%jjs, 1)
      nws_slave = SIZE(mesh_slave%jjs, 1)
      ! Done removing gauss, FL, Mar. 23

      IF (nws_master > nws_slave) THEN
         WRITE(*, *) ' BUG in load_interface: nws_master > nws_slave '
         STOP
      END IF

      dim = SIZE(mesh_master%rr, 1)
      IF (dim>2) THEN
         WRITE(*, *) ' Dimension 3 not yet programmed '
         STOP
      END IF

      ALLOCATE(virgin_elem(mesh_slave%mes), list(dim), interface_mesh1(mesh_master%mes), &
           interface_mesh2(mesh_slave%mes), interface_jjs1(nws_master, mesh_master%mes), &
           interface_jjs2(nws_slave, mesh_slave%mes))

      virgin_elem = .TRUE.

      ms = 0
      DO ms1 = 1, mesh_master%mes
         IF(MINVAL(ABS(list_inter - mesh_master%sides(ms1))) /= 0) CYCLE !not on interface
         r_norm = SUM(ABS(mesh_master%rr(:, mesh_master%jjs(1, ms1)) - &
              &    mesh_master%rr(:, mesh_master%jjs(2, ms1))))
         epsilon = eps_ref * r_norm
         okay = .FALSE.

         lp2 : DO ms2 = 1, mesh_slave%mes
            !IF(.NOT.virgin_elem(ms2)) CYCLE !element done!
            IF(MINVAL(ABS(list_inter - mesh_slave%sides(ms2))) /= 0) CYCLE !not on interface

            DO k = 0, dim - 1 !dim=2
               DO ns = 1, dim !dim=2
                  list(ns) = MODULO(ns - 1 + k, dim) + 1
               END DO

               IF (MAXVAL(ABS(mesh_master%rr(:, mesh_master%jjs(list, ms1))&
                    - mesh_slave%rr(:, mesh_slave%jjs(1:dim, ms2)))).GT.epsilon) CYCLE

               IF(.NOT.virgin_elem(ms2)) THEN
                  Okay = .TRUE.
                  CYCLE !already element done
               ENDIF

               m1 = mesh_master%neighs(ms1)
               m2 = mesh_slave%neighs(ms2)
               r_norm = SUM(ABS(mesh_master%rr(:, mesh_master%jj(1:3, m1)) &
                    - mesh_slave%rr(:, mesh_slave%jj(1:3, m2))))
               IF (r_norm .LE. 1d-9) THEN
                  CYCLE  ! two identical triangles
               END IF

               ms = ms + 1
               interface_mesh1(ms) = ms1
               interface_mesh2(ms) = ms2
               interface_jjs1(1:2, ms) = mesh_master%jjs(list, ms1)
               interface_jjs2(:, ms) = mesh_slave%jjs(:, ms2)
               IF (nws_master==3) THEN !P2 in two dimensions
                  interface_jjs1(3, ms) = mesh_master%jjs(3, ms1)
               END IF

               virgin_elem(ms2) = .FALSE.
               IF (.NOT.disjoint) virgin_elem(ms1) = .FALSE.
               okay = .TRUE.

               EXIT lp2
            END DO
         END DO lp2
         IF (.NOT.okay) THEN
            WRITE(*, *) ' BUG in load_interface: .NOT.okay'
            STOP
         END IF
      END DO

      mesh_interface%mes = ms

      ALLOCATE(mesh_interface%mesh1(ms), mesh_interface%mesh2(ms), &
           mesh_interface%jjs1(nws_master, ms), mesh_interface%jjs2(nws_slave, ms))

      mesh_interface%mesh1 = interface_mesh1(1:ms)
      mesh_interface%mesh2 = interface_mesh2(1:ms)
      mesh_interface%jjs1 = interface_jjs1(1:nws_master, 1:ms)
      mesh_interface%jjs2 = interface_jjs2(1:nws_slave, 1:ms)

      DEALLOCATE(virgin_elem, list, interface_mesh1, interface_mesh2, &
           interface_jjs1, interface_jjs2)

      ALLOCATE(virgin_elem(mesh_slave%mes_extra), list(dim), interface_mesh1(mesh_master%mes_extra), &
           interface_mesh2(mesh_slave%mes_extra), interface_jjs1(nws_master, mesh_master%mes_extra), &
           interface_jjs2(nws_slave, mesh_slave%mes_extra))

      virgin_elem = .TRUE.
      nw = 3
      ms = 0
      DO ms1 = 1, mesh_master%mes_extra

         IF(MINVAL(ABS(list_inter - mesh_master%sides_extra(ms1))) /= 0) CYCLE !not on interface

         cell_g = mesh_master%neighs_extra(ms1)
         DO m1 = 1, mesh_master%mextra !find associated extra cell
            IF (mesh_master%jcc_extra(m1) == cell_g) EXIT
         END DO
         DO n = 1, dim + 1 !===find side in cell
            IF (MINVAL(ABS(mesh_master%jj_extra(n, m1) - mesh_master%jjs_extra(:, ms1)))/=0) THEN
               k1 = n
               EXIT
            END IF
         ENDDO
         !==cell index of edge
         n1_ks = (/MODULO(k1, nw) + 1, MODULO(k1 + 1, nw) + 1/)

         r_norm = SUM(ABS(mesh_master%rrs_extra(:, n1_ks(1), ms1) - mesh_master%rrs_extra(:, n1_ks(2), ms1)))
         epsilon = eps_ref * r_norm
         okay = .FALSE.

         lp3 : DO ms2 = 1, mesh_slave%mes_extra
            IF(MINVAL(ABS(list_inter - mesh_slave%sides_extra(ms2))) /= 0) CYCLE !not on interface

            cell_g = mesh_slave%neighs_extra(ms2)
            DO m2 = 1, mesh_slave%mextra !find associated extra cell
               IF (mesh_slave%jcc_extra(m2) == cell_g) EXIT
            END DO

            DO n = 1, dim + 1 !===find side in cell
               IF (MINVAL(ABS(mesh_slave%jj_extra(n, m2) - mesh_slave%jjs_extra(:, ms2)))/=0) THEN
                  k2 = n
                  EXIT
               END IF
            ENDDO

            DO k = 0, dim - 1 !dim=2
               !==cell index of edge
               IF (k == 0) THEN
                  n2_ks = (/MODULO(k2, nw) + 1, MODULO(k2 + 1, nw) + 1/)
               ELSE
                  n2_ks = (/MODULO(k2 + 1, nw) + 1, MODULO(k2, nw) + 1/)
               END IF

               IF (MAXVAL(ABS(mesh_master%rrs_extra(:, n1_ks, ms1) - mesh_slave%rrs_extra(:, n2_ks, ms2))).GT.epsilon) THEN
                  CYCLE
               END IF

               IF(.NOT.virgin_elem(ms2)) THEN
                  Okay = .TRUE.
                  CYCLE !element already done
               ENDIF

               r_norm = SUM(ABS(mesh_master%rrs_extra(:, 1:3, ms1) - mesh_slave%rrs_extra(:, 1:3, ms2)))
               IF (r_norm .LE. 1d-9) CYCLE ! two identical triangles

               ms = ms + 1
               interface_mesh1(ms) = ms1
               interface_mesh2(ms) = ms2
               interface_jjs1(1:2, ms) = mesh_master%jjs(1:2, ms1)
               interface_jjs2(:, ms) = mesh_slave%jjs(:, ms2)
               IF (nws_master==3) THEN !P2 in two dimensions
                  interface_jjs1(3, ms) = mesh_master%jjs(3, ms1)
               END IF

               virgin_elem(ms2) = .FALSE.
               IF (.NOT.disjoint) virgin_elem(ms1) = .FALSE.
               okay = .TRUE.

               EXIT lp3
            END DO
         END DO lp3
         IF (.NOT.okay) THEN
            WRITE(*, *) ' BUG in load_interface extras: .NOT.okay'
            !STOP
         END IF
      END DO

      mesh_interface%mes_extra = ms
      !write(*, *) 'm1', ms, mesh_master%mes_extra, mesh_master%mextra
      !write(*, *) 'jcc', mesh_master%neighs_extra
      !write(*, *) 'side', mesh_master%sides_extra
      ALLOCATE(mesh_interface%mesh1_extra(ms), mesh_interface%mesh2_extra(ms), &
           mesh_interface%jjs1_extra(nws_master, ms), mesh_interface%jjs2_extra(nws_slave, ms))
      IF (ms > 0) THEN
         mesh_interface%mesh1_extra = interface_mesh1(1:ms)
         mesh_interface%mesh2_extra = interface_mesh2(1:ms)
         mesh_interface%jjs1_extra = interface_jjs1(1:nws_master, 1:ms)
         mesh_interface%jjs2_extra = interface_jjs2(1:nws_slave, 1:ms)
      END IF

      DEALLOCATE(virgin_elem, list, interface_mesh1, interface_mesh2, &
           interface_jjs1, interface_jjs2)

   END SUBROUTINE load_interface

END MODULE prep_interface
