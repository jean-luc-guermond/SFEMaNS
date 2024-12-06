MODULE refine_mesh
   USE mesh_tools
   PUBLIC :: create_iso_grid_distributed, refinement_iso_grid_distributed
   PRIVATE
CONTAINS

   SUBROUTINE create_iso_grid_distributed(mesh_p1, mesh, type_fe)
      !===jj(:, :)    nodes of the  volume_elements of the input grid
      !===jjs(:, :)    nodes of the surface_elements of the input grid
      !===rr(:, :)    cartesian coordinates of the nodes of the input grid
      !===m_op(:,:)   volume element opposite to each node
      !===neigh_el(:) volume element ajacent to the surface element
      !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
      !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
      !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
      USE def_type_mesh
      USE input_data
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p1, mesh
      INTEGER, INTENT(IN) :: type_fe
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: j_mid, jjs_mid
      INTEGER :: np, me, mes, nw, nws, kd, n, m, k, l, n_dof, dom_np
      INTEGER :: n1, n2, n3, n4, ms, n_start, n_end, n1_g, n2_g
      INTEGER :: n_k1, n_k2, m_op_k, kk, i, mm, ms_bord, p_e, p_c
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: r_mid
      INTEGER, DIMENSION(type_fe + 1) :: ns3
      REAL(KIND = 8), DIMENSION(2) :: rz
      REAL(KIND = 8), DIMENSION(type_fe + 1) :: scos
      REAL(KIND = 8) :: epsilon = 1.d-13, dist, d1, d2, s1, s2, s3, shalf, ref, scc, infinity, rescale
      INTEGER :: ns, ns1, index, nb_angle, f_dof, edge_g, edge_l, n_new_start, proc, nb_proc, edges, p, cell_g, cell_l
      INTEGER :: interface
      LOGICAL :: iso

      nw = SIZE(mesh_p1%jj, 1)   !===nodes in each volume element (3 in 2D)
      me = SIZE(mesh_p1%jj, 2)   !===number of cells
      kd = SIZE(mesh_p1%rr, 1)   !===space dimensions
      np = mesh_p1%np            !===number of P1 vertices connected to grid
      dom_np = mesh_p1%dom_np    !===number of P1 vertices attributed to proc
      mes = SIZE(mesh_p1%jjs, 2)
      nws = SIZE(mesh_p1%jjs, 1)
      f_dof = type_fe - 1
      nb_proc = SIZE(mesh_p1%domnp)

      IF (type_fe==1 .OR. mesh_p1%me == 0) THEN
         mesh%me = mesh_p1%me
         mesh%mes = mesh_p1%mes
         mesh%np = mesh_p1%np
         mesh%medge = mesh_p1%medge
         mesh%medges = mesh_p1%medges
         mesh%mextra = mesh_p1%mextra
         mesh%mes_extra = mesh_p1%mes_extra

         ALLOCATE(mesh%jj(nw, me))
         mesh%jj = mesh_p1%jj
         ALLOCATE(mesh%jjs(nws, mes))
         mesh%jjs = mesh_p1%jjs
         ALLOCATE(mesh%jjs_extra(nws, mesh%mes_extra))
         mesh%jjs_extra = mesh_p1%jjs_extra
         !ALLOCATE(mesh%iis(nws,mes))
         !mesh%iis = mesh_p1%iis
         ALLOCATE(mesh%jj_extra(nw, mesh%mextra))
         mesh%jj_extra = mesh_p1%jj_extra
         ALLOCATE(mesh%jce(nw, me))
         mesh%jce = mesh_p1%jce
         !ALLOCATE(mesh%jev(nw - 1, mesh%medge))
         !mesh%jev = mesh_p1%jev
         ALLOCATE(mesh%rr(kd, mesh%np))
         mesh%rr = mesh_p1%rr
         ALLOCATE(mesh%neigh(nw, mesh%me))
         mesh%neigh = mesh_p1%neigh
         ALLOCATE(mesh%sides(mesh%mes))
         mesh%sides = mesh_p1%sides
         ALLOCATE(mesh%neighs(mesh%mes))
         mesh%neighs = mesh_p1%neighs
         ALLOCATE(mesh%sides_extra(mesh%mes_extra))
         mesh%sides_extra = mesh_p1%sides_extra
         ALLOCATE(mesh%neighs_extra(mesh%mes_extra))
         mesh%neighs_extra = mesh_p1%neighs_extra
         ALLOCATE(mesh%rrs_extra(kd, nw, mesh%mes_extra))
         mesh%rrs_extra = mesh_p1%rrs_extra
         ALLOCATE(mesh%i_d(mesh%me))
         mesh%i_d = mesh_p1%i_d
         ALLOCATE(mesh%loc_to_glob(mesh%np))
         mesh%loc_to_glob = mesh_p1%loc_to_glob
         ALLOCATE(mesh%jcc_extra(mesh%mextra))
         mesh%jcc_extra = mesh_p1%jcc_extra

         mesh%dom_me = mesh_p1%dom_me
         mesh%dom_np = mesh_p1%dom_np
         mesh%dom_mes = mesh_p1%dom_mes
         ALLOCATE(mesh%disp(nb_proc + 1), mesh%domnp(nb_proc))
         mesh%domnp = mesh_p1%domnp
         mesh%disp = mesh_p1%disp
         ALLOCATE(mesh%discell(nb_proc + 1), mesh%domcell(nb_proc))
         mesh%domcell = mesh_p1%domcell
         mesh%discell = mesh_p1%discell
         ALLOCATE(mesh%disedge(nb_proc + 1), mesh%domedge(nb_proc))
         mesh%domedge = mesh_p1%domedge
         mesh%disedge = mesh_p1%disedge

         mesh%nis = mesh_p1%nis
         ALLOCATE(mesh%isolated_interfaces(mesh_p1%nis, 2))
         mesh%isolated_interfaces = mesh_p1%isolated_interfaces
         ALLOCATE(mesh%isolated_jjs(mesh_p1%nis))
         mesh%isolated_jjs = mesh_p1%isolated_jjs

         ALLOCATE(mesh%jce_extra(SIZE(mesh_p1%jce_extra, 1), SIZE(mesh_p1%jce_extra, 2)))
         mesh%jce_extra = mesh_p1%jce_extra
         ALLOCATE(mesh%jees(SIZE(mesh_p1%jees)))
         mesh%jees = mesh_p1%jees
         ALLOCATE(mesh%jecs(SIZE(mesh_p1%jecs)))
         mesh%jecs = mesh_p1%jecs

         mesh%gauss%n_w = 3
         mesh%gauss%n_ws = 2
         RETURN
      ELSE IF (type_fe==2) THEN
         mesh%me = mesh_p1%me
         mesh%mes = mesh_p1%mes
         mesh%np = mesh_p1%np + mesh_p1%medge + mesh_p1%medges
         mesh%medge = mesh_p1%medge
         mesh%medges = mesh_p1%medges
         mesh%mextra = mesh_p1%mextra
         mesh%mes_extra = mesh_p1%mes_extra

         ALLOCATE(mesh%jj(nw * (f_dof + 1), me))   !---->
         ALLOCATE(mesh%jjs(nws + f_dof, mes))   !---->
         ALLOCATE(mesh%jjs_extra(nws + f_dof, mesh%mes_extra))
         ALLOCATE(mesh%jj_extra(nw * (f_dof + 1), mesh%mextra)) !---->
         ALLOCATE(mesh%rrs_extra(kd, nw * (f_dof + 1), mesh%mes_extra))
         ALLOCATE(mesh%rr(kd, mesh%np))    !---->
         ALLOCATE(mesh%loc_to_glob(mesh%np)) !---->

         ALLOCATE(mesh%jce(nw, me))
         mesh%jce = mesh_p1%jce
         !ALLOCATE(mesh%jev(nw - 1, mesh%medge))
         !mesh%jev = mesh_p1%jev
         ALLOCATE(mesh%neigh(nw, mesh%me))
         mesh%neigh = mesh_p1%neigh
         ALLOCATE(mesh%sides(mesh%mes))
         mesh%sides = mesh_p1%sides
         ALLOCATE(mesh%neighs(mesh%mes))
         mesh%neighs = mesh_p1%neighs
         ALLOCATE(mesh%sides_extra(mesh%mes_extra))
         mesh%sides_extra = mesh_p1%sides_extra
         ALLOCATE(mesh%neighs_extra(mesh%mes_extra))
         mesh%neighs_extra = mesh_p1%neighs_extra
         ALLOCATE(mesh%i_d(mesh%me))
         mesh%i_d = mesh_p1%i_d
         ALLOCATE(mesh%jcc_extra(mesh%mextra))
         mesh%jcc_extra = mesh_p1%jcc_extra

         mesh%dom_me = mesh_p1%dom_me
         mesh%dom_np = mesh_p1%dom_np + mesh_p1%medge
         mesh%dom_mes = mesh_p1%dom_mes
         ALLOCATE(mesh%disp(nb_proc + 1), mesh%domnp(nb_proc))
         mesh%domnp = mesh_p1%domnp + mesh_p1%domedge
         mesh%disp = mesh_p1%disp + mesh_p1%disedge - 1
         ALLOCATE(mesh%discell(nb_proc + 1), mesh%domcell(nb_proc))
         mesh%domcell = mesh_p1%domcell
         mesh%discell = mesh_p1%discell
         ALLOCATE(mesh%disedge(nb_proc + 1), mesh%domedge(nb_proc))
         mesh%domedge = mesh_p1%domedge
         mesh%disedge = mesh_p1%disedge

         mesh%nis = mesh_p1%nis
         ALLOCATE(mesh%isolated_interfaces(mesh_p1%nis, 2))
         mesh%isolated_interfaces = mesh_p1%isolated_interfaces
         ALLOCATE(mesh%isolated_jjs(mesh_p1%nis))

         ALLOCATE(mesh%jce_extra(SIZE(mesh_p1%jce_extra, 1), SIZE(mesh_p1%jce_extra, 2)))
         mesh%jce_extra = mesh_p1%jce_extra
         ALLOCATE(mesh%jees(SIZE(mesh_p1%jees)))
         mesh%jees = mesh_p1%jees
         ALLOCATE(mesh%jecs(SIZE(mesh_p1%jecs)))
         mesh%jecs = mesh_p1%jecs

         mesh%gauss%n_w = 6
         mesh%gauss%n_ws = 3
      ELSE IF (type_fe==3) THEN
         mesh%me = mesh_p1%me
         mesh%mes = mesh_p1%mes
         mesh%np = mesh_p1%np + 2 * mesh_p1%medge + 2 * mesh_p1%medges + mesh_p1%me
         mesh%medge = mesh_p1%medge
         mesh%mextra = mesh_p1%mextra
         mesh%mes_extra = mesh_p1%mes_extra

         ALLOCATE(mesh%jj(nw * (f_dof + 1) + 1, me))   !----> done
         ALLOCATE(mesh%jjs(nws + f_dof, mes))   !---->
         ALLOCATE(mesh%jjs_extra(nws + f_dof, mesh%mextra))   !---->
         ALLOCATE(mesh%jj_extra(nw * (f_dof + 1) + 1, mesh%mextra)) !---->
         ALLOCATE(mesh%rr(kd, mesh%np))    !----> done
         ALLOCATE(mesh%rrs_extra(kd, nw * (f_dof + 1) + 1, mesh%mes_extra))
         ALLOCATE(mesh%loc_to_glob(mesh%np)) !----> done

         ALLOCATE(mesh%jce(nw, me))
         mesh%jce = mesh_p1%jce
         !ALLOCATE(mesh%jev(nw - 1, mesh%medge))
         !mesh%jev = mesh_p1%jev
         ALLOCATE(mesh%neigh(nw, mesh%me))
         mesh%neigh = mesh_p1%neigh
         ALLOCATE(mesh%sides(mesh%mes))
         mesh%sides = mesh_p1%sides
         ALLOCATE(mesh%neighs(mesh%mes))
         mesh%neighs = mesh_p1%neighs
         ALLOCATE(mesh%sides_extra(mesh%mes_extra))
         mesh%sides_extra = mesh_p1%sides_extra
         ALLOCATE(mesh%neighs_extra(mesh%mes_extra))
         mesh%neighs_extra = mesh_p1%neighs_extra
         ALLOCATE(mesh%i_d(mesh%me))
         mesh%i_d = mesh_p1%i_d
         ALLOCATE(mesh%jcc_extra(mesh%mextra))
         mesh%jcc_extra = mesh_p1%jcc_extra

         mesh%dom_me = mesh_p1%dom_me
         mesh%dom_np = mesh_p1%dom_np + 2 * mesh_p1%medge + mesh_p1%me
         mesh%dom_mes = mesh_p1%dom_mes
         ALLOCATE(mesh%disp(nb_proc + 1), mesh%domnp(nb_proc))
         mesh%domnp = mesh_p1%domnp + mesh_p1%domedge * 2 + mesh_p1%domcell
         mesh%disp = mesh_p1%disp + mesh_p1%disedge * 2 + mesh_p1%discell - 3
         ALLOCATE(mesh%discell(nb_proc + 1), mesh%domcell(nb_proc))
         mesh%domcell = mesh_p1%domcell
         mesh%discell = mesh_p1%discell
         ALLOCATE(mesh%disedge(nb_proc + 1), mesh%domedge(nb_proc))
         mesh%domedge = mesh_p1%domedge
         mesh%disedge = mesh_p1%disedge

         mesh%nis = mesh_p1%nis
         ALLOCATE(mesh%isolated_interfaces(mesh_p1%nis, 2))
         mesh%isolated_interfaces = mesh_p1%isolated_interfaces
         ALLOCATE(mesh%isolated_jjs(mesh_p1%nis))

         ALLOCATE(mesh%jce_extra(SIZE(mesh_p1%jce_extra, 1), SIZE(mesh_p1%jce_extra, 2)))
         mesh%jce_extra = mesh_p1%jce_extra
         ALLOCATE(mesh%jees(SIZE(mesh_p1%jees)))
         mesh%jees = mesh_p1%jees
         ALLOCATE(mesh%jecs(SIZE(mesh_p1%jecs)))
         mesh%jecs = mesh_p1%jecs

         mesh%gauss%n_w = 10
         mesh%gauss%n_ws = 4
      END IF

      nb_angle = 0
      ALLOCATE(virgin(mesh_p1%medge), j_mid(nw * f_dof, me), jjs_mid(f_dof, mes), r_mid(kd))

      IF (kd == 3) THEN
         WRITE(*, *) ' CREATE_GRID_Pk: 3D case not programmed yet !'
         STOP
      END IF

      DO proc = 1, nb_proc
         IF (mesh_p1%loc_to_glob(1) <= mesh_p1%disp(proc)) THEN
            EXIT
         END IF
      END DO

      !===GENERATION OF THE Pk GRID
      mesh%rr(:, 1:dom_np) = mesh_p1%rr(:, 1:dom_np)
      mesh%rr(:, mesh%dom_np + 1:mesh%dom_np + np - dom_np) = mesh_p1%rr(:, dom_np + 1:)
      mesh%jj(1:nw, :) = mesh_p1%jj
      mesh%jj_extra(1:nw, :) = mesh_p1%jj_extra
      mesh%jjs_extra(1:nws, :) = mesh_p1%jjs_extra
      mesh%rrs_extra(:, 1:nw, :) = mesh_p1%rrs_extra
      mesh%loc_to_glob(1:dom_np) = mesh_p1%loc_to_glob(1:dom_np) &
           + (mesh_p1%disedge(proc) - 1) * f_dof + (mesh_p1%discell(proc) - 1) * (f_dof - 1)
      mesh%isolated_jjs = mesh_p1%isolated_jjs &
           + (mesh_p1%disedge(proc) - 1) * f_dof + (mesh_p1%discell(proc) - 1) * (f_dof - 1)

      DO m = 1, mesh%me
         DO n = 1, nw
            IF (mesh_p1%jj(n, m) > mesh_p1%dom_np) THEN
               mesh%jj(n, m) = mesh_p1%jj(n, m) + mesh_p1%medge * f_dof + mesh_p1%me * (f_dof - 1)
            END IF
         END DO
      END DO

      DO m = 1, np - dom_np
         DO p = 1, nb_proc
            IF (mesh_p1%loc_to_glob(dom_np + m) < mesh_p1%disp(p + 1)) THEN
               EXIT
            END IF
         END DO
         mesh%loc_to_glob(mesh%dom_np + m) = mesh_p1%loc_to_glob(dom_np + m) &
              + (mesh_p1%disedge(p) - 1) * f_dof + (mesh_p1%discell(p) - 1) * (f_dof - 1)
      END DO

      DO m = 1, mesh_p1%mextra
         DO n = 1, nw
            DO p = 1, nb_proc
               IF (mesh_p1%jj_extra(n, m) < mesh_p1%disp(p + 1)) THEN
                  EXIT
               END IF
            END DO
            mesh%jj_extra(n, m) = mesh_p1%jj_extra(n, m) &
                 + (mesh_p1%disedge(p) - 1) * f_dof + (mesh_p1%discell(p) - 1) * (f_dof - 1)
         END DO
      END DO

      DO m = 1, mesh_p1%mes_extra
         DO n = 1, nws
            DO p = 1, nb_proc
               IF (mesh_p1%jjs_extra(n, m) < mesh_p1%disp(p + 1)) THEN
                  EXIT
               END IF
            END DO
            mesh%jjs_extra(n, m) = mesh_p1%jjs_extra(n, m) &
                 + (mesh_p1%disedge(p) - 1) * f_dof + (mesh_p1%discell(p) - 1) * (f_dof - 1)
         END DO
      END DO

      virgin = .TRUE.

      n_dof = 0
      DO m = 1, me !===loop on the elements
         DO k = 1, nw !===loop on the nodes (sides) of the element
            edge_g = mesh_p1%jce(k, m)
            edge_l = edge_g - mesh_p1%disedge(proc) + 1

            IF (edge_l <= 0) CYCLE

            n_new_start = (edge_l - 1) * f_dof + mesh_p1%dom_np

            m_op_k = mesh_p1%neigh(k, m)
            n_k1 = MODULO(k, nw) + 1
            n_k2 = MODULO(k + 1, nw) + 1
            n1 = mesh_p1%jj(n_k1, m)
            n2 = mesh_p1%jj(n_k2, m)

            IF (n_k1<n_k2) THEN !===Go from lowest global index to highest global index
               n_start = n1
               n_end = n2
            ELSE
               n_start = n2
               n_end = n1
            END IF

            iso = .FALSE.
            IF (m_op_k == 0) THEN  !===the side is on the boundary
               DO ms = 1, SIZE(mesh_p1%neighs) + 1
                  IF (ms == SIZE(mesh_p1%neighs) + 1) WRITE(*, *) &
                       'BUG in create_iso_grid: cell near boundary isnt in neighs', m_op_k, m
                  IF (mesh_p1%neighs(ms) == m) EXIT
               END DO
               CALL is_on_curved_interface(mesh_p1%sides(ms), iso, interface)

            END IF

            IF (virgin(edge_l)) THEN !===This side is new
               DO l = 1, f_dof
                  n_dof = n_dof + 1 !===New index created
                  j_mid((k - 1) * f_dof + l, m) = l + n_new_start
                  mesh%rr(:, l + n_new_start) = mesh_p1%rr(:, n_start) &
                       + l * (mesh_p1%rr(:, n_end) - mesh_p1%rr(:, n_start)) / type_fe
                  IF (iso) THEN
                     CALL rescale_to_curved_boundary(mesh%rr(:, l + n_new_start), interface)
                  END IF
                  mesh%loc_to_glob(l + n_new_start) = l + n_new_start + mesh%disp(proc) - 1
               END DO
            ELSE !===the side has been already considered
               mm = m_op_k
               DO i = 1, nw
                  IF (mesh_p1%neigh(i, mm) == m) THEN
                     kk = i
                     EXIT
                  END IF
               ENDDO
               DO l = 1, f_dof
                  j_mid((k - 1) * f_dof + l, m) = j_mid((kk - 1) * f_dof + l, mm) !===New index created
               END DO
            ENDIF
            virgin(edge_l) = .FALSE.
         ENDDO
      ENDDO

      IF (n_dof /= mesh_p1%medge * f_dof) THEN
         WRITE(*, *) 'BUG in create_iso_grid_distributed, n_dof /= mesh_p1%medge * f_dof'
         STOP
      END IF

      n_dof = 0
      DO edges = 1, mesh_p1%medges
         edge_g = mesh_p1%jees(edges)
         m = mesh_p1%jecs(edges)
         DO k = 1, nw
            IF (mesh_p1%jce(k, m) == edge_g) THEN
               EXIT
            END IF
         ENDDO

         DO p = 1, nb_proc
            IF (edge_g < mesh_p1%disedge(p + 1)) THEN
               EXIT
            END IF
         END DO

         n_k1 = MODULO(k, nw) + 1
         n_k2 = MODULO(k + 1, nw) + 1
         n1 = mesh_p1%jj(n_k1, m)
         n2 = mesh_p1%jj(n_k2, m)
         IF (n_k1<n_k2) THEN !===Go from lowest global index to highest global index
            n_start = n1
            n_end = n2
         ELSE
            n_start = n2
            n_end = n1
         END IF

         edge_l = edge_g - mesh_p1%disedge(p) + 1
         n_new_start = (edges - 1) * f_dof + mesh_p1%me * (f_dof - 1) + mesh_p1%medge * f_dof + mesh_p1%np

         DO l = 1, f_dof
            n_dof = n_dof + 1 !===New index created
            j_mid((k - 1) * f_dof + l, m) = l + n_new_start
            mesh%rr(:, n_new_start + l) = mesh_p1%rr(:, n_start) &
                 + l * (mesh_p1%rr(:, n_end) - mesh_p1%rr(:, n_start)) / type_fe
            mesh%loc_to_glob(n_new_start + l) = l + (edge_l - 1) * f_dof + mesh_p1%domnp(p) + mesh%disp(p) - 1
         END DO

      END DO

      IF (n_dof /= mesh_p1%medges * f_dof) THEN
         WRITE(*, *) 'BUG in create_iso_grid_distributed, n_dof /= mesh_p1%medge * f_dof'
         STOP
      END IF

      n_dof = 0
      !===connectivity array for iso grid
      DO m = 1, me
         DO n = 1, f_dof * nw
            mesh%jj(nw + n, m) = j_mid(n, m)
         END DO
         IF (type_fe==3) THEN
            n_dof = n_dof + 1
            n_new_start = m + mesh_p1%dom_np + mesh_p1%medge * 2
            mesh%jj(10, m) = n_new_start
            mesh%rr(:, n_new_start) = &
                 (mesh_p1%rr(:, mesh_p1%jj(1, m)) + mesh_p1%rr(:, mesh_p1%jj(2, m)) + mesh_p1%rr(:, mesh_p1%jj(3, m))) / 3
            mesh%loc_to_glob(n_new_start) = m + mesh_p1%medge * 2 + mesh_p1%dom_np + mesh%disp(proc) - 1
         END IF
      END DO

      IF (type_fe == 3 .and. n_dof /= mesh_p1%me) THEN
         WRITE(*, *) 'BUG in create_iso_grid_distributed, type_fe == 3 .and. n_dof /= mesh_p1%me'
         STOP
      END IF

      !==connectivity array the surface elements of the iso grid
      DO ms = 1, mes
         m = mesh_p1%neighs(ms)
         DO n = 1, kd + 1 !===Simplices
            IF (MINVAL(ABS(mesh_p1%jj(n, m) - mesh_p1%jjs(:, ms)))/=0) THEN
               kk = n
               EXIT
            END IF
         ENDDO
         DO l = 1, f_dof
            jjs_mid(l, ms) = j_mid((kk - 1) * f_dof + l, m) !===New index created
         END DO
      ENDDO
      mesh%jjs(1:nws, :) = mesh_p1%jjs
      DO i = 1, SIZE(mesh%jjs, 2)
         DO n = 1, nws
            IF (mesh%jjs(n, i) > mesh_p1%dom_np) THEN
               mesh%jjs(n, i) = mesh_p1%jjs(n, i) + mesh_p1%medge * f_dof + mesh_p1%me * (f_dof - 1)
            END IF
         END DO
      END DO

      mesh%jjs(nws + 1:, :) = jjs_mid

      DEALLOCATE(virgin, j_mid, jjs_mid, r_mid)

      !===new vertices on extra cells
      DO m = 1, mesh%mextra
         DO k = 1, nw !===loop on the nodes (sides) of the element
            edge_g = mesh_p1%jce_extra(k, m)
            cell_g = mesh_p1%jcc_extra(m)
            DO p_e = 1, nb_proc
               IF (edge_g < mesh_p1%disedge(p_e + 1)) THEN
                  EXIT
               END IF
            END DO

            DO p_c = 1, nb_proc
               IF (cell_g < mesh_p1%discell(p_c + 1)) THEN
                  EXIT
               END IF
            END DO
            edge_l = edge_g - mesh_p1%disedge(p_e) + 1
            cell_l = cell_g - mesh_p1%discell(p_c) + 1

            DO l = 1, f_dof
               mesh%jj_extra(nw + (k - 1) * f_dof + l, m) = l &
                    + (edge_l - 1) * f_dof + mesh_p1%domnp(p_e) + mesh%disp(p_e) - 1
            END DO

            IF (type_fe==3) THEN
               mesh%jj_extra(10, m) = cell_l + mesh_p1%domedge(p_c) * 2 + mesh_p1%domnp(p_c) + mesh%disp(p_c) - 1
            END IF
         END DO
      END DO

      !==connectivity array the surface elements of the iso grid for extras
      DO ms = 1, mesh%mes_extra
         iso = .FALSE.
         CALL is_on_curved_interface(mesh%sides_extra(ms), iso, interface)

         cell_g = mesh%neighs_extra(ms)
         DO m = 1, mesh%mextra !find associated extra cell
            IF (mesh_p1%jcc_extra(m) == cell_g) EXIT
         END DO
         DO n = 1, kd + 1 !===find side in cell
            IF (MINVAL(ABS(mesh%jj_extra(n, m) - mesh_p1%jjs_extra(:, ms)))/=0) THEN
               kk = n
               EXIT
            END IF
         ENDDO

         DO l = 1, f_dof
            mesh%jjs_extra(nws + l, ms) = mesh%jj_extra(nw + (kk - 1) * f_dof + l, m)
         END DO

         DO k = 1, kd + 1
            n_k1 = MODULO(k, nw) + 1
            n_k2 = MODULO(k + 1, nw) + 1
            IF (n_k1<n_k2) THEN !===Go from lowest index to highest index
               n_start = n_k1
               n_end = n_k2
            ELSE
               n_start = n_k2
               n_end = n_k1
            END IF

            DO l = 1, f_dof
               mesh%rrs_extra(:, nw + (k - 1) * f_dof + l, ms) = mesh_p1%rrs_extra(:, n_start, ms) &
                    + l * (mesh_p1%rrs_extra(:, n_end, ms) - mesh_p1%rrs_extra(:, n_start, ms)) / type_fe
               IF (iso) THEN
                  CALL rescale_to_curved_boundary(mesh%rrs_extra(:, nw + (k - 1) * f_dof + l, ms), interface)
               END IF
            END DO
         END DO

         IF (type_fe==3) THEN
            mesh%rrs_extra(:, 10, ms) = &
                 (mesh_p1%rrs_extra(:, 1, ms) + mesh_p1%rrs_extra(:, 2, ms) + mesh_p1%rrs_extra(:, 3, ms)) / 3
         END IF
      ENDDO
   END SUBROUTINE  create_iso_grid_distributed

   SUBROUTINE refinement_iso_grid_distributed(mesh_p1)
      !===jj(:, :)    nodes of the  volume_elements of the input grid
      !===jjs(:, :)    nodes of the surface_elements of the input grid
      !===rr(:, :)    cartesian coordinates of the nodes of the input grid
      !===m_op(:,:)   volume element opposite to each node
      !===neigh_el(:) volume element ajacent to the surface element
      !===jj_f(:, :)  nodes of the  volume_elements of the output p2 grid
      !===jjs_f(:, :)  nodes of the surface_elements of the output p2 grid
      !===rr_f(:, :)  cartesian coordinates of the nodes of the output p2 grid
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_p1, mesh
      LOGICAL, DIMENSION(:), ALLOCATABLE :: virgin
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: j_mid, jjs_mid
      INTEGER :: np, me, mes, nw, nws, kd, n, m, k, l, n_dof, dom_np, e_g, a, mextra
      INTEGER :: n1, n2, n3, n4, ms, n_start, n_end, n1_g, n2_g, neigh, k_neigh, n_kneigh1, n_kneigh2, swap
      INTEGER :: n_k1, n_k2, m_op_k, kk, i, mm, ms_bord, p_e, p_c, m_new, e_k, p_j, n_kks
      REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: r_mid
      INTEGER, DIMENSION(3) :: edges_g, edges_l, p_es
      INTEGER, DIMENSION(2) :: n_ks
      REAL(KIND = 8) :: epsilon = 1.d-13, dist, d1, d2, s1, s2, s3, shalf, ref, scc, infinity
      INTEGER :: ns, ns1, index, nb_angle, f_dof, edge_g, edge_l, n_new_start, proc, nb_proc, edges, p, cell_g, cell_l
      INTEGER :: m1, m2, interface, m_center, tab1, tab2, mes_int
      LOGICAL :: iso
      INTEGER, SAVE :: count = 0

      IF (mesh_p1%me == 0) THEN
         RETURN
      END IF

      nw = SIZE(mesh_p1%jj, 1)   !===nodes in each volume element (3 in 2D)
      me = SIZE(mesh_p1%jj, 2)   !===number of cells
      kd = SIZE(mesh_p1%rr, 1)   !===space dimensions
      np = mesh_p1%np            !===number of P1 vertices connected to grid
      dom_np = mesh_p1%dom_np    !===number of P1 vertices attributed to proc
      mes = SIZE(mesh_p1%jjs, 2)
      mes_int = SIZE(mesh_p1%jjs_int, 2)
      nws = SIZE(mesh_p1%jjs, 1)
      nb_proc = SIZE(mesh_p1%domnp)

      mesh%me = 4 * mesh_p1%me
      mesh%mes = 2 * mesh_p1%mes !---> something with news to take into account
      mesh%mes_int = 2 * mesh_p1%mes_int
      mesh%np = mesh_p1%np + mesh_p1%medge + mesh_p1%medges
      mesh%medge = 2 * mesh_p1%medge + 3 * mesh_p1%me
      mesh%medges = 2 * mesh_p1%medges
      ALLOCATE(mesh%jj(nw, mesh%me)) !--->done
      ALLOCATE(mesh%jjs(nws, mesh%mes))  !--->done
      ALLOCATE(mesh%rr(kd, mesh%np)) !--->done
      ALLOCATE(mesh%loc_to_glob(mesh%np)) !--->done

      ALLOCATE(mesh%jce(nw, mesh%me))  !--->done
      !ALLOCATE(mesh%jev(nw - 1, mesh%medge)) !----->never needed so not constructed
      ALLOCATE(mesh%jees(mesh%medges))  !--->done
      ALLOCATE(mesh%jecs(mesh%medges))  !--->done

      ALLOCATE(mesh%neigh(nw, mesh%me)) !--->done
      ALLOCATE(mesh%sides(mesh%mes))   !---> still don't know what this is ?? done ?
      ALLOCATE(mesh%neighs(mesh%mes))  !--->done
      ALLOCATE(mesh%i_d(mesh%me)) !--->done

      ALLOCATE(mesh%jjs_int(nws, mesh%mes_int))  !--->done
      ALLOCATE(mesh%neighs_int(2, mesh%mes_int)) !--->done
      ALLOCATE(mesh%sides_int(mesh%mes_int))

      mesh%dom_me = 4 * mesh_p1%dom_me
      mesh%dom_np = mesh_p1%dom_np + mesh_p1%medge
      mesh%dom_mes = 2 * mesh_p1%dom_mes
      ALLOCATE(mesh%disp(nb_proc + 1), mesh%domnp(nb_proc))
      mesh%domnp = mesh_p1%domnp + mesh_p1%domedge
      mesh%disp = mesh_p1%disp + mesh_p1%disedge - 1
      ALLOCATE(mesh%discell(nb_proc + 1), mesh%domcell(nb_proc))
      mesh%domcell = 4 * mesh_p1%domcell
      mesh%discell = 4 * mesh_p1%discell - 3
      ALLOCATE(mesh%disedge(nb_proc + 1), mesh%domedge(nb_proc))
      mesh%domedge = 2 * mesh_p1%domedge + 3 * mesh_p1%domcell
      mesh%disedge = 2 * mesh_p1%disedge + 3 * mesh_p1%discell - 4

      mesh%nis = mesh_p1%nis
      ALLOCATE(mesh%isolated_interfaces(mesh_p1%nis, 2))
      mesh%isolated_interfaces = mesh_p1%isolated_interfaces
      ALLOCATE(mesh%isolated_jjs(mesh_p1%nis))

      IF (kd == 3) THEN
         WRITE(*, *) ' CREATE_GRID_Pk: 3D case not programmed yet !'
         STOP
      END IF

      DO proc = 1, nb_proc
         IF (mesh_p1%loc_to_glob(1) <= mesh_p1%disp(proc))    EXIT
      END DO

      !===GENERATION OF THE Pk GRID
      mesh%rr(:, 1:dom_np) = mesh_p1%rr(:, 1:dom_np)
      mesh%rr(:, mesh%dom_np + 1:mesh%dom_np + np - dom_np) = mesh_p1%rr(:, dom_np + 1:)
      mesh%loc_to_glob(1:dom_np) = mesh_p1%loc_to_glob(1:dom_np) + mesh_p1%disedge(proc) - 1
      mesh%isolated_jjs = mesh_p1%isolated_jjs + mesh_p1%disedge(proc) - 1

      DO m = 1, np - dom_np
         DO p = 1, nb_proc
            IF (mesh_p1%loc_to_glob(dom_np + m) < mesh_p1%disp(p + 1)) THEN
               EXIT
            END IF
         END DO
         mesh%loc_to_glob(mesh%dom_np + m) = mesh_p1%loc_to_glob(dom_np + m) + mesh_p1%disedge(p) - 1
      END DO

      n_dof = 0
      DO m = 1, me !===loop on the elements unrefined mesh
         edges_g = mesh_p1%jce(:, m)
         edges_l = edges_g - mesh_p1%disedge(proc) + 1
         m_center = 4 * (m - 1) + 1
         !===Center cell
         mesh%i_d(m_center) = mesh_p1%i_d(m)
         DO i = 1, nw
            !===Creating the new points
            IF (edges_l(i) <=0) THEN
               DO p = 1, nb_proc
                  IF (edges_g(i) < mesh_p1%disedge(p + 1)) EXIT
               END DO
               DO ms = 1, mesh_p1%medges
                  IF (mesh_p1%jees(ms) == edges_g(i)) EXIT
               END DO
               mesh%jj(i, m_center) = mesh%dom_np + mesh_p1%np - mesh_p1%dom_np + ms
               mesh%loc_to_glob(mesh%jj(i, m_center)) = mesh%disp(p) - 1 + mesh_p1%domnp(p) + edges_g(i) - mesh_p1%disedge(p) + 1
            ELSE
               mesh%jj(i, m_center) = dom_np + edges_l(i)
               mesh%loc_to_glob(mesh%jj(i, m_center)) = mesh%disp(proc) - 1 + mesh%jj(i, m_center)
            END IF

            n1 = mesh_p1%jj(MODULO(i, nw) + 1, m)
            n2 = mesh_p1%jj(MODULO(i + 1, nw) + 1, m)
            mesh%rr(:, mesh%jj(i, m_center)) = (mesh_p1%rr(:, n2) + mesh_p1%rr(:, n1)) / 2.d0

            !===Setting up neighbours and edges
            mesh%neigh(i, m_center) = m_center + i
            mesh%i_d(mesh%neigh(i, m_center)) = mesh_p1%i_d(m)
            mesh%jce(i, m_center) = mesh%disedge(proc) - 1 + 2 * mesh_p1%medge + 3 * (m - 1) + i
         END DO


         !===Corner cells
         DO k = 1, nw
            m_new = mesh%neigh(k, m_center)

            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
            IF (n_ks(1)>n_ks(2)) THEN
               n_ks = (/n_ks(2), n_ks(1)/)
            END IF

            !===Setting the center cell as a neighbour
            mesh%neigh(1, m_new) = m_center
            mesh%jce(1, m_new) = mesh%jce(k, m_center)

            !===Adding the points using the center cell
            IF (mesh_p1%jj(k, m) > dom_np) THEN
               mesh%jj(1, m_new) = mesh_p1%jj(k, m) + mesh%dom_np - mesh_p1%dom_np
            ELSE
               mesh%jj(1, m_new) = mesh_p1%jj(k, m)
            END IF
            mesh%jj(2, m_new) = mesh%jj(n_ks(1), m_center)
            mesh%jj(3, m_new) = mesh%jj(n_ks(2), m_center)
            !===Adding the last edges and neighbours
            DO e_k = 1, 2

               !===Adding neighbours
               neigh = mesh_p1%neigh(n_ks(MODULO(e_k, 2) + 1), m)

               IF (neigh <= 0) THEN
                  mesh%neigh(e_k + 1, m_new) = neigh
               ELSE
                  DO k_neigh = 1, nw
                     IF (mesh_p1%neigh(k_neigh, neigh) == m) exit
                  END DO
                  n_kneigh1 = MODULO(k_neigh, nw) + 1
                  n_kneigh2 = MODULO(k_neigh + 1, nw) + 1
                  IF (n_kneigh1>n_kneigh2) THEN
                     swap = n_kneigh1
                     n_kneigh1 = n_kneigh2
                     n_kneigh2 = swap
                  END IF
                  IF (k < n_ks(e_k)) THEN
                     mesh%neigh(e_k + 1, m_new) = 4 * (neigh - 1) + 1 + n_kneigh1
                  ELSE
                     mesh%neigh(e_k + 1, m_new) = 4 * (neigh - 1) + 1 + n_kneigh2
                  END IF
               END IF

               !===Adding edge
               e_g = mesh_p1%jce(n_ks(MODULO(e_k, 2) + 1), m)
               IF (e_g < mesh_p1%disedge(proc)) THEN
                  DO p = 1, nb_proc
                     IF (e_g < mesh_p1%disedge(p + 1)) EXIT
                  END DO
                  DO ms = 1, mesh_p1%medges
                     IF (mesh_p1%jees(ms) == e_g) EXIT
                  END DO
               ELSE
                  p = proc
               END IF

               IF (k < n_ks(e_k)) THEN
                  mesh%jce(e_k + 1, m_new) = mesh%disedge(p) - 1 + e_g - mesh_p1%disedge(p) + 1
               ELSE
                  mesh%jce(e_k + 1, m_new) = mesh%disedge(p) - 1 + mesh_p1%domedge(p) + e_g - mesh_p1%disedge(p) + 1
               END IF

               IF (e_g < mesh_p1%disedge(proc)) THEN
                  IF (k < n_ks(e_k)) THEN
                     mesh%jees(ms) = mesh%jce(e_k + 1, m_new)
                     mesh%jecs(ms) = m_new
                  ELSE
                     mesh%jees(mesh_p1%medges + ms) = mesh%jce(e_k + 1, m_new)
                     mesh%jecs(mesh_p1%medges + ms) = m_new
                  END IF
               END IF
            END DO
         END DO
      END DO

      !===Surface elements
      DO ms = 1, mes
         m = mesh_p1%neighs(ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_p1%jj(k, m) - mesh_p1%jjs(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         mesh%jjs(1, ms) = mesh_p1%jj(n_ks(1), m)
         IF (mesh%jjs(1, ms) > mesh_p1%dom_np) mesh%jjs(1, ms) = mesh%jjs(1, ms) + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs(1, mes + ms) = mesh_p1%jj(n_ks(2), m)
         IF (mesh%jjs(1, mes + ms) > mesh_p1%dom_np) mesh%jjs(1, mes + ms) = mesh%jjs(1, mes + ms) &
              + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs(2, ms) = mesh%jj(k, 4 * (m - 1) + 1)
         mesh%jjs(2, mes + ms) = mesh%jj(k, 4 * (m - 1) + 1)
         mesh%neighs(ms) = mesh%neigh(n_ks(1), 4 * (m - 1) + 1)
         mesh%neighs(mes + ms) = mesh%neigh(n_ks(2), 4 * (m - 1) + 1)
         mesh%sides(ms) = mesh_p1%sides(ms)
         mesh%sides(mes + ms) = mesh_p1%sides(ms)

         CALL is_on_curved_interface(mesh_p1%sides(ms), iso, interface)
         IF (iso) THEN
            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         END IF
      ENDDO

      !===Internal surface elements
      DO ms = 1, mes_int
         m = mesh_p1%neighs_int(1, ms)
         !===Finding the corresponding side in the cell
         DO k = 1, nw
            IF (MINVAL(ABS(mesh_p1%jj(k, m) - mesh_p1%jjs_int(:, ms)))/=0) EXIT
         ENDDO

         n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         mesh%jjs_int(1, ms) = mesh_p1%jj(n_ks(1), m)
         IF (mesh%jjs_int(1, ms) > mesh_p1%dom_np) mesh%jjs_int(1, ms) = mesh%jjs_int(1, ms) + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs_int(1, mes_int + ms) = mesh_p1%jj(n_ks(2), m)
         IF (mesh%jjs_int(1, mes_int + ms) > mesh_p1%dom_np) mesh%jjs_int(1, mes_int + ms) = mesh%jjs_int(1, mes_int + ms) &
              + mesh%dom_np - mesh_p1%dom_np
         mesh%jjs_int(2, ms) = mesh%jj(k, 4 * (m - 1) + 1)
         mesh%jjs_int(2, mes_int + ms) = mesh%jj(k, 4 * (m - 1) + 1)
         mesh%neighs_int(1, ms) = mesh%neigh(n_ks(1), 4 * (m - 1) + 1)
         mesh%neighs_int(1, mes_int + ms) = mesh%neigh(n_ks(2), 4 * (m - 1) + 1)
         mesh%sides_int(ms) = mesh_p1%sides_int(ms)
         mesh%sides_int(mes_int + ms) = mesh_p1%sides_int(ms)

         CALL is_on_curved_interface(mesh_p1%sides_int(ms), iso, interface)
         IF (iso) THEN
            CALL rescale_to_curved_boundary(mesh%rr(:, mesh%jj(k, 4 * (m - 1) + 1)), interface)
         END IF
      ENDDO

      !===Counting number of new extra cells
      mesh%mextra = 0
      !===Need to rework that to do it the smart way and the update conditions when constructing cells
      !      DO m = 1, mesh_p1%mextra
      !         a = 0
      !         DO k = 1, nw
      !            IF (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(k, m) &
      !                 .and. mesh_p1%jce_extra(k, m) < mesh_p1%disedge(proc + 1) .and. a == 0) THEN
      !               a = 1
      !               mesh%mextra = mesh%mextra + 1
      !            END IF
      !         END DO
      !         DO k = 1, nw
      !            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
      !            IF (n_ks(1)>n_ks(2)) THEN
      !               n_ks = (/n_ks(2), n_ks(1)/)
      !            END IF
      !
      !            IF (mesh_p1%jj_extra(k, m) < mesh_p1%disp(proc + 1) .and. &
      !                 (mesh_p1%disp(proc) <= mesh_p1%jj_extra(k, m)  .or. &
      !                      (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(n_ks(1), m) &
      !                           .and. mesh_p1%jce_extra(n_ks(1), m) < mesh_p1%disedge(proc + 1)) .or.  &
      !                      (mesh_p1%disedge(proc) <= mesh_p1%jce_extra(n_ks(2), m) &
      !                           .and. mesh_p1%jce_extra(n_ks(2), m) < mesh_p1%disedge(proc + 1)))) THEN
      !               mesh%mextra = mesh%mextra + 1
      !            END IF
      !         END DO
      !      END DO
      !===In the mean time
      mesh%mextra = mesh_p1%mextra * 4
      ALLOCATE(mesh%jj_extra(nw, mesh%mextra)) !---->
      ALLOCATE(mesh%jce_extra(nw, mesh%mextra)) !---->
      ALLOCATE(mesh%jcc_extra(mesh%mextra)) !---->

      !===Constructing the extra cells
      mextra = 0
      DO m = 1, mesh_p1%mextra

         !center cell
         mextra = mextra + 1

         cell_g = mesh_p1%jcc_extra(m)
         DO p_c = 1, nb_proc
            IF (cell_g < mesh_p1%discell(p_c + 1)) EXIT
         END DO
         cell_l = cell_g - mesh_p1%discell(p_c) + 1
         mesh%jcc_extra(mextra) = 4 * (cell_l - 1) + 1 + mesh%discell(p_c) - 1

         DO i = 1, 3
            mesh%jce_extra(i, mextra) = mesh%disedge(p_c) - 1 + 2 * mesh_p1%domedge(p_c) + 3 * (cell_l - 1) + i
            e_g = mesh_p1%jce_extra(i, m)
            DO p = 1, nb_proc
               IF (e_g < mesh_p1%disedge(p + 1)) EXIT
            END DO
            mesh%jj_extra(i, mextra) = mesh%disp(p) - 1 + mesh_p1%domnp(p) + e_g - mesh_p1%disedge(p) + 1
         END DO

         !corner cells
         DO k = 1, nw
            n_ks = (/MODULO(k, nw) + 1, MODULO(k + 1, nw) + 1/)
            IF (n_ks(1)>n_ks(2)) THEN
               n_ks = (/n_ks(2), n_ks(1)/)
            END IF

            mextra = mextra + 1
            cell_g = mesh_p1%jcc_extra(m)
            DO p_c = 1, nb_proc
               IF (cell_g < mesh_p1%discell(p_c + 1))  EXIT
            END DO

            cell_l = cell_g - mesh_p1%discell(p_c) + 1
            mesh%jcc_extra(mextra) = mesh%discell(p_c) - 1 + 4 * (cell_l - 1) + 1 + k

            edges_g = mesh_p1%jce_extra(:, m)

            DO i = 1, nw
               DO p = 1, nb_proc
                  IF (edges_g(i) < mesh_p1%disedge(p + 1)) EXIT
               END DO
               p_es(i) = p
            END DO

            !===Adding the points
            DO p_j = 1, nb_proc
               IF (mesh_p1%jj_extra(k, m) < mesh_p1%disp(p_j + 1))  EXIT
            END DO
            mesh%jj_extra(1, mextra) = mesh_p1%jj_extra(k, m) + mesh_p1%disedge(p_j) - 1
            mesh%jj_extra(2, mextra) = mesh%disp(p_es(n_ks(1))) - 1 + &
                 mesh_p1%domnp(p_es(n_ks(1))) + edges_g(n_ks(1)) - mesh_p1%disedge(p_es(n_ks(1))) + 1
            mesh%jj_extra(3, mextra) = mesh%disp(p_es(n_ks(2))) - 1 + &
                 mesh_p1%domnp(p_es(n_ks(2))) + edges_g(n_ks(2)) - mesh_p1%disedge(p_es(n_ks(2))) + 1

            !===Adding the edges
            mesh%jce_extra(1, mextra) = mesh%disedge(p_c) - 1 + 2 * mesh_p1%domedge(p_c) + 3 * (cell_l - 1) + k

            DO e_k = 1, 2
               IF (k < n_ks(e_k)) THEN
                  mesh%jce_extra(e_k + 1, mextra) = mesh%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) - 1 + &
                       edges_g(n_ks(MODULO(e_k, 2) + 1)) - mesh_p1%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + 1
               ELSE
                  mesh%jce_extra(e_k + 1, mextra) = mesh%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) - 1 + &
                       mesh_p1%domedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + edges_g(n_ks(MODULO(e_k, 2) + 1)) &
                       - mesh_p1%disedge(p_es(n_ks(MODULO(e_k, 2) + 1))) + 1
               END IF
            END DO
         END DO
      END DO

      !===Constructing the extra cells at interfaces
      mesh%mes_extra = 2 * mesh_p1%mes_extra
      ALLOCATE(mesh%jjs_extra(nws, mesh%mes_extra))
      ALLOCATE(mesh%rrs_extra(2, nw, mesh%mes_extra))
      ALLOCATE(mesh%sides_extra(mesh%mes_extra), mesh%neighs_extra(mesh%mes_extra))

      mextra = 0
      DO m = 1, mesh_p1%mes_extra
         CALL is_on_curved_interface(mesh_p1%sides_extra(m), iso, interface)

         cell_g = mesh_p1%neighs_extra(m)
         DO m1 = 1, mesh_p1%mextra !find associated extra cell
            IF (mesh_p1%jcc_extra(m1) == cell_g) EXIT
         END DO

         DO p_c = 1, nb_proc
            IF (cell_g < mesh_p1%discell(p_c + 1))  EXIT
         END DO
         cell_l = cell_g - mesh_p1%discell(p_c) + 1

         DO n = 1, 3 !===find side in cell
            IF (MINVAL(ABS(mesh_p1%jj_extra(n, m1) - mesh_p1%jjs_extra(:, m)))/=0) THEN
               EXIT
            END IF
            IF (n == 3) write(*, *) 'BUG in refinement : didnt find face in extra cell for extra edge'
         ENDDO

         !==cell index of edge
         n_ks = (/MODULO(n, nw) + 1, MODULO(n + 1, nw) + 1/)
         IF (n_ks(1)>n_ks(2)) THEN
            n_ks = (/n_ks(2), n_ks(1)/)
         END IF

         DO k = 1, nws
            mextra = mextra + 1
            mesh%sides_extra(mextra) = mesh_p1%sides_extra(m)
            mesh%neighs_extra(mextra) = mesh%discell(p_c) - 1 + 4 * (cell_l - 1) + 1 + n_ks(k)
            DO m2 = 1, mesh%mextra !find associated extra cell
               IF (mesh%jcc_extra(m2) == mesh%neighs_extra(mextra)) EXIT
            END DO

            DO p_j = 1, nb_proc
               IF (mesh_p1%jj_extra(n_ks(k), m) < mesh_p1%disp(p_j + 1))  EXIT
            END DO

            mesh%jjs_extra(1, mextra) = mesh%jj_extra(1, m2)
            mesh%rrs_extra(:, 1, mextra) = mesh_p1%rrs_extra(:, n_ks(k), m)

            IF (n == 1) THEN
               tab1 = 2
               tab2 = 3
            ELSE IF (n == 3) THEN
               tab1 = 3
               tab2 = 2
            ELSE
               IF (k == 1) THEN
                  tab1 = 2
                  tab2 = 3
               ELSE
                  tab1 = 3
                  tab2 = 2
               END IF
            END IF
            mesh%jjs_extra(2, mextra) = mesh%jj_extra(tab1, m2)
            mesh%rrs_extra(:, tab1, mextra) = (mesh_p1%rrs_extra(:, n_ks(1), m) + mesh_p1%rrs_extra(:, n_ks(2), m)) / 2
            IF (iso) THEN
               CALL rescale_to_curved_boundary(mesh%rrs_extra(:, tab1, mextra), interface)
            END IF
            mesh%rrs_extra(:, tab2, mextra) = (mesh_p1%rrs_extra(:, n_ks(k), m) + mesh_p1%rrs_extra(:, n, m)) / 2
         END DO
      END DO


      !jjs_extra !(extra layer of cells not own by proc but with dofs own by proc)
      !rrs_extra  ! coordinates for cells at interfaces
      !sides_extra, neighs_extra !interfaces
      !mes_extra

      CALL free_mesh(mesh_p1)
      CALL copy_mesh(mesh, mesh_p1)
      CALL free_mesh(mesh)
   END SUBROUTINE refinement_iso_grid_distributed

   SUBROUTINE is_on_curved_interface(side, iso, interface)
      USE input_data
      INTEGER :: side, interface
      LOGICAL :: iso
      interface = -1
      iso = .FALSE.

      IF (inputs%nb_spherical + inputs%nb_curved > 0) THEN
         IF (MINVAL(ABS(side - inputs%list_spherical)) == 0 .OR. &
              MINVAL(ABS(side - inputs%list_curved)) == 0) THEN
            DO interface = 1, inputs%nb_spherical + inputs%nb_curved
               IF (interface <= inputs%nb_spherical) THEN
                  IF (side - inputs%list_spherical(interface) == 0) EXIT
               ELSE
                  IF (side - inputs%list_curved(interface - inputs%nb_spherical) == 0) EXIT
               END IF
            END DO
            iso = .TRUE.
         ELSE
            iso = .FALSE.
         END IF
      ELSE
         iso = .FALSE.
      END IF

   END SUBROUTINE is_on_curved_interface

   SUBROUTINE rescale_to_curved_boundary(rr, interface)
      USE input_data
      USE boundary
      REAL(KIND = 8), DIMENSION(2) :: rr, rr_ref
      INTEGER :: interface
      REAL(KIND = 8) :: rescale, pi = ACOS(-1.d0), theta
      IF (interface <= inputs%nb_spherical) THEN
         rr_ref = rr - inputs%origin_spherical(:, interface)
         rescale = inputs%radius_spherical(interface) / SQRT(SUM(rr_ref * rr_ref))
         rr = rr_ref * rescale + inputs%origin_spherical(:, interface)
      ELSE
         rr_ref = rr - inputs%origin_curved(:, interface - inputs%nb_spherical)
         theta = pi - pi / 2 * (1 + sgn(rr_ref(1))) * (1 - sgn(rr_ref(2) * rr_ref(2))) &
              - pi / 4 * (2 + sgn(rr_ref(1))) * sgn(rr_ref(2)) &
              - sgn(rr_ref(1) * rr_ref(2)) * ATAN((ABS(rr_ref(1)) - ABS(rr_ref(2))) / (ABS(rr_ref(1)) + ABS(rr_ref(2))))
         rescale = curved_boundary_radius(inputs%list_curved(interface - inputs%nb_spherical), theta) &
              / SQRT(SUM(rr_ref * rr_ref))
         rr = rr_ref * rescale + inputs%origin_curved(:, interface - inputs%nb_spherical)
      END IF

   END SUBROUTINE rescale_to_curved_boundary

   FUNCTION sgn(x) RESULT(out)
      REAL(KIND = 8) :: x, out
      IF (x > 0.d0) THEN
         out = 1.d0
      ELSE IF (x < 0.d0) THEN
         out = -1.d0
      ELSE
         out = 0.d0
      ENDIF
   END FUNCTION sgn

END MODULE refine_mesh