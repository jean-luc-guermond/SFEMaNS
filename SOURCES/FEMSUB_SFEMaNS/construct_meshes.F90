MODULE construct_meshes
   USE def_type_mesh
   USE HCT_mesh
   USE dir_nodes

   INTEGER :: max_nnz_per_row
CONTAINS

   SUBROUTINE construct_distributed_mesh
#include "petsc/finclude/petsc.h"
      USE petsc
      USE input_data
      USE prep_maill
      USE mod_gauss_points_2d
      USE metis_reorder_elements
      USE st_matrix
      USE dir_nodes_petsc
      USE my_util
      USE Dir_nodes
      USE sub_plot
      USE pressure_basis_change
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh_glob, mesh_p1, mesh_p1_r, HCT_mesh_p1
      INTEGER :: dom_np, n, m
      INTEGER, DIMENSION(:), ALLOCATABLE :: ranks
      REAL(KIND = 8) :: t1
      CHARACTER(len = 1) :: tit
      !TESTTTTTT
      REAL(KIND = 8) :: n_s(3, 4)
      INTEGER :: i, neigh, k
      !TESTTTTTT
      PetscErrorCode :: ierr
      PetscMPIInt    :: nb_proc, rank

      !===Number of procs
      CALL MPI_Comm_size(PETSC_COMM_WORLD, nb_proc, ierr)
      CALL MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)

      !===Prepare the grid
      t1 = user_time()

      CALL load_mesh_free_format_global_p1(inputs%directory, inputs%file_name, inputs%list_dom, &
           mesh_glob, inputs%if_mesh_formatted)
      IF (rank==0) WRITE(*, *) ' time load_mesh_free_format_iso', (user_time() - t1)

      !===Metis reorganizes the mesh
      t1 = user_time()
      CALL reorder_mesh(PETSC_COMM_WORLD, nb_proc, mesh_glob, mesh_p1)
      IF (rank==0) WRITE(*, *) ' reorder_mesh', (user_time() - t1)
      CALL deallocate_mesh(mesh_glob)

      !     ALLOCATE(ranks(mesh_p1%np))
      !     ranks = rank
      !     CALL plot_scalar_field(mesh_p1%jj, mesh_p1%rr, 1.d0 * ranks, 'mesh_p1_r0.plt')
      !     DEALLOCATE(ranks)
      IF (rank==0) WRITE(*, *) 'refinement domnp 0', mesh_p1%domnp
      DO n = 1, inputs%nb_refinement
         !===Create refined mesh
         t1 = user_time()
         CALL refinement_iso_grid_distributed(mesh_p1, mesh_p1_r)
         IF (rank==0) WRITE(*, *) 'refinement', (user_time() - t1)
         IF (rank==0) WRITE(*, *) 'refinement domnp', mesh_p1_r%domnp
         CALL deallocate_mesh(mesh_p1)
         CALL copy_mesh(mesh_p1_r, mesh_p1)
         CALL deallocate_mesh(mesh_p1_r)
         WRITE(tit, '(i1)') n
      END DO

      IF(inputs%if_pp_dg) THEN
         CALL HCT_iso_grid_distributed(mesh_p1, HCT_mesh_p1)
         CALL deallocate_mesh(mesh_p1)
         CALL copy_mesh(HCT_mesh_p1, mesh_p1)
         CALL deallocate_mesh(HCT_mesh_p1)
      END IF

      !===Create higher-order shape functions (P2 or P3)
      !===And setting up basis change for pressure
      t1 = user_time()
      CALL create_iso_grid_distributed(mesh_p1, uu_mesh, inputs%type_fe)
      IF (inputs%if_mesh_same) THEN
         CALL create_iso_grid_distributed(mesh_p1, pp_mesh, inputs%type_fe)
         CALL set_pressure_basis_change(inputs%type_fe, inputs%type_fe)
      ELSE
         CALL create_iso_grid_distributed(mesh_p1, pp_mesh, inputs%type_fe - 1)
         CALL set_pressure_basis_change(inputs%type_fe, inputs%type_fe - 1)
      END IF
      CALL deallocate_mesh(mesh_p1)



      !TESTTTTTTTTTTT
      IF (inputs%type_fe==3) THEN
         !         n_s(1, :) = (/2, 4, 5, 3/)
         !         n_s(2, :) = (/1, 6, 7, 3/)
         !         n_s(3, :) = (/1, 8, 9, 2/)
         !         DO m = 1, uu_mesh%me
         !            !WRITE(*,*) uu_mesh%jj(:,m)
         !            DO n = 1, 3
         !               !write(*,*) 'm = ', m, n, uu_mesh%jj(n_s(n,:),m)
         !               neigh = uu_mesh%neigh(n, m)
         !               IF (neigh.LE.0) CYCLE
         !               DO k = 1, 3
         !                  IF (uu_mesh%neigh(k, neigh) == m) EXIT
         !               END DO
         !               IF (ABS(uu_mesh%jce(n, m) - uu_mesh%jce(k, neigh)).NE.0) THEN
         !                  CALL error_petsc('ABS(uu_mesh%jce(n,m) - uu_mesh%jce(k,neigh)).NE.0')
         !               END IF
         !               i = SUM(ABS(uu_mesh%jj(n_s(n, :), m)) - uu_mesh%jj(n_s(k, :), neigh))
         !               IF (i.NE.0) THEN
         !                  CALL error_petsc('BUG')
         !               END IF
         !            END DO
         !         END DO
         !!$        WRITE(*,*) ' Coordinates'
         !!$        DO n = 1, uu_mesh%np
         !!$           WRITE(*,*) uu_mesh%rr(:,n)
         !!$        END DO
         !!$        WRITE(*,*) ' loc_to_glob'
         !!$        WRITE(*,*) uu_mesh%loc_to_glob
      END IF
      !TESTTTTTTTTTTTTTTTTT

      IF (rank==0) WRITE(*, *) 'P2/3 creation', (user_time() - t1)
      !IF (rank==0) WRITE(*, *) 'P2/3 creation domnp', uu_mesh%domnp
      !      do m = 1, mesh%me
      !         if(rank == 0 .and. minval(mesh%rr(1, mesh%jj(:, m))) == 0.d0) then
      !            write(*, *) 'r', mesh%rr(1, mesh%jj(:, m)), mesh%rr(2, mesh%jj(:, m))
      !         end if
      !      end do

      !===Create Sparsity pattern for the matrix (structure)
      t1 = user_time()
      CALL st_aij_csr_glob_block_with_extra_layer(PETSC_COMM_WORLD, 2, uu_mesh, LA_u)
      !TEST Scalar velocity matrix
      CALL st_aij_csr_glob_block_with_extra_layer(PETSC_COMM_WORLD, 1, uu_mesh, LA_u_scal)
      IF (rank==0) WRITE(*, *) ' st_aij_csr_glob_block', (user_time() - t1)

      dom_np = SIZE(LA_u%ia) - 1
      max_nnz_per_row = MAXVAL(LA_u%ia(1:dom_np) - LA_u%ia(0:dom_np - 1))

      !===Start Gauss points generation
      uu_mesh%edge_stab = .FALSE.
      t1 = user_time()
      CALL gauss_points_2d(uu_mesh, inputs%type_fe)
      IF (inputs%if_mesh_same) THEN
         CALL gauss_points_2d(pp_mesh, inputs%type_fe)
      ELSE
         CALL gauss_points_2d(pp_mesh, inputs%type_fe - 1)
      END IF
      IF (rank==0) WRITE(*, *) ' gauss_points_2d', (user_time() - t1)

      !===Dg approximation for pressure
      IF (inputs%if_pp_dg) THEN
         CALL create_dg_mesh_distributed(pp_mesh)
         !TODO isolated nodes at interfaces for dg
      END IF
      CALL st_aij_csr_glob_block_with_extra_layer(PETSC_COMM_WORLD, 1, pp_mesh, LA_p)

      !END IF
      !FIXME
      !ALLOCATE(mesh%iis(SIZE(mesh%jjs,1),mesh%mes))
      !CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
      !CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
      !mesh%nps = SIZE(mesh%j_s)
      !FIXME


      !===Boundary conditions on local domain
      CALL dirichlet_nodes_parallel(uu_mesh, inputs%Dir_list, js_d_loc_u)
      CALL dirichlet_nodes_parallel(pp_mesh, inputs%Dir_list, js_d_loc_p)

      !DEALLOCATE(js_D_loc_u)
      If (rank == 0) THEN
         ALLOCATE(js_D_loc_p(1))
         js_D_loc_p(1) = 1
      ELSE
         ALLOCATE(js_D_loc_p(0))
      End If

      OPEN(20, FILE = 'MESH_Velocity_Pressure', FORM = 'formatted')
      WRITE (20, *)  '# Cells  ', SUM(uu_mesh%domcell)
      WRITE (20, *)  '# Velocity grid points ', SUM(uu_mesh%domnp)
      WRITE (20, *)  '# Pressure grid points ', SUM(pp_mesh%domnp)
      CLOSE(20)

   END SUBROUTINE construct_distributed_mesh

   SUBROUTINE load_mesh_free_format_global_p1(dir, fil, list_dom, mesh, mesh_formatted)
      USE def_type_mesh
      USE chaine_caractere
      USE mod_gauss_points_2d
      USE Dir_nodes
      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
      TYPE(mesh_type) :: mesh_p1
      TYPE(mesh_type) :: mesh
      LOGICAL, INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
      INTEGER, DIMENSION(3) :: a_d
      INTEGER, DIMENSION(2) :: a_ds
      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els, &
           ancien_nd, ancien_el, ancien_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_el
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect
      LOGICAL :: test
      INTEGER :: mnouv, nnouv, i, nb_of_edges, n_a_d, nn
      INTEGER :: n, m, mop, edge, ms, msnouv, neighs1, neighs2
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh, nw_new, nws_new
      CHARACTER(len = 20) :: text
      CHARACTER(len = 2) :: truc
      INTEGER, PARAMETER :: type_fe = 1

      text = 'Mesh'
      DO n = 1, SIZE(list_dom)
         WRITE(truc, '(i2)') list_dom(n)
         text = TRIM(ADJUSTL(text)) // '_' // TRIM(ADJUSTL(truc))
      END DO

      d_end = last_c_leng (20, text)
      text = TRIM(ADJUSTL(text)) // '_FE_1'

      !WRITE (*, *) 'Loading mesh-file ...'
      IF (mesh_formatted) THEN
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'formatted')
      ELSE
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'unformatted')
      END IF
      OPEN(UNIT = 20, FILE = text, FORM = 'formatted', STATUS = 'unknown')

      !===READ GRID DATA AND ARRAY ALLOCATION
      !===Read P1 mesh
      IF (mesh_formatted) THEN
         READ  (30, *)  np, nw, me, nws, mes
      ELSE
         READ(30)  np, nw, me, nws, mes
      END IF

      IF (nw==3 .AND. nws==2) THEN !===Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed ', nw, nws
         STOP
      END IF

      ALLOCATE (jj_lect(nw, me), neigh_lect(nwneigh, me), i_d_lect(me))
      ALLOCATE (nouv_nd(np), ancien_nd(np), virgin_nd(np), &
           nouv_el(0:me), ancien_el(me), virgin_el(me))

      nouv_el = 0

      IF (mesh_formatted) THEN
         DO m = 1, me
            READ(30, *) jj_lect(:, m), neigh_lect(:, m), i_d_lect(m)
         END DO
      ELSE
         READ(30) jj_lect, neigh_lect, i_d_lect
      END IF

      !===Change enumeration
      virgin_nd = .TRUE.
      virgin_el = .TRUE.
      mnouv = 0
      nnouv = 0
      DO m = 1, me
         IF (MINVAL(ABS(list_dom - i_d_lect(m))) /= 0)  CYCLE !==cycle if cell m is not in the list_dom
         virgin_el(m) = .FALSE.
         mnouv = mnouv + 1  !==New element
         nouv_el(m) = mnouv
         ancien_el(mnouv) = m
         DO n = 1, nw; i = jj_lect(n, m)
         IF (virgin_nd(i)) THEN !==New point
            virgin_nd(i) = .FALSE.
            nnouv = nnouv + 1
            nouv_nd(i) = nnouv
            ancien_nd(nnouv) = i
         END IF
         END DO
      END DO
      mesh%me = mnouv
      mesh_p1%np = nnouv
      ALLOCATE(mesh_p1%jj(nw, mesh%me), mesh_p1%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))
      DO m = 1, mesh%me
         mesh_p1%jj(:, m) = nouv_nd(jj_lect(:, ancien_el(m)))
         mesh_p1%neigh(:, m) = nouv_el(neigh_lect(:, ancien_el(m)))
         mesh%i_d(m) = i_d_lect(ancien_el(m))
      END DO
      !===End change enumeration

      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))
      IF (mesh_formatted) THEN
         DO ms = 1, mes
            READ(30, *) jjs_lect(:, ms), neighs_lect(ms), sides_lect(ms)
         END DO
      ELSE
         READ(30) jjs_lect, neighs_lect, sides_lect
      END IF

      !===Change enumeration
      ALLOCATE(nouv_els(mes), ancien_els(mes))
      DEALLOCATE(virgin_el)
      ALLOCATE(virgin_el(mes))
      virgin_el = .TRUE.
      msnouv = 0
      DO ms = 1, mes
         neighs1 = neighs_lect(ms)
         DO n = 1, nw
            IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT
         END DO
         neighs2 = neigh_lect(n, neighs1)
         test = .FALSE.
         IF (MINVAL(ABS(list_dom - i_d_lect(neighs1))) == 0) THEN
            test = .TRUE.
         ELSE IF (neighs2 /= 0) THEN
            IF (MINVAL(ABS(list_dom - i_d_lect(neighs2))) == 0) THEN
               test = .TRUE.
               neighs_lect(ms) = neighs2 !===Change side
            END IF
         END IF
         IF (.NOT.test) CYCLE
         virgin_el(ms) = .FALSE.
         msnouv = msnouv + 1
         nouv_els(ms) = msnouv
         ancien_els(msnouv) = ms
      END DO
      mesh%mes = msnouv
      ALLOCATE (mesh_p1%jjs(nws, mesh%mes), mesh_p1%neighs(mesh%mes), &
           mesh%sides(mesh%mes))
      DO ms = 1, mesh%mes
         mesh_p1%jjs(:, ms) = nouv_nd(jjs_lect(:, ancien_els(ms)))
         mesh_p1%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
         mesh%sides(ms) = sides_lect(ancien_els(ms))
      END DO
      !===Check number of edges
      edge = 0
      DO m = 1, mesh%me
         DO n = 1, 3
            mop = mesh_p1%neigh(n, m)
            IF (mop==0) CYCLE !===Edge on boundary
            edge = edge + 1 !===New edge
         END DO
      END DO
      edge = edge / 2 !===Number of internal edges
      nb_of_edges = edge + mesh%mes
      IF (edge/=(3 * mesh%me - mesh%mes) / 2) THEN
         WRITE(*, *) ' BUG, edge/=(3*mesh%me + mesh%mes)/2'
         WRITE(*, *) ' edge ', edge, (3 * mesh%me - mesh%mes) / 2, MINVAL(mesh_p1%neigh), MAXVAL(mesh_p1%neigh)
         WRITE(*, *) ' mesh%mes ', mesh%mes, ' mesh%me ', mesh%me
         STOP
      END IF
      !===End check number of edges
      !===Change enumeration

      ALLOCATE(rr_lect(kd, np))
      IF (mesh_formatted) THEN
         DO n = 1, np
            READ(30, *) rr_lect(:, n)
         END DO
      ELSE
         READ(30) rr_lect
      END IF
      ALLOCATE(mesh_p1%rr(kd, mesh_p1%np))
      mesh_p1%rr = rr_lect(:, ancien_nd(1:mesh_p1%np))

      !===Make sure indexing in done with the lowest to highest index convention
      jj_lect(:, :mesh%me) = mesh_p1%jj
      jjs_lect(:, :mesh%mes) = mesh_p1%jjs
      neigh_lect(:, :mesh%me) = mesh_p1%neigh
      DO m = 1, mesh%me !===loop on the elements
         CALL tri_jlg(jj_lect(:, m), a_d, n_a_d)
         DO n = 1, nw
            i = mesh_p1%jj(n, m)
            DO nn = 1, nw
               IF (a_d(nn) == i) THEN
                  mesh_p1%neigh(nn, m) = neigh_lect(n, m)
                  EXIT
               END IF
            END DO
         END DO
         mesh_p1%jj(:, m) = a_d
      END DO
      DO ms = 1, mes !===loop on the elements
         CALL tri_jlg(jjs_lect(:, ms), a_ds, n_a_d)
         mesh_p1%jjs(:, ms) = a_ds
      END DO
      !===End make sure indexing in done with the lowest to highest index convention

      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_el, virgin_nd)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els, ancien_nd, ancien_el, ancien_els)
      !===End of P1 grid reading

      !===Prepare actual mesh (works in 2D only)
      IF (kd==3) THEN
         WRITE(*, *) 'k_d==3 not programmed yet'
         STOP
      END IF
      mesh%np = mesh_p1%np
      mesh%dom_np = mesh%np
      nw_new = (type_fe + 1) * (type_fe + 2) / 2 !===3
      nws_new = type_fe + 1 !===2
      ALLOCATE(mesh%jj(nw_new, mesh%me))
      ALLOCATE(mesh%neigh(nw, mesh%me))
      ALLOCATE(mesh%jjs(nws_new, mesh%mes))
      ALLOCATE(mesh%neighs(mesh%mes))
      ALLOCATE(mesh%rr(kd, mesh%np))

      ALLOCATE(mesh%disp(2), mesh%domnp(1))
      mesh%disp = (/1, mesh%np + 1/)
      mesh%domnp = (/mesh%np/)
      ALLOCATE(mesh%discell(2), mesh%domcell(1))
      mesh%discell = (/1, mesh%me + 1/)
      mesh%domcell = (/mesh%me/)

      ALLOCATE(mesh%loc_to_glob(mesh%np))
      DO n = 1, mesh%np
         mesh%loc_to_glob(n) = n
      END DO
      !===Create final mesh
      mesh%neigh = mesh_p1%neigh
      mesh%neighs = mesh_p1%neighs

      mesh%jj = mesh_p1%jj
      mesh%jjs = mesh_p1%jjs
      mesh%rr = mesh_p1%rr
      mesh%mextra = 0
      ALLOCATE(mesh%extra_jj(3, 0))
      ALLOCATE(mesh%extra_jce(3, 0))
      ALLOCATE(mesh%extra_jcc(0))

      DEALLOCATE(mesh_p1%jj, mesh_p1%jjs, mesh_p1%rr, mesh_p1%neigh, mesh_p1%neighs)
      !===End Prepare actual mesh

      mesh%nis = 0
      ALLOCATE(mesh%isolated_jjs(0), mesh%isolated_interfaces(0, 2))

      !===Prepare face structures (jce, jev)
      CALL prep_jce_jev(mesh)
      ALLOCATE(mesh%disedge(2), mesh%domedge(1))
      mesh%disedge = (/1, mesh%medge + 1/)
      mesh%domedge = (/mesh%medge/)
      mesh%medges = 0

      ALLOCATE(mesh%jees(0), mesh%jecs(0))
      WRITE (20, *)  'np ', mesh%np, 'me ', mesh%me, 'mes ', mesh%mes, 'nps ', mesh%nps

      mesh%edge_stab = .FALSE.
      mesh%mi = 0

      CLOSE(20)
      CLOSE(30)
   END SUBROUTINE load_mesh_free_format_global_p1

   SUBROUTINE neighbours(mesh, m, ni, neighs, nb_neighs)
      USE input_data
      IMPLICIT NONE
      TYPE(mesh_type), INTENT(IN) :: mesh
      INTEGER :: m, ni, nc, nn, np, k, k_start, nb_n, n_start, ki, kn, nb_neighs
      INTEGER, DIMENSION(2) :: n_ks, n_starts
      INTEGER, DIMENSION(:), ALLOCATABLE :: neighs
      LOGICAL :: side, break

      IF (MINVAL(ABS(mesh%neigh(:, m))) == 0) THEN
         nb_neighs = 1
         DEALLOCATE(neighs)
         ALLOCATE(neighs(nb_neighs))
         neighs(1) = m
         return
      END IF

      IF (ni <= 3) THEN
         n_ks = (/MODULO(ni, 3) + 1, MODULO(ni + 1, 3) + 1/)
         n_starts = mesh%neigh(n_ks, m)

         IF (MINVAL(ABS(n_starts)) == 0) THEN
            nb_neighs = 1
            DEALLOCATE(neighs)
            ALLOCATE(neighs(nb_neighs))
            neighs(1) = m
            return
         END IF

         nc = mesh%neigh(n_ks(1), m)
         IF (nc  == 0) THEN
            side = .true.
            nc = mesh%neigh(n_ks(2), m)
         ELSE
            side = .false.
         END IF

         np = m
         nb_neighs = 2
         break = .false.

         DO WHILE(nc /= m)

            DO ki = 1, 3
               IF (mesh%jj(ki, nc) == mesh%jj(ni, m)) EXIT
            END DO

            DO kn = 1, 3
               IF (mesh%neigh(kn, nc) == np) EXIT
            END DO

            IF (ki + kn == 3) THEN
               nn = mesh%neigh(3, nc)
            ELSE IF (ki + kn == 5) THEN
               nn = mesh%neigh(1, nc)
            ELSE IF (ki + kn == 4) THEN
               nn = mesh%neigh(2, nc)
            END IF

            IF (nn == 0) THEN
               IF (side) THEN
                  EXIT
               ELSE
                  side = .True.
                  nc = mesh%neigh(n_ks(2), m)
                  IF (nc  == 0) EXIT
                  np = m
                  nb_neighs = nb_neighs + 1
               END IF
            ELSE IF(nn == m) THEN
               EXIT
            ELSE
               nb_neighs = nb_neighs + 1
               np = nc
               nc = nn
            END IF

         END DO
         DEALLOCATE(neighs)
         ALLOCATE(neighs(nb_neighs))
         neighs(1) = m
         nc = mesh%neigh(n_ks(1), m)
         IF (nc  == 0) THEN
            side = .true.
            nc = mesh%neigh(n_ks(2), m)
         ELSE
            side = .false.
         END IF

         np = m
         nb_n = 2
         neighs(2) = nc
         break = .false.

         DO WHILE(nc /= m)

            DO ki = 1, 3
               IF (mesh%jj(ki, nc) == mesh%jj(ni, m)) EXIT
            END DO

            DO kn = 1, 3
               IF (mesh%neigh(kn, nc) == np) EXIT
            END DO

            IF (ki + kn == 3) THEN
               nn = mesh%neigh(3, nc)
            ELSE IF (ki + kn == 5) THEN
               nn = mesh%neigh(1, nc)
            ELSE IF (ki + kn == 4) THEN
               nn = mesh%neigh(2, nc)
            END IF

            IF (nn == 0) THEN
               IF (side) THEN
                  EXIT
               ELSE
                  side = .True.
                  nc = mesh%neigh(n_ks(2), m)
                  IF (nc  == 0) EXIT
                  np = m
                  nb_n = nb_n + 1
                  neighs(nb_n) = nc
               END IF
            ELSE IF(nn == m) THEN
               EXIT
            ELSE
               nb_n = nb_n + 1
               np = nc
               nc = nn
               neighs(nb_n) = nc
            END IF

         END DO
      ELSE IF (ni<=9) THEN
         nb_neighs = 2
         DEALLOCATE(neighs)
         ALLOCATE(neighs(nb_neighs))
         neighs(1) = m
         neighs(2) = mesh%neigh((ni - 4) / (inputs%type_fe - 1) + 1, m)
      ELSE
         nb_neighs = 1
         DEALLOCATE(neighs)
         ALLOCATE(neighs(nb_neighs))
         neighs(1) = m
      END IF

   END SUBROUTINE neighbours

END MODULE construct_meshes
