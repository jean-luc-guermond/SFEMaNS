!
!Authors Jean-Luc Guermond, Luigi Quarapelle, Copyrights 1994, 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE prep_maill

   IMPLICIT NONE

   PUBLIC :: load_mesh, load_mesh_formatted, load_mesh_free_format, &
        load_dg_mesh_free_format, load_mesh_free_format_ordered, prep_interfaces, &
        create_p3_mesh, incr_vrtx_indx_enumeration, incr_vrtx_indx_enumeration_for_interfaces, &
        create_iso_grid_distributed, prep_jce_jev, refinement_iso_grid_distributed
   PRIVATE

CONTAINS

   !------------------------------------------------------------------------------
   SUBROUTINE create_p3_mesh(p1_mesh, p2_mesh, p3_mesh, type_fe)
      USE def_type_mesh
      USE chaine_caractere
      USE Dir_nodes
      USE my_util
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type), INTENT(IN) :: p1_mesh
      TYPE(mesh_type), INTENT(IN) :: p2_mesh
      TYPE(mesh_type), INTENT(OUT) :: p3_mesh
      LOGICAL, DIMENSION(p2_mesh%np) :: p2_virgin
      INTEGER :: v_dof, f_dof, c_dof, n_dof, edge, i, h, k, l, m, ms, n, n1, n2, j1, j2, &
           nb_edges, nb_vertices, lmax, sgn, sgn_op, j1_op, j2_op, n1_op, n2_op, n_op, m_op, hh, hh_op, &
           n_mid, j_mid
      REAL(KIND = 8) :: s1, s2, s3, sh, sk
      REAL(KIND = 8), DIMENSION(2) :: r1, r2, r3
      !===Programmed for 2d Lagrage elements (type_fe can be arbitrary)
      v_dof = 3 !===Vertex dofs
      f_dof = type_fe - 1 !===Face dofs
      c_dof = (type_fe - 2) * (type_fe - 1) / 2 !===Cell dofs
      p3_mesh%gauss%n_w = (type_fe + 1) * (type_fe + 2) / 2
      p3_mesh%gauss%n_ws = type_fe + 1
      p3_mesh%gauss%k_d = 2

      p3_mesh%me = p1_mesh%me
      p3_mesh%mes = p1_mesh%mes
      !===Nb of edges
      edge = 0
      DO m = 1, p1_mesh%me
         DO n = 1, 3
            IF (p1_mesh%neigh(n, m)==0) CYCLE
            edge = edge + 1
         END DO
      END DO
      edge = edge / 2
      nb_edges = edge + p1_mesh%mes
      IF (edge /= (3 * p1_mesh%me - p1_mesh%mes) / 2) THEN
         CALL error_PETSC('BUG in create_p3_mesh, nb of edges is wrong')
      END IF
      !===Nb of P3 nodes
      nb_vertices = p2_mesh%np - nb_edges
      p3_mesh%np = nb_vertices + 2 * nb_edges + p1_mesh%me

      IF (nb_vertices.NE.p1_mesh%np) THEN
         CALL error_petsc(' BUG in create_p3_mesh, nb_vertices.NE.p1_mesh%np')
      END IF

      !===Allocate mesh structure
      ALLOCATE(p3_mesh%jj(p3_mesh%gauss%n_w, p3_mesh%me), p3_mesh%neigh(3, p3_mesh%me), p3_mesh%i_d(p3_mesh%me))
      ALLOCATE(p3_mesh%jjs(p3_mesh%gauss%n_ws, p3_mesh%mes), p3_mesh%neighs(p3_mesh%mes), p3_mesh%sides(p3_mesh%mes))
      ALLOCATE(p3_mesh%rr(p3_mesh%gauss%k_d, p3_mesh%np))

      !===Easy stuff
      p3_mesh%neigh = p1_mesh%neigh
      p3_mesh%neighs = p1_mesh%neighs
      p3_mesh%i_d = p1_mesh%i_d
      p3_mesh%sides = p1_mesh%sides

      !===Create vertices
      p3_mesh%jj(1:3, :) = p1_mesh%jj(1:3, :)
      IF (MAXVAL(p3_mesh%jj(1:3, :)).NE.nb_vertices) THEN
         write(*, *) MAXVAL(p1_mesh%jj(1:3, p1_mesh%me)), nb_vertices
         write(*, *) MAXVAL(p3_mesh%jj(1:3, p1_mesh%me)), nb_vertices
         CALL error_PETSC('BUG in create_p3_mesh, Vertex enumeration not done properly')
      END IF
      p3_mesh%rr(:, 1:nb_vertices) = p1_mesh%rr(:, 1:nb_vertices)

      !===Create face dofs
      p2_virgin = .TRUE.
      n_dof = nb_vertices
      DO m = 1, p3_mesh%me
         DO n = 1, v_dof !===n is the face number of the face to be populated
            n_mid = n + 3 !===n_mid is the mid point on the p2 mesh
            j_mid = p2_mesh%jj(n_mid, m)
            IF (.NOT.p2_virgin(j_mid)) THEN !===This edge has already been visited
               m_op = p1_mesh%neigh(n, m)
               DO n_op = 1, v_dof
                  IF (m == p1_mesh%neigh(n_op, m_op)) THEN
                     EXIT !===n_op is the index of the face I am looking for
                  END IF
               END DO
               n1 = MODULO(n, 3) + 1
               n2 = MODULO(n + 1, 3) + 1
               n1_op = MODULO(n_op, 3) + 1
               n2_op = MODULO(n_op + 1, 3) + 1
               j1 = p1_mesh%jj(n1, m)
               j1_op = p1_mesh%jj(n1_op, m_op)
               j2 = p1_mesh%jj(n2, m)
               j2_op = p1_mesh%jj(n2_op, m_op)
               IF (j1.NE.j1_op) THEN
                  i = n2
                  n2 = n1
                  n1 = i
               END IF
               IF (n2 - n1>0) THEN
                  sgn = 1
               ELSE
                  sgn = -1
               END IF
               IF (n2_op - n1_op>0) THEN
                  sgn_op = 1
               ELSE
                  sgn_op = -1
               END IF
               DO h = 1, f_dof
                  hh = h * sgn + ((1 - sgn) / 2) * (f_dof + 1)
                  hh_op = h * sgn_op + ((1 - sgn_op) / 2) * (f_dof + 1)
                  p3_mesh%jj(v_dof + (n - 1) * f_dof + hh, m) = p3_mesh%jj(v_dof + (n_op - 1) * f_dof + hh_op, m_op)
               END DO
            ELSE
               p2_virgin(j_mid) = .FALSE.
               n1 = MODULO(n, 3) + 1
               n2 = MODULO(n + 1, 3) + 1
               j1 = p2_mesh%jj(n1, m)
               j2 = p2_mesh%jj(n2, m)
               r1 = p2_mesh%rr(:, j1)
               r2 = p2_mesh%rr(:, j2)
               r3 = p2_mesh%rr(:, j_mid) !===middle node is P2 to have curved boundaries

               !===Create face dofs
               DO h = 1, f_dof
                  n_dof = n_dof + 1 !===New index created
                  p3_mesh%jj(v_dof + (n - 1) * f_dof + h, m) = n_dof
                  s1 = 0.d0
                  s2 = 1.d0
                  s3 = 0.5d0
                  IF (n1<n2) THEN
                     sh = h / DBLE(type_fe) !===Move from r1 to r2
                  ELSE
                     sh = 1 - h / DBLE(type_fe) !===Move from r2 to r1
                  END IF
                  p3_mesh%rr(:, n_dof) = r1(:) * (sh - s2) * (sh - s3) / ((s1 - s2) * (s1 - s3)) &
                       + r2(:) * (sh - s3) * (sh - s1) / ((s2 - s3) * (s2 - s1)) &
                       + r3(:) * (sh - s1) * (sh - s2) / ((s3 - s1) * (s3 - s2))
               END DO
            END IF
         END DO

         !===Create volume dofs
         r1 = p1_mesh%rr(:, p1_mesh%jj(1, m))
         r2 = p1_mesh%rr(:, p1_mesh%jj(2, m))
         r3 = p1_mesh%rr(:, p1_mesh%jj(3, m))
         l = v_dof + 3 * f_dof !===Initialize cell dofs
         lmax = type_fe - 2
         DO h = 1, lmax
            sh = h / DBLE(type_fe) !===move along y-axis
            DO k = 1, lmax - h + 1
               sk = k / DBLE(type_fe) !===move along x-axis
               n_dof = n_dof + 1
               p3_mesh%jj(p3_mesh%gauss%n_w, m) = n_dof !===Volume dof
               p3_mesh%rr(:, n_dof) = (1 - sh - sk) * r1 + sk * r2 + sh * r3
            END DO
         END DO
      END DO

      !===Create surface connectivity
      DO ms = 1, p1_mesh%mes
         m_op = p1_mesh%neighs(ms)
         j1 = p1_mesh%jjs(1, ms)
         j2 = p1_mesh%jjs(2, ms)
         DO n_op = 1, 3
            IF(j1==p1_mesh%jj(n_op, m_op) .OR. j2==p1_mesh%jj(n_op, m_op)) CYCLE
            !===n is the index of the face that is on the boundary
            p3_mesh%jjs(1:2, ms) = p1_mesh%jjs(1:2, ms) !===Vertex dofs
            n1_op = MODULO(n_op + 1, 3) + 1
            n2_op = MODULO(n_op + 2, 3) + 1
            j1_op = p1_mesh%jj(n1_op, m_op)
            IF (j1.NE.j1_op) THEN
               i = n1_op
               n1_op = n2_op
               n2_op = i
            END IF
            IF (n2_op>n1_op) THEN
               sgn_op = 1
            ELSE
               sgn_op = -1
            END IF
            DO h = 1, f_dof
               hh_op = h * sgn_op + ((1 - sgn_op) / 2) * (f_dof + 1)
               p3_mesh%jjs(2 + h, ms) = p3_mesh%jj(v_dof + (n_op - 1) * f_dof + hh_op, m_op) !===Face dofs
            END DO
         END DO
      END DO

      !===Take care of leftover
      ALLOCATE(p3_mesh%iis(p3_mesh%gauss%n_ws, p3_mesh%mes))
      CALL dirichlet_nodes(p3_mesh%jjs, SPREAD(1, 1, p3_mesh%mes), SPREAD(.TRUE., 1, 1), p3_mesh%j_s)
      CALL surf_nodes_i(p3_mesh%jjs, p3_mesh%j_s, p3_mesh%iis)
      p3_mesh%nps = SIZE(p3_mesh%j_s)

   END SUBROUTINE create_p3_mesh
   !===

   SUBROUTINE incr_vrtx_indx_enumeration(mesh, type_fe)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      INTEGER, INTENT(IN) :: type_fe
      INTEGER :: nmin(1), nmax(1), jjive(3)
      INTEGER :: m, ms, n
      !===JLG July 20, 2019, p3 mesh
      IF (type_fe<0) RETURN
      !===Use increasing vertex index enumeration
      DO m = 1, mesh%me
         nmin = minloc(mesh%jj(1:3, m))
         nmax = maxloc(mesh%jj(1:3, m))
         n = 6 / (nmin(1) * nmax(1))
         jjive(1) = mesh%jj(nmin(1), m)
         jjive(2) = mesh%jj(n, m)
         jjive(3) = mesh%jj(nmax(1), m)
         mesh%jj(1:3, m) = jjive
         !===Do it also for p2 (This is awkward. It should have been done well before.)
         IF (type_fe==2) THEN
            jjive(1) = mesh%jj(nmin(1) + 3, m)
            jjive(2) = mesh%jj(n + 3, m)
            jjive(3) = mesh%jj(nmax(1) + 3, m)
            mesh%jj(4:6, m) = jjive
         END IF
         !===Do it also for mesh%neigh
         jjive(1) = mesh%neigh(nmin(1), m)
         jjive(2) = mesh%neigh(n, m)
         jjive(3) = mesh%neigh(nmax(1), m)
         mesh%neigh(:, m) = jjive
      END DO
      DO ms = 1, mesh%mes
         nmin = minloc(mesh%jjs(1:2, ms))
         nmax = maxloc(mesh%jjs(1:2, ms))
         jjive(1) = mesh%jjs(nmin(1), ms)
         jjive(2) = mesh%jjs(nmax(1), ms)
         mesh%jjs(1:2, ms) = jjive(1:2)
         !===For p2 only (nws=3)
      END DO
      !===JLG July 20, 2019, p3 mesh
   END SUBROUTINE incr_vrtx_indx_enumeration

   SUBROUTINE incr_vrtx_indx_enumeration_for_interfaces(interface_XX, nws1, nws2)
      USE def_type_mesh
      INTEGER, INTENT(IN) :: nws1, nws2
      TYPE(interface_type) :: interface_XX
      INTEGER :: js1ive(nws1), js2ive(nws2)
      INTEGER :: ms
      DO ms = 1, interface_XX%mes
         IF (interface_XX%jjs1(2, ms)<interface_XX%jjs1(1, ms)) THEN
            js1ive = interface_XX%jjs1(:, ms)
            interface_XX%jjs1(2, ms) = js1ive(1)
            interface_XX%jjs1(1, ms) = js1ive(2)
            IF(nws1.GE.3) interface_XX%jjs1(nws1:3:-1, ms) = js1ive(3:nws1)
         END IF
         IF (interface_XX%jjs2(2, ms)<interface_XX%jjs2(1, ms)) THEN
            js2ive = interface_XX%jjs2(:, ms)
            interface_XX%jjs2(2, ms) = js2ive(1)
            interface_XX%jjs2(1, ms) = js2ive(2)
            IF(nws2.GE.3) interface_XX%jjs2(nws2:3:-1, ms) = js2ive(3:nws2)
         END IF
      END DO
   END SUBROUTINE incr_vrtx_indx_enumeration_for_interfaces

   SUBROUTINE load_dg_mesh_free_format(dir, fil, list_dom, list_inter, type_fe, &
        mesh, mesh_formatted)
      USE def_type_mesh
      USE chaine_caractere
      USE Dir_nodes
      USE my_util
      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom, list_inter
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type) :: mesh
      LOGICAL, INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect, stat
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_ms
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect
      LOGICAL :: t1, t2
      INTEGER :: mnouv, nnouv, i, dom
      INTEGER :: n, m, mop, ms, msnouv, neighs1, neighs2, mes_int
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh
      CHARACTER(len = 40) :: text
      CHARACTER(len = 2) :: truc

      mesh%gauss%k_d = 2 !===Space dimension (meridian section)

      text = 'Mesh'
      d_end = last_c_leng (20, text)
      DO n = 1, SIZE(list_dom)
         d_end = last_c_leng (20, text)
         WRITE(truc, '(i2)') list_dom(n)
         f_end = start_of_string (truc)
         text = text(1:d_end) // '_' // truc(f_end:)
      END DO

      d_end = last_c_leng (20, text)
      text = TRIM(ADJUSTL(text)) // '_FE_1'

      !WRITE (*,*) 'Loading mesh-file ...'
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

      IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==6 .AND. nws==3) THEN
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE IF (nw==10 .AND. nws==6) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed ', nw, nws
         STOP
      END IF

      ALLOCATE (jj_lect(nw, me), neigh_lect(nwneigh, me), i_d_lect(me))
      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))
      ALLOCATE(rr_lect(kd, np))

      IF (mesh_formatted) THEN
         DO m = 1, me
            READ(30, *) jj_lect(:, m), neigh_lect(:, m), i_d_lect(m)
         END DO
         DO ms = 1, mes
            READ(30, *) jjs_lect(:, ms), neighs_lect(ms), sides_lect(ms)
         END DO
         DO n = 1, np
            READ(30, *) rr_lect(:, n)
         END DO
      ELSE
         READ(30) jj_lect, neigh_lect, i_d_lect
         READ(30) jjs_lect, neighs_lect, sides_lect
         READ(30) rr_lect
      END IF

      !===Re-order at cell level to get proper edges orientation
      DO m = 1, me
         IF (jj_lect(1, m) >  jj_lect(2, m)) THEN
            CALL switch_cell_vertices(jj_lect, neigh_lect, m, 1, 2)
         END IF
         IF (jj_lect(2, m) >  jj_lect(3, m)) THEN
            CALL switch_cell_vertices(jj_lect, neigh_lect, m, 2, 3)
         END IF
         IF (jj_lect(1, m) >  jj_lect(2, m)) THEN
            CALL switch_cell_vertices(jj_lect, neigh_lect, m, 1, 2)
         END IF
      END DO

      DO ms = 1, mes
         IF (jjs_lect(1, ms) >  jjs_lect(2, ms)) THEN
            m = jjs_lect(1, ms)
            jjs_lect(1, ms) = jjs_lect(2, ms)
            jjs_lect(2, ms) = m
         END IF
      END DO

      ! Identify the status of faces
      ! stat = 1 (interface to be forgotten), stat = 2 (boundary), stat = 3 (real interface), stat = 4 (internal interface)
      ALLOCATE (stat(mes))
      DO ms = 1, mes
         neighs1 = neighs_lect(ms)
         IF (neighs1==0) THEN
            WRITE(*, *) ' BUG in prep_mesh, neighs1=0 '
            STOP
         END IF
         IF (MINVAL(ABS(i_d_lect(neighs1) - list_dom))==0) THEN
            t1 = .TRUE.
         ELSE
            t1 = .FALSE.
         END IF
         DO n = 1, nw
            IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT ! n not on the interface
         END DO
         neighs2 = neigh_lect(n, neighs1)
         IF (neighs2==0) THEN
            IF (t1) THEN
               stat(ms) = 2  ! face on the boundary of the domain of interest
            ELSE
               stat(ms) = 1  ! face does not touch the domain of interest
            END IF
            CYCLE
         END IF
         ! neighs2 /=0
         IF (MINVAL(ABS(i_d_lect(neighs2) - list_dom))==0) THEN
            t2 = .TRUE.
         ELSE
            t2 = .FALSE.
         END IF

         IF (t1) THEN
            IF (t2) THEN
               IF (SIZE(list_inter)==0) THEN
                  stat(ms) = 4 ! internal interface
               ELSE IF (MINVAL(ABS(sides_lect(ms) - list_inter))==0) THEN
                  stat(ms) = 3 ! real interface
               ELSE
                  stat(ms) = 4 ! internal interface
               END IF
            ELSE
               stat(ms) = 2 ! face at the boundary of the domain of interest
            END IF
         ELSE
            IF (t2) THEN
               stat(ms) = 2 ! face at the boundary of the domain of interest
            ELSE
               stat(ms) = 1 ! on an interface of no interest
            END IF
         END IF

      END DO

      !===Count number of cell/cells/borders within list_dom
      ALLOCATE (nouv_nd(np), nouv_els(mes), virgin_nd(np), virgin_ms(mes), nouv_el(me))
      nouv_nd = -1000
      virgin_nd = .TRUE.
      virgin_ms = .TRUE.
      mnouv = 0
      msnouv = 0
      nnouv = 0
      mes_int = 0
      DO dom = 1, SIZE(list_dom)
         DO m = 1, me ! Count new nodes from domain: i_d=dom
            IF (list_dom(dom) /= i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
            mnouv = mnouv + 1  ! Nouvel element
            nouv_el(m) = mnouv
            DO n = 1, nw; i = jj_lect(n, m)
            IF (virgin_nd(i)) THEN ! Nouveau point
               virgin_nd(i) = .FALSE.
               nnouv = nnouv + 1
            END IF
            END DO
         END DO

         DO ms = 1, mes
            IF (stat(ms) /= 1 .AND. stat(ms) /=2 .AND. stat(ms) /=3 .AND. stat(ms) /=4) THEN
               WRITE(*, *) ' BUG in prep_mesh, stat out of bounds '
               STOP
            END IF

            !I test if ms touches the current domain of interest: i_d = dom
            IF (stat(ms) == 1) CYCLE
            IF (stat(ms) == 4)  THEN ! if internal interface
               IF (virgin_ms(ms)) THEN
                  mes_int = mes_int + 1
                  virgin_ms(ms) = .FALSE.
               END IF
            ELSE
               neighs1 = neighs_lect(ms)
               DO n = 1, nw
                  IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT
                  ! exit when n is not on the interface
               END DO
               neighs2 = neigh_lect(n, neighs1)
               IF (i_d_lect(neighs1) /= list_dom(dom)) THEN
                  IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
                  IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
               END IF
               !End test if ms touches the domain of interest

               IF (virgin_ms(ms)) THEN !New interface
                  virgin_ms(ms) = .FALSE.
                  msnouv = msnouv + 1
               END IF
               IF (stat(ms) ==3) THEN
                  ! Nodes and sides on the interface are virgin again
                  virgin_nd(jjs_lect(:, ms)) = .TRUE. ! interface nodes are virgin again
                  virgin_ms(ms) = .TRUE.
               END IF
            END IF
         END DO
      END DO
      mesh%me = mnouv
      mesh%np = nnouv
      mesh%mes = msnouv
      mesh%mes_int = mes_int

      !===Re-ordering
      ALLOCATE(mesh%jj(nw, mesh%me), mesh%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))
      ALLOCATE(mesh%jjs(nws, mesh%mes), mesh%neighs(mesh%mes), mesh%sides(mesh%mes))
      ALLOCATE(mesh%jjs_int(nws, mesh%mes_int), mesh%neighs_int(2, mesh%mes_int), mesh%sides_int(mesh%mes_int))
      ALLOCATE(mesh%rr(kd, mesh%np))

      virgin_nd = .TRUE.
      virgin_ms = .TRUE.
      nnouv = 0
      msnouv = 0
      DO dom = 1, SIZE(list_dom)
         !Loop on me and get nouv_el and nouv_nd right
         DO m = 1, me
            IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
            DO n = 1, nw; i = jj_lect(n, m)
            IF (virgin_nd(i)) THEN ! Nouveau point
               virgin_nd(i) = .FALSE.
               nnouv = nnouv + 1
               nouv_nd(i) = nnouv
            END IF
            END DO
         END DO

         !Loop again on me and update
         DO m = 1, me
            IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
            DO n = 1, nw; i = jj_lect(n, m)
            IF (n .LE. nwneigh) THEN
               mop = neigh_lect(n, m)
               IF (mop .LE. 0) THEN
                  mesh%neigh(n, nouv_el(m)) = 0
               ELSE IF (MINVAL(ABS(list_dom - i_d_lect(mop))) == 0) THEN
                  mesh%neigh(n, nouv_el(m)) = nouv_el(mop)
               ELSE
                  mesh%neigh(n, nouv_el(m)) = 0
               END IF
            END IF
            mesh%rr(:, nouv_nd(i)) = rr_lect(:, i)
            END DO
            mesh%i_d(nouv_el(m)) = i_d_lect(m)
            mesh%jj(:, nouv_el(m)) = nouv_nd(jj_lect(:, m))
         END DO

         !Loop on mes and get neighs_lect and nouv_els right
         DO ms = 1, mes
            !I test if ms touches the current domain of interest: i_d = dom
            IF (stat(ms) == 1) CYCLE
            IF (stat(ms) == 4) CYCLE
            neighs1 = neighs_lect(ms)
            DO n = 1, nw
               IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT
               ! exit when n is not on the interface
            END DO
            neighs2 = neigh_lect(n, neighs1)
            IF (i_d_lect(neighs1) /= list_dom(dom)) THEN
               IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
               IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
            END IF
            !End test if ms touches the domain of interest

            IF (virgin_ms(ms)) THEN !New interface
               virgin_ms(ms) = .FALSE.
               msnouv = msnouv + 1
               nouv_els(ms) = msnouv
            END IF
            IF (stat(ms) ==3) THEN
               ! Nodes and sides on the interface are virgin again
               virgin_nd(jjs_lect(:, ms)) = .TRUE. ! interface nodes are virgin again
               virgin_ms(ms) = .TRUE.
            END IF

            ! Swapping problem
            neighs1 = neighs_lect(ms)
            IF ((ABS(i_d_lect(neighs1) - list_dom(dom)))==0) THEN
               t1 = .TRUE.
            ELSE
               t1 = .FALSE.
            END IF
            DO n = 1, nw
               IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT ! n not on the interface
            END DO
            neighs2 = neigh_lect(n, neighs1)
            IF (neighs2==0) THEN
               CYCLE
            END IF
            ! neighs2 /=0
            IF ((ABS(i_d_lect(neighs2) - list_dom(dom)))==0) THEN
               t2 = .TRUE.
            ELSE
               t2 = .FALSE.
            END IF
            IF (.NOT.t1 .AND. t2) THEN
               neighs_lect(ms) = neighs2 !get things right (swap neighs)
            END IF
         END DO

         mes_int = 0
         mesh%neighs_int = -1
         !Loop again on mes and update
         DO ms = 1, mes

            !I test if ms touches the current domain of interest: i_d = dom
            IF (stat(ms) == 1) CYCLE
            neighs1 = neighs_lect(ms)
            DO n = 1, nw
               IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT
               ! exit when n is not on the interface
            END DO
            neighs2 = neigh_lect(n, neighs1)
            IF (i_d_lect(neighs1) /= list_dom(dom)) THEN
               IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
               IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
            END IF
            !End test if ms touches the domain of interest
            IF (stat(ms) == 4) THEN
               mes_int = mes_int + 1
               mesh%jjs_int(:, mes_int) = nouv_nd(jjs_lect(:, ms))
               mesh%neighs_int(1, mes_int) = nouv_el(neighs_lect(ms))
               mesh%sides_int(mes_int) = sides_lect(ms)
            ELSE
               mesh%jjs(:, nouv_els(ms)) = nouv_nd(jjs_lect(:, ms))
               mesh%neighs(nouv_els(ms)) = nouv_el(neighs_lect(ms))
               mesh%sides(nouv_els(ms)) = sides_lect(ms)
               IF (stat(ms)==3) THEN ! side is an interface to be kept
                  neighs1 = neighs_lect(ms)
                  DO n = 1, nw
                     IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT ! n not on the interface
                  END DO
                  mesh%neigh(n, nouv_el(neighs1)) = 0
                  neighs2 = neigh_lect(n, neighs1)
                  DO n = 1, nw
                     IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs2))) /= 0) EXIT ! n not on the interface
                  END DO
                  mesh%neigh(n, nouv_el(neighs2)) = 0
               END IF
            END IF

         END DO

      END DO
      !===End reordring
      mesh%gauss%n_ws = SIZE(mesh%jjs, 1)
      mesh%gauss%n_w = SIZE(mesh%jj, 1)
      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_nd, virgin_ms, stat)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els)
      write(*,*) 'ok0', mesh%neighs_int(1, mes_int), mesh%neighs_int(2, mes_int), mes_int
      !===Prepare actual mesh (works in 2D only)
      IF (kd==3) THEN
         WRITE(*, *) 'k_d==3 not programmed yet'
         STOP
      END IF
      mesh%dom_np = mesh%np

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

      mesh%mextra = 0
      ALLOCATE(mesh%jj_extra(3, 0))
      ALLOCATE(mesh%jce_extra(3, 0))
      ALLOCATE(mesh%jcc_extra(0))

      mesh%mes_extra = 0
      ALLOCATE(mesh%neighs_extra(0))
      ALLOCATE(mesh%sides_extra(0))
      ALLOCATE(mesh%jjs_extra(2, 0))
      ALLOCATE(mesh%rrs_extra(2, 3, 0))

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

   END SUBROUTINE load_dg_mesh_free_format

   !------------------------------------------------------------------------------

   SUBROUTINE load_mesh_free_format_ordered(dir, fil, list_dom, type_fe, &
        mesh, mesh_formatted, edge_stab)

      USE def_type_mesh
      USE chaine_caractere
      USE Dir_nodes

      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type) :: mesh
      LOGICAL, INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
      LOGICAL, OPTIONAL, INTENT(IN) :: edge_stab

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els, &
           ancien_nd, ancien_el, ancien_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_el
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect

      LOGICAL :: test
      INTEGER :: mnouv, nnouv, i, dom
      INTEGER :: n, m, ms, msnouv, neighs1, neighs2
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh
      CHARACTER(len = 60) :: text
      CHARACTER(len = 2) :: truc
      text = 'Mesh'
      d_end = last_c_leng (20, text)
      DO n = 1, SIZE(list_dom)
         d_end = last_c_leng (20, text)
         WRITE(truc, '(i2)') list_dom(n)
         f_end = start_of_string (truc)
         text = text(1:d_end) // '_' // truc(f_end:)
      END DO

      d_end = last_c_leng (20, text)
      IF (type_fe==1) THEN
         text = text(1:d_end) // '_FE_1'
      ELSE
         text = text(1:d_end) // '_FE_2'
      END IF
      !WRITE (*,*) 'Loading mesh-file ...'
      !d_end = last_c_leng (64, dir)
      !f_end = last_c_leng (64, fil)
      IF (mesh_formatted) THEN
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'formatted')
      ELSE
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'unformatted')
      END IF
      OPEN(UNIT = 20, FILE = text, FORM = 'formatted', STATUS = 'unknown')

      !     READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

      IF (type_fe == 2) THEN    ! Skip P1 data if needed
         IF (mesh_formatted) THEN
            READ(30, *) np, nw, me, nws, mes
         ELSE
            READ(30) np, nw, me, nws, mes
         END IF

         IF (mesh_formatted) THEN
            DO m = 1, me
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF

         IF (mesh_formatted) THEN
            DO ms = 1, mes
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF

         IF (mesh_formatted) THEN
            DO n = 1, np
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF
      END IF

      IF (mesh_formatted) THEN
         READ  (30, *)  np, nw, me, nws, mes
      ELSE
         READ(30)  np, nw, me, nws, mes
      END IF

      IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==6 .AND. nws==3) THEN
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE IF (nw==10 .AND. nws==6) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed '
         WRITE(*, *) kd, nw, nws
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



      !---  Renumerotation------------------------------------------------------------
      virgin_nd = .TRUE.
      virgin_el = .TRUE.
      mnouv = 0
      nnouv = 0
      DO dom = 1, SIZE(list_dom)
         DO m = 1, me
            IF (ABS(list_dom(dom) - i_d_lect(m)) /= 0)  CYCLE ! i_d(m) dans la liste
            virgin_el(m) = .FALSE.
            mnouv = mnouv + 1   ! Nouvel element
            nouv_el(m) = mnouv
            ancien_el(mnouv) = m
            DO n = 1, nw; i = jj_lect(n, m)
            IF (virgin_nd(i)) THEN ! Nouveau point
               virgin_nd(i) = .FALSE.
               nnouv = nnouv + 1
               nouv_nd(i) = nnouv
               ancien_nd(nnouv) = i
            END IF
            END DO
         END DO
      END DO
      mesh%me = mnouv
      mesh%np = nnouv

      ALLOCATE(mesh%jj(nw, mesh%me), mesh%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))

      DO m = 1, mesh%me
         mesh%jj(:, m) = nouv_nd(jj_lect(:, ancien_el(m)))
         mesh%neigh(:, m) = nouv_el(neigh_lect(:, ancien_el(m)))
         mesh%i_d(m) = i_d_lect(ancien_el(m))
      END DO

      !---  Fin renumerotation--------------------------------------------------------
      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))

      IF (mesh_formatted) THEN
         DO ms = 1, mes
            READ(30, *) jjs_lect(:, ms), neighs_lect(ms), sides_lect(ms)
         END DO
      ELSE
         READ(30) jjs_lect, neighs_lect, sides_lect
      END IF

      !---  Renumerotation------------------------------------------------------------
      ALLOCATE(nouv_els(mes), ancien_els(mes))
      DEALLOCATE(virgin_el)
      ALLOCATE(virgin_el(mes))
      virgin_el = .TRUE.
      msnouv = 0
      DO dom = 1, SIZE(list_dom)
         DO ms = 1, mes
            neighs1 = neighs_lect(ms)
            DO n = 1, nw
               IF (MINVAL(ABS(jjs_lect(:, ms) - jj_lect(n, neighs1))) /= 0) EXIT
            END DO
            neighs2 = neigh_lect(n, neighs1)
            test = .FALSE.
            IF (ABS(list_dom(dom) - i_d_lect(neighs1)) == 0) THEN
               test = .TRUE.
            ELSE IF (neighs2 /= 0) THEN
               IF (ABS(list_dom(dom) - i_d_lect(neighs2)) == 0) THEN
                  test = .TRUE.
                  neighs_lect(ms) = neighs2 ! On change de cote
               END IF
            END IF

            IF (.NOT.test) CYCLE
            !11 June 2007
            IF (.NOT.virgin_el(ms)) CYCLE
            !11 June 2007
            virgin_el(ms) = .FALSE.
            msnouv = msnouv + 1
            nouv_els(ms) = msnouv
            ancien_els(msnouv) = ms
         END DO
      END DO
      mesh%mes = msnouv

      ALLOCATE (mesh%jjs(nws, mesh%mes), mesh%neighs(mesh%mes), &
           mesh%sides(mesh%mes))

      DO ms = 1, mesh%mes
         mesh%jjs(:, ms) = nouv_nd(jjs_lect(:, ancien_els(ms)))
         mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
         mesh%sides(ms) = sides_lect(ancien_els(ms))
      END DO

      !---  Fin renumerotation--------------------------------------------------------

      ALLOCATE(rr_lect(kd, np))

      IF (mesh_formatted) THEN
         DO n = 1, np
            READ(30, *) rr_lect(:, n)
         END DO
      ELSE
         READ(30) rr_lect
      END IF

      ALLOCATE(mesh%rr(kd, mesh%np))

      mesh%rr = rr_lect(:, ancien_nd(1:mesh%np))

      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_el, virgin_nd)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els, ancien_nd, ancien_el, ancien_els)
      !     END OF GRID READING --------------------------------------------------------

      ALLOCATE(mesh%iis(nws, mesh%mes))
      CALL dirichlet_nodes(mesh%jjs, SPREAD(1, 1, mesh%mes), SPREAD(.TRUE., 1, 1), mesh%j_s)
      CALL surf_nodes_i(mesh%jjs, mesh%j_s, mesh%iis)
      mesh%nps = SIZE(mesh%j_s)

      WRITE (20, *)  'np_lect', np, 'nw_lect', nw, 'nws_lect', nws, 'me_lect', me, &
           'mes_lect', mes
      WRITE (20, *)  'np ', mesh%np, 'me ', mesh%me, &
           'mes ', mesh%mes, 'nps ', mesh%nps
      WRITE (20, *) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
      !     CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
      !     mesh%jj, mesh%jjs, mesh%rr)
      IF (PRESENT(edge_stab)) THEN
         mesh%edge_stab = edge_stab
         IF (mesh%edge_stab) CALL prep_interfaces(mesh) ! JLG Fevrier 2010
      ELSE
         mesh%edge_stab = .FALSE.
      END IF

      CLOSE(20)
      CLOSE(30)

   END SUBROUTINE load_mesh_free_format_ordered

   !------------------------------------------------------------------------------

   SUBROUTINE load_mesh_free_format(dir, fil, list_dom, type_fe, mesh, mesh_formatted, edge_stab)

      USE def_type_mesh
      USE chaine_caractere
      USE Dir_nodes

      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type) :: mesh
      LOGICAL, INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
      LOGICAL, OPTIONAL, INTENT(IN) :: edge_stab

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els, &
           ancien_nd, ancien_el, ancien_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_el
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect

      LOGICAL :: test
      INTEGER :: mnouv, nnouv, i
      INTEGER :: n, m, ms, msnouv, neighs1, neighs2
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh
      CHARACTER(len = 60) :: text
      CHARACTER(len = 2) :: truc

      text = 'Mesh'
      d_end = last_c_leng (20, text)
      DO n = 1, SIZE(list_dom)
         d_end = last_c_leng (20, text)
         WRITE(truc, '(i2)') list_dom(n)
         f_end = start_of_string (truc)
         text = text(1:d_end) // '_' // truc(f_end:)
      END DO

      d_end = last_c_leng (20, text)
      IF (type_fe==1) THEN
         text = text(1:d_end) // '_FE_1'
      ELSE
         text = text(1:d_end) // '_FE_2'
      END IF

      !WRITE (*,*) 'Loading mesh-file ...'
      !d_end = last_c_leng (64, dir)
      !f_end = last_c_leng (64, fil)
      IF (mesh_formatted) THEN
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'formatted')
      ELSE
         OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'unformatted')
      END IF
      OPEN(UNIT = 20, FILE = text, FORM = 'formatted', STATUS = 'unknown')

      !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

      IF (type_fe == 2) THEN ! Skip P1 data if needed
         IF (mesh_formatted) THEN
            READ(30, *) np, nw, me, nws, mes
         ELSE
            READ(30) np, nw, me, nws, mes
         END IF

         IF (mesh_formatted) THEN
            DO m = 1, me
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF

         IF (mesh_formatted) THEN
            DO ms = 1, mes
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF

         IF (mesh_formatted) THEN
            DO n = 1, np
               READ(30, *)
            END DO
         ELSE
            READ(30)
         END IF
      END IF

      IF (mesh_formatted) THEN
         READ  (30, *)  np, nw, me, nws, mes
      ELSE
         READ(30)  np, nw, me, nws, mes
      END IF

      IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==6 .AND. nws==3) THEN
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE IF (nw==10 .AND. nws==6) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed '
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

      !---Renumerotation------------------------------------------------------------
      virgin_nd = .TRUE.
      virgin_el = .TRUE.
      mnouv = 0
      nnouv = 0

      DO m = 1, me
         IF (MINVAL(ABS(list_dom - i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
         virgin_el(m) = .FALSE.
         mnouv = mnouv + 1  ! Nouvel element
         nouv_el(m) = mnouv
         ancien_el(mnouv) = m
         DO n = 1, nw; i = jj_lect(n, m)
         IF (virgin_nd(i)) THEN ! Nouveau point
            virgin_nd(i) = .FALSE.
            nnouv = nnouv + 1
            nouv_nd(i) = nnouv
            ancien_nd(nnouv) = i
         END IF
         END DO
      END DO

      mesh%me = mnouv
      mesh%np = nnouv

      ALLOCATE(mesh%jj(nw, mesh%me), mesh%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))

      DO m = 1, mesh%me
         mesh%jj(:, m) = nouv_nd(jj_lect(:, ancien_el(m)))
         mesh%neigh(:, m) = nouv_el(neigh_lect(:, ancien_el(m)))
         mesh%i_d(m) = i_d_lect(ancien_el(m))
      END DO

      !---Fin renumerotation--------------------------------------------------------
      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))

      IF (mesh_formatted) THEN
         DO ms = 1, mes
            READ(30, *) jjs_lect(:, ms), neighs_lect(ms), sides_lect(ms)
         END DO
      ELSE
         READ(30) jjs_lect, neighs_lect, sides_lect
      END IF

      !---Renumerotation------------------------------------------------------------
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
               neighs_lect(ms) = neighs2 ! On change de cote
            END IF
         END IF

         IF (.NOT.test) CYCLE
         virgin_el(ms) = .FALSE.
         msnouv = msnouv + 1
         nouv_els(ms) = msnouv
         ancien_els(msnouv) = ms
      END DO

      mesh%mes = msnouv

      ALLOCATE (mesh%jjs(nws, mesh%mes), mesh%neighs(mesh%mes), &
           mesh%sides(mesh%mes))

      DO ms = 1, mesh%mes
         mesh%jjs(:, ms) = nouv_nd(jjs_lect(:, ancien_els(ms)))
         mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
         mesh%sides(ms) = sides_lect(ancien_els(ms))
      END DO

      !---Fin renumerotation--------------------------------------------------------

      ALLOCATE(rr_lect(kd, np))

      IF (mesh_formatted) THEN
         DO n = 1, np
            READ(30, *) rr_lect(:, n)
         END DO
      ELSE
         READ(30) rr_lect
      END IF

      ALLOCATE(mesh%rr(kd, mesh%np))

      mesh%rr = rr_lect(:, ancien_nd(1:mesh%np))

      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_el, virgin_nd)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els, ancien_nd, ancien_el, ancien_els)
      !  END OF GRID READING --------------------------------------------------------

      ALLOCATE(mesh%iis(nws, mesh%mes))
      CALL dirichlet_nodes(mesh%jjs, SPREAD(1, 1, mesh%mes), SPREAD(.TRUE., 1, 1), mesh%j_s)
      CALL surf_nodes_i(mesh%jjs, mesh%j_s, mesh%iis)
      mesh%nps = SIZE(mesh%j_s)

      WRITE (20, *)  'np_lect', np, 'nw_lect', nw, 'nws_lect', nws, 'me_lect', me, &
           'mes_lect', mes
      WRITE (20, *)  'np ', mesh%np, 'me ', mesh%me, &
           'mes ', mesh%mes, 'nps ', mesh%nps
      WRITE (20, *) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
      !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
      !                    mesh%jj, mesh%jjs, mesh%rr)
      IF (PRESENT(edge_stab)) THEN
         mesh%edge_stab = edge_stab
         IF (mesh%edge_stab) CALL prep_interfaces(mesh) ! JLG April 2009
      ELSE
         mesh%edge_stab = .FALSE.
      END IF

      CLOSE(20)
      CLOSE(30)

   END SUBROUTINE load_mesh_free_format

   !------------------------------------------------------------------------------

   SUBROUTINE load_mesh_formatted(dir, fil, list_dom, type_fe, mesh)

      USE def_type_mesh
      USE chaine_caractere
      USE Dir_nodes

      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type) :: mesh

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els, &
           ancien_nd, ancien_el, ancien_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_el
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect

      LOGICAL :: test
      INTEGER :: mnouv, nnouv, i
      INTEGER :: n, m, ms, msnouv, neighs1, neighs2
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh

      d_end = last_c_leng (64, dir)
      f_end = last_c_leng (64, fil)

      !WRITE (*,*) 'Loading mesh-file ...'

      OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'formatted')
      !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
      OPEN(UNIT = 20, FILE = 'error_mesh', FORM = 'formatted', STATUS = 'unknown')

      !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

      IF (type_fe == 2) THEN ! Skip P1 data if needed
         READ(30, *) np, nw, me, nws, mes
         !READ(30) np, nw, me, nws, mes
         DO m = 1, me
            READ(30, *)
         END DO
         !READ(30)
         DO ms = 1, mes
            READ(30, *)
         END DO
         !READ(30)
         DO n = 1, np
            READ(30, *)
         END DO
         !READ(30)
      END IF

      READ  (30, *)  np, nw, me, nws, mes
      !READ  (30)  np,  nw,  me,  nws,  mes

      IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==6 .AND. nws==3) THEN
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE IF (nw==10 .AND. nws==6) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed '
         STOP
      END IF

      ALLOCATE (jj_lect(nw, me), neigh_lect(nwneigh, me), i_d_lect(me))
      ALLOCATE (nouv_nd(np), ancien_nd(np), virgin_nd(np), &
           nouv_el(0:me), ancien_el(me), virgin_el(me))

      nouv_el = 0

      DO m = 1, me
         READ(30, *) jj_lect(:, m), neigh_lect(:, m), i_d_lect(m)
      END DO

      !READ(30) jj_lect, neigh_lect, i_d_lect

      !---Renumerotation------------------------------------------------------------
      virgin_nd = .TRUE.
      virgin_el = .TRUE.
      mnouv = 0
      nnouv = 0

      DO m = 1, me
         IF (MINVAL(ABS(list_dom - i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
         virgin_el(m) = .FALSE.
         mnouv = mnouv + 1  ! Nouvel element
         nouv_el(m) = mnouv
         ancien_el(mnouv) = m
         DO n = 1, nw; i = jj_lect(n, m)
         IF (virgin_nd(i)) THEN ! Nouveau point
            virgin_nd(i) = .FALSE.
            nnouv = nnouv + 1
            nouv_nd(i) = nnouv
            ancien_nd(nnouv) = i
         END IF
         END DO
      END DO

      mesh%me = mnouv
      mesh%np = nnouv

      ALLOCATE(mesh%jj(nw, mesh%me), mesh%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))

      DO m = 1, mesh%me
         mesh%jj(:, m) = nouv_nd(jj_lect(:, ancien_el(m)))
         mesh%neigh(:, m) = nouv_el(neigh_lect(:, ancien_el(m)))
         mesh%i_d(m) = i_d_lect(ancien_el(m))
      END DO

      !---Fin renumerotation--------------------------------------------------------
      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))

      DO ms = 1, mes
         READ(30, *) jjs_lect(:, ms), neighs_lect(ms), sides_lect(ms)
      END DO

      !READ(30) jjs_lect, neighs_lect, sides_lect

      !---Renumerotation------------------------------------------------------------
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
               neighs_lect(ms) = neighs2 ! On change de cote
            END IF
         END IF
         IF (.NOT.test) CYCLE
         virgin_el(ms) = .FALSE.
         msnouv = msnouv + 1
         nouv_els(ms) = msnouv
         ancien_els(msnouv) = ms
      END DO

      mesh%mes = msnouv

      ALLOCATE (mesh%jjs(nws, mesh%mes), mesh%neighs(mesh%mes), &
           mesh%sides(mesh%mes))

      DO ms = 1, mesh%mes
         mesh%jjs(:, ms) = nouv_nd(jjs_lect(:, ancien_els(ms)))
         mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
         mesh%sides(ms) = sides_lect(ancien_els(ms))
      END DO

      !---Fin renumerotation--------------------------------------------------------

      ALLOCATE(rr_lect(kd, np))

      DO n = 1, np
         READ(30, *) rr_lect(:, n)
      END DO
      !READ(30) rr_lect

      ALLOCATE(mesh%rr(kd, mesh%np))

      mesh%rr = rr_lect(:, ancien_nd(1:mesh%np))

      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_el, virgin_nd)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els, ancien_nd, ancien_el, ancien_els)
      !  END OF GRID READING --------------------------------------------------------

      ALLOCATE(mesh%iis(nws, mesh%mes))
      CALL dirichlet_nodes(mesh%jjs, SPREAD(1, 1, mesh%mes), SPREAD(.TRUE., 1, 1), mesh%j_s)
      CALL surf_nodes_i(mesh%jjs, mesh%j_s, mesh%iis)
      mesh%nps = SIZE(mesh%j_s)

      WRITE (20, *)  'np_lect', np, 'nw_lect', nw, 'nws_lect', nws, 'me_lect', me, &
           'mes_lect', mes
      WRITE (20, *)  'np ', mesh%np, 'me ', mesh%me, &
           'mes ', mesh%mes, 'nps ', mesh%nps
      WRITE (20, *) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
      !   CALL Gauss_gen(mesh%np, mesh%me, mesh%nps, mesh%mes, &
      !                    mesh%jj, mesh%jjs, mesh%rr)

      CLOSE(20)
      CLOSE(30)

   END SUBROUTINE load_mesh_formatted

   SUBROUTINE load_mesh(dir, fil, list_dom, type_fe, mesh)

      USE def_type_mesh
      USE chaine_caractere
      !USE gauss_points
      USE Dir_nodes

      IMPLICIT NONE
      CHARACTER(len = 200), INTENT(IN) :: dir, fil
      INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
      INTEGER, INTENT(IN) :: type_fe
      TYPE(mesh_type) :: mesh

      INTEGER, ALLOCATABLE, DIMENSION(:) :: nouv_nd, nouv_el, nouv_els, &
           ancien_nd, ancien_el, ancien_els
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: jj_lect, neigh_lect, jjs_lect
      INTEGER, ALLOCATABLE, DIMENSION(:) :: i_d_lect, sides_lect, neighs_lect
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: virgin_nd, virgin_el
      REAL(KIND = 8), ALLOCATABLE, DIMENSION(:, :) :: rr_lect

      LOGICAL :: test
      INTEGER :: mnouv, nnouv, i
      INTEGER :: n, m, ms, msnouv, neighs1, neighs2
      INTEGER :: d_end, f_end
      INTEGER :: np, nw, me, nws, mes, kd, nwneigh
      CHARACTER(len = 20) :: text
      CHARACTER(len = 2) :: truc

      text = 'Mesh'
      d_end = last_c_leng (20, text)
      DO n = 1, SIZE(list_dom)
         d_end = last_c_leng (20, text)
         WRITE(truc, '(i2)') list_dom(n)
         f_end = start_of_string (truc)
         text = text(1:d_end) // '_' // truc(f_end:)
      END DO

      d_end = last_c_leng (20, text)
      IF (type_fe==1) THEN
         text = text(1:d_end) // '_FE_1'
      ELSE
         text = text(1:d_end) // '_FE_2'
      END IF

      !WRITE (*,*) 'Loading mesh-file ...'
      d_end = last_c_leng (64, dir)
      f_end = last_c_leng (64, fil)
      !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
      OPEN(30, FILE = TRIM(ADJUSTL(dir)) // '/' // TRIM(ADJUSTL(fil)), FORM = 'unformatted')
      OPEN(UNIT = 20, FILE = text, FORM = 'formatted', STATUS = 'unknown')

      !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

      IF (type_fe == 2) THEN ! Skip P1 data if needed
         !READ(30,*) np, nw, me, nws, mes
         READ(30) np, nw, me, nws, mes

         !DO m = 1, me
         !   READ(30,*)
         !END DO
         READ(30)

         !DO ms = 1, mes
         !   READ(30,*)
         !END DO

         ! DO m = 1, mes
         !     READ(30)
         ! END DO
         READ(30)

         !DO n = 1, np
         !   READ(30,*)
         !END DO
         READ(30)
      END IF

      !READ  (30, *)  np,  nw,  me,  nws,  mes
      READ(30)  np, nw, me, nws, mes

      IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
         kd = 2; nwneigh = 3
      ELSE IF (nw==6 .AND. nws==3) THEN
         kd = 2; nwneigh = 3
      ELSE IF (nw==4 .AND. nws==3) THEN
         kd = 3; nwneigh = 4
      ELSE IF (nw==10 .AND. nws==6) THEN
         kd = 3; nwneigh = 4
      ELSE
         WRITE(*, *) ' Finite element not yet programmed '
         STOP
      END IF

      ALLOCATE (jj_lect(nw, me), neigh_lect(nwneigh, me), i_d_lect(me))
      ALLOCATE (nouv_nd(np), ancien_nd(np), virgin_nd(np), &
           nouv_el(0:me), ancien_el(me), virgin_el(me))

      nouv_el = 0

      !DO m = 1, me
      !   READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
      !END DO

      READ(30) jj_lect, neigh_lect, i_d_lect

      !---Renumerotation------------------------------------------------------------
      virgin_nd = .TRUE.
      virgin_el = .TRUE.
      mnouv = 0
      nnouv = 0

      DO m = 1, me
         IF (MINVAL(ABS(list_dom - i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
         virgin_el(m) = .FALSE.
         mnouv = mnouv + 1  ! Nouvel element
         nouv_el(m) = mnouv
         ancien_el(mnouv) = m
         DO n = 1, nw; i = jj_lect(n, m)
         IF (virgin_nd(i)) THEN ! Nouveau point
            virgin_nd(i) = .FALSE.
            nnouv = nnouv + 1
            nouv_nd(i) = nnouv
            ancien_nd(nnouv) = i
         END IF
         END DO
      END DO

      mesh%me = mnouv
      mesh%np = nnouv

      ALLOCATE(mesh%jj(nw, mesh%me), mesh%neigh(nwneigh, mesh%me), mesh%i_d(mesh%me))

      DO m = 1, mesh%me
         mesh%jj(:, m) = nouv_nd(jj_lect(:, ancien_el(m)))
         mesh%neigh(:, m) = nouv_el(neigh_lect(:, ancien_el(m)))
         mesh%i_d(m) = i_d_lect(ancien_el(m))
      END DO

      !---Fin renumerotation--------------------------------------------------------
      ALLOCATE (jjs_lect(nws, mes), neighs_lect(mes), sides_lect(mes))

      !DO ms = 1, mes
      !   READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
      !END DO


      !TEST
      !DO ms = 1, mes
      !READ(30) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
      !END DO
      !TEST
      READ(30) jjs_lect, neighs_lect, sides_lect

      !---Renumerotation------------------------------------------------------------
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
               neighs_lect(ms) = neighs2 ! On change de cote
            END IF
         END IF

         IF (.NOT.test) CYCLE
         virgin_el(ms) = .FALSE.
         msnouv = msnouv + 1
         nouv_els(ms) = msnouv
         ancien_els(msnouv) = ms
      END DO

      mesh%mes = msnouv

      ALLOCATE (mesh%jjs(nws, mesh%mes), mesh%neighs(mesh%mes), &
           mesh%sides(mesh%mes))

      DO ms = 1, mesh%mes
         mesh%jjs(:, ms) = nouv_nd(jjs_lect(:, ancien_els(ms)))
         mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
         mesh%sides(ms) = sides_lect(ancien_els(ms))
      END DO

      !---Fin renumerotation--------------------------------------------------------

      ALLOCATE(rr_lect(kd, np))

      !DO n = 1, np
      !   READ(30,*) rr_lect(:,n)
      !END DO

      READ(30) rr_lect

      ALLOCATE(mesh%rr(kd, mesh%np))

      mesh%rr = rr_lect(:, ancien_nd(1:mesh%np))

      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_el, virgin_nd)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els, ancien_nd, ancien_el, ancien_els)
      !  END OF GRID READING --------------------------------------------------------

      ALLOCATE(mesh%iis(nws, mesh%mes))
      CALL dirichlet_nodes(mesh%jjs, SPREAD(1, 1, mesh%mes), SPREAD(.TRUE., 1, 1), mesh%j_s)
      CALL surf_nodes_i(mesh%jjs, mesh%j_s, mesh%iis)
      mesh%nps = SIZE(mesh%j_s)

      WRITE (20, *)  'np_lect', np, 'nw_lect', nw, 'nws_lect', nws, 'me_lect', me, &
           'mes_lect', mes
      WRITE (20, *)  'np ', mesh%np, 'me ', mesh%me, &
           'mes ', mesh%mes, 'nps ', mesh%nps
      WRITE (20, *) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
      CLOSE(20)
      CLOSE(30)

   END SUBROUTINE load_mesh


   SUBROUTINE surf_nodes_i(jjs, j_s, iis)

      !  generation of the surface element connectivity matrix  iis
      !  based on the surface node numbering, starting from the
      !  connectivity matrix  jjs  of the surface elements according
      !  to the volume node numbering, and from the array  j_s  of
      !  the boundary nodes according to the volume node numbering

      IMPLICIT NONE

      INTEGER, DIMENSION(:, :), INTENT(IN) :: jjs
      INTEGER, DIMENSION(:), INTENT(IN) :: j_s
      INTEGER, DIMENSION(:, :), INTENT(OUT) :: iis

      INTEGER :: ms, ls, j, i

      DO ms = 1, SIZE(jjs, 2)
         DO ls = 1, SIZE(jjs, 1)
            j = jjs(ls, ms)
            DO i = 1, SIZE(j_s)
               IF (j == j_s(i))  iis(ls, ms) = i
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE surf_nodes_i

   SUBROUTINE prep_interfaces(mesh)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      LOGICAL, DIMENSION(mesh%me) :: virgin
      INTEGER :: m, mop, nw, me, l, lop, n, n1, n2, edge, nt, nws
      nw = SIZE(mesh%jj, 1)
      nws = SIZE(mesh%jjs, 1)
      me = mesh%me
      IF (SIZE(mesh%rr, 1)==2) THEN
         nt = 3
      ELSE
         WRITE(*, *) ' BUG: prep_interfaces, 3D not programmed yet '
         STOP
         nt = 4
      END IF

      virgin = .TRUE.
      edge = 0
      DO m = 1, me
         virgin(m) = .FALSE.
         DO n = 1, nt
            mop = mesh%neigh(n, m)
            IF (mop==0) CYCLE !Edge on boundary
            IF (.NOT.virgin(mop)) CYCLE !Edge already done
            edge = edge + 1 !New edge
         END DO
      END DO
      IF (SIZE(mesh%rr, 1)==2) THEN
         IF (edge/=(3 * mesh%me - mesh%mes) / 2) THEN
            WRITE(*, *) ' BUG in prep_interfaces, edge/=(3*mesh%me + mesh%mes)/2'
            WRITE(*, *) ' edge ', edge, (3 * mesh%me - mesh%mes) / 2
            WRITE(*, *) ' mesh%mes ', mesh%mes, ' mesh%me ', mesh%me
            STOP
         END IF
      END IF

      mesh%mi = edge
      ALLOCATE(mesh%jji(nw, 2, edge), mesh%jjsi(nws, edge), mesh%neighi(2, edge))

      edge = 0
      virgin = .TRUE.
      DO m = 1, me
         virgin(m) = .FALSE.
         DO n = 1, nt
            mop = mesh%neigh(n, m)
            IF (mop==0) CYCLE !Edge on boundary
            IF (.NOT.virgin(mop)) CYCLE !Edge already done
            edge = edge + 1 !New edge
            n1 = MODULO(n, nt) + 1
            n2 = MODULO(n + 1, nt) + 1 ! Works in 2D only
            mesh%jjsi(1, edge) = mesh%jj(n1, m)
            mesh%jjsi(2, edge) = mesh%jj(n2, m)
            IF (nws==3) THEN
               IF (n1 + n2==3) THEN
                  mesh%jjsi(3, edge) = mesh%jj(6, m)
               ELSE IF (n1 + n2==5) THEN
                  mesh%jjsi(3, edge) = mesh%jj(4, m)
               ELSE IF (n1 + n2==4) THEN
                  mesh%jjsi(3, edge) = mesh%jj(5, m)
               ELSE
                  WRITE(*, *) ' BUG prep_interfaces n1+n2 not correct '
                  STOP
               END IF
            END IF
            IF (m<mop) THEN
               l = 1 ! Side 1 on smallest cell index
               lop = 2
               mesh%neighi(1, edge) = m
               mesh%neighi(2, edge) = mop
            ELSE
               l = 2 ! Side 2 on largest cell index
               lop = 1
               mesh%neighi(1, edge) = mop
               mesh%neighi(2, edge) = m
            END IF
            DO n1 = 1, nw
               mesh%jji(n1, l, edge) = mesh%jj(n1, m)
               mesh%jji(n1, lop, edge) = mesh%jj(n1, mop)
            END DO
         END DO
      END DO

   END SUBROUTINE prep_interfaces


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
      write(*, *) '????', mesh_p1%neighs_int
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
         IF (mesh%jjs_int(1, mes_int + ms) > mesh_p1%dom_np) mesh%jjs_int(1, mes_int + ms) = mesh%jjs_int(1, mes + ms) &
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


   SUBROUTINE prep_jce_jev(mesh)
      USE def_type_mesh
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh
      LOGICAL, DIMENSION(mesh%me) :: virgin
      INTEGER :: m, mop, nw, me, l, lop, n, n1, n2, nmin, nmax, edge, nt, nws, f_dof, start, endf, nop
      LOGICAL :: test
      nw = SIZE(mesh%jj, 1)
      nws = SIZE(mesh%jjs, 1)
      me = mesh%me

      IF (SIZE(mesh%rr, 1)==2) THEN
         nt = 3
         f_dof = nws - 2
      ELSE
         WRITE(*, *) ' BUG: prep_jce_jev, 3D not programmed yet '
         STOP
         nt = 4
      END IF
      IF (nw.NE.nt) THEN
         WRITE(*, *) ' BUG in prep_jce_jev, nw.NE.nt'
         STOP
      END IF
      virgin = .TRUE.
      edge = 0
      DO m = 1, me
         virgin(m) = .FALSE.
         DO n = 1, nt
            mop = mesh%neigh(n, m)
            IF (mop>0) THEN
               IF (.NOT.virgin(mop)) CYCLE !Edge already done
            END IF
            edge = edge + 1 !New edge
         END DO
      END DO
      IF (SIZE(mesh%rr, 1)==2) THEN
         IF (edge/=(3 * mesh%me - mesh%mes) / 2 + mesh%mes) THEN
            WRITE(*, *) ' BUG in prep_jce_jev, edge/=(3*mesh%me - mesh%mes)/2+mesh%mes'
            WRITE(*, *) ' edge ', edge, (3 * mesh%me - mesh%mes) / 2 + mesh%mes
            WRITE(*, *) ' mesh%mes ', mesh%mes, ' mesh%me ', mesh%me
            STOP
         END IF
      END IF

      mesh%medge = edge
      !ALLOCATE(mesh%jev(nt - 1, mesh%medge))
      ALLOCATE(mesh%jce(nt, mesh%me))

      edge = 0
      virgin = .TRUE.
      DO m = 1, me
         virgin(m) = .FALSE.
         DO n = 1, nt
            mop = mesh%neigh(n, m)
            IF (mop>0) THEN
               IF (.NOT.virgin(mop)) THEN
                  n1 = MODULO(n, nt) + 1
                  n2 = MODULO(n + 1, nt) + 1
                  test = .false.
                  DO nop = 1, nw
                     IF (MIN(ABS(mesh%jj(nop, mop) - mesh%jj(n1, m)), ABS(mesh%jj(nop, mop) - mesh%jj(n2, m)))==0) THEN
                        CYCLE
                     ELSE
                        test = .true.
                        EXIT
                     END IF
                  END DO
                  IF (.NOT.test) THEN
                     WRITE(*, *) ' BUG in prep_jce_jev'
                     STOP
                  END IF
                  mesh%jce(n, m) = mesh%jce(nop, mop)
                  CYCLE !Edge already done
               END IF
            END IF
            edge = edge + 1 !New edge
            n1 = MODULO(n, nt) + 1
            n2 = MODULO(n + 1, nt) + 1 ! Works in 2D only
            nmin = MIN(n1, n2)
            nmax = MAX(n1, n2)
            !mesh%jev(1, edge) = mesh%jj(nmin, m)
            !mesh%jev(2, edge) = mesh%jj(nmax, m)
            mesh%jce(n, m) = edge
         END DO
      END DO
   END SUBROUTINE prep_jce_jev

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

   SUBROUTINE switch_cell_vertices(jj, neigh, m, i1, i2)
      INTEGER, DIMENSION(:, :), INTENT(INOUT) :: jj, neigh
      INTEGER :: m, i1, i2, save
      save = jj(i1, m)
      jj(i1, m) = jj(i2, m)
      jj(i2, m) = save

      save = neigh(i1, m)
      neigh(i1, m) = neigh(i2, m)
      neigh(i2, m) = save

   END SUBROUTINE switch_cell_vertices

   SUBROUTINE copy_mesh(mesh1, mesh2)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh1, mesh2

      mesh2%me = mesh1%me
      mesh2%mes = mesh1%mes
      mesh2%mes_int = mesh1%mes
      mesh2%np = mesh1%np
      mesh2%nps = mesh1%nps
      mesh2%mi = mesh1%mi
      mesh2%medge = mesh1%medge
      mesh2%medges = mesh1%medges
      mesh2%mextra = mesh1%mextra
      mesh2%mes_extra = mesh1%mes_extra

      mesh2%dom_me = mesh1%dom_me
      mesh2%dom_np = mesh1%dom_np
      mesh2%dom_mes = mesh1%dom_mes
      mesh2%nis = mesh1%nis

      ALLOCATE(mesh2%jj(SIZE(mesh1%jj, 1), SIZE(mesh1%jj, 2)))
      mesh2%jj = mesh1%jj
      ALLOCATE(mesh2%jjs(SIZE(mesh1%jjs, 1), SIZE(mesh1%jjs, 2)))
      mesh2%jjs = mesh1%jjs
      ALLOCATE(mesh2%rr(SIZE(mesh1%rr, 1), SIZE(mesh1%rr, 2)))
      mesh2%rr = mesh1%rr
      ALLOCATE(mesh2%loc_to_glob(SIZE(mesh1%loc_to_glob)))
      mesh2%loc_to_glob = mesh1%loc_to_glob
      ALLOCATE(mesh2%neigh(SIZE(mesh1%neigh, 1), SIZE(mesh1%neigh, 2)))
      mesh2%neigh = mesh1%neigh

      ALLOCATE(mesh2%sides(SIZE(mesh1%sides)))
      mesh2%sides = mesh1%sides
      ALLOCATE(mesh2%neighs(SIZE(mesh1%neighs)))
      mesh2%neighs = mesh1%neighs
      ALLOCATE(mesh2%i_d(SIZE(mesh1%i_d)))
      mesh2%i_d = mesh1%i_d

      ALLOCATE(mesh2%jjs_int(SIZE(mesh1%jjs_int, 1), SIZE(mesh1%jjs_int, 2)))
      mesh2%jjs_int = mesh1%jjs_int
      ALLOCATE(mesh2%sides_int(SIZE(mesh1%sides_int)))
      mesh2%sides_int = mesh1%sides_int
      ALLOCATE(mesh2%neighs_int(2,SIZE(mesh1%neighs_int, 2)))
      mesh2%neighs_int = mesh1%neighs_int

      ALLOCATE(mesh2%sides_extra(SIZE(mesh1%sides_extra)))
      mesh2%sides_extra = mesh1%sides_extra
      ALLOCATE(mesh2%neighs_extra(SIZE(mesh1%neighs_extra)))
      mesh2%neighs_extra = mesh1%neighs_extra

      ALLOCATE(mesh2%jce(SIZE(mesh1%jce, 1), SIZE(mesh1%jce, 2)))
      mesh2%jce = mesh1%jce
      ALLOCATE(mesh2%jees(SIZE(mesh1%jees)))
      mesh2%jees = mesh1%jees
      ALLOCATE(mesh2%jecs(SIZE(mesh1%jecs)))
      mesh2%jecs = mesh1%jecs

      ALLOCATE(mesh2%disp(SIZE(mesh1%disp)))
      mesh2%disp = mesh1%disp
      ALLOCATE(mesh2%domnp(SIZE(mesh1%domnp)))
      mesh2%domnp = mesh1%domnp
      ALLOCATE(mesh2%discell(SIZE(mesh1%disp)))
      mesh2%discell = mesh1%discell
      ALLOCATE(mesh2%domcell(SIZE(mesh1%domnp)))
      mesh2%domcell = mesh1%domcell
      ALLOCATE(mesh2%disedge(SIZE(mesh1%disp)))
      mesh2%disedge = mesh1%disedge
      ALLOCATE(mesh2%domedge(SIZE(mesh1%domnp)))
      mesh2%domedge = mesh1%domedge

      ALLOCATE(mesh2%jj_extra(SIZE(mesh1%jj_extra, 1), SIZE(mesh1%jj_extra, 2)))
      mesh2%jj_extra = mesh1%jj_extra
      ALLOCATE(mesh2%jce_extra(SIZE(mesh1%jce_extra, 1), SIZE(mesh1%jce_extra, 2)))
      mesh2%jce_extra = mesh1%jce_extra
      ALLOCATE(mesh2%jcc_extra(SIZE(mesh1%jcc_extra)))
      mesh2%jcc_extra = mesh1%jcc_extra

      ALLOCATE(mesh2%jjs_extra(SIZE(mesh1%jjs_extra, 1), SIZE(mesh1%jjs_extra, 2)))
      mesh2%jjs_extra = mesh1%jjs_extra
      ALLOCATE(mesh2%rrs_extra(SIZE(mesh1%rrs_extra, 1), SIZE(mesh1%rrs_extra, 2), SIZE(mesh1%rrs_extra, 3)))
      mesh2%rrs_extra = mesh1%rrs_extra

      ALLOCATE(mesh2%isolated_jjs(SIZE(mesh1%isolated_jjs)))
      mesh2%isolated_jjs = mesh1%isolated_jjs
      ALLOCATE(mesh2%isolated_interfaces(SIZE(mesh1%isolated_interfaces, 1), SIZE(mesh1%isolated_interfaces, 2)))
      mesh2%isolated_interfaces = mesh1%isolated_interfaces
   END SUBROUTINE copy_mesh

   SUBROUTINE free_mesh(mesh)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh

      DEALLOCATE(mesh%jj, mesh%i_d, mesh%loc_to_glob, mesh%rr, mesh%neigh)
      DEALLOCATE(mesh%jjs, mesh%sides, mesh%neighs)
      DEALLOCATE(mesh%jjs_int, mesh%sides_int, mesh%neighs_int)
      DEALLOCATE(mesh%disp, mesh%domnp, mesh%disedge, mesh%domedge, mesh%discell, mesh%domcell)
      DEALLOCATE(mesh%jce, mesh%jees, mesh%jecs)

      DEALLOCATE(mesh%jj_extra, mesh%jce_extra, mesh%jjs_extra, mesh%jcc_extra, mesh%rrs_extra)
      DEALLOCATE(mesh%sides_extra, mesh%neighs_extra) !interfaces
      DEALLOCATE(mesh%isolated_jjs, mesh%isolated_interfaces)

      !IF (mesh%edge_stab) THEN
      !   DEALLOCATE(mesh%iis)
      !   NULLIFY(mesh%jji)
      !   DEALLOCATE(mesh%jjsi)
      !   DEALLOCATE(mesh%neighi)
      !END IF

   END SUBROUTINE FREE_MESH
END MODULE prep_maill
