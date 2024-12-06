!
!Authors Jean-Luc Guermond, Luigi Quarapelle, Copyrights 1994, 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE load_mesh

   IMPLICIT NONE

   PUBLIC :: load_dg_mesh_free_format
   PRIVATE

CONTAINS

   SUBROUTINE load_dg_mesh_free_format(dir, fil, list_dom, list_inter, type_fe, &
        mesh, mesh_formatted)
      USE def_type_mesh
      USE chaine_caractere
      USE dir_nodes
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
            WRITE(*, *) ' BUG in load_mesh, neighs1=0 '
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
               WRITE(*, *) ' BUG in load_mesh, stat out of bounds '
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
      mesh%neighs_int = -1

      virgin_nd = .TRUE.
      virgin_ms = .TRUE.
      nnouv = 0
      msnouv = 0
      mes_int = 0
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
               IF (virgin_ms(ms)) THEN
                  mes_int = mes_int + 1
                  mesh%jjs_int(:, mes_int) = nouv_nd(jjs_lect(:, ms))
                  mesh%neighs_int(1, mes_int) = nouv_el(neighs_lect(ms))
                  mesh%sides_int(mes_int) = sides_lect(ms)
                  virgin_ms(ms) = .FALSE.
               END IF
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

      DO ms = 1, mesh%mes_int
         DO n = 1, 3
            IF (MINVAL(ABS(mesh%jjs_int(:, ms) - mesh%jj(n, mesh%neighs_int(1, ms)))) /= 0) EXIT ! n not on the interface
         END DO
         mesh%neighs_int(2, ms) = mesh%neigh(n, mesh%neighs_int(1, ms))
      END DO
      !===End reordring
      mesh%gauss%n_ws = SIZE(mesh%jjs, 1)
      mesh%gauss%n_w = SIZE(mesh%jj, 1)
      DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
      DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
      DEALLOCATE(rr_lect, virgin_nd, virgin_ms, stat)
      DEALLOCATE(nouv_nd, nouv_el, nouv_els)
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

END MODULE load_mesh
