!
!Authors Jean-Luc Guermond, Luigi Quarapelle, Copyrights 1994, 2005
!Revised June 2008, Jean-Luc Guermond
!
MODULE prep_maill

  IMPLICIT NONE

  PUBLIC :: load_mesh, load_mesh_formatted, load_mesh_free_format, &
       load_dg_mesh_free_format, load_mesh_free_format_ordered, prep_interfaces, &
       create_p3_mesh, incr_vrtx_indx_enumeration, incr_vrtx_indx_enumeration_for_interfaces
  PRIVATE

CONTAINS

  !------------------------------------------------------------------------------
  SUBROUTINE create_p3_mesh(p1_mesh,p2_mesh,p3_mesh,type_fe)
    USE def_type_mesh
    USE chaine_caractere
    USE Dir_nodes
    USE my_util
    IMPLICIT NONE
    INTEGER,          INTENT(IN)      :: type_fe
    TYPE(mesh_type),  INTENT(IN)      :: p1_mesh
    TYPE(mesh_type),  INTENT(IN)      :: p2_mesh
    TYPE(mesh_type),  INTENT(OUT)     :: p3_mesh
    LOGICAL, DIMENSION(p2_mesh%np)    :: p2_virgin
    INTEGER :: v_dof, f_dof, c_dof, n_dof, edge, i, h, k, l, m, ms, n, n1, n2, j1, j2, &
         nb_edges, nb_vertices, lmax, sgn, sgn_op, j1_op, j2_op, n1_op, n2_op, n_op, m_op, hh, hh_op, &
         n_mid, j_mid
    REAL(KIND=8) :: s1, s2, s3, sh, sk
    REAL(KIND=8), DIMENSION(2) :: r1, r2, r3
    !===Programmed for 2d Lagrage elements (type_fe can be arbitrary)
    v_dof = 3 !===Vertex dofs
    f_dof = type_fe-1 !===Face dofs
    c_dof = (type_fe-2)*(type_fe-1)/2 !===Cell dofs
    p3_mesh%gauss%n_w = (type_fe+1)*(type_fe+2)/2
    p3_mesh%gauss%n_ws = type_fe+1
    p3_mesh%gauss%k_d = 2

    p3_mesh%me  = p1_mesh%me
    p3_mesh%mes = p1_mesh%mes
    !===Nb of edges
    edge = 0
    DO m = 1, p1_mesh%me
       DO n = 1, 3
          IF (p1_mesh%neigh(n,m)==0) CYCLE
          edge = edge + 1
       END DO
    END DO
    edge = edge/2
    nb_edges = edge + p1_mesh%mes
    IF (edge /= (3*p1_mesh%me - p1_mesh%mes)/2) THEN
       CALL error_PETSC('BUG in create_p3_mesh, nb of edges is wrong')
    END IF
    !===Nb of P3 nodes
    nb_vertices = p2_mesh%np - nb_edges
    p3_mesh%np = nb_vertices + 2*nb_edges + p1_mesh%me

    IF (nb_vertices.NE.p1_mesh%np) THEN
       CALL error_petsc(' BUG in create_p3_mesh, nb_vertices.NE.p1_mesh%np')
    END IF

    !===Allocate mesh structure
    ALLOCATE(p3_mesh%jj(p3_mesh%gauss%n_w,p3_mesh%me), p3_mesh%neigh(3,p3_mesh%me), p3_mesh%i_d(p3_mesh%me))
    ALLOCATE(p3_mesh%jjs(p3_mesh%gauss%n_ws,p3_mesh%mes), p3_mesh%neighs(p3_mesh%mes), p3_mesh%sides(p3_mesh%mes))
    ALLOCATE(p3_mesh%rr(p3_mesh%gauss%k_d,p3_mesh%np))

    !===Easy stuff
    p3_mesh%neigh  = p1_mesh%neigh
    p3_mesh%neighs = p1_mesh%neighs
    p3_mesh%i_d    = p1_mesh%i_d
    p3_mesh%sides  = p1_mesh%sides

    !===Create vertices
    p3_mesh%jj(1:3,:) = p1_mesh%jj(1:3,:)
    IF (MAXVAL(p3_mesh%jj(1:3,:)).NE.nb_vertices) THEN
       write(*,*) MAXVAL(p1_mesh%jj(1:3,p1_mesh%me)), nb_vertices
       write(*,*) MAXVAL(p3_mesh%jj(1:3,p1_mesh%me)), nb_vertices
       CALL error_PETSC('BUG in create_p3_mesh, Vertex enumeration not done properly')
    END IF
    p3_mesh%rr(:,1:nb_vertices) = p1_mesh%rr(:,1:nb_vertices)

    !===Create face dofs
    p2_virgin = .TRUE.
    n_dof = nb_vertices
    DO m = 1, p3_mesh%me
       DO n = 1, v_dof !===n is the face number of the face to be populated
          n_mid = n+3 !===n_mid is the mid point on the p2 mesh
          j_mid = p2_mesh%jj(n_mid,m)
          IF (.NOT.p2_virgin(j_mid)) THEN !===This edge has already been visited
             m_op = p1_mesh%neigh(n,m)
             DO n_op = 1, v_dof
                IF (m == p1_mesh%neigh(n_op,m_op)) THEN
                   EXIT !===n_op is the index of the face I am looking for
                END IF
             END DO
             n1 = MODULO(n,3)+1
             n2 = MODULO(n+1,3)+1
             n1_op = MODULO(n_op,3)+1
             n2_op = MODULO(n_op+1,3)+1
             j1 = p1_mesh%jj(n1,m)
             j1_op = p1_mesh%jj(n1_op,m_op)
             j2 = p1_mesh%jj(n2,m)
             j2_op = p1_mesh%jj(n2_op,m_op)
             IF (j1.NE.j1_op) THEN
                i = n2
                n2 = n1
                n1 = i
             END IF
             IF (n2-n1>0) THEN
                sgn = 1
             ELSE
                sgn = -1
             END IF
             IF (n2_op-n1_op>0) THEN
                sgn_op = 1
             ELSE
                sgn_op = -1
             END IF
             DO h = 1, f_dof
                hh = h*sgn + ((1-sgn)/2)*(f_dof+1)
                hh_op = h*sgn_op + ((1-sgn_op)/2)*(f_dof+1)
                p3_mesh%jj(v_dof+(n-1)*f_dof+hh, m) = p3_mesh%jj(v_dof+(n_op-1)*f_dof+hh_op, m_op)
             END DO
          ELSE
             p2_virgin(j_mid) = .FALSE.
             n1 = MODULO(n,3)+1
             n2 = MODULO(n+1,3)+1
             j1 = p2_mesh%jj(n1,m)
             j2 = p2_mesh%jj(n2,m)
             r1 = p2_mesh%rr(:,j1)
             r2 = p2_mesh%rr(:,j2)
             r3 = p2_mesh%rr(:,j_mid) !===middle node is P2 to have curved boundaries

             !===Create face dofs
             DO h = 1, f_dof
                n_dof = n_dof + 1 !===New index created
                p3_mesh%jj(v_dof+(n-1)*f_dof+h, m) = n_dof
                s1 = 0.d0
                s2 = 1.d0
                s3 = 0.5d0
                IF (n1<n2) THEN
                   sh = h/DBLE(type_fe) !===Move from r1 to r2
                ELSE
                   sh = 1-h/DBLE(type_fe) !===Move from r2 to r1
                END IF
                p3_mesh%rr(:, n_dof) = r1(:)*(sh - s2)*(sh - s3)/((s1 - s2)*(s1 - s3)) &
                     + r2(:)*(sh - s3)*(sh - s1)/((s2 - s3)*(s2 - s1)) &
                     + r3(:)*(sh - s1)*(sh - s2)/((s3 - s1)*(s3 - s2))
             END DO
          END IF
       END DO

       !===Create volume dofs
       r1 = p1_mesh%rr(:,p1_mesh%jj(1,m))
       r2 = p1_mesh%rr(:,p1_mesh%jj(2,m))
       r3 = p1_mesh%rr(:,p1_mesh%jj(3,m))
       l = v_dof + 3*f_dof !===Initialize cell dofs
       lmax = type_fe-2
       DO h = 1, lmax
          sh = h/DBLE(type_fe) !===move along y-axis
          DO k = 1, lmax-h+1
             sk = k/DBLE(type_fe) !===move along x-axis
             n_dof = n_dof + 1
             p3_mesh%jj(p3_mesh%gauss%n_w,m) = n_dof !===Volume dof
             p3_mesh%rr(:,n_dof) = (1-sh-sk)*r1+sk*r2+sh*r3
          END DO
       END DO
    END DO

    !===Create surface connectivity
    DO ms = 1, p1_mesh%mes
       m_op  = p1_mesh%neighs(ms)
       j1 = p1_mesh%jjs(1,ms)
       j2 = p1_mesh%jjs(2,ms)
       DO n_op = 1, 3
          IF(j1==p1_mesh%jj(n_op,m_op) .OR. j2==p1_mesh%jj(n_op,m_op)) CYCLE
          !===n is the index of the face that is on the boundary
          p3_mesh%jjs(1:2,ms) = p1_mesh%jjs(1:2,ms) !===Vertex dofs
          n1_op = MODULO(n_op+1,3)+1
          n2_op = MODULO(n_op+2,3)+1
          j1_op = p1_mesh%jj(n1_op,m_op)
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
             hh_op = h*sgn_op + ((1-sgn_op)/2)*(f_dof+1)
             p3_mesh%jjs(2+h,ms) = p3_mesh%jj(v_dof+(n_op-1)*f_dof+hh_op, m_op) !===Face dofs
          END DO
       END DO
    END DO

    !===Take care of leftover
    ALLOCATE(p3_mesh%iis(p3_mesh%gauss%n_ws,p3_mesh%mes))
    CALL dirichlet_nodes(p3_mesh%jjs, SPREAD(1,1,p3_mesh%mes), SPREAD(.TRUE.,1,1), p3_mesh%j_s)
    CALL surf_nodes_i(p3_mesh%jjs, p3_mesh%j_s, p3_mesh%iis)
    p3_mesh%nps = SIZE(p3_mesh%j_s)

  END SUBROUTINE create_p3_mesh
  !===

  SUBROUTINE incr_vrtx_indx_enumeration(mesh,type_fe)
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
       nmin = minloc(mesh%jj(1:3,m))
       nmax = maxloc(mesh%jj(1:3,m))
       n = 6/(nmin(1)*nmax(1))
       jjive(1) = mesh%jj(nmin(1),m)
       jjive(2) = mesh%jj(n,m)
       jjive(3) = mesh%jj(nmax(1),m)
       mesh%jj(1:3,m) = jjive
       !===Do it also for p2 (This is awkward. It should have been done well before.)
       IF (type_fe==2) THEN
          jjive(1) = mesh%jj(nmin(1)+3,m)
          jjive(2) = mesh%jj(n+3,m)
          jjive(3) = mesh%jj(nmax(1)+3,m)
          mesh%jj(4:6,m) = jjive
       END IF
       !===Do it also for mesh%neigh
       jjive(1) = mesh%neigh(nmin(1),m)
       jjive(2) = mesh%neigh(n,m)
       jjive(3) = mesh%neigh(nmax(1),m)
       mesh%neigh(:,m) = jjive
    END DO
    DO ms = 1, mesh%mes
       nmin = minloc(mesh%jjs(1:2,ms))
       nmax = maxloc(mesh%jjs(1:2,ms))
       jjive(1) = mesh%jjs(nmin(1),ms)
       jjive(2) = mesh%jjs(nmax(1),ms)
       mesh%jjs(1:2,ms) = jjive(1:2)
       !===For p2 only (nws=3)
    END DO
    !===JLG July 20, 2019, p3 mesh
  END SUBROUTINE incr_vrtx_indx_enumeration

  SUBROUTINE incr_vrtx_indx_enumeration_for_interfaces(interface_XX,nws1,nws2)
    USE def_type_mesh
    INTEGER, INTENT(IN)   :: nws1, nws2
    TYPE(interface_type)  :: interface_XX
    INTEGER :: js1ive(nws1), js2ive(nws2)
    INTEGER :: ms
    DO ms = 1, interface_XX%mes
       IF (interface_XX%jjs1(2,ms)<interface_XX%jjs1(1,ms)) THEN
          js1ive = interface_XX%jjs1(:,ms)
          interface_XX%jjs1(2,ms) = js1ive(1)
          interface_XX%jjs1(1,ms) = js1ive(2)
          IF(nws1.GE.3) interface_XX%jjs1(nws1:3:-1,ms)=js1ive(3:nws1)
       END IF
       IF (interface_XX%jjs2(2,ms)<interface_XX%jjs2(1,ms)) THEN
          js2ive = interface_XX%jjs2(:,ms)
          interface_XX%jjs2(2,ms) = js2ive(1)
          interface_XX%jjs2(1,ms) = js2ive(2)
          IF(nws2.GE.3) interface_XX%jjs2(nws2:3:-1,ms)=js2ive(3:nws2)
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
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom, list_inter
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect, stat
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_ms
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect
    LOGICAL :: t1, t2
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, mop, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=40) :: text
    CHARACTER(len=2)  :: truc

    mesh%gauss%k_d = 2 !===Space dimension (meridian section)

    text = 'Mesh'
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    !WRITE (*,*) 'Loading mesh-file ...'
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------
    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
          DO m = 1, me
             READ(30,*)
          END DO
          DO ms = 1, mes
             READ(30,*)
          END DO
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30) np, nw, me, nws, mes
          READ(30)
          READ(30)
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ(30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
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
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))
    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
       END DO
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
       END DO
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) jj_lect, neigh_lect, i_d_lect
       READ(30) jjs_lect, neighs_lect, sides_lect
       READ(30) rr_lect
    END IF

    !---Renumerotation------------------------------------------------------------
    ! Identify the status of faces
    ! stat = 1 (interface to be forgotten), stat = 2 (boundary), stat = 3 (real interface)
    ALLOCATE (stat(mes))
    DO ms = 1, mes
       neighs1 = neighs_lect(ms)
       IF (neighs1==0) THEN
          WRITE(*,*) ' BUG in prep_mesh, neighs1=0 '
          STOP
       END IF
       IF (MINVAL(ABS(i_d_lect(neighs1) - list_dom))==0) THEN
          t1 = .TRUE.
       ELSE
          t1 = .FALSE.
       END IF
       DO n = 1, nw
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
       END DO
       neighs2 = neigh_lect(n,neighs1)
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
                stat(ms) = 1 ! no inteface to treat
             ELSE IF (MINVAL(ABS(sides_lect(ms)-list_inter))==0) THEN
                stat(ms) = 3 ! real interface
             ELSE
                stat(ms) = 1 ! interface to be forgotten
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

    ALLOCATE (nouv_nd(np),  nouv_els(mes), virgin_nd(np), virgin_ms(mes), nouv_el(me))
    nouv_nd = -1000
    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    mnouv = 0
    msnouv= 0
    nnouv = 0
    DO dom = 1, SIZE(list_dom)
       DO m = 1, me ! Count new nodes from domain: i_d=dom
          IF (list_dom(dom) /= i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          mnouv = mnouv + 1  ! Nouvel element
          nouv_el(m) = mnouv
          DO n = 1, nw; i = jj_lect(n,m)
             IF (virgin_nd(i)) THEN ! Nouveau point
                virgin_nd(i) = .FALSE.
                nnouv = nnouv + 1
             END IF
          END DO
       END DO

       DO ms = 1, mes
          IF (stat(ms) /= 1 .AND. stat(ms) /=2 .AND. stat(ms) /=3) THEN
             WRITE(*,*) ' BUG in prep_mesh, stat out of bounds '
             STOP
          END IF

          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
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
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
             virgin_ms(ms) = .TRUE.
          END IF
       END DO
    END DO
    mesh%me  = mnouv
    mesh%np  = nnouv
    mesh%mes = msnouv

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))
    ALLOCATE(mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes),  mesh%sides(mesh%mes))
    ALLOCATE(mesh%rr(kd,mesh%np))

    virgin_nd = .TRUE.
    virgin_ms = .TRUE.
    nnouv = 0
    msnouv = 0
    DO dom = 1, SIZE(list_dom)
       !Loop on me and get nouv_el and nouv_nd right
       DO m = 1, me
          IF (list_dom(dom) /=i_d_lect(m))  CYCLE ! i_d(m) pas dans la liste
          DO n = 1, nw; i = jj_lect(n,m)
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
          DO n = 1, nw; i = jj_lect(n,m)
             IF (n .LE. nwneigh) THEN
                mop = neigh_lect(n,m)
                IF (mop .LE. 0) THEN
                   mesh%neigh(n,nouv_el(m)) = 0
                ELSE IF (MINVAL(ABS(list_dom - i_d_lect(mop))) == 0) THEN
                   mesh%neigh(n,nouv_el(m)) = nouv_el(mop)
                ELSE
                   mesh%neigh(n,nouv_el(m)) = 0
                END IF
             END IF
             mesh%rr(:,nouv_nd(i)) = rr_lect(:,i)
          END DO
          mesh%i_d(nouv_el(m)) = i_d_lect(m)
          mesh%jj(:,nouv_el(m)) = nouv_nd(jj_lect(:,m))
       END DO

       !Loop on mes and get neighs_lect and nouv_els right
       DO ms = 1, mes
          !I test if ms touches the current domain of interest: i_d = dom
          IF (stat(ms) == 1) CYCLE
          neighs1 = neighs_lect(ms)
          DO n = 1, nw
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
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
             virgin_nd(jjs_lect(:,ms)) = .TRUE. ! interface nodes are virgin again
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
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
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
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
             ! exit when n is not on the interface
          END DO
          neighs2 = neigh_lect(n,neighs1)
          IF (i_d_lect(neighs1) /= list_dom(dom)) THEN
             IF (neighs2 == 0) CYCLE ! face on the boundary and does not touch dom yet
             IF (i_d_lect(neighs2) /= list_dom(dom)) CYCLE ! dom is on neither sides
          END IF
          !End test if ms touches the domain of interest

          mesh%jjs(:,nouv_els(ms))  = nouv_nd(jjs_lect(:,ms))
          mesh%neighs(nouv_els(ms)) = nouv_el(neighs_lect(ms))
          mesh%sides(nouv_els(ms))  = sides_lect(ms)
          IF (stat(ms)==3) THEN ! side is an interface to be kept
             neighs1 = neighs_lect(ms)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs1)) = 0
             neighs2 = neigh_lect(n,neighs1)
             DO n = 1, nw
                IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs2))) /= 0) EXIT ! n not on the interface
             END DO
             mesh%neigh(n,nouv_el(neighs2)) = 0
          END IF
       END DO

    END DO
    !===End reordring

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_nd, virgin_ms, stat)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps
    WRITE (20,*) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
    CLOSE(20)
    CLOSE(30)

    mesh%edge_stab = .FALSE.

  END SUBROUTINE load_dg_mesh_free_format

  !------------------------------------------------------------------------------

  SUBROUTINE load_mesh_free_format_ordered(dir, fil, list_dom, type_fe, &
       mesh, mesh_formatted, edge_stab)

    USE def_type_mesh
    USE chaine_caractere
    USE Dir_nodes

    IMPLICIT NONE
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    LOGICAL, OPTIONAL,     INTENT(IN) :: edge_stab

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i, dom
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=60)                 :: text
    CHARACTER(len=2)                  :: truc
    text = 'Mesh'
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF
    !WRITE (*,*) 'Loading mesh-file ...'
    !d_end = last_c_leng (64, dir)
    !f_end = last_c_leng (64, fil)
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !     READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN    ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
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
       WRITE(*,*) ' Finite element not yet programmed '
       WRITE(*,*) kd, nw, nws
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
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
          IF (ABS(list_dom(dom)-i_d_lect(m)) /= 0)  CYCLE ! i_d(m) dans la liste
          virgin_el(m) = .FALSE.
          mnouv = mnouv + 1   ! Nouvel element
          nouv_el(m) = mnouv
          ancien_el(mnouv) = m
          DO n = 1, nw; i = jj_lect(n,m)
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

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---  Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
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
             IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
          END DO
          neighs2 = neigh_lect(n,neighs1)
          test = .FALSE.
          IF (ABS(list_dom(dom)-i_d_lect(neighs1)) == 0) THEN
             test=.TRUE.
          ELSE IF (neighs2 /= 0) THEN
             IF (ABS(list_dom(dom)-i_d_lect(neighs2)) == 0) THEN
                test=.TRUE.
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

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---  Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !     END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps
    WRITE (20,*) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
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
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh
    LOGICAL,               INTENT(IN) :: mesh_formatted !formatted <=> mesh_formatted=.true.
    LOGICAL, OPTIONAL,     INTENT(IN) :: edge_stab

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=60)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    !WRITE (*,*) 'Loading mesh-file ...'
    !d_end = last_c_leng (64, dir)
    !f_end = last_c_leng (64, fil)
    IF (mesh_formatted) THEN
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    ELSE
       OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    END IF
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       IF (mesh_formatted) THEN
          READ(30,*) np, nw, me, nws, mes
       ELSE
          READ(30) np, nw, me, nws, mes
       END IF

       IF (mesh_formatted) THEN
          DO m = 1, me
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO ms = 1, mes
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF

       IF (mesh_formatted) THEN
          DO n = 1, np
             READ(30,*)
          END DO
       ELSE
          READ(30)
       END IF
    END IF

    IF (mesh_formatted) THEN
       READ  (30, *)  np,  nw,  me,  nws,  mes
    ELSE
       READ(30)  np,  nw,  me,  nws,  mes
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
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    IF (mesh_formatted) THEN
       DO m = 1, me
          READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
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
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
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

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    IF (mesh_formatted) THEN
       DO ms = 1, mes
          READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
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
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
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

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    IF (mesh_formatted) THEN
       DO n = 1, np
          READ(30,*) rr_lect(:,n)
       END DO
    ELSE
       READ(30) rr_lect
    END IF

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps
    WRITE (20,*) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
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
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh

    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)

    !WRITE (*,*) 'Loading mesh-file ...'

    OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='formatted')
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='unformatted')
    OPEN(UNIT=20,FILE='error_mesh', FORM='formatted',STATUS='unknown')

    !  READ GRID DATA AND ARRAY ALLOCATION ----------------------------------------

    IF (type_fe == 2) THEN ! Skip P1 data if needed
       READ(30,*) np, nw, me, nws, mes
       !READ(30) np, nw, me, nws, mes
       DO m = 1, me
          READ(30,*)
       END DO
       !READ(30)
       DO ms = 1, mes
          READ(30,*)
       END DO
       !READ(30)
       DO n = 1, np
          READ(30,*)
       END DO
       !READ(30)
    END IF

    READ  (30, *)  np,  nw,  me,  nws,  mes
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
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
         nouv_el(0:me), ancien_el(me), virgin_el(me))

    nouv_el = 0

    DO m = 1, me
       READ(30,*) jj_lect(:,m), neigh_lect(:,m), i_d_lect(m)
    END DO

    !READ(30) jj_lect, neigh_lect, i_d_lect

    !---Renumerotation------------------------------------------------------------
    virgin_nd = .TRUE.
    virgin_el = .TRUE.
    mnouv = 0
    nnouv = 0

    DO m = 1, me
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
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

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

    DO ms = 1, mes
       READ(30,*) jjs_lect(:,ms), neighs_lect(ms), sides_lect(ms)
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
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
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

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    DO n = 1, np
       READ(30,*) rr_lect(:,n)
    END DO
    !READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps
    WRITE (20,*) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
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
    CHARACTER(len=200),     INTENT(IN) :: dir, fil
    INTEGER, DIMENSION(:), INTENT(IN) :: list_dom
    INTEGER,               INTENT(IN) :: type_fe
    TYPE(mesh_type)                   :: mesh

    INTEGER, ALLOCATABLE, DIMENSION(:)   :: nouv_nd, nouv_el, nouv_els, &
         ancien_nd, ancien_el, ancien_els
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: jj_lect, neigh_lect, jjs_lect
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: i_d_lect, sides_lect, neighs_lect
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: virgin_nd, virgin_el
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: rr_lect

    LOGICAL :: test
    INTEGER :: mnouv, nnouv, i
    INTEGER :: n, m, ms, msnouv, neighs1, neighs2
    INTEGER :: d_end, f_end
    INTEGER :: np, nw, me, nws, mes, kd, nwneigh
    CHARACTER(len=20)                 :: text
    CHARACTER(len=2)                  :: truc

    text = 'Mesh'
    d_end = last_c_leng (20, text)
    DO n = 1, SIZE(list_dom)
       d_end = last_c_leng (20, text)
       WRITE(truc,'(i2)') list_dom(n)
       f_end = start_of_string (truc)
       text = text(1:d_end)//'_'//truc(f_end:)
    END DO

    d_end = last_c_leng (20, text)
    IF (type_fe==1) THEN
       text = text(1:d_end)//'_FE_1'
    ELSE
       text = text(1:d_end)//'_FE_2'
    END IF

    !WRITE (*,*) 'Loading mesh-file ...'
    d_end = last_c_leng (64, dir)
    f_end = last_c_leng (64, fil)
    !OPEN(30,FILE=dir(1:d_end)//'/'//fil(1:f_end),FORM='formatted')
    OPEN(30,FILE=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fil)),FORM='unformatted')
    OPEN(UNIT=20,FILE=text, FORM='formatted', STATUS='unknown')

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
    READ(30)  np,  nw,  me,  nws,  mes

    IF (nw==3 .AND. nws==2) THEN ! Decide about space dimension
       kd = 2; nwneigh = 3
    ELSE IF (nw==6 .AND. nws==3) THEN
       kd = 2; nwneigh = 3
    ELSE IF (nw==4 .AND. nws==3) THEN
       kd = 3; nwneigh = 4
    ELSE IF (nw==10 .AND. nws==6) THEN
       kd = 3; nwneigh = 4
    ELSE
       WRITE(*,*) ' Finite element not yet programmed '
       STOP
    END IF

    ALLOCATE (jj_lect(nw,me),neigh_lect(nwneigh,me),i_d_lect(me))
    ALLOCATE (nouv_nd(np),   ancien_nd(np), virgin_nd(np), &
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
       IF (MINVAL(ABS(list_dom-i_d_lect(m))) /= 0)  CYCLE ! i_d(m) dans la liste
       virgin_el(m) = .FALSE.
       mnouv = mnouv + 1  ! Nouvel element
       nouv_el(m) = mnouv
       ancien_el(mnouv) = m
       DO n = 1, nw; i = jj_lect(n,m)
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

    ALLOCATE(mesh%jj(nw,mesh%me), mesh%neigh(nwneigh,mesh%me), mesh%i_d(mesh%me))

    DO m = 1, mesh%me
       mesh%jj(:,m) = nouv_nd(jj_lect(:,ancien_el(m)))
       mesh%neigh(:,m) = nouv_el(neigh_lect(:,ancien_el(m)))
       mesh%i_d(m) = i_d_lect(ancien_el(m))
    END DO

    !---Fin renumerotation--------------------------------------------------------
    ALLOCATE (jjs_lect(nws,mes), neighs_lect(mes), sides_lect(mes))

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
          IF (MINVAL(ABS(jjs_lect(:,ms)-jj_lect(n,neighs1))) /= 0) EXIT
       END DO
       neighs2 = neigh_lect(n,neighs1)
       test = .FALSE.
       IF (MINVAL(ABS(list_dom-i_d_lect(neighs1))) == 0) THEN
          test=.TRUE.
       ELSE IF (neighs2 /= 0) THEN
          IF (MINVAL(ABS(list_dom-i_d_lect(neighs2))) == 0) THEN
             test=.TRUE.
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

    ALLOCATE (mesh%jjs(nws,mesh%mes), mesh%neighs(mesh%mes), &
         mesh%sides(mesh%mes))

    DO ms = 1, mesh%mes
       mesh%jjs(:,ms) = nouv_nd(jjs_lect(:,ancien_els(ms)))
       mesh%neighs(ms) = nouv_el(neighs_lect(ancien_els(ms)))
       mesh%sides(ms)  = sides_lect(ancien_els(ms))
    END DO

    !---Fin renumerotation--------------------------------------------------------

    ALLOCATE(rr_lect(kd,np))

    !DO n = 1, np
    !   READ(30,*) rr_lect(:,n)
    !END DO

    READ(30) rr_lect

    ALLOCATE(mesh%rr(kd,mesh%np))

    mesh%rr = rr_lect(:,ancien_nd(1:mesh%np))

    DEALLOCATE(jj_lect, neigh_lect, i_d_lect)
    DEALLOCATE(jjs_lect, neighs_lect, sides_lect)
    DEALLOCATE(rr_lect, virgin_el, virgin_nd)
    DEALLOCATE(nouv_nd,nouv_el,nouv_els,ancien_nd,ancien_el,ancien_els)
    !  END OF GRID READING --------------------------------------------------------

    ALLOCATE(mesh%iis(nws,mesh%mes))
    CALL dirichlet_nodes(mesh%jjs, SPREAD(1,1,mesh%mes), SPREAD(.TRUE.,1,1), mesh%j_s)
    CALL surf_nodes_i(mesh%jjs, mesh%j_s,  mesh%iis)
    mesh%nps = SIZE(mesh%j_s)

    WRITE (20,*)  'np_lect',np,'nw_lect',nw,'nws_lect',nws,'me_lect',me,&
         'mes_lect',mes
    WRITE (20,*)  'np ', mesh%np, 'me ',  mesh%me, &
         'mes ',mesh%mes,'nps ', mesh%nps
    WRITE (20,*) 'MAXVAL(sides)', MAXVAL(mesh%sides), 'MAXVAL(i_d)', MAXVAL(mesh%i_d)
    CLOSE(20)
    CLOSE(30)

  END SUBROUTINE load_mesh


  SUBROUTINE surf_nodes_i(jjs, j_s,  iis)

    !  generation of the surface element connectivity matrix  iis
    !  based on the surface node numbering, starting from the
    !  connectivity matrix  jjs  of the surface elements according
    !  to the volume node numbering, and from the array  j_s  of
    !  the boundary nodes according to the volume node numbering

    IMPLICIT NONE

    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jjs
    INTEGER, DIMENSION(:),   INTENT(IN)  :: j_s
    INTEGER, DIMENSION(:,:), INTENT(OUT) :: iis

    INTEGER :: ms, ls, j, i

    DO ms = 1, SIZE(jjs,2)
       DO ls = 1, SIZE(jjs,1)
          j = jjs(ls,ms)
          DO i = 1, SIZE(j_s)
             IF ( j == j_s(i) )  iis(ls,ms) = i
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE surf_nodes_i

  SUBROUTINE prep_interfaces(mesh)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                   :: mesh
    LOGICAL, DIMENSION(mesh%me) :: virgin
    INTEGER :: m, mop, nw, me, l, lop, n, n1, n2, edge, nt,nws
    nw = SIZE(mesh%jj,1)
    nws = SIZE(mesh%jjs,1)
    me = mesh%me
    IF (SIZE(mesh%rr,1)==2) THEN
       nt=3
    ELSE
       WRITE(*,*) ' BUG: prep_interfaces, 3D not programmed yet '
       STOP
       nt=4
    END IF

    virgin = .TRUE.
    edge = 0
    DO m = 1, me
       virgin(m) = .FALSE.
       DO n = 1, nt
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE !Edge on boundary
          IF (.NOT.virgin(mop)) CYCLE !Edge already done
          edge = edge + 1 !New edge
       END DO
    END DO
    IF (SIZE(mesh%rr,1)==2) THEN
       IF (edge/=(3*mesh%me - mesh%mes)/2) THEN
          WRITE(*,*) ' BUG in prep_interfaces, edge/=(3*mesh%me + mesh%mes)/2'
          WRITE(*,*) ' edge ', edge, (3*mesh%me - mesh%mes)/2
          WRITE(*,*) ' mesh%mes ', mesh%mes, ' mesh%me ',mesh%me
          STOP
       END IF
    END IF

    mesh%mi = edge
    ALLOCATE(mesh%jji(nw,2,edge), mesh%jjsi(nws,edge),mesh%neighi(2,edge) )

    edge = 0
    virgin = .TRUE.
    DO m = 1, me
       virgin(m) = .FALSE.
       DO n = 1, nt
          mop = mesh%neigh(n,m)
          IF (mop==0) CYCLE !Edge on boundary
          IF (.NOT.virgin(mop)) CYCLE !Edge already done
          edge = edge + 1 !New edge
          n1 = MODULO(n,nt) + 1
          n2 = MODULO(n+1,nt) + 1 ! Works in 2D only
          mesh%jjsi(1,edge) = mesh%jj(n1,m)
          mesh%jjsi(2,edge) = mesh%jj(n2,m)
          IF (nws==3) THEN
             IF (n1+n2==3) THEN
                mesh%jjsi(3,edge) = mesh%jj(6,m)
             ELSE IF (n1+n2==5) THEN
                mesh%jjsi(3,edge) = mesh%jj(4,m)
             ELSE IF (n1+n2==4) THEN
                mesh%jjsi(3,edge) = mesh%jj(5,m)
             ELSE
                WRITE(*,*) ' BUG prep_interfaces n1+n2 not correct '
                STOP
             END IF
          END IF
          IF (m<mop) THEN
             l = 1 ! Side 1 on smallest cell index
             lop = 2
             mesh%neighi(1,edge) = m
             mesh%neighi(2,edge) = mop
          ELSE
             l = 2 ! Side 2 on largest cell index
             lop = 1
             mesh%neighi(1,edge) = mop
             mesh%neighi(2,edge) = m
          END IF
          DO n1 = 1, nw
             mesh%jji(n1,l,edge) = mesh%jj(n1,m)
             mesh%jji(n1,lop,edge) = mesh%jj(n1,mop)
          END DO
       END DO
    END DO

  END SUBROUTINE prep_interfaces

END MODULE prep_maill
