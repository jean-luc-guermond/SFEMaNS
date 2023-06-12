!
!Authors Jean-Luc Guermond, Copyrights 2000
!
MODULE prep_mesh_interface

  IMPLICIT NONE

  PUBLIC :: load_interface, load_mesh_interface
  PRIVATE

CONTAINS

  !------------------------------------------------------------------------------


  SUBROUTINE load_interface(mesh_master, mesh_slave, list_inter, mesh_INTERFACE, disjoint)

    USE def_type_mesh

    TYPE(mesh_type),           INTENT(IN)  :: mesh_master, mesh_slave
    INTEGER, DIMENSION(:),     INTENT(IN)  :: list_inter
    TYPE(interface_type),      INTENT(OUT) ::  mesh_INTERFACE
    LOGICAL,                   INTENT(IN)  :: disjoint

    INTEGER                              :: dim, ms, ms1, ms2, ns, k, nws_master, nws_slave, m1, m2
    INTEGER, DIMENSION(:),   ALLOCATABLE :: list, interface_mesh1, interface_mesh2
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: interface_jjs1, interface_jjs2
    REAL(KIND=8)                         :: eps_ref=1.d-7, r_norm, epsilon
    LOGICAL                              :: okay
    LOGICAL, DIMENSION(:),   ALLOCATABLE :: virgin_elem

    ! Removing gauss, FL, Mar. 23
!!$    nws_master = mesh_master%gauss%n_ws
!!$    nws_slave  = mesh_slave%gauss%n_ws

    nws_master = SIZE(mesh_master%jjs,1)
    nws_slave  = SIZE(mesh_slave%jjs,1)
    ! Done removing gauss, FL, Mar. 23

    IF (nws_master > nws_slave) THEN
       WRITE(*,*) ' BUG in load_interface: nws_master > nws_slave '
       STOP
    END IF

    dim = SIZE(mesh_master%rr,1)
    IF (dim>2) THEN
       WRITE(*,*) ' Dimension 3 not yet programmed '
       STOP
    END IF

    ALLOCATE(virgin_elem(mesh_slave%me), list(dim), interface_mesh1(mesh_master%mes), &
         interface_mesh2(mesh_slave%mes), interface_jjs1(nws_master,mesh_master%mes), &
         interface_jjs2(nws_slave,mesh_slave%mes) )

    virgin_elem = .TRUE.

    ms = 0
    DO ms1 = 1, mesh_master%mes
       IF(MINVAL(ABS(list_inter-mesh_master%sides(ms1))) /= 0) CYCLE !not on interface
       r_norm = SUM(ABS(mesh_master%rr(:,mesh_master%jjs(1,ms1))- &
            &    mesh_master%rr(:,mesh_master%jjs(2,ms1))))
       epsilon = eps_ref*r_norm
       okay = .FALSE.

       lp2: DO ms2 = 1, mesh_slave%mes
          !IF(.NOT.virgin_elem(ms2)) CYCLE !element done!
          IF(MINVAL(ABS(list_inter-mesh_slave%sides(ms2))) /= 0) CYCLE !not on interface


          DO k = 0, dim-1 !dim=2
             DO ns = 1, dim !dim=2
                list(ns) = MODULO(ns-1+k,dim) + 1
             END DO

             IF (MAXVAL(ABS(mesh_master%rr(:,mesh_master%jjs(list,ms1))&
                  -mesh_slave%rr(:,mesh_slave%jjs(1:dim,ms2)))).GT.epsilon) CYCLE

             IF(.NOT.virgin_elem(ms2)) THEN
                Okay = .TRUE.
                CYCLE !already element done
             ENDIF

             m1  = mesh_master%neighs(ms1)
             m2  = mesh_slave%neighs(ms2)
             r_norm = SUM(ABS(mesh_master%rr(:,mesh_master%jj(1:3,m1)) &
                  - mesh_slave%rr(:,mesh_slave%jj(1:3,m2))))
             IF ( r_norm .LE. 1d-9) THEN
                CYCLE  ! two identical triangles
             END IF

             ms = ms + 1
             interface_mesh1(ms) = ms1
             interface_mesh2(ms) = ms2
             interface_jjs1(1:2,ms) = mesh_master%jjs(list,ms1)
             interface_jjs2(:,ms) = mesh_slave%jjs(:,ms2)
             IF (nws_master==3) THEN !P2 in two dimensions
                interface_jjs1(3,ms) = mesh_master%jjs(3,ms1)
             END IF

             virgin_elem(ms2) =.FALSE.
             IF (.NOT.disjoint) virgin_elem(ms1) =.FALSE.
             okay = .TRUE.

             EXIT lp2
          END DO
       END DO lp2
       IF (.NOT.okay) THEN
          WRITE(*,*) ' BUG in load_interface: .NOT.okay'
          STOP
       END IF
    END DO

    mesh_interface%mes = ms

    ALLOCATE( mesh_interface%mesh1(ms), mesh_interface%mesh2(ms),&
         mesh_interface%jjs1(nws_master,ms), mesh_interface%jjs2(nws_slave,ms))


    mesh_interface%mesh1 = interface_mesh1(1:ms)
    mesh_interface%mesh2 = interface_mesh2(1:ms)
    mesh_interface%jjs1 = interface_jjs1(1:nws_master,1:ms)
    mesh_interface%jjs2 = interface_jjs2(1:nws_slave,1:ms)


    DEALLOCATE(virgin_elem, list, interface_mesh1, interface_mesh2, &
         interface_jjs1,interface_jjs2)

  END SUBROUTINE load_interface

  !-----------------------------------------------------------------------

  SUBROUTINE load_mesh_interface(mesh_master, mesh_slave, list_inter,  mesh_INTERFACE)

    USE def_type_mesh

    TYPE(mesh_type),           INTENT(IN)  :: mesh_master, mesh_slave
    INTEGER, DIMENSION(:),     INTENT(IN)  :: list_inter
    TYPE(mesh_type_interface), INTENT(OUT) ::  mesh_INTERFACE

    INTEGER                              :: dim, ms1, ms2, ns, k, m, n, nn, &
         nws_master, nws_slave, nw_slave, n_count
    INTEGER, DIMENSION(:),   ALLOCATABLE :: list, interface_slave_elem_temp, &
         node_master
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: interface_master_node_temp
    REAL(KIND=8)                         :: eps_ref=1.d-7, r_norm, epsilon
    LOGICAL                              :: okay
    LOGICAL, DIMENSION(:),   ALLOCATABLE :: virgin_elem, virgin_node, virgin_node_inter

    nws_master = mesh_master%gauss%n_ws
    nws_slave  = mesh_slave%gauss%n_ws
    nw_slave  = mesh_slave%gauss%n_w

    IF (nws_master > nws_slave) THEN
       WRITE(*,*) ' BUG in load_mesh_interface: nws_master > nws_slave'
       STOP
    END IF

    dim = SIZE(mesh_master%rr,1)
    IF (dim>2) THEN
       WRITE(*,*) ' Dimension 3 not yet programmed '
       STOP
    END IF

    ALLOCATE(virgin_elem(mesh_slave%me), list(dim), virgin_node(mesh_slave%np), &
         node_master(mesh_slave%np), interface_slave_elem_temp(mesh_slave%me),&
         interface_master_node_temp(nw_slave,mesh_slave%me), virgin_node_inter(mesh_slave%nps))

    !==List of slave nodes on interface
    virgin_node_inter = .TRUE.
    DO ms2 = 1, mesh_slave%mes
       IF(MINVAL(ABS(list_inter-mesh_slave%sides(ms2))) /= 0) CYCLE !not on interface
       virgin_node_inter(mesh_slave%iis(:,ms2)) = .FALSE.
    END DO
    n_count = 0
    DO n = 1, mesh_slave%nps
       IF (virgin_node_inter(n)) CYCLE
       n_count = n_count + 1
    END DO
    ALLOCATE( mesh_interface%list_slave_node(n_count))
    n_count = 0
    DO n = 1, mesh_slave%nps
       IF (virgin_node_inter(n)) CYCLE
       n_count = n_count + 1
       mesh_interface%list_slave_node(n_count) = mesh_slave%j_s(n)
    END DO
    !==End list of slave nodes on interface

    virgin_elem = .TRUE.
    virgin_node = .TRUE.

    DO ms1 = 1, mesh_master%mes
       IF(MINVAL(ABS(list_inter-mesh_master%sides(ms1))) /= 0) CYCLE !not on interface
       r_norm = SUM(ABS(mesh_master%rr(:,mesh_master%jjs(1,ms1))- &
            &    mesh_master%rr(:,mesh_master%jjs(2,ms1))))
       epsilon = eps_ref*r_norm
       okay = .FALSE.
       lp2: DO ms2 = 1, mesh_slave%mes
          IF(.NOT.virgin_elem(ms2)) CYCLE !element done
          IF(MINVAL(ABS(list_inter-mesh_slave%sides(ms2))) /= 0) CYCLE !not on interface

          DO k = 0, dim-1 !dim=2
             DO ns = 1, dim !dim=2
                list(ns) = MODULO(ns-1+k,dim) + 1
             END DO

             IF (MAXVAL(ABS(mesh_master%rr(:,mesh_master%jjs(list,ms1))&
                  -mesh_slave%rr(:,mesh_slave%jjs(1:dim,ms2)))).GT.epsilon) CYCLE

             node_master(mesh_slave%jjs(:,ms2)) = mesh_master%jjs(list,ms1)
             IF (nws_master==3) THEN !P2 in two dimensions
                node_master(mesh_slave%jjs(3,ms2)) = mesh_master%jjs(3,ms1)
             END IF
             virgin_node(mesh_slave%jjs(1:nws_master,ms2)) = .FALSE.

             virgin_elem(ms2) =.FALSE.
             okay = .TRUE.

             EXIT lp2
          END DO
       END DO lp2
       IF (.NOT.okay) THEN
          WRITE(*,*) ' BUG in load_mesh_interface: .NOT.okay'
          STOP
       END IF
    END DO

    interface_master_node_temp = -1
    mesh_interface%me = 0
    DO m = 1, mesh_slave%me
       nn = 0
       DO n = 1, nw_slave
          IF(virgin_node(mesh_slave%jj(n,m))) CYCLE
          nn = nn + 1
          interface_master_node_temp(n,m) = node_master(mesh_slave%jj(n,m))
       END DO
       IF (nn/=0) THEN
          mesh_interface%me =  mesh_interface%me + 1
          interface_slave_elem_temp( mesh_interface%me) = m

       END IF
    END DO

    ALLOCATE( mesh_interface%slave_elem( mesh_interface%me))
    ALLOCATE( mesh_interface%master_node(nw_slave, mesh_interface%me), &
         mesh_interface%slave_node (nw_slave, mesh_interface%me))

    mesh_interface%slave_elem = interface_slave_elem_temp(1: mesh_interface%me)

    mesh_interface%slave_node = mesh_slave%jj(:, mesh_interface%slave_elem)

    mesh_interface%master_node = interface_master_node_temp(:, mesh_interface%slave_elem)

    DEALLOCATE(virgin_elem, list, virgin_node, virgin_node_inter, &
         node_master, interface_slave_elem_temp,&
         interface_master_node_temp)

  END SUBROUTINE load_mesh_interface

END MODULE prep_mesh_interface
