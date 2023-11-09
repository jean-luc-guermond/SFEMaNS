MODULE metis_sfemans
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC :: part_mesh_mhd, free_mesh, extract_mesh, free_interface, part_mesh_mhd_bis, part_mesh_M_T_H_phi
  PRIVATE
  REAL(KIND=8) :: epsilon = 1.d-10
!!$ Dummy for metis...
  INTEGER      :: METIS_NOPTIONS=40, METIS_OPTION_NUMBERING=18
  !==Go see in metis.h the actual values assigned to METIS_NOPTIONS, METIS_OPTION_NUMBERING
!!$ Dummy for metis...
CONTAINS
  !===Convention: Domain NS \subset domain Temp \subset Domain H
  SUBROUTINE part_mesh_M_T_H_phi(nb_proc, list_conc, list_u_in, list_T_in, list_h_in, &
       list_phi ,mesh,list_of_interfaces,part,my_periodic)
    USE def_type_mesh
    USE my_util
    USE sub_plot
    USE periodic
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    INTEGER, DIMENSION(mesh%me)              :: part
    INTEGER, DIMENSION(:)                    :: list_of_interfaces
    INTEGER, DIMENSION(:)                    :: list_conc, list_u_in, list_T_in
    INTEGER, DIMENSION(:)                    :: list_h_in, list_phi
    TYPE(periodic_data), OPTIONAL            :: my_periodic

    LOGICAL, DIMENSION(mesh%mes)             :: virgins
    INTEGER, DIMENSION(3,mesh%me)            :: neigh_new
    INTEGER, DIMENSION(5)                    :: opts
    INTEGER, DIMENSION(SIZE(mesh%jjs,1))     :: i_loc
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: xadj_conc, xadj_u, xadj_T, xadj_h, xadj_phi
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: list_u, list_T, list_h
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: xind_conc, xind_u, xind_T, xind_h, xind_phi
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: vwgt, adjwgt
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: conc2glob, u2glob, T2glob, h2glob, phi2glob
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: part_conc, part_u, part_T, part_h, part_phi
    INTEGER, DIMENSION(1)                    :: jm_loc
    INTEGER, DIMENSION(mesh%np,3)            :: per_pts
    INTEGER, DIMENSION(mesh%me)              :: glob2loc
    INTEGER, DIMENSION(mesh%np)              :: indicator
    INTEGER, DIMENSION(3)                    :: j_loc
    INTEGER :: nb_neigh, edge, m, ms, n, nb, numflag, p, wgtflag, j, &
         ns, nws, msop, nsop, proc, iop, mop, s2, k, me_conc, me_u, me_T, me_h, me_phi, idm
    REAL(KIND=8)                             :: err
    LOGICAL               :: test
    !===(JLG) Feb 20, 2019. Petsc developpers decide to use REAL(KIND=4) to interface with metis
    !REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tpwgts
    !REAL(KIND=8), DIMENSION(1)               :: ubvec
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE  :: tpwgts
    REAL(KIND=4), DIMENSION(1)               :: ubvec
    REAL(KIND=4)                             :: one_K4=1.0
    !===(JLG)Feb 20, 2019.
    INTEGER, DIMENSION(METIS_NOPTIONS)       :: metis_opt
    PetscMPIInt    :: nb_proc
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    IF (nb_proc==1) THEN
       part = 1
       RETURN
    END IF
    glob2loc = 0

    !===Create list_u = list_u_in \ list_conc
    nb = 0
    DO j = 1, SIZE(list_u_in)
       IF (MINVAL(ABS(list_u_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
    END DO
    ALLOCATE(list_u(nb))
    nb = 0
    DO j = 1, SIZE(list_u_in)
       IF (MINVAL(ABS(list_u_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
       list_u(nb) = list_u_in(j)
    END DO
    !==Done with list_u

    !===Create list_T = list_T_in \ list_u_in
    nb = 0
    DO j = 1, SIZE(list_T_in)
       IF (MINVAL(ABS(list_T_in(j)-list_u_in))==0) CYCLE
       IF (MINVAL(ABS(list_T_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
    END DO
    ALLOCATE(list_T(nb))
    nb = 0
    DO j = 1, SIZE(list_T_in)
       IF (MINVAL(ABS(list_T_in(j)-list_u_in))==0) CYCLE
       IF (MINVAL(ABS(list_T_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
       list_T(nb) = list_T_in(j)
    END DO
    !==Done with list_T

    !===Create list_h = list_h_in \ list_T_in
    nb = 0
    DO j = 1, SIZE(list_h_in)
       IF (MINVAL(ABS(list_h_in(j)-list_T_in))==0) CYCLE
       IF (MINVAL(ABS(list_h_in(j)-list_u_in))==0) CYCLE
       IF (MINVAL(ABS(list_h_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
    END DO
    ALLOCATE(list_h(nb))
    nb = 0
    DO j = 1, SIZE(list_h_in)
       IF (MINVAL(ABS(list_h_in(j)-list_T_in))==0) CYCLE
       IF (MINVAL(ABS(list_h_in(j)-list_u_in))==0) CYCLE
       IF (MINVAL(ABS(list_h_in(j)-list_conc))==0) CYCLE
       nb = nb + 1
       list_h(nb) = list_h_in(j)
    END DO
    !==Done with list_h

    !===Create neigh_new for interfaces
    nws = SIZE( mesh%jjs,1)
    neigh_new = mesh%neigh
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_M_T_H_phi, .NOT.test ')
          END IF
          DO n = 1, 3
             IF (neigh_new(n,mesh%neighs(msop))==0) THEN
                neigh_new(n,mesh%neighs(msop)) = mesh%neighs(ms)
             END IF
             IF (neigh_new(n,mesh%neighs(ms))==0) THEN
                neigh_new(n,mesh%neighs(ms)) = mesh%neighs(msop)
             END IF
          END DO
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
       END DO
    END IF
    !===End Create neigh_new for interfaces

    !===Create neigh_new for periodic faces
    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) == 0) THEN
                jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
                s2 = my_periodic%list_periodic(2,jm_loc(1))
                test = .FALSE.
                DO msop = 1, mesh%mes
                   IF (mesh%sides(msop) /= s2) CYCLE

                   err = 0.d0
                   DO ns = 1, SIZE(my_periodic%vect_e,1)
                      err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                           +my_periodic%vect_e(ns,jm_loc(1))))
                   END DO

                   IF (err .LE. epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   CALL error_Petsc('BUG in part_mesh_M_T_H_phi, mop not found')
                END IF
                mop = mesh%neighs(msop)
                DO n = 1, 3
                   IF (neigh_new(n,m) == 0) THEN
                      neigh_new(n,m) = mop
                   END IF
                   IF (neigh_new(n,mop) == 0) THEN
                      neigh_new(n,mop) = m
                   END IF
                END DO
             END IF
          END DO
       END IF
    END IF
    !===End Create neigh_new for periodic faces

    !===Create glob2loc and conc2glob, u2glob, T2glob, h2glob, phi2glob
    me_conc = 0
    me_u = 0
    me_T = 0
    me_h = 0
    me_phi = 0
    DO m = 1, mesh%me
       idm = mesh%i_d(m)
       IF (MINVAL(ABS(idm-list_conc))==0) THEN
          me_conc = me_conc + 1
       ELSE IF (MINVAL(ABS(idm-list_u))==0) THEN
          me_u = me_u + 1
       ELSE IF (MINVAL(ABS(idm-list_T))==0) THEN
          me_T = me_T + 1
       ELSE IF (MINVAL(ABS(idm-list_h))==0) THEN
          me_h = me_h + 1
       ELSE IF (MINVAL(ABS(idm-list_phi))==0) THEN
          me_phi = me_phi + 1
       ELSE
          CALL error_Petsc('BUG in part_mesh_M_T_H_phi : element not in the mesh')
       END IF
    END DO
    ALLOCATE(conc2glob(me_conc), u2glob(me_u), T2glob(me_T), h2glob(me_h), phi2glob(me_phi))
    me_conc = 0
    me_u = 0
    me_T = 0
    me_h = 0
    me_phi = 0
    DO m = 1, mesh%me
       idm = mesh%i_d(m)
       IF (MINVAL(ABS(idm-list_conc))==0) THEN
          me_conc = me_conc + 1
          conc2glob(me_conc) = m
          glob2loc(m) = me_conc
       ELSE IF (MINVAL(ABS(idm-list_u))==0) THEN
          me_u = me_u + 1
          u2glob(me_u) = m
          glob2loc(m) = me_u
       ELSE IF (MINVAL(ABS(idm-list_T))==0) THEN
          me_T = me_T + 1
          T2glob(me_T) = m
          glob2loc(m) = me_T
       ELSE IF (MINVAL(ABS(idm-list_h))==0) THEN
          me_h = me_h + 1
          h2glob(me_h) = m
          glob2loc(m) = me_h
       ELSE IF (MINVAL(ABS(idm-list_phi))==0) THEN
          me_phi = me_phi + 1
          phi2glob(me_phi) = m
          glob2loc(m) = me_phi
       ELSE
          CALL error_Petsc('BUG in part_mesh_M_T_H_phi: element not in the mesh')
       END IF
    END DO
    !===End Create glob2loc and u2glob, T2glob, h2glob, phi2glob, conc2glob

    !===Create the connectivity arrays Xind and Xadj based on neigh (for Metis)
    nb_neigh = SIZE(mesh%neigh,1)
    ALLOCATE(xind_conc(me_conc+1), xind_u(me_u+1), xind_T(me_T+1), xind_h(me_h+1), xind_phi(me_phi+1))
    xind_conc(1) = 1
    DO k = 1, me_conc
       m = conc2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_conc))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_conc(k+1) = xind_conc(k) + nb
    END DO
    xind_u(1) = 1
    DO k = 1, me_u
       m = u2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_u))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_u(k+1) = xind_u(k) + nb
    END DO
    xind_T(1) = 1
    DO k = 1, me_T
       m = T2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_T))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_T(k+1) = xind_T(k) + nb
    END DO
    xind_h(1) = 1
    DO k = 1, me_h
       m = h2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_h))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_h(k+1) = xind_h(k) + nb
    END DO
    xind_phi(1) = 1
    DO k = 1, me_phi
       m = phi2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_phi))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_phi(k+1) = xind_phi(k) + nb
    END DO

    ALLOCATE(xadj_conc(xind_conc(me_conc+1)-1))
    ALLOCATE(xadj_u(xind_u(me_u+1)-1))
    ALLOCATE(xadj_T(xind_T(me_T+1)-1))
    ALLOCATE(xadj_h(xind_h(me_h+1)-1))
    ALLOCATE(xadj_phi(xind_phi(me_phi+1)-1))
    p = 0
    DO k = 1, me_conc
       m = conc2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_conc))/=0) CYCLE
          p = p + 1
          xadj_conc(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_conc(me_conc+1)-1) THEN
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, p/=xind_conc(me_conc+1)-1')
    END IF
    p = 0
    DO k = 1, me_u
       m = u2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_u))/=0) CYCLE
          p = p + 1
          xadj_u(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_u(me_u+1)-1) THEN
       CALL error_Petsc('BUG in  part_mesh_M_T_H_phi, p/=xind_u(me_u+1)-1')
    END IF
    p = 0
    DO k = 1, me_T
       m = T2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_T))/=0) CYCLE
          p = p + 1
          xadj_T(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_T(me_T+1)-1) THEN
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, p/=xind_T(me_T+1)-1')
    END IF
    p = 0
    DO k = 1, me_h
       m = h2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_h))/=0) CYCLE
          p = p + 1
          xadj_h(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_h(me_h+1)-1) THEN
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, p/=xind_h(me_h+1)-1')
    END IF
    p = 0
    DO k = 1, me_phi
       m = phi2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_phi))/=0) CYCLE
          p = p + 1
          xadj_phi(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_phi(me_phi+1)-1) THEN
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, p/=xind_phi(me_phi+1)-1')
    END IF
    !===End Create the connectivity arrays Xind and Xadj based on neigh (for Metis)

    !===Create partitions
    opts    = 0
    numflag = 1
    wgtflag = 2
    ALLOCATE(tpwgts(nb_proc))
    tpwgts=one_K4/nb_proc
    CALL METIS_SetDefaultOptions(metis_opt)
    metis_opt(METIS_OPTION_NUMBERING)=1
    ubvec=1.01
    IF (me_conc /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_conc), adjwgt(SIZE(xadj_conc)), part_conc(me_conc))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_conc, 1, xind_conc, xadj_conc, vwgt, vwgt, adjwgt, nb_proc,tpwgts, &
            ubvec, metis_opt, edge, part_conc)
    END IF
    IF (me_u /= 0) THEN
       ALLOCATE(vwgt(me_u), adjwgt(SIZE(xadj_u)), part_u(me_u))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_u, 1, xind_u, xadj_u, vwgt, vwgt, adjwgt, nb_proc, tpwgts, &
            ubvec, metis_opt, edge, part_u)
    END IF
    IF (me_T /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_T), adjwgt(SIZE(xadj_T)), part_T(me_T))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_T, 1, xind_T, xadj_T, vwgt, vwgt, adjwgt, nb_proc, tpwgts, &
            ubvec, metis_opt, edge, part_T)
    END IF
    IF (me_h /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_h), adjwgt(SIZE(xadj_h)), part_h(me_h))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_h, 1, xind_h, xadj_h, vwgt, vwgt, adjwgt, nb_proc,tpwgts, &
            ubvec, metis_opt, edge, part_h)
    END IF
    IF (me_phi /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_phi), adjwgt(SIZE(xadj_phi)), part_phi(me_phi))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_phi, 1, xind_phi,xadj_phi,vwgt, vwgt, adjwgt, nb_proc,tpwgts, &
            ubvec, metis_opt, edge, part_phi)
    END IF
    !===End Create partitions

    !===Create global partition 'part'
    part = -1
    IF (me_conc/=0) THEN
       part(conc2glob(:)) = part_conc
    END IF
    IF (me_u/=0) THEN
       part(u2glob(:)) = part_u
    END IF
    IF (me_T/=0) THEN
       part(T2glob(:)) = part_T
    END IF
    IF (me_h/=0) THEN
       part(h2glob(:)) = part_h
    END IF
    IF (me_phi/=0) THEN
       part(phi2glob(:)) = part_phi
    END IF
    IF (MINVAL(part)==-1) THEN
       CALL error_Petsc('BUG in part_mesh_mhd_bis, MINVAL(part) == -1')
    END IF
    !===End Create global partition 'part'

    !===Create parts and modify part
    !===Search on the boundary whether ms is on a cut.
    IF (SIZE(mesh%jj,1)/=3) THEN
       write(*,*) 'SIZE(mesh%jj,1)', SIZE(mesh%jj,1)
       CALL error_Petsc('BUG in part_mesh_M_T_H_phi, SIZE(mesh%jj,1)/=3')
    END IF
    indicator = -1
    nws = SIZE( mesh%jjs,1)
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_M_T_H_phi, .NOT.test ')
          END IF
          IF (part(mesh%neighs(ms)) == part(mesh%neighs(msop))) CYCLE !==ms is an internal cut
          proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
          part(mesh%neighs(ms)) = proc !make sure interface are internal
          part(mesh%neighs(msop)) = proc !make sure interface are internal
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
          indicator(mesh%jjs(:,ms)) = proc
          indicator(mesh%jjs(:,msop)) = proc
       END DO
    END IF
    !===Fix the partition so that all the cells having one vertex on an
    !===interface belong to the same processor as those sharing this vertices and
    !===having two vertices on the interface (JLG + DCQ July 22 2015)
    DO m = 1, mesh%me
       j_loc = mesh%jj(:,m)
       n = MAXVAL(indicator(j_loc))
       IF (n == -1) CYCLE
       IF (indicator(j_loc(1))*indicator(j_loc(2))*indicator(j_loc(3))<0) CYCLE
       part(m) = n
    END DO
    !===End create parts and modify part

    !===Move the two elements with one periodic face on same processor
    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN
          DO j = 1, mesh%np
             per_pts(j,1) = j
          END DO
          per_pts(:,2:3) = 0
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF ((MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /=0) .AND. &
                  (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(2,:))) /=0) ) CYCLE
             DO ns = 1, SIZE(mesh%jjs,1)
                j = mesh%jjs(ns,ms)
                per_pts(j,2) = m
                DO msop = 1, mesh%mes
                   IF (MINVAL(ABS(mesh%sides(msop)-my_periodic%list_periodic(:,:))) /=0 ) CYCLE
                   IF (msop == ms) CYCLE
                   DO nsop = 1, SIZE(mesh%jjs,1)
                      IF (mesh%jjs(nsop,msop)==j) THEN
                         per_pts(j,3) = mesh%neighs(msop)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          CALL reassign_per_pts(mesh, part, per_pts)
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /= 0) CYCLE
             jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
             s2 = my_periodic%list_periodic(2,jm_loc(1))
             test = .FALSE.
             DO msop = 1, mesh%mes
                IF (mesh%sides(msop) /= s2) CYCLE
                err = 0.d0
                DO ns = 1, SIZE(my_periodic%vect_e,1)
                   err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                        +my_periodic%vect_e(ns,jm_loc(1))))
                END DO
                IF (err .LE. epsilon) THEN
                   test = .TRUE.
                   EXIT
                END IF
             END DO
             IF (.NOT.test) THEN
                CALL error_Petsc('BUG in part_mesh_M_T_H_phi, mop not found')
             END IF
             IF (part(mesh%neighs(ms)) /= part(mesh%neighs(msop))) THEN !==ms is an internal cut
                proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
                part(mesh%neighs(ms)) = proc !make sure interface are internal
                part(mesh%neighs(msop)) = proc !make sure interface are internal
             END IF
          END DO
       END IF
    END IF
    !===End Move the two elements with one periodic face on same processor

!!$ WARNING, FL 1/2/13 : TO BE ADDED IF NEEDED
    !================================================
    IF (rank==0) THEN
       CALL plot_const_p1_label(mesh%jj, mesh%rr, 1.d0*part, 'dd.plt')
    END IF
    !================================================
!!$ WARNING, FL 1/2/13 : TO BE ADDED IF NEEDED

    DEALLOCATE(vwgt,adjwgt)
    IF (ALLOCATED(xadj_conc)) DEALLOCATE(xadj_conc)
    IF (ALLOCATED(xadj_u)) DEALLOCATE(xadj_u)
    IF (ALLOCATED(xadj_T)) DEALLOCATE(xadj_T)
    IF (ALLOCATED(xadj_h)) DEALLOCATE(xadj_h)
    IF (ALLOCATED(xadj_phi)) DEALLOCATE(xadj_phi)
    IF (ALLOCATED(list_T)) DEALLOCATE(list_T)
    IF (ALLOCATED(list_h)) DEALLOCATE(list_h)
    IF (ALLOCATED(list_u)) DEALLOCATE(list_u)
    IF (ALLOCATED(xind_conc)) DEALLOCATE(xind_conc)
    IF (ALLOCATED(xind_u)) DEALLOCATE(xind_u)
    IF (ALLOCATED(xind_T)) DEALLOCATE(xind_T)
    IF (ALLOCATED(xind_h)) DEALLOCATE(xind_h)
    IF (ALLOCATED(xind_phi)) DEALLOCATE(xind_phi)
    IF (ALLOCATED(conc2glob)) DEALLOCATE(conc2glob)
    IF (ALLOCATED(u2glob)) DEALLOCATE(u2glob)
    IF (ALLOCATED(T2glob)) DEALLOCATE(T2glob)
    IF (ALLOCATED(h2glob)) DEALLOCATE(h2glob)
    IF (ALLOCATED(phi2glob)) DEALLOCATE(phi2glob)
    IF (ALLOCATED(part_conc)) DEALLOCATE(part_conc)
    IF (ALLOCATED(part_u)) DEALLOCATE(part_u)
    IF (ALLOCATED(part_T)) DEALLOCATE(part_T)
    IF (ALLOCATED(part_h)) DEALLOCATE(part_h)
    IF (ALLOCATED(part_phi)) DEALLOCATE(part_phi)

    DEALLOCATE(tpwgts)

  END SUBROUTINE part_mesh_M_T_H_phi

  SUBROUTINE part_mesh_mhd(nb_proc,vwgt,mesh,list_of_interfaces,part,my_periodic)
    USE def_type_mesh
    USE my_util
    USE sub_plot
    USE periodic
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    INTEGER, DIMENSION(mesh%me+1)            :: xind
    INTEGER, DIMENSION(mesh%me)              :: vwgt, part
    INTEGER, DIMENSION(:)                    :: list_of_interfaces
    TYPE(periodic_data), OPTIONAL            :: my_periodic

    LOGICAL, DIMENSION(mesh%mes)             :: virgins
    INTEGER, DIMENSION(3,mesh%me)            :: neigh_new
    INTEGER, DIMENSION(5)                    :: opts
    INTEGER, DIMENSION(SIZE(mesh%jjs,1))     :: i_loc
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: xadj, adjwgt
    INTEGER, DIMENSION(1)                    :: jm_loc
    INTEGER, DIMENSION(mesh%np,3)            :: per_pts
    INTEGER :: nb_neigh, edge, m, ms, n, nb, numflag, p, wgtflag, j, &
         ns, nws, msop, nsop, proc, iop, mop, s2
    REAL(KIND=8)                             :: err
    !===(JLG) Feb 20, 2019. Petsc developpers decide to use REAL(KIND=4) to interface with metis
    !REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tpwgts
    !REAL(KIND=8), DIMENSION(1)               :: ubvec
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE  :: tpwgts
    REAL(KIND=4), DIMENSION(1)               :: ubvec
    REAL(KIND=4)                             :: one_K4=1.0
    !===(JLG)Feb 20, 2019.
    INTEGER, DIMENSION(METIS_NOPTIONS)       :: metis_opt
    LOGICAL               :: test
    PetscMPIInt    :: nb_proc
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    IF (nb_proc==1) THEN
       part = 1
       RETURN
    END IF
    nws = SIZE( mesh%jjs,1)
    neigh_new = mesh%neigh
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_mhd, .NOT.test ')
          END IF
          DO n = 1, 3
             IF (neigh_new(n,mesh%neighs(msop))==0) THEN
                neigh_new(n,mesh%neighs(msop)) = mesh%neighs(ms)
             END IF
             IF (neigh_new(n,mesh%neighs(ms))==0) THEN
                neigh_new(n,mesh%neighs(ms)) = mesh%neighs(msop)
             END IF
          END DO
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
       END DO
    END IF

    !=================TEST PERIODIC==================
    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) == 0) THEN
                jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
                s2 = my_periodic%list_periodic(2,jm_loc(1))
                test = .FALSE.
                DO msop = 1, mesh%mes
                   IF (mesh%sides(msop) /= s2) CYCLE

                   err = 0.d0
                   DO ns = 1, SIZE(my_periodic%vect_e,1)
                      err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                           +my_periodic%vect_e(ns,jm_loc(1))))
                   END DO

                   IF (err .LE. epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   CALL error_Petsc('BUG, mop non trouve')
                END IF
                mop = mesh%neighs(msop)
                DO n = 1, 3
                   IF (neigh_new(n,m) == 0) THEN
                      neigh_new(n,m) = mop
                   END IF
                   IF (neigh_new(n,mop) == 0) THEN
                      neigh_new(n,mop) = m
                   END IF
                END DO
             END IF
          END DO
       END IF
    END IF
    !================================================


    ! Create the connectivity array based on neigh
    nb_neigh = SIZE(mesh%neigh,1)
    xind(1) = 1
    DO m = 1, mesh%me
       nb = 0
       DO n = 1, nb_neigh
          IF (neigh_new(n,m)==0) CYCLE
          nb = nb + 1
       END DO
       xind(m+1) = xind(m) + nb
    END DO
    ALLOCATE(xadj(xind(mesh%me+1)-1))
    p = 0
    DO m = 1, mesh%me
       DO n = 1, nb_neigh
          IF (neigh_new(n,m)==0) CYCLE
          p = p + 1
          xadj(p) = neigh_new(n,m)
       END DO
    END DO
    IF (p/=xind(mesh%me+1)-1) THEN
       CALL error_Petsc('BUG, p/=xind(mesh%me+1)-1')
    END IF
    ! End create the connectivity array based on neigh

    ALLOCATE(adjwgt(SIZE(xadj)))
    opts    = 0
    !vwgt    = 1
    adjwgt  = 1
    numflag = 1 ! Fortran numbering of processors
    wgtflag = 2
    ALLOCATE(tpwgts(nb_proc))
    tpwgts=one_K4/nb_proc
    CALL METIS_SetDefaultOptions(metis_opt)
    metis_opt(METIS_OPTION_NUMBERING)=1
    ubvec=1.01
    CALL METIS_PartGraphRecursive(mesh%me, 1, xind,xadj,vwgt, vwgt, adjwgt, nb_proc,tpwgts , ubvec, metis_opt, edge, part)
!!$    CALL METIS_PartGraphRecursive(mesh%me,xind,xadj,vwgt,adjwgt,wgtflag, numflag, nb_proc, opts, edge, part)

    ! Create parts and modify part
    !==Search on the boundary whether ms is on a cut.
    nws = SIZE( mesh%jjs,1)
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in create_local_mesh, .NOT.test ')
          END IF
          IF (part(mesh%neighs(ms)) == part(mesh%neighs(msop))) CYCLE !==ms is an internal cut
          proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
          part(mesh%neighs(ms)) = proc !make sure interface are internal
          part(mesh%neighs(msop)) = proc !make sure interface are internal
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
       END DO
    END IF
    ! End create parts and modify part

    !=================TEST PERIODIC==================
    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN


          DO j = 1, mesh%np
             per_pts(j,1) = j
          END DO
          per_pts(:,2:3) = 0
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF ((MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /=0) .AND. &
                  (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(2,:))) /=0) ) CYCLE
             DO ns = 1, SIZE(mesh%jjs,1)
                j = mesh%jjs(ns,ms)
                per_pts(j,2) = m
                DO msop = 1, mesh%mes
                   IF (MINVAL(ABS(mesh%sides(msop)-my_periodic%list_periodic(:,:))) /=0 ) CYCLE
                   IF (msop == ms) CYCLE
                   DO nsop = 1, SIZE(mesh%jjs,1)
                      IF (mesh%jjs(nsop,msop)==j) THEN
                         per_pts(j,3) = mesh%neighs(msop)
                      END IF
                   END DO
                END DO
             END DO
          END DO
!!$          WRITE(*,*) per_pts(:,2)
          CALL reassign_per_pts(mesh, part, per_pts)

          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /= 0) CYCLE
             jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
             s2 = my_periodic%list_periodic(2,jm_loc(1))
             test = .FALSE.
             DO msop = 1, mesh%mes
                IF (mesh%sides(msop) /= s2) CYCLE

                err = 0.d0
                DO ns = 1, SIZE(my_periodic%vect_e,1)
                   err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                        +my_periodic%vect_e(ns,jm_loc(1))))
                END DO

                IF (err .LE. epsilon) THEN
                   test = .TRUE.
                   EXIT
                END IF
             END DO
             IF (.NOT.test) THEN
                CALL error_Petsc('BUG, mop non trouve')
             END IF
             IF (part(mesh%neighs(ms)) /= part(mesh%neighs(msop))) THEN !==ms is an internal cut
                proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
                part(mesh%neighs(ms)) = proc !make sure interface are internal
                part(mesh%neighs(msop)) = proc !make sure interface are internal
             END IF
          END DO

       END IF
    END IF

!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    !================================================
    CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
    IF (rank==0) THEN
       CALL plot_const_p1_label(mesh%jj, mesh%rr, 1.d0*part, 'dd.plt')
    END IF
    !================================================
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    DEALLOCATE(xadj,adjwgt)
    DEALLOCATE(tpwgts)

  END SUBROUTINE part_mesh_mhd


  SUBROUTINE part_mesh_mhd_bis(nb_proc,list_u, list_h_in, list_phi ,mesh,list_of_interfaces,part,my_periodic)
    USE def_type_mesh
    USE my_util
    USE sub_plot
    USE periodic
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    INTEGER, DIMENSION(mesh%me)              :: part
    INTEGER, DIMENSION(:)                    :: list_of_interfaces
    INTEGER, DIMENSION(:)                    :: list_u, list_h_in, list_phi
    TYPE(periodic_data), OPTIONAL            :: my_periodic

    LOGICAL, DIMENSION(mesh%mes)             :: virgins
    INTEGER, DIMENSION(3,mesh%me)            :: neigh_new
    INTEGER, DIMENSION(5)                    :: opts
    INTEGER, DIMENSION(SIZE(mesh%jjs,1))     :: i_loc
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: xadj_u, xadj_h, xadj_phi, list_h
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: xind_u, xind_h, xind_phi
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: vwgt, adjwgt
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: u2glob, h2glob, phi2glob
    INTEGER, DIMENSION(:),   ALLOCATABLE     :: part_u, part_h, part_phi
    INTEGER, DIMENSION(1)                    :: jm_loc
!!$    INTEGER, DIMENSION(1)                    :: jloc
    INTEGER, DIMENSION(mesh%np,3)            :: per_pts
    INTEGER, DIMENSION(mesh%me)              :: glob2loc
    INTEGER, DIMENSION(mesh%np)              :: indicator
    INTEGER, DIMENSION(3)                    :: j_loc
    INTEGER :: nb_neigh, edge, m, ms, n, nb, numflag, p, wgtflag, j, &
         ns, nws, msop, nsop, proc, iop, mop, s2, k, me_u, me_h, me_phi, idm
    REAL(KIND=8)                             :: err
    LOGICAL               :: test
    !===(JLG) Feb 20, 2019. Petsc developpers decide to use REAL(KIND=4) to interface with metis
    !REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tpwgts
    !REAL(KIND=8), DIMENSION(1)               :: ubvec
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE  :: tpwgts
    REAL(KIND=4), DIMENSION(1)               :: ubvec
    REAL(KIND=4)                             :: one_K4=1.0
    !===(JLG)Feb 20, 2019.
    INTEGER, DIMENSION(METIS_NOPTIONS)       :: metis_opt
    PetscMPIInt    :: nb_proc
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    IF (nb_proc==1) THEN
       part = 1
       RETURN
    END IF
    glob2loc = 0
    !Create list_h = list_h_in \ list_ns
    nb = 0
    DO j = 1, SIZE(list_h_in)
       IF (MINVAL(ABS(list_h_in(j)-list_u))==0) CYCLE
       nb = nb + 1
    END DO
    ALLOCATE(list_h(nb))
    nb = 0
    DO j = 1, SIZE(list_h_in)
       IF (MINVAL(ABS(list_h_in(j)-list_u))==0) CYCLE
       nb = nb + 1
       list_h(nb) = list_h_in(j)
    END DO

    !Done with list_h

    nws = SIZE( mesh%jjs,1)
    neigh_new = mesh%neigh
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in part_mesh_mhd, .NOT.test ')
          END IF
          DO n = 1, 3
             IF (neigh_new(n,mesh%neighs(msop))==0) THEN
                neigh_new(n,mesh%neighs(msop)) = mesh%neighs(ms)
             END IF
             IF (neigh_new(n,mesh%neighs(ms))==0) THEN
                neigh_new(n,mesh%neighs(ms)) = mesh%neighs(msop)
             END IF
          END DO
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
       END DO
    END IF

    !=================TEST PERIODIC==================

    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) == 0) THEN
                jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
                s2 = my_periodic%list_periodic(2,jm_loc(1))
                test = .FALSE.
                DO msop = 1, mesh%mes
                   IF (mesh%sides(msop) /= s2) CYCLE

                   err = 0.d0
                   DO ns = 1, SIZE(my_periodic%vect_e,1)
                      err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                           +my_periodic%vect_e(ns,jm_loc(1))))
                   END DO

                   IF (err .LE. epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   CALL error_Petsc('BUG, mop non trouve')
                END IF
                mop = mesh%neighs(msop)
                DO n = 1, 3
                   IF (neigh_new(n,m) == 0) THEN
                      neigh_new(n,m) = mop
                   END IF
                   IF (neigh_new(n,mop) == 0) THEN
                      neigh_new(n,mop) = m
                   END IF
                END DO
             END IF
          END DO
       END IF
    END IF

    !================================================
    me_u = 0
    me_h = 0
    me_phi = 0
    DO m = 1, mesh%me
       idm = mesh%i_d(m)
       IF (MINVAL(ABS(idm-list_u))==0) THEN
          me_u = me_u + 1
       ELSE IF (MINVAL(ABS(idm-list_h))==0) THEN
          me_h = me_h + 1
       ELSE IF (MINVAL(ABS(idm-list_phi))==0) THEN
          me_phi = me_phi + 1
       ELSE
          CALL error_Petsc('BUG in part_mesh_mhd_bis : element not in the mesh')
       END IF
    END DO
    ALLOCATE(u2glob(me_u), h2glob(me_h),phi2glob(me_phi))
    me_u = 0
    me_h = 0
    me_phi = 0
    DO m = 1, mesh%me
       idm = mesh%i_d(m)
       IF (MINVAL(ABS(idm-list_u))==0) THEN
          me_u = me_u + 1
          u2glob(me_u) = m
          glob2loc(m) = me_u
       ELSE IF (MINVAL(ABS(idm-list_h))==0) THEN
          me_h = me_h + 1
          h2glob(me_h) = m
          glob2loc(m) = me_h
       ELSE IF (MINVAL(ABS(idm-list_phi))==0) THEN
          me_phi = me_phi + 1
          phi2glob(me_phi) = m
          glob2loc(m) = me_phi
       ELSE
          CALL error_Petsc('BUG in part_mesh_mhd_bis : element not in the mesh')
       END IF
    END DO
    ! Create the connectivity arrays based on neigh
    nb_neigh = SIZE(mesh%neigh,1)
    ALLOCATE(xind_u(me_u+1), xind_h(me_h+1), xind_phi(me_phi+1))
    xind_u(1) = 1
    DO k = 1, me_u
       m = u2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_u))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_u(k+1) = xind_u(k) + nb
    END DO
    xind_h(1) = 1
    DO k = 1, me_h
       m = h2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_h))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_h(k+1) = xind_h(k) + nb
    END DO
    xind_phi(1) = 1
    DO k = 1, me_phi
       m = phi2glob(k)
       nb = 0
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_phi))/=0) CYCLE
          nb = nb + 1
       END DO
       xind_phi(k+1) = xind_phi(k) + nb
    END DO

    ALLOCATE(xadj_u(xind_u(me_u+1)-1))
    ALLOCATE(xadj_h(xind_h(me_h+1)-1))
    ALLOCATE(xadj_phi(xind_phi(me_phi+1)-1))
    p = 0
    DO k = 1, me_u
       m = u2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_u))/=0) CYCLE
!!$          jloc = MINLOC(ABS(u2glob-mop))
          p = p + 1
!!$          xadj_u(p) = jloc(1)
          xadj_u(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_u(me_u+1)-1) THEN
       CALL error_Petsc('BUG, p/=xind_u(me_u+1)-1')
    END IF
    p = 0
    DO k = 1, me_h
       m = h2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_h))/=0) CYCLE
!!$          jloc = MINLOC(ABS(h2glob-mop))
          p = p + 1
!!$          xadj_h(p) = jloc(1)
          xadj_h(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_h(me_h+1)-1) THEN
       CALL error_Petsc('BUG, p/=xind_h(me_h+1)-1')
    END IF
    p = 0
    DO k = 1, me_phi
       m = phi2glob(k)
       DO n = 1, nb_neigh
          mop = neigh_new(n,m)
          IF (mop==0) CYCLE
          IF (MINVAL(ABS(mesh%i_d(mop)-list_phi))/=0) CYCLE
!!$          jloc = MINLOC(ABS(phi2glob-mop))
          p = p + 1
!!$          xadj_phi(p) = jloc(1)
          xadj_phi(p) = glob2loc(mop)
       END DO
    END DO
    IF (p/=xind_phi(me_phi+1)-1) THEN
       CALL error_Petsc('BUG, p/=xind_phi(me_phi+1)-1')
    END IF

    opts    = 0
    numflag = 1
    wgtflag = 2
    ! Create partitions
    ALLOCATE(tpwgts(nb_proc))
    tpwgts=one_K4/nb_proc
    CALL METIS_SetDefaultOptions(metis_opt)
    metis_opt(METIS_OPTION_NUMBERING)=1
    ubvec=1.01
    IF (me_u /= 0) THEN
       ALLOCATE(vwgt(me_u), adjwgt(SIZE(xadj_u)), part_u(me_u))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_u, 1, xind_u,xadj_u,vwgt, vwgt, adjwgt, nb_proc,tpwgts , ubvec, metis_opt, edge, part_u)
!!$       CALL METIS_PartGraphRecursive(me_u,xind_u,xadj_u,vwgt,adjwgt,wgtflag, numflag, nb_proc, opts, edge, part_u)
    END IF
    IF (me_h /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_h), adjwgt(SIZE(xadj_h)), part_h(me_h))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_h, 1, xind_h,xadj_h,vwgt, vwgt, adjwgt, nb_proc,tpwgts , ubvec, metis_opt, edge, part_h)
!!$       CALL METIS_PartGraphRecursive(me_h,xind_h,xadj_h,vwgt,adjwgt,wgtflag, numflag, nb_proc, opts, edge, part_h)
    END IF
    IF (me_phi /= 0) THEN
       IF (ALLOCATED(vwgt)) THEN
          DEALLOCATE(vwgt, adjwgt)
       END IF
       ALLOCATE(vwgt(me_phi), adjwgt(SIZE(xadj_phi)), part_phi(me_phi))
       vwgt   = 1
       adjwgt = 1
       CALL METIS_PartGraphRecursive(me_phi, 1, xind_phi,xadj_phi,vwgt, vwgt, adjwgt, nb_proc,tpwgts, &
            ubvec, metis_opt, edge, part_phi)
!!$       CALL METIS_PartGraphRecursive(me_phi,xind_phi,xadj_phi,vwgt,adjwgt,wgtflag, numflag, nb_proc, opts, edge, part_phi)
    END IF
!!$!TESTTT
!!$    do m = 1,  mesh%me
!!$       if (part_phi(m) ==1) THEN
!!$          part_phi(m) =2
!!$       ELSE if (part_phi(m) ==2) then
!!$          part_phi(m) =1
!!$      ELSE if (part_phi(m) ==3) then
!!$          part_phi(m) =4
!!$     ELSE if (part_phi(m) ==4) then
!!$          part_phi(m) =3
!!$     END IF
!!$   end do
!!$!TESTTT
    ! Create global partition 'part'
    part = -1
    IF (me_u/=0) THEN
       part(u2glob(:)) = part_u
    END IF
    IF (me_h/=0) THEN
       part(h2glob(:)) = part_h
    END IF
    IF (me_phi/=0) THEN
       part(phi2glob(:)) = part_phi
    END IF
    IF (MINVAL(part)==-1) THEN
       CALL error_Petsc('BUG in part_mesh_mhd_bis, MINVAL(part) == -1')
    END IF
    ! End Create global partition

    ! Create parts and modify part
    !==Search on the boundary whether ms is on a cut.
    IF (SIZE(mesh%jj,1)/=3) THEN
       write(*,*) 'SIZE(mesh%jj,1)', SIZE(mesh%jj,1)
       CALL error_Petsc('BUG in part_mesh_mhd_bis, SIZE(mesh%jj,1)/=3')
    END IF
    indicator = -1
    nws = SIZE( mesh%jjs,1)
    IF (SIZE(list_of_interfaces)/=0) THEN
       virgins = .TRUE.
       DO ms = 1, mesh%mes
          IF (.NOT.virgins(ms)) CYCLE
          IF (MINVAL(ABS(mesh%sides(ms)-list_of_interfaces))/=0) CYCLE !==ms not on a cut
          i_loc = mesh%jjs(:,ms)
          DO msop = 1, mesh%mes
             IF (msop==ms .OR. .NOT.virgins(msop)) CYCLE
             IF (MINVAL(ABS(mesh%sides(msop)-list_of_interfaces))/=0) CYCLE !==msop not on a cut
             DO ns = 1, nws
                test = .FALSE.
                DO nsop = 1, nws
                   iop = mesh%jjs(nsop,msop)
                   IF (MAXVAL(ABS(mesh%rr(:,i_loc(ns))-mesh%rr(:,iop))).LT.epsilon) THEN
                      test = .TRUE.
                      EXIT
                   END IF
                END DO
                IF (.NOT.test) THEN
                   EXIT !==This msop does not coincide with ms
                END IF
             END DO
             IF (test) EXIT
          END DO
          IF (.NOT.test) THEN
             CALL error_Petsc('BUG in create_local_mesh, .NOT.test ')
          END IF
          IF (part(mesh%neighs(ms)) == part(mesh%neighs(msop))) CYCLE !==ms is an internal cut
          proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
          part(mesh%neighs(ms)) = proc !make sure interface are internal
          part(mesh%neighs(msop)) = proc !make sure interface are internal
          virgins(ms) = .FALSE.
          virgins(msop) = .FALSE.
          indicator(mesh%jjs(:,ms)) = proc
          indicator(mesh%jjs(:,msop)) = proc
       END DO
    END IF
    !===Fix the partition so that all the cells having one vertex on an
    !===interface belong to the same processor as those sharing this vertices and
    !===having two vertices on the interface (JLG + DCL July 22 2015)
    DO m = 1, mesh%me
       j_loc = mesh%jj(:,m)
       n = MAXVAL(indicator(j_loc))
       IF (n == -1) CYCLE
       IF (indicator(j_loc(1))*indicator(j_loc(2))*indicator(j_loc(3))<0) CYCLE
       part(m) = n
    END DO
    ! End create parts and modify part

    !=================TEST PERIODIC==================
    IF (PRESENT(my_periodic)) THEN
       IF (my_periodic%nb_periodic_pairs/=0) THEN


          DO j = 1, mesh%np
             per_pts(j,1) = j
          END DO
          per_pts(:,2:3) = 0
          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF ((MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /=0) .AND. &
                  (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(2,:))) /=0) ) CYCLE
             DO ns = 1, SIZE(mesh%jjs,1)
                j = mesh%jjs(ns,ms)
                per_pts(j,2) = m
                DO msop = 1, mesh%mes
                   IF (MINVAL(ABS(mesh%sides(msop)-my_periodic%list_periodic(:,:))) /=0 ) CYCLE
                   IF (msop == ms) CYCLE
                   DO nsop = 1, SIZE(mesh%jjs,1)
                      IF (mesh%jjs(nsop,msop)==j) THEN
                         per_pts(j,3) = mesh%neighs(msop)
                      END IF
                   END DO
                END DO
             END DO
          END DO
!!$          WRITE(*,*) per_pts(:,2)
          CALL reassign_per_pts(mesh, part, per_pts)

          DO ms = 1, mesh%mes
             m = mesh%neighs(ms)
             IF (MINVAL(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:))) /= 0) CYCLE
             jm_loc = MINLOC(ABS(mesh%sides(ms)-my_periodic%list_periodic(1,:)))
             s2 = my_periodic%list_periodic(2,jm_loc(1))
             test = .FALSE.
             DO msop = 1, mesh%mes
                IF (mesh%sides(msop) /= s2) CYCLE

                err = 0.d0
                DO ns = 1, SIZE(my_periodic%vect_e,1)
                   err = err + ABS(SUM(mesh%rr(ns,mesh%jjs(:,ms))-mesh%rr(ns,mesh%jjs(:,msop)) &
                        +my_periodic%vect_e(ns,jm_loc(1))))
                END DO

                IF (err .LE. epsilon) THEN
                   test = .TRUE.
                   EXIT
                END IF
             END DO
             IF (.NOT.test) THEN
                CALL error_Petsc('BUG, mop non trouve')
             END IF
             IF (part(mesh%neighs(ms)) /= part(mesh%neighs(msop))) THEN !==ms is an internal cut
                proc = MIN(part(mesh%neighs(ms)),part(mesh%neighs(msop)))
                part(mesh%neighs(ms)) = proc !make sure interface are internal
                part(mesh%neighs(msop)) = proc !make sure interface are internal
             END IF
          END DO

       END IF
    END IF

!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED
    !================================================
    CALL MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
    IF (rank==0) THEN
       CALL plot_const_p1_label(mesh%jj, mesh%rr, 1.d0*part, 'dd.plt')
    END IF
    !================================================
!!$ WARNING, FL 1/2/13 : TO BE ADDED IN NEEDED

    DEALLOCATE(vwgt,adjwgt)
    IF (ALLOCATED(xadj_u)) DEALLOCATE(xadj_u)
    IF (ALLOCATED(xadj_h)) DEALLOCATE(xadj_h)
    IF (ALLOCATED(xadj_phi)) DEALLOCATE(xadj_phi)
    IF (ALLOCATED(list_h)) DEALLOCATE(list_h)
    IF (ALLOCATED(xind_u)) DEALLOCATE(xind_u)
    IF (ALLOCATED(xind_h)) DEALLOCATE(xind_h)
    IF (ALLOCATED(xind_phi)) DEALLOCATE(xind_phi)
    IF (ALLOCATED(u2glob)) DEALLOCATE(u2glob)
    IF (ALLOCATED(h2glob)) DEALLOCATE(h2glob)
    IF (ALLOCATED(phi2glob)) DEALLOCATE(phi2glob)
    IF (ALLOCATED(part_u)) DEALLOCATE(part_u)
    IF (ALLOCATED(part_h)) DEALLOCATE(part_h)
    IF (ALLOCATED(part_phi)) DEALLOCATE(part_phi)

    DEALLOCATE(tpwgts)

  END SUBROUTINE part_mesh_mhd_bis

  SUBROUTINE extract_mesh(communicator,nb_proc,mesh_glob,part,list_dom,mesh,mesh_loc)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh_glob, mesh, mesh_loc
    INTEGER, DIMENSION(:)                    :: part, list_dom
    INTEGER, DIMENSION(mesh_glob%me)         :: bat
    INTEGER, DIMENSION(mesh_glob%np)         :: i_old_to_new
    INTEGER, DIMENSION(mesh_glob%mes)        :: parts
    INTEGER, DIMENSION(nb_proc)              :: nblmt_per_proc, start, displ
    INTEGER, DIMENSION(2)                    :: np_loc, me_loc, mes_loc
    INTEGER, DIMENSION(:), ALLOCATABLE       :: list_m, tab, tabs
    INTEGER                                  :: nb_proc, ms, i, index, m, mop, n
    PetscErrorCode :: ierr
    PetscMPIInt    :: rank
    MPI_Comm       :: communicator
    CALL MPI_Comm_rank(communicator,rank,ierr)

    ! Create parts
    parts = part(mesh_glob%neighs)
    ! End create parts

    ! Create list_m
    i = 0
    DO m = 1, mesh_glob%me
       IF (MINVAL(ABS(list_dom-mesh_glob%i_d(m)))/=0) CYCLE
       i = i + 1
    END DO
    mesh%me = i
    ALLOCATE (list_m(mesh%me))
    i = 0
    DO m = 1, mesh_glob%me
       IF (MINVAL(ABS(list_dom-mesh_glob%i_d(m)))/=0) CYCLE
       i = i + 1
       list_m(i) = m
    END DO
    !End create list_m

    ! Count elements on processors
    nblmt_per_proc = 0
    DO i = 1, mesh%me
       m = list_m(i)
       nblmt_per_proc(part(m)) = nblmt_per_proc(part(m)) + 1
    END DO
    start(1) = 0
    DO n = 2, nb_proc
       start(n) = start(n-1) + nblmt_per_proc(n-1)
    END DO
    me_loc(1) = start(rank+1) + 1
    me_loc(2) = start(rank+1) + nblmt_per_proc(rank+1)
    displ = start
    ! End count elements on processors

    ! Re-order elements
    ALLOCATE(tab(mesh%me))
    bat = 0
    DO i = 1, mesh%me
       m = list_m(i)
       start(part(m)) = start(part(m)) + 1
       tab(start(part(m))) = m
       bat(m) = start(part(m))
    END DO
    ! Re-order elements

    ! Create mesh%jj
    mesh%gauss%n_w = SIZE(mesh_glob%jj,1)
    ALLOCATE(mesh%jj(SIZE(mesh_glob%jj,1),mesh%me))
    i_old_to_new = 0
    index = 0
    DO m = 1, mesh%me
       DO n = 1, SIZE(mesh_glob%jj,1)
          i = mesh_glob%jj(n,tab(m))
          IF (i_old_to_new(i)/=0) THEN
             mesh%jj(n,m) = i_old_to_new(i)
          ELSE
             index = index + 1
             i_old_to_new(i) = index
             mesh%jj(n,m) = i_old_to_new(i)
          END IF
       END DO
    END DO
    ! End Create mesh%jj

    ! Create mesh%rr
    mesh%np = index
    ALLOCATE(mesh%rr(2,mesh%np))
    DO i = 1, mesh_glob%np
       IF (i_old_to_new(i)==0) CYCLE
       mesh%rr(:,i_old_to_new(i)) = mesh_glob%rr(:,i)
    END DO
    !End Create mesh%rr

    ! Create mesh%neigh
    ALLOCATE(mesh%neigh(3,mesh%me))
    DO m = 1, mesh%me
       DO n = 1, 3
          mop = mesh_glob%neigh(n,tab(m))
          IF (mop==0) THEN
             mesh%neigh(n,m) = 0
          ELSE
             mesh%neigh(n,m) = bat(mop)
          END IF
       END DO
    END DO
    ! End  Create mesh%neigh

    ! Create mesh%i_d
    ALLOCATE(mesh%i_d(mesh%me))
    mesh%i_d =  mesh_glob%i_d(tab)
    ! End mesh%i_d

    ! Create np_loc
    IF (displ(rank+1)/=0) THEN
       np_loc(1) = MAXVAL(mesh%jj(:,1:displ(rank+1))) + 1
    ELSE
       np_loc(1) = 1
    END IF
    np_loc(2) = np_loc(1) - 1
    IF (me_loc(1).LE.me_loc(2)) THEN
       np_loc(2) = MAXVAL(mesh%jj(:,me_loc(1):me_loc(2)))
    END IF
    IF (np_loc(2) .LT. np_loc(1)-1) THEN
       np_loc(2) = np_loc(1) - 1
    END IF
    ! End create np_loc

    ! Create mes_loc
    nblmt_per_proc=0
    DO ms = 1, mesh_glob%mes
       IF (MINVAL(ABS(list_dom-mesh_glob%i_d(mesh_glob%neighs(ms))))/=0) CYCLE
       n = parts(ms)
       nblmt_per_proc(n) =  nblmt_per_proc(n) + 1
    END DO
    start(1) = 0
    DO n = 2, nb_proc
       start(n) = start(n-1) + nblmt_per_proc(n-1)
    END DO
    mes_loc(1) = start(rank+1) + 1
    mes_loc(2) = start(rank+1) + nblmt_per_proc(rank+1)
    mesh%mes = SUM(nblmt_per_proc)
    ! End create mes_loc

    ! Create tabs and sbat
    ALLOCATE(tabs(mesh%mes))
    DO ms = 1, mesh_glob%mes
       IF (MINVAL(ABS(list_dom-mesh_glob%i_d(mesh_glob%neighs(ms))))/=0) CYCLE
       start(parts(ms)) = start(parts(ms)) + 1
       tabs(start(parts(ms))) = ms
    END DO
    ! End create tabs and sbat

    ! Create neighs
    ALLOCATE(mesh%neighs(mesh%mes))
    mesh%neighs = bat(mesh_glob%neighs(tabs))
    ! End create neighs

    ! Re-order sides
    ALLOCATE(mesh%sides(mesh%mes))
    mesh%sides = mesh_glob%sides(tabs)
    ! End re-order sides

    ! Re-order jjs
    mesh%gauss%n_ws = SIZE(mesh_glob%jjs,1)
    ALLOCATE(mesh%jjs(SIZE(mesh_glob%jjs,1),mesh%mes))

    DO n = 1, SIZE(mesh%jjs,1)
       mesh%jjs(n,:) = i_old_to_new(mesh_glob%jjs(n,tabs))
    END DO
    ! End re-order jjs

    !==We create the local mesh now
    mesh%edge_stab = .FALSE.
    CALL create_local_mesh(mesh, mesh_loc, me_loc, mes_loc, np_loc)

    DEALLOCATE(list_m, tab, tabs)

  END SUBROUTINE extract_mesh

  SUBROUTINE create_local_mesh(mesh, mesh_loc, me_loc, mes_loc, np_loc)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh, mesh_loc
    INTEGER, DIMENSION(2),   INTENT(IN)      :: me_loc, mes_loc, np_loc
    INTEGER, DIMENSION(mesh%me)              :: m_glob_to_loc
    INTEGER, DIMENSION(:), ALLOCATABLE       :: m_loc_to_glob
    INTEGER, DIMENSION(mesh%np)              :: glob_to_loc,loc_to_glob
    LOGICAL, DIMENSION(mesh%np)              :: virgin
    INTEGER :: dim, nws, nw, m, ms, mop, ns, msup, minf, dof, &
         dom_me, nwc, dom_mes, dom_np, n, i
    LOGICAL :: test

    dim = SIZE(mesh%rr,1)
    nws = SIZE(mesh%jjs,1)
    nw = SIZE(mesh%jj,1)
    nwc = SIZE(mesh%neigh,1)
    !==Test if one proc only
    IF (me_loc(2) - me_loc(1) + 1==mesh%me) THEN
       mesh_loc%me = mesh%me
       mesh_loc%np = mesh%np
       mesh_loc%mes = mesh%mes
       mesh_loc%dom_me = mesh%me
       mesh_loc%dom_np = mesh%np
       mesh_loc%dom_mes = mesh%mes
       mesh_loc%gauss%n_w = nw
       ALLOCATE(mesh_loc%jj(nw,mesh%me))
       mesh_loc%jj = mesh%jj
       ALLOCATE(mesh_loc%loc_to_glob(mesh%np))
       DO n = 1, mesh%np
          mesh_loc%loc_to_glob(n) = n
       END DO
       ALLOCATE(mesh_loc%rr(dim,mesh%np))
       mesh_loc%rr=mesh%rr
       ALLOCATE(mesh_loc%neigh(nwc,mesh%me))
       mesh_loc%neigh = mesh%neigh
       ALLOCATE(mesh_loc%i_d(mesh%me))
       mesh_loc%i_d =  mesh%i_d
       ALLOCATE(mesh_loc%neighs(mesh_loc%mes))
       mesh_loc%neighs = mesh%neighs
       ALLOCATE(mesh_loc%sides(mesh_loc%mes))
       mesh_loc%sides = mesh%sides
       mesh_loc%gauss%n_ws = nws
       ALLOCATE(mesh_loc%jjs(nws,mesh_loc%mes))
       mesh_loc%jjs = mesh%jjs
       RETURN
    END IF
    !==End test if one proc only


    !==Create the new mesh
    dom_me = me_loc(2) - me_loc(1) + 1
    dom_mes = mes_loc(2) - mes_loc(1) + 1
    dom_np = np_loc(2) - np_loc(1) + 1
    mesh_loc%me = dom_me
    mesh_loc%mes = dom_mes
    mesh_loc%dom_me = dom_me
    mesh_loc%dom_np = dom_np
    mesh_loc%dom_mes = dom_mes

    IF (dom_me==0) THEN
       ALLOCATE(mesh_loc%jj(0,0))
       ALLOCATE(mesh_loc%loc_to_glob(0))
       ALLOCATE(mesh_loc%rr(0,0))
       ALLOCATE(mesh_loc%neigh(0,0))
       ALLOCATE(mesh_loc%i_d(0))
       ALLOCATE(mesh_loc%neighs(0))
       ALLOCATE(mesh_loc%sides(0))
       ALLOCATE(mesh_loc%jjs(0,0))
       mesh_loc%gauss%n_w = 0
       mesh_loc%gauss%n_ws = 0
       RETURN
    ELSE IF (dom_me<0) THEN
       CALL error_Petsc('BUG in create_local_mesh, dom_me<0 ')
    ELSE
       mesh_loc%gauss%n_w = nw
       mesh_loc%gauss%n_ws = nws
    END IF

    !==Re-order jj
    virgin = .TRUE.
    dof = 0
    DO m = me_loc(1), me_loc(2)
       DO n = 1, nw
          i = mesh%jj(n,m)
          IF(.NOT.virgin(i) .OR. i.GE.np_loc(1)) CYCLE
          virgin(i) = .FALSE.
          dof = dof + 1
       END DO
    END DO
    ALLOCATE(mesh_loc%jj(nw,mesh_loc%me))

    ALLOCATE(m_loc_to_glob(mesh_loc%me))
    m_glob_to_loc = 0
    virgin = .TRUE.
    dof = dom_np
    DO m = me_loc(1), me_loc(2)
       DO n = 1, nw
          i = mesh%jj(n,m)
          IF(virgin(i)) THEN
             virgin(i) = .FALSE.
             IF (i .LT. np_loc(1)) THEN
                dof = dof + 1
                glob_to_loc(i) = dof
                loc_to_glob(dof) = i
             ELSE
                glob_to_loc(i) = i-np_loc(1) + 1
                loc_to_glob(i-np_loc(1) + 1) = i
             END IF
          END IF
       END DO
       m_loc_to_glob(m-me_loc(1)+1) = m
       m_glob_to_loc(m) = m-me_loc(1)+1
    END DO

    DO n = 1, nw
       mesh_loc%jj(n,1:dom_me) = glob_to_loc(mesh%jj(n,me_loc(1):me_loc(2)))
    END DO
    !==End re-order jj

    !==Create mesh%loc_to_glob
    IF (MAXVAL(mesh_loc%jj)/=dof) THEN
       CALL error_Petsc('BUG in create_local_mesh, mesh_loc%jj)/=dof')
    END IF
    mesh_loc%np = dof
    ALLOCATE(mesh_loc%loc_to_glob(mesh_loc%np))
    mesh_loc%loc_to_glob = loc_to_glob(1:mesh_loc%np)
    !==End create mesh%loc_to_glob

    !==Re-order rr
    ALLOCATE(mesh_loc%rr(dim,mesh_loc%np))
    DO  n = 1, mesh_loc%np
       mesh_loc%rr(:,n) = mesh%rr(:,mesh_loc%loc_to_glob(n))
    END DO
    !==End re-order rr

    !==Re-order neigh
    ALLOCATE(mesh_loc%neigh(nwc,mesh_loc%me))
    msup = MAXVAL(m_loc_to_glob)
    minf = MINVAL(m_loc_to_glob)
    DO m = 1, mesh_loc%me
       DO n = 1, nwc
          mop = mesh%neigh(n,m_loc_to_glob(m))
          IF (mop==0) THEN
             mesh_loc%neigh(n,m) = 0
          ELSE IF(mop<minf .OR. mop>msup) THEN
             !mesh_loc%neigh(n,m) = mop
             mesh_loc%neigh(n,m) = -mop !JLG AUG 13 2015
          ELSE
             mesh_loc%neigh(n,m) = m_glob_to_loc(mop)
          END IF
       END DO
    END DO
    !==End re-order neigh

    !==Re-order i_d
    ALLOCATE(mesh_loc%i_d(mesh_loc%me))
    mesh_loc%i_d =  mesh%i_d(m_loc_to_glob)
    !==End re-order i_d

    !==Re-order neighs
    ALLOCATE(mesh_loc%neighs(mesh_loc%mes))
    mesh_loc%neighs = m_glob_to_loc(mesh%neighs(mes_loc(1):mes_loc(2)))
    !==End re-order neighs


    !==Re-order sides
    ALLOCATE(mesh_loc%sides(mesh_loc%mes))
    mesh_loc%sides = mesh%sides(mes_loc(1):mes_loc(2))
    !==End re-order sides

    !==Re-order jjs
    ALLOCATE(mesh_loc%jjs(nws,mesh_loc%mes))
    DO ns = 1, nws
       mesh_loc%jjs(ns,:) = glob_to_loc(mesh%jjs(ns,mes_loc(1):mes_loc(2)))
    END DO
    !==End re-order jjs

    !TEST

    DO ms = 1, mesh_loc%mes
       m = mesh_loc%neighs(ms)
       DO ns = 1, nws
          test = .TRUE.
          DO n = 1, nw
             IF (MAXVAL(ABS(mesh_loc%rr(:,mesh_loc%jj(n,m))-mesh_loc%rr(:,mesh_loc%jjs(ns,ms)))) .LT. 1.d-10) THEN
                test = .FALSE.
             END IF
          END DO
          IF (test) THEN
             WRITE(*,*) 'bug in create local mesh, non consistent numbering'
          END IF
       END DO
    END DO


    !TEST

    DEALLOCATE(m_loc_to_glob)

  END SUBROUTINE create_local_mesh

  SUBROUTINE free_mesh(mesh)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh

    DEALLOCATE(mesh%jj)
    DEALLOCATE(mesh%jjs)
    DEALLOCATE(mesh%rr)
    DEALLOCATE(mesh%neigh)
    DEALLOCATE(mesh%sides)
    DEALLOCATE(mesh%neighs)
    DEALLOCATE(mesh%i_d)

    NULLIFY(mesh%loc_to_glob)
    NULLIFY(mesh%disp)
    NULLIFY(mesh%domnp)
    NULLIFY(mesh%j_s)

    IF (mesh%edge_stab) THEN
       DEALLOCATE(mesh%iis)
       NULLIFY(mesh%jji)
       DEALLOCATE(mesh%jjsi)
       DEALLOCATE(mesh%neighi)
    END IF

    mesh%dom_me = 0
    mesh%dom_np = 0
    mesh%dom_mes = 0
    mesh%me = 0
    mesh%mes = 0
    mesh%np = 0
    mesh%nps = 0
    mesh%mi =0
    mesh%edge_stab = .FALSE.

  END SUBROUTINE FREE_MESH

  SUBROUTINE free_interface(interf)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(interface_type) :: interf

    interf%mes = 0
    DEALLOCATE(interf%mesh1)
    DEALLOCATE(interf%mesh2)
    DEALLOCATE(interf%jjs1)
    DEALLOCATE(interf%jjs2)
  END SUBROUTINE free_interface

  SUBROUTINE reassign_per_pts(mesh, partition, list_pts)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE

    TYPE(mesh_type)            ,   INTENT(IN)    :: mesh
    INTEGER, DIMENSION(mesh%me),   INTENT(INOUT) :: partition
    INTEGER, DIMENSION(mesh%np,3), INTENT(IN)    :: list_pts

    INTEGER                 :: i, j_loc, proc_min, index, i_loc, m, mop, n, proc1, proc2
    INTEGER, DIMENSION(50)  :: list_elmts
    LOGICAL                 :: okay


    list_elmts = 0
    DO i = 1, mesh%np
       IF (list_pts(i,2)==0) CYCLE
       j_loc = list_pts(i,1)
       list_elmts=0
       index = 1
       list_elmts(index) = list_pts(i,2)
       okay = .TRUE.
       DO WHILE (okay)
          m=list_elmts(index)
          okay = .FALSE.
          i_loc=index
          DO n = 1, 3
             mop = mesh%neigh(n, m)
             IF (mop == 0) CYCLE
             IF (MINVAL(ABS(mesh%jj(:,mop)-j_loc)) /=0) CYCLE
             IF (MINVAL(ABS(mop-list_elmts))==0) CYCLE
             okay=.TRUE.
             i_loc = i_loc + 1
             IF (i_loc-index==2) THEN
                CALL error_Petsc('BUG in reassign_per_pts, how is that possible?')
             END IF
             list_elmts(i_loc) = mop
          END DO
          index = i_loc
       END DO
!!$       WRITE(*,*) i, list_elmts(1:index)
       IF (list_pts(i,3) == 0) THEN  ! point au bord du bord periodique, ou sur une arete
          proc_min = MINVAL(partition(list_elmts(1:index)))
          partition(list_elmts(1)) = proc_min
       ELSE ! deux elements du bord periodique touchent le point
          IF (list_elmts(index) /= list_pts(i,3)) THEN
             CALL error_Petsc('BUG in reassign_per_pts, wrong element')
          END IF
          proc1 = partition(list_elmts(1))
          proc2 = partition(list_elmts(2))
          partition(list_elmts(2:index-1)) = MIN(proc1,proc2)
       END IF
    END DO

  END SUBROUTINE reassign_per_pts


END MODULE metis_sfemans
