!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE st_matrix

  PUBLIC :: st_aij_csr_glob_block, extract, create_my_ghost, st_aij_csr, tri_jlg, st_aij_csr_loc_block, &
       st_csr, st_csr_block, st_csr_bloc
  PRIVATE
#include "petsc/finclude/petsc.h"
CONTAINS

  !=========================================================================

  SUBROUTINE create_my_ghost(mesh,LA,ifrom)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    type(petsc_csr_LA)                       :: LA
    INTEGER, DIMENSION(:), POINTER           :: ifrom
    INTEGER                                  :: kmax, nifrom, start, fin, k
    kmax = SIZE(LA%loc_to_glob,1)
    nifrom = mesh%np-mesh%dom_np
    ALLOCATE(ifrom(kmax*nifrom))
    IF (nifrom/=0) THEN
       DO k = 1, kmax
          start = (k-1)*nifrom+1
          fin = start + nifrom - 1
          ifrom(start:fin)= LA%loc_to_glob(k,mesh%dom_np+1:)-1
       END DO
    END IF
  END SUBROUTINE create_my_ghost

  SUBROUTINE extract(xghost,ks,ke,LA,phi)
#include "petsc/finclude/petscvec.h"
    use petsc
    USE def_type_mesh
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:), INTENT(OUT):: phi
    type(petsc_csr_LA)                     :: LA
    INTEGER                                :: ks, ke
    INTEGER :: k, start, fin, nbghost, s, f
    Vec            :: xghost
    PetscErrorCode :: ierr
    PetscScalar, POINTER :: x_loc(:)
    CALL VecGetArrayF90(xghost, x_loc, ierr)
    DO k = ks, ke
       start = SUM(LA%dom_np(1:k-1)) + 1
       fin = start + LA%dom_np(k) - 1
       s = SUM(LA%np(ks:k-1)) + 1
       f = s + LA%dom_np(k) - 1
       phi(s:f) = x_loc(start:fin)
       nbghost =  LA%np(k) -  LA%dom_np(k)
       start = SUM(LA%dom_np) + SUM(LA%np(1:k-1)-LA%dom_np(1:k-1)) + 1
       fin =  start + nbghost - 1
       s = SUM(LA%np(ks:k-1))+LA%dom_np(k)+1
       f = s + nbghost - 1
       phi(s:f) = x_loc(start:fin)
    END DO
    CALL VecRestoreArrayF90(xghost, x_loc, ierr)
  END SUBROUTINE extract

  SUBROUTINE block_index(communicator,kmax,mesh,loc_to_glob_LA)
    USE def_type_mesh
    USE my_util
#include "petsc/finclude/petsc.h"
    use petsc
    IMPLICIT NONE
    TYPE(mesh_type)                           :: mesh
    INTEGER,                 INTENT(IN)       :: kmax
    INTEGER, DIMENSION(:,:), POINTER          :: loc_to_glob_LA
    INTEGER, DIMENSION(:),   POINTER          :: dom_np, disp
    INTEGER                                   :: code, nb_procs, rank
    INTEGER                                   :: i, p, n, k, i_loc, proc, iglob
    MPI_Comm       :: communicator

    CALL MPI_COMM_SIZE(communicator,nb_procs,code)
    CALL MPI_COMM_RANK(communicator,rank,code)
    ALLOCATE(dom_np(nb_procs) ,disp(nb_procs+1))
    CALL MPI_ALLGATHER(mesh%dom_np, 1, MPI_INTEGER, dom_np, 1, &
         MPI_INTEGER, communicator ,code)
    disp(1)=1
    DO n = 1, nb_procs
       disp(n+1) = disp(n) + dom_np(n)
    END DO
    IF (ASSOCIATED(mesh%disp)) THEN
       NULLIFY(mesh%disp)
    END IF
    IF (ASSOCIATED(mesh%domnp)) THEN
       NULLIFY(mesh%domnp)
    END IF
    ALLOCATE(mesh%disp(nb_procs+1))
    ALLOCATE(mesh%domnp(nb_procs))
    mesh%disp = disp
    mesh%domnp=dom_np

    ALLOCATE(loc_to_glob_LA(kmax,mesh%np))
    proc = rank + 1

    DO i = 1, mesh%dom_np
       DO k = 1, kmax
          loc_to_glob_LA(k,i) = kmax*(disp(proc)-1)+(k-1)*dom_np(proc)+i
       END DO
    END DO

!!$!TEST
!!$    DO i = 1, mesh%dom_np
!!$       iglob = mesh%loc_to_glob(i)
!!$       DO p = 2, nb_procs+1
!!$          IF (disp(p) > iglob) THEN
!!$             proc = p - 1
!!$             EXIT
!!$          END IF
!!$       END DO
!!$       IF (rank+1/=proc) THEN
!!$          write(*,*) 'BUG2', rank+1, proc
!!$          STOP
!!$       END IF
!!$
!!$       DO k = 2, kmax
!!$          IF (loc_to_glob_LA(k,i) - dom_np(proc) /= loc_to_glob_LA(k-1,i)) THEN
!!$             write(*,*) ' BUG1 ', rank
!!$             stop
!!$          END IF
!!$       END DO
!!$    END DO
!!$!TEST

    DO i = mesh%dom_np + 1, mesh%np
       iglob = mesh%loc_to_glob(i)
       DO p = 2, nb_procs+1
          IF (disp(p) > iglob) THEN
             proc = p - 1
             EXIT
          END IF
       END DO
       i_loc = iglob - disp(proc) + 1
       DO k = 1, kmax
          loc_to_glob_LA(k,i) = kmax*(disp(proc)-1)+(k-1)*dom_np(proc)+i_loc
       END DO
    END DO

!!$!TEST
!!$    DO i = 1, mesh%np
!!$       iglob = mesh%loc_to_glob(i)
!!$       DO p = 2, nb_procs+1
!!$          IF (disp(p) > iglob) THEN
!!$             proc = p - 1
!!$             EXIT
!!$          END IF
!!$       END DO
!!$       DO k = 2, kmax
!!$          IF (loc_to_glob_LA(k,i) - dom_np(proc) /= loc_to_glob_LA(k-1,i)) THEN
!!$             write(*,*) ' BUG ', rank
!!$             stop
!!$          END IF
!!$       END DO
!!$    END DO
!!$!TEST

    DEALLOCATE(dom_np,disp)

  END SUBROUTINE block_index

  SUBROUTINE st_aij_csr_glob_block(communicator,kmax,mesh_glob,mesh,LA, opt_per)
    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns
    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),         INTENT(IN)       :: mesh_glob,mesh
    INTEGER,                 INTENT(IN)       :: kmax
    TYPE(periodic_type), OPTIONAL, INTENT(IN) :: opt_per
    TYPE(petsc_csr_LA),       INTENT(OUT)     :: LA
    INTEGER :: nparm=200
    INTEGER :: me, nw, nmax, np, knp, ki, kj, k, njt, i1, i2
    INTEGER :: m, ni, nj, i, j, n_a_d, iloc, jloc, jglob, nb_procs, p, proc=-1
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    INTEGER, DIMENSION(SIZE(mesh_glob%jj,1))  :: jj_loc
    INTEGER, DIMENSION(:), ALLOCATABLE   :: per_loc
    LOGICAL :: out
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator

    CALL block_index(communicator,kmax,mesh,LA%loc_to_glob)
    nw = SIZE(mesh%jj, 1)
    me = mesh%dom_me
    np = mesh%dom_np
    knp = kmax*np
    nb_procs = SIZE(mesh%domnp)

    LA%kmax = kmax
    ALLOCATE(LA%dom_np(kmax),LA%np(kmax))
    LA%dom_np(:) = mesh%dom_np
    LA%np(:) = mesh%np

    IF (np==0) THEN
       ALLOCATE(LA%ia(0:0),LA%ja(0))
       LA%ia(0) = 0
       RETURN
    END IF

    ALLOCATE (ja_work(knp,kmax*nparm), a_d(kmax*nparm), nja(knp))
    ALLOCATE (per_loc(knp))
    ja_work = 0
    nja = 1
    DO ki = 1, kmax
       DO i = 1, np
          ja_work((ki-1)*np+i,1) = LA%loc_to_glob(ki,i)
       END DO
    END DO

    DO m = 1, mesh_glob%me
       jj_loc = mesh_glob%jj(:,m)
       IF (MAXVAL(jj_loc)< mesh%loc_to_glob(1) .OR. MINVAL(jj_loc)> mesh%loc_to_glob(1) + mesh%np -1) CYCLE
       DO ni = 1, nw
          iloc = jj_loc(ni)-mesh%loc_to_glob(1)+1
          IF (iloc<1 .OR. iloc>np) CYCLE
          DO ki = 1, kmax
             i = iloc + (ki-1)*np
             DO nj = 1, nw
                jglob = jj_loc(nj)
                IF (jglob< mesh%loc_to_glob(1) .OR. jglob> mesh%loc_to_glob(2)) THEN
                   DO p = 2, nb_procs+1
                      IF (mesh%disp(p) > jglob) THEN
                         proc = p - 1
                         EXIT
                      END IF
                   END DO
                   out = .TRUE.
                   jloc = jglob - mesh%disp(proc) + 1
                ELSE
                   out = .FALSE.
                   jloc = jglob - mesh%loc_to_glob(1) + 1
                END IF
                DO kj = 1, kmax
                   IF (out) THEN
                      j = kmax*(mesh%disp(proc)-1)+(kj-1)*mesh%domnp(proc)+jloc
                   ELSE
                      j = LA%loc_to_glob(kj,jloc)
                   END IF

                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
                      nja(i) = nja(i) + 1
                      ja_work(i,nja(i)) = j
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO

    IF (PRESENT(opt_per)) THEN
       DO k = 1, opt_per%n_bord
          DO i = 1, SIZE(opt_per%list(k)%DIL)
             per_loc=0
             i1 = opt_per%list(k)%DIL(i)
             i2 = opt_per%perlist(k)%DIL(i)
             njt = nja(i1)+nja(i2)
             IF (njt > kmax*nparm) THEN
                CALL error_Petsc('BUG in st_aij_glob_block, SIZE(ja) not large enough')
             END IF
             per_loc(1:nja(i1)) = ja_work(i1,1:nja(i1))
             per_loc(nja(i1)+1:nja(i1)+nja(i2)) = ja_work(i2,1:nja(i2))
             nja(i1) = njt
             nja(i2) = njt
             ja_work(i1,1:njt) = per_loc(1:njt)
             ja_work(i2,1:njt) = per_loc(1:njt)
          END DO
       END DO
    END IF

    IF (MAXVAL(nja)>nparm) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, knp
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(LA%ia(0:knp),LA%ja(0:nmax-1))
    LA%ia(0) = 0
    DO i = 1, knp
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
          STOP
       END IF
       LA%ia(i) = LA%ia(i-1) + nja(i)
       LA%ja(LA%ia(i-1):LA%ia(i)-1) = a_d(1:nja(i))-1
    END DO
    DEALLOCATE (ja_work, nja, a_d )
    DEALLOCATE (per_loc)
  END SUBROUTINE st_aij_csr_glob_block

!!$  SUBROUTINE st_aij_csr_glob_block_without_per(communicator,kmax,mesh_glob,mesh,LA)
!!$    !  input coefficient structure of the matrix and
!!$    !  perform the ordering of the unknowns
!!$    !  jj(nodes_per_element, number_of_elements)
!!$    !                  --->  node number in the grid
!!$    USE def_type_mesh
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),         INTENT(IN)  :: mesh_glob,mesh
!!$    INTEGER,                 INTENT(IN)  :: kmax
!!$    TYPE(petsc_csr_LA),       INTENT(OUT) :: LA
!!$    INTEGER :: nparm=200
!!$    INTEGER :: me, nw, nmax, np, knp, ki, kj
!!$    INTEGER :: m, ni, nj, i, j, n_a_d, iloc, jloc, jglob, nb_procs, p, proc=-1
!!$    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
!!$    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
!!$    INTEGER, DIMENSION(SIZE(mesh_glob%jj,1))  :: jj_loc
!!$    LOGICAL :: out
!!$    !#include "petsc/finclude/petsc.h"
!!$    MPI_Comm       :: communicator
!!$
!!$    CALL block_index(communicator,kmax,mesh,LA%loc_to_glob)
!!$    nw = SIZE(mesh%jj, 1)
!!$    me = mesh%dom_me
!!$    np = mesh%dom_np
!!$    knp = kmax*np
!!$    nb_procs = SIZE(mesh%domnp)
!!$
!!$    LA%kmax = kmax
!!$    ALLOCATE(LA%dom_np(kmax),LA%np(kmax))
!!$    LA%dom_np(:) = mesh%dom_np
!!$    LA%np(:) = mesh%np
!!$
!!$    IF (np==0) THEN
!!$       ALLOCATE(LA%ia(0:0),LA%ja(0))
!!$       LA%ia(0) = 0
!!$       RETURN
!!$    END IF
!!$
!!$    ALLOCATE (ja_work(knp,kmax*nparm), a_d(kmax*nparm), nja(knp))
!!$    ja_work = 0
!!$    nja = 1
!!$    DO ki = 1, kmax
!!$       DO i = 1, np
!!$          ja_work((ki-1)*np+i,1) = LA%loc_to_glob(ki,i)
!!$       END DO
!!$    END DO
!!$
!!$    DO m = 1, mesh_glob%me
!!$       jj_loc = mesh_glob%jj(:,m)
!!$       IF (MAXVAL(jj_loc)< mesh%loc_to_glob(1) .OR. MINVAL(jj_loc)> mesh%loc_to_glob(1) + mesh%np -1) CYCLE
!!$       DO ni = 1, nw
!!$          iloc = jj_loc(ni)-mesh%loc_to_glob(1)+1
!!$          IF (iloc<1 .OR. iloc>np) CYCLE
!!$          DO ki = 1, kmax
!!$             i = iloc + (ki-1)*np
!!$             DO nj = 1, nw
!!$                jglob = jj_loc(nj)
!!$                IF (jglob< mesh%loc_to_glob(1) .OR. jglob> mesh%loc_to_glob(2)) THEN
!!$                   DO p = 2, nb_procs+1
!!$                      IF (mesh%disp(p) > jglob) THEN
!!$                         proc = p - 1
!!$                         EXIT
!!$                      END IF
!!$                   END DO
!!$                   out = .TRUE.
!!$                   jloc = jglob - mesh%disp(proc) + 1
!!$                ELSE
!!$                   out = .FALSE.
!!$                   jloc = jglob - mesh%loc_to_glob(1) + 1
!!$                END IF
!!$                DO kj = 1, kmax
!!$                   IF (out) THEN
!!$                      j = kmax*(mesh%disp(proc)-1)+(kj-1)*mesh%domnp(proc)+jloc
!!$                   ELSE
!!$                      j = LA%loc_to_glob(kj,jloc)
!!$                   END IF
!!$
!!$                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
!!$                      nja(i) = nja(i) + 1
!!$                      ja_work(i,nja(i)) = j
!!$                   END IF
!!$                END DO
!!$             END DO
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF (MAXVAL(nja)>nparm) THEN
!!$       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
!!$       STOP
!!$    END IF
!!$
!!$    nmax = 0
!!$    DO i = 1, knp
!!$       nmax = nmax + nja(i)
!!$    END DO
!!$    ALLOCATE(LA%ia(0:knp),LA%ja(0:nmax-1))
!!$    LA%ia(0) = 0
!!$    DO i = 1, knp
!!$       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
!!$       IF (n_a_d /= nja(i)) THEN
!!$          WRITE(*,*) ' BUG : st_p1_CSR'
!!$          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
!!$          STOP
!!$       END IF
!!$       LA%ia(i) = LA%ia(i-1) + nja(i)
!!$       LA%ja(LA%ia(i-1):LA%ia(i)-1) = a_d(1:nja(i))-1
!!$    END DO
!!$    DEALLOCATE (ja_work, nja, a_d )
!!$  END SUBROUTINE st_aij_csr_glob_block_without_per

  SUBROUTINE st_aij_csr_loc_block(communicator,kmax,mesh,LA)
    !  input coefficient structure of the matrix and`
    !  perform the ordering of the unknowns
    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),         INTENT(IN)  :: mesh
    INTEGER,                 INTENT(IN)  :: kmax
    TYPE(petsc_csr_LA),       INTENT(OUT) :: LA
    INTEGER :: nparm=90
    INTEGER :: me, nw, nmax, np, knp, ki, kj
    INTEGER :: m, ni, nj, i, j, n_a_d, iloc
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
    INTEGER, DIMENSION(SIZE(mesh%jj,1))  :: j_loc
    !#include "petsc/finclude/petsc.h"
    MPI_Comm       :: communicator

    CALL block_index(communicator,kmax,mesh,LA%loc_to_glob)
    nw = SIZE(mesh%jj, 1)
    me = mesh%dom_me
    np = mesh%dom_np
    knp = kmax*np
    ALLOCATE (ja_work(knp,kmax*nparm), a_d(kmax*nparm), nja(knp))
    ja_work = 0
    nja = 1
    DO ki = 1, kmax
       DO i = 1, np
          ja_work((ki-1)*np+i,1) = LA%loc_to_glob(ki,i)
       END DO
    END DO

    DO m = 1, me
       j_loc = mesh%jj(:,m)
       DO ki = 1, kmax
          DO ni = 1, nw
             iloc = j_loc(ni)
             IF (iloc>np) CYCLE
             i = iloc + (ki-1)*np
             DO kj = 1, kmax
                DO nj = 1, nw
                   j = LA%loc_to_glob(kj,j_loc(nj))
                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
                      nja(i) = nja(i) + 1
                      ja_work(i,nja(i)) = j
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO

    IF (MAXVAL(nja)>nparm) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, knp
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(LA%ia(0:knp),LA%ja(0:nmax-1))
    LA%ia(0) = 0
    DO i = 1, knp
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
          STOP
       END IF
       LA%ia(i) = LA%ia(i-1) + nja(i)
       LA%ja(LA%ia(i-1):LA%ia(i)-1) = a_d(1:nja(i))-1
    END DO
    DEALLOCATE (ja_work, nja, a_d )
  END SUBROUTINE st_aij_csr_loc_block

  SUBROUTINE st_aij_csr(jj, aij)
    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns
    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(aij_type),          INTENT(OUT) :: aij
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
    INTEGER :: nparm=90
    INTEGER :: me, nw, nmax, np
    INTEGER :: m, ni, nj, i, j, n_a_d
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)
    ALLOCATE (ja_work(np,nparm), aij%ia(np+1), a_d(nparm), nja(np))
    ja_work = 0
    nja = 1
    DO i = 1, np
       ja_work(i,1) = i
    END DO

    DO m = 1, me
       DO ni = 1, nw
          i = jj(ni,m)
          DO nj = 1, nw
             j = jj(nj,m)
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
                nja(i) = nja(i) + 1
                ja_work(i,nja(i)) = j
             END IF
          END DO
       END DO
    END DO

    IF (MAXVAL(nja)>nparm) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, np
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(aij%ja(nmax))
    aij%ia(1) = 1
    DO i = 1, np
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
          STOP
       END IF
       aij%ia(i+1) = aij%ia(i) + nja(i)
       aij%ja(aij%ia(i):aij%ia(i+1)-1) = a_d(1:nja(i))
    END DO
    DEALLOCATE ( ja_work, nja, a_d )
  END SUBROUTINE st_aij_csr

  SUBROUTINE st_csr_block(n_b,in,out)
    !SUBROUTINE st_csr_block(in,vv_rt_loc_to_glob_LA,out)
    !
    ! Builds the CSR structure of a parallel block matrix
    ! from its scalar counterpart.
    !
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(aij_type),          INTENT(IN)  :: in
    TYPE(aij_type),          INTENT(OUT) :: out
    INTEGER,                 INTENT(IN)  :: n_b  ! Number of blocs
    INTEGER :: ki, kj, i, ib, jb, p, pb, bloc_size, np, nnz

    np = SIZE(in%ia) - 1
    nnz = in%ia(np+1) - in%ia(1)
    ALLOCATE (out%ia(n_b*np+1), out%ja(n_b*n_b*nnz))

    bloc_size = SIZE(in%ia) - 1

    out%ia(1) = in%ia(1)

    DO ki = 1, n_b
       DO i = 2, bloc_size + 1
          ib = i + (ki-1)*bloc_size
          out%ia(ib) = out%ia(ib-1) + n_b*(in%ia(i) - in%ia(i-1))
       END DO
    END DO

    DO ki = 1, n_b
       DO i = 1, bloc_size
          ib = i + (ki-1)*bloc_size

          DO kj = 1,  n_b
             DO p = in%ia(i), in%ia(i+1) - 1
                jb = in%ja(p) + (kj-1)*bloc_size

                pb = out%ia(ib)  +  p - in%ia(i)  +  (kj-1)*(in%ia(i+1)-in%ia(i))
                out%ja(pb) = jb

             END DO
          END DO

       END DO
    END DO
  END SUBROUTINE ST_CSR_BLOCK

  SUBROUTINE st_csr_bloc(ia,ja,ia_b,ja_b,n_b)
    !
    ! Builds the CSR structure of a 3D vector matrix
    ! from its scalar counterpart.
    !
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN)  :: ia, ja
    INTEGER, DIMENSION(:), POINTER     :: ia_b, ja_b
    INTEGER,               INTENT(IN)  :: n_b  ! Number of blocs
    INTEGER :: ki, kj, i, ib, jb, p, pb, bloc_size, np, nnz

    np = SIZE(ia) - 1
    nnz = ia(np+1) - ia(1)
    ALLOCATE (ia_b(n_b*np+1), ja_b(n_b*n_b*nnz))

    bloc_size = SIZE(ia) - 1

    ia_b(1) = ia(1)

    DO ki = 1, n_b
       DO i = 2, bloc_size + 1
          ib = i + (ki-1)*bloc_size
          ia_b(ib) = ia_b(ib-1) + n_b*(ia(i) - ia(i-1))
       END DO
    END DO

    DO ki = 1, n_b
       DO i = 1, bloc_size
          ib = i + (ki-1)*bloc_size

          DO kj = 1,  n_b
             DO p = ia(i), ia(i+1) - 1
                jb = ja(p) + (kj-1)*bloc_size

                pb = ia_b(ib)  +  p - ia(i)  +  (kj-1)*(ia(i+1)-ia(i))
                ja_b(pb) = jb

             END DO
          END DO
       END DO
    END DO

  END SUBROUTINE st_csr_bloc

  SUBROUTINE st_csr(jj, ia, ja)
    !  input coefficient structure of the matrix and
    !  perform the ordering of the unknowns
    !  jj(nodes_per_element, number_of_elements)
    !                  --->  node number in the grid
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), INTENT(IN)  :: jj
    INTEGER, DIMENSION(:),   POINTER     :: ia, ja
    INTEGER :: nparm=90
    INTEGER :: me, nw, nmax, np
    INTEGER :: m, ni, nj, i, j, n_a_d
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d

    nw = SIZE(jj, 1)
    me = SIZE(jj, 2)
    np = MAXVAL(jj)

    ALLOCATE (ja_work(np,nparm), ia(np+1), a_d(nparm), nja(np))

    ja_work = 0
    nja = 1
    DO i = 1, np
       ja_work(i,1) = i
    END DO

    DO m = 1, me
       DO ni = 1, nw
          i = jj(ni,m)

          DO nj = 1, nw
             j = jj(nj,m)
             IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
                nja(i) = nja(i) + 1
                ja_work(i,nja(i)) = j
             END IF
          END DO


       END DO
    END DO

    IF (MAXVAL(nja)>nparm) THEN
       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
       STOP
    END IF

    nmax = 0
    DO i = 1, np
       nmax = nmax + nja(i)
    END DO
    ALLOCATE(ja(nmax))

    ia(1) = 1
    DO i = 1, np
       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
       IF (n_a_d /= nja(i)) THEN
          WRITE(*,*) ' BUG : st_p1_CSR'
          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
          STOP
       END IF
       ia(i+1) = ia(i) + nja(i)
       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
    END DO


    DEALLOCATE ( ja_work, nja, a_d )

  END SUBROUTINE st_csr

!!$  SUBROUTINE st_csr_edge_stab(jji, ia, ja)
!!$    !  input coefficient structure of the matrix and
!!$    !  perform the ordering of the unknowns
!!$    !  jj(nodes_per_element, number_of_elements)
!!$    !                  --->  node number in the grid
!!$    IMPLICIT NONE
!!$    INTEGER, DIMENSION(:,:,:), INTENT(IN)  :: jji
!!$    INTEGER, DIMENSION(:),   POINTER     :: ia, ja
!!$    INTEGER :: nparm=180
!!$    INTEGER :: mi, nw, nmax, np
!!$    INTEGER :: ms, ni, nj, i, j, n_a_d, cotei, cotej
!!$    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ja_work
!!$    INTEGER, DIMENSION(:),   ALLOCATABLE :: nja, a_d
!!$
!!$    nw = SIZE(jji, 1)
!!$    mi = SIZE(jji, 3)
!!$    np = MAXVAL(jji)
!!$
!!$    ALLOCATE (ja_work(np,nparm), ia(np+1), a_d(nparm), nja(np))
!!$
!!$    ja_work = 0
!!$    nja = 1
!!$    DO i = 1, np
!!$       ja_work(i,1) = i
!!$    END DO
!!$
!!$    DO ms = 1, mi
!!$       DO cotei = 1, 2
!!$          DO ni = 1, nw
!!$             i = jji(ni,cotei,ms)
!!$             DO cotej = 1, 2
!!$                DO nj = 1, nw
!!$                   j = jji(nj,cotej,ms)
!!$                   IF (MINVAL(ABS(ja_work(i,1:nja(i))-j)) /= 0) THEN
!!$                      nja(i) = nja(i) + 1
!!$                      ja_work(i,nja(i)) = j
!!$                   END IF
!!$                END DO
!!$             END DO
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    IF (MAXVAL(nja)>nparm) THEN
!!$       WRITE(*,*) 'ST_SPARSEKIT: dimension de ja doit etre >= ',nparm
!!$       STOP
!!$    END IF
!!$
!!$    nmax = 0
!!$    DO i = 1, np
!!$       nmax = nmax + nja(i)
!!$    END DO
!!$    ALLOCATE(ja(nmax))
!!$
!!$    ia(1) = 1
!!$    DO i = 1, np
!!$       CALL tri_jlg (ja_work(i,1:nja(i)), a_d, n_a_d)
!!$       IF (n_a_d /= nja(i)) THEN
!!$          WRITE(*,*) ' BUG : st_p1_CSR'
!!$          WRITE(*,*) 'n_a_d ', n_a_d, 'nja(i)', nja(i)
!!$          STOP
!!$       END IF
!!$       ia(i+1) = ia(i) + nja(i)
!!$       ja(ia(i):ia(i+1)-1) = a_d(1:nja(i))
!!$    END DO
!!$
!!$    DEALLOCATE ( ja_work, nja, a_d )
!!$
!!$  END SUBROUTINE st_csr_edge_stab

  SUBROUTINE tri_jlg (a,  a_d, n_a_d)
    !  sort in ascending order of the integer array  a  and generation
    !  of the integer array  a_d  whose first  n_a_d  leading entries
    !  contain different values in ascending order, while all the
    !  remaining entries are set to zero
    !  sorting by Shell's method.
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(INOUT) :: a
    INTEGER, DIMENSION(:), INTENT(OUT)   :: a_d
    INTEGER,               INTENT(OUT)   :: n_a_d
    INTEGER :: n, na, inc, i, j, k, ia

    na = SIZE(a)

    !  sort phase

    IF (na == 0) THEN
       n_a_d = 0
       RETURN
    ENDIF

    inc = 1
    DO WHILE (inc <= na)
       inc = inc * 3
       inc = inc + 1
    ENDDO

    DO WHILE (inc > 1)
       inc = inc/3
       DO i = inc + 1, na
          ia = a(i)
          j = i
          DO WHILE (a(j-inc) > ia)
             a(j) = a(j-inc)
             j = j - inc
             IF (j <= inc) EXIT
          ENDDO
          a(j) = ia
       ENDDO
    ENDDO

    !  compression phase

    n = 1
    a_d(n) = a(1)
    DO k = 2, na
       IF (a(k) > a(k-1)) THEN
          n = n + 1
          a_d(n) = a(k)
          !TEST JUIN 13 2008
       ELSE
          WRITE(*,*) 'We have a problem in the compression phase of tri_jlg', a(k), a(k-1)
          !TEST JUIN 13 2008
       ENDIF
    ENDDO

    n_a_d = n

    a_d(n_a_d + 1 : na) = 0

  END SUBROUTINE tri_jlg

END MODULE st_matrix
