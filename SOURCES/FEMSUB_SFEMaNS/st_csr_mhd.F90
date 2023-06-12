MODULE st_csr_mhd
  PUBLIC :: st_scr_maxwell_mu_H_p_phi

  PRIVATE
#include "petsc/finclude/petsc.h"

CONTAINS

  SUBROUTINE st_scr_maxwell_mu_H_p_phi(communicator, H_mesh_glob, H_mesh, pmag_mesh_glob, pmag_mesh, &
       phi_mesh_glob, phi_mesh, interface_glob, interface_H_mu_glob, &
       LA_H, LA_pmag, LA_phi, LA_mhd, opt_per)
    USE my_util
    USE def_type_mesh
    USE st_matrix
    IMPLICIT NONE
    TYPE(mesh_type),           INTENT(IN) :: H_mesh_glob, pmag_mesh_glob, phi_mesh_glob
    TYPE(mesh_type),           INTENT(IN) :: H_mesh, pmag_mesh, phi_mesh
    TYPE(interface_type),      INTENT(IN) :: interface_glob, interface_H_mu_glob
    TYPE(periodic_type), OPTIONAL, INTENT(IN) :: opt_per
    TYPE(petsc_csr_LA),        INTENT(OUT):: LA_H, LA_pmag, LA_phi, LA_mhd
    INTEGER, PARAMETER                    :: param=200
    INTEGER, DIMENSION(param)             :: a_d
    INTEGER, DIMENSION(SIZE(H_mesh_glob%jj,1)) :: jj_loc, jj_loc1, jj_loc2
    INTEGER, DIMENSION(:,:),   POINTER    :: ja_work
    INTEGER, DIMENSION(:),     POINTER    :: nja_glob
    INTEGER                               :: np_glob, np_H, np_pmag, np_phi, np_m, code, rank, dp, n_a_d, &
         iglob, jglob, iloc, jloc, ni, nj, ci, m, mi, mj, m1, m2, i ,j, k, ki, kj, &
         jmin, jmax, p, nnz, p_i, p_f, ms, nb_procs, proc
    LOGICAL                               :: out
    INTEGER, DIMENSION(:), ALLOCATABLE    :: per_loc
    INTEGER                               :: njt, i1, i2
    MPI_Comm       :: communicator

    CALL MPI_COMM_RANK(communicator,rank,code)
    CALL MPI_COMM_SIZE(communicator,nb_procs,code)


    CALL st_aij_csr_glob_block(communicator,3,H_mesh_glob,H_mesh,LA_H)
    CALL st_aij_csr_glob_block(communicator,1,pmag_mesh_glob,pmag_mesh,LA_pmag)
    CALL st_aij_csr_glob_block(communicator,1,phi_mesh_glob,phi_mesh,LA_phi)

    np_H = 3*H_mesh%dom_np
    np_pmag = pmag_mesh%dom_np
    np_phi = phi_mesh%dom_np
    np_glob = np_H + np_pmag + np_phi
    np_m = H_mesh%dom_np


    ! Modify LA_H, LA_pmag, LA_phi
    dp = 3*(H_mesh%disp(rank+1)-1)+(pmag_mesh%disp(rank+1)-1)+(phi_mesh%disp(rank+1)-1)

!!$    WRITE(*,*) rank, 'H deb', MINVAL(LA_H%loc_to_glob(1,1:H_mesh%dom_np)), MAXVAL(LA_H%loc_to_glob(3,1:H_mesh%dom_np))
!!$    WRITE(*,*) rank, 'pmag deb', MINVAL(LA_pmag%loc_to_glob(1,1:pmag_mesh%dom_np)), MAXVAL(LA_pmag%loc_to_glob(1,1:pmag_mesh%dom_np))
!!$    WRITE(*,*) rank, 'phi deb', MINVAL(LA_phi%loc_to_glob(1,1:phi_mesh%dom_np)), MAXVAL(LA_phi%loc_to_glob(1,1:phi_mesh%dom_np))

!!$    DO i = 1, H_mesh%dom_np
!!$       DO k = 1, 3
!!$          LA_H%loc_to_glob(k,i) = LA_H%loc_to_glob(k,i)-3*(H_mesh%disp(rank+1)-1) + dp !+ (k-1)*np_m
!!$       END DO
!!$    END DO
!!$    DO i = 1, pmag_mesh%dom_np
!!$       LA_pmag%loc_to_glob(1,i) = LA_pmag%loc_to_glob(1,i) - (pmag_mesh%disp(rank+1)-1) +dp + np_H
!!$    END DO
!!$    DO i = 1, phi_mesh%dom_np
!!$       LA_phi%loc_to_glob(1,i) = LA_phi%loc_to_glob(1,i) - (phi_mesh%disp(rank+1)-1) + dp + np_H + np_pmag
!!$    END DO


    DO i = 1, H_mesh%np
       DO k = 1, 3
          CALL search_index_block(H_mesh, LA_H%loc_to_glob(k,i),rank, 3,  nb_procs, iloc, proc, out)
          LA_H%loc_to_glob(k,i) = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1)+ iloc
       END DO
    END DO

    DO i = 1, pmag_mesh%np
       CALL search_index_block(pmag_mesh, LA_pmag%loc_to_glob(1,i),rank, 1, nb_procs, iloc, proc, out)
       LA_pmag%loc_to_glob(1,i) =3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
            + 3*H_mesh%domnp(proc) + iloc
    END DO
    DO i = 1, phi_mesh%np
       CALL search_index_block(phi_mesh, LA_phi%loc_to_glob(1,i), rank, 1,  nb_procs, iloc, proc, out)
       LA_phi%loc_to_glob(1,i) =3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1)  &
            + 3*H_mesh%domnp(proc) +pmag_mesh%domnp(proc) + iloc
    END DO

    ! End Modify LA_H, LA_pmag, LA_phi

    ALLOCATE(ja_work(np_glob,param),nja_glob(np_glob))
    ALLOCATE(per_loc(2*param))
    nja_glob = 0

    ! Block HxH
    IF (H_mesh%me /=0) THEN
       DO m = 1, H_mesh_glob%me
          DO ni = 1, SIZE(H_mesh%jj,1)
             iglob = H_mesh_glob%jj(ni,m)
             IF (iglob<H_mesh%loc_to_glob(1) .OR. iglob>H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) CYCLE
             DO nj = 1, SIZE(H_mesh%jj,1)
                jglob = H_mesh_glob%jj(nj,m)
                CALL search_index(H_mesh,jglob,nb_procs,jloc,proc,out)
                DO kj = 1, 3
                   IF (out) THEN
                      j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) + (kj-1)*H_mesh%domnp(proc) &
                           + jloc
                   ELSE
                      j = LA_H%loc_to_glob(kj,jloc)
                   END IF
                   DO ki = 1, 3
                      i = iglob - H_mesh%loc_to_glob(1)+1 + (ki-1)*np_m
                      IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) /= 0) THEN
                         nja_glob(i) = nja_glob(i) + 1
                         ja_work(i,nja_glob(i)) = j
                      END IF
                   END DO
                END DO

             END DO
          END DO
       END DO
    END IF
    ! Block HxH

    ! Block pxp
    DO i = 1, np_pmag
       iglob =  i + np_H
       nnz = LA_pmag%ia(i) - LA_pmag%ia(i-1) !LA is C indexed
       p_i = LA_pmag%ia(i-1)
       p_f = LA_pmag%ia(i)-1
       DO p = p_i, p_f
          jglob = LA_pmag%ja(p)+1
          CALL search_index_block(pmag_mesh, jglob, rank, 1, nb_procs,jloc,proc,out)
          IF (out) THEN
             j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) + 3*H_mesh%domnp(proc) + jloc
          ELSE
             j = LA_pmag%loc_to_glob(1,jloc)
          END IF
          ja_work(iglob,p-p_i+1) = j
       END DO
       nja_glob(iglob) = nnz
    END DO
    ! End Block pxp

    ! Block phixphi
    DO i = 1, np_phi
       iglob =  i + np_H + np_pmag
       nnz = LA_phi%ia(i) - LA_phi%ia(i-1) !LA is C indexed
       p_i = LA_phi%ia(i-1)
       p_f = LA_phi%ia(i)-1
       DO p = p_i, p_f
          jglob = LA_phi%ja(p)+1
          CALL search_index_block(phi_mesh,jglob, rank, 1, nb_procs,jloc,proc,out)
          IF (out) THEN
             j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
                  + 3*H_mesh%domnp(proc) + pmag_mesh%domnp(proc) + jloc
          ELSE
             j = LA_phi%loc_to_glob(1,jloc)
          END IF
          ja_work(iglob,p-p_i+1) = j
       END DO
       nja_glob(iglob) = nnz
    END DO
    ! End Block phixphi

    ! Cleanup
    NULLIFY(LA_H%ia,LA_H%ja,LA_pmag%ia,LA_pmag%ja,LA_phi%ia,LA_phi%ja)
    ! Cleanup

    ! Block Hxp and pxH
    IF (H_mesh%me /=0) THEN
       DO m = 1, H_mesh_glob%me
          jj_loc = H_mesh_glob%jj(:,m)
          IF (MAXVAL(jj_loc)<H_mesh%loc_to_glob(1) .OR. MINVAL(jj_loc)>H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) CYCLE
          DO ni = 1, SIZE(H_mesh%jj,1)
             iglob = jj_loc(ni)

             DO nj = 1, SIZE(pmag_mesh%jj,1)
                jglob = pmag_mesh_glob%jj(nj,m)
                CALL search_index(pmag_mesh,jglob,nb_procs,jloc,proc,out)
                IF (out) THEN
                   j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) + 3*H_mesh%domnp(proc) + jloc
                ELSE
                   j = LA_pmag%loc_to_glob(1,jloc)
                END IF

                IF (iglob.GE.H_mesh%loc_to_glob(1) .AND. iglob.LE.H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) THEN
                   i = iglob - H_mesh%loc_to_glob(1)+1
                   IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) /= 0) THEN
                      nja_glob(i) = nja_glob(i) + 1
                      nja_glob(i+np_m) = nja_glob(i+np_m) + 1
                      nja_glob(i+2*np_m) = nja_glob(i+2*np_m) + 1
                      ja_work(i,nja_glob(i)) = j
                      ja_work(i+np_m,nja_glob(i+np_m)) = j
                      ja_work(i+2*np_m,nja_glob(i+2*np_m)) = j
                   END IF
                END IF

                IF (.NOT.out) THEN
                   i = np_H + jloc
                   CALL search_index(H_mesh,iglob,nb_procs,iloc,proc,out)
                   DO k = 1, 3
                      IF (out) THEN
                         j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
                              + (k-1)*H_mesh%domnp(proc) + iloc
                      ELSE
                         j = LA_H%loc_to_glob(k,iloc)
                      END IF
                      IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) /= 0) THEN
                         nja_glob(i) = nja_glob(i) + 1
                         ja_work(i,nja_glob(i)) = j
                      END IF
                   END DO
                END IF

             END DO
          END DO
       END DO
    END IF
    ! End Block Hxp and pxH

    ! Interface_H_mu
    IF (H_mesh%me /=0) THEN
       DO ms = 1, interface_H_mu_glob%mes
          m1 = H_mesh_glob%neighs(interface_H_mu_glob%mesh1(ms))
          m2 = H_mesh_glob%neighs(interface_H_mu_glob%mesh2(ms))

          jj_loc1 = H_mesh_glob%jj(:,m1)
          jj_loc2 = H_mesh_glob%jj(:,m2)
          jmin = MIN(MINVAL(jj_loc1),MINVAL(jj_loc2))
          jmax = MAX(MAXVAL(jj_loc1),MAXVAL(jj_loc2))

          IF (jmax<H_mesh%loc_to_glob(1) .OR. jmin>H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) CYCLE

          DO ci = 1, 2
             IF (ci==1) THEN
                mi = m1
                mj = m2
             ELSE
                mi = m2
                mj = m1
             END IF

             DO ni = 1, SIZE(H_mesh%jj,1)
                iglob = H_mesh_glob%jj(ni,mi)
                IF (iglob < H_mesh%loc_to_glob(1) .OR. iglob > H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) CYCLE

                DO nj = 1, SIZE(H_mesh%jj,1)
                   jglob = H_mesh_glob%jj(nj,mj)
                   CALL search_index(H_mesh,jglob,nb_procs,jloc,proc,out)
                   DO kj = 1, 3
                      IF (out) THEN
                         j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
                              + (kj-1)*H_mesh%domnp(proc) + jloc
                      ELSE
                         j = LA_H%loc_to_glob(kj,jloc)
                      END IF
                      DO ki = 1, 3
                         i = iglob - H_mesh%loc_to_glob(1) + 1 + (ki-1)*np_m
                         IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) == 0) CYCLE
                         nja_glob(i) = nja_glob(i) + 1
                         ja_work(i,nja_glob(i)) = j
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END IF
    ! End Interface_H_mu

    ! Interface_H_phi
    IF (H_mesh%me*phi_mesh%me /=0) THEN
       DO ms = 1, interface_glob%mes
          m1 = H_mesh_glob%neighs(interface_glob%mesh1(ms))
          m2 = phi_mesh_glob%neighs(interface_glob%mesh2(ms))

          DO ni = 1, SIZE(H_mesh%jj,1)
             iglob = H_mesh_glob%jj(ni,m1)
             IF (iglob < H_mesh%loc_to_glob(1) .OR. iglob > H_mesh%loc_to_glob(1) + H_mesh%dom_np -1) CYCLE

             DO nj = 1, SIZE(phi_mesh%jj,1)
                jglob =  phi_mesh_glob%jj(nj,m2)
                CALL search_index(phi_mesh,jglob,nb_procs,jloc,proc,out)
                IF (out) THEN
                   j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
                        + 3*H_mesh%domnp(proc) + pmag_mesh%domnp(proc) + jloc
                ELSE
                   j = LA_phi%loc_to_glob(1,jloc)
                END IF
                DO ki = 1, 3
                   i = iglob - H_mesh%loc_to_glob(1) + 1 + (ki-1)*np_m
                   IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) == 0) CYCLE
                   nja_glob(i) = nja_glob(i) + 1
                   ja_work(i,nja_glob(i)) = j
                END DO
             END DO
          END DO

          DO ni = 1, SIZE(phi_mesh%jj,1)
             iglob = phi_mesh_glob%jj(ni,m2)
             IF (iglob < phi_mesh%loc_to_glob(1) .OR. iglob > phi_mesh%loc_to_glob(1) + phi_mesh%dom_np -1) CYCLE
             i = iglob - phi_mesh%loc_to_glob(1) + 1 + np_H + np_pmag
             DO nj = 1, SIZE(H_mesh%jj,1)
                jglob = H_mesh_glob%jj(nj,m1)
                CALL search_index(H_mesh,jglob,nb_procs,jloc,proc,out)
                DO kj = 1, 3
                   IF (out) THEN
                      j = 3*(H_mesh%disp(proc)-1)+(pmag_mesh%disp(proc)-1)+(phi_mesh%disp(proc)-1) &
                           + (kj-1)*H_mesh%domnp(proc) + jloc
                   ELSE
                      j = LA_H%loc_to_glob(kj,jloc)
                   END IF
                   IF (MINVAL(ABS(ja_work(i,1:nja_glob(i))-j)) == 0) CYCLE
                   nja_glob(i) = nja_glob(i) + 1
                   ja_work(i,nja_glob(i)) = j
                END DO
             END DO
          END DO
       END DO
    END IF
    ! End Interface_H_phi

    IF (PRESENT(opt_per)) THEN
       DO k = 1, opt_per%n_bord
          DO i = 1, SIZE(opt_per%list(k)%DIL)
             per_loc=0
             i1 = opt_per%list(k)%DIL(i)
             i2 = opt_per%perlist(k)%DIL(i)
             njt = nja_glob(i1)+nja_glob(i2)
             IF (njt > param) THEN
                CALL error_Petsc('BUG in st_aij_glob_block, SIZE(ja) not large enough')
             END IF
             per_loc(1:nja_glob(i1)) = ja_work(i1,1:nja_glob(i1))
             per_loc(nja_glob(i1)+1:njt) = ja_work(i2,1:nja_glob(i2))
             nja_glob(i1) = njt
             nja_glob(i2) = njt
             ja_work(i1,1:njt) = per_loc(1:njt)
             ja_work(i2,1:njt) = per_loc(1:njt)
          END DO
       END DO
    END IF



    IF (MAXVAL(nja_glob)>param) THEN
       CALL error_Petsc('ST_SPARSEKIT: dimension of ja must be >= nparm')
    END IF


    nnz = SUM(nja_glob)
    ALLOCATE(LA_mhd%ia(0:np_glob),LA_mhd%ja(0:nnz-1))
    LA_mhd%ia(0) = 0
    DO i = 1, np_glob
       CALL tri_jlg (ja_work(i,1:nja_glob(i)), a_d, n_a_d)
       IF (n_a_d /= nja_glob(i)) THEN
          CALL error_Petsc(' BUG in st_scr_maxwell_mu_H_p_phi')
       END IF
       LA_mhd%ia(i) = LA_mhd%ia(i-1) + nja_glob(i)
       LA_mhd%ja(LA_mhd%ia(i-1):LA_mhd%ia(i)-1) = a_d(1:nja_glob(i))-1
    END DO

    ALLOCATE(LA_mhd%dom_np(5),LA_mhd%np(5))

    LA_mhd%dom_np(1:3) = H_mesh%dom_np
    LA_mhd%np(1:3) = H_mesh%np
    LA_mhd%dom_np(4) = pmag_mesh%dom_np
    LA_mhd%np(4) = pmag_mesh%np
    LA_mhd%dom_np(5) = phi_mesh%dom_np
    LA_mhd%np(5) = phi_mesh%np

    IF (ASSOCIATED(LA_mhd%loc_to_glob)) DEALLOCATE(LA_mhd%loc_to_glob)
    ALLOCATE(LA_mhd%loc_to_glob(1,np_glob))
    np_H = 3*H_mesh%dom_np
    np_pmag = pmag_mesh%dom_np
    np_phi = phi_mesh%dom_np
    np_m = H_mesh%dom_np
    DO i = 1, np_m
       DO j = 1, 3
          LA_mhd%loc_to_glob(1,(j-1)*np_m+i) = LA_H%loc_to_glob(j,i)
       END DO
    END DO
    DO i = 1, np_pmag
       LA_mhd%loc_to_glob(1,np_H+i) = LA_pmag%loc_to_glob(1,i)
    END DO
    DO i = 1, np_phi
       LA_mhd%loc_to_glob(1,np_H+np_pmag+i) = LA_phi%loc_to_glob(1,i)
    END DO

    DEALLOCATE (ja_work, nja_glob)
    DEALLOCATE (per_loc)

  END SUBROUTINE st_scr_maxwell_mu_H_p_phi

!!$  SUBROUTINE search_index_copy(mesh,jglob,nb_procs,jloc,proc,out)
!!$    USE def_type_mesh
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),           INTENT(IN) :: mesh
!!$    INTEGER ,                  INTENT(IN) :: jglob, nb_procs
!!$    INTEGER ,                  INTENT(OUT):: jloc, proc
!!$    LOGICAL,                   INTENT(OUT):: out
!!$    INTEGER                               :: p
!!$
!!$    IF (jglob< mesh%loc_to_glob(1) .OR. jglob> mesh%loc_to_glob(2)) THEN
!!$       DO p = 2, nb_procs+1
!!$          IF (mesh%disp(p) > jglob) THEN
!!$             proc = p - 1
!!$             EXIT
!!$          END IF
!!$       END DO
!!$       out = .TRUE.
!!$       jloc = jglob - mesh%disp(proc) + 1
!!$    ELSE
!!$       out = .FALSE.
!!$       jloc = jglob - mesh%loc_to_glob(1) + 1
!!$    END IF
!!$  END SUBROUTINE search_index_copy

  SUBROUTINE search_index(mesh,jglob,nb_procs,jloc,proc,out)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),           INTENT(IN) :: mesh
    INTEGER ,                  INTENT(IN) :: jglob, nb_procs
    INTEGER ,                  INTENT(OUT):: jloc, proc
    LOGICAL,                   INTENT(OUT):: out
    INTEGER                               :: p

    IF (mesh%dom_np == 0) THEN
       DO p = 2, nb_procs+1
          IF (mesh%disp(p) > jglob) THEN
             proc = p - 1
             EXIT
          END IF
       END DO
       out = .TRUE.
       jloc = jglob - mesh%disp(proc) + 1
    ELSE IF (jglob< mesh%loc_to_glob(1) .OR. jglob> mesh%loc_to_glob(mesh%dom_np)) THEN
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
  END SUBROUTINE search_index

  SUBROUTINE search_index_block(mesh, jglob, rank, kmax, nb_procs, jloc, proc, out)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type),           INTENT(IN)   :: mesh
    INTEGER,                   INTENT(IN)   :: jglob,rank, kmax, nb_procs
    INTEGER,                   INTENT(OUT)  :: jloc, proc
    LOGICAL,                   INTENT(OUT)  :: out
    INTEGER                                 :: p


    out = .TRUE.
    DO p = 2, nb_procs
       IF (kmax*(mesh%disp(p)-1) .GE.  jglob) THEN
          proc = p-1
          out = .FALSE.
          EXIT
       END IF
    END DO
    IF (out) proc = nb_procs

    jloc = jglob - kmax*(mesh%disp(proc)-1)
    out = (proc /= rank+1)

  END SUBROUTINE SEARCH_INDEX_BLOCK


END MODULE st_csr_mhd
