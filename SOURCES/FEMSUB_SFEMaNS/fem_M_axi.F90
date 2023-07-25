!
!
!Authors: Jean-Luc Guermond, Raphael Laguerre, Luigi Quartapelle, Copyrights 1996, 2000, 2004
!Revised, Jean-Luc Guermond, Franky Luddens, January, 2011
!
MODULE fem_M_axi
  USE def_type_mesh
#include "petsc/finclude/petsc.h"
  USE petsc
CONTAINS


  SUBROUTINE qs_diff_mass_vect_M (type_op, LA, mesh, visco, mass, stab, mode, matrix)
    !=================================================
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    REAL(KIND=8),                 INTENT(IN)               :: visco, mass, stab
    INTEGER,                      INTENT(IN)               :: type_op, mode
    TYPE(petsc_csr_LA)                                     :: LA
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(:),   ALLOCATABLE              :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc, b_loc
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE              :: mat_loc
    REAL(KIND=8)                                           :: ray, eps1, eps2
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, k_max, n_w, ix, jx, ki, kj
    REAL(KIND=8) :: viscolm, xij
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)


    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 ! 2x2 Structure
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !2x2 Structure
    ELSEIF (type_op == 3) THEN
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 ! Scalar Structure
    ELSE
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1
       CALL error_Petsc('BUG in qs_diff_mass_vect_M, type_op not correct')
    ENDIF

    IF (k_max/=SIZE(LA%loc_to_glob,1)) THEN
       CALL error_Petsc('BUG in qs_diff_mass_vect_petsc_M, k_max/=SIZE(LA%loc_to_glob,1)')
    END IF

    n_w = mesh%gauss%n_w
    DO l = 1, mesh%gauss%l_G
       DO ni = 1, n_w
          DO nj = 1, n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    ALLOCATE(mat_loc(k_max*n_w,k_max*n_w),idxm(k_max*n_w),idxn(k_max*n_w))
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)

       a_loc = 0.d0
       b_loc = 0.d0
       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          viscolm  = (visco + stab*mesh%hloc(m))*mesh%gauss%rj(l,m)
          DO nj = 1, n_w
             DO ni = 1, n_w

                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))

                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij    &
                     + viscolm*eps1*wwprod(ni,nj,l)/ray &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
                b_loc(ni,nj) = b_loc(ni,nj) + eps2*viscolm*2*mode*wwprod(ni,nj,l)/ray
             END DO
          END DO
       END DO

       mat_loc = 0.d0
       DO ki= 1, k_max
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxm(ix) = iglob-1
             DO kj = 1, k_max
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   idxn(jx) = jglob-1
                   IF (ki==kj) THEN
                      mat_loc(ix,jx) = mat_loc(ix,jx) + a_loc(ni,nj)
                   ELSE
                      mat_loc(ix,jx) = mat_loc(ix,jx) + b_loc(ni,nj)
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, k_max*n_w, idxm, k_max*n_w, idxn, mat_loc, ADD_VALUES, ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

    DEALLOCATE(mat_loc,idxm,idxn)

  END SUBROUTINE qs_diff_mass_vect_M

  SUBROUTINE qs_diff_mass_scal_M (mesh, LA, visco, mass, stab, mode, matrix)
    !=================================================
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    TYPE(petsc_csr_LA)                                     :: LA
    REAL(KIND=8),                 INTENT(IN)               :: visco, mass, stab
    INTEGER,                      INTENT(IN)               :: mode
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc
    REAL(KIND=8)                                           :: ray
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, n_w
    REAL(KIND=8) :: viscolm, xij
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    n_w = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       a_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          viscolm  = (visco + stab*mesh%hloc(m))*mesh%gauss%rj(l,m)

          DO nj = 1, n_w
             DO ni = 1, n_w
                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij  &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix,n_w,idxm,n_w,idxn,a_loc,ADD_VALUES,ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_diff_mass_scal_M


  SUBROUTINE qs_diff_mass_scal_M_variant(mesh, LA, heat_capa, visco, mass, temp_list_robin_sides, &
       convection_coeff, stab, mode, matrix)
    USE my_util
    !=================================================
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    TYPE(petsc_csr_LA)                                     :: LA
    REAL(KIND=8),                 INTENT(IN)               :: mass, stab
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)               :: heat_capa, visco, convection_coeff ! MODIFICATION: Robin coeff = convection coeff h for temperature case
    INTEGER, DIMENSION(:),   INTENT(IN)                    :: temp_list_robin_sides ! MODIFICATION: robin sides added to build the matrix
    INTEGER,                      INTENT(IN)               :: mode
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc
    REAL(KIND=8)                                           :: ray
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, n_w, n_ws, ms, ls, ib ! MODIFICATION: for Robin
    REAL(KIND=8) :: viscolm, xij, x ! MODIFICATION: for Robin
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,mesh%gauss%n_ws) :: h_phii_phij ! MODIFICATION: terms for Robin
    INTEGER, DIMENSION(1:1)                                  :: coeff_index ! MODIFICATION: index of robin coefficient in the list

    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    n_w = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       a_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          viscolm  = (visco(m) + stab*mesh%hloc(m))*mesh%gauss%rj(l,m)
          DO nj = 1, n_w
             DO ni = 1, n_w
                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij  &
                     + heat_capa(m) * mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix,n_w,idxm,n_w,idxn,a_loc,ADD_VALUES,ierr)
    ENDDO

    !===Robin conditions ! MODIFICATION: Robin: addition of the term int_(partial Omega) h*u*v, with h the convection coefficient

    IF (SIZE(temp_list_robin_sides) > 0) THEN

       n_ws = mesh%gauss%n_ws
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(temp_list_robin_sides - mesh%sides(ms))) > 0) CYCLE
          h_phii_phij = 0d0
          coeff_index = MINLOC(ABS(temp_list_robin_sides - mesh%sides(ms)))
          DO ls = 1, mesh%gauss%l_Gs
             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))* mesh%gauss%wws(:,ls))
             x = convection_coeff(coeff_index(1)) * ray * mesh%gauss%rjs(ls,ms)
             DO ni=1, n_ws
                DO nj=1, n_ws
                   h_phii_phij(ni,nj) = h_phii_phij(ni,nj) + &
                        x * mesh%gauss%wws(ni,ls) * mesh%gauss%wws(nj,ls)
                ENDDO
             ENDDO
          ENDDO
          DO ni = 1, n_ws
             i =  mesh%jjs(ni,ms)
             ib = LA%loc_to_glob(1,i)
             idxn(ni) = ib - 1
          END DO
          CALL MatSetValues(matrix, n_ws, idxn(1:n_ws), n_ws, idxn(1:n_ws), &
               h_phii_phij, ADD_VALUES, ierr)
       END DO

    END IF

    !===End Robin conditions

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_diff_mass_scal_M_variant

  SUBROUTINE qs_mass_vect_3x3(LA, mesh, mass, matrix)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 3np * 3np, pour trois composantes
    ! (V1,V4,V5), (V2,V3,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1, 4 et 5
    !                       2 pour les composantes 2, 3 et 6
    !------------------------------------------------------------------------
    USE my_util
    USE input_data
    IMPLICIT NONE

    REAL(KIND=8),                 INTENT(IN)    :: mass
    TYPE(petsc_csr_la)                          :: LA
    TYPE(mesh_type)                             :: mesh
    INTEGER ::l, m, ni, nj, i, j, ki, kj, n_w
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij
    REAL(KIND=8) :: ray
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w,3*mesh%gauss%n_w)   :: mat_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                    :: idxn, jdxn
    INTEGER                                                      :: ix, jx, iglob, jglob
    INTEGER,      DIMENSION(mesh%gauss%n_w)                      :: jj_loc
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    n_w  = mesh%gauss%n_w


    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, mesh%me
       jj_loc = mesh%jj(:,m)

       aij = 0.d0
       DO l = 1, mesh%gauss%l_G
          ray = 0.d0
          DO ni = 1, mesh%gauss%n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          DO nj = 1, mesh%gauss%n_w; j = jj_loc(nj)
             DO ni = 1, mesh%gauss%n_w;  i = jj_loc(ni)
                aij(ni,nj) =  aij(ni,nj) + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m)
             ENDDO
          ENDDO
       ENDDO

       mat_loc = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxn(ix) = iglob-1
             DO kj = 1, 3
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   jdxn(jx) = jglob-1
                   IF (ki==kj) THEN
                      mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)
    END DO

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_mass_vect_3x3


  SUBROUTINE qs_diff_mass_vect_3x3(type_op, LA, mesh, visco, mass, stab, stab_bdy_ns, i_mode, mode, matrix)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 3np * 3np, pour trois composantes
    ! (V1,V4,V5), (V2,V3,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1, 4 et 5
    !                       2 pour les composantes 2, 3 et 6
    !------------------------------------------------------------------------
    USE my_util
    USE input_data
    IMPLICIT NONE
    INTEGER     ,                 INTENT(IN)    :: type_op, mode, i_mode
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass, stab, stab_bdy_ns
    TYPE(petsc_csr_la)                          :: LA
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER :: k, l, m, ni, nj, i, j, np, ki, kj, k_max, ls, ms, n_w, n_ws
    REAL(KIND=8) :: xij, viscolm
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij, cij, dij, eij, fij
    REAL(KIND=8) :: ray, eps1, eps2, z, hm1, y, x, stab_normal
    REAL(KIND=8) :: two = 2.d0, coeff_stab_normal=10.d0
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w,3*mesh%gauss%n_w)   :: mat_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                    :: idxn, jdxn
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_ws,2*mesh%gauss%n_ws) :: mat_loc_s
    INTEGER,      DIMENSION(2*mesh%gauss%n_ws)                   :: idxn_s, jdxn_s
    INTEGER                                                      :: ix, jx, iglob, jglob
    INTEGER,      DIMENSION(mesh%gauss%n_w)                      :: jj_loc
    REAL(KIND=8) ::  hm, viscomode, hh
    !=====DEUX INSERTIONS VALEURS A FAIRE....
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    np   = SIZE(mesh%rr,2)
    n_w  = mesh%gauss%n_w
    n_ws = mesh%gauss%n_ws

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 3) THEN
       !cas du laplacien scalaire
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 !Structure scalaire
    ELSE
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1
       CALL error_petsc('probleme de type d''operateur')
    ENDIF


    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, mesh%me
       jj_loc = mesh%jj(:,m)

       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          hh=mesh%hloc(m)
          hm=MIN(mesh%hm(i_mode),hh)!WRONG choice
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)
          viscolm  = (visco + stab*hh)*mesh%gauss%rj(l,m) ! We add the constant stabilization here
          viscomode = (visco + stab*hm)*mesh%gauss%rj(l,m)

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)

             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.d0
                DO k = 1, mesh%gauss%k_d
                   xij =  xij + mesh%gauss%dw(k,nj,l,m) * mesh%gauss%dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                z = ray * viscolm* xij    &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                cij(ni,nj) =  cij(ni,nj) + z
                aij(ni,nj) =  aij(ni,nj) + z + viscomode*eps1*wwprod(ni,nj,l)/ray
                !blocs couplant
                bij(ni,nj) = bij(ni,nj) + eps2*viscomode*2*mode*wwprod(ni,nj,l)/ray
             ENDDO
          ENDDO

       ENDDO

       mat_loc = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxn(ix) = iglob-1
             DO kj = 1, 3
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   jdxn(jx) = jglob-1
                   IF ((ki .LT. 3) .AND. (kj .LT. 3)) THEN
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
                      ELSE
                         mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
                      END IF
                   ELSE ! ki=3 OR kj=3
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

       !==Calcul de visco (grad u)T . (grad v)
       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       dij = 0.d0
       eij = 0.d0
       fij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm  = visco*mesh%gauss%rj(l,m)*ray

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)
             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)
                aij(ni,nj) = aij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(1,ni,l,m)*mesh%gauss%dw(1,nj,l,m) + wwprod(ni,nj,l)/ray**2)
                bij(ni,nj) = bij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%ww(ni,l)*mesh%gauss%dw(1,nj,l,m)/ray+eps2*mode*wwprod(ni,nj,l)/ray**2)
                cij(ni,nj) = cij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(1,nj,l,m))
                dij(ni,nj) = dij(ni,nj) &
                     + viscolm*(-mesh%gauss%dw(1,ni,l,m)*mesh%gauss%ww(nj,l)/ray+(mode/ray)**2*wwprod(ni,nj,l) &
                     -mesh%gauss%dw(1,nj,l,m)*mesh%gauss%ww(ni,l)/ray)
                eij(ni,nj) = eij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%dw(2,ni,l,m)*mesh%gauss%ww(nj,l)/ray)
                fij(ni,nj) = fij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m))
             END DO
          END DO
       END DO

       mat_loc=0.d0
       idxn=0
       jdxn=0
       DO ni = 1, n_w
          DO ki = 1, 3
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w + ni
             idxn(ix) = iglob-1
          END DO
       END DO
       jdxn=idxn

       DO ni = 1, n_w
          DO nj = 1, n_w
             !=== Line i 1 (Vr)
             ix = ni
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)

             !=== Line i 2 (Vt)
             ix = ni+n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + dij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(ni,nj)

             !=== Line i 3 (Vz)
             ix = ni+2*n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(nj,ni)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + fij(ni,nj)
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)
    ENDDO
    !== Fin du Calcul de visco (grad u)T . (grad v)

    IF (inputs%vv_nb_dirichlet_normal_velocity>0) THEN
       !===Surface terms
       stab_normal = coeff_stab_normal*(1.d0+visco)
       !===Normalization is not done properly: (1.d0+visco)
       !===To be revisited some day.
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)- inputs%vv_list_dirichlet_normal_velocity_sides)).NE.0) CYCLE
          aij = 0.d0
          bij = 0.d0
          cij = 0.d0
          dij = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             !--------On calcule le rayon du point gauss
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
             IF (ray.LT.1.d-10) CYCLE
             hm1 = stab_bdy_ns/SUM(mesh%gauss%rjs(:,ms))
             x = two*stab_normal*hm1*mesh%gauss%rjs(ls,ms)*ray
             z = two*mesh%gauss%rjs(ls,ms)*ray*visco

             DO ni = 1, mesh%gauss%n_ws
                DO nj = 1, mesh%gauss%n_ws
                   y = x * mesh%gauss%wws(ni,ls)*mesh%gauss%wws(nj,ls)
                   aij(ni,nj) = aij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   bij(ni,nj) = bij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms) &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                   cij(ni,nj)  = cij(ni,nj)  + y *mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%rnorms(1,ls,ms) &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   dij(ni,nj) = dij(ni,nj) + y *mesh%gauss%rnorms(2,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                END DO
             END DO
          END DO

          !++++++++++++++++
          !=== In the following loops, ki=1 for Vr and ki=2 for Vz
          !===   ==> 2ki-1 = 1 for Vr, 2ki-1 = 3 for Vz   (LA%loc_to_glob(2*ki-1,i))
          mat_loc_s = 0.d0
          idxn_s    = 0
          jdxn_s    = 0
          DO ki = 1, 2
             DO ni = 1, n_ws
                i = mesh%jjs(ni,ms)
                iglob = LA%loc_to_glob(2*ki-1,i)
                ix = (ki-1)*n_ws+ni
                idxn_s(ix) = iglob-1
                DO kj = 1, 2
                   DO nj = 1, n_ws
                      j = mesh%jjs(nj,ms)
                      jglob = LA%loc_to_glob(2*kj-1,j)
                      jx = (kj-1)*n_ws+nj
                      jdxn_s(jx) = jglob-1
                      IF ((ki == 1) .AND. (kj == 1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + aij(ni,nj) + aij(nj,ni)
                      ELSE IF ((ki == 1) .AND. (kj==2)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + bij(ni,nj) + cij(nj,ni)
                      ELSE IF ((ki == 2) .AND. (kj==1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + cij(ni,nj) + bij(nj,ni)
                      ELSE
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + dij(ni,nj) + dij(nj,ni)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          CALL MatSetValues(matrix, 2*n_ws, idxn_s, 2*n_ws, jdxn_s, mat_loc_s, ADD_VALUES, ierr)
       END DO
    END IF

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_diff_mass_vect_3x3

  SUBROUTINE qs_diff_mass_vect_3x3_divpenal(type_op, LA, mesh, visco, mass, stab, stab_bdy_ns, i_mode, mode, matrix)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 3np * 3np, pour trois composantes
    ! (V1,V4,V5), (V2,V3,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1, 4 et 5
    !                       2 pour les composantes 2, 3 et 6
    !------------------------------------------------------------------------
    USE my_util
    USE input_data
    IMPLICIT NONE


    INTEGER     ,                 INTENT(IN)    :: type_op, mode, i_mode
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass, stab, stab_bdy_ns
    TYPE(petsc_csr_la)                          :: LA
    TYPE(mesh_type), TARGET                     :: mesh

!!$    INTEGER,      DIMENSION(:,:), POINTER       :: jj
!!$    INTEGER,                      POINTER       :: me
!!$    REAL(KIND=8), DIMENSION(:,:) , POINTER      :: rr

    INTEGER :: k, l, m, ni, nj, i, j, np, ki, kj, k_max, ls, ms, n_w, n_ws
    REAL(KIND=8) :: xij, viscolm, div_penal
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij, cij, dij, eij, fij
    REAL(KIND=8) :: ray, eps1, eps2, z, hm1, y, x, stab_normal
    REAL(KIND=8) :: two = 2.d0, coeff_stab_normal=10.d0
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w,3*mesh%gauss%n_w)   :: mat_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                    :: idxn, jdxn
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_ws,2*mesh%gauss%n_ws) :: mat_loc_s
    INTEGER,      DIMENSION(2*mesh%gauss%n_ws)                   :: idxn_s, jdxn_s
    INTEGER                                                      :: ix, jx, iglob, jglob
    INTEGER,      DIMENSION(mesh%gauss%n_w)                      :: jj_loc
    REAL(KIND=8) :: viscomode, hm
    !=====DEUX INSERTIONS VALEURS A FAIRE....
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    np   = SIZE(mesh%rr,2)
    n_w  = mesh%gauss%n_w
    n_ws = mesh%gauss%n_ws

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 3) THEN
       !cas du laplacien scalaire
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 !Structure scalaire
    ELSE
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1
       CALL error_petsc('probleme de type d''operateur')
    ENDIF


    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, mesh%me
       jj_loc = mesh%jj(:,m)

       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm  = (visco + stab*mesh%hloc(m))*mesh%gauss%rj(l,m) ! We add the constant stabilization here
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)
          hm=MIN(mesh%hm(i_mode),mesh%hloc(m))!WRONG choice
          viscomode = (visco + stab*hm)*mesh%gauss%rj(l,m)

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)

             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.d0
                DO k = 1, mesh%gauss%k_d
                   xij =  xij + mesh%gauss%dw(k,nj,l,m) * mesh%gauss%dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                z = ray * viscolm* xij    &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                cij(ni,nj) =  cij(ni,nj) + z
                aij(ni,nj) =  aij(ni,nj) + z + viscomode*eps1*wwprod(ni,nj,l)/ray
                !blocs couplant
                bij(ni,nj) = bij(ni,nj) + eps2*viscomode*2*mode*wwprod(ni,nj,l)/ray
             ENDDO
          ENDDO

       ENDDO



       !++++++++++++++++
!!$       DO ki= 1, k_max
!!$          DO ni = 1, mesh%gauss%n_w
!!$             i = mesh%jj(ni, m)
!!$             ib = i + (ki-1)*np
!!$             DO kj = 1, k_max
!!$                DO nj = 1, mesh%gauss%n_w
!!$                   j = mesh%jj(nj, m)
!!$                   jb = j + (kj-1)*np
!!$                   DO p = ia(ib),  ia(ib+1) - 1
!!$                      IF (ja(p) == jb) THEN
!!$                         IF (ki==kj) THEN
!!$                            a0(p) = a0(p) + aij(ni,nj)
!!$                         ELSE
!!$                            a0(p) = a0(p) + bij(ni,nj)
!!$                         END IF
!!$                         EXIT
!!$                      ENDIF
!!$                   END DO
!!$                END DO
!!$             END DO
!!$          END DO
!!$       END DO
!!$
!!$       IF (type_op /= 3) THEN
!!$          ki= 3
!!$          DO ni = 1, mesh%gauss%n_w
!!$             i = mesh%jj(ni, m)
!!$             ib = i + (ki-1)*np
!!$             kj = 3
!!$             DO nj = 1, mesh%gauss%n_w
!!$                j = mesh%jj(nj, m)
!!$                jb = j + (kj-1)*np
!!$                DO p = ia(ib),  ia(ib+1) - 1
!!$                   IF (ja(p) == jb) THEN
!!$                      a0(p) = a0(p) + cij(ni,nj)
!!$                      EXIT
!!$                   ENDIF
!!$                END DO
!!$             END DO
!!$          END DO
!!$       END IF

       mat_loc = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxn(ix) = iglob-1
             DO kj = 1, 3
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   jdxn(jx) = jglob-1
                   IF ((ki .LT. 3) .AND. (kj .LT. 3)) THEN
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
                      ELSE
                         mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
                      END IF
                   ELSE ! ki=3 OR kj=3
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

       !+++---------------------
!!$       IF (type_op==3) CYCLE
       !==Calcul de visco (grad u)T . (grad v)
       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       dij = 0.d0
       eij = 0.d0
       fij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm  = visco*mesh%gauss%rj(l,m)*ray
          !LC 2017/01/27
          !div_penal = inputs%div_stab_in_ns*mesh%gauss%rj(l,m)*ray
          div_penal = inputs%div_stab_in_ns/inputs%Re*mesh%gauss%rj(l,m)*ray

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)
             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)
                aij(ni,nj) = aij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(1,ni,l,m)*mesh%gauss%dw(1,nj,l,m) + wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(mesh%gauss%dw(1,nj,l,m) &
                     + mesh%gauss%ww(nj,l)/ray)
                bij(ni,nj) = bij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%ww(ni,l)*mesh%gauss%dw(1,nj,l,m)/ray+eps2*mode*wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(eps2*(mode/ray)*mesh%gauss%ww(nj,l))
                cij(ni,nj) = cij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(1,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*mesh%gauss%dw(2,nj,l,m)
                dij(ni,nj) = dij(ni,nj) &
                     + viscolm*(-mesh%gauss%dw(1,ni,l,m)*mesh%gauss%ww(nj,l)/ray+(mode/ray)**2*wwprod(ni,nj,l) &
                     -mesh%gauss%dw(1,nj,l,m)*mesh%gauss%ww(ni,l)/ray) &
                     + div_penal*wwprod(ni,nj,l)*(mode/ray)**2
                eij(ni,nj) = eij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%dw(2,ni,l,m)*mesh%gauss%ww(nj,l)/ray) &
                     + div_penal*eps2*(mode/ray)*mesh%gauss%ww(ni,l)*mesh%gauss%dw(2,nj,l,m)
                fij(ni,nj) = fij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m))
             END DO
          END DO
       END DO
       !aij = 0.d0
       !bij = 0.d0
       !cij = 0.d0
       !dij = 0.d0
       !eij = 0.d0
       !fij = 0.d0
       !++++++++++++++++
       mat_loc=0.d0
       idxn=0
       jdxn=0
       DO ni = 1, n_w
          DO ki = 1, 3
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w + ni
             idxn(ix) = iglob-1
          END DO
       END DO
       jdxn=idxn

       DO ni = 1, n_w
          DO nj = 1, n_w
             !=== Line i 1 (Vr)
             ix = ni
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)

             !=== Line i 2 (Vt)
             ix = ni+n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + dij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(ni,nj)

             !=== Line i 3 (Vz)
             ix = ni+2*n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(nj,ni)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + fij(ni,nj)
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

!!$       DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)
!!$          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)
!!$             ib = i
!!$             jb = j
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + aij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb= j + np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + bij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb = j + 2*np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + cij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             ib = i + np
!!$             jb = j
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + bij(nj,ni)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb = j + np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + dij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb = j + 2*np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + eij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             ib = i + 2*np
!!$             jb = j
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + cij(nj,ni)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb = j + np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + eij(nj,ni)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$             jb = j + 2*np
!!$             DO p = ia(ib),  ia(ib+1) - 1
!!$                IF (ja(p) == jb) THEN
!!$                   a0(p) = a0(p) + fij(ni,nj)
!!$                   EXIT
!!$                ENDIF
!!$             END DO
!!$          END DO
!!$       END DO
    ENDDO
    !== Fin du Calcul de visco (grad u)T . (grad v)
    !++++++++++++++++------

    IF (inputs%vv_nb_dirichlet_normal_velocity>0) THEN
       !===Surface terms
       stab_normal = coeff_stab_normal*(1.d0+visco)
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)- inputs%vv_list_dirichlet_normal_velocity_sides)).NE.0) CYCLE
          aij = 0.d0
          bij = 0.d0
          cij = 0.d0
          dij = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             !--------On calcule le rayon du point gauss
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
             IF (ray.LT.1.d-10) CYCLE
             hm1 = stab_bdy_ns/SUM(mesh%gauss%rjs(:,ms))
             x = two*stab_normal*hm1*mesh%gauss%rjs(ls,ms)*ray
             z = two*mesh%gauss%rjs(ls,ms)*ray*visco

             DO ni = 1, mesh%gauss%n_ws
                DO nj = 1, mesh%gauss%n_ws
                   y = x * mesh%gauss%wws(ni,ls)*mesh%gauss%wws(nj,ls)
                   aij(ni,nj) = aij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   bij(ni,nj) = bij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms) &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                   cij(ni,nj)  = cij(ni,nj)  + y *mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%rnorms(1,ls,ms) &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   dij(ni,nj) = dij(ni,nj) + y *mesh%gauss%rnorms(2,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                END DO
             END DO
          END DO

          !++++++++++++++++
          !=== In the following loops, ki=1 for Vr and ki=2 for Vz
          !===   ==> 2ki-1 = 1 for Vr, 2ki-1 = 3 for Vz   (LA%loc_to_glob(2*ki-1,i))
          mat_loc_s = 0.d0
          idxn_s    = 0
          jdxn_s    = 0
          DO ki = 1, 2
             DO ni = 1, n_ws
                i = mesh%jjs(ni,ms)
                iglob = LA%loc_to_glob(2*ki-1,i)
                ix = (ki-1)*n_ws+ni
                idxn_s(ix) = iglob-1
                DO kj = 1, 2
                   DO nj = 1, n_ws
                      j = mesh%jjs(nj,ms)
                      jglob = LA%loc_to_glob(2*kj-1,j)
                      jx = (kj-1)*n_ws+nj
                      jdxn_s(jx) = jglob-1
                      IF ((ki == 1) .AND. (kj == 1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + aij(ni,nj) + aij(nj,ni)
                      ELSE IF ((ki == 1) .AND. (kj==2)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + bij(ni,nj) + cij(nj,ni)
                      ELSE IF ((ki == 2) .AND. (kj==1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + cij(ni,nj) + bij(nj,ni)
                      ELSE
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + dij(ni,nj) + dij(nj,ni)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          CALL MatSetValues(matrix, 2*n_ws, idxn_s, 2*n_ws, jdxn_s, mat_loc_s, ADD_VALUES, ierr)

!!$       DO ki= 1, 3, 2
!!$          DO ni = 1, mesh%gauss%n_ws
!!$             i = mesh%jjs(ni, ms)
!!$             ib = i + (ki-1)*np
!!$             DO kj = 1, 3, 2
!!$                DO nj = 1, mesh%gauss%n_ws
!!$                   j = mesh%jjs(nj, ms)
!!$                   jb = j + (kj-1)*np
!!$                   DO p = ia(ib),  ia(ib+1) - 1
!!$                      IF (ja(p) == jb) THEN
!!$                         IF      (ki == 1 .AND. kj == 1) THEN
!!$                            a0(p) = a0(p) + aij(ni,nj) + aij(nj,ni)
!!$                         ELSE IF (ki == 1 .AND. kj == 3) THEN
!!$                            a0(p) = a0(p) + bij(ni,nj) + cij(nj,ni)
!!$                         ELSE IF (ki == 3 .AND. kj == 1) THEN
!!$                            a0(p) = a0(p) + cij(ni,nj) + bij(nj,ni)
!!$                         ELSE
!!$                            a0(p) = a0(p) + dij(ni,nj) + dij(nj,ni)
!!$                         END IF
!!$                         EXIT
!!$                      ENDIF
!!$                   END DO
!!$                END DO
!!$             END DO
!!$          END DO
!!$       END DO
          !++++++++++++++++----

       END DO
    END IF

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)


  END SUBROUTINE qs_diff_mass_vect_3x3_divpenal


  SUBROUTINE qs_diff_mass_vect_3x3_divpenal_art_comp(type_op, LA, mesh, visco, mass, stab, stab_bdy_ns, &
       stab_art_comp, i_mode, mode, matrix)
    !========================================================================
    ! laplacien vectoriel qui renvoie une matrice 3np * 3np, pour trois composantes
    ! (V1,V4,V5), (V2,V3,V6)
    ! calcul de l'operateur mass-visco(lap-eps1*1/r^2) et eps2*2m*visco/r^2
    ! eps1 et eps2 sont des parametres qui definissent le type d'operateur
    ! type de l'operateur : 1 pour les composantes 1, 4 et 5
    !                       2 pour les composantes 2, 3 et 6
    !------------------------------------------------------------------------
    USE my_util
    USE input_data
    IMPLICIT NONE


    INTEGER     ,                 INTENT(IN)    :: type_op, mode, i_mode
    REAL(KIND=8),                 INTENT(IN)    :: visco, mass, stab, stab_bdy_ns, stab_art_comp
    TYPE(petsc_csr_la)                          :: LA
    TYPE(mesh_type), TARGET                     :: mesh
    INTEGER :: k, l, m, ni, nj, i, j, np, ki, kj, k_max, ls, ms, n_w, n_ws
    REAL(KIND=8) :: xij, viscolm, div_penal
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: aij, bij, cij, dij, eij, fij
    REAL(KIND=8) :: ray, eps1, eps2, z, hm1, y, x, stab_normal
    REAL(KIND=8) :: two = 2.d0, coeff_stab_normal=10.d0
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w,3*mesh%gauss%n_w)   :: mat_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                    :: idxn, jdxn
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_ws,2*mesh%gauss%n_ws) :: mat_loc_s
    INTEGER,      DIMENSION(2*mesh%gauss%n_ws)                   :: idxn_s, jdxn_s
    INTEGER                                                      :: ix, jx, iglob, jglob
    INTEGER,      DIMENSION(mesh%gauss%n_w)                      :: jj_loc
    REAL(KIND=8) :: viscomode, hh, hm
    !=====DEUX INSERTIONS VALEURS A FAIRE....
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr

    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    np   = SIZE(mesh%rr,2)
    n_w  = mesh%gauss%n_w
    n_ws = mesh%gauss%n_ws

    IF (type_op == 1) THEN
       eps1 = 1.d0
       eps2 = 1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 2) THEN
       eps1 = 1.d0
       eps2 = -1.d0
       k_max = 2 !Structure vectorielle
    ELSEIF (type_op == 3) THEN
       !cas du laplacien scalaire
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1 !Structure scalaire
    ELSE
       eps1 = 0.d0
       eps2 = 0.d0
       k_max = 1
       CALL error_petsc('probleme de type d''operateur')
    ENDIF


    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    DO m = 1, mesh%me
       jj_loc = mesh%jj(:,m)

       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          hh=mesh%hloc(m)
          hm=MIN(mesh%hm(i_mode),hh)!WRONG choice
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)
          viscolm   = (visco + stab*hh)*mesh%gauss%rj(l,m) ! We add the constant stabilization here
          viscomode = (visco + stab*hm)*mesh%gauss%rj(l,m)

          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)

             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)

                !grad(u).grad(v) en r et z
                xij = 0.d0
                DO k = 1, mesh%gauss%k_d
                   xij =  xij + mesh%gauss%dw(k,nj,l,m) * mesh%gauss%dw(k,ni,l,m)
                END DO

                !blocs diagonaux
                z = ray * viscolm* xij    &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                cij(ni,nj) =  cij(ni,nj) + z
                aij(ni,nj) =  aij(ni,nj) + z + viscomode*eps1*wwprod(ni,nj,l)/ray
                !blocs couplant
                bij(ni,nj) = bij(ni,nj) + eps2*viscomode*2*mode*wwprod(ni,nj,l)/ray
             ENDDO
          ENDDO

       ENDDO

       mat_loc = 0.d0
       DO ki= 1, 3
          DO ni = 1, n_w
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w+ni
             idxn(ix) = iglob-1
             DO kj = 1, 3
                DO nj = 1, n_w
                   j = jj_loc(nj)
                   jglob = LA%loc_to_glob(kj,j)
                   jx = (kj-1)*n_w+nj
                   jdxn(jx) = jglob-1
                   IF ((ki .LT. 3) .AND. (kj .LT. 3)) THEN
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
                      ELSE
                         mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
                      END IF
                   ELSE ! ki=3 OR kj=3
                      IF (ki==kj) THEN
                         mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)
                      END IF
                   END IF
                END DO
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)

       !+++---------------------
       !==Calcul de visco (grad u)T . (grad v)
       aij = 0.d0
       bij = 0.d0
       cij = 0.d0
       dij = 0.d0
       eij = 0.d0
       fij = 0.d0
       DO l = 1, mesh%gauss%l_G

          !--------On calcule le rayon du point gauss
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni,m)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          viscolm   = visco*mesh%gauss%rj(l,m)*ray
          div_penal = stab_art_comp*mesh%gauss%rj(l,m)*ray
          DO nj = 1, mesh%gauss%n_w; j = mesh%jj(nj, m)
             DO ni = 1, mesh%gauss%n_w;  i = mesh%jj(ni, m)
                aij(ni,nj) = aij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(1,ni,l,m)*mesh%gauss%dw(1,nj,l,m) + wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(mesh%gauss%dw(1,nj,l,m) &
                     + mesh%gauss%ww(nj,l)/ray)
                bij(ni,nj) = bij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%ww(ni,l)*mesh%gauss%dw(1,nj,l,m)/ray+eps2*mode*wwprod(ni,nj,l)/ray**2) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*(eps2*(mode/ray)*mesh%gauss%ww(nj,l))
                cij(ni,nj) = cij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(1,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(1,ni,l,m) + mesh%gauss%ww(ni,l)/ray)*mesh%gauss%dw(2,nj,l,m)
                dij(ni,nj) = dij(ni,nj) &
                     + viscolm*(-mesh%gauss%dw(1,ni,l,m)*mesh%gauss%ww(nj,l)/ray+(mode/ray)**2*wwprod(ni,nj,l) &
                     -mesh%gauss%dw(1,nj,l,m)*mesh%gauss%ww(ni,l)/ray) &
                     + div_penal*wwprod(ni,nj,l)*(mode/ray)**2
                eij(ni,nj) = eij(ni,nj) &
                     + viscolm*(-eps2*mode*mesh%gauss%dw(2,ni,l,m)*mesh%gauss%ww(nj,l)/ray) &
                     + div_penal*eps2*(mode/ray)*mesh%gauss%ww(ni,l)*mesh%gauss%dw(2,nj,l,m)
                fij(ni,nj) = fij(ni,nj) &
                     + viscolm*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m)) &
                     + div_penal*(mesh%gauss%dw(2,ni,l,m)*mesh%gauss%dw(2,nj,l,m))
             END DO
          END DO
       END DO

       mat_loc=0.d0
       idxn=0
       jdxn=0
       DO ni = 1, n_w
          DO ki = 1, 3
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(ki,i)
             ix = (ki-1)*n_w + ni
             idxn(ix) = iglob-1
          END DO
       END DO
       jdxn=idxn

       DO ni = 1, n_w
          DO nj = 1, n_w
             !=== Line i 1 (Vr)
             ix = ni
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + aij(ni,nj)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(ni,nj)

             !=== Line i 2 (Vt)
             ix = ni+n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + bij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + dij(ni,nj)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(ni,nj)

             !=== Line i 3 (Vz)
             ix = ni+2*n_w
             !=== Column j 1 (Vr)
             jx = nj
             mat_loc(ix,jx) = mat_loc(ix,jx) + cij(nj,ni)
             !=== Column j 2 (Vt)
             jx = nj+n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + eij(nj,ni)
             !=== Column j 3 (Vz)
             jx = nj+2*n_w
             mat_loc(ix,jx) = mat_loc(ix,jx) + fij(ni,nj)
          END DO
       END DO
       CALL MatSetValues(matrix, 3*n_w, idxn, 3*n_w, jdxn, mat_loc, ADD_VALUES, ierr)
    ENDDO
    !== Fin du Calcul de visco (grad u)T . (grad v)
    !++++++++++++++++------

    IF (inputs%vv_nb_dirichlet_normal_velocity>0) THEN
       !===Surface terms
       stab_normal = coeff_stab_normal*(1.d0+visco)
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(mesh%sides(ms)- inputs%vv_list_dirichlet_normal_velocity_sides)).NE.0) CYCLE
          aij = 0.d0
          bij = 0.d0
          cij = 0.d0
          dij = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             !--------On calcule le rayon du point gauss
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))*mesh%gauss%wws(:,ls))
             IF (ray.LT.1.d-10) CYCLE
             hm1 = stab_bdy_ns/SUM(mesh%gauss%rjs(:,ms))
             x = two*stab_normal*hm1*mesh%gauss%rjs(ls,ms)*ray
             z = two*mesh%gauss%rjs(ls,ms)*ray*visco

             DO ni = 1, mesh%gauss%n_ws
                DO nj = 1, mesh%gauss%n_ws
                   y = x * mesh%gauss%wws(ni,ls)*mesh%gauss%wws(nj,ls)
                   aij(ni,nj) = aij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   bij(ni,nj) = bij(ni,nj) + y *mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms) &
                        - z*mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                   cij(ni,nj)  = cij(ni,nj)  + y *mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%rnorms(1,ls,ms) &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)**2*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(2,nj,ls,ms))
                   dij(ni,nj) = dij(ni,nj) + y *mesh%gauss%rnorms(2,ls,ms)**2 &
                        - z*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%wws(ni,ls) &
                        *(mesh%gauss%rnorms(1,ls,ms)*mesh%gauss%rnorms(2,ls,ms)*mesh%gauss%dw_s(1,nj,ls,ms) &
                        + mesh%gauss%rnorms(2,ls,ms)**2*mesh%gauss%dw_s(2,nj,ls,ms))
                END DO
             END DO
          END DO

          !++++++++++++++++
          !=== In the following loops, ki=1 for Vr and ki=2 for Vz
          !===   ==> 2ki-1 = 1 for Vr, 2ki-1 = 3 for Vz   (LA%loc_to_glob(2*ki-1,i))
          mat_loc_s = 0.d0
          idxn_s    = 0
          jdxn_s    = 0
          DO ki = 1, 2
             DO ni = 1, n_ws
                i = mesh%jjs(ni,ms)
                iglob = LA%loc_to_glob(2*ki-1,i)
                ix = (ki-1)*n_ws+ni
                idxn_s(ix) = iglob-1
                DO kj = 1, 2
                   DO nj = 1, n_ws
                      j = mesh%jjs(nj,ms)
                      jglob = LA%loc_to_glob(2*kj-1,j)
                      jx = (kj-1)*n_ws+nj
                      jdxn_s(jx) = jglob-1
                      IF ((ki == 1) .AND. (kj == 1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + aij(ni,nj) + aij(nj,ni)
                      ELSE IF ((ki == 1) .AND. (kj==2)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + bij(ni,nj) + cij(nj,ni)
                      ELSE IF ((ki == 2) .AND. (kj==1)) THEN
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + cij(ni,nj) + bij(nj,ni)
                      ELSE
                         mat_loc_s(ix,jx) = mat_loc_s(ix,jx) + dij(ni,nj) + dij(nj,ni)
                      END IF
                   END DO
                END DO
             END DO
          END DO
          CALL MatSetValues(matrix, 2*n_ws, idxn_s, 2*n_ws, jdxn_s, mat_loc_s, ADD_VALUES, ierr)
       END DO
    END IF

    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)


  END SUBROUTINE qs_diff_mass_vect_3x3_divpenal_art_comp

  SUBROUTINE qs_00_M (mesh, alpha, LA, matrix)
    !=================================================
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    REAL(KIND=8),                 INTENT(IN)               :: alpha
    TYPE(petsc_csr_LA)                                     :: LA
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: mat_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    INTEGER :: m, l, ni, nj, i, j, iglob, jglob
    REAL(KIND=8) :: al, x, ray
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)

    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       mat_loc = 0.d0
       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO

          al = alpha *mesh%gauss%rj(l,m)*ray
          DO nj = 1, mesh%gauss%n_w;
             j = jj_loc(nj)
             jglob = LA%loc_to_glob(1,j)
             idxn(nj) = jglob-1
             DO ni = 1, mesh%gauss%n_w;
                i = jj_loc(ni)
                iglob = LA%loc_to_glob(1,i)
                idxm(ni) = iglob-1
                x = al*mesh%gauss%ww(nj,l)*mesh%gauss%ww(ni,l)
                mat_loc(ni,nj) = mat_loc(ni,nj) + x
             ENDDO
          ENDDO
       ENDDO
       CALL MatSetValues(matrix, mesh%gauss%n_w, idxm, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)
  END SUBROUTINE qs_00_M

  SUBROUTINE qs_11_M (mesh, visco, LA, matrix)
    !=================================================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    REAL(KIND=8),                 INTENT(IN)    :: visco
    TYPE(petsc_csr_LA)                          :: LA
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: mat_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: idxm, idxn
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: jj_loc
    INTEGER :: m, l, ni, nj, i, j, iglob, jglob
    REAL(KIND=8) :: al, x, ray
    !#include "petsc/finclude/petsc.h"
    Mat            :: matrix
    PetscErrorCode :: ierr
    CALL MatZeroEntries (matrix, ierr)

    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       mat_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1, mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, mesh%gauss%n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          al = visco * mesh%gauss%rj(l,m)
          DO nj = 1, mesh%gauss%n_w;
             DO ni = 1, mesh%gauss%n_w;
                x =al*SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                mat_loc(ni,nj) = mat_loc(ni,nj) + x
             ENDDO
          ENDDO
       ENDDO
       CALL MatSetValues(matrix, mesh%gauss%n_w, idxm, mesh%gauss%n_w, idxn, mat_loc, ADD_VALUES, ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)
  END SUBROUTINE qs_11_M

  SUBROUTINE qs_diff_mass_scal_M_level (mesh, LA, visco, mass, stab, i_mode, mode, matrix)
    !=================================================
    USE input_data
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    TYPE(petsc_csr_LA)                                     :: LA
    REAL(KIND=8),                 INTENT(IN)               :: visco, mass, stab
    INTEGER,                      INTENT(IN)               :: mode, i_mode
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc
    REAL(KIND=8)                                           :: ray
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, n_w
    REAL(KIND=8) :: viscolm, xij, hm, viscomode, hh
    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    n_w = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       a_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          hh=mesh%hloc(m)
          hm=MIN(mesh%hm(i_mode),hh)!WRONG choice
          !hm=0.5d0/inputs%m_max
          !hm=mesh%hm(i_mode) !(JLG April 7 2017)
          viscolm  = (visco + stab*hh)*mesh%gauss%rj(l,m)
          viscomode = (visco + stab*hm)*mesh%gauss%rj(l,m)
          DO nj = 1, n_w
             DO ni = 1, n_w
                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij  &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscomode*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix,n_w,idxm,n_w,idxn,a_loc,ADD_VALUES,ierr)
    ENDDO
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_diff_mass_scal_M_level

  SUBROUTINE qs_diff_mass_scal_M_conc(mesh, LA, visco, mass, conc_list_robin_sides, &
       convection_coeff_conc_lhs, stab, mode, matrix)
    USE my_util
    !=================================================
    IMPLICIT NONE
    TYPE(mesh_type),              INTENT(IN)               :: mesh
    TYPE(petsc_csr_LA)                                     :: LA
    REAL(KIND=8),                 INTENT(IN)               :: mass, stab
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)               :: visco, convection_coeff_conc_lhs  ! Robin coeff = convection coeff h for temperature case
    INTEGER, DIMENSION(:),   INTENT(IN)                    :: conc_list_robin_sides ! Robin sides added to build the matrix
    INTEGER,                      INTENT(IN)               :: mode
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: jj_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm, idxn
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w,mesh%gauss%l_G) :: wwprod
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%n_w) :: a_loc
    REAL(KIND=8)                                           :: ray
    INTEGER      :: m, l, ni, nj, i, j, iglob, jglob, n_w, n_ws, ms, ls, ib ! for Robin
    REAL(KIND=8) :: viscolm, xij, x ! for Robin
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,mesh%gauss%n_ws) :: h_phii_phij ! terms for Robin
    INTEGER, DIMENSION(1:1)                                  :: coeff_index ! index of robin coefficient in the list

    !#include "petsc/finclude/petsc.h"
    Mat                                         :: matrix
    PetscErrorCode                              :: ierr
    CALL MatZeroEntries (matrix, ierr)
    CALL MatSetOption (matrix, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)

    DO l = 1, mesh%gauss%l_G
       DO ni = 1, mesh%gauss%n_w
          DO nj = 1, mesh%gauss%n_w
             wwprod(ni,nj,l) = mesh%gauss%ww(ni,l)*mesh%gauss%ww(nj,l)
          END DO
       END DO
    END DO

    n_w = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       a_loc = 0.d0

       DO nj = 1, mesh%gauss%n_w;
          j = jj_loc(nj)
          jglob = LA%loc_to_glob(1,j)
          idxn(nj) = jglob-1
          DO ni = 1, mesh%gauss%n_w;
             i = jj_loc(ni)
             iglob = LA%loc_to_glob(1,i)
             idxm(ni) = iglob-1
          END DO
       END DO

       DO l = 1,  mesh%gauss%l_G
          !Compute radius of Gauss point
          ray = 0
          DO ni = 1, n_w;  i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          viscolm  = (visco(m) + stab*mesh%hloc(m))*mesh%gauss%rj(l,m)
          DO nj = 1, n_w
             DO ni = 1, n_w
                !grad(u).grad(v) w.r.t. r and z
                xij = SUM(mesh%gauss%dw(:,nj,l,m)*mesh%gauss%dw(:,ni,l,m))
                !start diagonal block
                a_loc(ni,nj) =  a_loc(ni,nj) + ray * viscolm* xij  &
                     + mass*ray*wwprod(ni,nj,l)*mesh%gauss%rj(l,m) &
                     + viscolm*mode**2*wwprod(ni,nj,l)/ray
                !end diagonal block
             END DO
          END DO
       END DO
       CALL MatSetValues(matrix,n_w,idxm,n_w,idxn,a_loc,ADD_VALUES,ierr)
    ENDDO

    !===Robin conditions ! MODIFICATION: Robin: addition of the term int_(partial Omega) h*u*v
    IF (SIZE(conc_list_robin_sides) > 0) THEN

       n_ws = mesh%gauss%n_ws
       DO ms = 1, mesh%mes
          IF (MINVAL(ABS(conc_list_robin_sides - mesh%sides(ms))) > 0) CYCLE
          h_phii_phij = 0d0
          coeff_index = MINLOC(ABS(conc_list_robin_sides - mesh%sides(ms)))
          DO ls = 1, mesh%gauss%l_Gs
             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,mesh%jjs(:,ms))* mesh%gauss%wws(:,ls))
             x = convection_coeff_conc_lhs(coeff_index(1)) * ray * mesh%gauss%rjs(ls,ms)
             DO ni=1, n_ws
                DO nj=1, n_ws
                   h_phii_phij(ni,nj) = h_phii_phij(ni,nj) + &
                        x * mesh%gauss%wws(ni,ls) * mesh%gauss%wws(nj,ls)
                ENDDO
             ENDDO
          ENDDO
          DO ni = 1, n_ws
             i =  mesh%jjs(ni,ms)
             ib = LA%loc_to_glob(1,i)
             idxn(ni) = ib - 1
          END DO
          CALL MatSetValues(matrix, n_ws, idxn(1:n_ws), n_ws, idxn(1:n_ws), &
               h_phii_phij, ADD_VALUES, ierr)
       END DO
    END IF
    !===End Robin conditions
    CALL MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY,ierr)
    CALL MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY,ierr)

  END SUBROUTINE qs_diff_mass_scal_M_conc

END MODULE fem_M_axi
