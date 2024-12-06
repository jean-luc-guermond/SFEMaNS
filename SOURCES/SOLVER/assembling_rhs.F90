MODULE rhs_para_assembling

CONTAINS

  SUBROUTINE rhs_3x3(mesh, vv_3_LA, mode, rhs_in, vb_145, vb_236, opt_tensor, opt_tensor_scaln_bdy)
    USE def_type_mesh
    USE my_util
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                                                     :: mesh
    TYPE(petsc_csr_LA)                                                  :: vv_3_LA
    INTEGER,                                                 INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(mesh%dom_me*mesh%gauss%l_G,6),   INTENT(IN) :: rhs_in
    REAL(KIND=8), DIMENSION(3,mesh%dom_me*mesh%gauss%l_G,6), INTENT(IN), OPTIONAL :: opt_tensor
    REAL(KIND=8), DIMENSION(:,:),                            INTENT(IN), OPTIONAL :: opt_tensor_scaln_bdy
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w)                           :: v145_loc, v236_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                           :: idxm
    INTEGER,      DIMENSION(3*mesh%gauss%n_ws)                          :: idxms
    INTEGER,      DIMENSION(mesh%gauss%n_w)                             :: j_loc
    INTEGER,      DIMENSION(mesh%gauss%n_ws)                            :: js_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)              :: dw_loc
    REAL(KIND=8), DIMENSION(3,6)                                        :: tensor_loc
    REAL(KIND=8), DIMENSION(6)                                          :: smb, f_tensor
    INTEGER      :: m, l, i, ni, index, nw, ix, ki, iglob, type, k, ms, ls, indexs, nws
    REAL(KIND=8) :: ray
    Vec            :: vb_145, vb_236
    PetscErrorCode :: ierr

    CALL VecSet(vb_145, 0.d0, ierr)
    CALL VecSet(vb_236, 0.d0, ierr)

    nw = mesh%gauss%n_w
    index = 0
    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)

       DO ki = 1, 3
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm(ix) = iglob-1
          END DO
       END DO

       v145_loc = 0.d0
       v236_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index + 1
          dw_loc = mesh%gauss%dw(:,:,l,m)

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          IF (PRESENT(opt_tensor)) THEN
             DO TYPE = 1, 6
                DO k = 1, 3
                   tensor_loc(k,TYPE) = opt_tensor(k,index,TYPE)
                END DO
             END DO

             f_tensor(1) = (-mode*tensor_loc(1,4) + tensor_loc(2,3))/ray
             f_tensor(2) = (mode*tensor_loc(1,3) + tensor_loc(2,4))/ray
             f_tensor(3) = (-tensor_loc(1,3) - mode*tensor_loc(2,4))/ray
             f_tensor(4) = (-tensor_loc(1,4) + mode*tensor_loc(2,3))/ray
             f_tensor(5) = -tensor_loc(3,4)*mode/ray
             f_tensor(6) = tensor_loc(3,3)*mode/ray
          ELSE
             f_tensor = 0.d0
          END IF

          smb = (rhs_in(index,:)+f_tensor)*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(1)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(2)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(4)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(3)
          END DO

          ki = 3
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(5)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(6)
          END DO

          IF (PRESENT(opt_tensor)) THEN
             ki = 1
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(1,1) + dw_loc(2,ni)*tensor_loc(1,5))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(1,2) + dw_loc(2,ni)*tensor_loc(1,6))*ray*mesh%gauss%rj(l,m)
             END DO
             ki = 2
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(2,2) + dw_loc(2,ni)*tensor_loc(2,6))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(2,1) + dw_loc(2,ni)*tensor_loc(2,5))*ray*mesh%gauss%rj(l,m)
             END DO
             ki = 3
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(3,1) + dw_loc(2,ni)*tensor_loc(3,5))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(3,2) + dw_loc(2,ni)*tensor_loc(3,6))*ray*mesh%gauss%rj(l,m)
             END DO
          END IF
       ENDDO !l_G

       CALL VecSetValues(vb_145, 3*nw, idxm, v145_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_236, 3*nw, idxm, v236_loc, ADD_VALUES, ierr)

    ENDDO !dom_me

    IF (PRESENT(opt_tensor_scaln_bdy)) THEN
       nws = mesh%gauss%n_ws
       indexs = 0
       DO ms = 1, mesh%dom_mes
          js_loc = mesh%jjs(:,ms)

          DO ki = 1, 3
             DO ni = 1, nws
                i = js_loc(ni)
                iglob = vv_3_LA%loc_to_glob(ki,i)
                ix = (ki-1)*nws+ni
                idxms(ix) = iglob-1
             END DO
          END DO

          v145_loc = 0.d0
          v236_loc = 0.d0
          DO ls = 1, mesh%gauss%l_Gs
             indexs  = indexs + 1

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,js_loc)*mesh%gauss%wws(:,ls))

             !===Don't integrate on r=0
             IF (ray.LT.1.d-10) CYCLE

             smb(:) = opt_tensor_scaln_bdy(indexs,:)*ray*mesh%gauss%rjs(ls,ms)

             ki = 1
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145_loc(ix) = v145_loc(ix) + mesh%gauss%wws(ni,ls)*smb(1)
                v236_loc(ix) = v236_loc(ix) + mesh%gauss%wws(ni,ls)*smb(2)
             END DO

             ki = 2
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145_loc(ix) = v145_loc(ix) + mesh%gauss%wws(ni,ls)*smb(4)
                v236_loc(ix) = v236_loc(ix) + mesh%gauss%wws(ni,ls)*smb(3)
             END DO

             ki = 3
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145_loc(ix) = v145_loc(ix) + mesh%gauss%wws(ni,ls)*smb(5)
                v236_loc(ix) = v236_loc(ix) + mesh%gauss%wws(ni,ls)*smb(6)
             END DO


          END DO
          CALL VecSetValues(vb_145, 3*nws, idxms, v145_loc, ADD_VALUES, ierr)
          CALL VecSetValues(vb_236, 3*nws, idxms, v236_loc, ADD_VALUES, ierr)
       END DO
    END IF

    CALL VecAssemblyBegin(vb_145,ierr)
    CALL VecAssemblyEnd(vb_145,ierr)
    CALL VecAssemblyBegin(vb_236,ierr)
    CALL VecAssemblyEnd(vb_236,ierr)
  END SUBROUTINE rhs_3x3


  SUBROUTINE rhs_3x3_art_comp(mesh, vv_3_LA, mode, rhs_in, vb_145, vb_236, opt_tensor, opt_grad_div)
    USE def_type_mesh
    USE my_util
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                                                     :: mesh
    TYPE(petsc_csr_LA)                                                  :: vv_3_LA
    INTEGER,                                                 INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(mesh%dom_me*mesh%gauss%l_G,6),   INTENT(IN) :: rhs_in
    REAL(KIND=8), DIMENSION(3,mesh%dom_me*mesh%gauss%l_G,6), INTENT(IN), OPTIONAL :: opt_tensor
    REAL(KIND=8), DIMENSION(mesh%dom_me*mesh%gauss%l_G,2),   INTENT(IN), OPTIONAL :: opt_grad_div
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w)                           :: v145_loc, v236_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                           :: idxm
    INTEGER,      DIMENSION(mesh%gauss%n_w)                             :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)              :: dw_loc
    REAL(KIND=8), DIMENSION(3,6)                                        :: tensor_loc
    REAL(KIND=8), DIMENSION(2)                                          :: div_loc
    REAL(KIND=8), DIMENSION(6)                                          :: smb, f_tensor, f_div
    INTEGER      :: m, l, i, ni, index, nw, ix, ki, iglob, type, k
    REAL(KIND=8) :: ray
    Vec            :: vb_145, vb_236
    PetscErrorCode :: ierr

    CALL VecSet(vb_145, 0.d0, ierr)
    CALL VecSet(vb_236, 0.d0, ierr)

    nw = mesh%gauss%n_w
    index = 0
    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)

       DO ki = 1, 3
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm(ix) = iglob-1
          END DO
       END DO

       v145_loc = 0.d0
       v236_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index + 1
          dw_loc = mesh%gauss%dw(:,:,l,m)

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          IF (PRESENT(opt_tensor)) THEN
             DO TYPE = 1, 6
                DO k = 1, 3
                   tensor_loc(k,TYPE) = opt_tensor(k,index,TYPE)
                END DO
             END DO

             f_tensor(1) = (-mode*tensor_loc(1,4) + tensor_loc(2,3))/ray
             f_tensor(2) = (mode*tensor_loc(1,3) + tensor_loc(2,4))/ray
             f_tensor(3) = (-tensor_loc(1,3) - mode*tensor_loc(2,4))/ray
             f_tensor(4) = (-tensor_loc(1,4) + mode*tensor_loc(2,3))/ray
             f_tensor(5) = -tensor_loc(3,4)*mode/ray
             f_tensor(6) = tensor_loc(3,3)*mode/ray
          ELSE
             f_tensor = 0.d0
          END IF

          IF (PRESENT(opt_grad_div)) THEN
             DO TYPE = 1, 2
                div_loc(TYPE) = opt_grad_div(index,TYPE)
             END DO

             f_div(1) = div_loc(1)/ray
             f_div(2) = div_loc(2)/ray
             f_div(3) = -mode*div_loc(2)/ray
             f_div(4) = mode*div_loc(1)/ray
             f_div(5) = 0.d0
             f_div(6) = 0.d0
          ELSE
             f_div = 0.d0
          END IF

          smb = (rhs_in(index,:)+f_tensor+f_div)*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(1)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(2)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(4)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(3)
          END DO

          ki = 3
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v145_loc(ix) = v145_loc(ix) + mesh%gauss%ww(ni,l)*smb(5)
             v236_loc(ix) = v236_loc(ix) + mesh%gauss%ww(ni,l)*smb(6)
          END DO

          IF (PRESENT(opt_tensor)) THEN
             ki = 1
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(1,1) + dw_loc(2,ni)*tensor_loc(1,5))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(1,2) + dw_loc(2,ni)*tensor_loc(1,6))*ray*mesh%gauss%rj(l,m)
             END DO
             ki = 2
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(2,2) + dw_loc(2,ni)*tensor_loc(2,6))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(2,1) + dw_loc(2,ni)*tensor_loc(2,5))*ray*mesh%gauss%rj(l,m)
             END DO
             ki = 3
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(3,1) + dw_loc(2,ni)*tensor_loc(3,5))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*tensor_loc(3,2) + dw_loc(2,ni)*tensor_loc(3,6))*ray*mesh%gauss%rj(l,m)
             END DO
          END IF

          IF (PRESENT(opt_grad_div)) THEN
             ki = 1
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(1,ni)*div_loc(1))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(1,ni)*div_loc(2))*ray*mesh%gauss%rj(l,m)
             END DO
             ki = 2 !nothing to do
             ki = 3
             DO ni = 1, nw
                ix = (ki-1)*nw+ni
                v145_loc(ix) = v145_loc(ix) &
                     + (dw_loc(2,ni)*div_loc(1))*ray*mesh%gauss%rj(l,m)
                v236_loc(ix) = v236_loc(ix) &
                     + (dw_loc(2,ni)*div_loc(2))*ray*mesh%gauss%rj(l,m)
             END DO
          END IF

       ENDDO !l_G

       CALL VecSetValues(vb_145, 3*nw, idxm, v145_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_236, 3*nw, idxm, v236_loc, ADD_VALUES, ierr)

    ENDDO !dom_me


    CALL VecAssemblyBegin(vb_145,ierr)
    CALL VecAssemblyEnd(vb_145,ierr)
    CALL VecAssemblyBegin(vb_236,ierr)
    CALL VecAssemblyEnd(vb_236,ierr)
  END SUBROUTINE rhs_3x3_art_comp


END MODULE rhs_para_assembling
