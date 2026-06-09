MODULE assembling_rhs

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
    INTEGER      :: m, l, i, ni, index, nw, ix, ki, iglob, TYPE_VEC, k, ms, ls, indexs, nws
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
             DO TYPE_VEC = 1, 6
                DO k = 1, 3
                   tensor_loc(k,TYPE_VEC) = opt_tensor(k,index,TYPE_VEC)
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
    INTEGER      :: m, l, i, ni, index, nw, ix, ki, iglob, TYPE_VEC, k
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
             DO TYPE_VEC = 1, 6
                DO k = 1, 3
                   tensor_loc(k,TYPE_VEC) = opt_tensor(k,index,TYPE_VEC)
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
             DO TYPE_VEC = 1, 2
                div_loc(TYPE_VEC) = opt_grad_div(index,TYPE_VEC)
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

  SUBROUTINE rhs_3x3_weak_BC(mesh, vv_3_LA, mode, visco, stab_bdy_ns_Dirichlet, rhs_in, BC_in, &
             vb_145, vb_236, opt_tensor)
    USE def_type_mesh
    USE my_util
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                                                     :: mesh
    TYPE(petsc_csr_LA)                                                  :: vv_3_LA
    INTEGER,                                                 INTENT(IN) :: mode
    REAL(KIND=8),                                            INTENT(IN) :: visco, stab_bdy_ns_Dirichlet
    REAL(KIND=8), DIMENSION(mesh%dom_me*mesh%gauss%l_G,6),   INTENT(IN) :: rhs_in
    REAL(KIND=8), DIMENSION(mesh%dom_mes*mesh%gauss%l_Gs,6), INTENT(IN) :: BC_in
    REAL(KIND=8), DIMENSION(3,mesh%dom_me*mesh%gauss%l_G,6), INTENT(IN), OPTIONAL :: opt_tensor
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_w)                           :: v145_loc, v236_loc
    REAL(KIND=8), DIMENSION(3*mesh%gauss%n_ws)                          :: v145s_loc, v236s_loc
    INTEGER,      DIMENSION(3*mesh%gauss%n_w)                           :: idxm
    INTEGER,      DIMENSION(3*mesh%gauss%n_ws)                          :: idxms
    INTEGER,      DIMENSION(mesh%gauss%n_w)                             :: j_loc
    INTEGER,      DIMENSION(mesh%gauss%n_ws)                            :: js_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w)              :: dw_loc
    REAL(KIND=8), DIMENSION(3,6)                                        :: tensor_loc
    REAL(KIND=8), DIMENSION(6)                                          :: smb, f_tensor
    REAL(KIND=8), DIMENSION(6)                                          :: smb_norm, smb_tang
    INTEGER      :: m, l, i, ni, index, nw, ix, ki, iglob, TYPE_VEC, k, ms, ls, indexs, nws
    REAL(KIND=8) :: ray, z, hm1, x, stab_normal_bdy, stab_tangent_bdy, n1, n2
    Vec            :: vb_145, vb_236
    PetscErrorCode :: ierr

    CALL VecSet(vb_145, 0.d0, ierr)
    CALL VecSet(vb_236, 0.d0, ierr)

    !===Assemble volume terms
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
             DO TYPE_VEC = 1, 6
                DO k = 1, 3
                   tensor_loc(k,TYPE_VEC) = opt_tensor(k,index,TYPE_VEC)
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

    !===Assemble boundary terms
    ! Define stab_normal (always enforced) and stab_tangent
    stab_normal_bdy  = 1.d0 !1.d0 + visco
    stab_tangent_bdy = 1.d0 !visco

    nws = mesh%gauss%n_ws
    indexs = 0
    DO ms = 1, mesh%dom_mes
       js_loc = mesh%jjs(:,ms)
       !hm1 = inputs%stab_bdy_ns_Dirichlet/SUM(mesh%gauss%rjs(:,ms))
       hm1 = stab_bdy_ns_Dirichlet/SUM(mesh%gauss%rjs(:,ms))

       DO ki = 1, 3
          DO ni = 1, nws
             i = js_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nws+ni
             idxms(ix) = iglob-1
          END DO
       END DO
       v145s_loc = 0.d0
       v236s_loc = 0.d0

       IF (MINVAL(ABS(mesh%sides(ms) - inputs%vv_list_dirichlet_sides(1)%DIL))==0) THEN
          DO ls = 1, mesh%gauss%l_Gs
             indexs = indexs + 1

             !===Compute radius of Gauss point
             ray = SUM(mesh%rr(1,js_loc)*mesh%gauss%wws(:,ls))

             !===Don't integrate on r=0
             IF (ray.LT.1.d-10) CYCLE

             x = stab_normal_bdy*hm1*mesh%gauss%rjs(ls,ms)*ray
             z = stab_tangent_bdy*hm1*mesh%gauss%rjs(ls,ms)*ray
             n1 = mesh%gauss%rnorms(1,ls,ms)
             n2 = mesh%gauss%rnorms(2,ls,ms)

             !Normal terms: (BC_in . n) * (v . n)
             smb_norm(1) = x*(n1**2*BC_in(indexs,1) + n1*n2*BC_in(indexs,5))
             smb_norm(2) = x*(n1**2*BC_in(indexs,2) + n1*n2*BC_in(indexs,6))
             smb_norm(3) = 0.d0
             smb_norm(4) = 0.d0
             smb_norm(5) = x*(n1*n2*BC_in(indexs,1) +  n2**2*BC_in(indexs,5))
             smb_norm(6) = x*(n1*n2*BC_in(indexs,2) +  n2**2*BC_in(indexs,6))
             !Tangential terms: (BC_in x n ) . (v x n)
             smb_tang(1) = z*(n2**2*BC_in(indexs,1) - n1*n2*BC_in(indexs,5))
             smb_tang(2) = z*(n2**2*BC_in(indexs,2) - n1*n2*BC_in(indexs,6))
             smb_tang(3) = z*(n1**2+n2**2)*BC_in(indexs,3)
             smb_tang(4) = z*(n1**2+n2**2)*BC_in(indexs,4)
             smb_tang(5) = z*(-n1*n2*BC_in(indexs,1) + n1**2*BC_in(indexs,5))
             smb_tang(6) = z*(-n1*n2*BC_in(indexs,2) + n1**2*BC_in(indexs,6))

             ki = 1
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145s_loc(ix) = v145s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(1)+smb_tang(1))
                v236s_loc(ix) = v236s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(2)+smb_tang(2))
             END DO

             ki = 2
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145s_loc(ix) = v145s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(4)+smb_tang(4))
                v236s_loc(ix) = v236s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(3)+smb_tang(3))
             END DO

             ki = 3
             DO ni = 1, nws
                ix = (ki-1)*nws+ni
                v145s_loc(ix) = v145s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(5)+smb_tang(5))
                v236s_loc(ix) = v236s_loc(ix) + mesh%gauss%wws(ni,ls)*(smb_norm(6)+smb_tang(6))
             END DO

          END DO !l_Gs
          CALL VecSetValues(vb_145, 3*nws, idxms, v145s_loc, ADD_VALUES, ierr)
          CALL VecSetValues(vb_236, 3*nws, idxms, v236s_loc, ADD_VALUES, ierr)
       ELSE !No dirichlet face so we only update indexs
          DO ls = 1, mesh%gauss%l_Gs
             indexs  = indexs + 1
          END DO
       END IF
    END DO !dom_mes

    CALL VecAssemblyBegin(vb_145,ierr)
    CALL VecAssemblyEnd(vb_145,ierr)
    CALL VecAssemblyBegin(vb_236,ierr)
    CALL VecAssemblyEnd(vb_236,ierr)
  END SUBROUTINE rhs_3x3_weak_BC

  SUBROUTINE rhs_pp_consistent_pressure(type_fe_velocity, vv_mesh, pp_mesh, pp_1_LA, mode, visco, mass, &
             ff_gauss, curl_vel_gauss, error_vel_gauss, vel_exact_gauss, pb_1, pb_2)
    USE def_type_mesh
    USE basis_change
    USE my_util
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type)                          :: vv_mesh, pp_mesh
    TYPE(petsc_csr_LA)                       :: pp_1_LA
    INTEGER,                      INTENT(IN)  :: type_fe_velocity
    INTEGER,                      INTENT(IN) :: mode
    REAL(KIND=8),                 INTENT(IN) :: visco, mass
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: ff_gauss !=source_NS + 1/dt*un_m1 - curl(un)xun
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: curl_vel_gauss  !=curl(un)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: error_vel_gauss !=u_exact-un
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: vel_exact_gauss !=u_exact
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w)                   :: p1_loc, p2_loc
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_ws)                  :: p1s_loc, p2s_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_w)                   :: idxm
    INTEGER,      DIMENSION(pp_mesh%gauss%n_ws)                  :: idxms
    INTEGER,      DIMENSION(vv_mesh%gauss%n_w)                   :: ju_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_w)                   :: jp_loc
    INTEGER,      DIMENSION(vv_mesh%gauss%n_ws)                  :: jus_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_ws)                  :: jps_loc
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dwu_loc
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%k_d,pp_mesh%gauss%n_w) :: dwp_loc
    REAL(KIND=8), DIMENSION(vv_mesh%gauss%k_d,vv_mesh%gauss%n_w) :: dwus_loc
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%k_d,pp_mesh%gauss%n_w) :: dwps_loc
    INTEGER      :: m, l, n, i, ni, index, iglob, ms, ls, indexs
    INTEGER      :: nwu, nwus, nwp, nwps
    REAL(KIND=8) :: ray, z, hm1, x, y, stab_normal_bdy, stab_tangent_bdy, n1, n2
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w,vv_mesh%gauss%l_G)   ::  wp
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_ws,vv_mesh%gauss%l_Gs) ::  wps
    !===Change of basis arrays in save status
    REAL(KIND=8), DIMENSION(:,:), POINTER, SAVE :: aij_p1p2, aij_p2p3
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
    REAL(KIND=8), DIMENSION(:,:), POINTER, SAVE :: aij_p1p2_s, aij_p2p3_s
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij_s
    LOGICAL :: once_p1p2=.TRUE., once_p2p3=.TRUE.
    Vec            :: pb_1, pb_2
    PetscErrorCode :: ierr

    !===Construct aij_p2, aij_p3
    IF (type_fe_velocity==2) THEN
       IF (once_p1p2) THEN
          CALL p1_p2(aij_p1p2)
          CALL p1_p2_s(aij_p1p2_s)
          once_p1p2=.FALSE.
       END IF
       aij => aij_p1p2
       aij_s => aij_p1p2_s
    ELSE IF (type_fe_velocity==3) THEN
       IF (once_p2p3) THEN
          CALL p2_p3(aij_p2p3)
          CALL p2_p3_s(aij_p2p3_s)
          once_p2p3=.FALSE.
       END IF
       aij => aij_p2p3
       aij_s => aij_p2p3_s
    ELSE
       CALL error_petsc('qs_01_div_hybrid_generic, type_fe_velocity not correct')
    END IF

    CALL VecSet(pb_1, 0.d0, ierr)
    CALL VecSet(pb_2, 0.d0, ierr)

    nwu = vv_mesh%gauss%n_w
    nwp = pp_mesh%gauss%n_w
    nwus = vv_mesh%gauss%n_ws
    nwps = pp_mesh%gauss%n_ws

    !===Construct wp and wps
    DO l = 1, vv_mesh%gauss%l_G
       DO n = 1, nwp
          wp(n,l) = SUM(aij(n,:)*vv_mesh%gauss%ww(:,l))
       END DO
    END DO
    DO l = 1, vv_mesh%gauss%l_Gs
       DO n = 1, nwps
          wps(n,l) = SUM(aij_s(n,:)*vv_mesh%gauss%wws(:,l))
       END DO
    END DO
    !===End Construct wp and wps

    !===Assemble volume terms
    index = 0
    DO m = 1, vv_mesh%dom_me
       ju_loc = vv_mesh%jj(:,m)
       jp_loc = pp_mesh%jj(:,m)

       DO ni = 1, nwp
          i = jp_loc(ni)
          iglob = pp_1_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       p1_loc = 0.d0
       p2_loc = 0.d0
       DO l = 1, vv_mesh%gauss%l_G
          index  = index + 1
          dwu_loc = vv_mesh%gauss%dw(:,:,l,m)
          DO ni = 1, nwp
             dwp_loc(1,ni) = SUM(aij(ni,:)*dwu_loc(1,:))
             dwp_loc(2,ni) = SUM(aij(ni,:)*dwu_loc(2,:))
          END DO

          !===Compute radius of Gauss point
          ray = SUM(vv_mesh%rr(1,ju_loc)*vv_mesh%gauss%ww(:,l))

          !===Compute ff_gauss * GRAD(q)
          x = vv_mesh%gauss%rj(l,m)*ray
          DO ni = 1, nwp
             p1_loc(ni) = p1_loc(ni) + x * (ff_gauss(index,1)*dwp_loc(1,ni) &
                             - ff_gauss(index,4)*wp(ni,l)*(mode/ray) &
                             + ff_gauss(index,5)*dwp_loc(2,ni))
             p2_loc(ni) = p2_loc(ni) + x * (ff_gauss(index,2)*dwp_loc(1,ni) &
                             + ff_gauss(index,3)*wp(ni,l)*(mode/ray) &
                             + ff_gauss(index,6)*dwp_loc(2,ni))
          END DO
       ENDDO !l_G

       CALL VecSetValues(pb_1, nwp, idxm, p1_loc, ADD_VALUES, ierr)
       CALL VecSetValues(pb_2, nwp, idxm, p2_loc, ADD_VALUES, ierr)

    ENDDO !dom_me

    !===Assemble boundary terms
    ! Define stab_normal (always enforced) and stab_tangent
    stab_normal_bdy  = 1.d0 !1.d0 + visco
    stab_tangent_bdy = 1.d0 !visco

    indexs = 0
    DO ms = 1, vv_mesh%dom_mes
       m = vv_mesh%neighs(ms)
       ju_loc  = vv_mesh%jj(:,m) !for term with derivative of test function
       jus_loc = vv_mesh%jjs(:,ms)
       jp_loc  = pp_mesh%jj(:,m) !for term with derivative of test function
       jps_loc = pp_mesh%jjs(:,ms)
       hm1 = inputs%stab_bdy_ns_Dirichlet/SUM(vv_mesh%gauss%rjs(:,ms))

       !===For terms without derivative of test function
       DO ni = 1, nwps
          i = jps_loc(ni)
          iglob = pp_1_LA%loc_to_glob(1,i)
          idxms(ni) = iglob-1
       END DO
       p1s_loc = 0.d0
       p2s_loc = 0.d0

       !===For terms with derivative of test function
       DO ni = 1, nwp
          i = jp_loc(ni)
          iglob = pp_1_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO
       p1_loc = 0.d0
       p2_loc = 0.d0

       DO ls = 1, vv_mesh%gauss%l_Gs
          indexs = indexs + 1
          dwus_loc = vv_mesh%gauss%dw_s(:,:,ls,ms)
          DO ni = 1, nwp
             dwps_loc(1,ni) = SUM(aij(ni,:)*dwus_loc(1,:))
             dwps_loc(2,ni) = SUM(aij(ni,:)*dwus_loc(2,:))
          END DO

          !===Compute radius of Gauss point
          ray = SUM(vv_mesh%rr(1,jus_loc)*vv_mesh%gauss%wws(:,ls))

          !===Don't integrate on r=0
          IF (ray.LT.1.d-10) CYCLE

          !===Get normal vector on Gauss point
          n1 = vv_mesh%gauss%rnorms(1,ls,ms)
          n2 = vv_mesh%gauss%rnorms(2,ls,ms)

          x = vv_mesh%gauss%rjs(ls,ms)*ray
          y = stab_normal_bdy*hm1*x
          z = stab_tangent_bdy*hm1*x

          !===Compute -visco * CURL(u).(GRAD(q) x n) on all boundary faces
          DO ni = 1, nwps
             p1s_loc(ni) = p1s_loc(ni) &
                      - visco*x*( &
                         - n2*curl_vel_gauss(indexs,2)*wps(ni,ls)*(mode/ray) &
                         + n1*curl_vel_gauss(indexs,6)*wps(ni,ls)*(mode/ray) &
                      )
             p2s_loc(ni) = p2s_loc(ni) &
                      - visco*x*( &
                         + n2*curl_vel_gauss(indexs,1)*wps(ni,ls)*(mode/ray) &
                         - n1*curl_vel_gauss(indexs,5)*wps(ni,ls)*(mode/ray) &
                      )
          END DO

          DO ni = 1, nwp
             p1_loc(ni) = p1_loc(ni) &
                    - visco*x*( &
                         curl_vel_gauss(indexs,3)*(n1*dwps_loc(2,ni) - n2*dwps_loc(1,ni)) &
                    )
             p2_loc(ni) = p2_loc(ni) &
                    - visco*x*( &
                         curl_vel_gauss(indexs,4)*(n1*dwps_loc(2,ni) - n2*dwps_loc(1,ni)) &
                    )
          END DO

          !===Add contribution specific to Dirichlet faces
!TEST LC DEBUG CHECK IF NEEDED (add extra END IF later if needed)
          !IF (SIZE(inputs%vv_list_dirichlet_sides(1)%DIL) > 0) THEN
          !   IF (MINVAL(ABS(mesh%sides(ms) - inputs%vv_list_dirichlet_sides(1)%DIL))<1) THEN
!TEST LC DEBUG CHECK IF NEEDED
          IF (MINVAL(ABS(vv_mesh%sides(ms) - inputs%vv_list_dirichlet_sides(1)%DIL))==0) THEN
             DO ni = 1, nwps
                p1s_loc(ni) = p1s_loc(ni) &
                     !===Compute - mass * (vel_exact . n) q
                     - mass * x * wps(ni,ls)* &
                       (vel_exact_gauss(indexs,1)*n1 + vel_exact_gauss(indexs,5)*n2) &
                     !===Compute Tangent term: stab_tangent/h * (error_vel x n).(GRAD(q) x n)
                     + z * (-mode)*(n2**2+n1**2)*error_vel_gauss(indexs,4)*wps(ni,ls)/ray

                p2s_loc(ni) = p2s_loc(ni) &
                     !===Compute - mass * (vel_exact . n) q
                     - mass * x * wps(ni,ls)* &
                       (vel_exact_gauss(indexs,2)*n1 + vel_exact_gauss(indexs,6)*n2) &
                     !===Compute Tangent term: stab_tangent/h * (error_vel x n).(GRAD(q) x n)
                     + z * mode*(n2**2+n1**2)*error_vel_gauss(indexs,3)*wps(ni,ls)/ray
             END DO

             DO ni = 1, nwp
                p1_loc(ni) = p1_loc(ni) &
                     !===Compute Normal term: stab_normal/h * (error_vel . n) * (GRAD(q) . n)
                     + y * (n1*error_vel_gauss(indexs,1)+n2*error_vel_gauss(indexs,5)) &
                             *(n1*dwps_loc(1,ni) + n2*dwps_loc(2,ni)) &
                     !===Compute Tangent term: stab_tangent/h * (error_vel x n).(GRAD(q) x n)
                     + z * (n1*error_vel_gauss(indexs,5) - n2*error_vel_gauss(indexs,1)) &
                               *(n1*dwps_loc(2,ni) - n2*dwps_loc(1,ni))

                p2_loc(ni) = p2_loc(ni) &
                     !===Compute Normal term: stab_normal/h * (error_vel . n) * (GRAD(q) . n)
                     + y * (n1*error_vel_gauss(indexs,2)+n2*error_vel_gauss(indexs,6)) &
                            *(n1*dwps_loc(1,ni) + n2*dwps_loc(2,ni)) &
                     !===Compute Tangent term: stab_tangent/h * (error_vel x n).(GRAD(q) x n)
                     + z * (n1*error_vel_gauss(indexs,6) - n2*error_vel_gauss(indexs,2)) &
                               *(n1*dwps_loc(2,ni) - n2*dwps_loc(1,ni))
             END DO
          END IF
          !TEST LC DEBUG: add extra END IF if needed (see comment above)
       END DO !l_Gs 
       !Contribution without derivative in r and z
       CALL VecSetValues(pb_1, nwps, idxms, p1s_loc, ADD_VALUES, ierr)
       CALL VecSetValues(pb_2, nwps, idxms, p2s_loc, ADD_VALUES, ierr)
       !Contribution with derivative in r or z
       CALL VecSetValues(pb_1, nwp, idxm, p1_loc, ADD_VALUES, ierr)
       CALL VecSetValues(pb_2, nwp, idxm, p2_loc, ADD_VALUES, ierr)
    END DO !dom_mes

    CALL VecAssemblyBegin(pb_1,ierr)
    CALL VecAssemblyEnd(pb_1,ierr)
    CALL VecAssemblyBegin(pb_2,ierr)
    CALL VecAssemblyEnd(pb_2,ierr)

  END SUBROUTINE rhs_pp_consistent_pressure

END MODULE assembling_rhs
