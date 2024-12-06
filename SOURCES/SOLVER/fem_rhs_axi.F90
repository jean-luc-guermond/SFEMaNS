MODULE fem_rhs_axi
#include "petsc/finclude/petsc.h"
  USE petsc
  USE my_util
  USE input_data
CONTAINS

  SUBROUTINE qs_ns_stab_new(mesh,vv_1_LA,vv_2_LA,mode,ff,vel_tot,V1m,vit,P,dudt,phalf,nlhalf,dt,&
       vb_2_14,vb_2_23,vb_1_5,vb_1_6,rotv_v, vel_tot_max)
    !=================================
    !RHS for Navier-Stokes
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    !USE sub_plot
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_1_LA,vv_2_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v, dudt, nlhalf
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: vel_tot      !(noeud)
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, vit     !V(noeud, type)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P, phalf
    REAL(KIND=8),                               INTENT(IN) :: dt
    INTEGER,                                    INTENT(IN) :: mode
    !TEST LES LC October 28, 2014
    REAL(KIND=8)                 :: vel_tot_max
    !TEST LES LC October 28, 2014
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, fv, mult
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: vitloc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%dom_me)                   :: visc_plot
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8), DIMENSION(6)   :: visc1, visc2
    REAL(KIND=8)   :: ray, div1, div2, cfl, vloc, normal_vit
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: h
    REAL(KIND=8),                            SAVE :: coeff_ed_st, R_eff, &
         visc_eff, surf, nu_loc, coeff_visc_ordre_un
    LOGICAL,                                 SAVE :: once = .TRUE.
    REAL(KIND=8), DIMENSION(6,mesh%dom_me) :: viscosity
    REAL(KIND=8), DIMENSION(6) :: norm_vit
    REAL(KIND=8) ::  type_fe
    INTEGER :: m, l , i , ni , index, index2, TYPE, k
    INTEGER :: ms, nw, ix, ki, iglob
!!$ WARNING FL Variables removed (used for edge_stab)
!!$    INTEGER :: cotei, ls
!!$    REAL(KIND=8) :: dul, ed_st, h2
!!$    INTEGER, DIMENSION(mesh%gauss%n_w,2) :: jji_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2) :: u0loci, uloci
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%l_Gs, 2) :: dwni_loc
!!$ WARNING FL Variables removed (used for edge_stab)
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_2_14,vb_2_23,vb_1_5,vb_1_6
    PetscErrorCode :: ierr

    CALL VecSet(vb_2_14, 0.d0, ierr)
    CALL VecSet(vb_2_23, 0.d0, ierr)
    CALL VecSet(vb_1_5, 0.d0, ierr)
    CALL VecSet(vb_1_6, 0.d0, ierr)

    IF (once) THEN
       once =.FALSE.

       IF (.NOT.inputs%LES) THEN
          inputs%LES_coeff1=0.d0
          inputs%LES_coeff2=0.d0
          inputs%LES_coeff3=0.d0
          inputs%LES_coeff4=0.d0
       END IF

       surf = 0.d0
       DO m = 1, mesh%dom_me
          DO l = 1, mesh%gauss%l_G
             ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
             surf = surf + mesh%gauss%rj(l,m)*ray
          END DO
       END DO

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE
          type_fe = 2
       END IF
       coeff_ed_st         = inputs%LES_coeff3*0.02d0/type_fe
       coeff_visc_ordre_un = inputs%LES_coeff4
       IF (mesh%edge_stab) THEN
          ALLOCATE(h(mesh%mi))
          DO ms = 1, mesh%mi
             h(ms)  = SUM(mesh%gauss%rji(:,ms))/type_fe
          END DO
       END IF
    END IF

    !ATTENTION: JLG Jan 25 2010
    !ATTENTION: inputs%LES_coeff1 is assumed to be of the order of the convective velocity
    !that simplifies the semi-implicit treatment of the LES viscosity
!!$    normal_vit = MAXVAL(vel_tot)
    !TEST LES LC October 28, 2014
    normal_vit = vel_tot_max
    !TEST LES LC October 28, 2014
    DO TYPE = 1, 6
       norm_vit(TYPE) = SUM(ABS(vit(:,TYPE)))/mesh%np + 1.d-14
    END DO
    !ATTENTION: JLG Jan 25 2010
    R_eff = 0.d0
    cfl = 0
    index = 0
    index2 = 0
    IF (inputs%LES) THEN
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          vloc = MAXVAL(vel_tot(j_loc))
          cfl = MAX(vloc*dt/mesh%hloc(m),cfl)
          visc1 = 0
          visc2 = 0
          DO l = 1, mesh%gauss%l_G
             index2  = index2 +1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !--------radius of gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
             !fs Source term
             !ft Time derivative
             !fp Pressure gradient
             !--------calcul de la premiere composante
             ft(1) = SUM(dudt(j_loc,1) *mesh%gauss%ww(:,l))
             fp(1) = SUM(phalf(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante
             ft(2) = SUM(dudt(j_loc,2) * mesh%gauss%ww(:,l))
             fp(2) = SUM(phalf(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante
             ft(3) = SUM(dudt(j_loc,3) *mesh%gauss%ww(:,l))
             fp(3) = SUM(phalf(j_loc,2)*mesh%gauss%ww(:,l))*mode/ray
             !--------calcul de la quatrieme composante
             ft(4) = SUM(dudt(j_loc,4) *mesh%gauss%ww(:,l))
             fp(4) = -SUM(phalf(j_loc,1)*mesh%gauss%ww(:,l))*mode/ray
             !--------calcul de la cinquieme composante
             ft(5) = SUM(dudt(j_loc,5) *mesh%gauss%ww(:,l))
             fp(5) = SUM(phalf(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante
             ft(6) = SUM(dudt(j_loc,6) *mesh%gauss%ww(:,l))
             fp(6) = SUM(phalf(j_loc,2)*dw_loc(2,:))
             !-------calcul du second membre pour le terme nonlineaire------------------------

             visc1 = MAX(visc1,ABS(ft+fp+nlhalf(index2,:)))

             !--------Calcul du gradient de la vitesse
             DO TYPE = 1, 6
                DO k = 1 ,2
                   grad(k,TYPE,l) = SUM(vit(j_loc,TYPE)*dw_loc(k,:))
                END DO
                vitloc(TYPE,l) = SUM(vit(j_loc,TYPE)*mesh%gauss%ww(:,l))
             END DO

             !--------Calcul de la divergence
             div1 = ABS(grad(1,1,l) + vitloc(1,l)/ray + grad(2,5,l) + vitloc(4,l)*mode/ray)
             div2 = ABS(grad(1,2,l) + vitloc(2,l)/ray + grad(2,6,l) - vitloc(3,l)*mode/ray)
             visc2(1) = MAX(visc2(1),div1)
             visc2(2) = MAX(visc2(2),div2)
          END DO
          visc2(4) = visc2(1); visc2(5) = visc2(1)
          visc2(3) = visc2(2); visc2(6) = visc2(2)

          nu_loc = 0.d0
          DO TYPE = 1, 6
             !visc1(type) = MAX(visc1(type)/normal_vit,visc2(type))
             !visc1(type) = MAX(visc1(type),visc2(type)*normal_vit)/MAXVAL(ABS(vitloc(type,:)))
             visc1(TYPE) = MAX(visc1(TYPE),2*visc2(TYPE)*normal_vit)/norm_vit(TYPE)
             visc_eff = mesh%hloc(m)*MIN(coeff_visc_ordre_un*normal_vit,inputs%LES_coeff2*mesh%hloc(m)*visc1(TYPE))
             nu_loc = nu_loc + visc_eff
             !======Semi-implicit version==========
             viscosity(TYPE,m) = inputs%LES_coeff1*mesh%hloc(m) - visc_eff
             !======Semi-implicit version==========
          END DO
          R_eff = R_eff + (nu_loc + dt**2*mesh%hloc(m))*mesh%hloc(m)**2*ray
          visc_plot(m) = (nu_loc/6)/(coeff_visc_ordre_un*mesh%hloc(m)*normal_vit)
       END DO
    ELSE
       cfl = 0.d0
       viscosity = 0.d0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          vloc = maxval(vel_tot(j_loc))
          cfl = max(vloc*dt/mesh%hloc(m),cfl)
       END DO
    END IF


    !DO type = 1, 6
    !   CALL average(mesh,viscosity(type,:))
    !END DO

    nw = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       mult(:)= viscosity(:,m)
       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_1_LA%loc_to_glob(1,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_2_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)
          DO TYPE = 1, 6
             DO k = 1 ,2
                grad(k,TYPE,l) = SUM(vit(j_loc,TYPE)*dw_loc(k,:))
             END DO
             vitloc(TYPE,l) = SUM(vit(j_loc,TYPE)*mesh%gauss%ww(:,l))
          END DO

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fv(1) = ((mode*vitloc(1,l)+vitloc(4,l))*mode +mode*vitloc(4,l)+vitloc(1,l))/ray**2
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2) = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:))
          fv(2) = ((mode*vitloc(2,l)-vitloc(3,l))*mode -mode*vitloc(3,l)+vitloc(2,l))/ray**2
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3) = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3) = -SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          fv(3) = (-mode*vitloc(2,l)+vitloc(3,l) +(mode*vitloc(3,l)-vitloc(2,l))*mode)/ray**2
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4) = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4) = SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          fv(4) =  (mode*vitloc(1,l)+vitloc(4,l) +(mode*vitloc(4,l)+vitloc(1,l))*mode)/ray**2
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5) = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:))
          fv(5) =  vitloc(5,l)*(mode/ray)**2
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6) = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:))
          fv(6) = vitloc(6,l)*(mode/ray)**2
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          fv = mult*fv

          smb =  (ft+fp+fs+fv-rotv_v(index,:))*ray*mesh%gauss%rj(l,m)
          DO TYPE = 1, 6
             grad(:,TYPE,l) =  mult(TYPE)*grad(:,TYPE,l)*ray*mesh%gauss%rj(l,m)
          END DO

!!$        DO j=1,6
!!$           DO ni = 1, mesh%gauss%n_w
!!$              u0(j_loc(ni),j) = u0(j_loc(ni),j) +  mesh%gauss%ww(ni,l)*smb(j) + SUM(dw_loc(:,ni)*grad(:,j,l))
!!$           ENDDO
!!$        ENDDO

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(1) + SUM(dw_loc(:,ni)*grad(:,1,l))
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(2) + SUM(dw_loc(:,ni)*grad(:,2,l))
             v5_loc(ix)  = v5_loc(ix)  + mesh%gauss%ww(ni,l)*smb(5) + SUM(dw_loc(:,ni)*grad(:,5,l))
             v6_loc(ix)  = v6_loc(ix)  + mesh%gauss%ww(ni,l)*smb(6) + SUM(dw_loc(:,ni)*grad(:,6,l))
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(4) + SUM(dw_loc(:,ni)*grad(:,4,l))
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(3) + SUM(dw_loc(:,ni)*grad(:,3,l))
          END DO

       ENDDO

       CALL VecSetValues(vb_2_14, 2*nw, idxm_2, v14_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_2_23, 2*nw, idxm_2, v23_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_5,    nw, idxm_1,  v5_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_6,    nw, idxm_1,  v6_loc, ADD_VALUES, ierr)
    ENDDO

    !WRITE(*,*) ' CFL = ', cfl, ' R_eff for mode', mode, normal_vit*6*surf/R_eff
    !IF (mode==0) THEN
    !CALL plot_const_p1_label(mesh%jj, mesh%rr, visc_plot, 'tttt_0.plt')
    !END IF

    IF (inputs%LES) THEN
       IF (mesh%edge_stab) THEN
          CALL error_Petsc('BUG in qs_ns_stab_new: terms with edge_stab not yet assembled')
!!$     ed_st = normal_vit*coeff_ed_st
!!$     DO ms = 1, mesh%mi
!!$        dwni_loc = mesh%gauss%dwni(:,:,:,ms)
!!$        jji_loc = mesh%jji(:,:,ms)
!!$        h2 = -ed_st*h(ms)**2
!!$
!!$        j_loc(1:mesh%gauss%n_ws) = mesh%jjsi(1:mesh%gauss%n_ws,ms)
!!$        DO TYPE = 1, 6
!!$           DO cotei = 1, 2
!!$              uloci(:,cotei) = vit(jji_loc(:,cotei),TYPE)
!!$           END DO
!!$
!!$           u0loci = 0.d0
!!$           DO ls = 1, mesh%gauss%l_Gs
!!$              !===Compute radius of Gauss point at face center
!!$              ray = mesh%rr(1,j_loc(3))
!!$
!!$              !--------Calcul du saut de la derivee normale de la vitesse
!!$              dul = SUM(dwni_loc(:,ls,:)*uloci)*mesh%gauss%rji(ls,ms)*h2*ray
!!$              DO cotei = 1, 2
!!$                 DO ni = 1, mesh%gauss%n_w
!!$                    u0loci(ni, cotei) =  u0loci(ni, cotei) + dwni_loc(ni,ls,cotei)*dul
!!$                 END DO
!!$              END DO
!!$           END DO
!!$
!!$           DO cotei = 1, 2
!!$              DO ni = 1, mesh%gauss%n_w
!!$                 u0(jji_loc(ni,cotei),TYPE) = u0(jji_loc(ni,cotei),TYPE) + u0loci(ni,cotei)
!!$              END DO
!!$           END DO
!!$        END DO

       END IF
    END IF

    CALL VecAssemblyBegin(vb_2_14,ierr)
    CALL VecAssemblyEnd(vb_2_14,ierr)
    CALL VecAssemblyBegin(vb_2_23,ierr)
    CALL VecAssemblyEnd(vb_2_23,ierr)
    CALL VecAssemblyBegin(vb_1_5,ierr)
    CALL VecAssemblyEnd(vb_1_5,ierr)
    CALL VecAssemblyBegin(vb_1_6,ierr)
    CALL VecAssemblyEnd(vb_1_6,ierr)

  END SUBROUTINE qs_ns_stab_new

  SUBROUTINE qs_ns_stab_3x3(mesh,vv_3_LA,mode,ff,vel_tot,V1m,vit,P,dudt,phalf,nlhalf,dt,&
       vb_145, vb_236,rotv_v, vel_tot_max)
    !=================================
    !RHS for Navier-Stokes
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    USE rhs_para_assembling
    !USE sub_plot
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_3_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotv_v, dudt, nlhalf
    REAL(KIND=8), DIMENSION(:),                 INTENT(IN) :: vel_tot      !(noeud)
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, vit     !V(noeud, type)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P, phalf
    REAL(KIND=8),                               INTENT(IN) :: dt
    INTEGER,                                    INTENT(IN) :: mode
    !TEST LES LC October 28, 2014
    REAL(KIND=8)                 :: vel_tot_max
    !TEST LES LC October 28, 2014
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, fv, mult
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: vitloc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(mesh%dom_me)                   :: visc_plot
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8), DIMENSION(mesh%dom_me*mesh%gauss%l_G,6)  :: rhs_gauss
    REAL(KIND=8), DIMENSION(6)   :: visc1, visc2
    REAL(KIND=8)   :: ray, div1, div2, cfl, vloc, normal_vit

    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), SAVE :: h
    REAL(KIND=8),                            SAVE :: coeff_ed_st, R_eff, &
         visc_eff, surf, nu_loc, coeff_visc_ordre_un
    LOGICAL,                                 SAVE :: once = .TRUE.
    REAL(KIND=8), DIMENSION(6,mesh%dom_me) :: viscosity
    REAL(KIND=8), DIMENSION(6) :: norm_vit
    REAL(KIND=8) ::  type_fe
    INTEGER :: m, l , i , ni , index, index2, TYPE, k
    INTEGER :: ms, nw, ix, ki, iglob
!!$ WARNING FL Variables removed (used for edge_stab)
!!$    INTEGER :: cotei, ls
!!$    REAL(KIND=8) :: dul, ed_st, h2
!!$    INTEGER, DIMENSION(mesh%gauss%n_w,2) :: jji_loc
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,2) :: u0loci, uloci
!!$    REAL(KIND=8), DIMENSION(mesh%gauss%n_w,mesh%gauss%l_Gs, 2) :: dwni_loc
!!$ WARNING FL Variables removed (used for edge_stab)
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_145, vb_236
    !PetscErrorCode :: ierr

    !CALL VecSet(vb_3_145, 0.d0, ierr)
    !CALL VecSet(vb_3_236, 0.d0, ierr)

    IF (once) THEN
       once =.FALSE.

       IF (.NOT.inputs%LES) THEN
          inputs%LES_coeff1=0.d0
          inputs%LES_coeff2=0.d0
          inputs%LES_coeff3=0.d0
          inputs%LES_coeff4=0.d0
       END IF

       surf = 0.d0
       DO m = 1, mesh%dom_me
          DO l = 1, mesh%gauss%l_G
             ray = SUM(mesh%rr(1,mesh%jj(:,m))*mesh%gauss%ww(:,l))
             surf = surf + mesh%gauss%rj(l,m)*ray
          END DO
       END DO

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE
          type_fe = 2
       END IF
       coeff_ed_st         = inputs%LES_coeff3*0.02d0/type_fe
       coeff_visc_ordre_un = inputs%LES_coeff4
       IF (mesh%edge_stab) THEN
          ALLOCATE(h(mesh%mi))
          DO ms = 1, mesh%mi
             h(ms)  = SUM(mesh%gauss%rji(:,ms))/type_fe
          END DO
       END IF
    END IF

    !ATTENTION: JLG Jan 25 2010
    !ATTENTION: inputs%LES_coeff1 is assumed to be of the order of the convective velocity
    !that simplifies the semi-implicit treatment of the LES viscosity
!!$    normal_vit = MAXVAL(vel_tot)
    !TEST LES LC October 28, 2014
    normal_vit = vel_tot_max
    !TEST LES LC October 28, 2014
    DO TYPE = 1, 6
       norm_vit(TYPE) = SUM(ABS(vit(:,TYPE)))/mesh%np + 1.d-14
    END DO
    !ATTENTION: JLG Jan 25 2010
    R_eff = 0.d0
    cfl = 0
    index = 0
    index2 = 0
    IF (inputs%LES) THEN
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          vloc = MAXVAL(vel_tot(j_loc))
          cfl = MAX(vloc*dt/mesh%hloc(m),cfl)
          visc1 = 0
          visc2 = 0
          DO l = 1, mesh%gauss%l_G
             index2  = index2 +1
             dw_loc = mesh%gauss%dw(:,:,l,m)

             !--------radius of gauss point
             ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))
             !fs Source term
             !ft Time derivative
             !fp Pressure gradient
             !--------calcul de la premiere composante
             ft(1) = SUM(dudt(j_loc,1) *mesh%gauss%ww(:,l))
             fp(1) = SUM(phalf(j_loc,1)*dw_loc(1,:))
             !--------calcul de la seconde composante
             ft(2) = SUM(dudt(j_loc,2) * mesh%gauss%ww(:,l))
             fp(2) = SUM(phalf(j_loc,2)*dw_loc(1,:))
             !--------calcul de la troisieme composante
             ft(3) = SUM(dudt(j_loc,3) *mesh%gauss%ww(:,l))
             fp(3) = SUM(phalf(j_loc,2)*mesh%gauss%ww(:,l))*mode/ray
             !--------calcul de la quatrieme composante
             ft(4) = SUM(dudt(j_loc,4) *mesh%gauss%ww(:,l))
             fp(4) = -SUM(phalf(j_loc,1)*mesh%gauss%ww(:,l))*mode/ray
             !--------calcul de la cinquieme composante
             ft(5) = SUM(dudt(j_loc,5) *mesh%gauss%ww(:,l))
             fp(5) = SUM(phalf(j_loc,1)*dw_loc(2,:))
             !--------calcul de la sixieme composante
             ft(6) = SUM(dudt(j_loc,6) *mesh%gauss%ww(:,l))
             fp(6) = SUM(phalf(j_loc,2)*dw_loc(2,:))
             !-------calcul du second membre pour le terme nonlineaire------------------------

             visc1= MAX(visc1,ABS(ft+fp+nlhalf(index2,:)))

             !--------Calcul du gradient de la vitesse
             DO TYPE = 1, 6
                DO k = 1 ,2
                   grad(k,TYPE,l) = SUM(vit(j_loc,TYPE)*dw_loc(k,:))
                END DO
                vitloc(TYPE,l) = SUM(vit(j_loc,TYPE)*mesh%gauss%ww(:,l))
             END DO

             !--------Calcul de la divergence
             div1 = ABS(grad(1,1,l) + vitloc(1,l)/ray + grad(2,5,l) + vitloc(4,l)*mode/ray)
             div2 = ABS(grad(1,2,l) + vitloc(2,l)/ray + grad(2,6,l) - vitloc(3,l)*mode/ray)
             visc2(1) = MAX(visc2(1),div1)
             visc2(2) = MAX(visc2(2),div2)
          END DO
          visc2(4) = visc2(1); visc2(5) = visc2(1)
          visc2(3) = visc2(2); visc2(6) = visc2(2)

          nu_loc = 0.d0
          DO TYPE = 1, 6
             !visc1(type) = MAX(visc1(type)/normal_vit,visc2(type))
             !visc1(type) = MAX(visc1(type),visc2(type)*normal_vit)/MAXVAL(ABS(vitloc(type,:)))
             visc1(TYPE) = MAX(visc1(TYPE),2*visc2(TYPE)*normal_vit)/norm_vit(TYPE)
             visc_eff = mesh%hloc(m)*MIN(coeff_visc_ordre_un*normal_vit,inputs%LES_coeff2*mesh%hloc(m)*visc1(TYPE))
             nu_loc = nu_loc + visc_eff
             !======Semi-implicit version==========
             viscosity(TYPE,m) = inputs%LES_coeff1*mesh%hloc(m) - visc_eff
             !======Semi-implicit version==========
          END DO
          R_eff = R_eff + (nu_loc + dt**2*mesh%hloc(m))*mesh%hloc(m)**2*ray
          visc_plot(m) = (nu_loc/6)/(coeff_visc_ordre_un*mesh%hloc(m)*normal_vit)
       END DO
    ELSE
       cfl = 0.d0
       viscosity = 0.d0
       DO m = 1, mesh%dom_me
          j_loc = mesh%jj(:,m)
          vloc = maxval(vel_tot(j_loc))
          cfl = max(vloc*dt/mesh%hloc(m),cfl)
       END DO
    END IF


    !DO type = 1, 6
    !   CALL average(mesh,viscosity(type,:))
    !END DO

    nw = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       mult(:)= viscosity(:,m)
       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_3_LA%loc_to_glob(3,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)
          DO TYPE = 1, 6
             DO k = 1 ,2
                grad(k,TYPE,l) = SUM(vit(j_loc,TYPE)*dw_loc(k,:))
             END DO
             vitloc(TYPE,l) = SUM(vit(j_loc,TYPE)*mesh%gauss%ww(:,l))
          END DO

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fv(1) = ((mode*vitloc(1,l)+vitloc(4,l))*mode +mode*vitloc(4,l)+vitloc(1,l))/ray**2
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2) = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:))
          fv(2) = ((mode*vitloc(2,l)-vitloc(3,l))*mode -mode*vitloc(3,l)+vitloc(2,l))/ray**2
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3) = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3) = -SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          fv(3) = (-mode*vitloc(2,l)+vitloc(3,l) +(mode*vitloc(3,l)-vitloc(2,l))*mode)/ray**2
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4) = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4) = SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          fv(4) =  (mode*vitloc(1,l)+vitloc(4,l) +(mode*vitloc(4,l)+vitloc(1,l))*mode)/ray**2
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5) = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:))
          fv(5) =  vitloc(5,l)*(mode/ray)**2
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6) = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:))
          fv(6) = vitloc(6,l)*(mode/ray)**2
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          fv = mult*fv

          rhs_gauss(index, :) =  (ft+fp+fs+fv-rotv_v(index,:))

       ENDDO

    ENDDO

    CALL rhs_3x3(mesh, vv_3_LA, mode, rhs_gauss, vb_145, vb_236)

  END SUBROUTINE qs_ns_stab_3x3

  SUBROUTINE qs_01_div_hybrid_generic (type_fe_velocity,vv_mesh, pp_mesh, pp_LA, mode, gg, pb_1, pb_2)
    !================================================
    USE def_type_mesh
    USE basis_change
    IMPLICIT NONE
    TYPE(mesh_type)                             :: vv_mesh, pp_mesh
    TYPE(petsc_csr_LA)                          :: pp_LA
    INTEGER,                        INTENT(IN)  :: type_fe_velocity
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg
    INTEGER     ,                   INTENT(IN)  :: mode
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, vv_mesh%gauss%l_G) ::  w_c
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w)  :: v1_loc, v2_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_w)  :: idxm
    INTEGER :: m, l, n, i, nw, nwc, ni, iglob
    REAL(KIND=8), DIMENSION(3,2) :: f
    REAL(KIND=8)   :: ray
    INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
    INTEGER, DIMENSION(pp_mesh%gauss%n_w) :: jc_loc
    !===JLG July 20th, 2019
    !===Change of basis arrays in save status
    REAL(KIND=8), DIMENSION(:,:), POINTER, SAVE :: aij_p1p2, aij_p2p3
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
    LOGICAL :: once_p1p2=.TRUE., once_p2p3=.TRUE.
    Vec            :: pb_1, pb_2
    PetscErrorCode :: ierr

    !===JLG July 20th, 2019
    !===Construct aij_p2, aij_p3
    IF (type_fe_velocity==2) THEN
       IF (once_p1p2) THEN
          CALL p1_p2(aij_p1p2)
          once_p1p2=.FALSE.
       END IF
       aij => aij_p1p2
    ELSE IF (type_fe_velocity==3) THEN
       IF (once_p2p3) THEN
          CALL p2_p3(aij_p2p3)
          once_p2p3=.FALSE.
       END IF
       aij => aij_p2p3
    ELSE
       CALL error_petsc('qs_01_div_hybrid_generic, type_fe_velocity not correct')
    END IF
    !===JLG July 20th, 2019

    CALL VecSet(pb_1, 0.d0, ierr)
    CALL VecSet(pb_2, 0.d0, ierr)

    !===JLG July 20th, 2019
    !===Construct w_c
    nw = vv_mesh%gauss%n_w
    nwc = pp_mesh%gauss%n_w
    DO l = 1, vv_mesh%gauss%l_G
       DO n = 1, nwc
          w_c(n,l) = SUM(aij(n,:)*vv_mesh%gauss%ww(:,l))
       END DO
    END DO
    !===JLG July 20th, 2019

    DO m = 1, vv_mesh%dom_me
       j_loc(:) = vv_mesh%jj(:,m)
       jc_loc(:)= pp_mesh%jj(:,m)

       DO ni = 1, nwc
          i = jc_loc(ni)
          iglob = pp_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v1_loc = 0.d0
       v2_loc = 0.d0
       DO l = 1, vv_mesh%gauss%l_G

          !===radius of velocity Gauss point
          ray = 0
          DO n = 1, nw;  i = vv_mesh%jj(n,m)
             ray = ray + vv_mesh%rr(1,i)*vv_mesh%gauss%ww(n,l)
          END DO
          !----------------------

          !===Compute divergence on cosines
          f(1,1) = (ray*SUM(gg(j_loc,1)*vv_mesh%gauss%dw(1,:,l,m)) &
               + SUM(gg(j_loc,1)*vv_mesh%gauss%ww(:,l)))
          f(2,1) = mode*SUM(gg(j_loc,4)*vv_mesh%gauss%ww(:,l))
          f(3,1) =      SUM(gg(j_loc,5)*vv_mesh%gauss%dw(2,:,l,m)) * ray

          !===Compute divergence on sines
          f(1,2)  = (ray*SUM(gg(j_loc,2)*vv_mesh%gauss%dw(1,:,l,m)) &
               + SUM(gg(j_loc,2)*vv_mesh%gauss%ww(:,l)))
          f(2,2)  =-mode*SUM(gg(j_loc,3)*vv_mesh%gauss%ww(:,l))
          f(3,2)  =      SUM(gg(j_loc,6)*vv_mesh%gauss%dw(2,:,l,m)) * ray

          f = f *vv_mesh%gauss%rj(l,m)

          x = f(1,1)+f(2,1)+f(3,1)
          DO ni = 1, nwc
             v1_loc(ni) = v1_loc(ni) + w_c(ni,l)*x
          END DO

          x = f(1,2)+f(2,2)+f(3,2)
          DO ni = 1, nwc
             v2_loc(ni) = v2_loc(ni) + w_c(ni,l)*x
          END DO

       ENDDO
       CALL VecSetValues(pb_1, nwc, idxm, v1_loc, ADD_VALUES, ierr)
       CALL VecSetValues(pb_2, nwc, idxm, v2_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(pb_1,ierr)
    CALL VecAssemblyEnd(pb_1,ierr)

    CALL VecAssemblyBegin(pb_2,ierr)
    CALL VecAssemblyEnd(pb_2,ierr)

  END SUBROUTINE qs_01_div_hybrid_generic

  SUBROUTINE inject_generic(type_fe_velocity, jj_c, jj_f, pp_c, pp_f)
    USE basis_change
    IMPLICIT NONE
    INTEGER,                      INTENT(IN)  :: type_fe_velocity
    INTEGER,      DIMENSION(:,:), INTENT(IN)  :: jj_c, jj_f
    REAL(KIND=8), DIMENSION(:),   INTENT(IN)  :: pp_c
    REAL(KIND=8), DIMENSION(:),   INTENT(OUT) :: pp_f
    INTEGER:: m, n
    !===Change of basis arrays in save status
    REAL(KIND=8), DIMENSION(:,:), POINTER, SAVE :: aij_p1p2, aij_p2p3
    REAL(KIND=8), DIMENSION(:,:), POINTER :: aij
    LOGICAL :: once_p1p2=.TRUE., once_p2p3=.TRUE.

    !===JLG July 20th, 2019
    !===Construct aij_p2, aij_p3
    IF (type_fe_velocity==2) THEN
       IF (once_p1p2) THEN
          CALL p1_p2(aij_p1p2)
          once_p1p2=.FALSE.
       END IF
       aij => aij_p1p2
    ELSE IF (type_fe_velocity==3) THEN
       IF (once_p2p3) THEN
          CALL p2_p3(aij_p2p3)
          once_p2p3=.FALSE.
       END IF
       aij => aij_p2p3
    ELSE
       CALL error_petsc('inject_generic, type_fe_velocity not correct')
    END IF
    !===JLG July 20th, 2019

    DO m = 1, size(jj_f,2)
       DO n = 1, size(jj_f,1)
          pp_f(jj_f(n,m)) = SUM(aij(:,n)*pp_c(jj_c(:,m)))
       END DO
    END DO

  END SUBROUTINE inject_generic

  SUBROUTINE qs_01_div_hybrid_2006 (vv_mesh, pp_mesh, pp_LA, mode, gg, pb_1, pb_2)
    !================================================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                             :: vv_mesh, pp_mesh
    TYPE(petsc_csr_LA)                          :: pp_LA
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN)  :: gg
    INTEGER     ,                   INTENT(IN)  :: mode
    REAL(KIND=8) :: x
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w, vv_mesh%gauss%l_G) ::  w_c
    REAL(KIND=8), DIMENSION(pp_mesh%gauss%n_w)  :: v1_loc, v2_loc
    INTEGER,      DIMENSION(pp_mesh%gauss%n_w)  :: idxm
    INTEGER :: m, l, n, i, nw, nwc, ni, iglob
    REAL(KIND=8), DIMENSION(3,2) :: f
    REAL(KIND=8)   :: ray
    INTEGER, DIMENSION(vv_mesh%gauss%n_w) :: j_loc
    INTEGER, DIMENSION(pp_mesh%gauss%n_w) :: jc_loc
    !#include "petsc/finclude/petsc.h"
    Vec            :: pb_1, pb_2
    PetscErrorCode :: ierr

    !===Test for correctness
    IF (inputs%type_fe_velocity.NE.2) THEN
       call error_petsc('BUG qs_01_div_hybrid_2006, inputs%type_fe_velocity.NE.2')
    END IF

    CALL VecSet(pb_1, 0.d0, ierr)
    CALL VecSet(pb_2, 0.d0, ierr)

    nw = vv_mesh%gauss%n_w
    DO l = 1, vv_mesh%gauss%l_G
       !P1/P2
       w_c(1,l) = vv_mesh%gauss%ww(1,l) + 0.5*(vv_mesh%gauss%ww(nw-1,l) + vv_mesh%gauss%ww(nw,l))
       w_c(2,l) = vv_mesh%gauss%ww(2,l) + 0.5*(vv_mesh%gauss%ww(nw,l) + vv_mesh%gauss%ww(4,l))
       w_c(3,l) = vv_mesh%gauss%ww(3,l) + 0.5*(vv_mesh%gauss%ww(4,l) + vv_mesh%gauss%ww(nw-1,l))
    END DO

    nwc = pp_mesh%gauss%n_w
    DO m = 1, vv_mesh%dom_me
       j_loc(:) = vv_mesh%jj(:,m)
       jc_loc(:)= pp_mesh%jj(:,m)

       DO ni = 1, nwc
          i = jc_loc(ni)
          iglob = pp_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v1_loc = 0.d0
       v2_loc = 0.d0
       DO l = 1, vv_mesh%gauss%l_G

          !--------radius of P2 Gauss point
          ray = 0
          DO n = 1, nw;  i = vv_mesh%jj(n,m)
             ray = ray + vv_mesh%rr(1,i)*vv_mesh%gauss%ww(n,l)
          END DO
          !----------------------

          !calcul de la divergence sur les cos
          f(1,1) = (ray*SUM(gg(j_loc,1)*vv_mesh%gauss%dw(1,:,l,m)) &
               + SUM(gg(j_loc,1)*vv_mesh%gauss%ww(:,l)))
          f(2,1) = mode*SUM(gg(j_loc,4)*vv_mesh%gauss%ww(:,l))
          f(3,1) =      SUM(gg(j_loc,5)*vv_mesh%gauss%dw(2,:,l,m)) * ray

          !calcul de la divergence sur les sin
          f(1,2)  = (ray*SUM(gg(j_loc,2)*vv_mesh%gauss%dw(1,:,l,m)) &
               + SUM(gg(j_loc,2)*vv_mesh%gauss%ww(:,l)))
          f(2,2)  =-mode*SUM(gg(j_loc,3)*vv_mesh%gauss%ww(:,l))
          f(3,2)  =      SUM(gg(j_loc,6)*vv_mesh%gauss%dw(2,:,l,m)) * ray

          f = f *vv_mesh%gauss%rj(l,m)

!!$          DO j=1,2
!!$             x = f(1,j)+f(2,j)+f(3,j)
!!$             DO n=1, nwc
!!$                u0_c(jj_c(n,m),j) = u0_c(jj_c(n,m),j) + w_c(n,l)*x
!!$             ENDDO
!!$          ENDDO

          x = f(1,1)+f(2,1)+f(3,1)
          DO ni = 1, nwc
             v1_loc(ni) = v1_loc(ni) + w_c(ni,l)*x
          END DO

          x = f(1,2)+f(2,2)+f(3,2)
          DO ni = 1, nwc
             v2_loc(ni) = v2_loc(ni) + w_c(ni,l)*x
          END DO

       ENDDO
       CALL VecSetValues(pb_1, nwc, idxm, v1_loc, ADD_VALUES, ierr)
       CALL VecSetValues(pb_2, nwc, idxm, v2_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(pb_1,ierr)
    CALL VecAssemblyEnd(pb_1,ierr)

    CALL VecAssemblyBegin(pb_2,ierr)
    CALL VecAssemblyEnd(pb_2,ierr)

  END SUBROUTINE qs_01_div_hybrid_2006

  SUBROUTINE qs_00 (mesh, LA, ff, vect)
    !=================================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    type(petsc_csr_LA)                          :: LA
    REAL(KIND=8), DIMENSION(:), INTENT(IN)      :: ff
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: ff_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: jj_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: v_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: idxm
    INTEGER ::  i, m, l, ni, iglob
    REAL(KIND=8) :: fl, ray
    !#include "petsc/finclude/petsc.h"
    Vec                                         :: vect
    PetscErrorCode                              :: ierr

    CALL VecSet(vect, 0.d0, ierr)

    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       ff_loc = ff(jj_loc)
       DO ni = 1, mesh%gauss%n_w
          i = mesh%jj(ni,m)
          iglob = LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          ray = 0
          DO ni = 1, mesh%gauss%n_w
             i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          fl = SUM(ff_loc*mesh%gauss%ww(:,l))*mesh%gauss%rj(l,m)*ray
          DO ni = 1,  mesh%gauss%n_w
             v_loc(ni) = v_loc(ni) +  mesh%gauss%ww(ni,l) * fl
          END DO
       ENDDO
       CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(vect,ierr)
    CALL VecAssemblyEnd(vect,ierr)
  END SUBROUTINE qs_00

  SUBROUTINE qs_00_gauss (mesh, LA, heat_capa, ff, ff_gauss, vect)
    !=================================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    type(petsc_csr_LA)                          :: LA
    REAL(KIND=8), DIMENSION(:), INTENT(IN)      :: heat_capa
    REAL(KIND=8), DIMENSION(:), INTENT(IN)      :: ff, ff_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: ff_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: jj_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: v_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: idxm
    INTEGER ::  i, m, l, ni, iglob, index
    REAL(KIND=8) :: fl, ray
    !#include "petsc/finclude/petsc.h"
    Vec                                         :: vect
    PetscErrorCode                              :: ierr

    CALL VecSet(vect, 0.d0, ierr)

    index = 0
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       ff_loc = ff(jj_loc)
       DO ni = 1, mesh%gauss%n_w
          i = mesh%jj(ni,m)
          iglob = LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  =index + 1
          ray = 0
          DO ni = 1, mesh%gauss%n_w
             i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          fl = (heat_capa(m) * SUM(ff_loc*mesh%gauss%ww(:,l)) + ff_gauss(index))*mesh%gauss%rj(l,m)*ray
          DO ni = 1,  mesh%gauss%n_w
             v_loc(ni) = v_loc(ni) +  mesh%gauss%ww(ni,l) * fl
          END DO
       ENDDO
       CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(vect,ierr)
    CALL VecAssemblyEnd(vect,ierr)
  END SUBROUTINE qs_00_gauss

  SUBROUTINE qs_00_gauss_conc (mesh, LA, ff, ff_gauss, vect)
    !=================================
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type), TARGET                     :: mesh
    type(petsc_csr_LA)                          :: LA
    REAL(KIND=8), DIMENSION(:), INTENT(IN)      :: ff, ff_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: ff_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: jj_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)     :: v_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)     :: idxm
    INTEGER ::  i, m, l, ni, iglob, index
    REAL(KIND=8) :: fl, ray
    !#include "petsc/finclude/petsc.h"
    Vec                                         :: vect
    PetscErrorCode                              :: ierr

    CALL VecSet(vect, 0.d0, ierr)

    index = 0
    DO m = 1, mesh%dom_me
       jj_loc = mesh%jj(:,m)
       ff_loc = ff(jj_loc)
       DO ni = 1, mesh%gauss%n_w
          i = mesh%jj(ni,m)
          iglob = LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO

       v_loc = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  =index + 1
          ray = 0
          DO ni = 1, mesh%gauss%n_w
             i = jj_loc(ni)
             ray = ray + mesh%rr(1,i)*mesh%gauss%ww(ni,l)
          END DO
          fl = ( SUM(ff_loc*mesh%gauss%ww(:,l)) + ff_gauss(index))*mesh%gauss%rj(l,m)*ray
          DO ni = 1,  mesh%gauss%n_w
             v_loc(ni) = v_loc(ni) +  mesh%gauss%ww(ni,l) * fl
          END DO
       ENDDO
       CALL VecSetValues(vect, mesh%gauss%n_w, idxm, v_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(vect,ierr)
    CALL VecAssemblyEnd(vect,ierr)
  END SUBROUTINE qs_00_gauss_conc

  SUBROUTINE qs_ns_momentum_stab_new(mesh,vv_1_LA,vv_2_LA,mode,ff,V1m,P,&
       vb_2_14,vb_2_23,vb_1_5,vb_1_6, rotb_b, tensor, tensor_surface_gauss,&
       stab, momentum, momentum_LES, visc_grad_vel, visc_entro)
    !=================================
    !RHS for Navier-Stokes with momentum
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    !USE sub_plot
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_1_LA,vv_2_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotb_b
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor               !(node)
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor_surface_gauss !(gauss)
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, momentum, momentum_LES   !V(node, type)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8),                               INTENT(IN) :: stab
    INTEGER,                                    INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: visc_grad_vel        !V(r/th/z, gauss, type)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: visc_entro
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, mult
    REAL(KIND=8), DIMENSION(6)                             :: fnl, fvgm, fvgm_LES, fvgu
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad_mom
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: momloc, momloc_LES
    REAL(KIND=8), DIMENSION(3,6,mesh%gauss%l_G)            :: tensor_loc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8)   :: ray
    REAL(KIND=8), DIMENSION(mesh%dom_me)                   :: stab_loc
    LOGICAL,                                 SAVE :: once = .TRUE.
    REAL(KIND=8),                            SAVE :: type_fe
    INTEGER :: m, l , i , ni , index, TYPE, k
    INTEGER :: nw, ix, ki, iglob
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_2_14,vb_2_23,vb_1_5,vb_1_6
    PetscErrorCode :: ierr

    CALL VecSet(vb_2_14, 0.d0, ierr)
    CALL VecSet(vb_2_23, 0.d0, ierr)
    CALL VecSet(vb_1_5, 0.d0, ierr)
    CALL VecSet(vb_1_6, 0.d0, ierr)

    IF (once) THEN
       once =.FALSE.

       IF (.NOT.inputs%LES) THEN
          inputs%LES_coeff1=0.d0
       END IF

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE
          type_fe = 2
       END IF
    END IF ! end once

    IF (inputs%LES) THEN
       DO m = 1, mesh%dom_me
          stab_loc(m) = stab + inputs%LES_coeff1*mesh%hloc(m)
       END DO
    ELSE
       stab_loc(:) = stab
    END IF

    index = 0
    nw = mesh%gauss%n_w

    DO m = 1, mesh%dom_me
       mult(:)=  - visc_entro(m,1) !visc_entro is the same in type 1 or 2

       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_1_LA%loc_to_glob(1,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_2_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)
          DO TYPE = 1, 6
             DO k = 1, 3
                tensor_loc(k,TYPE,l) = SUM(tensor(k,j_loc,TYPE)*mesh%gauss%ww(:,l)) &
                     + tensor_surface_gauss(k,index,TYPE)
             END DO
          END DO

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1)  = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1)  = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1)  = -SUM(P(j_loc,1)*dw_loc(1,:))
          fnl(1) = (-mode*tensor_loc(1,4,l) + tensor_loc(2,3,l))/ray
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2)  = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2)  = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2)  = -SUM(P(j_loc,2)*dw_loc(1,:))
          fnl(2) = (mode*tensor_loc(1,3,l) + tensor_loc(2,4,l))/ray
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3)  = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3)  = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3)  = -SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          fnl(3) = (-tensor_loc(1,3,l) - mode*tensor_loc(2,4,l))/ray
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4)  = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4)  = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4)  = SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          fnl(4) = (-tensor_loc(1,4,l) + mode*tensor_loc(2,3,l))/ray
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5)  = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5)  = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5)  = -SUM(P(j_loc,1)*dw_loc(2,:))
          fnl(5) = -tensor_loc(3,4,l)*mode/ray
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6)  = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6)  = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6)  = -SUM(P(j_loc,2)*dw_loc(2,:))
          fnl(6) = tensor_loc(3,3,l)*mode/ray
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          IF (inputs%if_level_set) THEN
             DO TYPE = 1, 6

                DO k = 1 ,2
                   grad_mom(k,TYPE,l) = stab_loc(m)* SUM(momentum(j_loc,TYPE)*dw_loc(k,:)) &
                        + mult(TYPE)*SUM(momentum_LES(j_loc,TYPE)*dw_loc(k,:))
                END DO

                momloc(TYPE,l) = SUM(momentum(j_loc,TYPE)*mesh%gauss%ww(:,l))
                momloc_LES(TYPE,l) = SUM(momentum_LES(j_loc,TYPE)*mesh%gauss%ww(:,l))

             END DO
             fvgm(1) = ((mode*momloc(1,l)+momloc(4,l))*mode + mode*momloc(4,l)+momloc(1,l))/ray**2
             fvgm(2) = ((mode*momloc(2,l)-momloc(3,l))*mode - mode*momloc(3,l)+momloc(2,l))/ray**2
             fvgm(3) = (-mode*momloc(2,l)+momloc(3,l) + (mode*momloc(3,l)-momloc(2,l))*mode)/ray**2
             fvgm(4) = (mode*momloc(1,l)+momloc(4,l) + (mode*momloc(4,l)+momloc(1,l))*mode)/ray**2
             fvgm(5) = momloc(5,l)*(mode/ray)**2
             fvgm(6) = momloc(6,l)*(mode/ray)**2

             fvgm_LES(1) = ((mode*momloc_LES(1,l)+momloc_LES(4,l))*mode + mode*momloc_LES(4,l)+momloc_LES(1,l))/ray**2
             fvgm_LES(2) = ((mode*momloc_LES(2,l)-momloc_LES(3,l))*mode - mode*momloc_LES(3,l)+momloc_LES(2,l))/ray**2
             fvgm_LES(3) = (-mode*momloc_LES(2,l)+momloc_LES(3,l) + (mode*momloc_LES(3,l)-momloc_LES(2,l))*mode)/ray**2
             fvgm_LES(4) = (mode*momloc_LES(1,l)+momloc_LES(4,l) + (mode*momloc_LES(4,l)+momloc_LES(1,l))*mode)/ray**2
             fvgm_LES(5) = momloc_LES(5,l)*(mode/ray)**2
             fvgm_LES(6) = momloc_LES(6,l)*(mode/ray)**2

             fvgm = stab_loc(m)*fvgm + mult*fvgm_LES

             fvgu(1) = (-visc_grad_vel(1,index,4)*mode + visc_grad_vel(2,index,3))/ray
             fvgu(2) = (visc_grad_vel(1,index,3)*mode + visc_grad_vel(2,index,4))/ray
             fvgu(3) = (-visc_grad_vel(1,index,3) - visc_grad_vel(2,index,4)*mode)/ray
             fvgu(4) = (-visc_grad_vel(1,index,4) + visc_grad_vel(2,index,3)*mode)/ray
             fvgu(5) = -visc_grad_vel(3,index,4)*mode/ray
             fvgu(6) = visc_grad_vel(3,index,3)*mode/ray

             DO TYPE = 1, 6
                grad_mom(:,TYPE,l) = grad_mom(:,TYPE,l)*ray*mesh%gauss%rj(l,m)
             END DO
             tensor_loc(:,:,l) = tensor_loc(:,:,l) + visc_grad_vel(:,index,:)

          ELSE
             grad_mom = 0.d0
             fvgm     = 0.d0
             fvgu     = 0.d0
          END IF

          ! if NOT level_set then :
          ! rhs = (BDF2(n&n_m1 terms) - Grad_pressure + source_in_ns + precession + lorentz)*test_function
          !      + tensor_explicit:Grad(test_function)   (tensor by FFT)

          ! if level_set then :
          !rhs = (BDF2(n&n_m1 terms) - Grad_pressure + source_in_ns + precession + lorentz)*test_function
          !      + (tensor_explicit + visc_grad_vel):Grad(test_function)   (tensor and visc_grad_vel by FFT)
          !      + (LES + stab*grad_mom):GRAD(test_function)

          smb =  (ft+fp+fs+rotb_b(index,:)+fnl+fvgm+fvgu)*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(1) + SUM(dw_loc(:,ni)*grad_mom(:,1,l)) &
                  + (dw_loc(1,ni)*tensor_loc(1,1,l) + dw_loc(2,ni)*tensor_loc(1,5,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(2) + SUM(dw_loc(:,ni)*grad_mom(:,2,l)) &
                  + (dw_loc(1,ni)*tensor_loc(1,2,l) + dw_loc(2,ni)*tensor_loc(1,6,l))*ray*mesh%gauss%rj(l,m)
             v5_loc(ix)  = v5_loc(ix)  + mesh%gauss%ww(ni,l)*smb(5) + SUM(dw_loc(:,ni)*grad_mom(:,5,l)) &
                  + (dw_loc(1,ni)*tensor_loc(3,1,l) + dw_loc(2,ni)*tensor_loc(3,5,l))*ray*mesh%gauss%rj(l,m)
             v6_loc(ix)  = v6_loc(ix)  + mesh%gauss%ww(ni,l)*smb(6) + SUM(dw_loc(:,ni)*grad_mom(:,6,l)) &
                  + (dw_loc(1,ni)*tensor_loc(3,2,l) + dw_loc(2,ni)*tensor_loc(3,6,l))*ray*mesh%gauss%rj(l,m)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(4) + SUM(dw_loc(:,ni)*grad_mom(:,4,l)) &
                  + (dw_loc(1,ni)*tensor_loc(2,2,l) + dw_loc(2,ni)*tensor_loc(2,6,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(3) + SUM(dw_loc(:,ni)*grad_mom(:,3,l)) &
                  + (dw_loc(1,ni)*tensor_loc(2,1,l) + dw_loc(2,ni)*tensor_loc(2,5,l))*ray*mesh%gauss%rj(l,m)
          END DO
       ENDDO

       CALL VecSetValues(vb_2_14, 2*nw, idxm_2, v14_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_2_23, 2*nw, idxm_2, v23_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_5,    nw, idxm_1,  v5_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_6,    nw, idxm_1,  v6_loc, ADD_VALUES, ierr)
    ENDDO

    IF (inputs%LES) THEN
       IF (mesh%edge_stab) THEN
          CALL error_Petsc('BUG in qs_ns_stab_new: terms with edge_stab not yet assembled')
       END IF
    END IF

    CALL VecAssemblyBegin(vb_2_14,ierr)
    CALL VecAssemblyEnd(vb_2_14,ierr)
    CALL VecAssemblyBegin(vb_2_23,ierr)
    CALL VecAssemblyEnd(vb_2_23,ierr)
    CALL VecAssemblyBegin(vb_1_5,ierr)
    CALL VecAssemblyEnd(vb_1_5,ierr)
    CALL VecAssemblyBegin(vb_1_6,ierr)
    CALL VecAssemblyEnd(vb_1_6,ierr)

  END SUBROUTINE qs_ns_momentum_stab_new

  SUBROUTINE qs_ns_momentum_stab_3x3(mesh,vv_3_LA,mode,ff,V1m,P,&
       vb_3_145,vb_3_236, rotb_b, tensor, tensor_surface_gauss,&
       stab, momentum, visc_grad_vel)
    !=================================
    !RHS for Navier-Stokes
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    !USE sub_plot
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_3_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff, rotb_b
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor               !(node)
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor_surface_gauss !(gauss)
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m, momentum    !V(noeud, type)
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8),                               INTENT(IN) :: stab
    INTEGER,                                    INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: visc_grad_vel        !V(r/th/z, gauss, type)
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, smb, fnl
    REAL(KIND=8), DIMENSION(6)                             :: fvgm, fvgu, fvgmT
    REAL(KIND=8), DIMENSION(2,6,mesh%gauss%l_G)            :: grad_mom, grad_T_mom
    REAL(KIND=8), DIMENSION(6,mesh%gauss%l_G)              :: momloc
    REAL(KIND=8), DIMENSION(3,6,mesh%gauss%l_G)            :: tensor_loc
    INTEGER, DIMENSION(mesh%gauss%n_w)                     :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8)   :: ray
    REAL(KIND=8), DIMENSION(mesh%dom_me)                   :: stab_loc
    LOGICAL,                                 SAVE :: once = .TRUE.
    REAL(KIND=8),                            SAVE :: type_fe
    INTEGER :: m, l , i , ni , index, TYPE, k
    INTEGER :: nw, ix, ki, iglob
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_3_145, vb_3_236
    PetscErrorCode :: ierr

    CALL VecSet(vb_3_145, 0.d0, ierr)
    CALL VecSet(vb_3_236, 0.d0, ierr)

    IF (once) THEN
       once =.FALSE.

       IF (.NOT.inputs%LES) THEN
          inputs%LES_coeff1=0.d0
       END IF

       IF (mesh%gauss%n_w==3) THEN
          type_fe = 1
       ELSE
          type_fe = 2
       END IF
    END IF !end once

    stab_loc(:) = stab

    nw = mesh%gauss%n_w
    index = 0

    DO m = 1, mesh%dom_me

       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_3_LA%loc_to_glob(3,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)
          DO TYPE = 1, 6
             DO k = 1, 3
                tensor_loc(k,TYPE,l) = SUM(tensor(k,j_loc,TYPE)*mesh%gauss%ww(:,l)) &
                     + tensor_surface_gauss(k,index,TYPE)
             END DO
          END DO

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !-calcul des seconds membres pr les termes de forçage, temporels et de gradP
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1) = -SUM(P(j_loc,1)*dw_loc(1,:))
          fnl(1) = (-mode*tensor_loc(1,4,l) + tensor_loc(2,3,l))/ray
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2) = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2) = -SUM(P(j_loc,2)*dw_loc(1,:))
          fnl(2) = (mode*tensor_loc(1,3,l) + tensor_loc(2,4,l))/ray
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3) = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3) = -SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          fnl(3) = (-tensor_loc(1,3,l) - mode*tensor_loc(2,4,l))/ray
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4) = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4) = SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          fnl(4) = (-tensor_loc(1,4,l) + mode*tensor_loc(2,3,l))/ray
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5) = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5) = -SUM(P(j_loc,1)*dw_loc(2,:))
          fnl(5) = -tensor_loc(3,4,l)*mode/ray
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6) = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6) = -SUM(P(j_loc,2)*dw_loc(2,:))
          fnl(6) = tensor_loc(3,3,l)*mode/ray
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          IF (inputs%if_level_set) THEN
             DO TYPE = 1, 6
                DO k = 1 ,2
                   grad_mom(k,TYPE,l) = stab_loc(m)*SUM(momentum(j_loc,TYPE)*dw_loc(k,:))
                END DO

                momloc(TYPE,l) = SUM(momentum(j_loc,TYPE)*mesh%gauss%ww(:,l))
             END DO

             fvgm(1) = ((mode*momloc(1,l)+momloc(4,l))*mode + mode*momloc(4,l)+momloc(1,l))/ray**2
             fvgm(2) = ((mode*momloc(2,l)-momloc(3,l))*mode - mode*momloc(3,l)+momloc(2,l))/ray**2
             fvgm(3) = (-mode*momloc(2,l)+momloc(3,l) + (mode*momloc(3,l)-momloc(2,l))*mode)/ray**2
             fvgm(4) = (mode*momloc(1,l)+momloc(4,l) + (mode*momloc(4,l)+momloc(1,l))*mode)/ray**2
             fvgm(5) = momloc(5,l)*(mode/ray)**2
             fvgm(6) = momloc(6,l)*(mode/ray)**2

             !===Compute Grad_T_mom . r or z derivative of test function
             grad_T_mom(1,1,l) = SUM(momentum(j_loc,1)*dw_loc(1,:))
             grad_T_mom(2,1,l) = SUM(momentum(j_loc,5)*dw_loc(1,:))
             grad_T_mom(1,2,l) = SUM(momentum(j_loc,2)*dw_loc(1,:))
             grad_T_mom(2,2,l) = SUM(momentum(j_loc,6)*dw_loc(1,:))
             grad_T_mom(1,3,l) = (mode*momloc(2,l) - momloc(3,l))/ray
             grad_T_mom(2,3,l) = mode*momloc(6,l)/ray
             grad_T_mom(1,4,l) = (-mode*momloc(1,l) - momloc(4,l))/ray
             grad_T_mom(2,4,l) = -mode*momloc(5,l)/ray
             grad_T_mom(1,5,l) = SUM(momentum(j_loc,1)*dw_loc(2,:))
             grad_T_mom(2,5,l) = SUM(momentum(j_loc,5)*dw_loc(2,:))
             grad_T_mom(1,6,l) = SUM(momentum(j_loc,2)*dw_loc(2,:))
             grad_T_mom(2,6,l) = SUM(momentum(j_loc,6)*dw_loc(2,:))

             !===Compute Grad_T_mom . NO (r or z derivative) of test function
             fvgmT(1) = -mode*SUM(momentum(j_loc,4)*dw_loc(1,:))/ray &
                  + (mode*momloc(4,l)+momloc(1,l))/ray**2
             fvgmT(2) = mode*SUM(momentum(j_loc,3)*dw_loc(1,:))/ray &
                  + (-mode*momloc(3,l)+momloc(2,l))/ray**2
             fvgmT(3) = -SUM(momentum(j_loc,3)*dw_loc(1,:))/ray &
                  + mode*(mode*momloc(3,l)-momloc(2,l))/ray**2
             fvgmT(4) = -SUM(momentum(j_loc,4)*dw_loc(1,:))/ray &
                  + mode*(mode*momloc(4,l)+momloc(1,l))/ray**2
             fvgmT(5) = -mode*SUM(momentum(j_loc,4)*dw_loc(2,:))/ray
             fvgmT(6) = mode*SUM(momentum(j_loc,3)*dw_loc(2,:))/ray

             !===update grad_mom
             grad_mom(:,:,l) = grad_mom(:,:,l) + stab*grad_T_mom(:,:,l)
             fvgm     = stab_loc(m)*fvgm + stab*fvgmT

             fvgu(1) = (-visc_grad_vel(1,index,4)*mode + visc_grad_vel(2,index,3))/ray
             fvgu(2) = (visc_grad_vel(1,index,3)*mode + visc_grad_vel(2,index,4))/ray
             fvgu(3) = (-visc_grad_vel(1,index,3) - visc_grad_vel(2,index,4)*mode)/ray
             fvgu(4) = (-visc_grad_vel(1,index,4) + visc_grad_vel(2,index,3)*mode)/ray
             fvgu(5) = -visc_grad_vel(3,index,4)*mode/ray
             fvgu(6) = visc_grad_vel(3,index,3)*mode/ray

             DO TYPE = 1, 6
                grad_mom(:,TYPE,l) = grad_mom(:,TYPE,l)*ray*mesh%gauss%rj(l,m)
             END DO
             tensor_loc(:,:,l) = tensor_loc(:,:,l) + visc_grad_vel(:,index,:)

          ELSE
             grad_mom = 0.d0
             fvgm     = 0.d0
             fvgu     = 0.d0
          END IF

          ! if NOT level_set then :
          ! rhs = (BDF2(n&n_m1 terms) - Grad_pressure + source_in_ns + precession + lorentz)*test_function
          !      + tensor_explicit:Grad(test_function)   (tensor by FFT)

          ! if level_set then :
          !rhs = (BDF2(n&n_m1 terms) - Grad_pressure + source_in_ns + precession + lorentz)*test_function
          !      + (tensor_explicit + visc_grad_vel):Grad(test_function)   (tensor and visc_grad_vel by FFT)
          !      + (LES + stab*grad_mom):GRAD(test_function)

          smb =  (ft+fp+fs+rotb_b(index,:)+fnl+fvgm+fvgu)*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(1) + SUM(dw_loc(:,ni)*grad_mom(:,1,l)) &
                  + (dw_loc(1,ni)*tensor_loc(1,1,l) + dw_loc(2,ni)*tensor_loc(1,5,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(2) + SUM(dw_loc(:,ni)*grad_mom(:,2,l)) &
                  + (dw_loc(1,ni)*tensor_loc(1,2,l) + dw_loc(2,ni)*tensor_loc(1,6,l))*ray*mesh%gauss%rj(l,m)
             v5_loc(ix)  = v5_loc(ix)  + mesh%gauss%ww(ni,l)*smb(5) + SUM(dw_loc(:,ni)*grad_mom(:,5,l)) &
                  + (dw_loc(1,ni)*tensor_loc(3,1,l) + dw_loc(2,ni)*tensor_loc(3,5,l))*ray*mesh%gauss%rj(l,m)
             v6_loc(ix)  = v6_loc(ix)  + mesh%gauss%ww(ni,l)*smb(6) + SUM(dw_loc(:,ni)*grad_mom(:,6,l)) &
                  + (dw_loc(1,ni)*tensor_loc(3,2,l) + dw_loc(2,ni)*tensor_loc(3,6,l))*ray*mesh%gauss%rj(l,m)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(4) + SUM(dw_loc(:,ni)*grad_mom(:,4,l)) &
                  + (dw_loc(1,ni)*tensor_loc(2,2,l) + dw_loc(2,ni)*tensor_loc(2,6,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(3) + SUM(dw_loc(:,ni)*grad_mom(:,3,l)) &
                  + (dw_loc(1,ni)*tensor_loc(2,1,l) + dw_loc(2,ni)*tensor_loc(2,5,l))*ray*mesh%gauss%rj(l,m)
          END DO

       ENDDO

       CALL VecSetValues(vb_3_145, 2*nw, idxm_2, v14_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_236, 2*nw, idxm_2, v23_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_145,   nw, idxm_1,  v5_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_236,   nw, idxm_1,  v6_loc, ADD_VALUES, ierr)
    ENDDO

    IF (inputs%LES) THEN
       IF (mesh%edge_stab) THEN
          CALL error_Petsc('BUG in qs_ns_stab_new: terms with edge_stab not yet assembled')
       END IF
    END IF

    CALL VecAssemblyBegin(vb_3_145,ierr)
    CALL VecAssemblyEnd(vb_3_145,ierr)
    CALL VecAssemblyBegin(vb_3_236,ierr)
    CALL VecAssemblyEnd(vb_3_236,ierr)

  END SUBROUTINE qs_ns_momentum_stab_3x3

  SUBROUTINE qs_ns_mom_compute_residual_LES(mesh,vv_1_LA,vv_2_LA,mode,ff,V1m,P,rotb_b, &
       tensor,tensor_gauss,vb_2_14,vb_2_23,vb_1_5,vb_1_6)
    !=================================
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_1_LA,vv_2_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m   !V(noeud, type)
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor, tensor_gauss
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: rotb_b
    INTEGER,                                    INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, ftensor, smb
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(3,6,mesh%gauss%l_G)            :: tensor_loc
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8)   :: ray
    INTEGER        :: m, l , i , ni , index, TYPE, k
    INTEGER        :: nw, ix, ki, iglob
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_2_14,vb_2_23,vb_1_5,vb_1_6
    PetscErrorCode :: ierr

    CALL VecSet(vb_2_14, 0.d0, ierr)
    CALL VecSet(vb_2_23, 0.d0, ierr)
    CALL VecSet(vb_1_5, 0.d0, ierr)
    CALL VecSet(vb_1_6, 0.d0, ierr)

    index = 0
    nw = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_1_LA%loc_to_glob(1,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_2_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !===Compute tensor on gauss points:
          DO TYPE = 1, 6
             DO k = 1, 3
                tensor_loc(k,TYPE,l) = -SUM(tensor(k,j_loc,TYPE)*mesh%gauss%ww(:,l)) &
                     + tensor_gauss(k,index,TYPE)
             END DO
          END DO

          !===Compute residual part of source term, time derivative, pressure gradient without r&z derivatives
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1) = SUM(P(j_loc,1)*dw_loc(1,:))
          ftensor(1) = -mode/ray*tensor_loc(1,4,l) + tensor_loc(2,3,l)/ray
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2) = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2) = SUM(P(j_loc,2)*dw_loc(1,:))
          ftensor(2) = mode/ray*tensor_loc(1,3,l) + tensor_loc(2,4,l)/ray
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3) = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3) = SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          ftensor(3) = -mode/ray*tensor_loc(2,4,l) - tensor_loc(1,3,l)/ray
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4) = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4) = -SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          ftensor(4) = mode/ray*tensor_loc(2,3,l) - tensor_loc(1,4,l)/ray
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5) = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5) = SUM(P(j_loc,1)*dw_loc(2,:))
          ftensor(5) = -mode/ray*tensor_loc(3,4,l)
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6) = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6) = SUM(P(j_loc,2)*dw_loc(2,:))
          ftensor(6) = mode/ray*tensor_loc(3,3,l)
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp-fs+ftensor-rotb_b(index,:))*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(1) + &
                  (dw_loc(1,ni)*tensor_loc(1,1,l) + dw_loc(2,ni)*tensor_loc(1,5,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(2) + &
                  (dw_loc(1,ni)*tensor_loc(1,2,l) + dw_loc(2,ni)*tensor_loc(1,6,l))*ray*mesh%gauss%rj(l,m)
             v5_loc(ix)  = v5_loc(ix)  + mesh%gauss%ww(ni,l)*smb(5) + &
                  (dw_loc(1,ni)*tensor_loc(3,1,l) + dw_loc(2,ni)*tensor_loc(3,5,l))*ray*mesh%gauss%rj(l,m)
             v6_loc(ix)  = v6_loc(ix)  + mesh%gauss%ww(ni,l)*smb(6) + &
                  (dw_loc(1,ni)*tensor_loc(3,2,l) + dw_loc(2,ni)*tensor_loc(3,6,l))*ray*mesh%gauss%rj(l,m)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(4) + &
                  (dw_loc(1,ni)*tensor_loc(2,2,l) + dw_loc(2,ni)*tensor_loc(2,6,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(3) + &
                  (dw_loc(1,ni)*tensor_loc(2,1,l) + dw_loc(2,ni)*tensor_loc(2,5,l))*ray*mesh%gauss%rj(l,m)
          END DO
       ENDDO

       CALL VecSetValues(vb_2_14, 2*nw, idxm_2, v14_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_2_23, 2*nw, idxm_2, v23_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_5,    nw, idxm_1,  v5_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_1_6,    nw, idxm_1,  v6_loc, ADD_VALUES, ierr)
    ENDDO

    CALL VecAssemblyBegin(vb_2_14,ierr)
    CALL VecAssemblyEnd(vb_2_14,ierr)
    CALL VecAssemblyBegin(vb_2_23,ierr)
    CALL VecAssemblyEnd(vb_2_23,ierr)
    CALL VecAssemblyBegin(vb_1_5,ierr)
    CALL VecAssemblyEnd(vb_1_5,ierr)
    CALL VecAssemblyBegin(vb_1_6,ierr)
    CALL VecAssemblyEnd(vb_1_6,ierr)

  END SUBROUTINE qs_ns_mom_compute_residual_LES

  SUBROUTINE qs_ns_mom_compute_residual_LES_3x3(mesh,vv_3_LA,mode,ff,V1m,P,rotb_b, &
       tensor,tensor_gauss,vb_3_145,vb_3_236)
    !=================================
    USE def_type_mesh
    USE chaine_caractere
    USE my_util
    USE input_data
    IMPLICIT NONE
    TYPE(petsc_csr_LA)                                     :: vv_3_LA
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: ff
    TYPE(mesh_type), TARGET                                :: mesh
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: V1m   !V(noeud, type)
    REAL(KIND=8), DIMENSION(:,:,:),             INTENT(IN) :: tensor, tensor_gauss
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: P
    REAL(KIND=8), DIMENSION(:,:),               INTENT(IN) :: rotb_b
    INTEGER,                                    INTENT(IN) :: mode
    REAL(KIND=8), DIMENSION(6)                             :: fs , ft , fp, ftensor, smb
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: j_loc
    REAL(KIND=8), DIMENSION(mesh%gauss%k_d,mesh%gauss%n_w) :: dw_loc
    REAL(KIND=8), DIMENSION(3,6,mesh%gauss%l_G)            :: tensor_loc
    REAL(KIND=8), DIMENSION(2*mesh%gauss%n_w)              :: v14_loc, v23_loc
    INTEGER,      DIMENSION(2*mesh%gauss%n_w)              :: idxm_2
    REAL(KIND=8), DIMENSION(mesh%gauss%n_w)                :: v5_loc, v6_loc
    INTEGER,      DIMENSION(mesh%gauss%n_w)                :: idxm_1
    REAL(KIND=8)   :: ray
    INTEGER        :: m, l , i , ni , index, TYPE, k
    INTEGER        :: nw, ix, ki, iglob
    !#include "petsc/finclude/petsc.h"
    Vec            :: vb_3_145,vb_3_236
    PetscErrorCode :: ierr

    CALL VecSet(vb_3_145, 0.d0, ierr)
    CALL VecSet(vb_3_236, 0.d0, ierr)

    index = 0
    nw = mesh%gauss%n_w
    DO m = 1, mesh%dom_me
       j_loc = mesh%jj(:,m)

       DO ni = 1, nw
          i = j_loc(ni)
          iglob = vv_3_LA%loc_to_glob(3,i)
          idxm_1(ni) = iglob-1
       END DO

       DO ki = 1, 2
          DO ni = 1, nw
             i = j_loc(ni)
             iglob = vv_3_LA%loc_to_glob(ki,i)
             ix = (ki-1)*nw+ni
             idxm_2(ix) = iglob-1
          END DO
       END DO

       v14_loc = 0.d0
       v23_loc = 0.d0
       v5_loc  = 0.d0
       v6_loc  = 0.d0
       DO l = 1, mesh%gauss%l_G
          index  = index +1
          dw_loc = mesh%gauss%dw(:,:,l,m)

          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,j_loc)*mesh%gauss%ww(:,l))

          !===Compute tensor on gauss points:
          DO TYPE = 1, 6
             DO k = 1, 3
                tensor_loc(k,TYPE,l) = -SUM(tensor(k,j_loc,TYPE)*mesh%gauss%ww(:,l)) &
                     + tensor_gauss(k,index,TYPE)
             END DO
          END DO

          !===Compute residual part of source term, time derivative, pressure gradient without r&z derivatives
          !--------calcul de la premiere composante : u0(1,:) <--> f(r,m,c)
          fs(1) = SUM(ff(j_loc,1) * mesh%gauss%ww(:,l))
          ft(1) = SUM(V1m(j_loc,1) * mesh%gauss%ww(:,l))
          fp(1) = SUM(P(j_loc,1)*dw_loc(1,:))
          ftensor(1) = -mode/ray*tensor_loc(1,4,l) + tensor_loc(2,3,l)/ray
          !--------calcul de la seconde composante : u0(2,:) <--> f(r,m,s)
          fs(2) = SUM(ff(j_loc,2) * mesh%gauss%ww(:,l))
          ft(2) = SUM(V1m(j_loc,2) * mesh%gauss%ww(:,l))
          fp(2) = SUM(P(j_loc,2)*dw_loc(1,:))
          ftensor(2) = mode/ray*tensor_loc(1,3,l) + tensor_loc(2,4,l)/ray
          !--------calcul de la troisieme composante : u0(3,:) <--> f(th,m,c)
          fs(3) = SUM(ff(j_loc,3) * mesh%gauss%ww(:,l))
          ft(3) = SUM(V1m(j_loc,3) * mesh%gauss%ww(:,l))
          fp(3) = SUM(P(j_loc,2)*mesh%gauss%ww(:,l))/ray*mode
          ftensor(3) = -mode/ray*tensor_loc(2,4,l) - tensor_loc(1,3,l)/ray
          !--------calcul de la quatrieme composante : u0(4,:) <--> f(th,m,s)
          fs(4) = SUM(ff(j_loc,4) * mesh%gauss%ww(:,l))
          ft(4) = SUM(V1m(j_loc,4) * mesh%gauss%ww(:,l))
          fp(4) = -SUM(P(j_loc,1)*mesh%gauss%ww(:,l))/ray *mode
          ftensor(4) = mode/ray*tensor_loc(2,3,l) - tensor_loc(1,4,l)/ray
          !--------calcul de la cinquieme composante : u0(5,:) <--> f(z,m,c)
          fs(5) = SUM(ff(j_loc,5) * mesh%gauss%ww(:,l))
          ft(5) = SUM(V1m(j_loc,5) * mesh%gauss%ww(:,l))
          fp(5) = SUM(P(j_loc,1)*dw_loc(2,:))
          ftensor(5) = -mode/ray*tensor_loc(3,4,l)
          !--------calcul de la sixieme composante : u0(6,:) <--> f(z,m,s)
          fs(6) = SUM(ff(j_loc,6) * mesh%gauss%ww(:,l))
          ft(6) = SUM(V1m(j_loc,6) * mesh%gauss%ww(:,l))
          fp(6) = SUM(P(j_loc,2)*dw_loc(2,:))
          ftensor(6) = mode/ray*tensor_loc(3,3,l)
          !-------calcul du second membre pour le terme nonlineaire------------------------
          !-------------------------------------------------------------------------------

          smb =  (ft+fp-fs+ftensor-rotb_b(index,:))*ray*mesh%gauss%rj(l,m)

          ki = 1
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(1) + &
                  (dw_loc(1,ni)*tensor_loc(1,1,l) + dw_loc(2,ni)*tensor_loc(1,5,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(2) + &
                  (dw_loc(1,ni)*tensor_loc(1,2,l) + dw_loc(2,ni)*tensor_loc(1,6,l))*ray*mesh%gauss%rj(l,m)
             v5_loc(ix)  = v5_loc(ix)  + mesh%gauss%ww(ni,l)*smb(5) + &
                  (dw_loc(1,ni)*tensor_loc(3,1,l) + dw_loc(2,ni)*tensor_loc(3,5,l))*ray*mesh%gauss%rj(l,m)
             v6_loc(ix)  = v6_loc(ix)  + mesh%gauss%ww(ni,l)*smb(6) + &
                  (dw_loc(1,ni)*tensor_loc(3,2,l) + dw_loc(2,ni)*tensor_loc(3,6,l))*ray*mesh%gauss%rj(l,m)
          END DO

          ki = 2
          DO ni = 1, nw
             ix = (ki-1)*nw+ni
             v14_loc(ix) = v14_loc(ix) + mesh%gauss%ww(ni,l)*smb(4) + &
                  (dw_loc(1,ni)*tensor_loc(2,2,l) + dw_loc(2,ni)*tensor_loc(2,6,l))*ray*mesh%gauss%rj(l,m)
             v23_loc(ix) = v23_loc(ix) + mesh%gauss%ww(ni,l)*smb(3) + &
                  (dw_loc(1,ni)*tensor_loc(2,1,l) + dw_loc(2,ni)*tensor_loc(2,5,l))*ray*mesh%gauss%rj(l,m)
          END DO
       ENDDO

       CALL VecSetValues(vb_3_145, 2*nw, idxm_2, v14_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_236, 2*nw, idxm_2, v23_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_145,   nw, idxm_1,  v5_loc, ADD_VALUES, ierr)
       CALL VecSetValues(vb_3_236,   nw, idxm_1,  v6_loc, ADD_VALUES, ierr)
    ENDDO

    CALL VecAssemblyBegin(vb_3_145,ierr)
    CALL VecAssemblyEnd(vb_3_145,ierr)
    CALL VecAssemblyBegin(vb_3_236,ierr)
    CALL VecAssemblyEnd(vb_3_236,ierr)

  END SUBROUTINE qs_ns_mom_compute_residual_LES_3x3

  SUBROUTINE qs_00_gauss_surface(mesh, vv_1_LA, temp_list_robin_sides, convection_coeff, &
       exterior_temperature, cb) ! MODIFICATION: implementation of the term int_(partial Omega) h*Text*v, with h the convection coefficient
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    TYPE(petsc_csr_LA)                       :: vv_1_LA
    INTEGER     , DIMENSION(:), INTENT(IN)   :: temp_list_robin_sides
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: convection_coeff, exterior_temperature
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws) :: c_loc
    INTEGER,      DIMENSION(mesh%gauss%n_ws) :: idxm
    INTEGER                                  :: nws, ms, ls, ni, i, iglob
    REAL(KIND=8)                             :: ray, x
    INTEGER, DIMENSION(1:1)                  :: coeff_index
    !#include "petsc/finclude/petsc.h"
    Vec                                      :: cb
    PetscErrorCode                           :: ierr

    nws = mesh%gauss%n_ws

    DO ms = 1, mesh%mes
       IF (MINVAL(ABS(temp_list_robin_sides - mesh%sides(ms))) > 0) CYCLE
       c_loc = 0d0
       coeff_index = MINLOC(ABS(temp_list_robin_sides - mesh%sides(ms)))
       DO ls = 1, mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,mesh%jjs(:,ms)) * mesh%gauss%wws(:,ls))
          x = convection_coeff(coeff_index(1)) * exterior_temperature(coeff_index(1)) * &
               ray * mesh%gauss%rjs(ls,ms)
          DO ni = 1, nws
             c_loc(ni) = c_loc(ni) + x * mesh%gauss%wws(ni,ls)
          END DO
       ENDDO
       DO ni = 1, nws
          i = mesh%jjs(ni,ms)
          iglob = vv_1_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO
       CALL VecSetValues(cb, nws, idxm, c_loc, ADD_VALUES, ierr)
    ENDDO

    CALL VecAssemblyBegin(cb,ierr)
    CALL VecAssemblyEnd(cb,ierr)

  END SUBROUTINE qs_00_gauss_surface

  SUBROUTINE qs_00_gauss_surface_conc(mesh, vv_1_LA, conc_list_robin_sides, convection_coeff_conc_rhs, &
       exterior_concentration, cb) 
    !MODIFICATION: implementation of the term int_(partial Omega) h*Text*v, with h the convection coefficient
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    TYPE(petsc_csr_LA)                       :: vv_1_LA
    INTEGER     , DIMENSION(:), INTENT(IN)   :: conc_list_robin_sides
    REAL(KIND=8), DIMENSION(:), INTENT(IN)   :: convection_coeff_conc_rhs, exterior_concentration
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws) :: c_loc
    INTEGER,      DIMENSION(mesh%gauss%n_ws) :: idxm
    INTEGER                                  :: nws, ms, ls, ni, i, iglob
    REAL(KIND=8)                             :: ray, x
    INTEGER, DIMENSION(1:1)                  :: coeff_index

    !#include "petsc/finclude/petsc.h"
    Vec                                      :: cb
    PetscErrorCode                           :: ierr

    nws = mesh%gauss%n_ws

    DO ms = 1, mesh%mes
       IF (MINVAL(ABS(conc_list_robin_sides - mesh%sides(ms))) > 0) CYCLE
       c_loc = 0.d0
       coeff_index = MINLOC(ABS(conc_list_robin_sides - mesh%sides(ms)))
       DO ls = 1, mesh%gauss%l_Gs
          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,mesh%jjs(:,ms)) * mesh%gauss%wws(:,ls))
          x = convection_coeff_conc_rhs(coeff_index(1)) * exterior_concentration(coeff_index(1)) * &
               ray * mesh%gauss%rjs(ls,ms)
          DO ni = 1, nws
             c_loc(ni) = c_loc(ni) + x * mesh%gauss%wws(ni,ls)
          END DO
       ENDDO
       DO ni = 1, nws
          i = mesh%jjs(ni,ms)
          iglob = vv_1_LA%loc_to_glob(1,i)
          idxm(ni) = iglob-1
       END DO
       CALL VecSetValues(cb, nws, idxm, c_loc, ADD_VALUES, ierr)
    ENDDO
    CALL VecAssemblyBegin(cb,ierr)
    CALL VecAssemblyEnd(cb,ierr)
  END SUBROUTINE qs_00_gauss_surface_conc

  SUBROUTINE qs_00_gauss_H_conc(mesh, j_H_to_conc, conc_1_LA, cb_1, cb_2)
    USE def_type_mesh
    USE gauss_points
    USE boundary
    USE input_data
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN)   :: mesh
    REAL(KIND=8), DIMENSION(:,:),          INTENT(IN)   :: j_H_to_conc
    REAL(KIND=8), DIMENSION(4)                          :: j_H_to_conc_gauss
    REAL(KIND=8), DIMENSION(mesh%gauss%n_ws,2)          :: v_loc
    REAL(KIND=8), DIMENSION(2)                          :: normi
    REAL(KIND=8) ::  ray, stab, y
    INTEGER      :: i, ni, ms, ls, nws, iglob
    INTEGER, DIMENSION(mesh%gauss%n_ws)    :: idxn
    TYPE(petsc_csr_LA)                     :: conc_1_LA
    PetscErrorCode          :: ierr
    Vec                     :: cb_1, cb_2

    CALL gauss(mesh)
    nws = mesh%gauss%n_ws

    DO ms = 1, mesh%mes
       IF (MINVAL(ABS(inputs%list_inter_rot_h_jump - mesh%sides(ms))) > 0) CYCLE
       v_loc=0.d0
       stab = -inputs%MA/inputs%faraday_cst
       DO ls = 1, l_Gs
          !===Compute radius of Gauss point
          ray = SUM(mesh%rr(1,mesh%jjs(:,ms))* mesh%gauss%wws(:,ls))
          !===Get normal vector on Gauss point
          normi = rnorms(:,ls,ms)
          !===Compute Current on Gauss point
          j_H_to_conc_gauss(1) = SUM(j_H_to_conc(mesh%jjs(:,ms),1)*mesh%gauss%wws(:,ls))
          j_H_to_conc_gauss(2) = SUM(j_H_to_conc(mesh%jjs(:,ms),2)*mesh%gauss%wws(:,ls))
          j_H_to_conc_gauss(3) = SUM(j_H_to_conc(mesh%jjs(:,ms),5)*mesh%gauss%wws(:,ls))
          j_H_to_conc_gauss(4) = SUM(j_H_to_conc(mesh%jjs(:,ms),6)*mesh%gauss%wws(:,ls))
          !===Compute j dot n
          DO ni = 1, nws
             y = ray*rjs(ls,ms)*stab*wws(ni,ls)
             v_loc(ni,1) = v_loc(ni,1) + ( j_H_to_conc_gauss(1)*normi(1)+j_H_to_conc_gauss(3)*normi(2))*y
             v_loc(ni,2) = v_loc(ni,2) + (j_H_to_conc_gauss(2)*normi(1)+j_H_to_conc_gauss(4)*normi(2))*y
          END DO
       END DO ! end loop on ls
       DO ni = 1, nws
          i = mesh%jjs(ni,ms)
          iglob = conc_1_LA%loc_to_glob(1,i)
          idxn(ni) = iglob-1
       END DO
       CALL VecSetValues(cb_1, nws, idxn, v_loc(:,1), ADD_VALUES, ierr)
       CALL VecSetValues(cb_2, nws, idxn, v_loc(:,2), ADD_VALUES, ierr)
    END DO ! end loop on ms

    CALL VecAssemblyBegin(cb_1,ierr)
    CALL VecAssemblyEnd(cb_1,ierr)
    CALL VecAssemblyBegin(cb_2,ierr)
    CALL VecAssemblyEnd(cb_2,ierr)

  END SUBROUTINE qs_00_gauss_H_conc

END MODULE fem_rhs_axi
