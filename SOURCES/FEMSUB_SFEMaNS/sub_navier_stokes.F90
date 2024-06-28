MODULE update_navier_stokes
  PUBLIC :: navier_stokes_decouple
  PRIVATE

CONTAINS

  SUBROUTINE navier_stokes_decouple(comm_one_d_ns,time, vv_3_LA, pp_1_LA, &
       list_mode, pp_mesh, vv_mesh, incpn_m1, incpn, pn_m1, pn, un_m1, un, &
       vvz_per, pp_per, Hn_p2, Bn_p2, density_m2, density_m1, density, visco_dyn, tempn, concn,&
       level_set_m1, level_set, visc_entro_level, level_set_reg, visc_LES)

    USE def_type_mesh
    USE periodic
    USE input_data
    USE subroutine_ns_with_u
    USE subroutine_ns_with_m
    USE subroutine_ns_with_m_art_comp
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: vv_mesh, pp_mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(vv_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: un_m1, un
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode)), INTENT(IN)     :: density_m2, density_m1, density
    REAL(KIND=8), DIMENSION(vv_mesh%np,2,SIZE(list_mode)), INTENT(IN)     :: visco_dyn
    REAL(KIND=8), DIMENSION(:,:,:),                        INTENT(IN)     :: tempn
    REAL(KIND=8), DIMENSION(:,:,:),                        INTENT(IN)     :: concn
    REAL(KIND=8), DIMENSION(pp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT)  :: pn_m1, pn, incpn_m1, incpn
    REAL(KIND=8), DIMENSION(:,:,:),                        INTENT(IN)     :: Hn_p2
    REAL(KIND=8), DIMENSION(:,:,:),                        INTENT(IN)     :: Bn_p2
    REAL(KIND=8), DIMENSION(:,:,:,:),                      INTENT(IN)     :: level_set_m1, level_set
    REAL(KIND=8), DIMENSION(:,:,:,:),                      INTENT(IN)     :: level_set_reg
    REAL(KIND=8), DIMENSION(:,:,:,:),                      INTENT(INOUT)  :: visc_LES
    REAL(KIND=8),                   INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(:,:),   INTENT(OUT)    :: visc_entro_level
    TYPE(periodic_type),            INTENT(IN)     :: vvz_per, pp_per
    TYPE(petsc_csr_LA)                             :: vv_3_LA, pp_1_LA
#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_ns


    IF (inputs%if_navier_stokes_art_comp) THEN !Artificial compression for multiphase Navier-Stokes equations
       IF (.NOT.inputs%if_navier_stokes_with_u) THEN
          CALL BDF1_art_comp_with_m(comm_one_d_ns,time,vv_3_LA, pp_1_LA, vvz_per, pp_per, &
               inputs%dt, inputs%Re, list_mode, pp_mesh, vv_mesh, &
               pn_m1, pn, un_m1, un, Hn_p2, Bn_p2, tempn, concn, density_m2, density_m1, density,&
               visco_dyn, level_set, visc_entro_level, level_set_reg, visc_LES)
       ELSE
          CALL error_petsc('Artificial compressibility with velocity not programmed yet')
       END IF
    ELSE !Projection-correction scheme for Navier-Stokes equations
       IF (inputs%if_navier_stokes_with_u) THEN
          CALL BDF2_ns_stress_bc_with_u(comm_one_d_ns,time,vv_3_LA, pp_1_LA, vvz_per, pp_per, &
               inputs%dt, inputs%Re, list_mode, pp_mesh, vv_mesh, incpn_m1, incpn, &
               pn_m1, pn, un_m1, un, Hn_p2, Bn_p2, density, tempn, concn, visc_LES)
       ELSE
          CALL three_level_ns_tensor_sym_with_m(comm_one_d_ns,time,vv_3_LA, pp_1_LA, &
               inputs%dt, inputs%Re, list_mode, pp_mesh, vv_mesh, incpn_m1, incpn, &
               pn_m1, pn, un_m1, un, Hn_p2, Bn_p2, density_m2, density_m1, density, visco_dyn, tempn, concn, &
               level_set_m1, level_set, visc_entro_level, level_set_reg, visc_LES)
       END IF
    END IF
  END SUBROUTINE navier_stokes_decouple

END MODULE update_navier_stokes
