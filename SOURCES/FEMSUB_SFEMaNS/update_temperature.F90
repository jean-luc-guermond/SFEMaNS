MODULE update_temperature
  PUBLIC :: temperature_decouple
  PRIVATE

CONTAINS

  SUBROUTINE temperature_decouple(comm_one_d_temp, time, temp_1_LA, &
       list_mode, temp_mesh, tempn_m1, tempn, &
       v_to_energy, H_to_energy, pdt_H_to_energy, &
       vol_heat_capacity_field, temperature_diffusivity_field, temp_per, &
       heat_density_ns_m1, heat_density_ns, heat_density_ns_p1, heat_diffusivity_ns, &
       jj_v_to_temp)
    USE def_type_mesh
    USE periodic
    USE input_data
    USE subroutine_temperature_with_T
    USE subroutine_temperature_with_e
    USE my_util

    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: temp_mesh
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(temp_mesh%np,2,SIZE(list_mode)), INTENT(INOUT)  :: tempn_m1, tempn
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: v_to_energy, H_to_energy, pdt_H_to_energy
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: vol_heat_capacity_field
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: temperature_diffusivity_field
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: heat_density_ns_m1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: heat_density_ns
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: heat_density_ns_p1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: heat_diffusivity_ns
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: jj_v_to_temp
    REAL(KIND=8),                   INTENT(IN)     :: time
    TYPE(periodic_type),            INTENT(IN)     :: temp_per
    TYPE(petsc_csr_LA)                             :: temp_1_LA
#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d_temp


    IF (inputs%if_temperature_with_T) THEN
       CALL three_level_temperature_with_T(comm_one_d_temp, time, temp_1_LA, inputs%dt, list_mode, &
            temp_mesh, tempn_m1, tempn, v_to_energy, H_to_energy, pdt_H_to_energy, &
            vol_heat_capacity_field, temperature_diffusivity_field, &
            inputs%my_par_temperature, inputs%temperature_list_dirichlet_sides, &
            inputs%temperature_list_robin_sides, inputs%convection_coeff, &
            inputs%exterior_temperature, temp_per)
    ELSE
       CALL three_level_temperature_with_e(comm_one_d_temp, time, temp_1_LA, inputs%dt, list_mode, &
            temp_mesh, tempn_m1, tempn, v_to_energy, H_to_energy, pdt_H_to_energy, &
            vol_heat_capacity_field, temperature_diffusivity_field, &
            inputs%my_par_temperature, inputs%temperature_list_dirichlet_sides, &
            inputs%temperature_list_robin_sides, inputs%convection_coeff, &
            inputs%exterior_temperature, temp_per, heat_density_ns_m1, heat_density_ns, &
            heat_density_ns_p1, heat_diffusivity_ns, jj_v_to_temp)
    END IF

  END SUBROUTINE temperature_decouple

END MODULE update_temperature
