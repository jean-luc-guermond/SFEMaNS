module boundary
  use interface_boundary
  use point_to_boundary_generic

CONTAINS
  SUBROUTINE assign_boundary
    USE input_data
    USE my_util
    IMPLICIT NONE
    INTEGER :: i
    IF (inputs%test_de_convergence) THEN
       i=inputs%numero_du_test_debug
    ELSE
       i=0
    END IF
    SELECT case(i)
    case(0)
       init_velocity_pressure => init_velocity_pressure_generic
       init_temperature => init_temperature_generic
       init_concentration  => init_concentration_generic
       init_level_set =>  init_level_set_generic
       source_in_NS_momentum =>  source_in_NS_momentum_generic
       source_in_temperature =>  source_in_temperature_generic
       source_in_level_set =>  source_in_level_set_generic
       vv_exact =>  vv_exact_generic
       imposed_velocity_by_penalty =>  imposed_velocity_by_penalty_generic
       pp_exact =>  pp_exact_generic
       temperature_exact=>  temperature_exact_generic
       concentration_exact => concentration_exact_generic
       level_set_exact =>  level_set_exact_generic
       penal_in_real_space =>  penal_in_real_space_generic
       extension_velocity =>  extension_velocity_generic
       extension_temperature =>  extension_temperature_generic
       extension_concentration =>  extension_concentration_generic
       Vexact =>  Vexact_generic
       H_B_quasi_static =>  H_B_quasi_static_generic
       Hexact =>  Hexact_generic
       Phiexact =>  Phiexact_generic
       Jexact_gauss =>  Jexact_gauss_generic
       Eexact_gauss =>  Eexact_gauss_generic
       init_maxwell =>  init_maxwell_generic
       mu_bar_in_fourier_space =>  mu_bar_in_fourier_space_generic
       grad_mu_bar_in_fourier_space=>  grad_mu_bar_in_fourier_space_generic
       mu_in_real_space=>  mu_in_real_space_generic
       sigma_bar_in_fourier_space =>  sigma_bar_in_fourier_space_generic
       chi_coeff_law => chi_coeff_law_generic
       T_dchi_dT_coeff_law => T_dchi_dT_coeff_law_generic
       nu_tilde_law => nu_tilde_law_generic
       rot_H_jump_interface => rot_H_jump_interface_generic
       Derivative_of_potential_from_rhoLi => Derivative_of_potential_from_rhoLi_generic
       molar_fraction_from_concentration => molar_fraction_from_concentration_generic

    CASE DEFAULT
       CALL error_petsc(' BUG in assign_boundary, wrong test number')
    end SELECT
  END SUBROUTINE assign_boundary
end module boundary
