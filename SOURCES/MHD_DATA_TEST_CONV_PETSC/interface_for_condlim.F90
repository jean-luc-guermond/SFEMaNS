module interface_boundary
  use abstract_interface
  PROCEDURE(sub_init_velocity_pressure), POINTER :: init_velocity_pressure
  PROCEDURE(sub_init_temperature), POINTER :: init_temperature
  PROCEDURE(sub_init_concentration), POINTER :: init_concentration
  PROCEDURE(sub_init_level_set), POINTER :: init_level_set
  PROCEDURE(sub_source_in_NS_momentum), POINTER :: source_in_NS_momentum
  PROCEDURE(sub_source_in_temperature), POINTER :: source_in_temperature
  PROCEDURE(sub_source_in_level_set), POINTER :: source_in_level_set
  PROCEDURE(sub_vv_exact), POINTER :: vv_exact
  PROCEDURE(sub_imposed_velocity_by_penalty), POINTER :: imposed_velocity_by_penalty
  PROCEDURE(sub_pp_exact), POINTER :: pp_exact
  PROCEDURE(sub_temperature_exact), POINTER :: temperature_exact
  PROCEDURE(sub_concentration_exact), POINTER :: concentration_exact
  PROCEDURE(sub_level_set_exact), POINTER :: level_set_exact
  PROCEDURE(sub_penal_in_real_space), POINTER :: penal_in_real_space
  PROCEDURE(sub_extension_velocity), POINTER :: extension_velocity
  PROCEDURE(sub_extension_temperature), POINTER :: extension_temperature
  PROCEDURE(sub_extension_concentration), POINTER :: extension_concentration
  PROCEDURE(sub_Vexact), POINTER :: Vexact
  PROCEDURE(sub_H_B_quasi_static), POINTER :: H_B_quasi_static
  PROCEDURE(sub_Hexact), POINTER :: Hexact
  PROCEDURE(sub_Phiexact), POINTER :: Phiexact
  PROCEDURE(sub_Jexact_gauss), POINTER :: Jexact_gauss
  PROCEDURE(sub_Eexact_gauss), POINTER :: Eexact_gauss
  PROCEDURE(sub_init_maxwell), POINTER :: init_maxwell
  PROCEDURE(sub_mu_bar_in_fourier_space), POINTER :: mu_bar_in_fourier_space
  PROCEDURE(sub_grad_mu_bar_in_fourier_space), POINTER :: grad_mu_bar_in_fourier_space
  PROCEDURE(sub_mu_in_real_space), POINTER :: mu_in_real_space
  PROCEDURE(sub_sigma_bar_in_fourier_space), POINTER :: sigma_bar_in_fourier_space
  PROCEDURE(sub_chi_coeff_law), POINTER :: chi_coeff_law
  PROCEDURE(sub_T_dchi_dT_coeff_law), POINTER :: T_dchi_dT_coeff_law
  PROCEDURE(sub_nu_tilde_law), POINTER :: nu_tilde_law
  PROCEDURE(sub_rot_H_jump_interface), POINTER :: rot_H_jump_interface
  PROCEDURE(sub_Derivative_of_potential_from_rhoLi), POINTER :: Derivative_of_potential_from_rhoLi
  PROCEDURE(sub_molar_fraction_from_concentration), POINTER :: molar_fraction_from_concentration
END module interface_boundary

module point_to_boundary_generic
  USE boundary_generic_module, ONLY : init_velocity_pressure_generic => init_velocity_pressure, &
       init_temperature_generic => init_temperature, &
       init_concentration_generic => init_concentration, &
       init_level_set_generic =>  init_level_set, &
       source_in_NS_momentum_generic =>  source_in_NS_momentum, &
       source_in_temperature_generic =>  source_in_temperature, &
       source_in_level_set_generic =>  source_in_level_set, &
       vv_exact_generic =>  vv_exact, &
       imposed_velocity_by_penalty_generic =>  imposed_velocity_by_penalty, &
       pp_exact_generic =>  pp_exact, &
       temperature_exact_generic =>  temperature_exact, &
       concentration_exact_generic => concentration_exact, &
       level_set_exact_generic =>  level_set_exact, &
       penal_in_real_space_generic =>  penal_in_real_space, &
       extension_velocity_generic =>  extension_velocity, &
       extension_temperature_generic =>  extension_temperature, &
       extension_concentration_generic =>  extension_concentration, &
       Vexact_generic =>  Vexact, &
       H_B_quasi_static_generic =>  H_B_quasi_static, &
       Hexact_generic =>  Hexact, &
       Phiexact_generic =>  Phiexact, &
       Jexact_gauss_generic =>  Jexact_gauss, &
       Eexact_gauss_generic =>  Eexact_gauss, &
       init_maxwell_generic =>  init_maxwell, &
       mu_bar_in_fourier_space_generic =>  mu_bar_in_fourier_space, &
       grad_mu_bar_in_fourier_space_generic =>  grad_mu_bar_in_fourier_space, &
       mu_in_real_space_generic =>  mu_in_real_space, &
       sigma_bar_in_fourier_space_generic =>  sigma_bar_in_fourier_space, &
       chi_coeff_law_generic => chi_coeff_law, &
       T_dchi_dT_coeff_law_generic => T_dchi_dT_coeff_law, &
       nu_tilde_law_generic => nu_tilde_law, &
       rot_H_jump_interface_generic => rot_H_jump_interface, &
       Derivative_of_potential_from_rhoLi_generic => Derivative_of_potential_from_rhoLi, &
       molar_fraction_from_concentration_generic => molar_fraction_from_concentration
END module point_to_boundary_generic
