MODULE boundary_generic_module
  IMPLICIT NONE
  INTERFACE
     MODULE SUBROUTINE init_velocity_pressure(mesh_f, mesh_c, time, &
          dt, list_mode, un_m1, un, pn_m1, pn, phin_m1, phin)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh_f, mesh_c
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
     END SUBROUTINE init_velocity_pressure

     MODULE SUBROUTINE init_temperature(mesh, time, dt, list_mode, tempn_m1, tempn)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: tempn_m1, tempn
     END SUBROUTINE init_temperature

     MODULE SUBROUTINE init_concentration(mesh, time, dt, list_mode, concn_m1, concn)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: concn_m1, concn
     END SUBROUTINE init_concentration


     MODULE SUBROUTINE init_level_set(pp_mesh, time, &
          dt, list_mode, level_set_m1, level_set)
       USE def_type_mesh
       TYPE(mesh_type)                              :: pp_mesh
       REAL(KIND=8),                     INTENT(OUT):: time
       REAL(KIND=8),                     INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),       INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT):: level_set, level_set_m1
     END SUBROUTINE init_level_set

     MODULE FUNCTION source_in_NS_momentum(TYPE, rr, mode, i, time, Re, ty, &
          density, tempn, concn) RESULT(vv)
       INTEGER     ,                             INTENT(IN) :: TYPE
       REAL(KIND=8), DIMENSION(:,:),             INTENT(IN) :: rr
       INTEGER     ,                             INTENT(IN) :: mode, i
       REAL(KIND=8),                             INTENT(IN) :: time
       REAL(KIND=8),                             INTENT(IN) :: Re
       CHARACTER(LEN=2),                         INTENT(IN) :: ty
       REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: density
       REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: tempn
       REAL(KIND=8), DIMENSION(:,:,:),           INTENT(IN) :: concn
       REAL(KIND=8), DIMENSION(SIZE(rr,2))                  :: vv
     END FUNCTION source_in_NS_momentum


     MODULE FUNCTION source_in_temperature(TYPE, rr, m, t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION source_in_temperature

     MODULE FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m, interface_nb
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION source_in_level_set

     MODULE FUNCTION vv_exact(TYPE,rr,m,t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION vv_exact

     MODULE FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
     END FUNCTION imposed_velocity_by_penalty

     MODULE FUNCTION pp_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION pp_exact

     MODULE FUNCTION temperature_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION temperature_exact

     MODULE FUNCTION concentration_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION concentration_exact


     MODULE FUNCTION level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m, interface_nb
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION level_set_exact

     MODULE FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
       REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
       INTEGER,                    INTENT(IN)     :: nb_angles
       INTEGER,                    INTENT(IN)     :: nb, ne
       REAL(KIND=8),               INTENT(IN)     :: time
       REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
     END FUNCTION penal_in_real_space

     MODULE FUNCTION extension_velocity(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv
     END FUNCTION extension_velocity

     MODULE FUNCTION extension_temperature(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv
     END FUNCTION extension_temperature

     MODULE FUNCTION extension_concentration(TYPE, vv_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: vv_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(vv_Mesh%np)               :: vv
     END FUNCTION extension_concentration

     MODULE FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
       USE def_type_mesh
       TYPE(mesh_type),                       INTENT(IN) :: H_mesh
       INTEGER,                               INTENT(IN) :: m
       REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv
     END FUNCTION Vexact

     MODULE FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv)
       CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
     END FUNCTION H_B_quasi_static

     MODULE FUNCTION Hexact(H_mesh, TYPE, rr, m, mu_H_field, t) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION Hexact

     MODULE FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION Phiexact

     MODULE FUNCTION Jexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t, &
          mesh_id, opt_B_ext) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
       INTEGER     ,                        INTENT(IN)   :: mesh_id
       REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext
       REAL(KIND=8)                                      :: vv
     END FUNCTION Jexact_gauss

     MODULE FUNCTION Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
       REAL(KIND=8)                                      :: vv
     END FUNCTION Eexact_gauss

     MODULE SUBROUTINE init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
          list_mode, Hn1, Hn, phin1, phin)
       USE def_type_mesh
       TYPE(mesh_type)                            :: H_mesh, phi_mesh
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu_H_field
       REAL(KIND=8),                   INTENT(IN) :: mu_phi
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: Hn, Hn1
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: phin, phin1
     END SUBROUTINE init_maxwell

     MODULE FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       INTEGER,     INTENT(IN)                    :: nb, ne
       REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
       INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids
       REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
     END FUNCTION mu_bar_in_fourier_space

     MODULE FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
       REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
       INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
       REAL(KIND=8),DIMENSION(2)            :: vv
     END FUNCTION grad_mu_bar_in_fourier_space

     MODULE FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
       INTEGER, INTENT(IN)                        :: nb_angles
       INTEGER, INTENT(IN)                        :: nb, ne
       REAL(KIND=8), INTENT(IN)                   :: time
       REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
     END FUNCTION mu_in_real_space

     MODULE FUNCTION sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv
     END FUNCTION sigma_bar_in_fourier_space

     MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION chi_coeff_law

     MODULE FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION T_dchi_dT_coeff_law

     MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION nu_tilde_law

     MODULE FUNCTION rot_H_jump_interface(mesh,rr,list_mode) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type)                                   :: mesh
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,      DIMENSION(:),          INTENT(IN)   :: list_mode
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6,SIZE(list_mode)) :: vv
     END FUNCTION rot_H_jump_interface

     MODULE FUNCTION Derivative_of_potential_from_rhoLi(delta_rhoLi_phys) RESULT(vv)
       REAL(KIND=8) :: delta_rhoLi_phys
       REAL(KIND=8) :: vv
     END FUNCTION Derivative_of_potential_from_rhoLi

     MODULE FUNCTION molar_fraction_from_concentration(delta_rhoLi_phys) RESULT(vv)
       REAL(KIND=8) :: delta_rhoLi_phys
       REAL(KIND=8) :: vv
     END FUNCTION molar_fraction_from_concentration

     MODULE FUNCTION curved_boundary_radius(interface,theta) RESULT(vv)
        INTEGER, INTENT(IN)  :: interface
        REAL(KIND=8), INTENT(IN)   :: theta
        REAL(KIND = 8) :: vv
     END FUNCTION curved_boundary_radius
  END INTERFACE
END MODULE boundary_generic_module
