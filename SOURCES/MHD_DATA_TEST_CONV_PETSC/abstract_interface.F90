MODULE abstract_interface
  ABSTRACT INTERFACE
     SUBROUTINE sub_init_velocity_pressure(mesh_f, mesh_c, time, &
          dt, list_mode, un_m1, un, pn_m1, pn, phin_m1, phin)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh_f, mesh_c
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
     END SUBROUTINE sub_init_velocity_pressure

     SUBROUTINE sub_init_temperature(mesh, time, dt, list_mode, tempn_m1, tempn)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: tempn_m1, tempn
     END SUBROUTINE sub_init_temperature

     SUBROUTINE sub_init_concentration(mesh, time, dt, list_mode, concn_m1, concn)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8),                   INTENT(OUT):: time
       REAL(KIND=8),                   INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: concn_m1, concn
     END SUBROUTINE sub_init_concentration


     SUBROUTINE sub_init_level_set(pp_mesh, time, &
          dt, list_mode, level_set_m1, level_set)
       USE def_type_mesh
       TYPE(mesh_type)                              :: pp_mesh
       REAL(KIND=8),                     INTENT(OUT):: time
       REAL(KIND=8),                     INTENT(IN) :: dt
       INTEGER,      DIMENSION(:),       INTENT(IN) :: list_mode
       REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT):: level_set, level_set_m1
     END SUBROUTINE sub_init_level_set

     FUNCTION sub_source_in_NS_momentum(TYPE, rr, mode, i, time, Re, ty, &
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
     END FUNCTION sub_source_in_NS_momentum

     FUNCTION sub_source_in_temperature(TYPE, rr, m, t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_source_in_temperature

     FUNCTION sub_source_in_level_set(interface_nb,TYPE, rr, m, t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m, interface_nb
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_source_in_level_set

     FUNCTION sub_vv_exact(TYPE,rr,m,t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_vv_exact

     FUNCTION sub_imposed_velocity_by_penalty(rr,t) RESULT(vv)
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
     END FUNCTION sub_imposed_velocity_by_penalty

     FUNCTION sub_pp_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_pp_exact

     FUNCTION sub_temperature_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_temperature_exact

     FUNCTION sub_concentration_exact(TYPE,rr,m,t) RESULT (vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_concentration_exact


     FUNCTION sub_level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m, interface_nb
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_level_set_exact

     FUNCTION sub_penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type)                            :: mesh
       REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
       REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
       INTEGER,                    INTENT(IN)     :: nb_angles
       INTEGER,                    INTENT(IN)     :: nb, ne
       REAL(KIND=8),               INTENT(IN)     :: time
       REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
     END FUNCTION sub_penal_in_real_space

     FUNCTION sub_extension_velocity(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv
     END FUNCTION sub_extension_velocity

     FUNCTION sub_extension_temperature(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv
     END FUNCTION sub_extension_temperature

     FUNCTION sub_extension_concentration(TYPE, vv_mesh, mode, t, n_start) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: vv_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
       INTEGER,                             INTENT(IN)   :: mode
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(vv_Mesh%np)                :: vv
     END FUNCTION sub_extension_concentration

     FUNCTION sub_Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
       USE def_type_mesh
       TYPE(mesh_type),                       INTENT(IN) :: H_mesh
       INTEGER,                               INTENT(IN) :: m
       REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv
     END FUNCTION sub_Vexact

     FUNCTION sub_H_B_quasi_static(char_h_b, rr, m) RESULT(vv)
       CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv
     END FUNCTION sub_H_B_quasi_static

     FUNCTION sub_Hexact(H_mesh, TYPE, rr, m, mu_H_field, t) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: t
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_Hexact

     FUNCTION sub_Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
       REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
     END FUNCTION sub_Phiexact

     FUNCTION sub_Jexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t, &
          mesh_id, opt_B_ext) RESULT(vv)
       INTEGER     ,                        INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
       INTEGER     ,                        INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
       INTEGER     ,                        INTENT(IN)   :: mesh_id
       REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext
       REAL(KIND=8)                                      :: vv
     END FUNCTION sub_Jexact_gauss

     FUNCTION sub_Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
       INTEGER,                             INTENT(IN)   :: TYPE
       REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
       INTEGER,                             INTENT(IN)   :: m
       REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
       REAL(KIND=8)                                      :: vv
     END FUNCTION sub_Eexact_gauss

     SUBROUTINE sub_init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
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
     END SUBROUTINE sub_init_maxwell

     FUNCTION sub_mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       INTEGER,     INTENT(IN)                    :: nb, ne
       REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
       INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids
       REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
     END FUNCTION sub_mu_bar_in_fourier_space

     FUNCTION sub_grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
       REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
       INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
       REAL(KIND=8),DIMENSION(2)            :: vv
     END FUNCTION sub_grad_mu_bar_in_fourier_space

     FUNCTION sub_mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
       INTEGER, INTENT(IN)                        :: nb_angles
       INTEGER, INTENT(IN)                        :: nb, ne
       REAL(KIND=8), INTENT(IN)                   :: time
       REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv
     END FUNCTION sub_mu_in_real_space

     FUNCTION sub_sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type), INTENT(IN)                :: H_mesh
       REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv
     END FUNCTION sub_sigma_bar_in_fourier_space

     FUNCTION sub_chi_coeff_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION sub_chi_coeff_law

     FUNCTION sub_T_dchi_dT_coeff_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION sub_T_dchi_dT_coeff_law

     FUNCTION sub_nu_tilde_law(temp) RESULT(vv)
       REAL(KIND=8) :: temp
       REAL(KIND=8) :: vv
     END FUNCTION sub_nu_tilde_law

     FUNCTION sub_rot_H_jump_interface(mesh,rr,list_mode) RESULT(vv)
       USE def_type_mesh
       TYPE(mesh_type)                                   :: mesh
       REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
       INTEGER,      DIMENSION(:),          INTENT(IN)   :: list_mode
       REAL(KIND=8), DIMENSION(SIZE(rr,2),6,SIZE(list_mode)) :: vv
     END FUNCTION sub_rot_H_jump_interface

     FUNCTION sub_Derivative_of_potential_from_rhoLi(delta_rhoLi_phys) RESULT(vv)
       REAL(KIND=8) :: delta_rhoLi_phys
       REAL(KIND=8) :: vv
     END FUNCTION sub_Derivative_of_potential_from_rhoLi

     FUNCTION sub_molar_fraction_from_concentration(delta_rhoLi_phys) RESULT(vv)
       REAL(KIND=8) :: delta_rhoLi_phys
       REAL(KIND=8) :: vv
     END FUNCTION  sub_molar_fraction_from_concentration

  END INTERFACE
END MODULE abstract_interface
