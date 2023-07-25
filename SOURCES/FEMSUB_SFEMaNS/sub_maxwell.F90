MODULE update_maxwell
  PUBLIC:: maxwell_decouple
  PRIVATE

CONTAINS

  SUBROUTINE maxwell_decouple(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, stab_jump_h, sigma_in, &
       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H, conc_to_H)

    USE def_type_mesh
    USE periodic
    USE input_data
    USE update_maxwell_with_B
    USE update_maxwell_with_H
    USE update_maxwell_mxs_with_H
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: H_mesh, phi_mesh, pmag_mesh
    TYPE(interface_type),           INTENT(IN)     :: interface_H_phi, interface_H_mu
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: vel
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Hn, Hn1
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Bn, Bn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: phin, phin1
    REAL(KIND=8), DIMENSION(3),     INTENT(IN)     :: stab_in
    REAL(KIND=8),                   INTENT(IN)     :: stab_jump_h
    REAL(KIND=8),                   INTENT(IN)     :: R_fourier
    INTEGER,                        INTENT(IN)     :: index_fourier
    REAL(KIND=8),                   INTENT(IN)     :: mu_phi, time, dt, Rem
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: sigma_in, mu_H_field
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: one_over_sigma_ns
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: jj_v_to_H
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: conc_to_H
    TYPE(periodic_type),            INTENT(IN)     :: H_phi_per
    TYPE(petsc_csr_LA)                             :: LA_H, LA_pmag, LA_phi, LA_mhd
#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d

    IF (inputs%type_pb=='mhs') THEN
       IF (inputs%if_maxwell_with_H) THEN
          CALL maxwell_mxs_with_H(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
               interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, stab_jump_h, sigma_in, &
               R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
               H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H, conc_to_H)
       ELSE
          CALL  error_petsc('Maxwell_with_B for typ_pb=mhs not programmed yet')
       ENDIF
    ELSE
       IF (inputs%if_maxwell_with_H) THEN
          CALL maxwell_decouple_with_H(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
               interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, sigma_in, &
               R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
               H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H)
       ELSE
          CALL maxwell_decouple_with_B(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
               interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, sigma_in, &
               R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
               H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H)
       END IF
    ENDIF

  END SUBROUTINE maxwell_decouple

END MODULE update_maxwell
