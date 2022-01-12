MODULE update_maxwell
  PUBLIC:: maxwell_decouple
  PRIVATE

CONTAINS

  SUBROUTINE maxwell_decouple(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, sigma_in, &
       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
!TEST LC-CN 15/12/2016
!!$       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, sigma_ns, jj_v_to_H)
       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H)
!TEST LC-CN 15/12/2016

    USE def_type_mesh
    USE periodic
    USE input_data
    USE update_maxwell_with_B
    USE update_maxwell_with_H
    IMPLICIT NONE
    TYPE(mesh_type),                INTENT(IN)     :: H_mesh, phi_mesh, pmag_mesh
    TYPE(interface_type),           INTENT(IN)     :: interface_H_phi, interface_H_mu
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: vel  
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Hn, Hn1
    REAL(KIND=8), DIMENSION(H_mesh%np,6,SIZE(list_mode)), INTENT(INOUT)  :: Bn, Bn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(INOUT)  :: phin, phin1
    REAL(KIND=8), DIMENSION(3),     INTENT(IN)     :: stab_in 
    REAL(KIND=8),                   INTENT(IN)     :: R_fourier
    INTEGER,                        INTENT(IN)     :: index_fourier
    REAL(KIND=8),                   INTENT(IN)     :: mu_phi, time, dt, Rem     
    REAL(KIND=8), DIMENSION(:),     INTENT(IN)     :: sigma_in, mu_H_field
! TEST
!TEST LC-CN 15/12/2016
!!$    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: sigma_ns
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)     :: one_over_sigma_ns
!TEST LC-CN 15/12/2016
    INTEGER,      DIMENSION(:),     INTENT(IN)     :: jj_v_to_H
! TEST
    !jan 29/JLG+FL/Forget about it/We replace it by H_p_phi_per/Feb 2 2010
    TYPE(periodic_type),            INTENT(IN)     :: H_phi_per
    !jan 29/JLG+FL/Forget about it/We replace it by H_p_phi_per/Feb 2 2010
    TYPE(petsc_csr_LA)                             :: LA_H, LA_pmag, LA_phi, LA_mhd
#include "petsc/finclude/petsc.h"
    MPI_Comm, DIMENSION(:), POINTER  :: comm_one_d

    IF (inputs%if_maxwell_with_H) THEN
       CALL maxwell_decouple_with_H(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, sigma_in, &
       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
!TEST LC-CN 15/12/2016
!!$       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, sigma_ns, jj_v_to_H)
       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H)
!TEST LC-CN 15/12/2016
!!$       CALL maxwell_decouple_with_H(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
!!$       interface_H_mu, Hn, phin, Hn1, phin1, vel, stab_in, sigma_in, &
!!$       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
!!$       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd)
    ELSE
! TEST
       CALL maxwell_decouple_with_B(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
       interface_H_mu, Hn, Bn, phin, Hn1, Bn1, phin1, vel, stab_in, sigma_in, &
       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
!TEST LC-CN 15/12/2016
!!$       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, sigma_ns, jj_v_to_H)
       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd, one_over_sigma_ns, jj_v_to_H)
!TEST LC-CN 15/12/2016
!!$       CALL maxwell_decouple_with_B(comm_one_d, H_mesh, pmag_mesh, phi_mesh, interface_H_phi, &
!!$       interface_H_mu, Hn, phin, Hn1, phin1, vel, stab_in, sigma_in, &
!!$       R_fourier, index_fourier, mu_H_field, mu_phi, time, dt, Rem, list_mode, &
!!$       H_phi_per, LA_H, LA_pmag, LA_phi, LA_mhd)
! TEST
    END IF

  END SUBROUTINE maxwell_decouple

END MODULE update_maxwell
