SUBMODULE (boundary_generic_module) boundary_generic
  USE my_util
  USE def_type_mesh
  USE input_data
  USE bessel
  USE user_data

  REAL(KIND=8), PARAMETER   :: pi = 3.14159265358979323846d0
  REAL(KIND=8), PARAMETER   :: A_plus=2.d0
  REAL(KIND=8), PARAMETER   :: A_minus=1.d0
  REAL(KIND=8), PARAMETER   :: rmax=0.5d0

  !===Dummy variables to avoid warning
  REAL(KIND=8) :: rd
  INTEGER      :: nd
  CHARACTER(LEN=2)  :: cd2
  CHARACTER(LEN=1)  :: cd1
  !===Dummy variables to avoid warning
CONTAINS
  !===============================================================================
  !                       Boundary conditions for Navier-Stokes
  !===============================================================================

  !===Initialize velocity, pressure
  MODULE SUBROUTINE init_velocity_pressure(mesh_f, mesh_c, time, dt, list_mode, &
       un_m1, un, pn_m1, pn, phin_m1, phin)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh_f, mesh_c
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: un_m1, un
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: pn_m1, pn, phin_m1, phin
    INTEGER                                    :: mode, i, j
    REAL(KIND=8), DIMENSION(mesh_c%np)         :: pn_m2

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 6
          !===velocity
          un_m1(:,j,i) = vv_exact(j,mesh_f%rr,mode,time-dt)
          un   (:,j,i) = vv_exact(j,mesh_f%rr,mode,time)
       END DO
       DO j = 1, 2
          !===pressure
          pn_m2(:)       = pp_exact(j,mesh_c%rr,mode,time-2*dt)
          pn_m1  (:,j,i) = pp_exact(j,mesh_c%rr,mode,time-dt)
          pn     (:,j,i) = pp_exact(j,mesh_c%rr,mode,time)
          phin_m1(:,j,i) = pn_m1(:,j,i) - pn_m2(:)
          phin   (:,j,i) = pn   (:,j,i) - pn_m1(:,j,i)
       ENDDO
    ENDDO
  END SUBROUTINE init_velocity_pressure

  !===Initialize temperature
  MODULE SUBROUTINE init_temperature(mesh, time, dt, list_mode, tempn_m1, tempn)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: tempn_m1, tempn
    INTEGER                                    :: mode, i, j

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 2
          tempn_m1(:,j,i) = temperature_exact(j, mesh%rr, mode, time-dt)
          tempn   (:,j,i) = temperature_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
  END SUBROUTINE init_temperature

  MODULE SUBROUTINE init_concentration(mesh, time, dt, list_mode, concn_m1, concn)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: concn_m1, concn
    INTEGER                                    :: mode, i, j

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 2
          concn_m1(:,j,i) = concentration_exact(j, mesh%rr, mode, time-dt)
          concn   (:,j,i) = concentration_exact(j, mesh%rr, mode, time)
       ENDDO
    ENDDO
  END SUBROUTINE init_concentration

  !===Initialize level_set
  MODULE SUBROUTINE init_level_set(pp_mesh, time, &
       dt, list_mode, level_set_m1, level_set)
    IMPLICIT NONE
    TYPE(mesh_type)                              :: pp_mesh
    REAL(KIND=8),                     INTENT(OUT):: time
    REAL(KIND=8),                     INTENT(IN) :: dt
    INTEGER,      DIMENSION(:),       INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:,:), INTENT(OUT):: level_set, level_set_m1
    INTEGER                                      :: mode, i, j, n

    time = 0.d0
    DO i= 1, SIZE(list_mode)
       mode = list_mode(i)
       DO j = 1, 2
          !===level_set
          DO n = 1, inputs%nb_fluid -1
             level_set_m1(n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time-dt)
             level_set   (n,:,j,i) = level_set_exact(n,j,pp_mesh%rr,mode,time)
          END DO
       END DO
    END DO
  END SUBROUTINE init_level_set

  !===Source in momemtum equation. Always called.
  MODULE FUNCTION source_in_NS_momentum(TYPE, rr, mode, i, time, Re, ty, density, tempn, concn) RESULT(vv)
    IMPLICIT NONE
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

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=mode; nd=i; rd=time; rd=Re; cd2=ty
    nd=SIZE(density,1); nd=SIZE(tempn,1); nd=SIZE(concn,1)
    !===Dummy variables to avoid warning
  END FUNCTION source_in_NS_momentum

  !===Extra source in temperature equation. Always called.
  MODULE FUNCTION source_in_temperature(TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('source_in_temperature: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION source_in_temperature

  MODULE FUNCTION source_in_level_set(interface_nb,TYPE, rr, m, t)RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('source_in_level_set: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; nd=interface_nb; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION source_in_level_set

  !===Velocity for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize velocity in: init_velocity_pressure_temperature
  MODULE FUNCTION vv_exact(TYPE,rr,m,t) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION vv_exact

  !===Solid velocity imposed when using penalty technique
  !===Defined in Fourier space on mode 0 only.
  MODULE FUNCTION imposed_velocity_by_penalty(rr,t) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=SIZE(rr,1); rd=t
    !===Dummy variables to avoid warning
  END FUNCTION imposed_velocity_by_penalty

  !===Pressure for boundary conditions in Navier-Stokes.
  !===Can be used also to initialize pressure in: init_velocity_pressure_temperature
  !===Use this routine for outflow BCs only.
  !===CAUTION: Do not enfore BCs on pressure where normal component
  !            of velocity is prescribed.
  MODULE FUNCTION pp_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION pp_exact

  !===Temperature for boundary conditions in temperature equation.
  MODULE FUNCTION temperature_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('temperature_exact: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION temperature_exact

  !===Concentration for boundary conditions in concentration equation.
  MODULE FUNCTION concentration_exact(TYPE,rr,m,t) RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
     !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t
    !===Dummy variables to avoid warning
 END FUNCTION concentration_exact

  MODULE FUNCTION level_set_exact(interface_nb,TYPE,rr,m,t)  RESULT (vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m, interface_nb
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    CALL error_petsc('level_set_exact: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; nd=interface_nb; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION level_set_exact

  !===Penalty coefficient (if needed)
  !===This coefficient is equal to zero in subdomain
  !===where penalty is applied
  MODULE FUNCTION penal_in_real_space(mesh,rr_gauss,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN)   :: rr_gauss
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER,                    INTENT(IN)     :: nb_angles
    INTEGER,                    INTENT(IN)     :: nb, ne
    REAL(KIND=8),               INTENT(IN)     :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv

    vv = 1.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=mesh%np; nd=SIZE(rr_gauss,1); nd=SIZE(angles); nd=nb_angles; nd=nb; nd=ne; rd=time
    !===Dummy variables to avoid warning
  END FUNCTION penal_in_real_space

  !===Extension of the velocity field in the solid.
  !===Used only if Maxwell domain for H is larger than
  !===domain for Navier-Stokes velocity, or if there is
  !===no Navier-Stokes domain (i.e. mxw or mxx mode)
  MODULE FUNCTION extension_velocity(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_velocity

  MODULE FUNCTION extension_temperature(TYPE, H_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(H_Mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_temperature

  MODULE FUNCTION extension_concentration(TYPE, vv_mesh, mode, t, n_start) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: vv_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE, n_start
    INTEGER,                             INTENT(IN)   :: mode
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(vv_mesh%np)                :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=vv_mesh%np; nd=TYPE; nd=n_start; nd=mode; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION extension_concentration

  !===============================================================================
  !                       Boundary conditions for Maxwell
  !===============================================================================
  !===Velocity used in the induction equation.
  !===Used only in mxw mode
  MODULE FUNCTION Vexact(m, H_mesh) RESULT(vv)  !Set uniquement a l'induction
    IMPLICIT NONE
    TYPE(mesh_type),                       INTENT(IN) :: H_mesh
    INTEGER,                               INTENT(IN) :: m
    REAL(KIND=8), DIMENSION(H_mesh%np,6)              :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=m
    !===Dummy variables to avoid warning
  END FUNCTION Vexact

  MODULE FUNCTION H_B_quasi_static(char_h_b, rr, m) RESULT(vv)
    IMPLICIT NONE
    CHARACTER(LEN=1),                    INTENT(IN)   :: char_h_b
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6)             :: vv

    IF (inputs%if_quasi_static_approx) THEN
       vv = 0.d0
    ELSE
       CALL error_petsc('H_B_quasi_static should not be called')
    END IF
    RETURN
    !===Dummy variables to avoid warning
    cd1=char_h_b; nd=SIZE(rr,1); nd=m
    !===Dummy variables to avoid warning
  END FUNCTION H_B_quasi_static


  MODULE FUNCTION Hexact(H_mesh, TYPE, rr, m, mu_H_field, t) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type),                     INTENT(IN)   :: H_mesh
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: t
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: mu_H_field
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv
    INTEGER                                           :: n
    REAL(KIND=8), DIMENSION(5,SIZE(rr,2))             :: J
    REAL(KIND=8)                                      :: k_adim, k

    vv = 0.d0

    k_adim = 5.1356223d0
    k = k_adim/rmax
    CALL bessel_J(4,k*rr(1,:),J)

    IF (m==2 .AND. TYPE==1) THEN
       DO n = 1, SIZE(rr,2)
          vv(n) =  0.5d0*(A_plus*J(m+2,n)+A_minus*J(m,n))*EXP(k*rr(2,n))
       ENDDO
    END IF
    IF (m==2 .AND. TYPE==4) THEN
       DO n = 1, SIZE(rr,2)
          vv(n) =  0.5d0*(A_plus*J(m+2,n)-A_minus*J(m,n))*EXP(k*rr(2,n))
       ENDDO
    END IF
    IF (m==2 .AND. TYPE==5) THEN
       DO n = 1, SIZE(rr,2)
          vv(n) =  0.5d0*(A_minus-A_plus)*J(m+1,n)*EXP(k*rr(2,n))
       ENDDO
    END IF
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=TYPE; nd=SIZE(rr,1); nd=m; rd=t; nd=SIZE(mu_H_field)
    !===Dummy variables to avoid warning
    RETURN
  END FUNCTION Hexact

  !===Scalar potential for boundary conditions in the Maxwell equations.
  MODULE FUNCTION Phiexact(TYPE, rr, m, mu_phi,t) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, t
    REAL(KIND=8), DIMENSION(SIZE(rr,2))               :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION Phiexact

  !===Current in Ohm's law. Curl(H) = sigma(E + uxB) + current
  MODULE FUNCTION Jexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t, mesh_id, opt_B_ext) RESULT(vv)
    IMPLICIT NONE
    INTEGER     ,                        INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER     ,                        INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    INTEGER     ,                        INTENT(IN)   :: mesh_id
    REAL(KIND=8), DIMENSION(6), OPTIONAL,INTENT(IN)   :: opt_B_ext
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=sigma; rd=mu_H; rd=t; nd=mesh_id
    IF (PRESENT(opt_B_ext)) nd=SIZE(opt_B_ext)
    !===Dummy variables to avoid warning
  END FUNCTION Jexact_gauss

  !===Electric field for Neumann BC (cf. doc)
  MODULE FUNCTION Eexact_gauss(TYPE, rr, m, mu_phi, sigma, mu_H, t) RESULT(vv)
    IMPLICIT NONE
    INTEGER,                             INTENT(IN)   :: TYPE
    REAL(KIND=8), DIMENSION(:),          INTENT(IN)   :: rr
    INTEGER,                             INTENT(IN)   :: m
    REAL(KIND=8),                        INTENT(IN)   :: mu_phi, sigma, mu_H, t
    REAL(KIND=8)                                      :: vv

    vv = 0.d0
    RETURN
    !===Dummy variables to avoid warning
    nd=TYPE; nd=SIZE(rr,1); nd=m; rd=mu_phi; rd=sigma; rd=mu_H; rd=t
    !===Dummy variables to avoid warning
  END FUNCTION Eexact_gauss

  !===Initialization of magnetic field and scalar potential (if present)
  MODULE SUBROUTINE init_maxwell(H_mesh, phi_mesh, time, dt, mu_H_field, mu_phi, &
       list_mode, Hn1, Hn, phin1, phin)
    IMPLICIT NONE
    TYPE(mesh_type)                            :: H_mesh, phi_mesh
    REAL(KIND=8),                   INTENT(OUT):: time
    REAL(KIND=8),                   INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(:),     INTENT(IN) :: mu_H_field
    REAL(KIND=8),                   INTENT(IN) :: mu_phi
    INTEGER,      DIMENSION(:),     INTENT(IN) :: list_mode
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: Hn, Hn1
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT):: phin, phin1
    INTEGER                                    :: i, k

    time = -dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn1(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin1(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i) , mu_phi, time)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    time = time + dt
    DO k=1,6
       DO i=1, SIZE(list_mode)
          Hn(:,k,i) = Hexact(H_mesh,k, H_mesh%rr, list_mode(i), mu_H_field, time)
          IF (inputs%nb_dom_phi>0) THEN
             IF (k<3) THEN
                phin(:,k,i) = Phiexact(k, phi_mesh%rr, list_mode(i), mu_phi, time)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE init_maxwell

  !===Analytical permeability (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION mu_bar_in_fourier_space(H_mesh,nb,ne,pts,pts_ids) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(ne-nb+1)           :: vv
    INTEGER,     INTENT(IN)                    :: nb, ne
    REAL(KIND=8),DIMENSION(2,ne-nb+1),OPTIONAL :: pts
    INTEGER,     DIMENSION(ne-nb+1), OPTIONAL  :: pts_ids

    vv = 1.d0
    CALL error_petsc('mu_bar_in_fourier_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=nb; nd=ne
    IF (PRESENT(pts)) nd=SIZE(pts,1)
    IF (PRESENT(pts_ids)) nd=SIZE(pts_ids)
    !===Dummy variables to avoid warning
  END FUNCTION mu_bar_in_fourier_space

  !===Analytical mu_in_fourier_space (if needed)
  !===This function is not needed unless the flag
  !===     ===Use FEM Interpolation for magnetic permeability  (true/false)
  !===is activated and set to .FALSE. in the data data file. Default is .TRUE.
  MODULE FUNCTION grad_mu_bar_in_fourier_space(pt,pt_id) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8),DIMENSION(2), INTENT(in):: pt
    INTEGER,DIMENSION(1), INTENT(in)     :: pt_id
    REAL(KIND=8),DIMENSION(2)            :: vv

    vv = 0.d0
    CALL error_petsc('grad_mu_bar_in_fourier_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=SIZE(pt,1); nd=SIZE(pt_id)
    !===Dummy variables to avoid warning
  END FUNCTION grad_mu_bar_in_fourier_space

  !===Analytical permeability, mu in real space (if needed)
  MODULE FUNCTION mu_in_real_space(H_mesh,angles,nb_angles,nb,ne,time) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(:), INTENT(IN)     :: angles
    INTEGER, INTENT(IN)                        :: nb_angles
    INTEGER, INTENT(IN)                        :: nb, ne
    REAL(KIND=8), INTENT(IN)                   :: time
    REAL(KIND=8), DIMENSION(nb_angles,ne-nb+1) :: vv

    vv = 1.d0
    CALL error_petsc('mu_in_real_space: should not be called for this test')
    RETURN
    !===Dummy variables to avoid warning
    nd=H_mesh%np; nd=SIZE(angles); rd=time
    !===Dummy variables to avoid warning
  END FUNCTION mu_in_real_space

  !=== Only used for multiphase flow with variable electrical conductivity
  !=== This function is not needed unless the following line is in the data
  !=== ===Conductivity of fluid 0, fluid 1, ...
  !=== The above line should be followed by the value of the conductivity in the different fluids.
  !=== The result vv has to smaller than the effective conductivity at every nodes of the mesh.
  MODULE FUNCTION sigma_bar_in_fourier_space(H_mesh) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type), INTENT(IN)                :: H_mesh
    REAL(KIND=8), DIMENSION(SIZE(H_mesh%rr,2)) :: vv

    vv = 1.d0
    RETURN
  END FUNCTION sigma_bar_in_fourier_space

  !===Coefficient contaning the magnetic susceptibility for magnetic force in
  !ferrofluids:
  !===F = chi_coeff(T) * grad(H**2/2) (Kelvin form)
  !===or F = -H**2/2 * grad(chi_coeff(T)) (Helmholtz form)
  MODULE FUNCTION chi_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    RETURN
  END FUNCTION chi_coeff_law

  !===Coefficient contaning the temperature dependant factor in the
  !===pyromagnetic coefficient term of the temperature equation for ferrofluids:
  !===T * dchi/dT(T) * D/Dt(H**2/2)
  MODULE FUNCTION T_dchi_dT_coeff_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    RETURN
  END FUNCTION T_dchi_dT_coeff_law

  !===Kinematic viscosity's law of temperature
  MODULE FUNCTION nu_tilde_law(temp) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: temp
    REAL(KIND=8) :: vv

    vv = 0.d0*temp
    RETURN
  END FUNCTION nu_tilde_law

  MODULE FUNCTION rot_H_jump_interface(mesh,rr,list_mode) RESULT(vv)
    IMPLICIT NONE
    TYPE(mesh_type)                                   :: mesh
    REAL(KIND=8), DIMENSION(:,:),        INTENT(IN)   :: rr
    INTEGER,      DIMENSION(:),          INTENT(IN)   :: list_mode
    REAL(KIND=8), DIMENSION(SIZE(rr,2),6,SIZE(list_mode)) :: vv
    INTEGER                                           :: n, i, m
    REAL(KIND=8), DIMENSION(5,SIZE(rr,2))             :: J
    REAL(KIND=8) :: k_adim, k
    REAL(KIND=8) :: sigma1, sigma2, xx

    vv = 0.d0

    sigma1 = inputs%sigma(1)
    sigma2 = inputs%sigma(2)

    xx = (1.d0/sigma1-1.d0/sigma2)

    k_adim = 5.1356223d0
    k = k_adim/rmax
    CALL bessel_J(4,k*rr(1,:),J)

    DO i = 1, SIZE(list_mode)
       m = list_mode(i)

       IF (m==2) THEN
          DO n = 1, SIZE(rr,2)
             vv(n,1,i) = 0.d0
             vv(n,2,i) = 0.5d0*(J(m,n)-J(m+2,n))*EXP(k*rr(2,n))
             vv(n,3,i) = 0.5d0*(J(m+2,n)+J(m,n))*EXP(k*rr(2,n))
             vv(n,4,i) = 0.d0
             vv(n,5,i) = 0.d0
             vv(n,6,i) = J(m+1,n)*EXP(k*rr(2,n))
          ENDDO
          vv = xx*k/2.d0*(A_minus+A_plus)*vv
       END IF
    ENDDO
    RETURN
    !===Dummy variables to avoid warning
    nd=mesh%np
    !===Dummy variables to avoid warning
  END FUNCTION rot_H_jump_interface

  MODULE FUNCTION Derivative_of_potential_from_rhoLi(delta_rhoLi_phys) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: delta_rhoLi_phys
    REAL(KIND=8) :: vv

    vv = 0.d0*delta_rhoLi_phys
    RETURN
  END FUNCTION Derivative_of_potential_from_rhoLi

  MODULE FUNCTION molar_fraction_from_concentration(delta_rhoLi_phys) RESULT(vv)
    IMPLICIT NONE
    REAL(KIND=8) :: delta_rhoLi_phys
    REAL(KIND=8) :: vv

    vv = 0.d0*delta_rhoLi_phys
    RETURN
  END FUNCTION molar_fraction_from_concentration

  !====Local functions
  SUBROUTINE bessel_J(N_in,rr,J_out)
    IMPLICIT NONE
    INTEGER     ,                             INTENT(IN)  :: N_in
    REAL(KIND=8), DIMENSION(:),               INTENT(IN)  :: rr
    REAL(KIND=8), DIMENSION(:,:),             INTENT(OUT) :: J_out
    INTEGER                                               :: n, i

    DO n=1, SIZE(rr)
       J_out(1,n) = BESSJ0(rr(n))
       IF (N_in .EQ. 0) THEN
          CYCLE
       ELSE IF (N_in .EQ. 1) THEN
          J_out(2,n) = BESSJ1(rr(n))
       ELSE
          J_out(2,n) = BESSJ1(rr(n))
          DO i=2,N_in
             IF (rr(n) .LT. 1.d-12) THEN
                J_out(i+1,n) = 0.d0
             ELSE
                J_out(i+1,n) = 2*(i-1)*J_out(i,n)/(rr(n)+1.d-14) - J_out(i-1,n)
             END IF
          END DO
       ENDIF
    END DO
  END SUBROUTINE bessel_J

END SUBMODULE boundary_generic
