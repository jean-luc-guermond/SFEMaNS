MODULE projection_module
    USE boundary, ONLY:  extension_velocity, extension_concentration, extension_temperature
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    PUBLIC :: init_projection_connectivity_structures, &
              projection_velocity, &
              projection_concentration, &
              projection_mag_field, &
              projection_temperature

  !Connectivity structures----------------------------------------------------
  !October 7, 2008, JLG
    INTEGER, DIMENSION(:), ALLOCATABLE, TARGET, PUBLIC :: jj_v_to_H, jj_v_to_T, jj_T_to_H, jj_c_to_v, jj_c_to_H

CONTAINS

    SUBROUTINE init_projection_connectivity_structures(H_mesh, vv_mesh, temp_mesh, conc_mesh)
        USE input_data
        USE def_type_mesh
        USE SFEMaNS_logicals_module
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN) :: H_mesh, vv_mesh, temp_mesh, conc_mesh
        INTEGER :: m

    !===Create data structure jj_v_to_H=============================================
        IF (if_induction) THEN
            ALLOCATE(jj_v_to_H(H_mesh%np), SOURCE=-1)
            IF (if_momentum .OR. inputs%type_pb=='mxx') THEN
                DO m = 1, vv_mesh%me
                    jj_v_to_H(H_mesh%jj(:, m)) = vv_mesh%jj(1:H_mesh%gauss%n_w, m)
                END DO
            END IF
        END IF

    !===Create data structure jj_v_to_T==========================================
        IF (if_energy) THEN
            ALLOCATE(jj_v_to_T(temp_mesh%np), SOURCE=-1)
            IF (vv_mesh%me/=0) THEN
                DO m = 1, vv_mesh%me
                    jj_v_to_T(temp_mesh%jj(:, m)) = vv_mesh%jj(1:temp_mesh%gauss%n_w, m)
                END DO
            END IF
        END IF

    !===Create data structure jj_c_to_v==========================================
        IF (if_concentration) THEN
            ALLOCATE(jj_c_to_v(vv_mesh%np), SOURCE=-1)
            IF (conc_mesh%me/=0) THEN
                DO m = 1, conc_mesh%me
                    jj_c_to_v(vv_mesh%jj(:, m)) = conc_mesh%jj(1:vv_mesh%gauss%n_w, m)
                END DO
            END IF
        END IF

    !===Create data structure jj_T_to_H==========================================
        IF (if_induction) THEN
            ALLOCATE(jj_T_to_H(H_mesh%np), SOURCE=-1)
            IF (if_energy) THEN
                DO m = 1, temp_mesh%me
                    jj_T_to_H(H_mesh%jj(:, m)) = temp_mesh%jj(1:H_mesh%gauss%n_w, m)
                END DO
            END IF
        END IF

    !===Create data structure jj_c_to_H==========================================
        IF (if_induction) THEN
            ALLOCATE(jj_c_to_H(H_mesh%np), SOURCE=-1)
            IF (if_concentration) THEN
                DO m = 1, conc_mesh%me
                    jj_c_to_H(H_mesh%jj(:, m)) = conc_mesh%jj(1:H_mesh%gauss%n_w, m)
                END DO
            END IF
        END IF

    END SUBROUTINE init_projection_connectivity_structures

    !---------------------------------------------------------------------------
    SUBROUTINE projection_concentration(mesh, vn, which_projection, coupling_variable, list_mode, time)
        USE boundary
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN)                               :: mesh
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
        CHARACTER(LEN=*), INTENT(IN)                              :: which_projection
        INTEGER, DIMENSION(:), POINTER                            :: connectivity_structure
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
        INTEGER                                                   :: i,j,k
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8),        INTENT(IN) :: time

        !=== Point to correct connectivity array ===!
        SELECT CASE(which_projection)
        CASE('c_to_v')
            connectivity_structure => jj_c_to_v
        CASE('c_to_H')
            connectivity_structure => jj_c_to_H
        CASE DEFAULT
            CALL error_petsc("BUG in projection_concentration: wrong projection "//which_projection//", should be c_to_v or c_to_H.")
        END SELECT

        !=== Concentration smallest mesh ==> only extension is possible ===!
        IF (mesh%np>0) THEN
            DO i = 1, SIZE(list_mode)
                DO k = 1, 2 !===The user has to code extension_concentration
                    coupling_variable(:,k,i) = extension_concentration(k, mesh, list_mode(i),time,1)
                END DO
            END DO
        END IF
        IF (mesh%me /=0) THEN !===construction of the extended field
            DO j = 1, SIZE(connectivity_structure)
                IF (connectivity_structure(j) == -1) CYCLE
                coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
            END DO
        END IF
    END SUBROUTINE projection_concentration
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    SUBROUTINE projection_velocity(mesh, vn, which_projection, if_restriction, coupling_variable, list_mode, time, if_opt_zero_extension)
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN)                               :: mesh
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
        CHARACTER(LEN=*), INTENT(IN)                              :: which_projection
        INTEGER, DIMENSION(:), POINTER                            :: connectivity_structure
        LOGICAL, INTENT(IN)                                       :: if_restriction
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
        INTEGER                                                   :: i, j, k
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8),        INTENT(IN) :: time
        LOGICAL, OPTIONAL, INTENT(IN) :: if_opt_zero_extension
        LOGICAL :: if_zero_extension

        !=== Point to correct connectivity array ===!
        SELECT CASE(which_projection)
        CASE('c_to_v')
            connectivity_structure => jj_c_to_v
        CASE('v_to_T')
            connectivity_structure => jj_v_to_T
        CASE('v_to_H')
            connectivity_structure => jj_v_to_H
        CASE DEFAULT
            CALL error_petsc("BUG in projection_velocity: wrong projection "//which_projection//", should be c_to_v, v_to_T or v_to_H.")
        END SELECT

        !=== If restriction ===!
        IF (if_restriction) THEN
            IF (mesh%me /=0) THEN !===construction of the restricted velocity field
                DO j = 1, SIZE(connectivity_structure)
                    IF (connectivity_structure(j) == -1) CYCLE
                    coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
                END DO
            END IF
        !=== If extension ===!
        ELSE
            IF (mesh%np>0) THEN !===Check that ns is a subset of temp before constructing the extended vv field
                IF (PRESENT(if_opt_zero_extension)) THEN
                    if_zero_extension = if_opt_zero_extension
                ELSE
                    if_zero_extension = .FALSE.
                ENDIF
                coupling_variable = 0.d0
                IF (.NOT. if_zero_extension) THEN
                    DO i = 1, SIZE(list_mode)
                        DO k = 1, 6 !===The user has to code extension_vel
                            coupling_variable(:,k,i) = extension_velocity(k, mesh, list_mode(i), time, 1)
                        END DO
                    END DO
                ENDIF
            END IF
            IF (mesh%me /=0) THEN !===construction of the extended field
                DO j = 1, SIZE(connectivity_structure)
                    IF (connectivity_structure(j) == -1) CYCLE
                    coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
                END DO
            END IF
        END IF
    END SUBROUTINE projection_velocity
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    SUBROUTINE projection_temperature(mesh, vn, which_projection, if_restriction, coupling_variable, list_mode, time)
        USE boundary
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN)                               :: mesh
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
        CHARACTER(LEN=*), INTENT(IN)                              :: which_projection
        INTEGER, DIMENSION(:), POINTER                            :: connectivity_structure
        LOGICAL, INTENT(IN)                                       :: if_restriction
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
        INTEGER                                                   :: i,j,k
        INTEGER, DIMENSION(:), INTENT(IN) :: list_mode
        REAL(KIND = 8),        INTENT(IN) :: time

        !=== Point to correct connectivity array ===!
        SELECT CASE(which_projection)
        CASE('v_to_T')
            connectivity_structure => jj_v_to_T
        CASE('T_to_H')
            connectivity_structure => jj_T_to_H
        CASE DEFAULT
            CALL error_petsc("BUG in projection_temperature: wrong projection "//which_projection//", should be v_to_T or T_to_H.")
        END SELECT

        !=== if restriction ===!
        IF (if_restriction) THEN
            IF (mesh%me /=0) THEN !===construction of the restricted temperature field
                DO j = 1, SIZE(connectivity_structure)
                    IF (connectivity_structure(j) == -1) CYCLE
                    coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
                END DO
            END IF
        !=== if extension ===!
        ELSE
            IF (mesh%np>0) THEN !===Check that temp is a subset of H before constructing the extended vv field
                DO i = 1, SIZE(list_mode)
                    DO k = 1, 2 !===The user has to code extension_temperature
                        coupling_variable(:,k,i) = extension_temperature(k, mesh, list_mode(i),time,1)
                    END DO
                END DO
            END IF
            IF (mesh%me /=0) THEN !===construction of the extended field
                DO j = 1, SIZE(connectivity_structure)
                    IF (connectivity_structure(j) == -1) CYCLE
                    coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
                END DO
            END IF
        END IF

    END SUBROUTINE projection_temperature
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    SUBROUTINE projection_mag_field(mesh, H_mesh, vn, which_projection, if_restriction, coupling_variable)
        USE my_util
        IMPLICIT NONE
        TYPE(mesh_type), INTENT(IN)                               :: mesh, H_mesh
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN)                :: vn
        CHARACTER(LEN=*), INTENT(IN)                              :: which_projection
        INTEGER, DIMENSION(:), POINTER                            :: connectivity_structure
        LOGICAL, INTENT(IN)                                       :: if_restriction
        REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT)               :: coupling_variable
        INTEGER                                                   :: j, m
        
        !=== Point to correct connectivity array ===!
        SELECT CASE(which_projection)
        CASE('c_to_H')
            connectivity_structure => jj_c_to_H
        CASE('v_to_H')
            connectivity_structure => jj_v_to_H
        CASE('v_to_T')
            connectivity_structure => jj_v_to_T
        CASE('T_to_H')
            connectivity_structure => jj_T_to_H
        CASE DEFAULT
            CALL error_petsc("BUG in projection_mag_field: wrong projection "//which_projection//", should be c_to_H, v_to_T, v_to_H or T_to_H.")
        END SELECT

        DO j = 1, SIZE(connectivity_structure)
            IF (connectivity_structure(j) == -1) CYCLE
            IF (if_restriction) THEN
                coupling_variable(connectivity_structure(j),:,:) = vn(j,:,:)
            ELSE
                coupling_variable(j,:,:) = vn(connectivity_structure(j),:,:)
            END IF
        END DO

        !=== If restriction from P2 to P1 ===!
        IF (if_restriction) THEN
            IF (H_mesh%gauss%n_w/=mesh%gauss%n_w) THEN
                DO m = 1, mesh%me
                    coupling_variable(mesh%jj(4,m),:,:) = (coupling_variable(mesh%jj(2,m),:,:) &
                        + coupling_variable(mesh%jj(3,m),:,:)) / 2
                    coupling_variable(mesh%jj(5,m),:,:) = (coupling_variable(mesh%jj(3,m),:,:) &
                        + coupling_variable(mesh%jj(1,m),:,:)) / 2
                    coupling_variable(mesh%jj(6,m),:,:) = (coupling_variable(mesh%jj(1,m),:,:) &
                        + coupling_variable(mesh%jj(2,m),:,:)) / 2
                END DO
            END IF
        END IF

    END SUBROUTINE projection_mag_field
    !---------------------------------------------------------------------------

    
END MODULE projection_module