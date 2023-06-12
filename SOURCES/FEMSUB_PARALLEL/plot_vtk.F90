MODULE vtk_viz
#include "petsc/finclude/petsc.h"
  USE petsc
  PUBLIC ::  create_pvd_file, check_list,&
       make_vtu_file_arpack
  PUBLIC :: create_vtu_file_axi3D
  PUBLIC :: make_vtu_file_3D, create_xml_vtu_scal_file, create_xml_vtu_vect_file
  !PUBLIC :: make_vtu_file_scalar_2D
  !PUBLIC :: make_vtu_file_axi3D, create_vtk_file
  PRIVATE
CONTAINS

  SUBROUTINE check_list(communicator, file_list, check)
    IMPLICIT NONE
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: dummy_list
    INTEGER, DIMENSION(SIZE(file_list)) :: check_mylist
    INTEGER                             :: check, n, count
    !#include "petsc/finclude/petsc.h"
    MPI_Comm                            :: communicator
    PetscMPIInt                         :: rank, nb_procs
    PetscErrorCode                      :: ierr
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)

    CALL MPI_ALLGATHER(check, 1, MPI_INTEGER, check_mylist, 1, &
         MPI_INTEGER, communicator, ierr)

    count = 0
    DO n = 1, SIZE(file_list)
       IF (check_mylist(n)==0) CYCLE
       count = count + 1
    END DO
    ALLOCATE(dummy_list(count))
    count = 0
    DO n = 1, SIZE(file_list)
       IF (check_mylist(n)==0) CYCLE
       count = count + 1
       dummy_list(count) = file_list(n)
    END DO
    DEALLOCATE(file_list)
    ALLOCATE(file_list(count))
    file_list = dummy_list
  END SUBROUTINE check_list

  SUBROUTINE create_pvd_file(file_list, file_header, time_step, what)
    IMPLICIT NONE
    CHARACTER(*), DIMENSION(:), INTENT(IN) :: file_list
    CHARACTER(*),               INTENT(IN) :: file_header, what
    INTEGER,                    INTENT(IN) :: time_step
    INTEGER                                :: unit_file=789, j
    CHARACTER(len=5)                       :: tit, tit_part

    IF (what=='new') THEN
       OPEN (UNIT=unit_file, FILE=file_header//'.pvd', FORM = 'formatted', &
            ACCESS = 'append', STATUS = 'replace')
       WRITE(unit_file, '(A)') '<?xml version="1.0"?>'
       WRITE(unit_file, '(A)') '<VTKFile type="Collection" version="0.1"'// &
            ' byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
       WRITE(unit_file, '(A)') '<Collection>'
    ELSE
       OPEN (UNIT=unit_file, FILE=file_header//'.pvd', FORM = 'formatted', &
            ACCESS = 'append', STATUS = 'old')
       BACKSPACE(unit_file)
       BACKSPACE(unit_file)
    END IF
    WRITE(tit,'(I5)') time_step
    DO j = 1, SIZE(file_list)
       WRITE(tit_part,'(I5)') j
       WRITE(unit_file,'(A)') '<DataSet timestep="'//TRIM(ADJUSTL(tit))//&
            '" group="" part="'// &
            TRIM(ADJUSTL(tit_part))//'" file="./'//TRIM(ADJUSTL(file_list(j)))//&
            '.vtu'//'"/>'
    END DO
    WRITE(unit_file, '(A)') '</Collection>'
    WRITE(unit_file, '(A)') '</VTKFile>'

    CLOSE(unit_file) !Ecriture pour paraview
  END SUBROUTINE create_pvd_file

  SUBROUTINE create_xml_vtu_scal_file(field, mesh, file_name, field_name)
    USE def_type_mesh
    USE input_data
    USE zlib_base64
    IMPLICIT NONE
    TYPE(mesh_type)                        :: mesh
    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: field
    CHARACTER(*),               INTENT(IN) :: file_name, field_name
    INTEGER                                :: unit_file=789, m, n, type_cell
    REAL(KIND=4),    DIMENSION(3*mesh%np)              :: r4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%gauss%n_w*mesh%me) :: i4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%me)                :: i4_xml_field
    INTEGER(KIND=1), DIMENSION(mesh%me)                :: i1_xml_field
    CHARACTER(LEN=200)                         :: ascii_or_binary

    IF (SIZE(field)==0) RETURN

    IF (inputs%if_xml) THEN
       ascii_or_binary = 'binary'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<?xml version="1.0" ?>'
       WRITE(unit_file,'(A)', advance="no") '<VTKFile type="UnstructuredGrid" version="0.1" '
       WRITE(unit_file,'(A)') 'compressor="vtkZLibDataCompressor" byte_order="LittleEndian">'
    ELSE
       ascii_or_binary = 'ascii'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
            FORM = 'formatted', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
            ' byte_order="LittleEndian">'
    END IF

    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
         '" NumberOfCells="', mesh%me, '">'

    !===PointData Block
    WRITE(unit_file,'(A)') '<PointData Scalars="truc">'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'&
         //TRIM(ADJUSTL(field_name))//'" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       CALL write_compressed_block(unit_file,  REAL(field(1:mesh%np),4))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(e14.7)') field(n)
       ENDDO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</PointData>'
    !===End of PointData Block

    !===CellData Block
    WRITE(unit_file,'(A)') '<CellData>'
    WRITE(unit_file,'(A)') '</CellData>'
    !===End of CellData Block

    !===Points Block
    WRITE(unit_file,'(A)') '<Points>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
         'NumberOfComponents="3" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1) = REAL(mesh%rr(1,n),4)
          r4_threed_xml_field(3*(n-1)+2) = REAL(0.,4)
          r4_threed_xml_field(3*(n-1)+3) = REAL(mesh%rr(2,n),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') mesh%rr(1,n), 0.d0 , mesh%rr(3,n)
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Points>'
    !===End of Points Block

    !===Cells Block
    IF (mesh%gauss%n_w==3) THEN
       type_cell = 5
    ELSE IF (mesh%gauss%n_w==6) THEN
       type_cell = 22
    END IF
    WRITE(unit_file,'(A)') '<Cells>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="connectivity" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+1) = INT(mesh%jj(1,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+2) = INT(mesh%jj(2,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+3) = INT(mesh%jj(3,m)-1,4)
          IF (type_cell==22) THEN
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+4) = INT(mesh%jj(6,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+5) = INT(mesh%jj(4,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+6) = INT(mesh%jj(5,m)-1,4)
          END IF
       END DO
       CALL write_compressed_block(unit_file, i4_threed_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1
          IF (type_cell==22) THEN
             WRITE(unit_file,'(3(I8,1x))') mesh%jj(6,m)-1 , mesh%jj(4,m)-1 , mesh%jj(5,m)-1
          END IF
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'

    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="offsets" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_xml_field(m) = INT(mesh%gauss%n_w*m,4)
       END DO
       CALL write_compressed_block(unit_file, i4_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') m*mesh%gauss%n_w
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'

    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i1_xml_field(m) = INT(type_cell,1)
       END DO
       CALL write_compressed_block(unit_file, i1_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') type_cell
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Cells>'
    !===End of Cells Block

    WRITE(unit_file,'(A)') '</Piece>'
    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
    WRITE(unit_file,'(A)') '</VTKFile>'

    CLOSE(unit_file)
  END SUBROUTINE create_xml_vtu_scal_file

  SUBROUTINE create_xml_vtu_vect_file(field, mesh, file_name, opt_st)
    USE def_type_mesh
    USE input_data
    USE zlib_base64
    IMPLICIT NONE
    TYPE(mesh_type)                          :: mesh
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: field
    CHARACTER(*),                 INTENT(IN) :: file_name
    INTEGER                                  :: unit_file=789, m, n, type_cell
    CHARACTER(*), OPTIONAL,       INTENT(IN) :: opt_st
    CHARACTER(LEN=200)                       :: field_name
    REAL(KIND=4),    DIMENSION(3*mesh%np)              :: r4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%gauss%n_w*mesh%me) :: i4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%me)                :: i4_xml_field
    INTEGER(KIND=1), DIMENSION(mesh%me)                :: i1_xml_field
    CHARACTER(LEN=200)                         :: ascii_or_binary

    IF (PRESENT(opt_st)) THEN
       field_name=opt_st
    ELSE
       field_name='field'
    END IF

    IF (SIZE(field)==0) RETURN

    IF (inputs%if_xml) THEN
       ascii_or_binary = 'binary'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<?xml version="1.0" ?>'
       WRITE(unit_file,'(A)', advance="no") '<VTKFile type="UnstructuredGrid" version="0.1" '
       WRITE(unit_file,'(A)') 'compressor="vtkZLibDataCompressor" byte_order="LittleEndian">'
    ELSE
       ascii_or_binary = 'ascii'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
            FORM = 'formatted', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
            ' byte_order="LittleEndian">'
    END IF

    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
         '" NumberOfCells="', mesh%me, '">'

    !===PointData Block
    WRITE(unit_file,'(A)') '<PointData>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
         '_cos" format="'//TRIM(ADJUSTL(ascii_or_binary))//'" NumberOfComponents="3">'
    IF (inputs%if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1:3*n) = REAL(field(n,1:5:2),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') field(n,1), field(n,3), field(n,5) !===JLG+NORE Jan 18, 2022
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
         '_sin" format="'//TRIM(ADJUSTL(ascii_or_binary))//'" NumberOfComponents="3">'
    IF (inputs%if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1:3*n) = REAL(field(n,2:6:2),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') field(n,2), field(n,4), field(n,6)
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</PointData>'
    !===End of PointData Block

    !===CellData Block
    WRITE(unit_file,'(A)') '<CellData>'
    WRITE(unit_file,'(A)') '</CellData>'
    !===End of CellData Block

    !===Points Block
    WRITE(unit_file,'(A)') '<Points>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
         'NumberOfComponents="3" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1) = REAL(mesh%rr(1,n),4)
          r4_threed_xml_field(3*(n-1)+2) = REAL(0.,4)
          r4_threed_xml_field(3*(n-1)+3) = REAL(mesh%rr(2,n),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') mesh%rr(1,n), 0.d0 , mesh%rr(3,n)
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Points>'
    !===End of Points Block

    !===Cells Block
    IF (mesh%gauss%n_w==3) THEN
       type_cell = 5
    ELSE IF (mesh%gauss%n_w==6) THEN
       type_cell = 22
    END IF
    WRITE(unit_file,'(A)') '<Cells>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="connectivity" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+1) = INT(mesh%jj(1,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+2) = INT(mesh%jj(2,m)-1,4)
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+3) = INT(mesh%jj(3,m)-1,4)
          IF (type_cell==22) THEN
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+4) = INT(mesh%jj(6,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+5) = INT(mesh%jj(4,m)-1,4)
             i4_threed_xml_field(mesh%gauss%n_w*(m-1)+6) = INT(mesh%jj(5,m)-1,4)
          END IF
       END DO
       CALL write_compressed_block(unit_file, i4_threed_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1
          IF (type_cell==22) THEN
             WRITE(unit_file,'(3(I8,1x))') mesh%jj(6,m)-1 , mesh%jj(4,m)-1 , mesh%jj(5,m)-1
          END IF
       END DO
    END IF

    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="offsets" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_xml_field(m) = INT(mesh%gauss%n_w*m,4)
       END DO
       CALL write_compressed_block(unit_file, i4_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') m*mesh%gauss%n_w
       END DO
    END IF

    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i1_xml_field(m) = INT(type_cell,1)
       END DO
       CALL write_compressed_block(unit_file, i1_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') type_cell
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Cells>'
    ! End of Cells Block ---------------------------------------------------------

    WRITE(unit_file,'(A)') '</Piece>'
    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
    WRITE(unit_file,'(A)') '</VTKFile>'

    CLOSE(unit_file)
  END SUBROUTINE create_xml_vtu_vect_file



  SUBROUTINE make_vtu_file_arpack(communicator, mesh, header, field, field_name, what, num_vp)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),               INTENT(IN) :: mesh
    CHARACTER(*),                  INTENT(IN) :: header
    CHARACTER(*),                  INTENT(IN) :: field_name, what
    INTEGER,                       INTENT(IN) :: num_vp
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN) :: field
    CHARACTER(LEN=200), DIMENSION(1)          :: file_list
    CHARACTER(LEN=3)                          :: st_rank
    !#include "petsc/finclude/petsc.h"
    PetscErrorCode                            :: ierr
    PetscMPIInt                               :: rank, nb_procs
    MPI_Comm                                  :: communicator
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)

    WRITE(st_rank,'(I3)') num_vp
    file_list(1) = TRIM(header)//'_eigen_'//TRIM(ADJUSTL(st_rank))

    CALL create_pvd_file(file_list, TRIM(header), num_vp, TRIM(what))

    !=========TEST FL Feb. 11th, 2013
    IF (SIZE(field,2) == 6) THEN
       !CALL create_vtu_vect_file(field, mesh, TRIM(ADJUSTL(file_list(1))), field_name)
       CALL create_xml_vtu_vect_file(field, mesh, TRIM(ADJUSTL(file_list(1))), field_name)
    ELSE IF (SIZE(field,2) == 1) THEN
       CALL create_xml_vtu_scal_file(field(:,1), mesh, TRIM(ADJUSTL(file_list(1))), field_name)
    ELSE IF (SIZE(field,2) .GT. 0) THEN
       CALL create_xml_vtu_scal_file(field(:,1), mesh, TRIM(ADJUSTL(file_list(1))), field_name)
    ELSE
       CALL error_Petsc('Bug in make_vtu_file_arpack: field needs at least one component')
    END IF
    !=========TEST FL Feb. 11th, 2013
  END SUBROUTINE make_vtu_file_arpack

  SUBROUTINE create_vtu_file_axi3D(field, mesh, file_name, opt_st)
    USE def_type_mesh
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:,:), INTENT(IN) :: field
    CHARACTER(*),                   INTENT(IN) :: file_name
    CHARACTER(*), OPTIONAL,         INTENT(IN) :: opt_st
    CHARACTER(LEN=200)                         :: field_name
    REAL(KIND=8)                               :: theta, pi, dtheta
    INTEGER                                    :: unit_file=789, m, i, type_cell, &
         nb_angle, k

    IF (SIZE(field,2)==0) RETURN

    IF (PRESENT(opt_st)) THEN
       field_name=opt_st
    ELSE
       field_name='field'
    END IF

    nb_angle = SIZE(field,1)
    pi = ACOS(-1.d0)
    dtheta = 2*pi/nb_angle

    OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
         FORM = 'formatted', STATUS = 'unknown')

    WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
         ' byte_order="LittleEndian">'
    WRITE(unit_file,'(A)') '<UnstructuredGrid>'

    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', nb_angle*mesh%np, &
         '" NumberOfCells="', nb_angle*mesh%me, '">'
    ! PointData Block ------------------------------------------------------------
    WRITE(unit_file,'(A)') '<PointData Scalars="truc">'

    IF (SIZE(field,2)==1) THEN
       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
            '" format="ascii">'
       DO k = 1, nb_angle
          DO i = 1, mesh%np
             WRITE(unit_file,'(e14.7)') field(k, 1, i)
          ENDDO
       END DO
    ELSE
       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
            '" format="ascii" NumberOfComponents="3">'
       DO k = 1, nb_angle
          DO i = 1, mesh%np
             WRITE(unit_file,'(3(e14.7,x))') field(k, 1, i), field(k, 2, i), field(k, 3, i)
          ENDDO
       END DO
    END IF

    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</PointData>'
    ! End of PointData Block -----------------------------------------------------

    ! CellData Block -------------------------------------------------------------
    WRITE(unit_file,'(A)') '<CellData>'
    WRITE(unit_file,'(A)') '</CellData>'
    ! End of CellData Block ------------------------------------------------------

    ! Points Block ---------------------------------------------------------------
    WRITE(unit_file,'(A)') '<Points>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
         'NumberOfComponents="3" format="ascii">'
    DO k = 1, nb_angle
       theta = (k-1)*dtheta
       DO i=1, mesh%np
          WRITE(unit_file,'(3(e14.7,2x))') mesh%rr(1,i)*COS(theta), &
               mesh%rr(1,i)*SIN(theta), mesh%rr(2,i)
       ENDDO
    END DO
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Points>'
    ! End of Points Block --------------------------------------------------------

    ! Cells Block ----------------------------------------------------------------
    type_cell = 13
    WRITE(unit_file,'(A)') '<Cells>'
    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="connectivity" format="ascii">'
    DO k = 1, nb_angle-1
       DO m = 1, mesh%me
          WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1+(k-1)*mesh%np
          WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1+k*mesh%np
       END DO
    END DO
    k = nb_angle
    DO m = 1, mesh%me
       WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1+(k-1)*mesh%np
       WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1
    END DO
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="offsets" format="ascii">'
    DO m = 1, nb_angle*mesh%me
       WRITE(unit_file,'(I8)') 6*m
    END DO
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
    DO m = 1, nb_angle*mesh%me
       WRITE(unit_file,'(I8)') type_cell
    END DO
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Cells>'
    ! End of Cells Block ---------------------------------------------------------

    WRITE(unit_file,'(A)') '</Piece>'
    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
    WRITE(unit_file,'(A)') '</VTKFile>'

    CLOSE(unit_file)
  END SUBROUTINE create_vtu_file_axi3D

  SUBROUTINE make_vtu_file_3D(communicator, mesh, header, &
       field, field_name, what, opt_it)
    USE def_type_mesh
    USE my_util
    IMPLICIT NONE
    TYPE(mesh_type),               INTENT(IN) :: mesh
    CHARACTER(*),                  INTENT(IN) :: header
    CHARACTER(*),                  INTENT(IN) :: field_name, what
    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
    REAL(KIND=8), DIMENSION(:,:),  INTENT(IN) :: field
    INTEGER                                   :: j, it
    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
    CHARACTER(LEN=3)                          :: st_rank, st_it
    !#include "petsc/finclude/petsc.h"
    PetscErrorCode                            :: ierr
    PetscMPIInt                               :: rank, nb_procs
    MPI_Comm                                  :: communicator
    CALL MPI_Comm_rank(communicator, rank, ierr)
    CALL MPI_Comm_size(communicator, nb_procs, ierr)
    ALLOCATE(file_list(nb_procs))

    IF (PRESENT(opt_it)) THEN
       it = opt_it
       WRITE(st_it,'(I3)') it
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))//&
               '_it_'//TRIM(ADJUSTL(st_it))
       END DO
    ELSE
       DO j = 1, nb_procs
          WRITE(st_rank,'(I3)') j
          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))
       END DO
    END IF

    CALL check_list(communicator, file_list, mesh%np)
    IF (rank==0) THEN
       IF (PRESENT(opt_it)) THEN
          it = opt_it
       ELSE
          it = 1

       END IF
       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
    END IF

    !CALL create_vtu_file_3D(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), &
    !     opt_st=field_name)
    CALL create_xml_vtu_file_3D(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), &
         opt_st=field_name)

  END SUBROUTINE make_vtu_file_3D



  SUBROUTINE create_xml_vtu_file_3D(field, mesh, file_name, opt_st)
    USE def_type_mesh
    USE my_util
    USE input_data
    USE zlib_base64
    IMPLICIT NONE
    TYPE(mesh_type)                            :: mesh
    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: field
    CHARACTER(*),                   INTENT(IN) :: file_name
    CHARACTER(*), OPTIONAL,         INTENT(IN) :: opt_st
    CHARACTER(LEN=200)                         :: field_name
    INTEGER                                    :: unit_file=789, m, n, type_cell, stride
    REAL(KIND=4),    DIMENSION(3*mesh%np)              :: r4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%gauss%n_w*mesh%me) :: i4_threed_xml_field
    INTEGER(KIND=4), DIMENSION(mesh%me)                :: i4_xml_field
    INTEGER(KIND=1), DIMENSION(mesh%me)                :: i1_xml_field
    CHARACTER(LEN=200)                         :: ascii_or_binary

    IF (SIZE(field,2)==0) RETURN

    IF (SIZE(mesh%jj,1)==6) THEN
       type_cell = 13
       stride = 6
    ELSE IF (SIZE(mesh%jj,1)==15) THEN
       type_cell = 26
       stride = 15
    ELSE
       type_cell = 0
       stride = 0
       CALL error_petsc('Bug in create_vtu_file_3D: SIZE(mesh%jj,1) is wrong')
    END IF

    IF (PRESENT(opt_st)) THEN
       field_name=opt_st
    ELSE
       field_name='field'
    END IF

    IF (inputs%if_xml) THEN
       ascii_or_binary = 'binary'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<?xml version="1.0" ?>'
       WRITE(unit_file,'(A)', advance="no") '<VTKFile type="UnstructuredGrid" version="0.1" '
       WRITE(unit_file,'(A)') 'compressor="vtkZLibDataCompressor" byte_order="LittleEndian">'
    ELSE
       ascii_or_binary = 'ascii'
       OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
            FORM = 'formatted', STATUS = 'unknown')
       WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
            ' byte_order="LittleEndian">'
    END IF

    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
         '" NumberOfCells="', mesh%me, '">'

    !===PointData Block
    WRITE(unit_file,'(A)') '<PointData Scalars="truc">'
    IF (SIZE(field,1)==1) THEN
       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
            '" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
       IF (inputs%if_xml) THEN
          CALL write_compressed_block(unit_file, REAL(field(1,1:mesh%np),4))
       ELSE
          DO n = 1, mesh%np
             WRITE(unit_file,'(e14.7)') field(1,n)
          END DO
       END IF
    ELSE
       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
            '" format="'//TRIM(ADJUSTL(ascii_or_binary))//'" NumberOfComponents="3">'
       IF (inputs%if_xml) THEN
          DO n = 1, mesh%np
             r4_threed_xml_field(3*(n-1)+1:3*n) = REAL(field(1:3,n),4)
          END DO
          CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
       ELSE
          DO n = 1, mesh%np
             WRITE(unit_file,'(3(e14.7,x))') field(1,n), field(2,n), field(3,n)
          END DO
       END IF
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</PointData>'
    !===End of PointData Block

    !===CellData Block
    WRITE(unit_file,'(A)') '<CellData>'
    WRITE(unit_file,'(A)') '</CellData>'
    !===End of CellData Block

    !===Points Block
    WRITE(unit_file,'(A)') '<Points>'
    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
         'NumberOfComponents="3" format="'//TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO n = 1, mesh%np
          r4_threed_xml_field(3*(n-1)+1:3*n) = REAL(mesh%rr(1:3,n),4)
       END DO
       CALL write_compressed_block(unit_file, r4_threed_xml_field(1:3*mesh%np))
    ELSE
       DO n = 1, mesh%np
          WRITE(unit_file,'(3(e14.7,x))') mesh%rr(1,n), mesh%rr(2,n), mesh%rr(3,n)
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Points>'
    !===End of Points Block

    !===Cells Block
    WRITE(unit_file,'(A)') '<Cells>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="connectivity" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_threed_xml_field(mesh%gauss%n_w*(m-1)+1:mesh%gauss%n_w*m) = INT(mesh%jj(1:mesh%gauss%n_w,m)-1,4)
       END DO
       CALL write_compressed_block(unit_file, i4_threed_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(15(I8,1x))') mesh%jj(:,m)-1
       END DO
    END IF

    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="Int32" Name="offsets" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i4_xml_field(m) = INT(stride*m,4)
       END DO
       CALL write_compressed_block(unit_file, i4_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') stride*m
       END DO
    END IF

    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="'&
         //TRIM(ADJUSTL(ascii_or_binary))//'">'
    IF (inputs%if_xml) THEN
       DO m = 1, mesh%me
          i1_xml_field(m) = INT(type_cell,1)
       END DO
       CALL write_compressed_block(unit_file, i1_xml_field)
    ELSE
       DO m = 1, mesh%me
          WRITE(unit_file,'(I8)') type_cell
       END DO
    END IF
    WRITE(unit_file,'(A)') '</DataArray>'
    WRITE(unit_file,'(A)') '</Cells>'
    !===End of Cells Block

    WRITE(unit_file,'(A)') '</Piece>'
    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
    WRITE(unit_file,'(A)') '</VTKFile>'

    CLOSE(unit_file)
  END SUBROUTINE create_xml_vtu_file_3D

!!$  SUBROUTINE make_vtu_file_scalar_2D(comm, mesh, header, field, field_name, what, opt_it)
!!$    USE def_type_mesh
!!$    USE my_util
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),               INTENT(IN) :: mesh
!!$    CHARACTER(*),                  INTENT(IN) :: header
!!$    CHARACTER(*),                  INTENT(IN) :: field_name, what
!!$    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
!!$    REAL(KIND=8), DIMENSION(:),    INTENT(IN) :: field
!!$    INTEGER                                   :: j, it
!!$    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
!!$    CHARACTER(LEN=3)                          :: st_rank, st_it
!!$!#include "petsc/finclude/petsc.h"
!!$    PetscErrorCode                            :: ierr
!!$    PetscMPIInt                               :: rank, nb_procs
!!$    MPI_Comm                                  :: comm
!!$    CALL MPI_Comm_rank(comm, rank, ierr)
!!$    CALL MPI_Comm_size(comm, nb_procs, ierr)
!!$    ALLOCATE(file_list(nb_procs))
!!$    IF (PRESENT(opt_it)) THEN
!!$       it = opt_it
!!$       WRITE(st_it,'(I3)') it
!!$       DO j = 1, nb_procs
!!$          WRITE(st_rank,'(I3)') j
!!$          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))//&
!!$               '_it_'//TRIM(ADJUSTL(st_it))
!!$       END DO
!!$    ELSE
!!$       DO j = 1, nb_procs
!!$          WRITE(st_rank,'(I3)') j
!!$          file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))
!!$       END DO
!!$    END IF
!!$
!!$    CALL check_list(comm, file_list, mesh%np)
!!$    IF (rank==0) THEN
!!$       IF (PRESENT(opt_it)) THEN
!!$          it = opt_it
!!$       ELSE
!!$          it = 1
!!$       END IF
!!$       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
!!$    END IF
!!$    CALL create_xml_vtu_scal_file(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), TRIM(ADJUSTL(field_name)))
!!$  END SUBROUTINE make_vtu_file_scalar_2D

!!$  SUBROUTINE create_vtk_file(comm, field, mesh, file_name, opt_it)
!!$    USE def_type_mesh
!!$    USE chaine_caractere
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type)                        :: mesh
!!$    REAL(KIND=8), DIMENSION(:), INTENT(IN) :: field
!!$    CHARACTER(*),         INTENT(IN)       :: file_name
!!$    INTEGER, OPTIONAL                      :: opt_it
!!$    CHARACTER(len=4)                       :: tit
!!$    CHARACTER(len=3)                       :: st_it
!!$    INTEGER                                :: unit_file, l, lblank, start, i, j, nit
!!$    PetscErrorCode :: ierr
!!$    PetscMPIInt    :: rank, nb_procs
!!$    MPI_Comm       :: comm
!!$
!!$    CALL MPI_Comm_rank(comm,rank,ierr)
!!$    CALL MPI_Comm_size(comm,nb_procs,ierr)
!!$
!!$    IF (PRESENT(opt_it)) THEN
!!$       WRITE(st_it,'(i3)') opt_it
!!$       lblank = eval_blank(3,st_it)
!!$       DO l = 1, lblank - 1
!!$          st_it(l:l) = '0'
!!$       END DO
!!$       nit = opt_it
!!$    ELSE
!!$       st_it='001'
!!$       nit = 1
!!$    END IF
!!$
!!$    unit_file=rank+10
!!$    IF (rank==0) THEN
!!$       OPEN (UNIT=unit_file, FILE=file_name//'.visit', FORM = 'formatted', STATUS = 'unknown')
!!$       WRITE(unit_file,'(A,i7)') '!NBLOCKS ', nb_procs
!!$       DO j = 0, nb_procs-1
!!$          WRITE(tit,'(I4)') j
!!$          lblank = eval_blank(4,tit)
!!$          DO l = 1, lblank - 1
!!$             tit(l:l) = '0'
!!$          END DO
!!$          start = 4 - nb_digit() + 1
!!$          WRITE(unit_file,'(A)') file_name//'_'//tit(start:)//'_'//st_it//'.vtk'
!!$       END DO
!!$       CLOSE(unit_file) !Ecriture pour visit
!!$
!!$       OPEN (UNIT=unit_file, FILE=file_name//'.pvd', FORM = 'formatted', STATUS = 'unknown')
!!$       IF (nit == 1) THEN
!!$          WRITE(unit_file, '(A)') '<?xml version="1.0"?>'
!!$          WRITE(unit_file, '(A)') '<VTKFile type="Collection" version="0.1" '// &
!!$               'byte_order="BigEndian" compressor="vtkZLibDataCompressor">'
!!$          WRITE(unit_file, '(A)') '<Collection>'
!!$       ELSE
!!$          DO j = 1, 3+(nit-1)*nb_procs
!!$             READ(unit_file,*)
!!$          END DO
!!$       END IF
!!$       DO j = 0, nb_procs-1
!!$          WRITE(tit,'(I4)') j
!!$          lblank = eval_blank(4,tit)
!!$          DO l = 1, lblank - 1
!!$             tit(l:l) = '0'
!!$          END DO
!!$          start = 4 - nb_digit() + 1
!!$          WRITE(unit_file,'(A)') '<DataSet timestep="'//st_it//'" group="" '// &
!!$               'part="'//tit(2:)//'" file="./'//file_name//'_'//tit(start:)//'_'//st_it//'.vtu'//'"/>'
!!$       END DO
!!$       WRITE(unit_file, '(A)') '</Collection>'
!!$       WRITE(unit_file, '(A)') '</VTKFile>'
!!$
!!$       CLOSE(unit_file) !Ecriture pour paraview
!!$    END IF
!!$
!!$    IF (SIZE(field,1)==0) RETURN
!!$
!!$    WRITE(tit,'(I4)') rank
!!$    lblank = eval_blank(4,tit)
!!$    DO l = 1, lblank - 1
!!$       tit(l:l) = '0'
!!$    END DO
!!$    start = 4 - nb_digit() + 1
!!$
!!$    unit_file=rank+10
!!$    OPEN (UNIT=unit_file, FILE=file_name//'_'//tit(start:)//'_'//st_it//'.vtk', &
!!$         FORM = 'formatted', STATUS = 'unknown')
!!$    WRITE(unit_file,'(A)') '# vtk DataFile Version 3.0'
!!$    WRITE(unit_file,'(A)') 'vtk '//file_name//''
!!$    WRITE(unit_file,'(A)')'ASCII'
!!$    WRITE(unit_file,'(A)')'DATASET UNSTRUCTURED_GRID'
!!$
!!$    WRITE(unit_file,'(A,I7,A)')'POINTS ', mesh%np, ' float'
!!$    WRITE(*,*) 'points ...'
!!$    DO i=1, mesh%np
!!$       WRITE(unit_file,'(2(e14.7,2x),A)') mesh%rr(1,i), &
!!$            mesh%rr(2,i), ' 0.0 '
!!$    ENDDO
!!$    WRITE(*,*) 'cells ...'
!!$    IF (mesh%gauss%n_w==3) THEN
!!$       WRITE(unit_file,'(A,I7,I8)') 'CELLS ', mesh%me, mesh%me*4
!!$       DO i=1, mesh%me
!!$          WRITE(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(1,i)-1, mesh%jj(2,i)-1, mesh%jj(3,i)-1
!!$       ENDDO
!!$       WRITE(unit_file,'(A,I7)') 'CELL_TYPES ', mesh%me
!!$       DO i=1, mesh%me
!!$          WRITE(unit_file,'(A)') '5'
!!$       ENDDO
!!$    ELSE IF (mesh%gauss%n_w==6) THEN
!!$       WRITE(unit_file,'(A,I7,I8)') 'CELLS ', 4*mesh%me, 4*mesh%me*4
!!$       DO i=1, mesh%me
!!$          WRITE(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(1,i)-1, mesh%jj(6,i)-1, mesh%jj(5,i)-1
!!$          WRITE(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(2,i)-1, mesh%jj(4,i)-1, mesh%jj(6,i)-1
!!$          WRITE(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(3,i)-1, mesh%jj(5,i)-1, mesh%jj(4,i)-1
!!$          WRITE(unit_file,'(A,3(I8,1x))') '3 ',  mesh%jj(6,i)-1, mesh%jj(4,i)-1, mesh%jj(5,i)-1
!!$       ENDDO
!!$       WRITE(unit_file,'(A,I7)') 'CELL_TYPES ', 4*mesh%me
!!$
!!$       DO i=1, 4*mesh%me
!!$          WRITE(unit_file,'(A)') '5'
!!$       ENDDO
!!$    END IF
!!$
!!$    WRITE(*,*) 'data ...'
!!$    WRITE(unit_file,'(A,I7)') 'POINT_DATA ',mesh%np
!!$    WRITE(unit_file,'(A)') 'SCALARS scalars float 1'
!!$    WRITE(unit_file,'(A)') 'LOOKUP_TABLE default'
!!$    DO i=1, mesh%np
!!$       WRITE(unit_file,'(e14.7,2x)') field(i)
!!$    ENDDO
!!$    CLOSE(unit_file)
!!$  END SUBROUTINE create_vtk_file

!!$  SUBROUTINE create_vtu_vect_file(field, mesh, file_name, opt_st)
!!$    USE def_type_mesh
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type)                          :: mesh
!!$    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: field
!!$    CHARACTER(*),                 INTENT(IN) :: file_name
!!$    INTEGER                                  :: unit_file=789, m, i, type_cell
!!$    CHARACTER(*), OPTIONAL,       INTENT(IN) :: opt_st
!!$    CHARACTER(LEN=200)                       :: field_name
!!$
!!$    IF (PRESENT(opt_st)) THEN
!!$       field_name=opt_st
!!$    ELSE
!!$       field_name='field'
!!$    END IF
!!$
!!$    IF (SIZE(field)==0) RETURN
!!$    OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
!!$         FORM = 'formatted', STATUS = 'unknown')
!!$
!!$    WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
!!$         ' byte_order="LittleEndian">'
!!$    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
!!$
!!$    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
!!$         '" NumberOfCells="', mesh%me, '">'
!!$    ! PointData Block ------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<PointData>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
!!$         '_cos" format="ascii" NumberOfComponents="3">'
!!$    DO i = 1, mesh%np
!!$       WRITE(unit_file,'(e14.7)') field(i,1),  field(i,3),  field(i,5)
!!$    ENDDO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
!!$         '_sin" format="ascii" NumberOfComponents="3">'
!!$    DO i = 1, mesh%np
!!$       WRITE(unit_file,'(e14.7)') field(i,2),  field(i,4),  field(i,6)
!!$    ENDDO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</PointData>'
!!$    ! End of PointData Block -----------------------------------------------------
!!$
!!$    ! CellData Block -------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<CellData>'
!!$    WRITE(unit_file,'(A)') '</CellData>'
!!$    ! End of CellData Block ------------------------------------------------------
!!$
!!$    ! Points Block ---------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<Points>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
!!$         'NumberOfComponents="3" format="ascii">'
!!$    DO i=1, mesh%np
!!$       WRITE(unit_file,'(e14.7,A,e14.7)') mesh%rr(1,i), ' 0.0 ' , &
!!$            mesh%rr(2,i)
!!$    ENDDO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</Points>'
!!$    ! End of Points Block --------------------------------------------------------
!!$
!!$    ! Cells Block ----------------------------------------------------------------
!!$    IF (mesh%gauss%n_w==3) THEN
!!$       type_cell = 5
!!$    ELSE IF (mesh%gauss%n_w==6) THEN
!!$       type_cell = 22
!!$    END IF
!!$    WRITE(unit_file,'(A)') '<Cells>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="connectivity" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(3(I8,1x))') mesh%jj(1:3,m)-1
!!$       IF (type_cell==22) THEN
!!$          WRITE(unit_file,'(3(I8,1x))') mesh%jj(6,m)-1 , mesh%jj(4,m)-1 , mesh%jj(5,m)-1
!!$       END IF
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="offsets" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(I8)') m*mesh%gauss%n_w
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(I8)') type_cell
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</Cells>'
!!$    ! End of Cells Block ---------------------------------------------------------
!!$
!!$    WRITE(unit_file,'(A)') '</Piece>'
!!$    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
!!$    WRITE(unit_file,'(A)') '</VTKFile>'
!!$
!!$    CLOSE(unit_file)
!!$  END SUBROUTINE create_vtu_vect_file

!!$  SUBROUTINE make_vtu_file_axi3D(communicator, mesh, header, &
!!$       field, field_name, what, opt_it)
!!$    USE def_type_mesh
!!$    USE my_util
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type),               INTENT(IN) :: mesh
!!$    CHARACTER(*),                  INTENT(IN) :: header
!!$    CHARACTER(*),                  INTENT(IN) :: field_name, what
!!$    INTEGER, OPTIONAL,             INTENT(IN) :: opt_it
!!$    REAL(KIND=8), DIMENSION(:,:,:),INTENT(IN) :: field
!!$    INTEGER                                   :: j, it
!!$    CHARACTER(LEN=200), DIMENSION(:), POINTER :: file_list
!!$    CHARACTER(LEN=3)                          :: st_rank
!!$!#include "petsc/finclude/petsc.h"
!!$    PetscErrorCode                            :: ierr
!!$    PetscMPIInt                               :: rank, nb_procs
!!$    MPI_Comm                                  :: communicator
!!$    CALL MPI_Comm_rank(communicator, rank, ierr)
!!$    CALL MPI_Comm_size(communicator, nb_procs, ierr)
!!$    ALLOCATE(file_list(nb_procs))
!!$    DO j = 1, nb_procs
!!$       WRITE(st_rank,'(I3)') j
!!$       file_list(j) = TRIM(header)//'_proc_'//TRIM(ADJUSTL(st_rank))
!!$    END DO
!!$    CALL check_list(communicator, file_list, mesh%np)
!!$    IF (rank==0) THEN
!!$       IF (PRESENT(opt_it)) THEN
!!$          it = opt_it
!!$       ELSE
!!$          it = 1
!!$
!!$       END IF
!!$       CALL create_pvd_file(file_list, TRIM(header), it, TRIM(what))
!!$    END IF
!!$
!!$    CALL create_vtu_file_axi3D(field, mesh, TRIM(ADJUSTL(file_list(rank+1))), &
!!$         opt_st=field_name)
!!$  END SUBROUTINE make_vtu_file_axi3D

!!$  SUBROUTINE create_vtu_file_3D(field, mesh, file_name, opt_st)
!!$    USE def_type_mesh
!!$    USE my_util
!!$    IMPLICIT NONE
!!$    TYPE(mesh_type)                            :: mesh
!!$    REAL(KIND=8), DIMENSION(:,:),   INTENT(IN) :: field
!!$    CHARACTER(*),                   INTENT(IN) :: file_name
!!$    CHARACTER(*), OPTIONAL,         INTENT(IN) :: opt_st
!!$    CHARACTER(LEN=200)                         :: field_name
!!$    INTEGER                                    :: unit_file=789, m, n, type_cell, stride
!!$
!!$    IF (SIZE(field,2)==0) RETURN
!!$
!!$    IF (SIZE(mesh%jj,1)==6) THEN
!!$       type_cell = 13
!!$       stride = 6
!!$    ELSE IF (SIZE(mesh%jj,1)==15) THEN
!!$       type_cell = 26
!!$       stride = 15
!!$    ELSE
!!$       type_cell = 0
!!$       stride = 0
!!$       CALL error_petsc('Bug in create_vtu_file_3D: SIZE(mesh%jj,1) is wrong')
!!$    END IF
!!$
!!$    IF (PRESENT(opt_st)) THEN
!!$       field_name=opt_st
!!$    ELSE
!!$       field_name='field'
!!$    END IF
!!$
!!$    OPEN (UNIT=unit_file, FILE=file_name//'.vtu',&
!!$         FORM = 'formatted', STATUS = 'unknown')
!!$
!!$    WRITE(unit_file,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1"'// &
!!$         ' byte_order="LittleEndian">'
!!$    WRITE(unit_file,'(A)') '<UnstructuredGrid>'
!!$
!!$    WRITE(unit_file,'(A,I9,A,I9,A)') '<Piece NumberOfPoints="', mesh%np, &
!!$         '" NumberOfCells="', mesh%me, '">'
!!$    ! PointData Block ------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<PointData Scalars="truc">'
!!$
!!$    IF (SIZE(field,1)==1) THEN
!!$       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
!!$            '" format="ascii">'
!!$       DO n = 1, mesh%np
!!$          WRITE(unit_file,'(e14.7)') field(1,n)
!!$       END DO
!!$    ELSE
!!$       WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="'//TRIM(ADJUSTL(field_name))//&
!!$            '" format="ascii" NumberOfComponents="3">'
!!$       DO n = 1, mesh%np
!!$          WRITE(unit_file,'(3(e14.7,x))') field(1,n), field(2,n), field(3,n)
!!$       END DO
!!$    END IF
!!$
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</PointData>'
!!$    ! End of PointData Block -----------------------------------------------------
!!$
!!$    ! CellData Block -------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<CellData>'
!!$    WRITE(unit_file,'(A)') '</CellData>'
!!$    ! End of CellData Block ------------------------------------------------------
!!$
!!$    ! Points Block ---------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<Points>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Float32" Name="Points" '//&
!!$         'NumberOfComponents="3" format="ascii">'
!!$    DO n = 1, mesh%np
!!$       WRITE(unit_file,'(3(e14.7,x))') mesh%rr(1,n), mesh%rr(2,n), mesh%rr(3,n)
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</Points>'
!!$    ! End of Points Block --------------------------------------------------------
!!$
!!$    ! Cells Block ----------------------------------------------------------------
!!$    WRITE(unit_file,'(A)') '<Cells>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="connectivity" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(15(I8,1x))') mesh%jj(:,m)-1
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="Int64" Name="offsets" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(I8)') stride*m
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
!!$    DO m = 1, mesh%me
!!$       WRITE(unit_file,'(I8)') type_cell
!!$    END DO
!!$    WRITE(unit_file,'(A)') '</DataArray>'
!!$    WRITE(unit_file,'(A)') '</Cells>'
!!$    ! End of Cells Block ---------------------------------------------------------
!!$
!!$    WRITE(unit_file,'(A)') '</Piece>'
!!$    WRITE(unit_file,'(A)') '</UnstructuredGrid>'
!!$    WRITE(unit_file,'(A)') '</VTKFile>'
!!$
!!$    CLOSE(unit_file)
!!$  END SUBROUTINE create_vtu_file_3D

!!$  FUNCTION nb_digit() RESULT(dg)
!!$    USE my_util
!!$    IMPLICIT NONE
!!$    INTEGER :: code, nb_procs, dg
!!$!#include "petsc/finclude/petsc.h"
!!$    CALL MPI_Comm_size(PETSC_COMM_WORLD,nb_procs,code)
!!$
!!$    IF (nb_procs>9999) THEN
!!$       CALL error_Petsc('nb_procs>9999')
!!$    END IF
!!$
!!$    IF (nb_procs < 10) THEN
!!$       dg=1
!!$    ELSE IF (nb_procs < 99) THEN
!!$       dg=2
!!$    ELSE IF (nb_procs < 999) THEN
!!$       dg=3
!!$    ELSE IF (nb_procs < 9999) THEN
!!$       dg=4
!!$    ELSE
!!$       dg=0
!!$    END IF
!!$  END FUNCTION nb_digit
!!$
END MODULE vtk_viz
