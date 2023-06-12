MODULE input_data
  IMPLICIT NONE
  PUBLIC :: read_my_data
  TYPE my_data
     CHARACTER(len=200)             :: directory
     CHARACTER(len=200)             :: file_name
     LOGICAL                        :: if_mesh_formatted
     LOGICAL                        :: if_plot_vtu
     INTEGER                        :: nb_dom
     INTEGER, DIMENSION(:), POINTER :: list_dom
     INTEGER                        :: type_fe
     REAL(KIND=8)                   :: Tfinal
     REAL(KIND=8)                   :: CFL
     REAL(KIND=8)                   :: ce
     LOGICAL                        :: if_viscous
     INTEGER                        :: type_test
     REAL(KIND=8)                   :: dt
     INTEGER                        :: nb_Dir_bdy
     INTEGER, DIMENSION(:), POINTER :: Dir_list
     LOGICAL                        :: if_xml
  END type my_data
  TYPE(my_data), PUBLIC  :: inputs
  PRIVATE
CONTAINS
  SUBROUTINE read_my_data(data_fichier)
    USE chaine_caractere
    IMPLICIT NONE
    INTEGER, PARAMETER           :: in_unit=21
    CHARACTER(len=*), INTENT(IN) :: data_fichier
    LOGICAL                      :: test
    OPEN(UNIT = in_unit, FILE = data_fichier, FORM = 'formatted', STATUS = 'unknown')
    CALL read_until(in_unit, "===Name of directory for mesh file===")
    READ (in_unit,*) inputs%directory
    CALL read_until(in_unit, "===Name of mesh file===")
    READ (in_unit,*) inputs%file_name
    CALL read_until(in_unit, "===Is the mesh formatted? (True/False)===")
    READ (in_unit,*) inputs%if_mesh_formatted
    CALL read_until(in_unit, '===Number of subdomains in the mesh===')
    READ(21,*) inputs%nb_dom
    ALLOCATE(inputs%list_dom(inputs%nb_dom))
    CALL read_until(21, '===List of subdomain in the mesh===')
    READ(21,*) inputs%list_dom
    CALL read_until(21, '===Type of finite element===')
    READ(21,*) inputs%type_fe

    CALL read_until(in_unit, "===How many Dirichlet boundaries?===")
    READ (in_unit,*)  inputs%nb_Dir_bdy
    CALL read_until(in_unit, "===List of Dirichlet boundaries?===")
    ALLOCATE(inputs%Dir_list(inputs%nb_Dir_bdy))
    READ (in_unit,*) inputs%Dir_list


    !===Format for paraview============================!
    CALL read_until(21, '===Create vtu files? (true/false)===')
    READ(21,*) inputs%if_plot_vtu
    CALL find_string(21, '===Vtu files in xml format? (true=xml/false=ascii)===', test)
    IF (test) THEN
       READ (21,*) inputs%if_xml
    ELSE
       inputs%if_xml=.TRUE.
    END IF


    CLOSE(in_unit)
  END SUBROUTINE read_my_data
END MODULE input_data
