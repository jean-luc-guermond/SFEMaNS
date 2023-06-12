MODULE sfemans_tools

CONTAINS

  FUNCTION find_point(mesh,r,z)  RESULT(n)

    USE def_type_mesh

    IMPLICIT NONE
    TYPE(mesh_type),                 INTENT(IN) :: mesh !type de maillage
    REAL(KIND=8),                    INTENT(IN) :: r, z
    INTEGER                                     :: n
    REAL(KIND=8)                                :: hmax
    INTEGER, DIMENSION(1) :: jlg

    hmax = MAXVAL(mesh%hloc)

    jlg = MINLOC((mesh%rr(1,:)-r)**2 + (mesh%rr(2,:)-z)**2)
    n = jlg(1)

    IF (((mesh%rr(1,n)-r)**2 + (mesh%rr(2,n)-z)**2) .GT. hmax**2) THEN
       n = 0
    END IF

    IF (n .GT. mesh%dom_np) THEN
       n = 0
    END IF

    RETURN

  END FUNCTION find_point

END MODULE sfemans_tools
