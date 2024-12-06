MODULE create_comm
  PUBLIC :: create_cart_comm
  PRIVATE
CONTAINS
  SUBROUTINE create_cart_comm(ndim,comm_cart,comm_one_d,coord_cart)
#include "petsc/finclude/petsc.h"
    USE petsc
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN)  :: ndim
    INTEGER,               INTENT(OUT) :: comm_cart
    INTEGER, DIMENSION(:), POINTER     :: comm_one_d, coord_cart
    LOGICAL, DIMENSION(SIZE(ndim))     :: period, remain
    LOGICAL :: reorder
    INTEGER :: dim, nb_procs, code, n, i, rank

    !==Verification of compatibility==!
    CALL MPI_COMM_SIZE(PETSC_COMM_WORLD, nb_procs, code)
    dim = SIZE(ndim)
    n= 1
    DO i = 1, dim
       n = n*ndim(i)
    END DO
    IF (nb_procs/=n) THEN
       WRITE(*,*) ' CREATE_CART_COMM: Nb of procs not compatible with Cartesian decomposition'
       STOP
    END IF

    !==Create Cartesian communication==!
    ALLOCATE(comm_one_d(dim),coord_cart(dim))
    period   = .FALSE.
    reorder  = .FALSE.
    CALL MPI_CART_CREATE(PETSC_COMM_WORLD, dim, ndim, period, reorder, comm_cart, code)
    CALL MPI_COMM_RANK(comm_cart, rank, code)
    CALL MPI_CART_COORDS(comm_cart, rank, dim, coord_cart, code)

    !==Create one D communication for each dimension==!
    DO i = 1, dim
       remain = .FALSE.
       remain(i) = .TRUE.
       CALL MPI_CART_SUB(comm_cart, remain, comm_one_d(i), code)
    END DO
  END SUBROUTINE create_cart_comm

END MODULE create_comm
