MODULE zlib_base64
  USE iso_C_binding, ONLY : C_INT, C_PTR, c_loc
  IMPLICIT NONE

  INTERFACE
     SUBROUTINE write_compressed_block_(a, b, c) &
          BIND(C, name="write_compressed_block")
       IMPORT C_INT, C_PTR
       IMPLICIT NONE
       INTEGER(KIND=C_INT) :: a, c
       TYPE(C_PTR) :: b
     END SUBROUTINE write_compressed_block_
  END INTERFACE

  INTERFACE write_compressed_block
     MODULE PROCEDURE &
          write_compressed_block_I1, &
          write_compressed_block_I4, &
          write_compressed_block_R4, &
          write_compressed_block_R8
  END INTERFACE write_compressed_block

CONTAINS

  SUBROUTINE write_compressed_block_I1(unit, ptr)
#ifdef __INTEL_COMPILER
    USE ifposix
#endif
    IMPLICIT NONE
    INTEGER(KIND=4) :: unit, length, fd, ierror
    INTEGER(KIND=1), TARGET :: ptr(:)
    length = 1 * SIZE(ptr)
    FLUSH(unit)
#ifdef __INTEL_COMPILER
    CALL pxffileno(unit, fd, ierror)
#else
    fd = fnum(unit)
    ierror = 0 !===To avoid warning
#endif
    CALL write_compressed_block_(fd, c_LOC(ptr), length)
  END SUBROUTINE write_compressed_block_I1

  SUBROUTINE write_compressed_block_I4(unit, ptr)
#ifdef __INTEL_COMPILER
    USE ifposix
#endif
    IMPLICIT NONE
    INTEGER(KIND=4) :: unit, length, fd, ierror
    INTEGER(KIND=4), TARGET :: ptr(:)
    length = 4 * SIZE(ptr)
    FLUSH(unit)
#ifdef __INTEL_COMPILER
    CALL pxffileno(unit, fd, ierror)
#else
    fd = fnum(unit)
    ierror = 0 !===To avoid warning
#endif
    CALL write_compressed_block_(fd, c_LOC(ptr), length)
  END SUBROUTINE write_compressed_block_I4

  SUBROUTINE write_compressed_block_R4(unit, ptr)
#ifdef __INTEL_COMPILER
    USE ifposix
#endif
    IMPLICIT NONE
    INTEGER(KIND=4) :: unit, length, fd, ierror
    REAL(KIND=4), TARGET :: ptr(:)
    length = 4 * SIZE(ptr)
    FLUSH(unit)
#ifdef __INTEL_COMPILER
    CALL pxffileno(unit, fd, ierror)
#else
    fd = fnum(unit)
    ierror = 0 !===To avoid warning
#endif
    CALL write_compressed_block_(fd, c_LOC(ptr), length)
  END SUBROUTINE write_compressed_block_R4

  SUBROUTINE write_compressed_block_R8(unit, ptr)
#ifdef __INTEL_COMPILER
    USE ifposix
#endif
    IMPLICIT NONE
    INTEGER(KIND=4) :: unit, length, fd, ierror
    REAL(KIND=8), TARGET :: ptr(:)
    length = 8 * SIZE(ptr)
    FLUSH(unit)
#ifdef __INTEL_COMPILER
    CALL pxffileno(unit, fd, ierror)
#else
    fd = fnum(unit)
    ierror = 0 !===To avoid warning
#endif
    CALL write_compressed_block_(fd, c_LOC(ptr), length)
  END SUBROUTINE write_compressed_block_R8

END MODULE zlib_base64
