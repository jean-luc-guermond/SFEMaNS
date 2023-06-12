!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE chaine_caractere

  PUBLIC :: last_c_leng, last_of_string, start_of_string

CONTAINS

  FUNCTION last_c_leng (len_str, string) RESULT (leng)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .EQ. ' ' ) THEN
          leng = i-1; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_c_leng

  !========================================================================

  FUNCTION  eval_blank(len_str, string) RESULT (leng)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: len_str
    CHARACTER (LEN=len_str), INTENT(IN) :: string
    INTEGER :: leng

    INTEGER :: i

    leng = len_str

    DO i=1,len_str
       IF ( string(i:i) .NE. ' ' ) THEN
          leng = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION eval_blank

  !========================================================================

  FUNCTION start_of_string (string) RESULT (start)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: string
    INTEGER :: start

    INTEGER :: i

    start = 1

    DO i = 1, LEN(string)
       IF ( string(i:i) .NE. ' ' ) THEN
          start = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION start_of_string

  !========================================================================

  FUNCTION last_of_string (string) RESULT (last)

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT(IN) :: string
    INTEGER :: last

    INTEGER :: i

    last = 1

    DO i = LEN(string), 1, -1
       IF ( string(i:i) .NE. ' ' ) THEN
          last = i; EXIT
       ENDIF
    ENDDO

  END FUNCTION last_of_string
  !========================================================================

  SUBROUTINE read_until(unit, string, error)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    INTEGER                            :: d_end, d_start
    LOGICAL, OPTIONAL                  :: error
    IF (PRESENT(error)) error =.FALSE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       d_start = start_of_string(control)
       d_end =   last_of_string(control)
       IF (control(d_start:d_end)==string) RETURN
    END DO

    RETURN
11  WRITE(*,*) ' Error reading data file '; IF (PRESENT(error)) error=.TRUE.; RETURN
22  WRITE(*,*) ' Data string ',string,' not found '; IF (PRESENT(error)) error=.TRUE.; RETURN

  END SUBROUTINE read_until

  SUBROUTINE find_string(unit, string, okay)
    IMPLICIT NONE
    INTEGER, PARAMETER                 :: long_max=128
    INTEGER,                INTENT(IN) :: unit
    CHARACTER(LEN=*),       INTENT(IN) :: string
    CHARACTER(len=long_max)            :: control
    INTEGER                            :: d_end, d_start
    LOGICAL                            :: okay

    okay = .TRUE.
    REWIND(unit)
    DO WHILE (.TRUE.)
       READ(unit,'(64A)',ERR=11,END=22) control
       d_start = start_of_string(control)
       d_end =   last_of_string(control)
       IF (control(d_start:d_end)==string) RETURN
    END DO

11  WRITE(*,*) ' Erreur de lecture '; STOP
22  okay = .FALSE.; RETURN

  END SUBROUTINE find_string
  !========================================================================

END MODULE chaine_caractere
