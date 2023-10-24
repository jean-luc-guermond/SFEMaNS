MODULE user_data_module
  TYPE personalized_data
     !===I declare my own data here==================================================
     LOGICAL                                 :: if_my_stuff
     !.......Continue here ................................
     !===End I declare my own data here==============================================
  END TYPE personalized_data
END MODULE user_data_module

MODULE user_data
  USE user_data_module
  IMPLICIT NONE
  PUBLIC :: read_user_data
  TYPE(personalized_data), PUBLIC  :: user
  PRIVATE

CONTAINS

  SUBROUTINE read_user_data(data_file)
    USE my_util
    USE chaine_caractere
    IMPLICIT NONE
    CHARACTER(*),       INTENT(IN) :: data_file
    INTEGER                        :: unit_file=22
    LOGICAL                        :: test

    OPEN(UNIT=unit_file, FILE = data_file, FORM = 'formatted', STATUS = 'unknown')

    !===I add lines that the code SFEMaNS reads in the data file=========================
    CALL find_string(unit_file, '===Should I read my stuff? (true/false)', test)
    IF (test) THEN
       READ (unit_file, *) user%if_my_stuff
    ELSE
       user%if_my_stuff = .FALSE.
    END IF
    !.......Continue here ................................
    !===End I add lines that the code SFEMaNS reads in the data file=====================

    CLOSE(unit_file)
  END SUBROUTINE read_user_data

END MODULE user_data
