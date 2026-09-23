MODULE user_data_module
  TYPE personalized_data
     !===I declare my own data here==================================================
     LOGICAL                                 :: if_my_stuff
     LOGICAL :: if_nonlin_m, if_nonlin
     REAL(KIND=8) :: amp_nl, amp_nl_m
     INTEGER :: mode_ampli
     !.......Continue here ................................

     LOGICAL                                 :: if_solid_walls
     REAL(KIND=8)                            :: Omega_Int
     REAL(KIND=8)                            :: Omega_Out
     REAL(KIND=8)                            :: Bz_imposed  

     ! my information
     !Time avg info and restart
     LOGICAL                                 :: compute_time_average_info, read_time_average_info
     INTEGER                                 :: time_avg_freq_restart

     !Forts info
     INTEGER      :: nb_fort_mom, nb_fort_mom_avg
     INTEGER      :: nb_fort_mag, nb_fort_mag_avg
     INTEGER      :: nb_fort_conc, nb_fort_temp
     INTEGER, DIMENSION(:), ALLOCATABLE   :: list_fort_mom, list_fort_mom_avg
     INTEGER, DIMENSION(:), ALLOCATABLE   :: list_fort_mag, list_fort_mag_avg
     INTEGER, DIMENSION(:), ALLOCATABLE   :: list_fort_conc, list_fort_temp
     INTEGER :: pod_tub_slices
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
        !===Nonlinear restart================================================
        CALL find_string(unit_file, '===Is it a nonlinear restart on modes/=0 ?(true/false)', test)
        IF (test) THEN
            READ(unit_file, *) user%if_nonlin
            IF (user%if_nonlin) THEN
                CALL find_string(unit_file, '===Amplitude of amp_nl', test)
                READ(unit_file, *) user%amp_nl
            ELSE
                user%if_nonlin = .FALSE.
                user%amp_nl = 1.d0
            END IF
        ELSE
            user%if_nonlin = .FALSE.
            user%amp_nl = 1.d0
        END IF

        CALL find_string(unit_file, '===Is it a nonlinear restart on one mode?(true/false)', test)
        IF (test) THEN
            READ(unit_file, *) user%if_nonlin_m
            IF (user%if_nonlin_m) THEN
                CALL find_string(unit_file, '===Mode amplified', test)
                READ(unit_file, *) user%mode_ampli
                CALL find_string(unit_file, '===Amplitude of amp_nl_m', test)
                READ(unit_file, *) user%amp_nl_m
            ELSE
                user%if_nonlin_m = .FALSE.
                user%mode_ampli = 0
                user%amp_nl_m = 1.d0
            END IF
        ELSE
            user%if_nonlin_m = .FALSE.
            user%mode_ampli = 0
            user%amp_nl_m = 1.d0
        END IF
        !   WRITE(*,*) 'user%if_nonlin_m = ', user%if_nonlin_m
        !   WRITE(*,*) 'user%mode_ampli = ', user%mode_ampli
        !   WRITE(*,*) 'user%amp_nl_m = ', user%amp_nl_m

        !=== For MRI
        CALL find_string(unit_file, '===Inner angular velocity?', test)
        READ(unit_file,*) user%Omega_Int
        CALL find_string(unit_file, '===Outer angular velocity?', test)
        READ(unit_file,*) user%Omega_Out
        CALL find_string(unit_file, '===Imposed axial magnetic field?', test)
        READ (unit_file, *) user%Bz_imposed

        !===Read time average info ===============================================
        CALL find_string(unit_file, '===Compute time average info? (true/false)', test)
        IF (test) THEN
           READ (unit_file,*) user%compute_time_average_info
        ELSE
           user%compute_time_average_info=.FALSE.
        END IF

        CALL find_string(unit_file, '===Read time average info? (true/false)', test)
        IF (test) THEN
            READ (unit_file, *) user%read_time_average_info
        ELSE
            user%read_time_average_info = .FALSE.
        END IF

        CALL find_string(unit_file, '===Frequency to write time avg restart file', test)
        IF (test) THEN
            READ (unit_file, *) user%time_avg_freq_restart
        ELSE
            user%time_avg_freq_restart = 1d5
        END IF
!===Read forts info ===============================================
        ! Momentum
        CALL find_string(unit_file, '===How many forts for momentum? (in 11,88,98,888,889,998,999)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_mom
        ELSE
           user%nb_fort_mom=0
        END IF
        ALLOCATE(user%list_fort_mom(user%nb_fort_mom))
        IF (user%nb_fort_mom > 0) THEN
          CALL read_until(unit_file, '===List of forts for momentum', test)
          READ (unit_file,*) user%list_fort_mom
        END IF
        ! Avg Momentum
        CALL find_string(unit_file, '===How many forts for avg momentum? (in 86,97,99)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_mom_avg
        ELSE
           user%nb_fort_mom_avg=0
        END IF
        ALLOCATE(user%list_fort_mom_avg(user%nb_fort_mom_avg))
        IF (user%nb_fort_mom_avg > 0) THEN
          CALL read_until(unit_file, '===List of forts for momentum avg', test)
          READ (unit_file,*) user%list_fort_mom_avg
        END IF
    
        ! Induction
        CALL find_string(unit_file, '===How many forts for induction? (in 41,51,750,751,752,78)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_mag
        ELSE
           user%nb_fort_mag=0
        END IF
        ALLOCATE(user%list_fort_mag(user%nb_fort_mag))
        IF (user%nb_fort_mag > 0) THEN
          CALL read_until(unit_file, '===List of forts for induction', test)
          READ (unit_file,*) user%list_fort_mag
        END IF
        ! Avg Induction
        CALL find_string(unit_file, '===How many forts for avg induction? (in 52)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_mag_avg
        ELSE
           user%nb_fort_mag_avg=0
        END IF
        ALLOCATE(user%list_fort_mag_avg(user%nb_fort_mag_avg))
        IF (user%nb_fort_mag_avg > 0) THEN
          CALL read_until(unit_file, '===List of forts for induction avg', test)
          READ (unit_file,*) user%list_fort_mag_avg
        END IF
    
        ! Concentration
        CALL find_string(unit_file, '===How many forts for concentration? (in 199)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_conc
        ELSE
           user%nb_fort_conc=0
        END IF
        ALLOCATE(user%list_fort_conc(user%nb_fort_conc))
        IF (user%nb_fort_conc > 0) THEN
          CALL read_until(unit_file, '===List of forts for concentration', test)
          READ (unit_file,*) user%list_fort_conc
        END IF
    
        ! Temperature
        CALL find_string(unit_file, '===How many forts for temperature? (in 299)', test)
        IF (test) THEN
           READ (unit_file,*) user%nb_fort_temp
        ELSE
           user%nb_fort_temp=0
        END IF
        ALLOCATE(user%list_fort_temp(user%nb_fort_temp))
        IF (user%nb_fort_temp > 0) THEN
          CALL read_until(unit_file, '===List of forts for temperature', test)
          READ (unit_file,*) user%list_fort_temp
        END IF
    
        ! POD
        CALL find_string(unit_file, '===How many slices for Tub from POD? (fort.7XX)', test)
        IF (test) THEN
           READ (unit_file,*) user%pod_tub_slices
        ELSE
           user%pod_tub_slices=0
        END IF

    !===End I add lines that the code SFEMaNS reads in the data file=====================

    CLOSE(unit_file)
  END SUBROUTINE read_user_data

END MODULE user_data
