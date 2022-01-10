!
!Authors: Jean-Luc Guermond, Lugi Quartapelle, Copyright 1994
!
MODULE dyn_line
  TYPE dyn_int_line
     INTEGER, DIMENSION(:), POINTER :: DIL  ! POINTER attribute
  END TYPE dyn_int_line                     ! instead of ALLOCATABLE

  TYPE dyn_real_line
     REAL(KIND=8), DIMENSION(:), POINTER :: DRL  ! POINTER attribute
  END TYPE dyn_real_line                        ! instead of ALLOCATABLE

  TYPE dyn_real_array
     REAL(KIND=8), DIMENSION(:,:), POINTER :: DRA  ! POINTER attribute
  END TYPE dyn_real_array                       ! instead of ALLOCATABLE

  TYPE dyn_real_array_three
     REAL(KIND=8), DIMENSION(:,:,:), POINTER :: DRT  ! POINTER attribute
  END TYPE dyn_real_array_three                   ! instead of ALLOCATABLE
END MODULE dyn_line
