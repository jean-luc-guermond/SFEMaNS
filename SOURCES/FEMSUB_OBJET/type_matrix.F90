MODULE matrix_type
  TYPE matrice_bloc
     REAL(KIND=8), POINTER, DIMENSION(:) :: aa
     INTEGER,      POINTER, DIMENSION(:) :: ia, ja
  END TYPE matrice_bloc
END MODULE matrix_type
