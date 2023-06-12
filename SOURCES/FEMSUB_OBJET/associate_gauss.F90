!
!Authors: Jean-Luc Guermond, Copyright 2000
!
MODULE Gauss_points

  USE def_type_mesh

  PUBLIC
  INTEGER, PUBLIC :: k_d, n_w, l_G , n_ws, l_Gs
  REAL(KIND=8), DIMENSION(:, :),       POINTER :: ww
  REAL(KIND=8), DIMENSION(:, :),       POINTER :: wws
  REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw
  REAL(KIND=8), DIMENSION(:,    :, :), POINTER :: rnorms
  REAL(KIND=8), DIMENSION(:, :),       POINTER :: rj
  REAL(KIND=8), DIMENSION(:, :),       POINTER :: rjs
  REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dw_s
  REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dwps !special!
  REAL(KIND=8), DIMENSION(:, :, :, :), POINTER :: dws  !SPECIAL!


CONTAINS

  SUBROUTINE gauss(mesh)

    IMPLICIT NONE

    TYPE(mesh_type) :: mesh

    ww => mesh%gauss%ww
    wws => mesh%gauss%wws
    dw => mesh%gauss%dw
    rnorms => mesh%gauss%rnorms
    rj => mesh%gauss%rj
    rjs => mesh%gauss%rjs
    dw_s => mesh%gauss%dw_s
    dwps => mesh%gauss%dwps
    dws => mesh%gauss%dws

    k_d = SIZE(dw,1)
    n_w = SIZE(ww,1)
    l_G = SIZE(ww,2)
    n_ws = SIZE(wws,1)
    l_Gs = SIZE(wws,2)

  END SUBROUTINE gauss

END MODULE Gauss_points
