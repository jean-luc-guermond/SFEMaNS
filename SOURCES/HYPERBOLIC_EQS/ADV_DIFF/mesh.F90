MODULE mesh_handling
  USE def_type_mesh
  USE dir_nodes
  TYPE(mesh_type), PUBLIC                :: mesh
  INTEGER, POINTER, DIMENSION(:), PUBLIC :: js_D_loc, js_D_glob
  type(petsc_csr_LA), PUBLIC             :: LA
  INTEGER                                :: max_nnz_per_row
CONTAINS
  SUBROUTINE construct_mesh
#include "petsc/finclude/petsc.h"
    USE petsc
    USE input_data
    USE prep_maill
    USE mod_gauss_points_2d
    USE metis_reorder_elements
    USE st_matrix
    USE dir_nodes_petsc
    IMPLICIT NONE
    TYPE(mesh_type)               :: mesh_glob
    INTEGER                       :: dom_np
    PetscErrorCode :: ierr
    PetscMPIInt    :: nb_proc

    !===Number of procs
    CALL MPI_Comm_size(PETSC_COMM_WORLD,nb_proc,ierr)

    !===Prepare the grid
    CALL load_mesh_free_format_ordered(inputs%directory,inputs%file_name,inputs%list_dom,&
         inputs%type_fe,mesh_glob,inputs%if_mesh_formatted,edge_stab=.FALSE.)
    CALL incr_vrtx_indx_enumeration(mesh_glob,inputs%type_fe)

    !===Metis reorganizes the mesh
    CALL reorder_mesh(PETSC_COMM_WORLD,nb_proc,mesh_glob,mesh)

    !===Create Sparsity pattern for the matrix (structure)
    CALL st_aij_csr_glob_block(PETSC_COMM_WORLD,1,mesh_glob,mesh,LA)
    dom_np = SIZE(LA%ia) - 1
    max_nnz_per_row = MAXVAL(LA%ia(1:dom_np)-LA%ia(0:dom_np-1))

    !===Deallocate global mesh
    CALL free_mesh_after(mesh_glob)

    !===Start Gauss points generation
    mesh%edge_stab=.FALSE.
    CALL gauss_points_2d(mesh,inputs%type_fe)

    !===Boundary conditions on local domain
    CALL dirichlet_nodes_parallel(mesh, inputs%Dir_list, js_d_loc)
    ALLOCATE(js_D_glob(SIZE(js_D_loc)))
    js_D_glob = LA%loc_to_glob(1,js_D_loc)-1

  END SUBROUTINE construct_mesh

END MODULE MESH_HANDLING
