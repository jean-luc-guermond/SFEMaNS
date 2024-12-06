MODULE mesh_tools

   PUBLIC :: copy_mesh, free_mesh, clean_mesh

CONTAINS

   SUBROUTINE copy_mesh(mesh1, mesh2)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh1, mesh2

      mesh2%me = mesh1%me
      mesh2%mes = mesh1%mes
      mesh2%mes_int = mesh1%mes
      mesh2%np = mesh1%np
      mesh2%nps = mesh1%nps
      mesh2%mi = mesh1%mi
      mesh2%medge = mesh1%medge
      mesh2%medges = mesh1%medges
      mesh2%mextra = mesh1%mextra
      mesh2%mes_extra = mesh1%mes_extra

      mesh2%dom_me = mesh1%dom_me
      mesh2%dom_np = mesh1%dom_np
      mesh2%dom_mes = mesh1%dom_mes
      mesh2%nis = mesh1%nis

      ALLOCATE(mesh2%jj(SIZE(mesh1%jj, 1), SIZE(mesh1%jj, 2)))
      mesh2%jj = mesh1%jj
      ALLOCATE(mesh2%jjs(SIZE(mesh1%jjs, 1), SIZE(mesh1%jjs, 2)))
      mesh2%jjs = mesh1%jjs
      ALLOCATE(mesh2%rr(SIZE(mesh1%rr, 1), SIZE(mesh1%rr, 2)))
      mesh2%rr = mesh1%rr
      ALLOCATE(mesh2%loc_to_glob(SIZE(mesh1%loc_to_glob)))
      mesh2%loc_to_glob = mesh1%loc_to_glob
      ALLOCATE(mesh2%neigh(SIZE(mesh1%neigh, 1), SIZE(mesh1%neigh, 2)))
      mesh2%neigh = mesh1%neigh

      ALLOCATE(mesh2%sides(SIZE(mesh1%sides)))
      mesh2%sides = mesh1%sides
      ALLOCATE(mesh2%neighs(SIZE(mesh1%neighs)))
      mesh2%neighs = mesh1%neighs
      ALLOCATE(mesh2%i_d(SIZE(mesh1%i_d)))
      mesh2%i_d = mesh1%i_d

      ALLOCATE(mesh2%jjs_int(SIZE(mesh1%jjs_int, 1), SIZE(mesh1%jjs_int, 2)))
      mesh2%jjs_int = mesh1%jjs_int
      ALLOCATE(mesh2%sides_int(SIZE(mesh1%sides_int)))
      mesh2%sides_int = mesh1%sides_int
      ALLOCATE(mesh2%neighs_int(2, SIZE(mesh1%neighs_int, 2)))
      mesh2%neighs_int = mesh1%neighs_int

      ALLOCATE(mesh2%sides_extra(SIZE(mesh1%sides_extra)))
      mesh2%sides_extra = mesh1%sides_extra
      ALLOCATE(mesh2%neighs_extra(SIZE(mesh1%neighs_extra)))
      mesh2%neighs_extra = mesh1%neighs_extra

      ALLOCATE(mesh2%jce(SIZE(mesh1%jce, 1), SIZE(mesh1%jce, 2)))
      mesh2%jce = mesh1%jce
      ALLOCATE(mesh2%jees(SIZE(mesh1%jees)))
      mesh2%jees = mesh1%jees
      ALLOCATE(mesh2%jecs(SIZE(mesh1%jecs)))
      mesh2%jecs = mesh1%jecs

      ALLOCATE(mesh2%disp(SIZE(mesh1%disp)))
      mesh2%disp = mesh1%disp
      ALLOCATE(mesh2%domnp(SIZE(mesh1%domnp)))
      mesh2%domnp = mesh1%domnp
      ALLOCATE(mesh2%discell(SIZE(mesh1%disp)))
      mesh2%discell = mesh1%discell
      ALLOCATE(mesh2%domcell(SIZE(mesh1%domnp)))
      mesh2%domcell = mesh1%domcell
      ALLOCATE(mesh2%disedge(SIZE(mesh1%disp)))
      mesh2%disedge = mesh1%disedge
      ALLOCATE(mesh2%domedge(SIZE(mesh1%domnp)))
      mesh2%domedge = mesh1%domedge

      ALLOCATE(mesh2%jj_extra(SIZE(mesh1%jj_extra, 1), SIZE(mesh1%jj_extra, 2)))
      mesh2%jj_extra = mesh1%jj_extra
      ALLOCATE(mesh2%jce_extra(SIZE(mesh1%jce_extra, 1), SIZE(mesh1%jce_extra, 2)))
      mesh2%jce_extra = mesh1%jce_extra
      ALLOCATE(mesh2%jcc_extra(SIZE(mesh1%jcc_extra)))
      mesh2%jcc_extra = mesh1%jcc_extra

      ALLOCATE(mesh2%jjs_extra(SIZE(mesh1%jjs_extra, 1), SIZE(mesh1%jjs_extra, 2)))
      mesh2%jjs_extra = mesh1%jjs_extra
      ALLOCATE(mesh2%rrs_extra(SIZE(mesh1%rrs_extra, 1), SIZE(mesh1%rrs_extra, 2), SIZE(mesh1%rrs_extra, 3)))
      mesh2%rrs_extra = mesh1%rrs_extra

      ALLOCATE(mesh2%isolated_jjs(SIZE(mesh1%isolated_jjs)))
      mesh2%isolated_jjs = mesh1%isolated_jjs
      ALLOCATE(mesh2%isolated_interfaces(SIZE(mesh1%isolated_interfaces, 1), SIZE(mesh1%isolated_interfaces, 2)))
      mesh2%isolated_interfaces = mesh1%isolated_interfaces
   END SUBROUTINE copy_mesh

   SUBROUTINE free_mesh(mesh)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh

      DEALLOCATE(mesh%jj, mesh%i_d, mesh%loc_to_glob, mesh%rr, mesh%neigh)
      DEALLOCATE(mesh%jjs, mesh%sides, mesh%neighs)
      DEALLOCATE(mesh%jjs_int, mesh%sides_int, mesh%neighs_int)
      DEALLOCATE(mesh%disp, mesh%domnp, mesh%disedge, mesh%domedge, mesh%discell, mesh%domcell)
      DEALLOCATE(mesh%jce, mesh%jees, mesh%jecs)

      DEALLOCATE(mesh%jj_extra, mesh%jce_extra, mesh%jjs_extra, mesh%jcc_extra, mesh%rrs_extra)
      DEALLOCATE(mesh%sides_extra, mesh%neighs_extra) !interfaces
      DEALLOCATE(mesh%isolated_jjs, mesh%isolated_interfaces)

      !IF (mesh%edge_stab) THEN
      !   DEALLOCATE(mesh%iis)
      !   NULLIFY(mesh%jji)
      !   DEALLOCATE(mesh%jjsi)
      !   DEALLOCATE(mesh%neighi)
      !END IF

   END SUBROUTINE free_mesh


   SUBROUTINE clean_mesh(mesh)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(mesh_type) :: mesh

      DEALLOCATE(mesh%disp, mesh%domnp, mesh%disedge, mesh%domedge, mesh%discell, mesh%domcell)
      DEALLOCATE(mesh%jce, mesh%jees, mesh%jecs)

      DEALLOCATE(mesh%jj_extra, mesh%jce_extra, mesh%jjs_extra, mesh%jcc_extra, mesh%rrs_extra)
      DEALLOCATE(mesh%sides_extra, mesh%neighs_extra) !interfaces

      !IF (mesh%edge_stab) THEN
      !   DEALLOCATE(mesh%iis)
      !   NULLIFY(mesh%jji)
      !   DEALLOCATE(mesh%jjsi)
      !   DEALLOCATE(mesh%neighi)
      !END IF

   END SUBROUTINE clean_mesh

   SUBROUTINE free_interface(interf)
      USE def_type_mesh
      USE my_util
      IMPLICIT NONE
      TYPE(interface_type) :: interf

      interf%mes = 0
      DEALLOCATE(interf%mesh1)
      DEALLOCATE(interf%mesh2)
      DEALLOCATE(interf%jjs1)
      DEALLOCATE(interf%jjs2)
      DEALLOCATE(interf%mesh1_extra, interf%mesh2_extra, interf%jjs1_extra, interf%jjs2_extra)
   END SUBROUTINE free_interface


END MODULE mesh_tools