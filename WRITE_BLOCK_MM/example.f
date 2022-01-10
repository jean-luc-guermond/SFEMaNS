      program example
        !use iso_c_binding, only: c_loc
        use zlib_base64, only: write_compressed_block
        implicit none

        INTEGER(KIND=4) :: number_of_points
        INTEGER(KIND=4) :: number_of_cells
        REAL(KIND=4), POINTER :: points(:), solution(:)
        INTEGER(KIND=4), POINTER :: connectivity(:), offsets(:)
        INTEGER(KIND=1), POINTER :: types(:)
        INTEGER :: unit_file

        number_of_points = 16
        number_of_cells = 4

        allocate(points(3 * number_of_points))
        points = (/-1.,-1.,0.,0.,-1.,0.,-1.,0.,0.,0.,0.,0.,0.,-1.,0.,1.,-1.,0.,0.,0.,0.,1.,0.,0.,-1.,0.,0.,0.,0.,0.,-1.,1.,0.,&
          0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0./)

        allocate(connectivity(4 * number_of_cells))
        connectivity = (/0,1,3,2,4,5,7,6,8,9,11,10,12,13,15,14/)

        allocate(offsets(number_of_cells))
        offsets = (/4,8,12,16/)

        allocate(types(number_of_cells))
        types = (/9,9,9,9/)

        allocate(solution(number_of_points))
        solution = (/0.,0.,0.,0.375,0.,0.,0.375,0.,0.,0.375,0.,0., &
          0.375,0.,0.,0./)

        unit_file = 42
        open(unit = unit_file, file = "solution.vtk")

        write(42,'(A)') '<?xml version="1.0" ?>'
        write(42,'(A)', advance="no") '<VTKFile type="UnstructuredGrid" version="0.1" '
        write(42,'(A)') 'compressor="vtkZLibDataCompressor" byte_order="LittleEndian">'
        write(42,'(A)') '<UnstructuredGrid>'
        write(42,'(A, I0, A, I0, A)') '<Piece NumberOfPoints="', number_of_points, '" NumberOfCells="', number_of_cells, '">'
        write(42,'(A)') '<Points>'
        write(42,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="binary">'

        call write_compressed_block(unit_file, points)

        write(42,'(A)') '</DataArray>'
        write(42,'(A)') '</Points>'
        write(42,'(A)') '<Cells>'
        write(42,'(A)') '<DataArray type="Int32" Name="connectivity" format="binary">'

        call write_compressed_block(unit_file, connectivity)

        write(42,'(A)') '</DataArray>'
        write(42,'(A)') '<DataArray type="Int32" Name="offsets" format="binary">'

        call write_compressed_block(unit_file, offsets)

        write(42,'(A)') '</DataArray>'
        write(42,'(A)') '<DataArray type="UInt8" Name="types" format="binary">'

        call write_compressed_block(unit_file, types)

        write(42,'(A)') '</DataArray>'
        write(42,'(A)') '</Cells>'
        write(42,'(A)') '<PointData Scalars="scalars">'
        write(42,'(A)') '<DataArray type="Float32" Name="solution" format="binary">'

        call write_compressed_block(unit_file, solution)

        write(42,'(A)') '</DataArray>'
        write(42,'(A)') '</PointData>'
        write(42,'(A)') '</Piece>'
        write(42,'(A)') '</UnstructuredGrid>'
        write(42,'(A)') '</VTKFile>'

        close(unit_file)
      end program example
