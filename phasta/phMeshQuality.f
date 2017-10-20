      module core_mesh_quality
      use iso_c_binding
c
      interface
c
        subroutine core_measure_mesh(x1, x2, x3, numnp, minq)
     &    bind(C, NAME='core_measure_mesh')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &                               x1(numnp), x2(numnp), x3(numnp)
          real(c_double),intent(out) :: minq
        end subroutine
c
      end interface
c
      end module
