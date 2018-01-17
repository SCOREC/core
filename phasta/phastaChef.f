      module core_snap
      use iso_c_binding
      interface

        subroutine core_get_pos_on_surf (dx, dy, dz, numnp, px, py, pz) &
          bind(C, NAME='core_get_pos_on_surf')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) :: &
                                    dx(numnp), dy(numnp), dz(numnp)
          real(c_double),intent(out),dimension(:) :: &
                                    px(numnp), py(numnp), pz(numnp)
        end subroutine

        subroutine core_is_in_closure ( e_dim, e_tag, &
                                       t_dim, t_tag, answer ) &
          bind(C, NAME='core_is_in_closure')
        use iso_c_binding
          integer(c_int),intent(in),value :: e_dim, e_tag, t_dim, t_tag
          integer(c_int),intent(out) :: answer
        end subroutine

      end interface
      end module


      module core_mesh_quality
      use iso_c_binding
      interface

        subroutine core_measure_mesh(x1, x2, x3, numnp, minq) &
          bind(C, NAME='core_measure_mesh')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) :: &
                                     x1(numnp), x2(numnp), x3(numnp)
          real(c_double),intent(out) :: minq
        end subroutine

      end interface
      end module
