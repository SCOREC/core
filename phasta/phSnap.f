      module simmetrix_snap
      use iso_c_binding
c
      interface
c
        subroutine sim_get_pos_on_surf ( dx, dy, dz, numnp, px, py, pz )
     &    bind(C, NAME='sim_get_pos_on_surf')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &                               dx(numnp), dy(numnp), dz(numnp)
          real(c_double),intent(out),dimension(:) ::
     &                               px(numnp), py(numnp), pz(numnp)
        end subroutine
c
        subroutine sim_is_in_closure ( e_dim, e_tag,
     &                                 t_dim, t_tag, answer )
     &    bind(C, NAME='sim_is_in_closure')
        use iso_c_binding
          integer(c_int),intent(in),value :: e_dim, e_tag, t_dim, t_tag
          integer(c_int),intent(out) :: answer
        end subroutine
c
      end interface
c
      end module
