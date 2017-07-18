      module simmetrix_snap
      use iso_c_binding
c
      interface
        subroutine sim_get_pos_on_surf ( dx, dy, dz, numnp, px, py, pz )
     &    bind(C, NAME='sim_get_pos_on_surf')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(out),dimension(:) ::
     &                               dx(numnp), dy(numnp), dz(numnp)
          real(c_double),intent(in),dimension(:) ::
     &                               px(numnp), py(numnp), pz(numnp)
        end subroutine
      end interface
      end module
