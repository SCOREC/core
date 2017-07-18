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
c
c        subroutine get_model_velocity (v)
c     &    bind(C, NAME='get_model_velocity')
c        use iso_c_binding
c          real(c_double) :: v
c        end subroutine
      end interface
c
      end module

