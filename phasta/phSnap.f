      module simmetrix_snap
      use iso_c_binding
c
      interface
        subroutine sim_get_pos_on_surf ( dx, dy, dz, id )
     &    bind(C, NAME='sim_get_pos_on_surf')
        use iso_c_binding
          real(c_double), value :: dx, dy, dz
          integer(c_int), value :: id
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

