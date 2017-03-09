      module phiotimer
      use :: iso_c_binding
      public
        enum, bind(C)
          enumerator :: GEOMBC_READ, GEOMBC_WRITE,
     &                  RESTART_READ, RESTART_WRITE,
     &                  NUM_PHASTAIO_MODES
        end enum
      interface 
        subroutine phastaio_setfile(fileidx)
     &   bind(C, NAME='phastaio_setfile')
        use :: iso_c_binding
          integer(c_int), intent(in), value:: fileidx
        end subroutine
      end interface
      end module
