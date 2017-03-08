      module phiotimer
      use :: iso_c_binding
      public
        enum, bind(C)
          enumerator :: CHEF_GEOMBC, CHEF_RESTART,
     &                  PHASTA_GEOMBC, PHASTA_RESTART,
     &                  NUM_PHASTA_FILES
        end enum
      interface 
        subroutine phastaio_setfile(fileidx)
     &   bind(C, NAME='phastaio_setfile')
        use :: iso_c_binding
          integer(c_int), intent(in), value:: fileidx
        end subroutine
      end interface
      end module
