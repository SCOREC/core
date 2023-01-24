      module core_snap
      use iso_c_binding
c
      interface
c
        subroutine core_get_pos_on_surf (dx,dy,dz, numnp, f, px,py,pz)
     &    bind(C, NAME='core_get_pos_on_surf')
        use iso_c_binding
          integer(c_int),value :: numnp
          integer(c_int),intent(in),dimension(:) ::  f(numnp)
          real(c_double),intent(in),dimension(:) ::
     &                               dx(numnp), dy(numnp), dz(numnp)
          real(c_double),intent(out),dimension(:) ::
     &                               px(numnp), py(numnp), pz(numnp)
        end subroutine
c
        subroutine core_is_in_closure ( e_dim, e_tag,
     &                                  t_dim, t_tag, answer )
     &    bind(C, NAME='core_is_in_closure')
        use iso_c_binding
          integer(c_int),intent(in),value :: e_dim, e_tag, t_dim, t_tag
          integer(c_int),intent(out) :: answer
        end subroutine
c
        subroutine core_get_surf_normal (flag,f_tag, e_or_v_tag, par1,
     &           par2, n1, n2, n3)
     &    bind(C, NAME='core_get_surf_normal')
        use iso_c_binding
          integer(c_int),intent(in),value :: flag,f_tag, e_or_v_tag
          real(c_double),intent(in),value :: par1, par2
          real(c_double),intent(out) :: n1, n2, n3
        end subroutine
c
c
      end interface
c
      end module
c-------------------------------------------------------------------------
c
c-------------------------------------------------------------------------

      module core_ScalarSmoother
      use iso_c_binding

      interface

        subroutine core_ShockVelSmoother(x1, x2, x3, numnp, V_i, V_o)
     &    bind(C, NAME='core_ShockVelSmoother')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &         x1(numnp), x2(numnp), x3(numnp), V_i(numnp)
          real(c_double),intent(out) :: V_o(numnp)
        end subroutine
      end interface

      end module


c-------------------------------------------------------------------------
c
c-------------------------------------------------------------------------
      module core_mesh_quality
      use iso_c_binding
c
      interface
c
        subroutine core_measure_mesh(x1, x2, x3, numnp, minvq, minfq)
     &    bind(C, NAME='core_measure_mesh')
        use iso_c_binding
          integer(c_int),value :: numnp
          real(c_double),intent(in),dimension(:) ::
     &                               x1(numnp), x2(numnp), x3(numnp)
          real(c_double),intent(out) :: minvq, minfq
        end subroutine
c
      end interface
c
      end module
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      module core_rigid_body
      use iso_c_binding
c
      interface
c
        subroutine core_get_centroid(r_tag, ct)
     &    bind(C, NAME='core_get_centroid')
        use iso_c_binding
          integer(c_int),intent(in) :: r_tag
          real(c_double),intent(out),dimension(:) :: ct(3)
        end subroutine
c
        subroutine core_update_rbms(
     &                   tx, ty, tz,
     &                   ax, ay, az,
     &                   px, py, pz,
     &                   ag, sc, tags, numRbm)
     &    bind(C, NAME='core_update_rbms')
        use iso_c_binding
          integer(c_int),value :: numRbm
          real(c_double),intent(in),dimension(:) ::
     &                               tx(numRbm), ty(numRbm), tz(numRBM),
     &                               ax(numRbm), ay(numRbm), az(numRBM),
     &                               px(numRbm), py(numRbm), pz(numRBM),
     &                               ag(numRbm), sc(numRbm)
          integer(c_int),intent(in),dimension(:) :: tags(numRbm)
        end subroutine
c
      end interface
c
      end module

