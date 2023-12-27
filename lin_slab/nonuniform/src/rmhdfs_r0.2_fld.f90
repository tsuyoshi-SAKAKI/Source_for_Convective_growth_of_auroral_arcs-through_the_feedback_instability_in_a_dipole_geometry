MODULE RMHDFS_fld
!-------------------------------------------------------------------------------
!
!    Time integration of the reduced MHD equations
!
!-------------------------------------------------------------------------------

  use RMHDFS_header
  use RMHDFS_mpienv
  use RMHDFS_pssn,  only: pssn_derivative, pssn_brackets, pssn_divgrad
  use RMHDFS_bndry
  use LB_matrix,    only: LB_matrix_slvtridmrc, LB_matrix_slvcyctridmrc, &
                          LB_matrix_slvperitridmrc


  implicit none

  private

  public   fld_calfld, fld_lapl


CONTAINS


!--------------------------------------
  SUBROUTINE fld_calfld( omg, psi, dns, phi, cpr )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    if( bndry_type == "periodic" ) then

      call fld_calfld_peri( omg, psi, dns, phi, cpr )

    else if( bndry_type == "zerofix" ) then

      call fld_calfld_fix( omg, psi, dns, phi, cpr )

    end if


  END SUBROUTINE fld_calfld


!--------------------------------------
  SUBROUTINE fld_calfld_peri( omg, psi, dns, phi, cpr )
!
!   Periodic boundary at x = +- Lx for phi ans psi
!
!--------------------------------------
!     coompute phi and cpr

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    complex(kind=DP), dimension(0:ny,-nz:nz-1) :: phi_sum

    real(kind=DP) :: cef
    integer :: my, iz

    real(kind=DP), dimension(-nx:nx-1) :: mtxa
    real(kind=DP), dimension(-nx:nx-1) :: mtxb
    real(kind=DP), dimension(-nx:nx-1) :: mtxc


      cef   = 1._DP / ( xx(1) - xx(0) )**2

      mtxa(:) = cef
      mtxc(:) = cef


      call fld_lapl( psi, cpr, 2*nz )



      my = 0

        mtxb(:) = - 2._DP * cef - ky(my)**2

!$OMP parallel do private(iz) shared(my)
        do iz = -nz, nz-1
          call LB_matrix_slvperitridmrc( mtxa, mtxb, mtxc, omg(-nx:nx-1,my,iz), &
                                         2*nx, phi(-nx:nx-1,my,iz) )
          phi( nx,my,iz) = phi(-nx,my,iz)
        end do


      do my = 1, ny

        mtxb(:) = - 2._DP * cef - ky(my)**2

!$OMP parallel do private(iz) shared(my)
        do iz = -nz, nz-1
          call LB_matrix_slvcyctridmrc( mtxa, mtxb, mtxc, omg(-nx:nx-1,my,iz), &
                                        2*nx, phi(-nx:nx-1,my,iz) )
          phi( nx,my,iz) = phi(-nx,my,iz)
        end do

      end do


! ionospheric boundary
      if( rankz == 0 ) then

! phi -> j
        call iono_cpr( phi(:,:,-nz), dns(:,:,-nz), cpr(:,:,-nz) )

! j -> psi
        my = 0

          mtxb(:) = - 2._DP * cef - ky(my)**2

          call LB_matrix_slvperitridmrc( mtxa, mtxb, mtxc, cpr(-nx:nx-1,my,-nz), &
                                         2*nx, psi(-nx:nx-1,my,-nz) )
          psi( nx,my,-nz) = psi(-nx,my,-nz)


        do my = 1, ny

          mtxb(:) = - 2._DP * cef - ky(my)**2

          call LB_matrix_slvcyctridmrc( mtxa, mtxb, mtxc, cpr(-nx:nx-1,my,-nz), &
                                        2*nx, psi(-nx:nx-1,my,-nz) )
          psi( nx,my,-nz) = psi(-nx,my,-nz)
        end do


      end if



  END SUBROUTINE fld_calfld_peri


!--------------------------------------
  SUBROUTINE fld_calfld_fix( omg, psi, dns, phi, cpr )
!
!   Zero fixed boundary at x = +- Lx for phi ans psi
!
!      But, condision for psi may need to be changed...
!      open boundary with omg = cpr = 0 for x<-Lx or x>Lx 
!      should be used => See Birdsall and Langdon
!
!--------------------------------------
!     coompute phi and cpr

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    real(kind=DP) :: cef
    integer :: my, iz

    real(kind=DP), dimension(-nx+2:nx-1) :: mtxa
    real(kind=DP), dimension(-nx+1:nx-1) :: mtxb
    real(kind=DP), dimension(-nx+1:nx-2) :: mtxc



      cef   = 1._DP / ( xx(1) - xx(0) )**2

      mtxa(:) = cef
      mtxc(:) = cef

      call fld_lapl( psi, cpr, 2*nz )

      do my = 0, ny

        mtxb(:) = - 2._DP * cef - ky(my)**2

!$OMP parallel do private(iz) shared(my)
        do iz = -nz, nz-1

          call LB_matrix_slvtridmrc( mtxa, mtxb, mtxc, omg(-nx+1:nx-1,my,iz), &
                                     2*nx-1, phi(-nx+1:nx-1,my,iz) )

          phi(-nx,my,iz) = ( 0._DP, 0._DP )
          phi( nx,my,iz) = ( 0._DP, 0._DP )

        end do
      end do


! ionospheric boundary
      if( rankz == 0 ) then

! phi -> j
        call iono_cpr( phi(:,:,-nz), dns(:,:,-nz), cpr(:,:,-nz) )

! j -> psi
        do my = 0, ny

          mtxb(:) = - 2._DP * cef - ky(my)**2

          call LB_matrix_slvtridmrc( mtxa, mtxb, mtxc, cpr(-nx+1:nx-1,my,-nz), &
                                     2*nx-1, psi(-nx+1:nx-1,my,-nz) )

          psi(-nx,my,-nz) = ( 0._DP, 0._DP )
          psi( nx,my,-nz) = ( 0._DP, 0._DP )

        end do


      end if


  END SUBROUTINE fld_calfld_fix


!--------------------------------------
  SUBROUTINE iono_cpr( iphi, idns, icpr )
!--------------------------------------
!    compute the ionospheric potential

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny)  :: iphi, idns
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)  :: icpr

    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: ipsnb, idvgd, lapl_iphi, lapl_idns

    complex(kind=DP), &
      dimension(-nx:nx,0:ny) :: igrdnx, igrdny

    integer :: ix, my


      if( trim(calc_type) == "nonlinear" ) then
        call pssn_brackets( iphi, idns, ipsnb, 1 )
        call pssn_divgrad ( idns, iphi, idvgd, 1 )
      else
        ipsnb(:,:) = ( 0._DP, 0._DP )
        idvgd(:,:) = ( 0._DP, 0._DP )
      end if
      
      call fld_lapl( iphi, lapl_iphi, 1 )
      call fld_lapl( idns, lapl_idns, 1 )
      call pssn_derivative( idns, igrdnx, igrdny )

!$OMP parallel do collapse(2) private(my,ix)
      do my = 0, ny
        do ix = -nx, nx
          icpr(ix,my) = mp * idns0 * lapl_iphi(ix,my)                    &
                      - e0 * ( ( mp * cos(theta2) - mh * sin(theta2) )   &
                      * igrdnx(ix,my) + ( mp * sin(theta2) + mh          &
                      * cos(theta2) ) * igrdny(ix,my) )                  &
                      + dperp * lapl_idns(ix,my)                         & 
                      + mp * idvgd(ix,my) + ipsnb(ix,my)
        end do
      end do

     
  END SUBROUTINE iono_cpr


!--------------------------------------
  SUBROUTINE fld_lapl( wrk, lapl_wrk, nm )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: wrk
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: lapl_wrk
    integer, intent(in) :: nm


    if( bndry_type == "periodic" ) then

      call fld_lapl_peri( wrk, lapl_wrk, nm )

    else if( bndry_type == "zerofix" ) then

      call fld_lapl_fix( wrk, lapl_wrk, nm )

    end if


  END SUBROUTINE fld_lapl


!--------------------------------------
  SUBROUTINE fld_lapl_peri( wrk, lapl_wrk, nm )
!
!   Periodic boundary at x = +- Lx
!
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: wrk
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: lapl_wrk
    integer, intent(in) :: nm

    real(kind=DP) :: cef
    integer :: ix, my, iz


      cef   = 1._DP / ( xx(1) - xx(0) )**2

! 2nd-order in x
!$OMP parallel do private(iz,my,ix)
      do iz = 0, nm-1
        do my = 0, ny
          do ix = -nx+1, nx-1
            lapl_wrk(ix,my,iz)  = cef * ( wrk(ix+1,my,iz) - 2._DP*wrk(ix,my,iz) + wrk(ix-1,my,iz) ) &
                                - ky(my)**2 * wrk(ix,my,iz)
          end do
            lapl_wrk(-nx,my,iz) = cef * ( wrk(-nx+1,my,iz) - 2._DP*wrk(-nx,my,iz) + wrk(nx-1,my,iz) ) &
                                - ky(my)**2 * wrk(-nx,my,iz)
            lapl_wrk( nx,my,iz) = lapl_wrk(-nx,my,iz)
        end do
      end do

  END SUBROUTINE fld_lapl_peri


!--------------------------------------
  SUBROUTINE fld_lapl_fix( wrk, lapl_wrk, nm )
!
!   Zero fixed boundary at x = +- Lx
!
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: wrk
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,0:nm-1)  :: lapl_wrk
    integer, intent(in) :: nm

    real(kind=DP) :: cef
    integer :: ix, my, iz


      cef   = 1._DP / ( xx(1) - xx(0) )**2

! 2nd-order in x
!$OMP parallel do private(iz,my,ix)
      do iz = 0, nm-1
        do my = 0, ny
          do ix = -nx+1, nx-1
            lapl_wrk(ix,my,iz)  = cef * ( wrk(ix+1,my,iz) - 2._DP*wrk(ix,my,iz) + wrk(ix-1,my,iz) ) &
                                - ky(my)**2 * wrk(ix,my,iz)
          end do
            lapl_wrk(-nx,my,iz) = ( 0._DP, 0._DP )
            lapl_wrk( nx,my,iz) = ( 0._DP, 0._DP )
        end do
      end do

  END SUBROUTINE fld_lapl_fix


END MODULE RMHDFS_fld
