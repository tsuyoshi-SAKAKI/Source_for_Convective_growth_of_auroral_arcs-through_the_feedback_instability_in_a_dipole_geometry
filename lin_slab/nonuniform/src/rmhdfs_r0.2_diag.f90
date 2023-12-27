MODULE RMHDFS_diag
!-------------------------------------------------------------------------------
!
!    Diagnostics including entropy transfer analysis
!
!-------------------------------------------------------------------------------

  use RMHDFS_header
  use RMHDFS_mpienv
!  use RMHDFS_set,    only: set_close
  use RMHDFS_fld,    only: fld_calfld
  use RMHDFS_intgrl, only: intgrl_zeta
  use RMHDFS_pssn,   only: pssn_derivative_xy

  implicit none

  private

  integer, parameter :: nk = min( nx, ny )

  public   diag_cntrl


CONTAINS


!--------------------------------------
  SUBROUTINE diag_cntrl( omg, psi, dns, time, id )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: psi

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, dns

    real(kind=DP), intent(in) :: time

    integer, intent(in) :: id

    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr, lpo

    real(kind=DP), save :: tout_mag, tout_ion, tout_eng

    integer :: ix, my, iz


      call fld_calfld( omg, psi, dns, phi, cpr )


!       call dt_check2 ( time, phi, psi )     !SAKAKI 13Mar


      if( id == 0 ) then

        if ( dtout_mag /= 0._DP ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 0 )
        end if


        if ( dtout_ion /= 0._DP ) then 
          call wrt ( omg, psi, dns, phi, cpr, time, 1 )
        end if  
    
        call wrt ( omg, psi, dns, phi, cpr, time, 2 )

        tout_mag  = ( int( ( time + eps )/dtout_mag ) + 1 ) * dtout_mag
        tout_ion  = ( int( ( time + eps )/dtout_ion ) + 1 ) * dtout_ion
        tout_eng  = ( int( ( time + eps )/dtout_eng ) + 1 ) * dtout_eng


      else if( id == 1 ) then

        if ( time >= tout_mag - eps ) then
          tout_mag   = tout_mag + dtout_mag
          if ( dtout_mag /= 0._DP ) then 
            call wrt ( omg, psi, dns, phi, cpr, time, 0 )
            write( olog, * ) &
              " # Magnetospheric data output at time = ", time

          end if 
        end if

        if ( time >= tout_ion - eps ) then
          tout_ion   = tout_ion + dtout_ion
          if ( dtout_ion /= 0._DP ) then 
            call wrt ( omg, psi, dns, phi, cpr, time, 1 )
            call contnu ( omg, psi, dns, time )
            write( olog, * ) &
               " # Ionospheric data output at time = ", time
          end if

        end if

! ---

        if ( time >= tout_eng - eps ) then
          tout_eng   = tout_eng + dtout_eng
          call wrt ( omg, psi, dns, phi, cpr, time, 2 )
        end if


      else if( id == 2 ) then

          call contnu ( omg, psi, dns, time )

      end if


  END SUBROUTINE diag_cntrl


!--------------------------------------
  SUBROUTINE contnu ( omg, psi, dns, time )
!--------------------------------------

    implicit none

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    real(kind=DP), intent(in) :: time


    integer ::ix, my, iz


      rewind ocnt
      write( unit=ocnt ) time, dt, omg, psi, dns


  END SUBROUTINE contnu



!--------------------------------------
  SUBROUTINE wrt ( omg, psi, dns, phi, cpr, time, id )
!--------------------------------------

! --- arguments

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    real(kind=DP), intent(in) :: time

    integer, intent(in) :: id

    real(kind=DP), dimension(0:ny) :: phi_mode_y, psi_mode_y

    real(kind=DP) :: phi_total_e, psi_total_e

    real(kind=DP), dimension(0:ny) :: idns_mode_y, icpr_mode_y, iomg_mode_y
 
    real(kind=DP) :: idns_total_e, icpr_total_e, iomg_total_e

    integer        ::  ix, my, iz


      if( id == 0 ) then

          write( unit=omag ) time, omg, cpr 
!! debug
!        my = 1
!        write(unit=olog,*) "# time = ", time
!        do iz = -nz, nz-1
!          write(unit=olog,fmt="(1p, 6e15.7)") ( xx(ix), zz(iz), omg(ix,my,iz), cpr(ix,my,iz), ix=-nx,nx )
!          write(unit=olog,*)
!        end do
!! debug

      else if( id == 1 ) then

        if ( rankz == 0 ) then
!          write( unit=oion ) time, dns(:,:,-nz), omg(:,:,-nz), cpr(:,:,-nz)
! THW
          my = 0
          write(oion,*) "# time = ", time
          write(unit=oion,fmt="(1p, 7e15.7)") ( xx(ix), dns(ix,my,-nz), omg(ix,my,-nz), cpr(ix,my,-nz), ix=-nx,nx )
          write(oion,*)

!After          write( unit=oion ) time, xx(ix), dns(:,:,-nz)
! debug
!          write(olog,*) "# time = ", time
!          write(unit=olog,*) "# time = ", time
!Before          write(unit=olog,fmt="(1p, 7e15.7)") ( xx(ix), dns(ix,my,-nz), omg(ix,my,-nz), cpr(ix,my,-nz), ix=-nx,nx )
!          write(unit=olog,fmt="(1p, 8e15.7)") ( xx(ix), zz(nz-1), dns(ix,my,nz-1), omg(ix,my,nz-1), cpr(ix,my,nz-1), ix=-nx,nx )
!          write(olog,*)

!!new
!          write(otest,*) "# time = ", time
!!          write(unit=olog,*) "# time = ", time
!!Before          write(unit=olog,fmt="(1p, 7e15.7)") ( xx(ix), dns(ix,my,-nz), omg(ix,my,-nz), cpr(ix,my,-nz), ix=-nx,nx )
!          write(unit=otest,fmt="(1p, 8e15.7)") ( xx(ix), zz(nz-1), dns(ix,my,nz-1), omg(ix,my,nz-1), cpr(ix,my,nz-1), ix=-nx,nx )
!          write(otest,*)
!
!          write(unit=olog,*)
!!new
!          write(oll,*) "# time = ", time
!          write(unit=oll,fmt="(1p, 3e15.7)") ( xx(ix), zz(nz), ll(ix), ix=-nx,nx )
!          write(oll,*)
!
! debug
        end if

      else if( id == 2 ) then

        call mode_energy ( phi, phi_mode_y, phi_total_e )
        call mode_energy ( psi, psi_mode_y, psi_total_e )

!!!        call dt_check ( time, phi_total_e, psi_total_e )
!!!        call dt_check2 ( time, phi, psi )

        if ( rankz == 0 ) then
          call imode_energy ( dns(:,:,-nz), idns_mode_y, idns_total_e )
          call imode_energy ( cpr(:,:,-nz), icpr_mode_y, icpr_total_e )
          call imode_energy ( omg(:,:,-nz), iomg_mode_y, iomg_total_e )
        end if

        if ( rank == 0 ) then
          write( unit=omph, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, phi_total_e, phi_mode_y
          write( unit=omps, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, psi_total_e, psi_mode_y

          write( unit=oidn, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, idns_total_e, idns_mode_y
          write( unit=oicr, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, icpr_total_e, icpr_mode_y
          write( unit=oiog, fmt="(f15.8, SP, 2050ES24.15e3)" ) &
               time, iomg_total_e, iomg_mode_y
        end if

      end if


  END SUBROUTINE wrt


!--------------------------------------
  SUBROUTINE dt_check ( time, phi_total_e, psi_total_e )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    real(kind=DP), intent(in) :: phi_total_e, psi_total_e

    real(kind=DP) :: cfl, dx

      dx = xx(1) - xx(0)

      cfl = sqrt( max( phi_total_e, psi_total_e ) / ( 2._DP * lz ) ) * dt / dx

      if( cfl > 0.25_DP ) then
!!!      if( cfl > 0.1_DP ) then
        write(olog,fmt="(a,1p,e15.7,a,e15.7)") "### now, cfl = ", cfl, " at t = ", time
        write(olog,fmt="(a,1p,e15.7)")         "### then, dt is changed to dt = ", dt*0.5_DP
        dt   = dt * 0.5_DP
      end if


  END SUBROUTINE dt_check
 

!--------------------------------------
  SUBROUTINE dt_check2 ( time, phi, psi )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, psi
!!! real(kind=DP), dimension(-nz:nz-1)        :: dpara, valf, b0
    real(kind=DP), dimension(0:nny-1,0:nnx) :: vx, vy

    real(kind=DP) :: cfl, cfl_z, cfl_v, cfl_g, dx

    integer :: ix, iy, iz


!!!      dx = 2._DP * pi / kx(nx)
      dx = xx(1) - xx(0)

      cfl   = 0._DP
      cfl_z = 0._DP
      cfl_v = 0._DP

      do iz = -nz, nz-1
        cfl_z = max( cfl_z, abs(valf(iz)*dt/dpara(iz)) )
      end do

      do iz = -nz, nz-1
        call pssn_derivative_xy( phi(:,:,iz), vy, vx )
        do ix = 0, nnx-1
          do iy = 0, nny-1
            cfl_v = max( cfl_v, abs(vx(iy,ix)), abs(vx(iy,ix)) )
          end do
        end do
      end do
        cfl = max( cfl_z, cfl_v*dt/dx )

        call MPI_Allreduce( cfl, cfl_g, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )


      if( cfl_g > 0.7_DP ) then
!!!      if( cfl > 0.1_DP ) then
        write(olog,fmt="(a,1p,2e15.7,a,e15.7)") &
             "### now, cfl, cfl_g = ", cfl, cfl_g, " at t = ", time
        write(olog,fmt="(a,1p,e15.7)")         "### then, dt is changed to dt = ", dt*0.5_DP
        dt   = dt * 0.5_DP
      end if


  END SUBROUTINE dt_check2


!--------------------------------------
  SUBROUTINE mode_energy ( wrk, mode_y, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: wrk

    real(kind=DP), dimension(0:ny) :: mode_y
 
    real(kind=DP) :: total_e

! --- local variables

    real(kind=DP), dimension(:,:,:), allocatable :: wr3
    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  ix, my, iz

    integer, dimension(0:nk) :: ic

    real(kind=DP) :: dx

      dx = xx(1) - xx(0)


      allocate( wr3(-nx:nx,0:ny,-nz:nz-1) )
      allocate( wr2(-nx:nx,0:ny) )

      mode_y(:)  = 0._DP
      total_e    = 0._DP

      ic(:)      = 0


        do iz = -nz, nz-1
          do my = 0, ny
            do ix = -nx, nx
              wr3(ix,my,iz) = real( wrk(ix,my,iz) * conjg( wrk(ix,my,iz) ), kind=DP )
            end do
          end do
        end do


      call intgrl_zeta ( wr3, dpara, wr2 )


! --- ky-modes
        do my = 0, ny
          do ix = -nx, nx
            mode_y(my) = mode_y(my) + wr2(ix,my)*dx
          end do
        end do
    
    
! --- total in whole k-space
      do my = 0, ny
        total_e = total_e + mode_y(my)
      end do


      deallocate( wr3 )
      deallocate( wr2 )


  END SUBROUTINE mode_energy


!--------------------------------------
  SUBROUTINE imode_energy ( wrk, mode_y, total_e )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny) :: wrk

    real(kind=DP), dimension(0:ny) :: mode_y
 
    real(kind=DP) :: total_e

! --- local variables

    real(kind=DP), dimension(:,:),   allocatable :: wr2

    integer  ::  ix, my, iz

    integer, dimension(0:nk) :: ic

    real(kind=DP) :: dx

      dx = xx(1) - xx(0)


      allocate( wr2(-nx:nx,0:ny) )


      mode_y(:)  = 0._DP
      total_e    = 0._DP

      ic(:)      = 0

          do my = 0, ny
            do ix = -nx, nx
              wr2(ix,my) = real( wrk(ix,my) * conjg( wrk(ix,my) ), kind=DP )
            end do
          end do


! --- ky-modes
        do my = 0, ny
          do ix = -nx, nx
            mode_y(my) = mode_y(my) + wr2(ix,my)*dx
          end do
        end do
    

! --- total in whole k-space
      do my = 0, ny
        total_e = total_e + mode_y(my)
      end do


      deallocate( wr2 )


  END SUBROUTINE imode_energy


END MODULE RMHDFS_diag
