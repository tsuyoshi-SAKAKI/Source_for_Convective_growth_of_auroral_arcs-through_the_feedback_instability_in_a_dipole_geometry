MODULE RMHDFS_bndry
!-------------------------------------------------------------------------------
!
!    Some useful tools and tips
!
!      communication in z-direction
!      zero-fixed boundary at x=-lx and lx
!
!-------------------------------------------------------------------------------

  use RMHDFS_header
  use RMHDFS_mpienv

  implicit none

  private

  public   bndry_bound, bndry_xbuf


CONTAINS


!--------------------------------------
  SUBROUTINE bndry_bound ( ww )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nb:nz-1+nb)   :: ww

    complex(kind=DP), dimension(:,:,:), allocatable :: zb1, zb2

! --- local variables

    integer  ::  mx, my, ib
    integer  ::  slngz


      slngz  = (2*nx+1)*(ny+1)*nb


      allocate( zb1(-nx:nx,0:ny,nb*2) )
      allocate( zb2(-nx:nx,0:ny,nb*2) )


      do ib = 1, nb
        do my = 0, ny
          do mx = -nx, nx
            zb1(mx,my,ib)    = ww(mx,my,-nz   +ib-1)
            zb1(mx,my,nb+ib) = ww(mx,my, nz-nb+ib-1)
          end do
        end do
      end do

      call MPI_sendrecv( zb1(-nx,0,1),    slngz, MPI_DOUBLE_COMPLEX, izdn, 1, &
                         zb2(-nx,0,nb+1), slngz, MPI_DOUBLE_COMPLEX, izup, 1, &
                         MPI_COMM_WORLD, status, ierr_mpi )

      call MPI_sendrecv( zb1(-nx,0,nb+1), slngz, MPI_DOUBLE_COMPLEX, izup, 2, &
                         zb2(-nx,0,1),    slngz, MPI_DOUBLE_COMPLEX, izdn, 2, &
                         MPI_COMM_WORLD, status, ierr_mpi )


      if( rankz /= 0 ) then

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,-nz-nb+ib-1) = zb2(mx,my,ib)
            end do
          end do
        end do

      else

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,-nz-nb+ib-1) = ( 0._DP, 0._DP )
            end do
          end do
        end do

      end if


      if( rankz /= nprocz-1 ) then

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my, nz+ib-1) = zb2(mx,my,nb+ib)
            end do
          end do
        end do

      else

        do ib = 1, nb
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my, nz+ib-1) = ( 0._DP, 0._DP )
            end do
          end do
        end do

      end if



  END SUBROUTINE bndry_bound


!--------------------------------------
  SUBROUTINE bndry_xbuf ( ww, wu )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)   :: ww
    complex(kind=DP), intent(inout), &
      dimension(-nx-nb:nx+nb,0:ny,-nz:nz-1)   :: wu

    integer :: ix, my, iz, ib


      do iz = -nz, nz-1
        do my = 0, ny
          do ix = -nx, nx
            wu(ix,my,iz) = ww(ix,my,iz)
          end do
          do ib = 1, nb
            wu(-nx-ib,my,iz) = ( 0._DP, 0._DP )
            wu( nx+ib,my,iz) = ( 0._DP, 0._DP )
          end do
        end do
      end do


  END SUBROUTINE bndry_xbuf


END MODULE RMHDFS_bndry
