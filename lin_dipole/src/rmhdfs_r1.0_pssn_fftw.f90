MODULE RMHDFS_pssn
!-------------------------------------------------------------------------------
!
!    Calculation of Poisson bracket
!
!      created by Maeyama
!      modified by THW for aurora spectral code with FD in x
!
!      In this version, nnx=2*nx so that
!      ww(-nx:nx,my) = wwkk(my,0:nnx)
!
!      FFTW:   1D FFT
!      OpenMP: mx,ky-loop parallelization
!
!-------------------------------------------------------------------------------
  use RMHDFS_header
  use RMHDFS_mpienv

  implicit none
  include "fftw3.f"
  private

  public  pssn_brackets, pssn_divgrad, pssn_derivative, pssn_derivative_xy

  integer(kind=DP), save :: plan_backward_x, plan_forward_x
  integer(kind=DP), save :: plan_backward_y, plan_forward_y
  complex(kind=DP), dimension(0:nny/2,0:nnx) :: wwkk

  integer, save :: iflg
  data iflg / 0 /


 CONTAINS


SUBROUTINE pssn_brackets( fk, pk, pb, nm )
!-------------------------------------------------------------------------------
!
!    Calculate time-differential term
!
!-------------------------------------------------------------------------------
  use RMHDFS_header

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: pk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(out) :: pb
  integer, intent(in) :: nm

  real(kind=DP), dimension(0:nny-1,0:nnx) :: dfdx, dfdy, dpdx, dpdy, pbxy

  integer :: mx, my, ix, iy, im ,iz


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


!%%% Calculate Poisson bracket %%%
!%%% !$OMP parallel default(none) &
!$OMP parallel &
!$OMP shared(fk,pk,dfdx,dfdy,dpdx,dpdy,pbxy,pb,nm) &
!$OMP private(mx,my,im,ix,iy)
    do im = 0, nm-1

      call pssn_derivative_xy( fk(:,:,im), im, dfdx, dfdy )
      call pssn_derivative_xy( pk(:,:,im), im, dpdx, dpdy )

!==============================SAKAKI===============================
!$OMP do
      do ix = 0, nnx
        do iy = 0, nny-1
          pbxy(iy,ix) =  - dpdx(iy,ix) * dfdy(iy,ix) + dpdy(iy,ix) * dfdx(iy,ix) 
        end do
      end do
!$OMP end do
      call exb_fft_forward(pbxy, pb(:,:,im))
    end do
!$OMP end parallel
!==============================SAKAKI===============================

END SUBROUTINE pssn_brackets


SUBROUTINE pssn_derivative_xy( gk, iz, dgdx, dgdy )
!-------------------------------------------------------------------------------
!   Compute spatial derivative of g in x and y space
!-------------------------------------------------------------------------------

  use RMHDFS_header

  implicit none

  complex(kind=DP), dimension(-nx:nx,0:ny), intent(in)  :: gk
  integer, intent(in) :: iz

  real(kind=DP), dimension(0:nny-1,0:nnx), intent(out)  :: dgdx, dgdy

  complex(kind=DP), dimension(-nx:nx,0:ny) :: ikxg, ikyg

  integer :: mx, my, ix, iy


!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if


    call pssn_derivative( gk, iz, ikxg, ikyg )
    call exb_fft_backward( ikxg, dgdx )
    call exb_fft_backward( ikyg, dgdy )


END SUBROUTINE pssn_derivative_xy


SUBROUTINE pssn_derivative( gk, iz, ikxg, ikyg )
!-------------------------------------------------------------------------------
!   Compute spatial derivative of g in x and ky space
!-------------------------------------------------------------------------------
! finite ky

  use RMHDFS_header

  implicit none

  complex(kind=DP), dimension(-nx:nx,0:ny), intent(in)  :: gk
  integer, intent(in) :: iz

  complex(kind=DP), dimension(-nx:nx,0:ny), intent(out) :: ikxg, ikyg

  integer :: mx, my, ix, iy

  real(kind=DP) :: cef
 


      cef   = 1._DP / ( 12._DP * ( xi1(1) - xi1(0) ) )


!$OMP do
      do my = 0, ny

! spectral in xi2 but FD in xi1

        do ix = -nx+2, nx-2
          ikxg(ix,my) =   ( -gk(ix+2,my)+gk(ix-2,my) + 8._DP*(gk(ix+1,my)-gk(ix-1,my)) )*cef
        end do

        if( bndry_type == "periodic" ) then
          ix = -nx
          ikxg(ix,my) = ( -gk(ix+2,my)+gk(nx-2,my) + 8._DP*(gk(ix+1,my)-gk(nx-1,my)) )*cef
          ix = -nx+1
          ikxg(ix,my) = ( -gk(ix+2,my)+gk(nx-1,my) + 8._DP*(gk(ix+1,my)-gk(ix-1,my)) )*cef
          ix =  nx-1
          ikxg(ix,my) = ( -gk(-nx+1,my)+gk(ix-2,my) + 8._DP*(gk(ix+1,my)-gk(ix-1,my)) )*cef
          ikxg(nx,my) = ikxg(-nx,my)

        else if( bndry_type == "zerofix" ) then

          ix = -nx
          ikxg(ix,my) = ( 0._DP, 0._DP )

          ix = -nx+1
          ikxg(ix,my) = ( gk(ix+1,my) - gk(ix-1,my) ) * cef * 6._DP

          ix =  nx-1
          ikxg(ix,my) = ( gk(ix+1,my) - gk(ix-1,my) ) * cef * 6._DP

          ix =  nx
          ikxg(ix,my) = ( 0._DP, 0._DP )

        end if

        do ix = -nx, nx
          ikxg(ix,my) =  ikxg(ix,my)             / h1(ix,iz-nz)       ! SAKAKI
          ikyg(ix,my) =  ui * ky(my) * gk(ix,my) / h2(ix,iz-nz)       ! SAKAKI
        end do

      end do
!$OMP end do


END SUBROUTINE pssn_derivative


SUBROUTINE pssn_divgrad( fk, pk, dg, nm )
!-------------------------------------------------------------------------------
!
!    Calculate time-differential term
!
!-------------------------------------------------------------------------------
  use RMHDFS_header

  implicit none
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: fk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(in) :: pk
  complex(kind=DP), dimension(-nx:nx,0:ny,0:nm-1), intent(out) :: dg
  integer, intent(in) :: nm

  complex(kind=DP), dimension(-nx:nx,0:ny) :: ikxg, ikyg, ikxp, ikyp
  real(kind=DP), dimension(0:nny-1,0:nnx) :: fdpdx, fdpdy, dpdx, dpdy, reff
  integer :: mx, my, ix, iy, im

  real(kind=DP) :: cef

!%%% Initialize FFTW %%%
    if (iflg == 0) then
      iflg = 1
      call exb_fft_pre
    end if

      cef   = 1._DP / ( 12._DP * ( xi1(1) - xi1(0) ) )



!%%% Calculate Poisson bracket %%%
!%%% !$OMP parallel default(none) &
!$OMP parallel &
!$OMP shared(ky,pk,ikxg,ikyg,fdpdx,fdpdy,reff,dg,nm) &
!$OMP private(mx,my,im,ix,iy)
    do im = 0, nm-1

      call exb_fft_backward( fk(:,:,im), reff)

      call pssn_derivative ( pk(:,:,im), im, ikxp, ikyp )
      call exb_fft_backward( ikxp, dpdx )
      call exb_fft_backward( ikyp, dpdy )


!$OMP do
      do ix = 0, nnx
        do iy = 0, nny-1
          fdpdx(iy,ix) = reff(iy,ix) * dpdx(iy,ix) * cef_d1(ix-nx,im-nz) / h3(ix-nx,im-nz)
          fdpdy(iy,ix) = reff(iy,ix) * dpdy(iy,ix) * cef_d2(ix-nx,im-nz) / h3(ix-nx,im-nz)
!          fdpdx(iy,ix) = reff(iy,ix) * dpdx(iy,ix) * cef_d1(ix-nx,im-nz)
!          fdpdy(iy,ix) = reff(iy,ix) * dpdy(iy,ix) * cef_d2(ix-nx,im-nz)
        end do
      end do
!$OMP end do
      call exb_fft_forward(fdpdx, ikxg)
      call exb_fft_forward(fdpdy, ikyg)

!$OMP do
      do my = 0, ny

! spectral in y but FD in x

        do ix = -nx+2, nx-2
          dg(ix,my,im) = ( -ikxg(ix+2,my)+ikxg(ix-2,my) + 8._DP*(ikxg(ix+1,my)-ikxg(ix-1,my)) )*cef
        end do

        if( bndry_type == "periodic" ) then

          ix = -nx
          dg(ix,my,im) = ( -ikxg(ix+2,my)+ikxg(nx-2,my) + 8._DP*(ikxg(ix+1,my)-ikxg(nx-1,my)) )*cef
          ix = -nx+1
          dg(ix,my,im) = ( -ikxg(ix+2,my)+ikxg(nx-1,my) + 8._DP*(ikxg(ix+1,my)-ikxg(ix-1,my)) )*cef
          ix =  nx-1
          dg(ix,my,im) = ( -ikxg(-nx+1,my)+ikxg(ix-2,my) + 8._DP*(ikxg(ix+1,my)-ikxg(ix-1,my)) )*cef
          ix =  nx
          dg(ix,my,im) = dg(-nx,my,im)


        else if( bndry_type == "zerofix" ) then


          ix = -nx
          dg(ix,my,im) = ( 0._DP, 0._DP )

          ix = -nx+1
          dg(ix,my,im) = ( ikxg(ix+1,my) - ikxg(ix-1,my) ) * cef * 6._DP

          ix =  nx-1
          dg(ix,my,im) = ( ikxg(ix+1,my) - ikxg(ix-1,my) ) * cef * 6._DP

          ix =  nx
          dg(ix,my,im) = ( 0._DP, 0._DP )


        end if

        do ix = -nx, nx
          dg(ix,my,im) =  ( dg(ix,my,im) + ui * ky(my) * ikyg(ix,my) ) / jcb(ix,im-nz) * h3(ix,im-nz) ! only used in ionosheric dns
        end do


      end do
!$OMP end do

    end do
!$OMP end parallel


END SUBROUTINE pssn_divgrad


SUBROUTINE exb_fft_pre
!--------------------------------------
!  Initialization of FFT
  complex(kind=DP), dimension(0:nnx-1) :: wkx1, wkx2
!  complex(kind=DP), dimension(0:ny/2) :: wky1
  complex(kind=DP), dimension(0:nny/2) :: wky1
  real(kind=DP), dimension(0:nny-1) :: wky2

    call dfftw_plan_dft_1d(plan_backward_x,  &
                           nnx,               &
                           wkx1, wkx2,       &
                           FFTW_BACKWARD,    &
                           FFTW_MEASURE)
    call dfftw_plan_dft_c2r_1d(plan_backward_y,  &
                               nny,               &
                               wky1, wky2,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_r2c_1d(plan_forward_y,   &
                               nny,               &
                               wky2, wky1,       &
                               FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_forward_x,   &
                           nnx,               &
                           wkx2, wkx1,       &
                           FFTW_FORWARD,     &
                           FFTW_MEASURE)

END SUBROUTINE exb_fft_pre
  

SUBROUTINE exb_fft_backward(ww, wwxy)
!--------------------------------------
!  Execution of FFT
  complex(kind=DP), dimension(-nx:nx,0:ny), intent(in) :: ww
  real(kind=DP), dimension(0:nny-1,0:nnx), intent(out)    :: wwxy

  complex(kind=DP), dimension(0:nnx-1) :: w1, w2
  integer :: mx, my

    w1(:) = (0._DP, 0._DP)
!$OMP workshare
    wwkk(:,:) = (0._DP, 0._DP)
!$OMP end workshare
!$OMP do
! spectral in xi2 but FD in xi1
    do my = 0, ny
      wwkk(my,0:nnx) = ww(-nx:nx,my)
    end do
!$OMP end do

!$OMP do
!    do mx = 0, nnx-1
    do mx = 0, nnx 
      call dfftw_execute_dft_c2r(plan_backward_y, wwkk(:,mx), wwxy(:,mx))
    end do
!$OMP end do

END SUBROUTINE exb_fft_backward


SUBROUTINE exb_fft_forward(wwxy, ww)
!--------------------------------------
!  Execution of FFT
  real(kind=DP), dimension(0:nny-1,0:nnx), intent(in)      :: wwxy
  complex(kind=DP), dimension(-nx:nx,0:ny), intent(out) :: ww

  complex(kind=DP), dimension(0:nnx-1) :: w1, w2
  real(kind=DP) :: ceff_norm
  integer :: mx, my

! spectral in xi2 but FD in xi1
    ceff_norm = 1._DP / real(nny, kind=DP)

!$OMP do
    do mx = 0, nnx
      call dfftw_execute_dft_r2c(plan_forward_y, wwxy(:,mx), wwkk(:,mx))
    end do
!$OMP end do
!$OMP do
! spectral in xi2 but FD in xi1

    do my = 0, ny         
      ww(-nx:nx,my) = wwkk(my,0:nnx) * ceff_norm
    end do

!$OMP end do
    my = 0  
      do mx = -nx, nx
      ww(mx,my) = cmplx( real(ww(mx,my)), 0._DP )
      end do

END SUBROUTINE exb_fft_forward


END MODULE RMHDFS_pssn
