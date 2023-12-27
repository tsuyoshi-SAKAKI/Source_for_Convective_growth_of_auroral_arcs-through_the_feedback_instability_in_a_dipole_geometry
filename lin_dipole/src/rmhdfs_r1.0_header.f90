MODULE RMHDFS_header
!-------------------------------------------------------------------------------
!
!    Nonlinear fluxtube reduced MHD FD and spectral code RMHDFS
!      - Finite-Difference in x and spectral in y
!
!      Header for feedback instability
!
!      RMHDF r0.1 ( T.-H.Watanabe, April 2019)
!
!      Coding style is based on GKV-plus r2.1
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)

!--------------------------------------
!  Dimension size (grid numbers)
!--------------------------------------

  integer, parameter :: nxw = 1280, nyw = 1
  integer, parameter :: global_nz = 64
  integer, parameter :: nb = 2

!--------------------------------------
!  Data distribution for MPI
!--------------------------------------

  integer, parameter :: nprocz = 32

!--------------------------------------
!  Parameters for variable sizes
!--------------------------------------

  integer, parameter :: nx = nxw, ny = (nyw/3)*2  

  integer, parameter :: nz = global_nz / nprocz

  integer, parameter :: nxyz  = (2*nx+1)*(ny+1)*(2*nz),    &
                        nxy   = (2*nx+1)*(ny+1)     

  integer, parameter :: nnx = nxw*2, nny = nyw*2

!--------------------------------------
!  Constants
!--------------------------------------

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, &
                                 twopi = pi * 2._DP,         &
                                 eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )
 
  real(kind=DP), parameter :: alp = 1000._DP
  real(kind=DP), parameter :: theta0 = pi/9._DP      ! footpoint of the box center
  real(kind=DP), parameter :: h2_ref = sin(theta0)

!--------------------------------------
!  Parameters for time
!--------------------------------------

  real(kind=DP) :: e_limit                           ! elapsed time limit of a job
  real(kind=DP) :: tend                              ! end time
  real(kind=DP) :: dtout_mag, dtout_ion, dtout_eng   ! time-spacing for output


!--------------------------------------
!  Configuration parameters to be 
!    initialized in init subroutine
!--------------------------------------
  real(kind=DP), dimension(-nx:nx , -nz:nz-1)  ::  ig11, ig22, ig33, g33, cef_d1, cef_d2, jcb, ig13
  real(kind=DP), dimension(-nx:nx , -nz:nz-1)  ::  h1, h2, h3 , dcx, dcz, radius_grid,theta_grid
  real(kind=DP), dimension(-nz:nz-1)           ::  h1f, h2f, h3f
  real(kind=DP), dimension(0:nnx,0:2*nz-1)         ::  pb_cef

  real(kind=DP), dimension(-nx:nx) :: vaa, ll

  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ksq, ksqi

  real(kind=DP), dimension(-nx:nx)          :: kx, xi1, kxx, slp_icr
  real(kind=DP), dimension(0:ny)            :: ky
  real(kind=DP), dimension(-nz:nz-1)        :: xi3
  real(kind=DP), dimension(-nz:nz-1)        :: valf
  real(kind=DP), dimension(-nx:nx,-nz:nz-1) :: dpara

  real(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: b0

  complex(kind=DP), dimension(0:ny)       :: ck
  integer, dimension(0:ny)                :: dj

! --- box size
  real(kind=DP) :: lx, ly, lz,llz

! --- basic parameters for reduced MHD
  real(kind=DP) :: nu,   &    ! viscosity
                   eta,  &    ! resitivity
                   s_hat      ! shear parameter

! --- parameters for M-I coupling
  real(kind=DP) :: e0,    &    ! convection electric field
                   idns0, &    ! ionospheric number density
                   mp,    &    ! Pedersen mobility normalized by the ExB mobility
                   mh,    &    ! Hall mobility normalized by the ExB mobility
                   dperp, &    ! diffusion coefficient
                   alpha, &    ! normalized recombination rate : alpha * n0
                   theta2,&
                   delta

  real(kind=DP) :: dt


!--------------------------------------
!  Type of calculation
!--------------------------------------

  character(9)  :: calc_type, fd_type, fd_filt
  character(9)  :: bndry_type

!--------------------------------------
!  Parameters for numerical settings
!--------------------------------------

  integer :: inum
  integer :: loop

! --- unit numbers for I/O
  integer, parameter :: inml =  9,  &
                        olog = 10,  &
                        ogrid = 40, &
                        icnt = 20,  &
                        omag = 30,  &
!                        oion = 31,  &
                        oiodn = 31,  &
                        oiocr = 32,  &
                        oioog = 33,  &
                        ocnt = 50,  &
                        omph = 60,  &
                        omps = 61,  &
                        oidn = 62,  &
                        oicr = 63,  &
                        oiog = 64


END MODULE RMHDFS_header
