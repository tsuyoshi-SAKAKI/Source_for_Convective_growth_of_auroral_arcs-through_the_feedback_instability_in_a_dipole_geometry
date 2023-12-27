MODULE dsp_header

  implicit none

  integer, parameter :: DP = selected_real_kind(14)
  integer, parameter :: QP = selected_real_kind(28)
  real(kind=DP), parameter  :: pi = 3.141592653589793_DP


  real(kind=DP), parameter  :: zeta_0 = pi /180._DP * ( 90._DP - &
                                        70._DP) ! theta0
!  real(kind=DP), parameter  :: zeta_0 = pi /180._DP * ( 90._DP - &
!                                        69.08_DP) ! [0:800]
!  real(kind=DP), parameter  :: zeta_0 = pi /180._DP *( 90._DP -  &
!                                        68.385_DP) ! [400:1200]
!  real(kind=DP), parameter  :: zeta_0 = pi /180._DP * ( 90._DP -  &
!                                        67.761_DP) ! [800:1600]

  real(kind=DP)             :: theta0 = zeta_0*180._DP/pi 

  real(kind=DP), parameter :: ell = (1._DP/(2._DP*sqrt(3._DP)*(sin(zeta_0))**2))*  &
                             (sqrt(3._DP)*cos(zeta_0)*sqrt(3._DP*(cos(zeta_0))**2+1._DP)&
                             +log(sqrt(3._DP)*cos(zeta_0)+sqrt(3._DP*(cos(zeta_0))**2+1._DP)))

 
  integer :: olog = 6

  public

END MODULE dsp_header
