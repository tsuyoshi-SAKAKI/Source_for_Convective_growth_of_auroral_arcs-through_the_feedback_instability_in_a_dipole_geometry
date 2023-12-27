MODULE RMHDFS_feedback
!---------------------------------------------------------------------
!     analysis of feedback instability
!           based on the correct treatment of the diffusion term
!
!       Extended for the dipole coordinates
!
!                            by tomo  Jul. 23, 2009.
!---------------------------------------------------------------------

  use RMHDFS_header

  implicit none

  private

!------------
! old parameters
!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet
!------------

  real(kind=DP) :: zl, zet0, dzet

  public feedback_linana, scale_cal


CONTAINS
 
!---------------------------
  SUBROUTINE feedback_linana( kx, ky, zom, egphi, egpsi, egdns, nzg )
!---------------------------

    integer          ::   nzg

    complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
    complex(kind=DP) ::  zom, egdns
    real(kind=DP)    :: kx, ky

    complex(kind=DP) :: zom0

      zl   = llz
      zet0 = 0._DP
      dzet = (1._DP - zet0) / real( nzg, kind=DP )     
!      zom0 = ( 2.d-3, 1.d-4 ) * 1000.d0/zl
      zom0 = ( 2.2d-1, 2.5d-2 )

      call solv( zom0, kx, ky, nzg, zom, egphi, egpsi, egdns )


    return

  END SUBROUTINE feedback_linana


!----------------------------------------------
!     solver in Newton method
!
!----------------------------------------------
  SUBROUTINE solv( zom0, kx, ky, nzg, zom, egphi, egpsi, egdns )

      integer          :: nzg
      real(kind=DP)    :: kx, ky
      complex(kind=DP) :: zom, zom0, egphi(0:nzg), egpsi(0:nzg), egdns

      complex(kind=DP) :: f, df, z, znew

      real(kind=DP)    :: eps, eps0

      integer      isw, i


      data     eps0     /  1.d-8 /

        z          = zom0
        znew       = zom0


        i          = 0


   10   continue
          i          = i + 1
          call  fun ( z,  f, 0, kx, ky, nzg, egphi, egpsi, egdns )
          call dfunr( z, df, 0, kx, ky, nzg, egphi, egpsi, egdns )

!          call dfuni( z, df, 0, 
!     .              kx, ky, nzg, egphi, egpsi, egdns )
!                               ----  these two routines give the
!                                     same result
!                                     (Cauchy-Riemann condition) 

          znew       = z - f / df
          eps        = abs( znew - z ) / abs( z )
          z          = znew
          if( i .ge. 30 ) goto 999
        if( eps .gt. eps0 )   goto 10

  999   continue


        zom = z

        call  fun( z,  f, 1, kx, ky, nzg, egphi, egpsi, egdns )

!
        write(6,*) '   *   solution is obtained ; ', z
        write(6,*) '   *   convergence check ; ', i, eps
!

      return

  END SUBROUTINE solv


!----------------------------------------------
!     function to be solved
!
!----------------------------------------------
  SUBROUTINE fun( z, f, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, f
      integer          :: iflg

      integer          :: nzg
      real(kind=DP)    :: kx, ky
      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns


!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet

      complex(kind=DP) :: zres
      real(kind=DP)    :: omgr, omgi, theta, kperp
      integer          :: nout, k


        omgr       =  dble(z)
        omgi       =  imag(z)
        kperp      =  sqrt( kx**2 + ky**2 )
        theta      =  atan2( ky, kx )

!
! magnetospheric response Phi/J in the slab geom.
!
!        zres       = - ( cdexp(ui*kpara*zl) + cdexp(-ui*kpara*zl) )
!     .               / ( cdexp(ui*kpara*zl) - cdexp(-ui*kpara*zl) )
!


        call response( z, theta, kperp, nzg, zres, egphi, egpsi )

        egdns = kperp**2 * egpsi(0) / ( ui * z - 2.d+0 * alpha )

        f   = ( z + ui * 2.d+0 * alpha ) * ( 1.d+0 + zres * mp * idns0 )      &
            - kperp * e0 * ( mp * cos(theta - theta2) + sin(theta - theta2) ) &
            + ui * kperp**2 * dperp



        if( iflg .eq. 1 ) then
           nout    = 6

           write(nout,630) '#------  '
           write(nout,630) '# kx, ky, omgr, omgi, f '
           write(nout,600) '#   ', kx, ky, z, f
           write(nout,630) '#------  '

  600        format(  a, 1p, 6e15.7 )
  610        format(  a, 1p, 2e15.7 )
  620        format( i5, 1p, 4e15.7 )
  630        format( a )
        end if


        return

  END SUBROUTINE fun


!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfunr( z, df, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, df
      real(kind=DP)    :: kx, ky
      integer          :: iflg, nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns


      complex(kind=DP) :: ff1, ff2, z1, z2
      real(kind=DP)    :: dx, d

      data     dx  /  1.d-10 /


        d          = dble(z) * dx
        z1         = z + dcmplx( d*0.5d0, 0.d0 )
        z2         = z - dcmplx( d*0.5d0, 0.d0 )


        call fun( z1, ff1, iflg, kx, ky, nzg, egphi, egpsi, egdns )
        call fun( z2, ff2, iflg, kx, ky, nzg, egphi, egpsi, egdns )


        df         = ( ff1 - ff2 ) / d


      return


  END SUBROUTINE dfunr


!----------------------------------------------
!     finite diffence function of f
!
!----------------------------------------------
  SUBROUTINE dfuni( z, df, iflg, kx, ky, nzg, egphi, egpsi, egdns )

      complex(kind=DP) :: z, df
      real(kind=DP)    :: kx, ky
      integer          :: iflg, nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
      complex(kind=DP) :: egdns

      complex(kind=DP) :: ff1, ff2, z1, z2
      real(kind=DP)    :: dx, d

      data     dx  /  1.d-7 /


        d          = imag(z) * dx
        z1         = z + dcmplx( 0.d0, d*0.5d0 )
        z2         = z - dcmplx( 0.d0, d*0.5d0 )


        call fun( z1, ff1, iflg, kx, ky, nzg, egphi, egpsi, egdns )
        call fun( z2, ff2, iflg, kx, ky, nzg, egphi, egpsi, egdns )


        df         = ( ff1 - ff2 ) / dcmplx( 0.d0, d )


      return


  END SUBROUTINE dfuni


!----------------------------------------------
!     magnetospheric response Phi/J
!
!----------------------------------------------
  SUBROUTINE response( z, theta, kperp, nzg, zres, egphi, egpsi )

      complex(kind=DP) :: z, zres
      real(kind=DP)    :: theta, kperp
      integer          :: nzg

      complex(kind=DP), dimension(0:nzg) :: egphi, egpsi

!      real*8    zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .          dmeta, dmnu, zet0, dzet
!      common /comprm/ zl, e0, idns0, mp, dperp, alpha, ogi2n,
!     .                dmeta, dmnu, zet0, dzet

      complex(kind=DP) :: psi, chi, phi, omg, bph
      real(kind=DP)    :: kperp_sq, zet, hph, hnu, dzds, sss, b00, valf
      real(kind=DP)    :: hph0, hnu0

      complex(kind=DP), dimension(4) :: chi_k, chi_w, bph_k, bph_w
      real(kind=DP),    dimension(4) :: ck(4), ci(4)

      real(kind=DP)    :: fc1, fc2
      integer          :: i, j, k, isw


!
! coefficients used in the RK method
        ck(1)   = 1.d0/6.d0
        ck(2)   = 1.d0/3.d0
        ck(3)   = 1.d0/3.d0
        ck(4)   = 1.d0/6.d0

        ci(1)   = 0.d0
        ci(2)   = 0.5d0
        ci(3)   = 0.5d0
        ci(4)   = 1.d0


        chi_k = ( 0.d0, 0.d0 )
        chi_w = ( 0.d0, 0.d0 )
        bph_k = ( 0.d0, 0.d0 )
        bph_w = ( 0.d0, 0.d0 )

! boundary condition on the equator
        chi   = ( 0.d0, 0.d0 )
        bph   = ( 1.d0, 0.d0 )


! THW 2022.01.25
        k     = 0
        zet   = zet0 + dzet*dble(k-nzg)
        call metric( zl, zet0, zet, hnu0, hph0, &
                     dzds, sss, b00, valf )

        k     = nzg
        zet   = zet0 + dzet*dble(k-nzg)
        call metric( zl, zet0, zet, hnu, hph, &
                     dzds, sss, b00, valf )

        kperp_sq = kperp**2                    &
                 * ( ( sin(theta) * hph0 / hph )**2   &
                   + ( cos(theta) * hnu0 / hnu )**2 )
!                 * ( ( cos(theta) / hph )**2   &
!                   + ( sin(theta) / hnu )**2 )
! THW 2022.01.25

        psi   = - chi / kperp_sq
        omg   = - kperp_sq * bph / b00

        egpsi(nzg) = psi * b00
        egphi(nzg) = bph


        do k = nzg, 1, -1

! Runge-Kutta starts

          do j = 1, 4
            if( j .eq. 1 ) then
              chi_w(j) = chi
              bph_w(j) = bph
            else
              chi_w(j) = chi + ci(j) * chi_k(j-1) * dzet
              bph_w(j) = bph + ci(j) * bph_k(j-1) * dzet
            end if

            if( j .eq. 1 ) then
! THW 2022.01.25
!              zet   = zet0 + dzet*dble(k)
              zet   = zet0 + dzet*dble(k-nzg)
! THW 2022.01.25
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b00, valf )
            else if ( j .eq. 4 ) then
! THW 2022.01.25
!              zet   = zet0 + dzet*dble(k)
              zet   = zet0 + dzet*dble(k-nzg-1)
! THW 2022.01.25
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b00, valf )
            else
! THW 2022.01.25
!              zet   = zet0 + dzet*dble(k)
              zet   = zet0 + dzet*(dble(k-nzg)-0.5d+0)
! THW 2022.01.25
              call metric( zl, zet0, zet, hnu, hph, &
                           dzds, sss, b00, valf )
            end if

! THW 2022.01.25
            kperp_sq = kperp**2                    &
                     * ( ( sin(theta) * hph0 / hph )**2   &
                       + ( cos(theta) * hnu0 / hnu )**2 )
!                 * ( ( cos(theta) / hph )**2   &
!                   + ( sin(theta) / hnu )**2 )
! THW 2022.01.25

            fc1  = - kperp_sq /( valf**2 * dzds * b00 )
            fc2  = - b00 / ( dzds * kperp_sq )

            chi_k(j) = ui * z * bph_w(j) * fc1
            bph_k(j) = ui * z * chi_w(j) * fc2
          end do

          do j = 1, 4
            chi   = chi + dzet * ck(j)*chi_k(j)
            bph   = bph + dzet * ck(j)*bph_k(j)
          end do

          psi  = - chi / kperp_sq
          omg  = - kperp_sq * bph / b00

          egpsi(k-1) = psi * b00
          egphi(k-1) = bph

        end do

        phi   = - omg / kperp**2

        zres = kperp**2 * phi / chi / &
               (2._DP * sqrt( (1._DP - sin(theta0)**2 ) &
             / (4._DP - 3._DP * sin(theta0)**2 )))


      return


  END SUBROUTINE response


!----------------------------------------------
!     metric
!
!----------------------------------------------

  SUBROUTINE metric( zl, zet0, zet, hnu, hph, dzds, sss, b00, valf )

      real(kind=DP) :: zl, zet0, zet
      real(kind=DP) :: hnu, hph, dzds, sss, b00, valf

      real(kind=DP) :: hnu0, hph0
      real(kind=DP) :: cn11,cn22,cn33,cn13,co33,bb,radius,the
      real(kind=DP), parameter  :: xi0 = (sin(theta0))**2


      call scale_cal(xi0,zet,cn11,cn22,cn33,cn13,co33,bb,radius,the)

        hnu   = 1._DP / sqrt(cn11)
        hph   = 1._DP / sqrt(cn22)
        dzds  = 1._DP / sqrt(co33)
        sss   = zet / dzds 
        b00   = bb
        valf  = 1.d0

      return


  END SUBROUTINE metric


!!! SUBROUTINE metric_slab( zl, zet0, zet, hnu, hph, dzds, sss, b00, valf )
!
!      real(kind=DP) :: zl, zet0, zet
!      real(kind=DP) :: hnu, hph, dzds, sss, b00, valf
!
!      real(kind=DP) :: hnu0, hph0, w10, w20, w1, w2, al
!
!
!        hph   = 1.d0
!        hnu   = 1.d0
!        dzds  = ( 0.5d0 * pi - zet0 ) / zl
!        sss   = zl * ( zet - zet0 ) / ( 0.5d0 * pi - zet0 )
!        b00    = 1.d0
!        valf  = 1.d0
!
!
!      return
!
!!!  END SUBROUTINE metric_slab


!----------------------------------------------
!     metric
!
!----------------------------------------------
  SUBROUTINE metric_dipole( zl, zet0, zet, hnu, hph, dzds, sss, b00, valf )

      real(kind=DP) :: zl, zet0, zet
      real(kind=DP) :: hnu, hph, dzds, sss, b00, valf

      real(kind=DP) :: hnu0, hph0, w10, w20, w1, w2, al


        w10   = sqrt( 3.d0 ) * cos( zet0 )
        w20   = sqrt( 1.d0 + w10**2 )
        hph0  = sin( zet0 )**3
        hnu0  = sin( zet0 )**3 / w20

        al    = zl * 2.d0 * sqrt( 3.d0 )  &
              / ( w10 * w20 + log( abs( w10 + w20 ) ) )

        w1    = sqrt( 3.d0 ) * cos( zet )
        w2    = sqrt( 1.d0 + w1**2 )

        hph   = sin( zet )**3      / hph0
        hnu   = sin( zet )**3 / w2 / hnu0
        dzds  = 1.d0 / ( al * sin( zet ) * w2 )
        sss   = al / ( 2.d0*sqrt(3.d0) )                     &
                   * ( w10 * w20 + log( abs( w10 + w20 ) )   &
                     - w1  * w2  - log( abs( w1  + w2  ) ) )
        b00    = 1.d0 / ( hnu * hph )
        valf  = 1.d0


      return


  END SUBROUTINE metric_dipole


 SUBROUTINE scale_cal(x1,x3,cn11,cn22,cn33,cn13,co33,bb,radius,the)

      real(kind=DP) , parameter  :: uc1 =2._DP**(7._DP/3._DP)*3**(-1._DP/3._DP)
      real(kind=DP) , parameter  :: uc2 = 2**(1._DP/3._DP)*3**(2._DP/3._DP)
      real(kind=DP) , parameter  :: lnalp = log(alp + sqrt(alp**2 + 1._DP))
      real(kind=DP)              :: zeta, dbu, gam
      real(kind=DP), intent(in)  :: x1, x3
      real(kind=DP), intent(out) :: cn11, cn22, cn33, cn13,co33, bb, radius,the
      real(kind=DP)              :: dxdm, capth, the0, usi
      real(kind=DP)              :: mu

           mu = sinh(lnalp * x3) / alp
           
           if (mu == 0._DP) then
            zeta = 0._DP
            gam  = 0._DP
            dbu  = 0._DP
            usi  = ( sin(pi/2._DP) ) ** 2

            else

            zeta = ( mu**2 * (1._DP - x1) ) / x1**4
            gam  = (9._DP * zeta + sqrt(3._DP) * sqrt(27._DP * zeta**2 + 256._DP * zeta**3 ))**(1._DP/3._DP)
            dbu  = - uc1 / gam + gam / uc2 / zeta
            usi  = -0.5_DP * sqrt(dbu) + 0.5_DP * sqrt(-dbu + 2._DP / zeta / sqrt(dbu))

           end if

            the    = asin(sqrt(usi))
            radius = usi /x1
            the0   = acos(sqrt(1._DP - x1)) 
            dxdm   = alp / lnalp / sqrt((alp * mu)**2 + 1._DP)
            capth  = 1._DP + 3._DP * (cos(the))**2

            cn11   = (1._DP / radius**4) * usi * capth 
            cn22   = 1._DP / radius**2 / usi
            cn33   = dxdm**2 * (1._DP / radius**6 / (1._DP - x1)**3) *            &
                          ((cos(the) * (1._DP + 3._DP*(1._DP - x1)))**2 / 4._DP   &
                             + usi * (1._DP - (1._DP /radius))**2)
            cn13   = -dxdm / radius**5 * x1 * cos(the) / 2._DP / (cos(the0))**3 * capth
            co33   = radius**6 * (1._DP - x1) / dxdm**2 / capth
            bb     = sqrt( 1._DP + 3._DP * (cos(the))**2 ) / radius**3            &
                                 / sqrt(1._DP + 3._DP * (cos(theta0))**2 ) 

 END SUBROUTINE scale_cal


END MODULE RMHDFS_feedback
