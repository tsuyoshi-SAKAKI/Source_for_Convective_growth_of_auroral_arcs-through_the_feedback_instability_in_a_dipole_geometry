MODULE RMHDFS_advnc
!-------------------------------------------------------------------------------
!
!    Time integration of the reduced MHD equations
!
!-------------------------------------------------------------------------------

  use RMHDFS_header
  use RMHDFS_mpienv
  use RMHDFS_bndry, only: bndry_bound
  use RMHDFS_fld,   only: fld_calfld, fld_lapl
  use RMHDFS_pssn,  only: pssn_brackets

  implicit none

  private

  public   advnc_rkgsteps


CONTAINS


!--------------------------------------
  SUBROUTINE advnc_rkgsteps( omg, psi, dns )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg, psi, dns


    complex(kind=DP), dimension(:,:,:), allocatable :: domg, qomg
    complex(kind=DP), dimension(:,:,:), allocatable :: dpsi, qpsi
    complex(kind=DP), dimension(:,:,:), allocatable :: ddns, qdns

    integer :: istep, mx, my, iz

      allocate( domg(-nx:nx,0:ny,-nz:nz-1), qomg(-nx:nx,0:ny,-nz:nz-1) ) 
      allocate( dpsi(-nx:nx,0:ny,-nz:nz-1), qpsi(-nx:nx,0:ny,-nz:nz-1) ) 
      allocate( ddns(-nx:nx,0:ny,-nz:nz-1), qdns(-nx:nx,0:ny,-nz:nz-1) ) 


      domg(:,:,:) = ( 0._DP, 0._DP )
      qomg(:,:,:) = ( 0._DP, 0._DP )
      dpsi(:,:,:) = ( 0._DP, 0._DP )
      qpsi(:,:,:) = ( 0._DP, 0._DP )
      ddns(:,:,:) = ( 0._DP, 0._DP )
      qdns(:,:,:) = ( 0._DP, 0._DP )


      do istep = 1, 4


        if( trim(fd_type) == "upwind"  ) then
          call caldlt_u( omg, psi, dns, domg, dpsi, ddns )
        else
          print *, "# selection of fd_type is invalid"
          stop
        end if


        call rkg   ( omg, psi, dns, domg, dpsi, ddns,       &
                                    qomg, qpsi, qdns, istep )


      end do 

      deallocate( domg, qomg ) 
      deallocate( dpsi, qpsi ) 
      deallocate( ddns, qdns ) 


  END SUBROUTINE advnc_rkgsteps


!--------------------------------------
  SUBROUTINE rkg( omg, psi, dns, domg, dpsi, ddns, qomg, qpsi, qdns, istep )
!--------------------------------------
!     Runge-Kutta-Gill

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg,  psi,  dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: domg, dpsi, ddns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: qomg, qpsi, qdns

    integer, intent(in) :: istep

    real(kind=DP) :: c1, c2, cq, c0
    integer :: mx, my, iz


      if      ( istep == 1 ) then
        c1   =  0.5_DP
        c2   = -1._DP
        cq   = -2._DP
        c0   =  1._DP
      else if ( istep == 2 ) then
        c1   =  1._DP - sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 3 ) then
        c1   =  1._DP + sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 4 ) then
        c1   =  1._DP / 6._DP
        c2   = -1._DP / 3._DP
        cq   =  0._DP
        c0   =  0._DP
      end if


!$OMP parallel do collapse(3) private(mx,my,iz)
      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            omg (mx,my,iz) = omg(mx,my,iz) &
                           + c1 * domg(mx,my,iz) + c2 * qomg(mx,my,iz)
            psi (mx,my,iz) = psi(mx,my,iz) &
                           + c1 * dpsi(mx,my,iz) + c2 * qpsi(mx,my,iz)
            dns (mx,my,iz) = dns(mx,my,iz) &
                           + c1 * ddns(mx,my,iz) + c2 * qdns(mx,my,iz)

            qomg(mx,my,iz) = cq * qomg(mx,my,iz) + c0 * domg(mx,my,iz)
            qpsi(mx,my,iz) = cq * qpsi(mx,my,iz) + c0 * dpsi(mx,my,iz)
            qdns(mx,my,iz) = cq * qdns(mx,my,iz) + c0 * ddns(mx,my,iz)
          end do
        end do
      end do


  END SUBROUTINE rkg


!--------------------------------------
  SUBROUTINE caldlt_u( omg, psi, dns, domg, dpsi, ddns )
!--------------------------------------
!     increment within a time step
!     with up-wind for the parallel derivatives
!     Note: this is valid only for the uniform V_A case

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: omg,  psi,  dns
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: domg, dpsi, ddns

    complex(kind=DP), dimension(:,:,:), allocatable :: phi, cpr, lpo
    complex(kind=DP), dimension(:,:), allocatable :: bcpr          !SAKAKIdebug
    complex(kind=DP), dimension(:,:,:), allocatable :: psnb1, psnb2, psnb3
    complex(kind=DP), dimension(:,:,:), allocatable :: dphidz, dcprdz, dpsidz

    complex(kind=DP), dimension(:,:,:), allocatable :: upg, dwg
    complex(kind=DP), dimension(:,:,:), allocatable :: dudz, dddz

!===============SAKAKI=======
    complex(kind=DP), dimension(:,:,:), allocatable :: upgv, dwgv
    complex(kind=DP), dimension(:,:,:), allocatable :: dudzv, dddzv
!===============SAKAKI=======

    real(kind=DP), dimension(-nz:nz-1) :: domg_sum, dpsi_sum, ddns_sum

    real(kind=DP), dimension(-nx:nx) :: a

    integer  ::  ix, my, iz


      allocate( phi(-nx:nx,0:ny,-nz:nz-1) )
      allocate( cpr(-nx:nx,0:ny,-nz:nz-1) )
      allocate( bcpr(-nx:nx,0:ny) ) 
      allocate( lpo(-nx:nx,0:ny,-nz:nz-1) )

      allocate( psnb1(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb2(-nx:nx,0:ny,-nz:nz-1) )
      allocate( psnb3(-nx:nx,0:ny,-nz:nz-1) )

      allocate( dphidz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dcprdz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dpsidz(-nx:nx,0:ny,-nz:nz-1) )

      allocate( upg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dwg(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dudz(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dddz(-nx:nx,0:ny,-nz:nz-1) )

!=============SAKAKI=======
      allocate( upgv(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dwgv(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dudzv(-nx:nx,0:ny,-nz:nz-1) )
      allocate( dddzv(-nx:nx,0:ny,-nz:nz-1) )
!=============SAKAKI=======

      call fld_calfld( omg, psi, dns, phi, cpr )

     
      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( phi, omg, psnb1, 2*nz )
        call pssn_brackets( psi, cpr, psnb2, 2*nz )
        call pssn_brackets( psi, phi, psnb3, 2*nz )

        if( rankz == 0 ) then
!$OMP parallel workshare
          psnb3(:,:,-nz) = ( 0._DP, 0._DP )
!$OMP end parallel workshare
        end if

      else

!$OMP parallel workshare
        psnb1(:,:,:) = ( 0._DP, 0._DP )
        psnb2(:,:,:) = ( 0._DP, 0._DP )
        psnb3(:,:,:) = ( 0._DP, 0._DP )
!$OMP end parallel workshare

      end if


!--- up-wind
!$OMP parallel workshare
        upg(:,:,:)   = ( phi(:,:,:) - psi(:,:,:) ) * b0(:,:,:)
        dwg(:,:,:)   = ( phi(:,:,:) + psi(:,:,:) ) * b0(:,:,:)

        upgv(:,:,:) = ( omg(:,:,:) - cpr(:,:,:) )
        dwgv(:,:,:) = ( omg(:,:,:) + cpr(:,:,:) )

!$OMP end parallel workshare

      call dudz_dddz( upg, dwg, dudz, dddz )

      call dudz_dddz( upgv, dwgv, dudzv, dddzv )

!$OMP parallel workshare
        dphidz(:,:,:) =    0.5_DP * ( dudz(:,:,:) + dddz(:,:,:) ) / b0(:,:,:)
        dcprdz(:,:,:) =  - 0.5_DP * ( dudzv(:,:,:) - dddzv(:,:,:) )
!$OMP end parallel workshare
!--- up-wind


      call fld_lapl( omg, jcb, cef_d1, cef_d2, lpo, 2*nz )


!$OMP parallel do collapse(3) private(ix,my,iz)
      do iz = -nz, nz-1
        do my = 0, ny
          do ix = -nx, nx
            domg(ix,my,iz) = (     -psnb1(ix,my,iz)                                                            &
                   + valf(iz)**2  *( psnb2(ix,my,iz) + vaa(ix)* dcprdz(ix,my,iz) )                             &
                   + nu  * lpo(ix,my,iz) ) * dt
            dpsi(ix,my,iz) = (                                                                                 & 
                   +               ( psnb3(ix,my,iz) + vaa(ix)* dphidz(ix,my,iz) )                             &
                   + eta * cpr(ix,my,iz) ) * dt
            ddns(ix,my,iz) = ( 0._DP, 0._DP )
          end do
        end do
      end do

! ionospheric boundary
      if( rankz == 0 ) then

        do my = 0, ny
          do ix = -nx, nx
            bcpr(ix,my) = cpr(ix,my,-nz) * b0(ix,my,-nz)
          end do
        end do

!        call iono_dns( dns(:,:,-nz), cpr(:,:,-nz), phi(:,:,-nz), ddns(:,:,-nz) )
        call iono_dns( dns(:,:,-nz), bcpr(:,:), phi(:,:,-nz), ddns(:,:,-nz) )

      end if


! --- reality condition
      my = 0
!$OMP parallel do collapse(2) private(ix,iz)
        do iz = -nz, nz-1
          do ix = -nx, nx
            domg(ix,my,iz) = dcmplx( real(domg(ix,my,iz)), 0._DP )
            dpsi(ix,my,iz) = dcmplx( real(dpsi(ix,my,iz)), 0._DP )
            ddns(ix,my,iz) = dcmplx( real(ddns(ix,my,iz)), 0._DP )
          end do
        end do
         

! --- filter
!      if ( trim(fd_filt) == "on" ) then
!        call zfilter ( domg,  1 )
!        call zfilter ( dpsi, -1 )
!      end if


!  zero-zero
      if( bndry_type == "periodic" ) then

!$OMP parallel workshare
        domg_sum(:) = 0._DP
        dpsi_sum(:) = 0._DP
        ddns_sum(:) = 0._DP
!$OMP end parallel workshare

!$OMP parallel do reduction(+:domg_sum,dpsi_sum,ddns_sum)
        do ix = -nx, nx-1
          domg_sum(:) = domg_sum(:) + real( domg(ix,0,:) )
          dpsi_sum(:) = dpsi_sum(:) + real( dpsi(ix,0,:) )
          ddns_sum(:) = ddns_sum(:) + real( ddns(ix,0,:) )
        end do

!$OMP parallel do private(ix)
        do ix = -nx, nx
          domg(ix,0,:) = domg(ix,0,:) - domg_sum(:) / real( 2*nx, kind=DP )
          dpsi(ix,0,:) = dpsi(ix,0,:) - dpsi_sum(:) / real( 2*nx, kind=DP )
          ddns(ix,0,:) = ddns(ix,0,:) - ddns_sum(:) / real( 2*nx, kind=DP )
        end do

      end if


      deallocate( phi, cpr, lpo )
      deallocate( bcpr )      ! SAKAKIdebug 
      deallocate( psnb1, psnb2, psnb3 )
      deallocate( dphidz, dcprdz, dpsidz )

      deallocate( upg, dwg )
      deallocate( dudz, dddz )


  END SUBROUTINE caldlt_u


!--------------------------------------
  SUBROUTINE dudz_dddz( upg, dwg, dudz, dddz )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: upg, dwg

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)  :: dudz, dddz

    integer  ::  mx, my, iz

! --- local variables

    complex(kind=DP), dimension(:,:,:), allocatable :: u_l, d_l

    real(kind=DP),dimension(-nx:nx,-nz:nz-1) :: cefz


!$OMP parallel do private(iz)
      do iz = -nz, nz-1
       do mx = -nx, nx
        cefz(mx,iz)   = 1._DP / ( 6._DP * dpara(mx,iz) )
       end do
      end do


      allocate( u_l(-nx:nx,0:ny,-nz-nb:nz-1+nb) )
      allocate( d_l(-nx:nx,0:ny,-nz-nb:nz-1+nb) )


!$OMP parallel do collapse(3) private(iz,my,mx)
      do iz = -nz, nz-1
        do my = 0, ny
          do mx = -nx, nx
            u_l(mx,my,iz) = upg(mx,my,iz)
            d_l(mx,my,iz) = dwg(mx,my,iz)
          end do
        end do
      end do


      call bndry_bound ( u_l )
      call bndry_bound ( d_l )


      if( rankz == nprocz-1 ) then

!$OMP parallel do collapse(2) private(my,mx)
          do my = 0, ny
            do mx = -nx, nx
              u_l(mx,my,nz  ) = ( 4._DP*u_l(mx,my,nz-1) - u_l(mx,my,nz-2) ) / 3._DP
              u_l(mx,my,nz+1) =   d_l(mx,my,nz-1)

              d_l(mx,my,nz  ) =   u_l(mx,my,nz  )
              d_l(mx,my,nz+1) =   u_l(mx,my,nz-1)
            end do
          end do

      end if


! --- gradient
!$OMP parallel do collapse(2) private(my,mx) 
      do mx = -nx, nx
        do my = 0, ny
          do iz = -nz, nz-1
            dudz(mx,my,iz) = cefz(mx,iz) * (   2._DP * u_l(mx,my,iz+1) &
                                          + 3._DP * u_l(mx,my,iz  ) &
                                          - 6._DP * u_l(mx,my,iz-1) &
                                          +         u_l(mx,my,iz-2) )

            dddz(mx,my,iz) = cefz(mx,iz) * ( - 2._DP * d_l(mx,my,iz-1) &
                                          - 3._DP * d_l(mx,my,iz  ) &
                                          + 6._DP * d_l(mx,my,iz+1) &
                                          -         d_l(mx,my,iz+2) )
          end do
        end do
      end do


      if( rankz == 0 ) then

!$OMP parallel do collapse(2) private(iz,my,mx)

        do mx = -nx, nx
          do my = 0, ny
            iz = -nz+1
            dudz(mx,my,iz) = ( u_l(mx,my,iz) - u_l(mx,my,iz-1) ) /  dpara(mx,iz) 

            iz = -nz
            dddz(mx,my,iz) = (     - d_l(mx,my,iz+2) &
                           + 4._DP * d_l(mx,my,iz+1) &
                           - 3._DP * d_l(mx,my,iz  ) &
                           ) / ( dpara(mx,iz) * 2._DP )
            dudz(mx,my,iz) = dudz(mx,my,iz+1)
          end do
        end do

      end if


  END SUBROUTINE dudz_dddz


!--------------------------------------
  SUBROUTINE zfilter ( vv, isw )
!--------------------------------------
!     z-derivative of f

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: vv

    integer, intent(in)  ::  isw

    complex(kind=DP), dimension(:,:,:), allocatable :: ww

    real(kind=DP) :: aaa

    integer  ::  mx, my, iz, im


      allocate( ww(-nx:nx,0:ny,-nz-nb:nz-1+nb) )


        aaa   = 1._DP


        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              ww(mx,my,iz) = vv(mx,my,iz)
            end do
          end do
        end do


        call bndry_bound ( ww )


        if( rankz == nprocz-1 ) then

          if( isw >= 0 ) then        ! symmetric

!$OMP parallel do collapse(2) private(my,mx)
            do my = 0, ny
              do mx = -nx, nx
                ww(mx,my,nz  ) = ( 4._DP*ww(mx,my,nz-1) - ww(mx,my,nz-2) ) / 3._DP
                ww(mx,my,nz+1) =   ww(mx,my,nz-1)
              end do
            end do

          else if( isw < 0 ) then    ! anti-symmetric

!$OMP parallel do collapse(2) private(my,mx)
            do my = 0, ny
              do mx = -nx, nx
                ww(mx,my,nz  ) = ( 0._DP, 0._DP)
                ww(mx,my,nz+1) = - ww(mx,my,nz-1)
              end do
            end do

          end if

        end if


        if( rankz /= 0 ) then

          do iz = -nz, nz-1
!$OMP parallel do collapse(2) private(my,mx)
            do mx = -nx, nx
              do my = 0, ny
                vv(mx,my,iz) =                                &
                          ( 1._DP - aaa )  * ww(mx,my,iz  )   &
                        + aaa * ( -          ww(mx,my,iz+2)   &
                                  +  4._DP * ww(mx,my,iz+1)   &
                                  + 10._DP * ww(mx,my,iz  )   &
                                  +  4._DP * ww(mx,my,iz-1)   &
                                  -          ww(mx,my,iz-2) ) &
                               / 16._DP
              end do
            end do
          end do

        else

          do iz = -nz+2, nz-1
!$OMP parallel do collapse(2) private(my,mx)
            do my = 0, ny
              do mx = -nx, nx
                vv(mx,my,iz) =                                &
                          ( 1._DP - aaa )  * ww(mx,my,iz  )   &
                        + aaa * ( -          ww(mx,my,iz+2)   &
                                  +  4._DP * ww(mx,my,iz+1)   &
                                  + 10._DP * ww(mx,my,iz  )   &
                                  +  4._DP * ww(mx,my,iz-1)   &
                                  -          ww(mx,my,iz-2) ) &
                               / 16._DP
              end do
            end do
          end do


        end if


       deallocate( ww )


  END SUBROUTINE zfilter


!--------------------------------------
  SUBROUTINE iono_dns( idns, icpr, iphi, didns )
!--------------------------------------
!     increment within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny)  :: idns, icpr, iphi, didns

    complex(kind=DP), &
      dimension(-nx:nx,0:ny)  :: ipsnb

    integer  ::  mx, my


      if( trim(calc_type) == "nonlinear" ) then

        call pssn_brackets( iphi, idns, ipsnb, 1 )

      else

        ipsnb(:,:) = ( 0._DP, 0._DP )

      end if

!$OMP parallel do collapse(2) private(my,mx)

      do my = 0, ny
        do mx = -nx, nx
         didns(mx,my) = ( slp_icr(mx) * icpr(mx,my) - 2.d+0 * alpha * idns(mx,my) &
                      - ipsnb(mx,my) ) * dt
        end do
      end do

      didns(:,0) = dcmplx( real(didns(:,0)), 0._DP )


  END SUBROUTINE iono_dns


END MODULE RMHDFS_advnc
