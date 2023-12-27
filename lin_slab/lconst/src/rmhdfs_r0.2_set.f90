MODULE RMHDFS_set
!-------------------------------------------------------------------------------
!
!    Set file I/O, and read parameters from namelist
!
!-------------------------------------------------------------------------------

  use RMHDFS_header
  use RMHDFS_mpienv
  use RMHDFS_math,     only: math_random
  use RMHDFS_fld,      only: fld_calfld, fld_lapl
  use RMHDFS_bndry,    only: bndry_bound
  use RMHDFS_feedback, only: feedback_linana

  implicit none

  private

  character(8)  :: equib_type

  public   set_init, set_close



CONTAINS


!--------------------------------------
  SUBROUTINE set_init( omg, psi, dns, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    real(kind=DP), intent(out) :: time


      call set_start
      call set_param


      if ( trim(bndry_type) == "periodic" .OR. &
           trim(bndry_type) == "zerofix" ) then

        if ( rank == 0 ) then
          write(olog,*) "# boundary type = ", bndry_type
        end if

      else

        if ( rank == 0 ) then
          write(*,*) "# error!! on namelist: bndry_type"
        end if
        call MPI_Finalize (ierr_mpi)
        stop

      end if


      if ( trim(equib_type) == "slab"   .OR. &
           trim(equib_type) == "dipole" .OR. &
           trim(equib_type) == "torus" ) then

        call set_cnfig

      else

        if ( rank == 0 ) then
          write(*,*) "set_cnfig_error!! on namelist: equib"
        end if
        call MPI_Finalize (ierr_mpi)
        stop

      end if


      call set_value( omg, psi, dns, time )


    return


  END SUBROUTINE set_init


!--------------------------------------
  SUBROUTINE set_start
!--------------------------------------

    character(128) :: memo
    character(512) :: f_log, f_hst, f_fld, f_cnt

    character(4)   :: crank
    character(3)   :: cold, cnew

    character(10)   :: cdate, ctime

    namelist /cmemo/ memo
    namelist /calct/ calc_type, fd_type, fd_filt
    namelist /run_n/ inum
    namelist /files/ f_log, f_hst, f_fld, f_cnt


      call date_and_time( cdate, ctime )


      read(inml,nml=cmemo)


      read(inml,nml=calct)


      inum = 1
      read(inml,nml=run_n)


      read(inml,nml=files)

      write( crank, fmt="(i4.4)" ) rank
      write( cold,  fmt="(i3.3)" ) inum-1
      write( cnew,  fmt="(i3.3)" ) inum


      open( olog, file=trim(f_log)//"g"//crank//".log."//cnew )

      if ( inum > 1 ) then
        open( icnt, file=trim(f_cnt)//"g"//crank//".cnt."//cold, &
              form="unformatted", status="old", action="read" )
      end if

      open( omag, file=trim(f_fld)//"g"//crank//".mag."//cnew, form="unformatted" )

      if( rankz == 0 ) then
         open( oion, file=trim(f_fld)//"g"//crank//".ion."//cnew )
      end if

      open( ocnt, file=trim(f_cnt)//"g"//crank//".cnt."//cnew, form="unformatted" )

      if( rank == 0 ) then
        open( omph, file=trim(f_hst)//"mph."//cnew )
        open( omps, file=trim(f_hst)//"mps."//cnew )
        open( oidn, file=trim(f_hst)//"idn."//cnew )
        open( oicr, file=trim(f_hst)//"icr."//cnew )
        open( oiog, file=trim(f_hst)//"iog."//cnew )
      end if


      write( olog, * ) "##### ", trim(memo), " #####"
      write( olog, * ) ""
      write( olog, * ) "# Date : ", cdate
      write( olog, * ) "# Time : ", ctime
      write( olog, * ) ""
      write( olog, * ) "# Type of calc.  : ", trim(calc_type)
      write( olog, * ) "#   finite diff. : ", trim(fd_type)
      write( olog, * ) "#   filter       : ", trim(fd_filt)
      write( olog, * ) ""
      write( olog, * ) "# Run number = ", inum
      write( olog, * ) ""


    return
  

  END SUBROUTINE set_start


!--------------------------------------
  SUBROUTINE set_close
!--------------------------------------
 
     close( olog )
     close( icnt )
     close( omag )
     close( oion )
     close( ocnt )

     if( rank == 0 ) then
       close( omph )
       close( omps )
     end if


  END SUBROUTINE set_close


!--------------------------------------
  SUBROUTINE set_param
!--------------------------------------

    namelist /runlm/ e_limit
    namelist /times/ tend, dtout_mag, dtout_ion, dtout_eng
    namelist /deltt/ dt
    namelist /equib/ equib_type
    namelist /bndry/ bndry_type


      e_limit   = 1._DP * 3600._DP - 300._DP

      read(inml,nml=runlm)


      tend      = 10.00_DP
      dtout_mag = 0.5_DP
      dtout_ion = 0.5_DP
      dtout_eng = 0.005_DP

      read(inml,nml=times)

      dt        = 1._DP

      read(inml,nml=deltt)
      read(inml,nml=equib)
      read(inml,nml=bndry)

      if( dtout_eng < dt ) dtout_eng = dt


        write( olog, * ) " # Numerical parameters "
        write( olog, * ) ""
        write( olog, * ) " # nxw, nyw  = ", nxw, nyw
        write( olog, * ) " # global_nz = ", global_nz
        write( olog, * ) ""
        write( olog, * ) " # nx, ny, nz   = ", nx, ny, nz
        write( olog, * ) ""
        write( olog, * ) " # e_limit      = ", e_limit
        write( olog, * ) " # tend         = ", tend
        write( olog, * ) " # dtout_mag = ", dtout_mag
        write( olog, * ) " # dtout_ion = ", dtout_ion
        write( olog, * ) " # dtout_eng = ", dtout_eng
        write( olog, * ) ""
        write( olog, * ) " # time step dt = ", dt
        write( olog, * ) ""

        write( olog, * ) " # MPI parallelization parameters "
        write( olog, * ) ""
        write( olog, * ) " # nproc , rank  = ", nproc , rank
        write( olog, * ) ""
        write( olog, * ) " # nprocz, rankz = ", nprocz, rankz
        write( olog, * ) ""
        write( olog, * ) " # izup, izdn    = ", izup, izdn
        write( olog, * ) ""


  END SUBROUTINE set_param


!--------------------------------------
  SUBROUTINE set_cnfig
!--------------------------------------

    real(kind=DP) :: kxmin, kymin, dx, dz
    real(kind=DP) :: lz_l, z0, z0_l

    integer       :: n_alp, n_tht, m_j
    real(kind=DP) :: del_c

    integer       :: ix, mx, my, iz, ierr_mpi


    namelist /physp/  nu,   &   ! viscosity
                      eta,  &   ! resistivity
                      s_hat     ! shear parameter

    namelist /ionop/  e0,    &   ! convection electric field
                      idns0, &   ! ionospheric number density
                      mp,    &   ! Pedersen mobility normalized by the ExB mobility
                      mh,    &   ! Hall mobility normalized by the ExB mobility
                      dperp, &   ! diffusion coefficient
                      alpha, &   ! normalized recombination rate : alpha * n0
                      theta2,&
                      delta 


    namelist /boxsize/ lx, ly, lz



      read(inml,nml=physp)

      read(inml,nml=ionop)


        write( olog, * ) " # Physical parameters"
        write( olog, * ) ""
        write( olog, * ) " # nu           = ", nu
        write( olog, * ) " # eta          = ", eta
        write( olog, * ) " # s_hat        = ", s_hat
        write( olog, * ) ""
        write( olog, * ) " # e0           = ", e0
        write( olog, * ) " # idns0        = ", idns0
        write( olog, * ) " # mp           = ", mp
        write( olog, * ) " # dperp        = ", dperp
        write( olog, * ) " # alpha        = ", alpha
        write( olog, * ) " # theta2       = ", theta2
        write( olog, * ) " # delta        = ", delta
        write( olog, * ) ""


! --- coordinate settings ---

        lx       =   0.10790697190940028_DP        ! total x-length
        ly       =   0.5_DP                        ! total y-length
        lz       =   10.790697190940028_DP * 0.5_DP  ! total z-length

      read(inml,nml=boxsize)

        dx       = lx / real( nx, kind=DP )

        kxmin    = pi / lx
        kymin    = pi / ly


        write( olog, * ) " # lx           = ", lx
        write( olog, * ) " # kymin        = ", kymin
        write( olog, * ) ""


        lz_l     = lz / real( nprocz, kind=DP )       ! local z-length


        z0       = 0._DP                    ! global lower boundary
        z0_l     = lz_l * real( rankz, kind=DP ) + z0
                                            ! local lower boundary

        dz       = lz_l / real( nz, kind=DP )

        do iz = -nz, nz-1
          zz(iz)   = dz * real( iz + nz, kind=DP ) + z0_l
        end do

          radius = (4._DP*lz*sqrt(3._DP)*(sin(pi/9._DP))**2)/                            &
                   (sqrt(3._DP)*cos(pi/9._DP)*sqrt(3._DP*(cos(pi/9._DP))**2+1._DP)       &
                   +log(sqrt(3._DP)*cos(pi/9._DP)+sqrt(3._DP*(cos(pi/9._DP))**2+1._DP))) 

          x_0 = (radius*pi)/9._DP

        do mx = -nx, nx
          xx(mx)   = dx    * real( mx, kind=DP )
          kx(mx)   = kxmin * real( mx, kind=DP )

          theta0(mx) = (xx(mx)+x_0)/radius
  
          ll(mx) = (radius/(2._DP*sqrt(3._DP)*(sin(theta0(mx)))**2))*                        &
                   (sqrt(3._DP)*cos(theta0(mx))*sqrt(3._DP*(cos(theta0(mx)))**2+1._DP)       &
                   +log(sqrt(3._DP)*cos(theta0(mx))+sqrt(3._DP*(cos(theta0(mx)))**2+1._DP))) 
!          ll(mx)  = 2._DP * lz
          
!          vaa(mx) = 2._DP*lz/ll(mx)
          vaa(mx) = 1._DP
        end do

          write(olog,*) "# Field line setting:"
          write(unit=olog,fmt="(1p, 3e15.7)") ( xx(ix), ll(ix), vaa(ix), ix=-nx,nx )
          write(olog,*)


        do my = 0, ny
          ky(my)   = kymin * real( my, kind=DP )
        end do


        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              ksq(mx,my,iz) = ( kx(mx) + s_hat * zz(iz) * ky(my) )**2 + ky(my)**2

              if ( mx == 0  .and.  my == 0 ) then
                ksqi(mx,my,iz) = 0._DP
              else
                ksqi(mx,my,iz) = 1._DP / ksq(mx,my,iz)
              end if

            end do
          end do
        end do


        write( olog, * ) " # Numerical parameters"
        write( olog, * ) ""
        write( olog, * ) " # lx, ly, lz   = ", lx, ly, lz
        write( olog, * ) " # lz,   z0     = ", lz, z0
        write( olog, * ) " # lz_l, z0_l   = ", lz_l, z0_l
        write( olog, * ) " # kxmin, kymin = ", kxmin, kymin
        write( olog, * ) " # kxmax, kymax = ", kxmin*nx, kymin*ny
        write( olog, * ) " # kperp_max    = ", sqrt((kxmin*nx)**2+(kymin*ny)**2)
        write( olog, * ) ""
        write( olog, * ) " # dx           = ", dx
        write( olog, * ) " # dz           = ", dz
        write( olog, * ) ""


! --- coordinate settings ---

 
! --- operator settings ---



!!! for the concentric and large-aspect-ratio model !!!
        if( trim(equib_type) == "slab" ) then

          do iz = -nz, nz-1
            b0(iz)    = 1._DP
            valf (iz) = 1._DP
            dpara(iz) = dz
          end do

        else

          write( olog, * ) " # wrong choice of the equilibrium "
          stop

        end if


! --- operator settings ---



  END SUBROUTINE set_cnfig


!--------------------------------------
  SUBROUTINE set_value( omg, psi, dns, time )
!--------------------------------------

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: omg, psi, dns

    real(kind=DP), intent(out) :: time

    real(kind=DP),    dimension(:),       allocatable :: rr
    integer :: rri1, rri2

    real(kind=DP)   :: iamp, pamp     ! initial perturbation amplitude
    integer         :: mmi, mmm
    integer, dimension(nx+1) :: mxi, myi ! mode number of the perturbation
    namelist /initd/ iamp, pamp, mmm, mxi, myi

    integer, parameter :: nzg = global_nz * 2


    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
    complex(kind=DP)  :: zom, egdns

    real(kind=DP) :: kpara_r, kpara_i, phase1, dtw


    integer :: input_status
    integer :: ix, mx, my, iz, global_iz


      omg(:,:,:)    = ( 0._DP, 0._DP )
      phi(:,:,:)    = ( 0._DP, 0._DP )
      psi(:,:,:)    = ( 0._DP, 0._DP )
      dns(:,:,:)    = ( 1._DP, 0._DP )


      if( rank == 0 ) then
        dns(:,:,-nz)    = ( 0._DP, 0._DP )   ! ionospheric density perturbation
      end if


      if( inum == 1 ) then

        time     = 0._DP

        mxi(:) = 0
        myi(:) = 0

        read(inml,nml=initd)


        write( olog, * ) " # iamp, pamp, = ", iamp, pamp
        write( olog, * ) " # mxi   = ", mxi(1:mmm)
        write( olog, * ) " # myi   = ", myi(1:mmm)
        write( olog, * ) ""


        do mmi = 1, mmm

        call feedback_linana( kx(mxi(mmi)), ky(myi(mmi)), zom, egphi, egpsi, egdns, nzg ) 

        write( olog, * ) ""
        write( olog, * ) " # mmi      = ", mmi
        write( olog, * ) " # egphi(0) = ", egphi(0)
        write( olog, * ) " # egpsi(0) = ", egpsi(0)
        write( olog, * ) " # egdns    = ", egdns
        write( olog, * ) " # zom      = ", zom
        write( olog, * ) ""

        kpara_r = real (zom)
        kpara_i = aimag(zom)


        write(olog,*) "# kx(mxi(mmi)), ky(myi(mmi)), omega =", kx(mxi(mmi)), ky(myi(mmi)), zom


          allocate( rr((2*nx+1)*(ny+1)) )

          call math_random ( rr )




            do my =   0, ny
              do mx = -nx, nx
                do ix = -nx, nx
                  phase1    = 2._DP * pi * rr(nx+1+mx+my*(2*nx+1)) + kx(mx)*xx(ix)
                  do iz = -nz, nz-1
                  global_iz = iz + nz + 2 * nz * rankz
!!! Ver Sep 2019
!                  omg(ix,my,iz) = omg(ix,my,iz)                  &
!                                +                  pamp / b0(iz) &
!                                * egphi(global_iz) * exp( ui*phase1 )
! May 29, 2020

!Before
!                  phi(ix,my,iz) = phi(ix,my,iz)                          &
!                                + ksqi(mx,my,iz) * pamp / b0(iz)         &
!                                * egphi(global_iz) * exp( ui*phase1 )
!                  psi(ix,my,iz) = psi(ix,my,iz)                          &
!                                + ksqi(mx,my,iz) * pamp / b0(iz)         &
!                                * egpsi(global_iz) * exp( ui*phase1 )
!After
                  phi(ix,my,iz) = phi(ix,my,iz)                                     &
                                + ksqi(mx,my,iz) * pamp / b0(iz)                    &
                                * exp( -(xx(ix))**2/( 2._DP * delta**2 ) )          &
                                * egphi(global_iz) * exp( ui*phase1 )
                  psi(ix,my,iz) = psi(ix,my,iz)                                     &
                                + ksqi(mx,my,iz) * pamp / b0(iz)                    &
                                * exp( -(xx(ix))**2/( 2._DP * delta**2 ) )          &
                                * egpsi(global_iz) * exp( ui*phase1 )
                  end do
                end do
              end do
            end do


                  write(olog,fmt="(a)") "# eigen function "
                  write(olog,fmt="(a)") "# zz, phi, psi "
            mx   = mxi(mmi)
            my   = myi(mmi)
              do ix = -nx, nx
                phase1    = 2._DP * pi * rr(nx+1+mx+my*(2*nx+1)) + kx(mx)*xx(ix)
                do iz = -nz, nz-1
                  global_iz = iz + nz + 2 * nz * rankz
!!! Ver Sep 2019
!!!                  omg(ix,my,iz) = omg(ix,my,iz) &
!!!                                - ksq(mx,my,iz) &
!!!                                * iamp / b0(iz) * egphi(global_iz) * exp( ui*phase1 )
! May 29, 2020

!Before
!                  phi(ix,my,iz) = phi(ix,my,iz)                                          &
!                                + iamp / b0(iz) * egphi(global_iz) * exp( ui*phase1 )
!                  psi(ix,my,iz) = psi(ix,my,iz)                                          &
!                                + iamp / b0(iz) * egpsi(global_iz) * exp( ui*phase1 )

!After
                  phi(ix,my,iz) = phi(ix,my,iz) + iamp / b0(iz)                     &
                                * exp( -(xx(ix))**2/( 2._DP * delta**2 ) ) &
                                * egphi(global_iz) * exp( ui*phase1 )
                  psi(ix,my,iz) = psi(ix,my,iz) + iamp / b0(iz)                     &
                                * exp( -(xx(ix))**2/( 2._DP * delta**2 ) ) &
                                * egpsi(global_iz) * exp( ui*phase1 )

!!!                  write(olog,fmt="(1p,5e15.7)") zz(iz), phi(mx,my,iz), psi(mx,my,iz)
                end do
!!!                  write(olog,*)

                if( rankz == 0 ) then
!Before
!                  dns(ix,my,-nz)= dns(ix,my,-nz) + iamp * egdns * exp( ui*phase1 ) 
!After
                  dns(ix,my,-nz) = dns(ix,my,-nz) + iamp                             &
                                 * exp( -(xx(ix))**2/( 2._DP * delta**2 ) ) &
                                 * egdns * exp( ui*phase1 ) 
!!!                  write(olog,fmt="(a,1p,3e15.7)") "# zz, dns ", zz(-nz), dns(mx,my,-nz)
!!!                  write(olog,*) "#"
                end if
              end do


          deallocate( rr )

        end do

! May 29, 2020
          call fld_lapl( phi, omg, 2*nz )

!!!          omg(:,:,:) = - ksq(:,:,:) * phi(:,:,:)

          call fld_calfld( omg, psi, dns, phi, cpr )


!--- imported from f0.45

      else


        time   = - 1._DP

        do 
          read( unit=icnt, iostat=input_status ) time, dtw, omg, psi, dns

          if ( input_status < 0 ) then
            write( olog, * ) &
               " # end of file of unit=30 is detected --> stop"
            call MPI_Finalize ( ierr_mpi )
            stop
          end if

          if ( input_status > 0 ) then
            write( olog, * ) &
               " # input error of unit=30 is detected --> stop"
            call MPI_Finalize ( ierr_mpi )
            stop
          end if

          if ( time > eps ) exit
        end do

          dt = min( dt, dtw )

          write( olog, * ) ""
          write( olog, * ) " # simulation is re-started at t = ", time
          write( olog, * ) " # with dt = ", dt
          write( olog, * ) ""

      end if


! --- reality condition
      my = 0
        do iz = -nz, nz-1
!!!          do mx = 1, nx
!!!            omg(mx,my,iz) = conjg( omg(-mx,-my,iz) )
!!!            psi(mx,my,iz) = conjg( psi(-mx,-my,iz) )
!!!            dns(mx,my,iz) = conjg( dns(-mx,-my,iz) )
!!!          end do
          do mx = -nx, nx
            omg(mx,my,iz) = cmplx( real(omg(mx,my,iz)), 0._DP )
            psi(mx,my,iz) = cmplx( real(psi(mx,my,iz)), 0._DP )
            dns(mx,my,iz) = cmplx( real(dns(mx,my,iz)), 0._DP )
          end do
        end do



  END SUBROUTINE set_value


END MODULE RMHDFS_set
