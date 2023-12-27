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
  use RMHDFS_feedback, only: feedback_linana, scale_cal
  use RMHDFS_diag,     only: diag_cntrl


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
    character(512) :: f_log, f_hst, f_fld, f_cnt, f_test, f_grid

    character(4)   :: crank
    character(3)   :: cold, cnew

    character(10)   :: cdate, ctime

    namelist /cmemo/ memo
    namelist /calct/ calc_type, fd_type, fd_filt
    namelist /run_n/ inum
    namelist /files/ f_log, f_hst, f_grid, f_fld, f_cnt


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
      open( ogrid, file=trim(f_grid)//"g"//crank//".grid."//cnew )
      

      if ( inum > 1 ) then
        open( icnt, file=trim(f_cnt)//"g"//crank//".cnt."//cold, &
              form="unformatted", status="old", action="read" )
      end if

      open( omag, file=trim(f_fld)//"g"//crank//".mag."//cnew, form="unformatted" )

      if( rankz == 0 ) then
!         open( oion, file=trim(f_fld)//"g"//crank//".ion."//cnew )
         open( oiodn, file=trim(f_fld)//"g"//crank//".idn."//cnew )
         open( oiocr, file=trim(f_fld)//"g"//crank//".icr."//cnew )
         open( oioog, file=trim(f_fld)//"g"//crank//".iog."//cnew )
      end if

      open( ocnt, file=trim(f_cnt)//"g"//crank//".cnt."//cnew, form="unformatted" )

      if( rank == 0 ) then
        open( omph, file=trim(f_hst)//"mph."//cnew )
        open( omps, file=trim(f_hst)//"mps."//cnew )
        open( oidn, file=trim(f_hst)//"idn."//cnew )
        open( oicr, file=trim(f_hst)//"icr."//cnew )
        open( oiog, file=trim(f_hst)//"iog."//cnew )
      end if


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
     close( ogrid )
     close( icnt )
     close( omag )
!     close( oion )
     close( oiodn )
     close( oiocr )
     close( oioog )
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

    real(kind=DP) :: kxmin, kymin, dx, dz, kxxmin
    real(kind=DP) :: lz_l, z0, z0_l, mu0, x0
    real(kind=DP) :: cn11, cn22, cn33, cn13, co33, bb, x1, x3, radius, the 
    integer       :: n_alp, n_tht, m_j
    real(kind=DP) :: del_c

    integer       :: ix, mx, my, iz, ierr_mpi
    
    real(kind=DP), parameter  ::  bound = pi/36._DP
    real(kind=DP), parameter  ::  theta_box_cntr = theta0
    real(kind=DP), parameter  ::  lnalp = log(alp + sqrt(alp**2 + 1._DP))

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

        x0  = (sin(theta_box_cntr))**2                                                    ! box center of xi1
        mu0 = -1._DP                                                                      ! see lysak u3 

        lx       = -(sin(theta_box_cntr - bound))**2 +(sin(theta_box_cntr + bound))**2    ! xi1-range
        ly       =   2._DP * pi                                                           ! xi2-range
        lz       = -log(alp * mu0 + sqrt(1._DP + (alp * mu0)**2)) / lnalp                 ! xi3-range
        llz = (1._DP/(2._DP*sqrt(3._DP)*(sin(theta0))**2))*                 &
                  (sqrt(3._DP)*cos(theta0)*sqrt(3._DP*(cos(theta0))**2+1._DP)       &
                   +log(sqrt(3._DP)*cos(theta0)+sqrt(3._DP*(cos(theta0))**2+1._DP)))      ! field line length

        write( olog, * ) " # lx            = ", lx
        write( olog, * ) " # ly            = ", ly
        write( olog, * ) " # lz            = ", lz
        write( olog, * ) " # llz           = ", llz
        write( olog, * ) ""


        dx    = lx / 2._DP / real( nx, kind=DP )
        kymin = 598.9260225878468_DP * sin(theta0) / 2._DP
        ly    = 2._DP * pi / kymin


        write( olog, * ) " # lx           = ", lx
        write( olog, * ) " # kymin        = ", kymin
        write( olog, * ) ""


        lz_l     = lz / real( nprocz, kind=DP )                    ! local z-length
        z0       = - lz                                            ! global lower boundary
        z0_l     = lz_l * real( rankz, kind=DP ) + z0              ! local lower boundary
        dz       = lz_l / real( nz, kind=DP ) / 2._DP        


   !-------------------- coordinate calculation ---------------------

       do iz = -nz, nz-1
        do mx = -nx, nx              

          x1 = dx * real( mx, kind=DP ) + x0                  !xi1
          x3 = dz * real( iz + nz, kind=DP ) + z0_l           !xi3
          xi1(mx) = x1
          xi3(iz) = x3

          vaa(mx) = 1._DP 


        call scale_cal(x1,x3,cn11,cn22,cn33,cn13,co33,bb,radius,the)

          ig11(mx,iz)   = cn11
          ig22(mx,iz)   = cn22
          ig33(mx,iz)   = cn33
          ig13(mx,iz)   = cn13
          h1(mx,iz)     = 1._DP / sqrt(cn11)
          h2(mx,iz)     = 1._DP / sqrt(cn22)
          h3(mx,iz)     = sqrt(co33)
          jcb(mx,iz)    = h1(mx,iz) * h2(mx,iz) * h3(mx,iz)    ! Jacobian
          cef_d1(mx,iz) = jcb(mx,iz) * ig11(mx,iz)             ! h2 h3/h1
          cef_d2(mx,iz) = jcb(mx,iz) * ig22(mx,iz)             ! h3 h1/h2
          slp_icr(mx)   = 2._DP * sqrt( (1._DP - xi1(mx) ) / (4._DP - 3._DP * xi1(mx) ) )

          radius_grid(mx,iz) = radius
          theta_grid(mx,iz)  = the

         do my = 0, ny
          b0(mx,my,iz) = bb
         end do

        end do
       end do

        do mx = -nx, nx              
          kxxmin = 2._DP * pi / (  ( asin(sqrt(xi1(nx)))-asin(sqrt(xi1(-nx))) )) 
          kxx(mx)  = kxxmin * real( mx, kind=DP )

          kxmin    = 2._DP * pi / ( xi1(nx) - xi1(-nx)  ) 
          kx(mx)   = kxmin * real( mx, kind=DP )
        end do


        do my = 0, ny
          ky(my)   = kymin * real( my, kind=DP )
        end do


        do iz = -nz, nz-1
          do my = 0, ny
            do mx = -nx, nx
              ksq(mx,my,iz) = ( kx(mx) + s_hat * xi3(iz) * ky(my) )**2 + ky(my)**2

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
            valf (iz) = 1._DP
          end do

          do iz = -nz, nz-1
           do ix = -nx, nx
            dpara(ix,iz) = dz * h3(ix,iz)
           end do
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
    real(kind=DP),    dimension(:,:),     allocatable :: ramp, phase 
    integer :: rri1, rri2

    real(kind=DP)   :: iamp, pamp               ! initial perturbation amplitude
    integer         :: mmi, mmm
    integer, dimension(nx+1) :: mxi, myi        ! mode number of the perturbation
    namelist /initd/ iamp, pamp, mmm, mxi, myi

    integer, parameter :: nzg = global_nz * 2


    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)  :: phi, cpr

    complex(kind=DP), dimension(0:nzg) :: egphi, egpsi
    complex(kind=DP)  :: zom, egdns

    real(kind=DP) :: kpara_r, kpara_i, phase1, dtw

    integer :: input_status
    integer :: ix, iy, mx, my, iz, global_iz

    real(kind=DP)    ::  xx0                    ! center of the initial wave packet

      omg(:,:,:)    = ( 0._DP, 0._DP )
      phi(:,:,:)    = ( 0._DP, 0._DP )
      psi(:,:,:)    = ( 0._DP, 0._DP )
      dns(:,:,:)    = ( 1._DP, 0._DP )

      cpr(:,:,:)    = ( 0._DP, 0._DP )

  
      xx0 = sin(theta0) **2                     ! center of the initial wave packet


      if( rank == 0 ) then
        dns(:,:,-nz)    = ( 0._DP, 0._DP )      ! ionospheric density perturbation
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

        call feedback_linana( 1630.3783_DP , ky(myi(mmi)) / sin(theta0) , zom, egphi, egpsi, egdns, nzg ) 

        write( olog, * ) ""
        write( olog, * ) " # mmi      = ", mmi
        write( olog, * ) " # egphi(0) = ", egphi(0)
        write( olog, * ) " # egpsi(0) = ", egpsi(0)
        write( olog, * ) " # egdns    = ", egdns
        write( olog, * ) ""
        write( olog, * ) " # zom      = ", zom 
        write( olog, * ) " # ky      = ", ky(myi(mmi)) 

        kpara_r = real (zom)
        kpara_i = aimag(zom)

        write( olog, * ) " # kpara_r  = ", kpara_r  
        write( olog, * ) " # kpara_i  = ", kpara_i  


          allocate( rr((2*nx+1)*(ny+1)) )


           rr(:) = 0._DP


      allocate( ramp(-nx:nx,0:ny) )
      allocate( phase(-nx:nx,0:ny) )
      do my = 0, ny
       call math_random ( ramp(-nx:nx,my) )
      end do


      do my = 0, ny
          do ix = -nx, nx
            phase(ix,my) = 2._DP*pi*ramp(ix,my)
               do iz = -nz, nz-1
                  phi(ix,my,iz) = phi(ix,my,iz)                     &
                                + pamp*ksqi(ix,my,iz)/b0(ix,my,iz)  &
                                * exp( ui*phase(ix,my) )
                  psi(ix,my,iz) = psi(ix,iy,iz)                     &
                                + pamp*ksqi(ix,iy,iz)/b0(ix,my,iz)  &
                                * exp( ui*phase(ix,my) ) 
               end do
                if( rank == 0 ) then
                  dns(ix,my,-nz) = dns(ix,my,-nz)                    &
                                 + pamp*ksqi(ix,my,-nz)/b0(ix,my,-nz)&
                                 * exp( ui*phase(ix,my) )
                end if
          end do 
      end do

      


      deallocate( ramp, phase )

                  write(olog,fmt="(a)") "# eigen function "
                  write(olog,fmt="(a)") "# xi1, xi3, phi, psi, dns "
            mx   = mxi(mmi)
            my   = myi(mmi)
              do ix = -nx, nx
                phase1    = 1630.38_DP*asin(sqrt(xi1(ix)))
                do iz = -nz, nz-1
                  global_iz = iz + nz + 2 * nz * rankz


                  phi(ix,my,iz) = phi(ix,my,iz) + iamp / b0(ix,my,iz)                     &
                                * exp( -(xi1(ix)-xx0)**2/( 2._DP * delta**2 ) ) &
                                * egphi(global_iz) * exp( ui*phase1 )
                  psi(ix,my,iz) = psi(ix,my,iz) + iamp / b0(ix,my,iz)                     &
                                * exp( -(xi1(ix)-xx0)**2/( 2._DP * delta**2 ) ) &
                                * egpsi(global_iz) * exp( ui*phase1 )
                end do

                if( rankz == 0 ) then

                  dns(ix,my,-nz)= dns(ix,my,-nz) + iamp * egdns * exp( ui*phase1 ) &
                                * exp( -(xi1(ix)-xx0)**2/( 2._DP * delta**2 ) ) &
                                / b0(ix,my,-nz)

                end if
              end do


          deallocate( rr )

        end do

         write(olog,*) "# feedback linear ana check"
         do iz = -nz,nz-1
           global_iz = iz + nz + 2 * nz * rankz
            write(olog,fmt="(1p,9e15.6)")   xi3(iz),  egphi(global_iz), egpsi(global_iz), phi(0,0,iz), psi(0,0,iz)
         end do
         write(olog,*) 



          call fld_lapl( phi, jcb, cef_d1, cef_d2, omg, 2*nz )
          call fld_lapl( psi, jcb, cef_d1, cef_d2, cpr, 2*nz )
          call fld_calfld( omg, psi, dns, phi, cpr )



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
          do mx = -nx, nx
            omg(mx,my,iz) = cmplx( real(omg(mx,my,iz)), 0._DP )
            phi(mx,my,iz) = cmplx( real(phi(mx,my,iz)), 0._DP )
            psi(mx,my,iz) = cmplx( real(psi(mx,my,iz)), 0._DP )
            cpr(mx,my,iz) = cmplx( real(cpr(mx,my,iz)), 0._DP )
            dns(mx,my,iz) = cmplx( real(dns(mx,my,iz)), 0._DP )
          end do
        end do



        my = 0
        do iz = -nz, nz-1
          do ix = -nx, nx
               write(ogrid,fmt="(1p,8e15.6)") radius_grid(ix,iz),theta_grid(ix,iz), b0(ix,my,iz), h1(ix,iz), h2(ix,iz)&
                                              ,h3(ix,iz), xi1(ix), xi3(iz)
          end do
            write(ogrid,*)
        end do



  END SUBROUTINE set_value


END MODULE RMHDFS_set
