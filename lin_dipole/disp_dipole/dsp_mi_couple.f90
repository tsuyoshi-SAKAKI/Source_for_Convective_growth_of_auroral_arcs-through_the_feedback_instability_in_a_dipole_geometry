MODULE dsp_mi_couple

  use dsp_header
  use dsp_slv
  use dsp_fbi
  use dsp_mag

  implicit none

  real(kind=DP)            :: theta = 0._DP

  private

  public mi_couple_find, mi_couple_func


CONTAINS


  SUBROUTINE mi_couple_find

    complex(kind=DP) :: omega, k_para

    complex(kind=DP) :: zin, zout
    real(kind=DP)    :: k_perp

    real(kind=DP)    :: dk, dth 
    integer          :: i, j, k, m

    integer, parameter :: iloop=500, jloop=4 
    integer, parameter :: nth = 1  
    complex(kind=DP), dimension (0:jloop,0:nth-1,-1:1) :: zsave
    character         :: length*20, latitude*20, filename*35


      zin    = cmplx( pi/ell, 0.0_DP ) 

      dk     = pi * 100._DP * 92.6724179453029_DP / real( iloop, kind=DP ) 
      dth    = pi * 0.5_DP  / real( nth, kind=DP )
 
      k_perp = dk

      write(length,*) int(ell*100)
      write(latitude,*) int(theta0*100)
      filename ="l"//trim(adjustl(length))//"_lat"//trim(adjustl(latitude))//'.d'
      open(olog, file=filename)

      call fbi_set_k_perp( k_perp )

      call fbi_init
      call mag_init
      call mi_couple_init


      do m = 0, nth-1

!      do k = 1, -1, -2 
      k = 1 

      do j = 0, jloop

        if( m == 0 ) then
          zin = cmplx( pi*(j+1)*k/ell, 0._DP, kind=DP )
        else
          zin = zsave(j,m-1,k)
        end if

        theta = dth * real( m, kind=DP ) + pi*0.5_DP
!        theta = dth * real( m, kind=DP ) + atan2(0.5_DP,-1._DP) ! finite ky

        call fbi_set_theta( theta )
        call mag_set_theta( theta )


        do i = iloop, 2, -1

          k_perp = dk*real(i,kind=DP)
          call fbi_set_k_perp( k_perp )
          call mag_set_k_perp( k_perp )

          call slv_newton( mi_couple_func, zin, zout )

          k_para = zout
          call fbi_func( omega )

          write(olog,fmt="(1p,6e15.7)") k_para, k_perp, theta, omega

          zin = zout

          if( i == iloop ) then
            zsave(j,m,k) = zout
          end if

        end do

          write(olog,*) "#j=" , j
          write(olog,fmt="(1p,5e15.7)")

      end do


      end do

     close(olog) 

  END SUBROUTINE mi_couple_find


  SUBROUTINE mi_couple_init

      write(olog,fmt="(a,1p,e15.7)")
      write(olog,fmt="(a,1p,e15.7)") "# ell   = ", ell

  END SUBROUTINE mi_couple_init


  FUNCTION mi_couple_func( zin )

    complex(kind=DP) :: mi_couple_func
    complex(kind=DP), intent(in) :: zin

    complex(kind=DP) :: k_para, omega_m, omega_i, kpl
    complex(kind=DP) :: phi_jpara

    complex(kind=DP) :: ui = ( 0._DP, 1._DP )

    real(kind=DP)    :: aomg, delt

    complex(kind=DP) :: k_para_sq


      omega_m = zin
      call mag_response( omega_m, phi_jpara )

      call fbi_set_phi_jpara( phi_jpara )

      call fbi_func( omega_i )


      mi_couple_func = omega_i / omega_m - 1._DP
       

  END FUNCTION mi_couple_func


END MODULE dsp_mi_couple
