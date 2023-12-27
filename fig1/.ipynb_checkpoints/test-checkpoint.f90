program kenkyuu
 implicit none
 real(8), parameter   ::   uc1 =2.d0**(7.d0/3.d0)*3.d0**(-1.d0/3.d0)  , uc2 = 2.d0**(1.d0/3.d0)*3.d0**(2.d0/3.d0)
 real(8), parameter   ::   pi = dacos(-1.d0)
 real(8), parameter   ::   bound = pi / 36.d0, alp = 1000.d0
 real(8), parameter   ::   ifalp = dlog(alp + dsqrt(alp ** 2.d0 + 1.d0)) 
 integer,parameter    ::   nx = 12, nz = 25, ny = 100
 real(8)              ::   z0, z0_l, mu0, x0, lz, lx, dx , dz,dy , y0
 real(8),dimension(-nx:nx , -nz:nz)  :: zeta, usi , dbu, gam, dxdm, capth , radius, theta
 real(8),dimension(-nx:nx ,0:ny ,-nz:nz)  :: ig11, ig22, ig33, g33, h1, h2, h3, dcx, dcy,dcz, b0, jcb
 real(8),dimension(-nx:nx)    ::      xx, xi1, theta0, ll
 real(8),dimension(0:ny)  :: xi2
 real(8),dimension(-nz:nz)  ::      zz, xi3, mu
 integer              ::   mx, iz, my

 y0 = 0.d0
 dy = 2.d0 * pi / real(ny, kind=8)


 x0 = (dsin(pi/9.d0))**2.d0                     ! xi1(=nu) の初期値
 mu0 = -dcos(pi/9.d0) / dsqrt(1.d0 - x0)     ! lysak の u3 
! x0 = (sin(pi/9.d0-bound))**2                     ! xi1(=nu) の初期値
! mu0 = -cos(pi/9.d0-bound) / sqrt(1.d0 - x0)     ! lysak の u3 
 lx       = -(dsin(pi/9.d0 - bound))**2.d0 +(dsin(pi/9.d0 + bound))**2.d0    ! total x-range nu方向
 lz       = -dlog(alp * mu0 + dsqrt(1.d0 + (alp * mu0)**2.d0)) / ifalp  ! total z-range 
! dx       = lx / 2.d0 / nx
 dx       = lx / 2.d0 / real(nx, kind=8)
! lz_l     = lz / nz
 z0       = - lz
 dz       = lz / 2.d0 / real(nz, kind=8)
!  dz      = lz / (2.d0 * nz)

 write(*,*)"#dx= " ,dx


 zz = 0.d0
 xx = 0.d0
 xi1 = 0.d0
 xi3 = 0.d0

!do mx = 0, nx
 do iz = -nz, nz
  do my = 0, ny
  do mx = -nx, nx 

  zz(iz) = dz * real(iz + nz, kind=8 ) + z0 
  xx(mx) = dx * real(mx , kind = 8) + x0
  xi2(my) = dy * real(my, kind=8 )+ y0 
 
!  xx(mx) = xi1(mx)
!  zz(iz) = xi3(iz)
   xi3(iz) = zz(iz)
   xi1(mx) = xx(mx)

  mu(iz) = dsinh(ifalp * xi3(iz)) / alp

  if (mu(iz) == 0.d0) then
  
  zeta(mx,iz) = 0.d0
  gam(mx,iz) = 0.d0
  dbu(mx,iz) = 0.d0
  usi(mx,iz) = (dsin(pi / 2.d0))**2.d0

  else

  zeta(mx,iz) = (mu(iz) ** 2.d0 * (1.d0 - xi1(mx))) / (xi1(mx) ** 4.d0 ) 
  gam(mx,iz) = (9.d0 * zeta(mx,iz) + dsqrt(3.d0) * dsqrt(27.d0 * zeta(mx,iz) ** 2.d0 + &
                256.d0 * zeta(mx,iz) ** 3.d0 )) ** (1.d0 / 3.d0)


  dbu(mx,iz) = - uc1 / gam(mx,iz) + gam(mx,iz) / uc2 / zeta(mx,iz)
  usi(mx,iz) = -0.5d0 * dsqrt(dbu(mx,iz)) + 0.5d0 *                                  &
                         dsqrt(-dbu(mx,iz) + 2.d0 / zeta(mx,iz) / dsqrt(dbu(mx,iz))) 

  end if

  theta(mx,iz) = dasin(dsqrt(usi(mx,iz))) 
  radius(mx,iz) = usi(mx,iz) /xi1(mx)
  theta0(mx) = dacos(dsqrt(1.d0 - xi1(mx)))  

  ll(mx) = (radius(mx,iz)/(2.d0*dsqrt(3.d0)*(dsin(theta0(mx)))**2.d0))*                 &
           (dsqrt(3.d0)*dcos(theta0(mx))*dsqrt(3.d0*(dcos(theta0(mx)))**2.d0+1.d0)       &
           +dlog(dsqrt(3.d0)*dcos(theta0(mx))+dsqrt(3.d0*(dcos(theta0(mx)))**2.d0+1.d0)))


  dxdm(mx,iz) = alp / ifalp / dsqrt((alp * mu(iz))**2.d0 + 1.d0)
  capth(mx,iz) = 1.d0 + 3.d0*(dcos(theta(mx,iz)))**2.d0
  ig11(mx,my,iz) = (1.d0 / radius(mx,iz) ** 4.d0) * usi(mx,iz) * capth(mx,iz)
  ig22(mx,my,iz) = 1.d0 / radius(mx,iz) **2.d0 / usi(mx,iz)
  ig33(mx,my,iz) = dxdm(mx,iz) **2.d0 * (1.d0 / radius(mx,iz)**6.d0 / (1.d0 - xi1(mx))**3.d0) *    &
                   ((dcos(theta(mx,iz)) * (1.d0 + 3.d0*(1.d0 - xi1(mx))))**2.d0 / 4.d0   &
                   + usi(mx,iz) * (1.d0 - (1.d0 /radius(mx,iz)))**2.d0)

  g33(mx,my,iz)  = radius(mx,iz) ** 6.d0 * (1.d0 - xi1(mx)) / dxdm(mx,iz) ** 2.d0 / capth(mx,iz)
 
  jcb(mx,my,iz) = radius(mx,iz) ** 6.d0 * (1.d0 - xi1(mx)) ** 0.5d0 / capth(mx,iz) / dxdm(mx,iz)

  b0(mx,my,iz) = 1.d0 / dsqrt( ig11(mx,my,iz) * ig22(mx,my,iz) )

  h1(mx,my,iz) = 1.d0 / dsqrt(ig11(mx,my,iz))
  h2(mx,my,iz) = 1.d0 / dsqrt(ig22(mx,my,iz))
  h3(mx,my,iz) = dsqrt(g33(mx,my,iz)) 
  
  dcx(mx,my,iz) = radius(mx,iz) * sin(theta(mx,iz)) * cos(xi2(my))
  dcy(mx,my,iz) = radius(mx,iz) * sin(theta(mx,iz)) * sin(xi2(my))
  dcz(mx,my,iz) = radius(mx,iz) * cos(theta(mx,iz))

!  cef_d1(mx,iz) = sqrt(ig22(mx,iz) * g33(mx,iz) / ig11(mx,iz) )
!  cef_d2(mx,iz) = sqrt(g33(mx,iz) * ig11(mx,iz) / ig22(mx,iz) )
!  cef_para(mx,my,iz) = 1.d0 / sqrt(g33(mx,iz)) * dxdm(mx,iz)

  open(10, file='test.dat')
  write(10,fmt="(1p,12e15.7)") &
                   xx(mx),xi2(my), zz(iz), h1(mx,my,iz), h2(mx,my,iz), h3(mx,my,iz),dcx(mx,my,iz), &
                  dcy(mx,my,iz), dcz(mx,my,iz), dxdm(mx,iz), radius(mx,iz), theta(mx,iz)
  end do
  write(10,*)
  end do
  write(10,*)
 end do
close(10)

end program kenkyuu
