program kenkyuu
 implicit none
 real(8), parameter   ::   uc1 =2.d0**(7.d0/3.d0)*3**(-1.d0/3.d0)  , uc2 = 2.d0**(1.d0/3.d0)*3.d0**(2.d0/3.d0)
 real(8), parameter   ::   pi = acos(-1.d0)
 real(8), parameter   ::   bound = pi / 36.d0, alp = 1000.d0
 real(8), parameter   ::   ifalp = log(alp + sqrt(alp ** 2 + 1.d0)) 
 integer,parameter    ::   nx = 10, nz = 40000, ny = 1
 real(8)              ::   z0, z0_l, mu0, x0, lz, lx, dx , dz,dy , y0
 real(8),dimension(-nx:nx , -nz:nz)  :: zeta, usi , dbu, gam, dxdm, capth , radius, theta
 real(8),dimension(-nx:nx ,0:ny ,-nz:nz)  :: ig11, ig22, ig33, g33, h1, h2, h3, dcx, dcy,dcz, b0, jcb
 real(8),dimension(-nx:nx)    ::      xx, xi1, theta0, ll
 real(8),dimension(0:ny)  :: xi2
 real(8),dimension(-nz:nz)  ::      zz, xi3, mu
 integer              ::   mx, iz, my

 y0 = 0.d0
 dy = 2.d0 * pi / real(ny, kind=8)


 x0  = (sin(pi/9.d0))**2      
 mu0 = -cos(pi/9.d0)
 lx  = -(sin(pi/9.d0 - bound))**2 +(sin(pi/9.d0 + bound))**2  
 lz  = -mu0
 dx  = lx / 2.d0 / real(nx, kind=8)
 z0  = - lz
 dz  = lz / 2.d0 / real(nz, kind=8)

! write(*,*)"#dx= " ,dx

 zz = 0.d0
 xx = 0.d0
 xi1 = 0.d0
 xi3 = 0.d0

 do iz = -nz, nz
  do my = 0, ny
  do mx = -nx, nx 

  zz(iz) = dz * real(iz + nz, kind=8 ) + z0 
  xx(mx) = dx * real(mx , kind = 8) + x0
  xi2(my) = dy * real(my, kind=8 )+ y0 
 
   xi3(iz) = zz(iz)
   xi1(mx) = xx(mx)

  mu(iz) = xi3(iz)

  if (mu(iz) == 0.d0) then
  
  zeta(mx,iz) = 0.d0
  gam(mx,iz) = 0.d0
  dbu(mx,iz) = 0.d0
  usi(mx,iz) = (sin(pi / 2.d0))**2

  else

  zeta(mx,iz) = mu(iz) ** 2 / xi1(mx) ** 4
  gam(mx,iz) = (9.d0 * zeta(mx,iz) + sqrt(3.d0) * sqrt(27.d0 * zeta(mx,iz) ** 2 + &
                256.d0 * zeta(mx,iz) ** 3 )) ** (1.d0 / 3.d0)


  dbu(mx,iz) = - uc1 / gam(mx,iz) + gam(mx,iz) / uc2 / zeta(mx,iz)
  usi(mx,iz) = -0.5d0 * sqrt(dbu(mx,iz)) + 0.5d0 *                                  &
                         sqrt(-dbu(mx,iz) + 2.d0 / zeta(mx,iz) / sqrt(dbu(mx,iz))) 

  end if

  theta(mx,iz) = asin(sqrt(usi(mx,iz))) 
  radius(mx,iz) = usi(mx,iz) /xi1(mx)
  theta0(mx) = acos(sqrt(1.d0 - xi1(mx)))  

  ll(mx) = (radius(mx,iz)/(2.d0*sqrt(3.d0)*(sin(theta0(mx)))**2))*                 &
           (sqrt(3.d0)*cos(theta0(mx))*sqrt(3.d0*(cos(theta0(mx)))**2+1.d0)       &
           +log(sqrt(3.d0)*cos(theta0(mx))+sqrt(3.d0*(cos(theta0(mx)))**2+1.d0)))


  dxdm(mx,iz) = alp / ifalp / sqrt((alp * mu(iz))**2 + 1.d0)
  capth(mx,iz) = 1.d0 + 3.d0*(cos(theta(mx,iz)))**2
  ig11(mx,my,iz) = (1.d0 / radius(mx,iz) ** 4) * usi(mx,iz) * capth(mx,iz)
  ig22(mx,my,iz) = 1.d0 / radius(mx,iz) **2 / usi(mx,iz)
   g33(mx,my,iz)  = radius(mx,iz) ** 6 / (1.d0 + 3.d0 * (cos(theta(mx,iz)) **2) )

  b0(mx,my,iz) = 1.d0 / sqrt( ig11(mx,my,iz) * ig22(mx,my,iz) )

  h1(mx,my,iz) = 1.d0 / sqrt(ig11(mx,my,iz))
  h2(mx,my,iz) = 1.d0 / sqrt(ig22(mx,my,iz))
  h3(mx,my,iz) = sqrt(g33(mx,my,iz)) 
 
  jcb(mx,my,iz) = 1.d0 / (h1(mx,my,iz)*h2(mx,my,iz)*h3(mx,my,iz)  )   
 
  dcx(mx,my,iz) = radius(mx,iz) * sin(theta(mx,iz)) * cos(xi2(my))
  dcy(mx,my,iz) = radius(mx,iz) * sin(theta(mx,iz)) * sin(xi2(my))
  dcz(mx,my,iz) = radius(mx,iz) * cos(theta(mx,iz))


  open(10, file='fig2_sd.dat')
  write(10,fmt="(1p,5e15.7)") &
                   xx(mx),xi2(my), zz(iz), h3(mx,my,iz), theta(mx,iz)
  end do
  write(10,*)
  end do
  write(10,*)
 end do
close(10)

end program kenkyuu
