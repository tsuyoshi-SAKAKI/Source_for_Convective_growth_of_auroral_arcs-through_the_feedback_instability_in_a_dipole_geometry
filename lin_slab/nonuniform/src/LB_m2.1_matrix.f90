MODULE LB_matrix

  implicit none

  integer, parameter :: DP = selected_real_kind(14)

  private

  public :: LB_matrix_slvtridmrr, LB_matrix_slvtridmrc, LB_matrix_slvtridmcc
  public :: LB_matrix_slvcyctridmrr, LB_matrix_slvcyctridmrc
  public :: LB_matrix_slvperitridmrc


CONTAINS

!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvtridmrr( a, b, c, y, n, x )
!-------------------------------------------------------------------------------
!  matrix : real
!  r.h.s. : real

  integer, intent(in) :: n
  real(kind=DP), intent(in),  dimension(n)   :: b, y
  real(kind=DP), intent(in),  dimension(n-1) :: a, c
  real(kind=DP), intent(out), dimension(n)   :: x

  real(kind=DP), dimension(n)   :: r, u
  integer :: i

    r(1) = b(1)
    u(1) = y(1)
    do i = 2, n
      r(i) = b(i) - a(i-1) * c(i-1) / r(i-1)
      u(i) = y(i) - a(i-1) * u(i-1) / r(i-1)
    end do

    x(n) = u(n) / r(n)
    do i = n-1, 1, -1
      x(i) = ( u(i) - c(i) * x(i+1) ) / r(i)
    end do

    return


  END SUBROUTINE LB_matrix_slvtridmrr


!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvtridmrc( a, b, c, y, n, x )
!-------------------------------------------------------------------------------
!  matrix : real
!  r.h.s. : complex

  integer, intent(in) :: n
  complex(kind=DP), intent(in),  dimension(n)   :: y
  real(kind=DP),    intent(in),  dimension(n)   :: b
  real(kind=DP),    intent(in),  dimension(n-1) :: a, c
  complex(kind=DP), intent(out), dimension(n)   :: x

  complex(kind=DP), dimension(n)   :: u
  real(kind=DP),    dimension(n)   :: r
  integer :: i

    r(1) = b(1)
    u(1) = y(1)
    do i = 2, n
      r(i) = b(i) - a(i-1) * c(i-1) / r(i-1)
      u(i) = y(i) - a(i-1) * u(i-1) / r(i-1)
    end do

    x(n) = u(n) / r(n)
    do i = n-1, 1, -1
      x(i) = ( u(i) - c(i) * x(i+1) ) / r(i)
    end do

!!!    print *, "# n = ", n

    return


  END SUBROUTINE LB_matrix_slvtridmrc


!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvtridmcc( a, b, c, y, n, x )
!-------------------------------------------------------------------------------
!  matrix : complex
!  r.h.s. : complex

  integer, intent(in) :: n
  complex(kind=DP), intent(in),  dimension(n)   :: b, y
  complex(kind=DP), intent(in),  dimension(n-1) :: a, c
  complex(kind=DP), intent(out), dimension(n)   :: x

  complex(kind=DP), dimension(n)   :: r, u
  integer :: i

    r(1) = b(1)
    u(1) = y(1)
    do i = 2, n
      r(i) = b(i) - a(i-1) * c(i-1) / r(i-1)
      u(i) = y(i) - a(i-1) * u(i-1) / r(i-1)
    end do

    x(n) = u(n) / r(n)
    do i = n-1, 1, -1
      x(i) = ( u(i) - c(i) * x(i+1) ) / r(i)
    end do

    return


  END SUBROUTINE LB_matrix_slvtridmcc


!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvperitridmrc( a, b, c, y, n, x )
!
!   to solve the tridiagonal matrix equation 
!   with the periodic boundary condition 
!   where solution of Ax = y with 
!            ( b1 c1 ... ... 0  )
!            ( a2 b2 c2  ... 0  )
!        A = ( 0  a3 b3 c3 . 0  )
!            ( ...              )
!            ( 0  ... ... an bn )
!   is used
!
!   Note the size of A is n-1 while the period is n
!   Condition of <y> = <x> = 0 is assumed
!-------------------------------------------------------------------------------
!  matrix : real
!  r.h.s. : complex

  integer, intent(in) :: n
  complex(kind=DP), intent(in),  dimension(n) :: y
  real(kind=DP),    intent(in),  dimension(n) :: a, b, c
  complex(kind=DP), intent(out), dimension(n) :: x

  complex(kind=DP) :: ave
  integer :: i


    call LB_matrix_slvtridmrc( a(2:n-1), b(1:n-1), c(2:n-1), y(1:n-1), n-1, x(1:n-1) )

    x(n) = (0.d0,0.d0)

    ave = (0.d0,0.d0)
    do i = 1, n
      ave = ave + x(i)
    end do
      ave = ave / dble( n )

    do i = 1, n
      x(i) = x(i) - ave
    end do

    return


  END SUBROUTINE LB_matrix_slvperitridmrc


!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvcyctridmrr( a, b, c, y, n, x )
!
!   to solve the cyclic tridiagonal matrix equation of Ax = y where 
!            ( b1 c1 ... ... a1 )
!            ( a2 b2 c2  ... 0  )
!        A = ( 0  a3 b3 c3 . 0  )
!            ( ...              )
!            ( cn ... ... an bn )
!
!   Note the index rule for ai is different from that of the tridiagonal matrix
!-------------------------------------------------------------------------------
!  matrix : real
!  r.h.s. : real

  integer, intent(in) :: n
  real(kind=DP), intent(in),  dimension(n) :: y
  real(kind=DP), intent(in),  dimension(n) :: a, b, c
  real(kind=DP), intent(out), dimension(n) :: x

  real(kind=DP), dimension(n)   :: p, q, s, t, xi, et
  integer :: i

    p(1)  = c(1) / b(1)
    q(1)  = a(1) / b(1)
    xi(1) = y(1) / b(1)
    s(1)  = c(n)
    t(1)  = b(n)
    et(1) = y(n)
    do i = 2, n-1
      p(i)  = c(i)                    / ( b(i) - a(i)*p(i-1) )
      q(i)  =-a(i)*q(i-1)             / ( b(i) - a(i)*p(i-1) )
      xi(i) = ( y(i) - a(i)*xi(i-1) ) / ( b(i) - a(i)*p(i-1) )

      s(i)  =         - s(i-1)*p (i-1)
      t(i)  = t (i-1) - s(i-1)*q (i-1)
      et(i) = et(i-1) - s(i-1)*xi(i-1)
    end do
      s (n-1) = s(n-1) + a(n)

    xi(n)   = ( et(n-1) -  xi(n-1)           *s(n-1) ) &
            / ( t (n-1) - ( p(n-1) + q(n-1) )*s(n-1) )
    x(n)    = xi(n)

    do i = n-1, 1, -1
      x(i) = -p(i)*x(i+1) - q(i)*x(n) + xi(i)
    end do

    return


  END SUBROUTINE LB_matrix_slvcyctridmrr


!-------------------------------------------------------------------------------
  SUBROUTINE LB_matrix_slvcyctridmrc( a, b, c, y, n, x )
!
!   to solve the cyclic tridiagonal matrix equation of Ax = y where 
!            ( b1 c1 ... ... a1 )
!            ( a2 b2 c2  ... 0  )
!        A = ( 0  a3 b3 c3 . 0  )
!            ( ...              )
!            ( cn ... ... an bn )
!
!   Note the index rule for ai is different from that of the tridiagonal matrix
!-------------------------------------------------------------------------------
!  matrix : real
!  r.h.s. : complex

  integer, intent(in) :: n
  complex(kind=DP), intent(in),  dimension(n) :: y
  real(kind=DP),    intent(in),  dimension(n) :: a, b, c
  complex(kind=DP), intent(out), dimension(n) :: x

  complex(kind=DP), dimension(n)   :: xi, et
  real(kind=DP),    dimension(n)   :: p, q, s, t
  integer :: i

    p(1)  = c(1) / b(1)
    q(1)  = a(1) / b(1)
    xi(1) = y(1) / b(1)
    s(1)  = c(n)
    t(1)  = b(n)
    et(1) = y(n)
    do i = 2, n-1
      p(i)  = c(i)                    / ( b(i) - a(i)*p(i-1) )
      q(i)  =-a(i)*q(i-1)             / ( b(i) - a(i)*p(i-1) )
      xi(i) = ( y(i) - a(i)*xi(i-1) ) / ( b(i) - a(i)*p(i-1) )

      s(i)  =         - s(i-1)*p (i-1)
      t(i)  = t (i-1) - s(i-1)*q (i-1)
      et(i) = et(i-1) - s(i-1)*xi(i-1)
    end do
      s (n-1) = s(n-1) + a(n)

    xi(n)   = ( et(n-1) -  xi(n-1)           *s(n-1) ) &
            / ( t (n-1) - ( p(n-1) + q(n-1) )*s(n-1) )
    x(n)    = xi(n)

    do i = n-1, 1, -1
      x(i) = -p(i)*x(i+1) - q(i)*x(n) + xi(i)
    end do

    return


  END SUBROUTINE LB_matrix_slvcyctridmrc

END MODULE LB_matrix
