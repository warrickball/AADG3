! Copyright 2018 Warrick Ball & Bill Chaplin

! This file is part of the AsteroFLAG Artificial Dataset Generator v3 (AADG3).

! AADG3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! AADG3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with AADG3.  If not, see <https://www.gnu.org/licenses/>.

module math
  ! Various mathematical constants and functions.

  use types, only: dp

  implicit none

  real(dp), parameter :: PI = 3.1415926535897932384626433_dp
  real(dp), parameter :: TWOPI = 2*PI
  
contains

  pure integer function factorial(n)
    integer, intent(in) :: n
    integer :: i

    factorial = 1
    do i = 1, n
       factorial = factorial*i
    end do

  end function factorial
  

  subroutine get_Plm(x, lmax, Plm)
    real(dp), intent(in) :: x
    integer, intent(in) :: lmax
    real(dp), intent(out) :: Plm(0:lmax, 0:lmax)

    real(dp) :: y
    integer :: l, m

    if (x == 1d0) then
       Plm = 0d0
       Plm(:,0) = 1d0
    else
       y = sqrt(1.0_dp - x**2)

       Plm = 0
       Plm(0, 0) = 1
       Plm(1, 0) = x
       Plm(1, 1) = -y

       do l = 1, lmax-1
          do m = 0, l
             Plm(l+1, m) = ((2*l+1)*x*Plm(l, m) - (l+m)*Plm(l-1, m))/real(l-m+1, dp)
          end do
          Plm(l+1, l+1) = (x*Plm(l+1, l) - (2*l+1)*Plm(l, l))/y
       end do
    end if
    
  end subroutine get_Plm
  

  subroutine get_Elm(inclination, lmax, Elm)
    real(dp), intent(in) :: inclination   ! in degrees
    integer, intent(in) :: lmax   ! maximum angular degree
    real(dp), intent(out) :: Elm(0:lmax, 0:lmax)

    real(dp) :: x
    real(dp) :: Plm(0:lmax, 0:lmax)
    integer :: l, m, i

    x = cos(inclination*TWOPI/360.0_dp)
    call get_Plm(x, lmax, Plm)
    do l = 0, lmax
       do m = 0, l
          ! this would be a simple version but factorials explode
          ! Elm(l,m) = real(factorial(l-m), dp)/real(factorial(l+m), dp)*Plm(l,m)**2
          ! note that (l-m)!/(l+m)! = 1/[(l+m) x (l+m-1) x ... x (l-m+1)]
          Elm(l, m) = 1.0_dp
          do i = l-m+1, l+m
             Elm(l,m) = Elm(l,m)/real(i, dp)
          end do
          Elm(l,m) = Elm(l,m)*Plm(l,m)**2
       end do
    end do
    
  end subroutine get_Elm
  
  
  ! this is much slower but also much simpler
  ! https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Polar_form

  ! "vector" version, internalizes loop to avoid overhead of
  ! looping over subroutine.  limited gain (~5%)
  subroutine randn(x)
    ! Populates the array ``x`` with normally-distributed random
    ! variates with zero mean and unit variance using the polar form
    ! of the Box-Müller transform.  The length of ``x`` must be even
    ! otherwise the subroutine will ``stop``.
    real(dp), intent(out) :: x(:)
    real(dp) :: w
    integer :: i, n

    n = size(x)

    ! let's just be VERY defensive for now
    if (mod(n, 2) == 1) then
       write(*,*) 'ERROR in math.randn: output array must have '
       write(*,*) '                     even number of elements'
       stop
    end if
    
    call random_number(x)  ! one initial call saves some time
    x = 2d0*x - 1d0

    do i = 1, n, 2
       w = x(i)*x(i) + x(i+1)*x(i+1)  ! was 2d0 to trigger do while
       do while (w >= 1.0)
          call random_number(x(i:i+1))
          ! x(i:i+1) = 2d0*x(i:i+1) - 1d0
          x(i) = 2d0*x(i) - 1d0
          x(i+1) = 2d0*x(i+1) - 1d0
          w = x(i)*x(i) + x(i+1)*x(i+1)
       end do

       w = sqrt(-2d0*log(w)/w)
       ! x(i:i+1) = x(i:i+1)*w
       x(i) = x(i)*w
       x(i+1) = x(i+1)*w
    end do

  end subroutine randn

  ! generator from Leva, "A Fast Normal Random Number Generator",
  ! Transactions on Mathematical Software, Vol. 18, No. 4, December 1992, 449-453.
  ! http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.544.5806&rep=rep1&type=pdf
  subroutine leva_randn(x)
    real(dp) :: x(:)
    integer :: i, n

    n = size(x)
    do i = 1, n
       x(i) = leva_one_randn()
    end do
  end subroutine leva_randn


  real(dp) function leva_one_randn()
    real(dp) :: u, v, x, y, Q
    real(dp), parameter :: s = 0.449871d0
    real(dp), parameter :: t = -0.386595d0
    real(dp), parameter :: a = 0.19600d0
    real(dp), parameter :: b = 0.25472d0

    do
       call random_number(u)
       call random_number(v)
       v = 1.7156d0*(v - 0.5d0)
       x = u - s
       y = abs(v) - t
       Q = x*x + y*(a*y - b*x)
       if (Q < 0.27597d0) then
          exit
       elseif (Q > 0.27846d0) then
          cycle
       elseif (v*v > -4d0*u*u*log(u)) then
          cycle
       end if
    end do

    leva_one_randn = v/u

  end function leva_one_randn

  subroutine lm_to_k(l, m, k)
    ! A one-to-one relation that maps a pair of integers (l, m), with
    ! ``0 <= m <= l``, into a continuous sequence.  Useful for looping
    ! over sequences of angular degrees and azimuthal orders,
    !
    ! This is the inverse of ``k_to_lm``.
    integer, intent(in) :: l, m
    integer, intent(out) :: k

    k = l*(l+1) + m + 1
  end subroutine lm_to_k

  subroutine k_to_lm(k, l, m)
    ! A one-to-one relation that maps a non-negative integer to pair
    ! of integers (l, m), with ``0 <= m <= l``, assuming zero maps to
    ! (0, 0).  Useful for looping over sequences of angular degrees
    ! and azimuthal orders,
    !
    ! This is the inverse of ``lm_to_k``.
    integer, intent(in) :: k
    integer, intent(out) :: l, m

    l = int(sqrt(real(k-1, dp)))
    m = k - 1 - l*(l+1)
  end subroutine k_to_lm
  
end module math
