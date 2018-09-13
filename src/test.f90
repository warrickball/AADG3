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

program test_AADG3
  use types

  implicit none

  call test_factorial
  call test_get_Plm
  call test_get_Elm
  call test_k_to_lm
  call test_lm_to_k
  call test_laplace_solution

contains

  subroutine test_factorial
    use math, only: factorial
    integer, parameter :: n = 3
    integer :: i, input(n), output(n)

    input = [1, 6, 12]
    output = [1, 720, 479001600]

    do i = 1, n
       if (factorial(input(i)) /= output(i)) then
          write(*,'(A24,A8)') 'test_factorial          ', 'FAILED'
          write(*,*) 'i, input(i) =', i, input(i)
          return
       end if
    end do

    write(*,'(A24,A8)') 'test_factorial          ', 'PASSED'
    
  end subroutine test_factorial

  subroutine test_get_Plm
    use math, only: get_Plm
    integer, parameter :: lmax = 3
    integer :: l, m
    real(dp) :: Plm(0:lmax,0:lmax), Plm_check(0:lmax, 0:lmax)
    real(dp) :: x, y

    x = 0.5_dp
    y = sqrt(1.0_dp - x**2)

    Plm_check = 0
    Plm_check(0, 0) = 1
    Plm_check(1, 0) = x
    Plm_check(1, 1) = -y
    Plm_check(2, 0) = (3.0_dp*x**2 - 1.0_dp)/2.0_dp
    Plm_check(2, 1) = -3.0_dp*x*y
    Plm_check(2, 2) = 3.0_dp*y**2
    Plm_check(3, 0) = (5.0_dp*x**3 - 3.0_dp*x)/2.0_dp
    Plm_check(3, 1) = -1.5_dp*(5.0_dp*x**2 - 1.0_dp)*y
    Plm_check(3, 2) = 15.0_dp*x*y**2
    Plm_check(3, 3) = -15.0_dp*y**3

    call get_Plm(x, lmax, Plm)

    do l = 0, lmax
       do m = 0, lmax
          if (.not. almost_equal(Plm(l, m), Plm_check(l, m))) then
             write(*,'(A24,A8)') 'test_get_Plm            ', 'FAILED'
             write(*,*) 'l, m =', l, m
             return
          end if
       end do
    end do

    write(*,'(A24,A8)') 'test_get_Plm            ', 'PASSED'
    
  end subroutine test_get_Plm
  

  subroutine test_get_Elm
    use math, only: get_Elm, PI
    integer, parameter :: lmax = 2
    real(dp) :: Elm(0:lmax,0:lmax), Elm_check(0:lmax,0:lmax)
    real(dp) :: inc_degrees, inc_radians
    integer :: l, m

    inc_degrees = 60
    inc_radians = inc_degrees*PI/180

    Elm_check = 0

    ! from Gizon & Solanki (2003), eq. (12) -- (16)
    Elm_check(0, 0) = 1
    Elm_check(1, 0) = cos(inc_radians)**2
    Elm_check(1, 1) = sin(inc_radians)**2/2.0_dp
    Elm_check(2, 0) = (3.0_dp*cos(inc_radians)**2 - 1.0_dp)**2/4.0_dp
    Elm_check(2, 1) = 0.375_dp*sin(2*inc_radians)**2
    Elm_check(2, 2) = 0.375_dp*sin(inc_radians)**4

    call get_Elm(inc_degrees, lmax, Elm)

    do l = 0, lmax
       do m = 0, l
          if (.not. almost_equal(Elm(l, m), Elm_check(l, m))) then
             write(*,'(A24,A8)') 'test_get_Elm            ', 'FAILED'
             write(*,*) 'l, m =', l, m
             return
          end if
       end do
    end do
    
    write(*,'(A24,A8)') 'test_get_Elm            ', 'PASSED'
    
  end subroutine test_get_Elm


  subroutine test_k_to_lm
    use math, only: k_to_lm
    integer, parameter :: lmax = 20
    integer :: l, m, k, ll, mm

    k = 1
    do l = 0, lmax
       do m = -l, l
          call k_to_lm(k, ll, mm)
          if (l /= ll) then
             write(*,'(A24,A8)') 'test_k_to_lm            ', 'FAILED'
                write(*,*) 'k, l, ll =', k, l, ll
             return
          end if
          k = k+1
       end do
    end do
    
    write(*,'(A24,A8)') 'test_k_to_lm            ', 'PASSED'

  end subroutine test_k_to_lm
  

  subroutine test_lm_to_k
    use math, only: lm_to_k
    integer, parameter :: lmax = 20
    integer :: l, m, k, kk
    k = 1
    do l = 0, lmax
       do m = -l, l
          call lm_to_k(l, m, kk)
          if (k /= kk) then
             write(*,'(A24,A8)') 'test_lm_to_k            ', 'FAILED'
             write(*,*) 'l, m, k, kk =', l, m, k, kk
             return
          end if
          k = k+1
       end do
    end do
    
    write(*,'(A24,A8)') 'test_lm_to_k            ', 'PASSED'
    
  end subroutine test_lm_to_k


  subroutine test_laplace_solution
    use core, only: laplace_solution
    use math, only: TWOPI

    integer :: n_cadences, n_relax
    real(dp), allocatable :: kick(:), v(:), v_test(:), i_cadences(:)
    real(dp) :: freq0, damp0, power, dnu, dwid, pcad, phicad
    integer :: i

    n_relax = 1000
    n_cadences = 10000 + n_relax

    allocate(kick(n_relax + n_cadences))
    allocate(i_cadences(n_cadences))
    allocate(v(n_cadences))
    allocate(v_test(n_cadences))
    i_cadences = [(i, i = 1, n_cadences)]

    kick = 0
    kick(1) = 1
    freq0 = 0.001_dp
    damp0 = 0
    power = 1
    dnu = 0
    dwid = 0
    pcad = 100
    phicad = 0

    call laplace_solution(n_cadences, n_relax, kick, freq0, damp0, &
         power, dnu, dwid, pcad, phicad, v)

    ! https://physics.stackexchange.com/questions/101129/harmonic-oscillator-driven-by-a-dirac-delta-like-force
    v_test = sqrt(2.0_dp)*cos(TWOPI*i_cadences*freq0)

    do i = 1, n_cadences
       if (.not. almost_equal(v(i), v_test(i))) then
          write(*,'(A24,A8)') 'test_laplace_solution   ', 'FAILED'
          write(*,*) i, v(i), v_test(i)
          deallocate(kick, v, v_test)
          return
       end if
    end do

    deallocate(kick, v, v_test)

    write(*,'(A24,A8)') 'test_laplace_solution   ', 'PASSED'
    
  end subroutine test_laplace_solution

  
  pure logical function almost_equal(x, y)
    real(dp), intent(in) :: x, y
    real(dp), parameter :: ATOL = 1e-12_dp

    almost_equal = (abs(x-y) < ATOL)
    
  end function almost_equal

end program test_AADG3
