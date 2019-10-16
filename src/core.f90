! Copyright 2018-2019 Warrick Ball & Bill Chaplin

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

module core
  ! The ``core`` module contains the main subroutines that AADG3 uses
  ! to generate the timeseries for an oscillation mode.

  use types, only: dp

  implicit none

contains

  subroutine generate_kicks(n_fine, cadence, tau, nsd, kicks)
    ! Generates an exponentially-damped random walk (i.e. an AR(1)
    ! process) that AADG3 uses to model granulation, for either the
    ! correlated or uncorrelated component of the driving term.
    
    use math, only: randn
    
    integer, intent(in) :: n_fine
    real(dp), intent(in) :: cadence, tau, nsd
    real(dp), intent(out) :: kicks(:)

    real(dp) :: v, vsum, cadfine_tau
    real(dp) :: da(n_fine)
    integer :: i, j, Nsteps

    Nsteps = size(kicks)
    cadfine_tau = cadence/n_fine/tau

    call randn(da)
    
    v = nsd*da(1)
    do i = 1, Nsteps
       call randn(da)

       vsum = v
       do j = 2, n_fine
          v = exp(-cadfine_tau)*v + nsd*da(j)
          vsum = vsum + v
       end do
       kicks(i) = vsum/n_fine
    end do

  end subroutine generate_kicks


  subroutine laplace_solution(n_cadences, n_relax, kick, freq0, damp0, &
       power, dnu, dwid, pcad, phicad, vout)
    ! Implements the solution of the Laplace transformation of the
    ! damped, driven, harmonic oscillator, which AADG3 uses to compute
    ! the timeseries of each oscillation mode.
    
    use math, only: TWOPI
    
    integer, intent(in) :: n_cadences, n_relax
    real(dp), intent(in) :: freq0, damp0, power
    real(dp), intent(in) :: dnu, dwid, pcad, phicad
    real(dp), intent(in) :: kick(n_cadences+n_relax)
    real(dp), intent(out) :: vout(n_cadences)

    real(dp) :: x, v, w1, c1
    real(dp) :: freq, damp, msv
    real(dp) :: asum(n_cadences)
    integer :: i
    real(dp) :: i_cadences(n_cadences)

    i_cadences = [(i, i = 1, n_cadences)]
    asum = sin(TWOPI*(i_cadences-phicad)/pcad)**4

    x = 0
    v = 0

    freq = freq0
    damp = damp0

    do i = 1, n_relax
       w1 = sqrt((TWOPI*freq)**2 - damp**2)
       c1 = kick(i) + damp*x + v
       v = ((c1 - damp*x)*cos(w1) - (w1*x + damp*c1/w1)*sin(w1))*exp(-damp)
       x = (x*cos(w1) + c1/w1*sin(w1))*exp(-damp)
    end do
    
    do i = 1, n_cadences
       w1 = sqrt((TWOPI*freq)**2 - damp**2)
       c1 = kick(n_relax + i) + damp*x + v
       v = ((c1 - damp*x)*cos(w1) - (w1*x + damp*c1/w1)*sin(w1))*exp(-damp)
       freq = freq0 + dnu*asum(i)
       damp = damp0*(1.0_dp + dwid*asum(i))
       x = (x*cos(w1) + c1/w1*sin(w1))*exp(-damp)
       vout(i) = v
    end do
    
    msv = sum(vout**2)/n_cadences
    vout = vout*(1.0_dp - 0.5_dp*dwid*asum)*sqrt(power/msv)
    
  end subroutine laplace_solution
  
end module core
