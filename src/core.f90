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

  use types, only: dp, mode

  implicit none

contains

  subroutine overtones(n_fine, n_relax, n_cadences, rho, tau_cad, nsdi, &
       period_cad, phase_cad, modes, add_granulation, vout)

    integer, intent(in) :: n_fine, n_relax, n_cadences
    real(dp), intent(in) :: rho, tau_cad, nsdi, period_cad, phase_cad
    type(mode), intent(in) :: modes(:)
    logical, intent(in) :: add_granulation
    real(dp), dimension(n_cadences), intent(out) :: vout

    integer :: i
    real(dp), dimension(n_relax+n_cadences) :: ckick, ukick, kickthis
    real(dp), dimension(n_cadences) :: vi

    vi = 0
    vout = 0

    ckick = 0
    call generate_kicks(n_fine, tau_cad, nsdi, ckick)

    ! loop over the overtones for this (l,|m|) value
    do i = 1, size(modes)
       ukick = 0
       call generate_kicks(n_fine, tau_cad, nsdi, ukick)

       ! Get shift factor for frequency:
       ! dnu = modes(i)% freq_shift  ! activity cycle frequency shift
       ! ... and factor for width:
       ! dwid = (sdnui/0.40)*0.25  ! activity cycle width shift

       !   dwid = (sdnui)/0.40)*0.25*(1.0-dabs((fc(i,ic)/(1000.0*nuac))-2.90d0)/1.1d0)
       !   if ((fc(i,ic)/nuac) < 1800.0.or.((fc(i,ic)/nuac) > 4000.0) dwid = 0.0

       kickthis = rho*ckick + ukick*sqrt(1.0_dp - rho**2)

       call laplace_solution(n_relax, n_cadences, kickthis, modes(i), &
            period_cad, phase_cad, vi)
       vout = vout + vi
    end do

    if (add_granulation) then
       vout = vout + ckick(n_relax+1:n_relax+n_cadences)
    end if

    return

  end subroutine overtones


  subroutine generate_kicks(n_fine, tau_cad, nsd, kicks)
    ! Generates an exponentially-damped random walk (i.e. an AR(1)
    ! process) that AADG3 uses to model granulation, for either the
    ! correlated or uncorrelated component of the driving term.
    
    use math, only: randn
    
    integer, intent(in) :: n_fine
    real(dp), intent(in) :: tau_cad, nsd
    real(dp), intent(out) :: kicks(:)

    real(dp) :: v, vsum, tau_fine_inv
    real(dp) :: da(n_fine)
    integer :: i, j

    tau_fine_inv = 1.0_dp/n_fine/tau_cad

    call randn(da)
    
    v = nsd*da(1)
    do i = 1, size(kicks)
       call randn(da)

       vsum = v
       do j = 2, n_fine
          v = exp(-tau_fine_inv)*v + nsd*da(j)
          vsum = vsum + v
       end do
       kicks(i) = vsum/n_fine
    end do

  end subroutine generate_kicks


  subroutine laplace_solution(n_relax, n_cadences, kick, m, &
       period_cad, phase_cad, vout)
    ! Implements the solution of the Laplace transformation of the
    ! damped, driven, harmonic oscillator, which AADG3 uses to compute
    ! the timeseries of each oscillation mode.
    
    use math, only: TWOPI
    
    integer, intent(in) :: n_relax, n_cadences
    type(mode), intent(in) :: m
    real(dp), intent(in) :: period_cad, phase_cad
    real(dp), intent(in) :: kick(n_cadences+n_relax)
    real(dp), intent(out) :: vout(n_cadences)

    real(dp) :: x, v, w1, c1
    real(dp) :: freq, damp, msv
    real(dp) :: asum(n_cadences)
    integer :: i
    real(dp) :: i_cadences(n_cadences)

    i_cadences = [(i, i = 1, n_cadences)]
    asum = sin(TWOPI*(i_cadences-phase_cad)/period_cad)**4

    x = 0
    v = 0

    freq = m% freq
    damp = m% damp

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
       freq = m% freq + m% freq_shift*asum(i)
       damp = m% damp + m% damp_shift*asum(i)
       x = (x*cos(w1) + c1/w1*sin(w1))*exp(-damp)
       vout(i) = v
    end do
    
    msv = sum(vout**2)/n_cadences
    vout = vout*(1.0_dp - 0.5_dp*m% damp_shift/m% damp*asum)*sqrt(m% power/msv)
    
  end subroutine laplace_solution
  
end module core
