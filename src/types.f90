! Copyright 2018-2020 Warrick Ball & Bill Chaplin

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

module types
  ! Defines common type definitions for the rest of AADG3.

  implicit none

  integer, parameter :: dp = kind(0.0d0)

  type mode
     ! Single object containing all the parameters for a single mode
     ! (from the simulation's perspective).
     
     ! All times are in cadences
     real(dp) :: freq
     real(dp) :: damp ! = PI*width
     real(dp) :: power
     real(dp) :: freq_shift = 0.0d0
     real(dp) :: damp_shift = 0.0d0
     real(dp) :: power_shift = 0.0d0
  end type mode

end module types
