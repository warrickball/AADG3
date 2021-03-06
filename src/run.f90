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

program AADG3

  use types, only: dp, mode
  use math, only: PI, TWOPI, k_to_lm
  use core, only: overtones
  
  implicit none

  integer :: ierr, i, l, m
  integer, parameter :: ntype = 16
  integer :: nc(ntype), user_seed
  integer :: n_cadences, n_relax, n_fine
  real(dp), allocatable :: vtotal(:), v(:)
  real(dp) :: freq_factor, tau_cad
  real(dp) :: rho, tau, sig
  real(dp) :: inclination, cycle_period, cycle_phase, sdnu(ntype)
  real(dp) :: period_cad, phase_cad, nsdsum, p(0:3), nuac
  real(dp) :: cadence, nsd(ntype), pvis(ntype), ptot
  type(mode), dimension(2000, ntype) :: modes
  character(50) :: input_filename, modes_filename, rotation_filename, output_filename
  character(50) :: output_fmt
  logical :: add_granulation
  integer :: iounit
  integer :: n_seed
  integer, allocatable :: seed(:)
  logical :: verbose

  namelist /controls/ user_seed, n_relax, n_cadences, n_fine, cadence, &
       sig, rho, tau, inclination, cycle_period, cycle_phase, nuac, sdnu, p, &
       add_granulation, modes_filename, rotation_filename, &
       output_filename, output_fmt, verbose

  include 'defaults.in'
  
  ierr = 0

  call parse_args

  if (verbose) then
     write(*,*) ''
     call show_version
     write(*,*) ''
  end if
  
  if (verbose) write(*,'(A)',advance='no') 'Checking command line arguments... '
  call check_args
  if (verbose) write(*,'(A)') 'Done.'

  freq_factor = 1d-6*cadence

  if (verbose) write(*,'(A)',advance='no') 'Initialising random number generator... '
  if (user_seed == 0) then
     call random_seed
  elseif (user_seed > 0) then
     call random_seed(size=n_seed)
     allocate(seed(n_seed))
     seed(1) = user_seed
     do i = 2, n_seed
        seed(i) = user_seed + mod(i**2 + 14897, 1234)
     end do
     call random_seed(put=seed)
     deallocate(seed)
  else
     if (verbose) write(*,*) ''
     write(*,*) 'ERROR in run.f90: user_seed must be >= 0'
     stop 1
  endif
  
  allocate(vtotal(n_cadences))
  
  call random_number(vtotal)  ! warms up RNG
  vtotal = 0
  if (verbose) write(*,'(A)') 'Done.'
  
  if (verbose) write(*,'(A)',advance='no') 'Loading mode data... '
  call get_modes
  if (verbose) write(*,'(A)') 'Done.'
  
  ! Next, determine the variances for the correlated noise arrays
  ! this is buggy if angular degrees are missing because correlated
  ! noise won't be added if ncomp == 0
  ptot = sum(pvis)

  tau_cad = tau/cadence

  ! Calculate variances:
  nsd = sqrt(1.0_dp/tau_cad/n_fine)*sig*sqrt(pvis/ptot)
  nsdsum = sum(nsd**2)
  
  ! convert stellar cycle period and phase into multiples of cadence
  period_cad = 2.0*cycle_period*86400.0*365.25/cadence
  phase_cad = 2.0*cycle_phase*86400.0*365.25/cadence

  ! main loop, to make overtones of each (l,m):
  if (verbose) write(*,'(A)') 'Computing oscillations... '
  if (verbose) write(*,'(3A6)') 'l', 'm', 'nc'
  
  !$OMP PARALLEL DO SCHEDULE(dynamic) PRIVATE(i,l,m,v) REDUCTION(+:vtotal)
  do i = ntype, 1, -1
     allocate(v(n_cadences))
     call k_to_lm(i, l, m)
     if (verbose) write(*,'(4I6,A4,I2)') l, m, nc(i), ntype+1-i, ' of ', ntype
     call overtones(n_fine, n_relax, n_cadences, rho, &
          tau_cad, nsd(i), period_cad, phase_cad, modes(1:nc(i), i), &
          add_granulation, v)
     vtotal = vtotal + v
     deallocate(v)
  end do
  !$OMP END PARALLEL DO
  
  if (verbose) write(*,'(A)') 'Finished computing oscillations.'

  ! Now output time series to disk:
  if (verbose) write(*,'(A)',advance='no') 'Saving output timeseries to '//trim(output_filename)//'... '
  ierr = 0
  open(newunit=iounit, file=output_filename, status='replace', iostat=ierr)
  if (ierr /= 0) then
     if (verbose) write(*,*) ''
     write(*,*) 'ERROR in AADG3 while opening '//trim(output_filename)
     stop
  end if
  write(iounit, output_fmt) (vtotal(i), i=1, n_cadences)
  close(iounit)
  if (verbose) write(*,'(A)') 'Done.'

  deallocate(vtotal)

contains
  
  subroutine get_modes
    use io, only: load_rotation, skip_comments
    use math, only: get_Elm, lm_to_k
    
    integer :: k, l, m, n
    real(dp) :: x, d1, d2, d3, d4
    real(dp) :: s(-2000:200,3,3)
    real(dp) :: Elm(0:3,0:3)
    integer :: iounit

    x = cos(inclination*TWOPI/360d0)

    ! fix visibilities of the individual components:
    ! Gizon & Solanki (2003), eq. (11)
    ! pvis(l,m) = V2(l)*E(l,m)
    ! k = l(l+1) + m + 1
    
    call get_Elm(inclination, 3, Elm)
    do l = 0, 3
       do m = -l, l
          call lm_to_k(l, m, k)
          pvis(k) = p(l)*Elm(l,abs(m))
       end do
    end do

    if (rotation_filename == '') then
       s = 0d0
    else
       call load_rotation(rotation_filename, lbound(s, 1), s)
    end if

    ierr = 0
    open(newunit=iounit, file=modes_filename, status='old', iostat=ierr)
    if (ierr /= 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in AADG3.get_modes while opening '//trim(modes_filename)
       stop
    end if

    nc = 0
    do
       call skip_comments(iounit, ierr)
       if (ierr /= 0) then
          ierr = 0
          exit
       end if
       
       read(iounit, *, iostat=ierr) l, n, d1, d2, d3, d4
       if (ierr /= 0) then
          if (verbose) write(*,*) ''
          write(*,*) 'WARNING in AADG3.get_modes: unexpected early exit while reading '//trim(modes_filename)
          exit
       end if

       do m = -l, l
          call lm_to_k(l, m, k)
          if ((pvis(k) > 1d-8) .and. (d3 > 1d-8)) then
             nc(k) = nc(k) + 1
             if (m == 0) then
                modes(nc(k), k)% freq = d1
             else
                modes(nc(k), k)% freq = d1 + m*s(n,l,abs(m))
             end if
             modes(nc(k), k)% freq = modes(nc(k), k)% freq*freq_factor
             modes(nc(k), k)% damp = d2*freq_factor*PI
             modes(nc(k), k)% power = pvis(k)*d3
             modes(nc(k), k)% freq_shift = d4*freq_factor*sdnu(k)
             modes(nc(k), k)% damp_shift = &
                  modes(nc(k), k)% damp*sdnu(k)/0.4d0*0.25d0
          end if
          
       end do
    end do

    if (ierr > 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in AADG3.get_modes while reading '//trim(modes_filename)
       stop
    end if
    
    close(iounit)

    return
    
  end subroutine get_modes

  
  subroutine parse_args
    character(80) :: arg
    integer :: i
    
    call get_command_argument(1, arg)
    if (arg == '-h' .or. arg == '--help') then
       call show_help
       stop
    else if (arg == '-V' .or. arg == '--version') then
       call show_version
       stop
    else if (arg == ' ') then
       call show_usage
       stop
    end if
    
    call get_command_argument(1, input_filename)
    
    ! Get basic parameters from the information file:
    open(newunit=iounit, file=input_filename, status='old', iostat=ierr)
    
    if (ierr /= 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in AADG3 while opening '//trim(input_filename)
       stop 1
    end if
    
    read(iounit, nml=controls, iostat=ierr)
    close(iounit)
    
    if (ierr /= 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in AADG3 while reading '//trim(input_filename)
       write(*,*) 'the following runtime error may be informative'
       open(newunit=iounit, file=input_filename, status='old', iostat=ierr)
       read(iounit, nml=controls)
       stop 1
    end if

    ! loop over remaining command line arguments to get override input
    ! values
    i = 2
    call get_command_argument(i, arg)
    do while (arg  /= ' ')
       ! write(*,*) i, arg
       
       ! integer options
       select case (arg)
       case ('--user_seed', '--user-seed')
          call get_int_arg(i, user_seed)
       case ('--n_relax', '--n-relax')
          call get_int_arg(i, n_relax)
       case ('--n_cadences', '--n-cadences')
          call get_int_arg(i, n_cadences)
       case ('--n_fine', '--n-fine')
          call get_int_arg(i, n_fine)
       ! float options
       case ('--cadence')
          call get_real_arg(i, cadence)
       case ('--sig')
          call get_real_arg(i, sig)
       case ('--rho')
          call get_real_arg(i, rho)
       case ('--tau')
          call get_real_arg(i, tau)
       case ('--inclination')
          call get_real_arg(i, inclination)
       case ('--cycle_period', '--cycle-period')
          call get_real_arg(i, cycle_period)
       case ('--cycle_phase', '--cycle-phase')
          call get_real_arg(i, cycle_phase)
       ! logical options
       case ('--add_granulation', '--add-granulation')
          add_granulation = .true.
       case ('--no-add_granulation', '--no-add-granulation')
          add_granulation = .false.
       case ('--verbose', '-v')
          verbose = .true.
       case ('--no-verbose', '--quiet', '-q')
          verbose = .false.
       ! string options
       case ('--modes_filename', '--modes-filename')
          i = i + 1
          call get_command_argument(i, modes_filename)
       case ('--rotation_filename', '--rotation-filename')
          i = i + 1
          call get_command_argument(i, rotation_filename)
       case ('--output_filename', '--output-filename', '-o')
          i = i + 1
          call get_command_argument(i, output_filename)
       case default
          if (verbose) write(*,*) ''
          write(*,*) 'ERROR in AADG3 while parsing command-line arguments'
          write(*,*) 'argument '//trim(arg)//' is not valid'
          stop 1
       end select

       i = i + 1
       call get_command_argument(i, arg)
    end do
    
  end subroutine parse_args


  subroutine get_int_arg(i, var)

    integer, intent(inout) :: i
    integer, intent(out) :: var
    character(80) :: arg
    
    i = i + 1
    call get_command_argument(i, arg)
    read(arg, *, iostat=ierr) var
    if (ierr /= 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in get_int_arg: could not parse string '//trim(arg)//' as integer'
       stop 1
    end if

  end subroutine get_int_arg
  

  subroutine get_real_arg(i, var)

    integer, intent(inout) :: i
    real(dp), intent(out) :: var
    character(80) :: arg

    i = i + 1
    call get_command_argument(i, arg)
    read(arg, *, iostat=ierr) var
    if (ierr /= 0) then
       if (verbose) write(*,*) ''
       write(*,*) 'ERROR in get_real_arg: could not parse string '//trim(arg)//' as real(dp)'
       stop 1
    end if

  end subroutine get_real_arg
  

  subroutine check_args

    call check_bounds_dp('inclination', inclination, lower=0d0, upper=90d0)
    call check_bounds_dp('rho', rho, lower=0d0, upper=1d0)

    call check_bounds_dp('cadence', cadence, lower=0d0)
    call check_bounds_dp('sig', sig, lower=0d0)
    call check_bounds_dp('tau', tau, lower=0d0)
    call check_bounds_dp('cycle_period', cycle_period, lower=0d0)
    call check_bounds_dp('cycle_phase', cycle_phase, lower=0d0)
    
    call check_bounds_int('n_fine', n_fine, lower=0)
    call check_bounds_int('n_relax', n_relax, lower=0)
    call check_bounds_int('n_cadences', n_cadences, lower=0)
    
  end subroutine check_args


  subroutine check_bounds_dp(name, val, lower, upper)
    
    character(*), intent(in) :: name
    real(dp), intent(in) :: val
    real(dp), intent(in), optional :: lower, upper

    if (present(lower)) then
       if (val < lower) then
          if (verbose) write(*,*) ''
          write(*,*) 'ERROR: '//trim(name)//' must be greater than', lower
          write(*,*) '       but got ', val
          stop 1
       end if
    end if

    if (present(upper)) then
       if (val > upper) then
          if (verbose) write(*,*) ''
          write(*,*) 'ERROR: '//trim(name)//' must be less than', upper
          write(*,*) '       but got ', val
          stop 1
       end if
    end if
    
  end subroutine check_bounds_dp


  subroutine check_bounds_int(name, val, lower, upper)
    
    character(*), intent(in) :: name
    integer, intent(in) :: val
    integer, intent(in), optional :: lower, upper

    if (present(lower)) then
       if (val < lower) then
          if (verbose) write(*,*) ''
          write(*,*) 'ERROR: '//trim(name)//' must be greater than', lower
          write(*,*) '       but got ', val
          stop 1
       end if
    end if

    if (present(upper)) then
       if (val > upper) then
          if (verbose) write(*,*) ''
          write(*,*) 'ERROR: '//trim(name)//' must be less than', upper
          write(*,*) '       but got ', val
          stop 1
       end if
    end if
    
  end subroutine check_bounds_int
  
  
  subroutine show_help

    write(*,*)
    write(*,*) "asteroFLAG Artificial Dataset Generator 3 (AADG3)"
    write(*,*) "Copyright 2018 Warrick Ball & Bill Chaplin"
    write(*,*) "https://warrickball.github.io/AADG3"
    write(*,*)
    call show_usage
    write(*,*)
    write(*,*) 'Initial options (must be first argument):'
    write(*,*) '  -h, --help            show this help and exit'
    write(*,*) '  -V, --version         show version and exit'
    write(*,*) ''
    write(*,*) 'You can override any of the options in the input file by adding the option '
    write(*,*) 'preceded by `--` and followed by the desired value.  e.g. the argument '
    write(*,*) '`--n-cadences 1000` would override the parameter `n_cadences` in the input file '
    write(*,*) 'and replace it with 1000.  For logical options, either set the parameter with '
    write(*,*) '`--name` or unset it with --no-name`. e.g. `--no-add-granulation`.'
    write(*,*)

  end subroutine show_help


  subroutine show_usage

    write(*,*) 'usage: AADG3 [-h] [-V] <input file> [options]'
    
  end subroutine show_usage

  subroutine show_version

    write(*,*) 'AADG3 v3.0.2'

  end subroutine show_version
  
end program AADG3
