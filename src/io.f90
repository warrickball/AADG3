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

module io

  use types, only: dp

  implicit none

contains

  subroutine load_modes(filename, n_modes, l, n, freq, width, power, cycshift)
    character(*), intent(in) :: filename
    integer, intent(out) :: n_modes
    integer, dimension(:), intent(out) :: l, n
    real(dp), dimension(:), intent(out) :: freq, width, power, cycshift

    integer :: i, ierr, iounit

    ierr = 0
    open(newunit=iounit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'ERROR in io.load_modes while opening ', filename
       stop
    end if

    i = 1
    do
       call skip_comments(iounit, ierr)
       if (ierr /= 0) then
          ierr = 0
          exit
       end if
       
       read(iounit, *, iostat=ierr) &
         l(i), n(i), freq(i), width(i), power(i), cycshift(i)
       if (ierr /= 0) exit
       i = i + 1
    end do
    n_modes = i - 1

    if (ierr > 0) then
       write(*,*) 'ERROR in io.load_modes while reading ', filename
       stop
    end if
    
    close(iounit)

  end subroutine load_modes
  

  subroutine load_rotation(filename, low, splitting)
    character(*), intent(in) :: filename
    integer, intent(in) :: low
    real(dp), intent(out) :: splitting(low:,:,:)
    integer :: n, l, m, ierr, iounit
    real(dp) :: s_nlm

    ierr = 0
    open(newunit=iounit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
       write(*,*) 'ERROR in io.load_rotation while opening ', filename
       stop
    end if
    
    do
       call skip_comments(iounit, ierr)
       if (ierr /= 0) then
          ierr = 0
          exit
       end if
       
       read(iounit, *, iostat=ierr) n, l, m, s_nlm
       if (ierr /= 0) then
          write(*,*) 'WARNING in io.load_rotation: unexpected early exit while reading ', filename
          exit
       end if

       splitting(n,l,m) = s_nlm
    end do
    
    if (ierr > 0) then
       write(*,*) 'ERROR in io.load_rotation while reading ', filename
       stop
    end if
    
    close(iounit)

  end subroutine load_rotation


  subroutine skip_comments(iounit, ierr)
    ! Given a file on unit ``iounit``, advances to the next line not
    ! starting with ``!`` or ``#``.
    integer, intent(in) :: iounit
    integer, intent(out) :: ierr

    character(132) :: comment_test

    ierr = 0

    do
       read(iounit, '(a)', iostat=ierr) comment_test
       if (ierr /= 0) then
          return
       elseif ((comment_test(1:1) == '!') .or. &
           (comment_test(1:1) == '#')) then
          cycle
       else
          backspace(iounit)
          return
       end if
    end do
    
  end subroutine skip_comments

end module io
