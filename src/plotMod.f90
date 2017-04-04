!--------------------------------------------------------------------------
!   Copyright 2012-2016 Lasse Lambrecht (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2015-2017 Andre Lamert (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of NEXD 3D.
!
!   NEXD 3D is free software: you can redistribute it and/or modify it
!   under the terms of the GNU General Public License as published by the
!   Free Software Foundation, either version 3 of the License, or (at your
!   option) any later version.
!
!   NEXD 3D is distributed in the hope that it will be useful, but WITHOUT
!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!   FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License v3.0
!   along with NEXD 3D. If not, see <http://www.gnu.org/licenses/>.
!--------------------------------------------------------------------------
module plotMod
	use constantsMod
  implicit none
! module to save some files to disk for plotting purpose

contains

  subroutine plotPoints2d(x,z,filename)
    implicit none
    real(kind=CUSTOM_REAL), dimension(:) :: x,z
    character(len=80) :: filename
    integer :: i

    open(unit=19,file=trim(filename))
    do i=1,size(x)
       write(19,*) x(i),z(i)
    end do
    close(19)
  end subroutine plotPoints2d

end module plotMod
