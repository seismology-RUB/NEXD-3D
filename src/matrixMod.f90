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
module matrixMod
! some service routines for matrices
	use constantsMod
 implicit none

contains

  subroutine invert(V)
    !invert matrix
    real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: V
    real(kind=CUSTOM_REAL), dimension(2*size(v(:,1))) :: work
    integer, dimension(size(v(:,1))) :: ipvt
    integer :: ierr, iwork
    integer :: N
    N=size(v(:,1))
    iwork=2*N
    !LU trafo
    if (CUSTOM_REAL==SIZE_REAL) then
        call sgetrf(N,N,V,N,ipvt,ierr)
        if (ierr/=0) write(*,*) "Error LU vdm2d ",ierr
        ! ivert Pr
        call sgetri(N,V,N,ipvt,work,iwork,ierr)
        if (ierr/=0) write(*,*) "Error invert vdm2d ",ierr
    else if (CUSTOM_REAL==SIZE_DOUBLE) then
        call dgetrf(N,N,V,N,ipvt,ierr)
        if (ierr/=0) write(*,*) "Error LU vdm2d DP",ierr
        ! ivert Pr
        call dgetri(N,V,N,ipvt,work,iwork,ierr)
        if (ierr/=0) write(*,*) "Error invert vdm2d DP",ierr
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        stop
    endif
  end subroutine invert

  subroutine matrixLU(V)
    !invert matrix
    real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: V
    integer, dimension(size(v(:,1))) :: ipvt
    integer :: ierr
    integer :: N
    N=size(v(:,1))
    !LU trafo
    call sgetrf(N,N,V,N,ipvt,ierr)
    if (ierr/=0) write(*,*) "Error LU vdm2d ",ierr
  end subroutine matrixLU

end module matrixMod
