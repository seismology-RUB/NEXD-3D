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

module dmatricesMod
  	use constantsMod
	use vandermondeMod
  implicit none

contains

  subroutine dmatrices3d(dr,ds,dt,r,s,t,v)
    ! INPUT: node coordinates reference tet r,s,t, van-der-monde-matrix v
    ! DOES: calculates the differentiation matrices on the simplex at (r,s,t)
    ! OUTPUT: differentiation matrices dr,ds,dt with respect to r,s,t
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
    real(kind=CUSTOM_REAL), dimension(size(v(:,1)),Np) :: dr,ds,dt,w,vr,vs,vt
    
    call gradVdm3D(vr,vs,vt,r,s,t)
    call invVdm3D(v,w,0)

    dr=matmul(vr,w)
    ds=matmul(vs,w)
    dt=matmul(vt,w)
    
  end subroutine dmatrices3d
end module dmatricesMod
