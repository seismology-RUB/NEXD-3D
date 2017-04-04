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

module geometricFactorsMod
! module the compute the metric elements for local mappings of the element
	use constantsMod
 implicit none

contains

  subroutine geometricFactors3d(rx,sx,tx,ry,sy,ty,rz,sz,tz,J,x,y,z,Dr,Ds,Dt)
    ! INPUT: differentation matrices Dr, Ds, Dt, global coordinates of points x,y,z
    ! DOES: compute the metric elements for local mappings of an element
    ! OUTPUT: geometric factors rx,sx,tx,ry,sy,ty,rz,sz,tz; Jacobian J
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y,z
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds,Dt
    real(kind=CUSTOM_REAL), dimension(size(x)) :: rx,sx,tx,ry,sy,ty,rz,sz,tz,J
    real(kind=CUSTOM_REAL), dimension(size(x)) :: xr,xs,xt,yr,ys,yt,zr,zs,zt
    integer :: i

    xr = matmul(Dr,x)               ! derivative of with respect to r
    xs = matmul(Ds,x)
    xt = matmul(Dt,x)
    yr = matmul(Dr,y)
    ys = matmul(Ds,y)
    yt = matmul(Dt,y)
    zr = matmul(Dr,z)
    zs = matmul(Ds,z)
    zt = matmul(Dt,z)

    J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt)  ! get Jacobian

    do i=1,size(J)
       if ( abs(J(i) - 0.0) < EPS) then
          write(*,*)  i, J(i)
          stop "Jacobian is zero in geometric factors"
       end if
    end do

    rx =  (ys*zt - zs*yt)/J     ! inverse mapping (derivative of r with respect to x)
    ry = -(xs*zt - zs*xt)/J
    rz =  (xs*yt - ys*xt)/J

    sx = -(yr*zt - zr*yt)/J
    sy =  (xr*zt - zr*xt)/J
    sz = -(xr*yt - yr*xt)/J

    tx =  (yr*zs - zr*ys)/J
    ty = -(xr*zs - zr*xs)/J
    tz =  (xr*ys - yr*xs)/J
  end subroutine geometricFactors3d
end module geometricFactorsMod
