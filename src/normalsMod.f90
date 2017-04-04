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
module normalsMod
! module to compute the outwart pointing normals at element faces and the surface jacobians
	use constantsMod
 implicit none

contains

  subroutine normals3d(nx,ny,nz,sj,Dr,Ds,Dt,x,y,z,fmask)
    ! INPUT: differentation matrices Dr,Ds,Dt, global coordinates x,y,z, surface mask fmask
    ! DOES: computes normals and their absolute value
    ! OUTPUT: normals nx,ny,nz and absolute value sj
     real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Dr,Ds,Dt
     real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y,z
     integer, dimension(:,:), intent(in) :: fmask
     real(kind=CUSTOM_REAL), dimension(4*NpF), intent(out) :: nx,ny,nz,sj
     real(kind=CUSTOM_REAL), dimension(Np) :: xr,yr,zr,xs,ys,zs,xt,yt,zt,J
     real(kind=CUSTOM_REAL), dimension(Np) :: rx,sx,tx,ry,sy,ty,rz,sz,tz
     real(kind=CUSTOM_REAL), dimension(4*NpF) :: frx,fry,frz,fsx,fsy,fsz,ftx,fty,ftz

     integer, dimension(4*NpF) :: fmaskv
     integer, dimension(NpF) :: fid1,fid2,fid3,fid4
     integer :: i,k,l
!
     k=1
     do l=1,4       ! for all surfaces of tet
        do i=1,NpF  ! for all points on surface
           fmaskv(k)=fmask(i,l) ! make vector with points on surface and 0 if specific point is not on surface
           if (l==1) then
              fid1(i)=k     ! vector with 4 blocks, first number 1 to NpF, rest zeros
           else if (l==2) then
              fid2(i)=k     ! vector with 4 blocks, first zeros, than numbers NpF+1 to 2*NpF, than zeros again
           else if (l==3) then
              fid3(i)=k     ! like above (in third block)
           else if (l==4) then
              fid4(i)=k     ! like above (in fourth block)
           end if
           k=k+1
        end do
     end do
           
    xr = matmul(Dr,x)
    xs = matmul(Ds,x)
    xt = matmul(Dt,x)
    yr = matmul(Dr,y)
    ys = matmul(Ds,y)
    yt = matmul(Dt,y)
    zr = matmul(Dr,z)
    zs = matmul(Ds,z)
    zt = matmul(Dt,z)

    J = xr*(ys*zt-zs*yt) - yr*(xs*zt-zs*xt) + zr*(xs*yt-ys*xt)

    do i=1,size(J)
       if ( abs(J(i) - 0.0) < EPS) then
          write(*,*)  i, J(i)
          stop "jacobian zero in normals"
       end if
    end do

    rx =  (ys*zt - zs*yt)/J
    ry = -(xs*zt - zs*xt)/J
    rz =  (xs*yt - ys*xt)/J

    sx = -(yr*zt - zr*yt)/J
    sy =  (xr*zt - zr*xt)/J
    sz = -(xr*yt - yr*xt)/J

    tx =  (yr*zs - zr*ys)/J
    ty = -(xr*zs - zr*xs)/J
    tz =  (xr*ys - yr*xs)/J


! interpolate geometric factors to face nodes
     frx = rx(fmaskv)
     fsx = sx(fmaskv)
     ftx = tx(fmaskv)
     fry = ry(fmaskv)
     fsy = sy(fmaskv)
     fty = ty(fmaskv)
     frz = rz(fmaskv)
     fsz = sz(fmaskv)
     ftz = tz(fmaskv)

! build normals (use derivatives at surface points as normals)

! face 1
     nx(fid1) = -ftx(fid1)
     ny(fid1) = -fty(fid1)
     nz(fid1) = -ftz(fid1)
! face 2
     nx(fid2) = -fsx(fid2)
     ny(fid2) = -fsy(fid2)
     nz(fid2) = -fsz(fid2)
! face 3
     nx(fid3) = frx(fid3)+fsx(fid3)+ftx(fid3)
     ny(fid3) = fry(fid3)+fsy(fid3)+fty(fid3)
     nz(fid3) = frz(fid3)+fsz(fid3)+ftz(fid3)

! face 4
     nx(fid4) = -frx(fid4)
     ny(fid4) = -fry(fid4)
     nz(fid4) = -frz(fid4)

! normalize
     sj = sqrt(nx*nx + ny*ny + nz*nz)
     nx = nx/sJ
     ny = ny/sJ
     nz = nz/sJ

     sJ = sJ * J(fmaskv(:))
  end subroutine normals3d
end module normalsMod
