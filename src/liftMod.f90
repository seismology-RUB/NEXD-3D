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
module liftMod
! module to compute the lifting from surf to vol
	use constantsMod
	use vandermondeMod
 implicit none

contains

  subroutine lift3D(lift,fmask,vdm,r,s,t)
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(Npf,Npf) :: v,m,inv
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: vdm
    integer, dimension(:,:), intent(in) :: fmask
    real(kind=CUSTOM_REAL), dimension(Np,4*NpF), intent(out) :: lift
    integer :: iface,i
    real(kind=CUSTOM_REAL), dimension(Np,4*NpF) :: emat,temp
    real(kind=CUSTOM_REAL), dimension(NpF) :: faceR, faceS
    integer, dimension(Npf) :: idr,idc,id
    emat(:,:) = 0

    do iface=1,4
       v=0
       inv=0
       m=0
       if (iface == 1) then
          id=fmask(:,1)
          faceR=r(id)
          faceS=s(id)
       else if (iface == 2) then
          id=fmask(:,2)
          faceR=r(id)
          faceS=t(id)
       else if (iface == 3) then
          id=fmask(:,3)
          faceR=s(id)
          faceS=t(id)
       else if (iface == 4) then
          id=fmask(:,4)
          faceR=s(id)
          faceS=t(id)
       end if

       call vdm2d(v,faceR,faceS)
       call invVdm2D(v,inv,0,Npf)

       m = matmul(transpose(inv),inv)

       idr=fmask(:,iface)
       do i=1,Npf
          idc(i)=((iface-1)*Npf+1)+(i-1)
       end do
       emat(idr,idc) = emat(idr,idc)+m
    end do

    temp = matmul(transpose(vdm),emat)
    lift = matmul(vdm,temp)
  end subroutine lift3D
end module liftMod
