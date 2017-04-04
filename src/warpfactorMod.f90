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
module warpfactorMod
! module to calculate the warp factors for interpolation in the triangle with quasi gll points
	use constantsMod
	use gllMod
 implicit none

contains

  function warpfactor(xnodes,xout)
    ! INPUT: node locations xnodes, locations of nodes to interpolate warp function at (xout)
    ! DOES: calculates warp functions at positions xout
    ! RETURNS: warp function
    integer :: p
    real(kind=CUSTOM_REAL), dimension(NGLL) :: xnodes,xeq
    real(kind=CUSTOM_REAL), dimension(Np) :: warpfactor,d, warp ,xout
    real(kind=CUSTOM_REAL) :: dx,prod,prod1,prod2,e,sf
    logical :: zerof
    integer :: i,j,k, zeroff

    p=NGLL-1

! create equispaced vector in same ordering as xnodes
    dx=2./(NGLL-1)
    j=NGLL
    do i=1,NGLL
       xeq(i)=((j-1.0)*dx)-1.0
       j=j-1
    end do
    d(:)=0.
    do k=1,Np
       do i=1,NGLL
          e=(xnodes(i)-xeq(i))
          prod1=1.
          prod2=1.
          do j=1,NGLL
             if (i/=j) then
                prod1=prod1*(xout(k)-xeq(j))
                prod2=prod2*(xeq(i)-xeq(j))
             end if
          end do
          prod=e*(prod1/prod2)
          d(k)=d(k)+prod
       end do
       warp(k)=d(k)

! scale warp factor
       zerof = abs(xout(k)) < (1.0-EPS)
       if (zerof) then
          zeroff=1.0
       else
          zeroff=0.0
       end if
       sf = 1.0-(zeroff*xout(k))**2
       warpfactor(k)=warp(k)/sf + warp(k)*(zeroff-1.0)
    end do
  end function warpfactor

  subroutine shiftwarp(alpha,x,y,L1,L2,L3)
    ! INPUT: warping parameters alpha,
    ! DOES: performs blend and warp technique
    ! RETURNS: node coordinates of a triangle x,y, equidistributed node coordinates, L1, L2, L3
    real(kind=CUSTOM_REAL), intent(in) :: alpha
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: x,y
    real(kind=CUSTOM_REAL), dimension(Np) :: L1,L2,L3, blend1, blend2, blend3
    real(kind=CUSTOM_REAL), dimension(Np) :: warpfactor1,warpfactor2,warpfactor3
    real(kind=CUSTOM_REAL), dimension(Np) :: warp1,warp2,warp3
    real(kind=CUSTOM_REAL), dimension(NGLL) :: xgll
    integer :: p

    p=NGLL-1

    ! get GLL nodes
    xgll=getGLL()

    ! compute blending
    blend1(:) = 4.0*L2(:)*L3(:)
    blend2(:) = 4.0*L1(:)*L3(:)
    blend3(:) = 4.0*L1(:)*L2(:)

    !get warpfactors
    warpfactor1 = warpfactor(xgll, L3(:)-L2(:))
    warpfactor2 = warpfactor(xgll, L1(:)-L3(:))
    warpfactor3 = warpfactor(xgll, L2(:)-L1(:))

    ! combine blend and warp
    warp1(:) = blend1(:)*warpfactor1(:)*(1.0+(alpha*L1(:))**2)
    warp2(:) = blend2(:)*warpfactor2(:)*(1.0+(alpha*L2(:))**2)
    warp3(:) = blend3(:)*warpfactor3(:)*(1.0+(alpha*L3(:))**2)

    ! accumulate deformations associated with each edge
    x(:) =  1.0*warp1(:) + cos(2.0*PI/3.0)*warp2(:) + cos(4.0*PI/3.0)*warp3(:)
    y(:) =  0.0*warp1(:) + sin(2.0*PI/3.0)*warp2(:) + sin(4.0*PI/3.0)*warp3(:)

  end subroutine shiftwarp
!
end module warpfactorMod
