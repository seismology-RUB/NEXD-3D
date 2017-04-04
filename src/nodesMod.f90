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
module nodesMod
! module to create 'warp and blend' nodes 3d following T.Warburton 2006
	use constantsMod
	use gllMod
	use warpfactorMod
 implicit none

contains

  subroutine nodes3D(alpha,x,y,z)

    ! INPUT: warping parameters alpha
    ! DOES: performs warp and blend technique suggested by Warburton: An explicit construction of interpolation nodes
    ! on the simplex, 2006
    ! RETURNS: node coordinates in tetrahedron x,y,z
    real(kind=CUSTOM_REAL), intent(in) :: alpha
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: x,y,z
    real(kind=CUSTOM_REAL), dimension(Np) :: L1,L2,L3,L4,r,s,t
    real(kind=CUSTOM_REAL), dimension(Np) :: La,Lb,Lc,Ld
    real(kind=CUSTOM_REAL), dimension(Np) :: warp1,warp2
    real(kind=CUSTOM_REAL), dimension(Np) :: blend,denom
    integer, dimension(3) :: bool
    integer :: p=order
    real(kind=CUSTOM_REAL) :: tol=1e-10
    integer :: sk,n,m,q,i,j
    real(kind=CUSTOM_REAL), dimension(3) :: v1,v2,v3,v4
    real(kind=CUSTOM_REAL), dimension(Np,3) :: XYZ,shift
    real(kind=CUSTOM_REAL), dimension(4,3) :: t1,t2

    !make eqidistant point in tet
    sk = 1
    do n=1,p+1
        do m=1,p+2-n
            do q=1,p+3-n-m
                r(sk) = -1.0 + (q-1.0)*2.0/p
                s(sk) = -1.0 + (m-1.0)*2.0/p
                t(sk) = -1.0 + (n-1.0)*2.0/p
                sk = sk+1
            end do
        end do
    end do

    L1(:) = (1.0+t(:))/2.0
    L2(:) = (1.0+s(:))/2.0
    L3(:) = -(1.0+r(:)+s(:)+t(:))/2.0
    L4(:) = (1.0+r(:))/2.0
    
    ! set vertices of the tetra
    v1 = (/-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v2 = (/ 1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v3 = (/ 0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v4 = (/ 0.0,            0.0,  3.0/sqrt(6.0) /)

    ! set orthogonal axis tengents on the faces
    t1(1,:) = v2(:) - v1(:)
    t1(2,:) = v2(:) - v1(:)
    t1(3,:) = v3(:) - v2(:)
    t1(4,:) = v3(:) - v1(:)
    t2(1,:) = v3(:) - 0.5*(v1(:)+v2(:))
    t2(2,:) = v4(:) - 0.5*(v1(:)+v2(:))
    t2(3,:) = v4(:) - 0.5*(v2(:)+v3(:))
    t2(4,:) = v4(:) - 0.5*(v1(:)+v3(:))

    ! normalise tangents
    do i=1,4
        t1(i,:) = t1(i,:)/sqrt(t1(i,1)**2 + t1(i,2)**2 + t1(i,3)**2)
        t2(i,:) = t2(i,:)/sqrt(t2(i,1)**2 + t2(i,2)**2 + t2(i,3)**2)
    end do
    
    ! warp and blend for each face
    shift(:,:) = 0.0
    
    do i=1,Np
        do j=1,3
            XYZ(i,j)=L3(i)*v1(j) + L4(i)*v2(j) + L2(i)*v3(j) + L1(i)*v4(j)
        end do
    end do
    
    do i=1,4
        if (i == 1) then
            La = L1
            Lb = L2
            Lc = L3
            Ld = L4
        else if( i == 2 ) then
            La = L2
            Lb = L1
            Lc = L3
            Ld = L4
        else if( i == 3 ) then
            La = L3
            Lb = L1
            Lc = L4
            Ld = L2
        else if( i == 4 ) then
            La = L4
            Lb = L1
            Lc = L3
            Ld = L2
        end if
        ! compute warp tangential to face
        call shiftwarp(alpha,warp1,warp2,Lb,Lc,Ld)

        ! volume blending
        blend(:) = Lb(:)*Lc(:)*Ld(:)
        denom(:) = (Lb(:)+0.5*La(:))*(Lc(:)+0.5*La(:))*(Ld(:)+0.5*La(:))
        do j=1,Np
            if (denom(j) > tol) then
                blend(j) = ( 1.0 + (alpha*La(j))**2 ) * blend(j)/denom(j)
            end if
        end do

        !compute warp and shift
        do j=1,3
            shift(:,j) = shift(:,j) + (blend(:) *warp1) * t1(i,j) + (blend(:) *warp2) * t2(i,j)
        end do

        ! fix face warp
        do j=1,Np
            bool(:)=0
            if (Lb(j) > tol) bool(1)=1
            if (Lc(j) > tol) bool(2)=1
            if (Ld(j) > tol) bool(2)=1
            if ( (La(j) < tol) .AND. ( sum(bool)<3 )) then
                do n=1,3
                    shift(j,n) = warp1(j) *t1(i,n) + warp2(j) *t2(i,n)
                end do
            end if
        end do
    end do
    XYZ = XYZ + shift
    x=XYZ(:,1)
    y=XYZ(:,2)
    z=XYZ(:,3)
    
  end subroutine nodes3D

end module nodesMod

