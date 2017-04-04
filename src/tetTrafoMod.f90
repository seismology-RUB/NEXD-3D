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
module tetTrafoMod
! module to transform between different tri and tet representations
	use matrixMod
	use constantsMod
 implicit none
!
contains
!
  subroutine rsToAb(r,s,a,b)
    ! transform from reference triange into reference quad
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: a,b
    integer :: i

    do i=1,size(r)
       if (abs(s(i) - 1.) > EPS) then
          a(i) = 2.0*(1.0+r(i))/(1.0-s(i))-1.0
       else
          a(i) = -1.0
       end if
    end do
    b(:)=s(:)
  end subroutine rsToAb

  subroutine xyToRs(x,y,r,s)
    ! transform from equilateral tri into reference tri
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: r,s
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y
    real(kind=CUSTOM_REAL), dimension(Np) :: L1,L2,L3
    !
    L1(:) = (sqrt(3.0)*y(:)+1.0)/3.0
    L2(:) = (-3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
    L3(:) = ( 3.0*x(:) - sqrt(3.0) * y(:) + 2.0)/6.0
    r(:) = -L2(:) + L3(:) - L1(:)
    s(:) = -L2(:) - L3(:) + L1(:)
  end subroutine xyToRs

  subroutine xyzToRst(x,y,z,r,s,t)
    ! INPUT: node coordinates in tetrahedron x,y,z
    ! DOES: transform from equilateral tetrahedron into reference tetrahedron
    ! OUTPUT: node coordinates in reference tetrahedron r,s,t
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x,y,z
    real(kind=CUSTOM_REAL), dimension(3) :: v1,v2,v3,v4
    real(kind=CUSTOM_REAL), dimension(3,3) :: A
    real(kind=CUSTOM_REAL), dimension(3,Np) :: B,C

    A(:,:) = 0
    B(:,:) = 0
    C(:,:) = 0

    ! set vertices of the tetra
    v1 = (/-1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v2 = (/ 1.0, -1.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v3 = (/ 0.0,  2.0/sqrt(3.0), -1.0/sqrt(6.0) /)
    v4 = (/ 0.0,            0.0,  3.0/sqrt(6.0) /)


    B(1,:) = x(:) - 0.5*(v2(1)+v3(1)+v4(1)-v1(1)) 
    B(2,:) = y(:) - 0.5*(v2(2)+v3(2)+v4(2)-v1(2)) 
    B(3,:) = z(:) - 0.5*(v2(3)+v3(3)+v4(3)-v1(3)) 
    
    A(:,1) = 0.5*(v2-v1) 
    A(:,2) = 0.5*(v3-v1)
    A(:,3) = 0.5*(v4-v1)

    call invert(A)
    C = matmul(A,B)

    r= C(1,:)
    s= C(2,:)
    t= C(3,:)
  end subroutine xyzToRst

  subroutine rstToAbc(r,s,t,a,b,c)
    ! INPUT: node coordinates in reference tet r,s,t
    ! DOES: transform from reference tet into reference quad
    ! OUTPUT: node coordinates in reference quad a,b,c
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: a,b,c
    integer :: i

    do i=1,size(r)
       if (abs(s(i) + t(i)) > EPS) then
          a(i) = 2.0*(1.0+r(i))/(-s(i)-t(i))-1.0
       else
          a(i) = -1.0
       end if
       if (abs(t(i) - 1.) > EPS) then
          b(i) = 2.0*(1.0+s(i))/(1.0-t(i))-1.0
       else
          b(i) = -1.0
       end if
    end do
    c(:)=t(:)
  end subroutine rstToAbc
end module tetTrafoMod
