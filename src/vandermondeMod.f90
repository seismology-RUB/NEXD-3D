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
module vandermondeMod
! module to create the vandermonde matrix
	use constantsMod
	use simplexMod
	use tetTrafoMod
	use jacobiMod
	use matrixMod
  implicit none

contains
  subroutine vdm2d(v,r,s)
    ! initialize the 2D vandermonde matrix
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s
    real(kind=CUSTOM_REAL), dimension(size(r),size(r)), intent(out) :: v
    real(kind=CUSTOM_REAL), dimension(size(r)) :: a,b
    integer :: i,j,sk
    !
    call rstoab(r,s,a,b)
    sk=1
    do i=0,NGLL-1
       do j=0,NGLL-1-i
          call simplex2DP(v(:,sk),a,b,i,j)
          sk=sk+1
       end do
    end do
  end subroutine vdm2d

  subroutine invVdm2D(v,w,trans,n)
    ! calulates v^-1' if trans==0 -> inverse if trans == 1 => inverse transponierte
    integer, intent(in) ::n
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
    real(kind=CUSTOM_REAL), dimension(size(v(:,1)),n) :: w
    integer, intent(in) :: trans

    w=v
    call invert(w)

    if (trans==1) then
       w=transpose(w)
    end if

  end subroutine invVdm2D

  subroutine vdm3d(v,r,s,t)
    ! INPUT: node coordinates in reference tetrahedron r,s,t
    ! DOES: initialize the 3D vandermonde matrix
    ! OUTPUT: van-der-monde matrix v
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(size(r),Np), intent(out) :: v
    real(kind=CUSTOM_REAL), dimension(size(r)) :: a,b,c
    integer :: i,j,k,sk
    !
    call rstToAbc(r,s,t,a,b,c)  ! transform into reference tet

    sk=1
    do i=0,NGLL-1
       do j=0,NGLL-1-i
          do k=0,NGLL-1-i-j
             call simplex3DP(v(:,sk),a,b,c,i,j,k)
             sk=sk+1
          end do
       end do
    end do
  end subroutine vdm3d
!
  subroutine invVdm3D(v,w,trans)
    ! INPUT: van-der-monde matrix v, flag trans
    ! DOES: calulates v^-1', if trans==0 -> inverse; if trans == 1 => inverse transpose
    ! OUTPUT: inverse (transpose) van-der-monde matrix w
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: v
    real(kind=CUSTOM_REAL), dimension(size(v(:,1)),Np) :: w
    integer, intent(in) :: trans

    w=v
    call invert(w)
    if (trans==1) then
        w=transpose(w)
    end if
  end subroutine invVdm3D
!
  subroutine gradVdm3D(v3dr,v3ds,v3dt,r,s,t)
    ! INPUT: node coordinates reference tet r,s,t
    ! DOES: build the gradient of the modal basis (i,j,k) at (r,s,t) in referenze tet
    ! OUTPUT: gradients to basis r,s,t (v3dr, v3ds, v3dt)
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: r,s,t
    real(kind=CUSTOM_REAL), dimension(size(r),Np), intent(out) :: v3dr,v3ds,v3dt
    real(kind=CUSTOM_REAL), dimension(size(r)) :: a,b,c
    integer :: i,j,k,sk
    !
    call rstToabc(r,s,t,a,b,c)
    
    sk=1
    do i=0,NGLL-1
       do j=0,NGLL-1-i
          do k=0,NGLL-1-i-j
             call gradSimplex3DP(v3dr(:,sk),v3ds(:,sk),v3dt(:,sk),a,b,c,i,j,k)
             sk=sk+1
          end do
       end do
    end do
  end subroutine gradVdm3D
end module vandermondeMod
