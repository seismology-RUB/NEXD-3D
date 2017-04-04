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
module simplexMod
! Evaluate orthonormal polynomial on simplex
	use constantsMod
	use jacobiMod
 implicit none

contains

  subroutine simplex2DP(Pr,a,b,i,j)
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: Pr
    real(kind=CUSTOM_REAL), dimension(:), intent(in):: a,b
    real(kind=CUSTOM_REAL), dimension(size(a)) :: h1,h2
    integer, intent(in) :: i,j
    real(kind=CUSTOM_REAL) :: zero,idvar1

    zero=0.0
    idvar1=2.0*i+1.0

    call jacobiP(h1,a,zero,zero,i)
    call jacobiP(h2,b,idvar1,zero,j)

    Pr(:) = sqrt(2.0)*h1*h2*(1.0-b)**i
  end subroutine simplex2DP

  subroutine simplex3DP(Pr,a,b,c,i,j,k)
    ! INPUT: node coordinates in reference quad a,b,c, orders i,j,k
    ! DOES: Calculate orthonormal polynomial on simplex at a,b,c
    ! RETURNS: orthonormal polynomial Pr
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: Pr
    real(kind=CUSTOM_REAL), dimension(:), intent(in):: a,b,c
    real(kind=CUSTOM_REAL), dimension(size(a)) :: h1,h2,h3
    integer, intent(in) :: i,j,k
    real(kind=CUSTOM_REAL) :: zero,idvar1,idvar2

    zero=0.0
    idvar1=2.0*i+1.0
    idvar2=2.0*(i+j)+2.0

    call jacobiP(h1,a,zero,zero,i)
    call jacobiP(h2,b,idvar1,zero,j)
    call jacobiP(h3,c,idvar2,zero,k)

    Pr(:) = 2*sqrt(2.0)*h1*h2*((1.0-b)**i)*h3*((1.0-c)**(i+j))
  end subroutine simplex3DP

  subroutine gradSimplex2DP(dPdr,dPds,a,b,id,jd)
    ! computes the differentiation matrices Dr and Ds
    real(kind=CUSTOM_REAL), dimension(:), intent(in) ::a,b
    integer, intent(in) :: id,jd
    real(kind=CUSTOM_REAL), dimension(size(a)), intent(out) :: dPdr,dPds
    real(kind=CUSTOM_REAL), dimension(size(a)) :: fa,dfa
    real(kind=CUSTOM_REAL), dimension(size(b)) :: gb,dgb
    real(kind=CUSTOM_REAL), dimension(size(a)) :: temp
    real(kind=CUSTOM_REAL) :: zero,idvar1

    zero=0.0
    idvar1=2.0*id+1.0

    call jacobiP(fa,a,zero,zero,id)
    call gradJacobiP(dfa,a,zero,zero,id)
    call jacobiP(gb,b,idvar1,zero,jd)
    call gradJacobiP(dgb,b,idvar1,zero,jd)
    
    dPdr(:)=dfa(:)*gb(:)
    if (id>0) then
       dPdr(:)=dPdr(:)*((0.5*(1.0-b(:)))**(id-1.0))
    end if

    dPds(:) = dfa(:)*(gb(:)*(0.5*(1.0+a(:))))
    if (id>0) then
       dPds(:)=dPds(:)*((0.5*(1.0-b(:)))**(id-1.0))
    end if

    temp(:) = dgb(:)*((0.5*(1.0-b(:)))**id)
    if (id>0) then
       temp(:)=temp(:)-0.5*id*gb(:)*((0.5*(1.0-b(:)))**(id-1.0))
    end if
    dPds(:) = dPds(:)+fa(:)*temp(:)

    !normalize
    dPdr(:) = 2.0**(id+0.5)*dPdr(:)
    dPds(:) = 2.0**(id+0.5)*dPds(:)
  end subroutine gradSimplex2DP
!
  subroutine gradSimplex3DP(dPdr,dPds,dPdt,a,b,c,id,jd,kd)
    ! computes the differentiation matrices Dr and Ds
    real(kind=CUSTOM_REAL), dimension(:), intent(in) ::a,b,c
    integer, intent(in) :: id,jd,kd
    real(kind=CUSTOM_REAL), dimension(size(a)), intent(out) :: dPdr,dPds,dPdt
    real(kind=CUSTOM_REAL), dimension(size(a)) :: fa,dfa
    real(kind=CUSTOM_REAL), dimension(size(b)) :: gb,dgb
    real(kind=CUSTOM_REAL), dimension(size(c)) :: hc,dhc
    real(kind=CUSTOM_REAL), dimension(size(a)) :: temp
    real(kind=CUSTOM_REAL) :: zero,idvar1,idvar2

    zero=0.0
    idvar1=2.0*id+1.0
    idvar2=2.0*(id+jd)+2.0

    call jacobiP(fa,a,zero,zero,id)
    call gradJacobiP(dfa,a,zero,zero,id)
    call jacobiP(gb,b,idvar1,zero,jd)
    call gradJacobiP(dgb,b,idvar1,zero,jd)
    call jacobiP(hc,c,idvar2,zero,kd)
    call gradJacobiP(dhc,c,idvar2,zero,kd)

    ! r-deriverate
    
    dPdr(:) = dfa(:)*(gb(:)*hc(:))
    if (id>0) then
       dPdr(:) = dPdr(:)*((0.5*(1.0-b(:)))**(id-1.0))
    end if
    if ((id+jd)>0) then
       dPdr(:) = dPdr(:)*((0.5*(1.0-c(:)))**(id+jd-1.0))
    end if
    
    ! s-deriverate

    dPds(:) = 0.5*(1.0+a(:))*dPdr(:)
    temp(:) = dgb(:)*((0.5*(1.0-b(:)))**id)
    if (id>0) then
       temp(:) = temp(:)+(-0.5*id)*(gb(:)*(0.5*(1.0-b(:)))**(id-1.0))
    end if
    if ((id+jd)>0) then
       temp(:) = temp(:)*((0.5*(1.0-c(:)))**(id+jd-1.0))
    end if
    temp(:) = fa(:)*temp(:)*hc(:)
    dPds(:) = dpds(:)+temp(:)
    
    ! t-deriverate
    
    dPdt(:) = 0.5 * (1.0+a(:)) * dPdr(:) + 0.5*(1.0+b(:))*temp(:)
    temp(:) = dhc(:) *((0.5*(1.0-c(:)))**(id+jd))
    if ((id+jd)>0) then
       temp(:) = temp(:) -0.5*(id+jd)*(hc*((0.5*(1.0-c(:)))**(id+jd-1.0)))
    end if
    temp(:) = fa(:) * (gb(:)*temp(:))
    temp(:) = temp(:) *((0.5*(1.0-b(:)))**id)

    dPdt(:) = dPdt(:) + temp(:)

    !normalize
    dPdr(:) = (2.0**(2.0*id+jd+1.5))*dPdr(:)
    dPds(:) = (2.0**(2.0*id+jd+1.5))*dPds(:)
    dPdt(:) = (2.0**(2.0*id+jd+1.5))*dPdt(:)
  end subroutine gradSimplex3DP
end module simplexMod
