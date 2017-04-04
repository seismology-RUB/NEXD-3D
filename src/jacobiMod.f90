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
module jacobiMod
! module to deal with jacobi polinomials
	use constantsMod
	use rosettaGammaMod
 implicit none

contains

  subroutine jacobiP(P,x,alpha,beta,N)
    ! INPUT: x,,alpha,beta,N: Parameters of Jacobi polynomials
    ! DOES: Calculates the Jacobi polynomials
    ! OUTPUT: P: Array with Jacobi polynomial
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x
    integer, intent(in) :: N
    real(kind=CUSTOM_REAL), dimension(N+1,size(x)) :: PL
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: P
    real(kind=CUSTOM_REAL) :: gamma0, gamma1, aold, anew, h1, bnew
    real(kind=CUSTOM_REAL), intent(in) :: alpha, beta
    integer :: dim, i
    
    dim=size(x)
    gamma0 = 2**(alpha+beta+1.0)/(alpha+beta+1.0)*lacz_gamma(real(alpha+1.0))*lacz_gamma(real(beta+1.0))/lacz_gamma(real(alpha+beta+1.0))
    PL(1,:) = 1.0/sqrt(gamma0)

    if (N == 0) then
        P=PL(1,:)
        return
    end if

    gamma1 = (alpha+1.0)*(beta+1.0)/(alpha+beta+3.0)*gamma0
    PL(2,:) = ((alpha+beta+2.0)*x(:)/2.0 + (alpha-beta)/2.0) / sqrt(gamma1)

    if (N == 1) then
        P=PL(N+1,:)
        return
    end if
    
    ! repeat value in recurrence
    aold = 2.0/(2.0+alpha+beta)*sqrt((alpha+1.0)*(beta+1.0)/(alpha+beta+3.0))
    
    ! forward recurrence using the symetry of the recurrence
    do i=1,N-1
        h1 = 2.0*i+alpha+beta
        anew = 2.0/(h1+2.0)*sqrt( (i+1.0) * (i+1.0+alpha+beta)*(i+1.0+alpha)*(i+1.0+beta)/(h1+1.0)/(h1+3.0) )
        bnew = - (alpha**2-beta**2)/h1/(h1+2.0)
        PL(i+2,:) = 1.0/anew*( -aold*PL(i,:) + (x(:)-bnew)*PL(i+1,:))
        aold=anew
    end do

    P=PL(N+1,:)
  end subroutine jacobiP

  subroutine gradJacobiP(dP,x,alpha,beta,N)
    ! INPUT: x,,alpha,beta,N: Parameters of Jacobi gradient
    ! DOES: compute grad of jacobi polynomials
    ! OUTPUT: P: Array with Jacobi gradient
    real(kind=CUSTOM_REAL), dimension(:), intent(out) :: dP
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: x
    real(kind=CUSTOM_REAL), intent(in) :: alpha, beta
    real(kind=CUSTOM_REAL), dimension(size(dp)) ::h1
    integer, intent(in) :: N
    if (N==0) then
        dP(:)=0
    else
        call jacobiP(h1,x,alpha+1,beta+1,N-1)
        dp(:) = sqrt(N*(N+alpha+beta+1.0))*h1(:)
    end if
  end subroutine gradJacobiP

end module jacobiMod
