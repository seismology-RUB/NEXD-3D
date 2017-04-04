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
module anelasticMod
    use constantsMod
    use linearSystemMod
    implicit none

contains

  subroutine anelasticPreCalc(f0,fr,wl,wk)
    ! For information and equations see Kaeser et al. (2006): An arbitrary high-order Discontinous Galerkin
    ! method for elastic waves on unstructured  meshes - III. Viscoelastic attenuation
    !INPUT: center and lower frequency of attenuation model f0 and fr
    !DOES:  calculates relaxation frequencies for general Maxwell body approach
    !RETURNS:   relaxation frequencies wl and wk
    integer :: m,c,i
    real(kind=CUSTOM_REAL) :: wmax,wmin,f0,fr                     ! maximum and minimum frequency
    real(Kind=CUSTOM_REAL), dimension(:) :: wl, wk

    m = 2*nMB-1                 ! nMB=number of Maxwell bodies, m=number of values for wk (frequencies) to obtain
                                ! constant quality factor (see Kaeser et al.)
    wmax = 2.5*pi*f0            ! maximum and minimum frequency of wk use f0 and fr from Parfile
    wmin = f0/sqrt(fr)

    if (nMB > 1) then               ! setup equidistant (log) frequencies in the requestet frequency band
        do i=1,m
            wk(i) = exp(log(wmin) + (i-1.0)/(2.0*(nMB-1)) * log(fr) )
        end do
    else
        wk(:)=wmax
    end if
    c=1
    do i=1,m
        if(mod(i,2)==1) then
            wl(c)=wk(i)           ! get the relaxation frequencies of the nMB mechanisms (see text between eq. 5 and 6 of Kaeser et al.)
            c=c+1
        end if
    end do

  end subroutine anelasticPreCalc

  subroutine calcAnelasticCoefficients(Qp,Qs,vp,vs,rho,f0,fr,wl,wk,ylambda,ymu,lambda_u,mu_u)
    ! For information and equations see Kaeser et al. (2006): An arbitrary high-order Discontinous Galerkin
    ! method for elastic waves on unstructured  meshes - III. Viscoelastic attenuation
    !INPUT: for every material: Quality factors Qp and Qs, velocities at center frequency vp and vs,
    !       density rho, center and lower frequency of attenuation model f0 and fr
    !DOES:  calculation of anelastic coefficients and parameters for a general Maxwell Body approach
    !RETURNS:   relaxation frequencies wl, anelastic coefficients lambda and ymu, unrelaxed Lamé parameters lambda_u and mu_u

    real(kind=CUSTOM_REAL) :: Qp,Qs,vp,vs,rho,f0,fr         ! see input
    real(kind=CUSTOM_REAL), dimension(:) :: ylambda,ymu     ! see returns
    real(kind=CUSTOM_REAL) :: lambda_u, mu_u
    real(Kind=CUSTOM_REAL), dimension(:) :: wl, wk

    integer :: i,j,m
    real(Kind=CUSTOM_REAL), dimension(:), allocatable :: invQp, invQs   ! used frequencies, inverse Qp and Qs at frequencies wk
    real(Kind=CUSTOM_REAL), dimension(:,:), allocatable :: Mp,Ms            ! Matrices for calculation of anelastic coefficients
    real(kind=CUSTOM_REAL), dimension(nMB) :: yp,ys
    real(kind=CUSTOM_REAL) :: theta1, theta2,R
    real(kind=CUSTOM_REAL) :: Mup,Mus


    ! Only necesssary to plot velocities and quality factors
    !integer :: dim=1000         ! number of datapoints to plot test functions for Qp and vp
    !real, dimension(1000) :: qtest, wtest
    !real :: temp1, temp2,vptest,vstest,vtemp,wtemp


!    set up matrix of the system where we want to obtain the anelastic constants for the wave speeds
    m = 2*nMB-1
    allocate(Mp(m,nMB))
    allocate(Ms(m,nMB))
    allocate(invQp(m), invQs(m))     ! prepare frequencies for determination and inverse values for Qp and Qs

    do i=1,m
        do j=1,nMB
             Mp(i,j) =  ( wl(j)*wk(i)+wl(j)**2/Qp ) / (wl(j)**2+wk(i)**2)  ! see eq. 6 in Kaeser et al.
             Ms(i,j) =  ( wl(j)*wk(i)+wl(j)**2/Qs ) / (wl(j)**2+wk(i)**2)
         end do
    end do

    invQp(:) = 1.0/Qp       ! define inverse quality factors at every used frequency wk
    invQs(:) = 1.0/Qs
    call solveLinearSystemQR(Mp,invQp,yp,m,nMB) ! solve the system invQp=Mp for yp (same for ys)
    call solveLinearSystemQR(Ms,invQs,ys,m,nMB)

    ! get the unrelaxed Lame parameters mu_u and lambda_u

    ! P waves
    theta1=1.0
    theta2=0.0
    do j=1,nMB
        theta1 = theta1 - yp(j)/(1+(f0/wl(j))**2)
        theta2 = theta2 + yp(j)*(f0/wl(j))/(1+(f0/wl(j))**2)
    end do
    R = sqrt(theta1**2+theta2**2)
    Mup = rho * vp**2 * (R + theta1)/(2*R**2)

    ! S waves
    theta1=1.0
    theta2=0.0
    do j=1,nMB
        theta1 = theta1 - ys(j)/(1+(f0/wl(j))**2)
        theta2 = theta2 + ys(j)*(f0/wl(j))/(1+(f0/wl(j))**2)
    end do
    R = sqrt(theta1**2+theta2**2)
    Mus = rho * vs**2 * (R + theta1)/(2*R**2)

    ! Lame parameters:
    mu_u = Mus
    lambda_u = Mup-2*Mu_u
    ! get the anelastic coefficients for the Lame parameters
    do j=1,nMB
         ylambda(j) = (1+ (2 * mu_u) / lambda_u) * yp(j) - (2 * mu_u)/lambda_u *ys(j) ! eq. 7 Kaeser et al.
         ymu(j)= ys(j)
    end do

    ! Here, the calculations are finished. The remaining parts provide two files to print Qp and vp as a function of frequency
    ! Can be commented out if not needed
    !do i=1,dim
    !    wtest(i)=(i/2)**2        ! frequencies at which Qp, vp and vs will be plotted
    !end do

    !open(unit=27,file="qtest",status='unknown')
    !do i=1,dim          ! calculate for every wtest Qp
    !    temp1=0
    !    temp2=0
    !    do j=1,nMB
    !        temp1=temp1 + ( wl(j)**2 * yp(j) ) / ( wl(j)**2 + wtest(i)**2 )
    !        temp2=temp2 + ( wl(j)*wtest(i) * yp(j) ) / ( wl(j)**2 + wtest(i)**2 )
    !    end do
    !    qtest(i) = (1.0-temp1)/temp2
    !    write(27,*) wtest(i),qtest(i)
    !end do
    !close(27)

    !open(unit=27,file="vptest",status='unknown')
    !vtemp=99999.
    !do i=1,dim
    !    itemp1=0
    !    itemp2=0
    !    do j=1,nMB
    !        itemp1=itemp1 + wl(j)*ylambda(j)/(wl(j)+cmplx(0.0,1.0)*wtest(i))
    !        itemp2=itemp2 + wl(j)*ymu(j)/(wl(j)+cmplx(0.0,1.0)*wtest(i))
    !    end do
    !    iltemp = lambda_u*(1-itemp1)
    !    imtemp = mu_u*(1-itemp2)
    !    vptest = real(sqrt((iltemp+2*imtemp)/rho))
    !    vstest = real(sqrt(imtemp/rho))
    !    write(27,*) wtest(i),vptest,vstest
    !    if (vstest<vtemp) then
    !        vtemp=vstest
    !        wtemp=wtest(i)
    !    endif
    !end do
    !write(*,*) 'vs_min: ',vtemp,' at w: ',wtemp
    !close(27)

    deallocate(Mp,Ms)
    deallocate(invQp, invQs)

  end subroutine calcAnelasticCoefficients
end module anelasticMod
