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
module waveMod
  ! calculate elastic fluxes
	use constantsMod
	use matrixMod
	use mpiMod
  implicit none

contains

  subroutine elasticFluxes(q,elasticfluxvar,f,g,h)
    ! INPUT: q: simulation parameters (stresses, velocities), lambda, mu: Lamé constants, rho: density
    ! DOES: calculate elastic fluxes
    ! RETURNS: f,g,h: elastic fluxes for x,y,z [1, eq. 2.97 to 2.99 multiplicated with q, compare eq. 2.129]
    real(kind=CUSTOM_REAL), dimension(4) :: elasticfluxvar
    real(kind=CUSTOM_REAL), dimension(:,:) :: f,g,h
    real(kind=CUSTOM_REAL), dimension(:,:) :: q

    f(:,1)= elasticfluxvar(1)*q(:,7)
    f(:,2)= elasticfluxvar(2)*q(:,7)
    f(:,3)= elasticfluxvar(2)*q(:,7)
    f(:,4)= elasticfluxvar(3)*q(:,8)
    f(:,6)= elasticfluxvar(3)*q(:,9)
    f(:,7)= elasticfluxvar(4)*q(:,1)
    f(:,8)= elasticfluxvar(4)*q(:,4)
    f(:,9)= elasticfluxvar(4)*q(:,6)

    g(:,1)= elasticfluxvar(2)*q(:,8)
    g(:,2)= elasticfluxvar(1)*q(:,8)
    g(:,3)= elasticfluxvar(2)*q(:,8)
    g(:,4)= elasticfluxvar(3)*q(:,7)
    g(:,5)= elasticfluxvar(3)*q(:,9)
    g(:,7)= elasticfluxvar(4)*q(:,4)
    g(:,8)= elasticfluxvar(4)*q(:,2)
    g(:,9)= elasticfluxvar(4)*q(:,5)

    h(:,1)= elasticfluxvar(2)*q(:,9)
    h(:,2)= elasticfluxvar(2)*q(:,9)
    h(:,3)= elasticfluxvar(1)*q(:,9)
    h(:,5)= elasticfluxvar(3)*q(:,8)
    h(:,6)= elasticfluxvar(3)*q(:,7)
    h(:,7)= elasticfluxvar(4)*q(:,6)
    h(:,8)= elasticfluxvar(4)*q(:,5)
    h(:,9)= elasticfluxvar(4)*q(:,3)
  end subroutine elasticFluxes
  
  subroutine anelasticFluxes(q,wl,f,g,h)
    ! INPUT: q: simulation parameters (stresses, velocities), wl: attenuation frequencies
    ! DOES: calculate anelastic fluxes
    ! RETURNS: f,g,h: anelastic fluxes for x,y,z [1, eq. 2.88 for 2D case, matrices are multiplicated by q]
  real(kind=CUSTOM_REAL), dimension(nMB) :: wl
  real(kind=CUSTOM_REAL), dimension(Np,3) :: temp1, temp2
  real(kind=CUSTOM_REAL), dimension(:,:) :: f,g,h
  real(kind=CUSTOM_REAL), dimension(:,:) :: q
  integer :: i
  do i=1,nMB
    temp1(:,:) = -0.5*wl(i)*q(:,7:9)
    temp2(:,:) = -wl(i)*q(:,7:9)

    f(:,(i-1)*6+1) = temp2(:,1)
    f(:,(i-1)*6+4) = temp1(:,2)
    f(:,(i-1)*6+6) = temp1(:,3)

    g(:,(i-1)*6+2) = temp2(:,2)
    g(:,(i-1)*6+4) = temp1(:,1)
    g(:,(i-1)*6+5) = temp1(:,3)

    h(:,(i-1)*6+3) = temp2(:,3)
    h(:,(i-1)*6+5) = temp1(:,2)
    h(:,(i-1)*6+6) = temp1(:,1)
  end do
    
  end subroutine anelasticFluxes    
    
  !
  function getAPA(vp,vs,rho,lambda,mu)
    ! INPUT: vp, vs: seismic velocities, rho: density, lambda, mu: Lamé constants
    ! DOES:  compute A + |A| for exact riemann fluxes
    ! RETURNS: A + |A| [1,eq. 2.55]
    real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
    real(kind=CUSTOM_REAL), dimension(9,9) :: getAPA
    getAPA=0.0
    getAPA(1,1) = vp
    getAPA(1,7) = -(lambda+2*mu)
    getAPA(2,1) = (lambda*vp)/(lambda+2*mu)
    getAPA(2,7) = -lambda
    getAPA(3,1) = (lambda*vp)/(lambda+2*mu)
    getAPA(3,7) = -lambda
    getAPA(4,4) = vs
    getAPA(4,8) = -mu
    getAPA(6,6) = vs
    getAPA(6,9) = -mu
    getAPA(7,1) = -1/rho
    getAPA(7,7) = vp
    getAPA(8,4) = -1/rho
    getAPA(8,8) = vs
    getAPA(9,6) = -1/rho
    getAPA(9,9) = vs

  end function getAPA
  !
  function getAMA(vp,vs,rho,lambda,mu)
    ! INPUT: vp, vs: seismic velocities, rho: density, lambda, mu: Lamé constants
    ! DOES:  compute A - |A| for exact riemann fluxes
    ! RETURNS: A - |A| [1,eq. 2.55]
    real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
    real(kind=CUSTOM_REAL), dimension(9,9) :: getAMA
    getAMA=0.0
    getAMA(1,1) = -vp
    getAMA(1,7) = -(lambda+2*mu)
    getAMA(2,1) = -(lambda*vp)/(lambda+2*mu)
    getAMA(2,7) = -lambda
    getAMA(3,1) = -(lambda*vp)/(lambda+2*mu)
    getAMA(3,7) = -lambda
    getAMA(4,4) = -vs
    getAMA(4,8) = -mu
    getAMA(6,6) = -vs
    getAMA(6,9) = -mu
    getAMA(7,1) = -1/rho
    getAMA(7,7) = -vp
    getAMA(8,4) = -1/rho
    getAMA(8,8) = -vs
    getAMA(9,6) = -1/rho
    getAMA(9,9) = -vs
  end function getAMA
  
  
  function getAnelasticAPA(vp,vs,rho,lambda,mu)
    ! INPUT: vp, vs: seismic velocities, rho: density, lambda, mu: Lamé constants
    ! DOES:  compute anelastiv part of A + |A| for exact riemann fluxes
    ! RETURNS: anelastic A + |A| [1,eq. 2.88 (without omega) and eq. 2.93 for 2D case]
    real(kind=CUSTOM_REAL) :: vp,vs,rho,lambda,mu
    real(kind=CUSTOM_REAL), dimension(6,9) :: getAnelasticAPA
    getAnelasticAPA=0.0
    getAnelasticAPA(1,1) = 1.0/(vp*rho)
    getAnelasticAPA(1,7) = -1.0
    getAnelasticAPA(4,4) = 1.0/(2*vs*rho)
    getAnelasticAPA(4,8) = -1.0/2.0
    getAnelasticAPA(6,6) = 1.0/(2*vs*rho)
    getAnelasticAPA(6,9) = -1.0/2.0
  end function getAnelasticAPA
!
  function getT(n,s,t)
    ! INPUT: n,s,t: normal components
    ! DOES: calculates rotation matrix
    ! RETURNS: getT: rotation matrix [1, eq: 2.131]
    real(kind=CUSTOM_REAL), dimension(:) :: n,s,t 
    real(kind=CUSTOM_REAL), dimension(9,9) :: getT
    getT=0.0
    getT(1,1) = n(1)*n(1)
    getT(1,2) = s(1)*s(1)
    getT(1,3) = t(1)*t(1)
    getT(1,4) = 2*n(1)*s(1)
    getT(1,5) = 2*s(1)*t(1)
    getT(1,6) = 2*n(1)*t(1)
    getT(2,1) = n(2)*n(2)
    getT(2,2) = s(2)*s(2)
    getT(2,3) = t(2)*t(2)
    getT(2,4) = 2*n(2)*s(2)
    getT(2,5) = 2*s(2)*t(2)
    getT(2,6) = 2*n(2)*t(2)
    getT(3,1) = n(3)*n(3)
    getT(3,2) = s(3)*s(3)
    getT(3,3) = t(3)*t(3)
    getT(3,4) = 2*n(3)*s(3)
    getT(3,5) = 2*s(3)*t(3)
    getT(3,6) = 2*n(3)*t(3)
    getT(4,1) = n(2)*n(1)
    getT(4,2) = s(2)*s(1)
    getT(4,3) = t(2)*t(1)
    getT(4,4) = n(2)*s(1)+n(1)*s(2)
    getT(4,5) = s(2)*t(1)+s(1)*t(2)
    getT(4,6) = n(2)*t(1)+n(1)*t(2)
    getT(5,1) = n(3)*n(2)
    getT(5,2) = s(3)*s(2)
    getT(5,3) = t(3)*t(2)
    getT(5,4) = n(3)*s(2)+n(2)*s(3)
    getT(5,5) = s(3)*t(2)+s(2)*t(3)
    getT(5,6) = n(3)*t(2)+n(2)*t(3)
    getT(6,1) = n(3)*n(1)
    getT(6,2) = s(3)*s(1)
    getT(6,3) = t(3)*t(1)
    getT(6,4) = n(3)*s(1)+n(1)*s(3)
    getT(6,5) = s(3)*t(1)+s(1)*t(3)
    getT(6,6) = n(3)*t(1)+n(1)*t(3)
    getT(7,7) = n(1)
    getT(7,8) = s(1)
    getT(7,9) = t(1)
    getT(8,7) = n(2)
    getT(8,8) = s(2)
    getT(8,9) = t(2)
    getT(9,7) = n(3)
    getT(9,8) = s(3)
    getT(9,9) = t(3)
  end function getT
!
  function getInvT(n,s,t)
    ! INPUT: n,s,t: normal components
    ! DOES: calculates inverse rotation matrix
    ! RETURNS: getinvT: rotation matrix [1, eq: 2.132]
    real(kind=CUSTOM_REAL), dimension(:) :: n,s,t 
    real(kind=CUSTOM_REAL), dimension(9,9) :: getInvT
    getInvT=0.0
    getInvT(1,1) = n(1)*n(1)
    getInvT(1,2) = n(2)*n(2)
    getInvT(1,3) = n(3)*n(3)
    getInvT(1,4) = 2*n(2)*n(1)
    getInvT(1,5) = 2*n(3)*n(2)
    getInvT(1,6) = 2*n(3)*n(1)
    getInvT(2,1) = s(1)*s(1)
    getInvT(2,2) = s(2)*s(2)
    getInvT(2,3) = s(3)*s(3)
    getInvT(2,4) = 2*s(2)*s(1)
    getInvT(2,5) = 2*s(3)*s(2)
    getInvT(2,6) = 2*s(3)*s(1)
    getInvT(3,1) = t(1)*t(1)
    getInvT(3,2) = t(2)*t(2)
    getInvT(3,3) = t(3)*t(3)
    getInvT(3,4) = 2*t(2)*t(1)
    getInvT(3,5) = 2*t(3)*t(2)
    getInvT(3,6) = 2*t(3)*t(1)
    getInvT(4,1) = n(1)*s(1)
    getInvT(4,2) = n(2)*s(2)
    getInvT(4,3) = n(3)*s(3)
    getInvT(4,4) = n(2)*s(1)+n(1)*s(2)
    getInvT(4,5) = n(3)*s(2)+n(2)*s(3)
    getInvT(4,6) = n(3)*s(1)+n(1)*s(3)
    getInvT(5,1) = s(1)*t(1)
    getInvT(5,2) = s(2)*t(2)
    getInvT(5,3) = s(3)*t(3)
    getInvT(5,4) = s(2)*t(1)+s(1)*t(2)
    getInvT(5,5) = s(3)*t(2)+s(2)*t(3)
    getInvT(5,6) = s(3)*t(1)+s(1)*t(3)
    getInvT(6,1) = n(1)*t(1)
    getInvT(6,2) = n(2)*t(2)
    getInvT(6,3) = n(3)*t(3)
    getInvT(6,4) = n(2)*t(1)+n(1)*t(2)
    getInvT(6,5) = n(3)*t(2)+n(2)*t(3)
    getInvT(6,6) = n(3)*t(1)+n(1)*t(3)
    getInvT(7,7) = n(1)
    getInvT(7,8) = n(2)
    getInvT(7,9) = n(3)
    getInvT(8,7) = s(1)
    getInvT(8,8) = s(2)
    getInvT(8,9) = s(3)
    getInvT(9,7) = t(1)
    getInvT(9,8) = t(2)
    getInvT(9,9) = t(3)
  end function getInvT
!
  subroutine computeExactRiemannSF(ie,myrank,flux,q,qi,neighbor,VT,VTfree,face, mpi_interface,&
       mpi_connection,mpi_ibool, mpi_ibt, mpi_icon,ibt,ibn)
    ! INPUT: see timeloop for describtion of parameters, nn, st, tt are the normal components
    ! DOES: computes the Riemann flux of an element [1, eq. 2.130 with modifications for free/absorbing boundaries]
    ! RETURNS: flux: Riemann flux
    integer, dimension(:) :: neighbor,face
    integer, dimension(:,:) :: ibt,ibn
    real(kind=CUSTOM_REAL), dimension(4*NpF,9) :: flux
    real(kind=CUSTOM_REAL), dimension(:,:) :: q
    real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
    integer, dimension(:,:) :: mpi_interface
    integer, dimension(:,:,:) :: mpi_connection
    integer, dimension(:) :: mpi_ibool
    integer, dimension(:,:) :: mpi_ibt
    integer, dimension(:) :: mpi_icon
    real(kind=CUSTOM_REAL), dimension(Npf,9) :: qtemp_n,qtemp_t
    real(kind=CUSTOM_REAL), dimension(4,9,9) :: VT,VTfree
    integer :: is,in,i,j,ie
    integer :: mpi_e, mpi_n
    integer :: myrank
    flux=0.0
    do is=1,4                               ! for all surfaces
        in = neighbor(is)                   ! get neigbor element
        if (in>0) then                      ! not a boundary
            qtemp_t=q(ibt(:,is),:)          ! copy simulation parameter for element and neighbor
            qtemp_n=q(ibn(:,is),:)
            do i=1,NpF
                qtemp_t(i,:) = matmul(VT(is,:,:),(qtemp_n(i,:)-qtemp_t(i,:)))       ! apply remaining parts of [1,eq. 2.130] WAS IST MIT DEM 1/2??? AL
            end do
            do j=1,9
                flux(((is-1)*NpF+1):(is*NpF),j) = qtemp_t(:,j)              ! make flux vector for all points of all surfaces
            end do
        ! boundaries
        ! "------------------------------------------------------------------"
        else if ((in==0).and.face(is) == -1) then                           ! free surface element, works as above, only VT is different
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                qtemp_t(i,:) = matmul(VTfree(is,:,:),qtemp_t(i,:))
            end do
            do j=1,9
                flux(((is-1)*NpF+1):(is*NpF),j) = qtemp_t(:,j)
            end do
        ! "------------------------------------------------------------------"
        else if ((in==0).and.face(is) == -2) then                           ! absorbing surface element, qtemp_t is negative for absorbtion
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                qtemp_t(i,:) = -matmul(VT(is,:,:),qtemp_t(i,:))
            end do
            do j=1,9
                flux(((is-1)*NpF+1):(is*NpF),j) = qtemp_t(:,j)
            end do
        else if (in == -1) then                                             ! MPI interface flux, works like above as well, just get information from mpi variables
            mpi_e=mpi_ibool(is)
            mpi_n=mpi_interface(4,is)
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
            end do
            do i=1,NpF
                qtemp_t(i,:) = matmul(VT(is,:,:),(qtemp_n(i,:)-qtemp_t(i,:)))
            end do
            do j=1,9
                flux(((is-1)*NpF+1):(is*Npf),j) = qtemp_t(:,j)
            end do
        end if
    end do

  end subroutine computeExactRiemannSF
  
  subroutine computeExactRiemannSFAnelastic(ie,myrank,flux,q,qi,neighbor,aVT,aVTfree,face, mpi_interface,&
       mpi_connection,mpi_ibool, mpi_ibt, mpi_icon,ibt,ibn)
    ! INPUT: see timeloop for describtion of parameters, nn, st, tt are the normal components
    ! DOES: computes the anelastic Riemann flux of an element [1, eq. 2.130 with modifications for free/absorbing boundaries]
    ! RETURNS: flux: anelastic Riemann flux
    integer, dimension(:) :: neighbor,face
    integer, dimension(:,:) :: ibt,ibn
    real(kind=CUSTOM_REAL), dimension(4*NpF,6) :: flux
    real(kind=CUSTOM_REAL), dimension(:,:) :: q
    real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
    integer, dimension(:,:) :: mpi_interface
    integer, dimension(:,:,:) :: mpi_connection
    integer, dimension(:) :: mpi_ibool
    integer, dimension(:,:) :: mpi_ibt
    integer, dimension(:) :: mpi_icon
    real(kind=CUSTOM_REAL), dimension(Npf,9) :: qtemp_n,qtemp_t
    real(kind=CUSTOM_REAL), dimension(NpF,6) :: atemp_t
    real(kind=CUSTOM_REAL), dimension(6,9) :: aVT,aVTfree
    integer :: is,in,i,j,ie
    integer :: mpi_e, mpi_n
    integer :: myrank
    flux=0.0
    do is=1,4                                   ! works like the elastic part
        in = neighbor(is)
        if (in>0) then
            qtemp_t=q(ibt(:,is),:)
            qtemp_n=q(ibn(:,is),:)
            do i=1,NpF
                atemp_t(i,:) = matmul(aVT,(qtemp_n(i,:)-qtemp_t(i,:)))
            end do
            do j=1,6
             flux(((is-1)*NpF+1):(is*NpF),j) = atemp_t(:,j)
            end do
        ! boundaries
        ! "------------------------------------------------------------------"
        else if ((in==0).and.(face(is) == -1)) then ! free
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                atemp_t(i,:) = matmul(aVTfree,qtemp_t(i,:))
            end do
            do j=1,6
                flux(((is-1)*NpF+1):(is*NpF),j) = atemp_t(:,j)
            end do
        else if ((in==0).and.face(is) == -2) then !absorb
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                atemp_t(i,:) = -matmul(aVT,qtemp_t(i,:))
            end do
            do j=1,6
                flux(((is-1)*NpF+1):(is*NpF),j) = atemp_t(:,j)
            end do
        else if (in == -1) then ! mpi interface
            mpi_e=mpi_ibool(is)
            mpi_n=mpi_interface(4,is)
            qtemp_t=q(ibt(:,is),:)
            do i=1,NpF
                qtemp_n(i,:)=qi(mpi_ibt(i,is),:,mpi_e,mpi_n)
            end do
            do i=1,NpF
                atemp_t(i,:) = matmul(aVT,(qtemp_n(i,:)-qtemp_t(i,:)))
            end do
            do j=1,6
                flux(((is-1)*NpF+1):(is*Npf),j) = atemp_t(:,j)
            end do
        end if
    end do

  end subroutine computeExactRiemannSFAnelastic
end module waveMod
