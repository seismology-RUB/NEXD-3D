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
module rungekuttaMod
    use constantsMod
    use typeMod
implicit none

contains

     subroutine RungeKuttaSteps(par,wrk,irk,nelem,ibool,q,qn,rq,deltat,theta,thetan,a_rq,fprime,gprime,hprime,fprimen,gprimen,hprimen,&
                            alphax,alphay,alphaz,ddx,ddy,ddz,kx,ky,kz,ftemp,gtemp,htemp,pmlcheck,onehalf,threefour,onefour,twothree,onethree)

         type(ParameterVar) :: par
         integer, dimension(:,:) :: ibool
         integer, dimension(Np) :: iv
         integer :: wrk, irk,ie,c,i,j,nelem
         real(kind=CUSTOM_REAL), dimension(:,:,:) :: theta, thetan
         real(kind=CUSTOM_REAL), dimension(:,:) :: q, qn, rq, a_rq,fprime,gprime, hprime,fprimen,gprimen,hprimen
         real(kind=CUSTOM_REAL), dimension(:) :: alphax,alphay,alphaz,ddx,ddy,ddz,kx,ky,kz
         real(kind=CUSTOM_REAL), dimension(Np,9) :: ftemp,gtemp,htemp
         real(kind=CUSTOM_REAL) :: deltat, onehalf,threefour,onefour,twothree,onethree
         logical, dimension(:) :: pmlcheck

         if (wrk==1) then            ! RK 2 and Euler
             q(:,:) = RK2(q(:,:),qn(:,:),rq(:,:),deltat,irk,onehalf)
             if(par%attenuation) then
                 c=1
                 do i=1,nMB
                     do j=1,6
                         theta(:,j,i) = RK2(theta(:,j,i),thetan(:,j,i),a_rq(:,c),deltat,irk,onehalf)
                         c=c+1
                     end do
                 end do
             end if
             do ie=1,nelem
                 if (pmlcheck(ie)) then                                              ! Euler step for PML fluxes [1, eq. 2.67]
                     iv=ibool(:,ie)
                     do i=1,9
                         do j=1,Np
                             fprime(iv(j),i) = RK2(fprime(iv(j),i),fprimen(iv(j),i),-alphax(iv(j))*fprime(iv(j),i) - ddx(iv(j))/kx(iv(j))*(fprime(iv(j),i)+ftemp(j,i)),deltat,irk,onehalf)
                             gprime(iv(j),i) = RK2(gprime(iv(j),i),gprimen(iv(j),i),-alphay(iv(j))*gprime(iv(j),i) - ddy(iv(j))/ky(iv(j))*(gprime(iv(j),i)+gtemp(j,i)),deltat,irk,onehalf)
                             hprime(iv(j),i) = RK2(hprime(iv(j),i),hprimen(iv(j),i),-alphaz(iv(j))*hprime(iv(j),i) - ddz(iv(j))/kz(iv(j))*(hprime(iv(j),i)+htemp(j,i)),deltat,irk,onehalf)
                         enddo
                     enddo
                 endif
             enddo
         elseif (wrk==2) then            ! RK 2 and Euler
             q(:,:) = RK3(q(:,:),qn(:,:),rq(:,:),deltat,irk,threefour,onefour,twothree,onethree)
             if(par%attenuation) then
                 c=1
                 do i=1,nMB
                     do j=1,6
                         theta(:,j,i) = RK3(theta(:,j,i),thetan(:,j,i),a_rq(:,c),deltat,irk,threefour,onefour,twothree,onethree)
                         c=c+1
                     end do
                 end do
             end if
             do ie=1,nelem
                 if (pmlcheck(ie)) then                                              ! Euler step for PML fluxes [1, eq. 2.67]
                     iv=ibool(:,ie)
                     do i=1,9
                         do j=1,Np
                             fprime(iv(j),i) = RK3(fprime(iv(j),i),fprimen(iv(j),i),-alphax(iv(j))*fprime(iv(j),i) - ddx(iv(j))/kx(iv(j))*(fprime(iv(j),i)+ftemp(j,i)),deltat,irk,threefour,onefour,twothree,onethree)
                             gprime(iv(j),i) = RK3(gprime(iv(j),i),gprimen(iv(j),i),-alphay(iv(j))*gprime(iv(j),i) - ddy(iv(j))/ky(iv(j))*(gprime(iv(j),i)+gtemp(j,i)),deltat,irk,threefour,onefour,twothree,onethree)
                             hprime(iv(j),i) = RK3(hprime(iv(j),i),hprimen(iv(j),i),-alphaz(iv(j))*hprime(iv(j),i) - ddz(iv(j))/kz(iv(j))*(hprime(iv(j),i)+htemp(j,i)),deltat,irk,threefour,onefour,twothree,onethree)
                         enddo
                     enddo
                 endif
             enddo
            endif
     end subroutine RungeKuttaSteps

     elemental function RK2(a,an,dadt,dt,irk,onehalf)

        real(kind=CUSTOM_REAL),intent(in) :: a,an, dadt
        real(kind=CUSTOM_REAL),intent(in) :: dt,onehalf
        real(kind=CUSTOM_REAL) :: RK2
        integer, intent(in) :: irk

        if (irk==1) then
            RK2=a+dadt*dt
        else
            RK2=onehalf*(an+a+dadt*dt)
        endif

     end function RK2

     elemental function RK3(a,an,dadt,dt,irk,threefour,onefour,twothree,onethree)

        real(kind=CUSTOM_REAL),intent(in) :: a,an, dadt
        real(kind=CUSTOM_REAL),intent(in) :: dt,threefour,onefour,twothree,onethree
        real(kind=CUSTOM_REAL) :: RK3
        integer, intent(in) :: irk

        if (irk==1) then
            RK3=a+dadt*dt
        elseif (irk==2) then
            RK3=threefour*an+onefour*(a+dadt*dt)
        else
            RK3=onethree*an+twothree*(a+dadt*dt)
        endif

     end function RK3


end module rungekuttaMod
