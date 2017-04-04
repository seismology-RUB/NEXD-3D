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
module pmlMod
  	use constantsMod
  	use typeMod
implicit none

contains

  subroutine dampingProfile(ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz,vx,vy,vz,vp,xmin,xmax,ymin,ymax,zmin,zmax,pml_delta,rc,amax,kmax,iv,pmlloc,i,myrank)
    ! INPUT: vx,vy,vz: coordinates of interpolation points, vp: p-wave velocity, xmin - zmax: coordinate boundaries of setup, pml_delta: thickness of PML,
    !        rc: reflection coefficient, amax, kmax: stretching factor variables, iv: interpolation points of element, pmlloc: flags to which PML (x,y,z)
    !        element belongs, i: element number, myrank: rank number
    ! DOES:  calculates damping profiles and stretching factor variables for Perfectly Matched Layers using eq. 6 from
    !        Martin, Komatitisch et al: A High-Order Time and Space Formulation of the Unsplit Perfectly Matched Layer for the Seismic Wave Equation Using
    !        Auxiliary Differential Equations (ADE-PML), 2010
    ! RETURNS: ddx,ddy,ddz: damping profiles, alphax - kz: stretching factor variables
    real(kind=CUSTOM_REAL), dimension(:) :: vx,vz,vy
    real(kind=CUSTOM_REAL), dimension(:) :: ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz
    real(kind=CUSTOM_REAL), dimension(Np) :: ar,br
    real(kind=CUSTOM_REAL) :: pml_delta,vp,rc,amax,kmax, xmin, xmax, ymin, ymax, zmin, zmax
    integer :: i,myrank
    integer, dimension(3) :: pmlloc
    integer, dimension(NP) :: iv

    if ( pmlloc(1) == -1) then !xmin (does the same for every if depending on specfic coordinate)
       ar=sqrt((xmin+pml_delta)**2)     ! get position where PML begings
       br=sqrt(vx(iv)**2)               ! get position if interpolation point
       br=abs(ar-br)                    ! get x distance to inner edge of PML
       ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2   ! equations from Martin, Komatitisch et al, 2010
       alphax(iv) = amax*(1-br/pml_delta)
       kx(iv) = 1+kmax*(br/pml_delta)**2
    else if ( pmlloc(1) == 1 ) then !xmax
       ar=sqrt((xmax-pml_delta)**2)
       br=sqrt(vx(iv)**2)
       br=abs(ar-br)
       ddx(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2
       alphax(iv) = amax*(1-br/pml_delta)
       kx(iv) = 1+kmax*(br/pml_delta)**2
    endif
    if ( pmlloc(2) == -1 ) then !ymin
       ar=sqrt((ymin+pml_delta)**2)
       br=sqrt(vy(iv)**2)
       br=abs(ar-br)
       ddy(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2
       alphay(iv) = amax*(1-br/pml_delta)
       ky(iv) = 1+kmax*(br/pml_delta)**2
    else if ( pmlloc(2) == 1 ) then !ymax
       ar=sqrt((ymax-pml_delta)**2)
       br=sqrt(vy(iv)**2)
       br=abs(ar-br)
       ddy(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2
       alphay(iv) = amax*(1-br/pml_delta)
       ky(iv) = 1+kmax*(br/pml_delta)**2
    endif
    if ( pmlloc(3) == -1 ) then !zmin
       ar=sqrt((zmin+pml_delta)**2)
       br=sqrt(vz(iv)**2)
       br=abs(ar-br)
       ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2
       alphaz(iv) = amax*(1-br/pml_delta)
       kz(iv) = 1+kmax*(br/pml_delta)**2
    else if ( pmlloc(3) == 1 ) then !zmax
       ar=sqrt((zmax-pml_delta)**2)
       br=sqrt(vz(iv)**2)
       br=abs(ar-br)
       ddz(iv) = -3.0 / (2*pml_delta) * vp * alog(real(rc))*( br/pml_delta )**2
       alphaz(iv) = amax*(1-br/pml_delta)
       kz(iv) = 1+kmax*(br/pml_delta)**2
    endif
  end subroutine dampingProfile

     subroutine PMLPreCalculation(mesh,par,f0,ddx,ddy,ddz,kx,ky,kz,alphax,alphay,alphaz,pmlloc,pmlcheck,myrank)

         type(meshVar) :: mesh
         type(ParameterVar) :: par
         real(kind=CUSTOM_REAL) :: f0
         real(kind=CUSTOM_REAL), dimension(:) :: ddx,ddy,ddz,kx,ky,kz,alphax,alphay,alphaz
         integer, dimension(:,:) :: pmlloc
         logical, dimension(:) :: pmlcheck
         integer :: myrank

         real(kind=CUSTOM_REAL) :: amax,xmean,ymean,zmean, pxy, pxz, pyx, pyz, pzx, pzy, ddx_temp, ddy_temp, ddz_temp
         integer, dimension(Np) :: iv
         integer :: i,j


         if (par%log.and.myrank==0) write(*,*) "use PML boundary conditions"
         if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
         amax=par%pml_afac*pi*f0             ! get maximum for alpha [2, p. 160]
         ddx=0.0                         ! predefine damping profiles and variables of PML stretching factor
         ddy=0.0
         ddz=0.0
         kx=1.0
         ky=1.0
         kz=1.0
         alphax=0.0
         alphay=0.0
         alphaz=0.0
         pmlloc=0

         do i=1, mesh%nelem
             iv=mesh%ibool(:,i)
             if (mesh%pml(i)>0) then
                 pmlcheck(i)=.true.      ! flag element as PML element
                 xmean=0.
                 ymean=0.
                 zmean=0.
                 do j=1,NP               ! get center of element to determine its position
                     xmean=xmean+mesh%vx(iv(j))/Np
                     ymean=ymean+mesh%vy(iv(j))/Np
                     zmean=zmean+mesh%vz(iv(j))/Np
                 enddo
                 if (xmean<mesh%pmlxmin+par%pml_delta) pmlloc(i,1)=-1     ! define to which PML (x,y,z) an element belongs (possibly multiple)
                 if (xmean>mesh%pmlxmax-par%pml_delta) pmlloc(i,1)=1
                 if (ymean<mesh%pmlymin+par%pml_delta) pmlloc(i,2)=-1
                 if (ymean>mesh%pmlymax-par%pml_delta) pmlloc(i,2)=1
                 if (zmean<mesh%pmlzmin+par%pml_delta) pmlloc(i,3)=-1
                 if (zmean>mesh%pmlzmax+par%pml_delta) pmlloc(i,3)=1
                 call dampingProfile(ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz,mesh%vx,mesh%vy,mesh%vz,mesh%vp(i),mesh%pmlxmin,mesh%pmlxmax,mesh%pmlymin,mesh%pmlymax,mesh%pmlzmin,mesh%pmlzmax,par%pml_delta,par%pml_rc,amax,par%pml_kmax,iv,pmlloc(i,:),i,myrank)    ! get damping profiles
             endif
         enddo
         pyx=0.3     ! here we use Multiaxial PML (i.e. the damping profile of x depends also on y and z)
         pzx=0.3     ! These factors describe the influence of y and z to damping profile of x
         pxy=0.3     ! See [3, p. 130]
         pzy=0.3
         pxz=0.3
         pyz=0.3
         do i=1,mesh%nglob
             if (mesh%vx(i)>mesh%pmlxmin+par%pml_delta.and.mesh%vx(i)<mesh%pmlxmax-par%pml_delta.and.mesh%vy(i)>mesh%pmlymin+par%pml_delta.and.mesh%vy(i)<mesh%pmlymax-par%pml_delta.and.mesh%vz(i)>mesh%pmlzmin+par%pml_delta) then
                 ddx(i)=0.0          ! set all points outside of PML to non-PML values (maybe outdated)
                 ddy(i)=0.0
                 ddz(i)=0.0
                 alphax(i)=0.0
                 alphay(i)=0.0
                 alphaz(i)=0.0
                 kx(i)=1.0
                 ky(i)=1.0
                 kz(i)=1.0
             endif
             ddx_temp=ddx(i)
             ddy_temp=ddy(i)
             ddz_temp=ddz(i)
             ddx(i)=ddx_temp+pyx*ddy_temp+pzx*ddz_temp       ! apply p factors
             ddy(i)=pxy*ddx_temp+ddy_temp+pzy*ddz_temp
             ddz(i)=pxz*ddx_temp+pyz*ddy_temp+ddz_temp
         enddo
     end subroutine PMLPreCalculation
end module pmlMod
