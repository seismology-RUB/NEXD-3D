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
module timeloopMod
    ! module to calculate the time loop for the DG3G program with a RK TVD time integration
    ! written by LL LL
    use constantsMod
    use meshMod
    use parameterMod
    use sourceReceiverMod
    use mpiMod
    use vandermondeMod
    use stfMod
    use vtkMod
    use plotMod
    use waveMod
    use logoMod
    use pmlMod
    use rungekuttaMod
    use typeMod
    use allocateMod
    implicit none
!
contains
    !
    subroutine timeloop3D(par, myrank)
        ! INPUT: paramter container par, number of rank
        ! DOES: performs calculation of NDG method in 3D for every rank
        ! RETURNS: Nothing
        ! REFERENCE PAPERS/THESES:
            ! 1: Lambrecht, L: Forward and inverse modelling of seismic waves for reconaissance in mechanized tunneling, Dissertation, 2015
            ! 2: Komatitsch, Martin: An unsplit convolutional perfectly matched layer improved at grazing incidence for the seismic wave equation, 2007
            ! 3: YiFeng Li: Development of numerical simualtion method for nonlinear elastodynamic: application to acoustic imaging of defect with the help
            !               of cavity chaotic transducer, Dissertation, 2011
            ! 4: Kaeser et al.: An arbitrary high-order Discontinous Galerkin method for elastic waves on unstructured meshes -
            !                   III: Viscoelastic attenuation, 2007
        type(parameterVar) :: par                                                   ! container for parameters from parfile
        type(meshVar) :: mesh                                                       ! container for mesh variables
        type(srcVar) ::src                                                          ! source variables
        type(recVar) ::rec                                                          ! receiver variables
        ! time variables
        real(kind=CUSTOM_REAL) :: deltat                                 ! timestep with cfl, timestep without cfl
        real(kind=CUSTOM_REAL), dimension(:),allocatable :: t0  ! t0 for source time function, amplitude of stf
        real(kind=CUSTOM_REAL) :: f0, f0temp                                     ! global highest center frequency of sources, local highest center frequency in every rank
        ! timer
        integer :: t1,t2,count_max,rate                                             ! real time variables for system_clock function
        real(kind=CUSTOM_REAL) :: dt1,dtges
        integer :: n_output                                                         ! timesteps per ouput
        ! mesh
         integer(kind=CUSTOM_REAL) :: elem_sum                                                                         ! number of global elements                                                     ! point numbers of inner points
        ! TEST
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: elasticfluxvar
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: anelasticvar
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: T,invT, aT, VT, VTfree, aVT, aVTfree
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: APA,aAPA                               ! matrices: A + |A|, A - |A|
        !PML
        integer, dimension(:,:), allocatable :: pmlloc                                              ! position of PML (x,y or z PML)
        logical, dimension(:), allocatable :: pmlcheck                                              ! shows if an element belongs to PML
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ddx,ddy, ddz, alphax, alphay, alphaz, kx, ky, kz   ! Parameters for PML
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: fprime, gprime, hprime, fprimen, gprimen, hprimen, fprimem, gprimem, hprimem ! PML fluxes
        real(kind=CUSTOM_REAL) :: avg_energy1, avg_energy2, sta_lta                ! Trigger parameters if PML gets instable
        integer :: timecrit, avg_window1, avg_window2, pmlswitch
        logical :: pmlcrit=.true.                                                                   ! to switch off trigger after it turned off PML
        !Energy
        real(kind=SIZE_DOUBLE), dimension(:), allocatable :: all_energy, energy_loc   ! global and local total, kinetic, potential energy
        real(kind=CUSTOM_REAL), dimension(Np,Np) ::invmass, vdmTinv                    ! van-der-monde matrix, its inverse, inverse mass matrix, inverse transponse vdm
        !output
        character(len=256) ::filename                            ! filename for input or output
        ! MPI
        integer :: myrank                                                   ! rank number
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send       ! variables from q to send to other ranks
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_rec        ! variables from q to be received from other ranks
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi       ! ??? AL
        ! src
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt   !??? AL
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf    ! 1: time axis of source time function, 2: source time function
        real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: srcArray    ! to turn stf parallel to srcangle_force
        real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: srcArrayM   ! to use moment tensor as source
        integer, dimension(:), allocatable :: srcelemv                          ! vector containing source numbering
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: Ms                   ! moment tensor variables
        ! rec
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt   ! ??? AL
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: recTemp
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: plotux,plotuy,plotuz ! plot variables for receiver
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: plotvx,plotvy,plotvz
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: plotax,plotay,plotaz
        real(kind=CUSTOM_REAL), dimension(1) :: r_v,s_v,t_v                 ! Temp variable for receiver coordinates
        ! fields
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uy,uz,ax,ay,az,uplot    ! displacement, acceleration, norm of displacement
        real(kind=CUSTOM_REAL) :: maxu                                                  ! maximum norm of displacement
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q,rq,qn,qm               ! vector of stresses and velocities, vector of their time derivaties, dummy vectors to save values
        real(kind=CUSTOM_REAL), dimension(Np) :: dFdr,dFds,dFdt                         ! differentation matrices times fluxes
        real(kind=CUSTOM_REAL), dimension(Np) :: dGdr,dGds,dGdt
        real(kind=CUSTOM_REAL), dimension(Np) :: dHdr,dHds,dHdt
        real(kind=CUSTOM_REAL), dimension(Npf*4,9) :: flux                              ! flux
        real(kind=CUSTOM_REAL) , dimension(:), allocatable :: stf                      ! time in simulation at specific timestep
        real(kind=CUSTOM_REAL) , dimension(9,9) :: free                                 ! matrix for calculation of fluxes at free surface
        ! attenuation
        real(kind=CUSTOM_REAL), dimension(Np,6*nMB) :: a_ftemp,a_gtemp,a_htemp          ! attenuation fluxes
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: theta,thetam,thetan    ! attenuation coefficients
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: a_rq                     ! extra term of rq for attenuation
        real(kind=CUSTOM_REAL), dimension(Npf*4,6) :: aflux                             ! extra flux term for attenuation
        ! tempvars and counter
        real(kind=8), dimension(:,:), allocatable :: e          ! energy dummy
        real(kind=8) :: temp1, temp2
        real(kind=CUSTOM_REAL), dimension(Np,9) :: ftemp,gtemp,htemp              ! fluxes
        integer, dimension(Np) :: iv                            ! numbers of interpolation points in global numbering
        integer :: i,j,r                                        ! counters
        integer :: iglob
        integer :: it,ie,is,c
        integer :: irk, nrk, wrk
        ! usefulls
        real(kind=CUSTOM_REAL) :: onehalf = 1./2.               ! to improve calculation speed
        real(kind=CUSTOM_REAL) :: onethree = 1./3.
        real(kind=CUSTOM_REAL) :: twothree = 2./3.
        real(kind=CUSTOM_REAL) :: onefour = 1./4.
        real(kind=CUSTOM_REAL) :: threefour = 3./4.

        call writeEmptyLineToScreen(myrank,par%nproc)
! load database from meshVar files
        call readMesh(myrank,mesh)
        call writeEmptyLineToScreen(myrank,par%nproc)
! load sources
        if (mesh%has_src) then
            call readSrc(myrank,src)
            write(*,*) "found source in rank ",myrank+1, "and element",src%srcelem
            f0temp=maxval(src%srcf0)
        else
            f0temp=0.0
        end if

        if (par%nproc > 1) then
            call maxval_real_all(f0temp,f0,CUSTOM_REAL)                                 ! get maximum central frequency of all sources
        else
            f0=f0temp
        endif

        call writeEmptyLineToScreen(myrank,par%nproc)

        if (mesh%has_rec) then
            call readRec(myrank,rec)
            write(*,*) "found receiver in rank ",myrank+1, "and element",rec%recelem
        end if
        if (par%nproc > 1) call sync_mpi()

        call allocateForwardVar(mesh,par,elasticfluxvar,APA,T,invT,VT,VTfree,ux,uy,uz,ax,ay,az,uplot,rQ,Q,Qn,Qm,&
                                    pmlcheck,all_energy,energy_loc,e,&
                                    stf,srcelemV)
        if (par%nproc > 1) call allocateMPIVar(mesh,q_send,q_rec,qi)
        if (par%attenuation) call allocateAttenuationVar(mesh,a_rQ,theta,thetan,thetam,anelasticvar,aVT,aVTfree,aT,aAPA)
        if (par%set_pml) call allocatePMLVar(mesh,fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem,pmlloc, &
                                                ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz)
        if (mesh%has_src) call allocateSrcVar(src,par,srcInt,t0,srcArray,srcArrayM,plotstf,Ms)
        if (mesh%has_rec) call allocateRecVar(rec,par,recInt,recTemp,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz)

        ! prepare variables
        if (mesh%has_rec) then
            plotux=0.0
            plotuy=0.0
            plotuz=0.0
            plotvx=0.0
            plotvy=0.0
            plotvz=0.0
            plotax=0.0
            plotay=0.0
            plotaz=0.0
        end if
        q=1e-24                         ! set start values for simulation parameters
        rQ=1e-24                        ! set stresses and velocities (and others) to very small values but not zero
        if (par%attenuation) a_rQ=1e-24
        if (par%attenuation)theta=1e-24
        Qn=1e-24
        ux=1e-24
        uy=1e-24
        uz=1e-24
        ax=1e-24
        ay=1e-24
        az=1e-24
        uplot=1e-24
        dtges = 0.0
        energy_loc=0.0
        all_energy=0.0
        deltat=mesh%dtfactor*par%cfl
        flux=1e-24
        if (par%attenuation) aflux=1e-24
        pmlcheck=.false.
        timecrit=2.4/f0/deltat
        n_output=10
        pmlswitch=0             ! needs to be assigned even when PML is switched off, would lead to segmentation fault while check at end of time loop

        free = 0.0                      ! define matrix for free surfaces
        free(1,1)= -2
        free(2,2)= 0
        free(3,3)= 0
        free(4,4)= -2
        free(5,5)= 0
        free(6,6)= -2
        free(7,7)=  0
        free(8,8)=  0
        free(9,9)=  0

        ftemp=0.0
        gtemp=0.0
        htemp=0.0

        if (par%attenuation) then
            a_ftemp=0.0
            a_gtemp=0.0
            a_htemp=0.0
        endif

        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log.and.myrank==0) write(*,*) "dt =",deltat
        if (par%timeint == 1) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) "use euler"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 1 !euler
            wrk = 1
        else if (par%timeint == 2) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) "use tvd runge kutta second order"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 2 ! rk2
            wrk = 1
        else if (par%timeint == 3) then
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            if (par%log.and.myrank==0) write(*,*) "use tvd runge kutta third order"
            if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
            nrk = 3 ! rk3
            wrk = 2
        end if

        if (par%set_pml) then
            call PMLPreCalculation(mesh,par,f0,ddx,ddy,ddz,kx,ky,kz,alphax,alphay,alphaz,pmlloc,pmlcheck,myrank)
            fprime=0.0
            gprime=0.0
            hprime=0.0
        endif

        do ie=1,mesh%nelem
            elasticfluxvar(ie,1)=-(mesh%lambdau(ie)+2*mesh%muu(ie))
            elasticfluxvar(ie,2)=-mesh%lambdau(ie)
            elasticfluxvar(ie,3)=-mesh%muu(ie)
            elasticfluxvar(ie,4)=-1.0/mesh%rho(ie)
            APA(ie,:,:)=getAPA(mesh%vpu(ie),mesh%vsu(ie),mesh%rho(ie),mesh%lambdau(ie),mesh%muu(ie))
            if (par%attenuation) aAPA(ie,:,:) = getAnelasticAPA(mesh%vpu(ie),mesh%vsu(ie),mesh%rho(ie),mesh%lambdau(ie),mesh%muu(ie))
            do is=1,4
                T(ie,is,:,:)      =getT(mesh%nnormal(:,is,ie),mesh%stangential(:,is,ie),mesh%ttangential(:,is,ie))
                invT(ie,is,:,:)   =getinvT(mesh%nnormal(:,is,ie),mesh%stangential(:,is,ie),mesh%ttangential(:,is,ie))
                if (par%attenuation) aT(ie,is,:,:)     = T(ie,is,1:6,1:6)
                VT(ie,is,:,:)     = matmul(T(ie,is,:,:),matmul(APA(ie,:,:),invT(ie,is,:,:)))
                VTfree(ie,is,:,:) = matmul(T(ie,is,:,:),matmul(matmul(APA(ie,:,:),free),invT(ie,is,:,:)))
                if (par%attenuation) aVT(ie,is,:,:)    = matmul(aT(ie,is,:,:),matmul(aAPA(ie,:,:),invT(ie,is,:,:)))
                if (par%attenuation) aVTfree(ie,is,:,:)= matmul(aT(ie,is,:,:),matmul(matmul(aAPA(ie,:,:),free),invT(ie,is,:,:)))
            enddo
        enddo
        if (par%attenuation) then
            do ie=1,mesh%nelem
                do i=1,nMB
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+1) = -mesh%lambdau(ie)*mesh%ylambda(i,ie)
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+2) = -2*mesh%muu(ie)*mesh%ymu(i,ie)
                    anelasticvar(3*nMB*(ie-1)+3*(i-1)+3) = anelasticvar(3*nMB*(ie-1)+3*(i-1)+1)+anelasticvar(3*nMB*(ie-1)+3*(i-1)+2)
                enddo
            enddo
        endif

        vdmTinv=mesh%vdm                 ! create inverse transpose van-der-monde matrix for source and receiver interpolation
        vdmTinv=transpose(vdmTinv)
        call invert(vdmTinv)

        srcelemv=0
        call prepareSources(mesh,src,par,srcInt,vdmTinv,t0,srcelemv,Ms,plotstf,deltat)
        if (par%nproc > 1) call sync_mpi()

        if (mesh%has_rec) then
            do i=1,rec%nrec                                 ! do for receivers the same than for sources
                r_v(1) = rec%recrst(1,i)
                s_v(1) = rec%recrst(2,i)
                t_v(1) = rec%recrst(3,i)
                call vdm3D(recTemp(:),r_v,s_v,t_v)
                recInt(:,i)=matmul(vdmTinv,recTemp(:))
            end do
        end if

        invmass=matmul(mesh%vdm,transpose(mesh%vdm))                      ! precalculate inverse mass matrix

        if (par%nproc > 1) then
            call sum_int(int(mesh%nelem,CUSTOM_REAL), elem_sum,CUSTOM_REAL)           ! sumup all elements from all ranks
        else
            elem_sum=mesh%nelem
        endif

        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log.and.myrank==0) write(*,*) "Start timeloop over" ,elem_sum," elements with dt: ", deltat," and nt : ",par%nt
        if (par%log.and.myrank==0) write(*,*) "--------------------------------------------------------------------------------------"
        ! ---------------------------------------------------------------------------------------------
        ! ------------------------------------timeloop-------------------------------------------------
        ! ---------------------------------------------------------------------------------------------



        do it=1,par%nt
            stf(it)=it*deltat                       ! get time in simulation
        enddo

        if (par%log.and.myrank==0) then
            call  system_clock(t1, rate, count_max) !get real start time of simulation to estimate total simulation time
        end if

        do it=1,par%nt
            if (mesh%has_src) then
                do i=1,src%nsrc
                    if (src%srctype(i)==0) then
                        ! take one element as acting point for point source
                        srcArray(:,1,i,it)=srcInt(:,i)*(src%srcangle_force(1,i)*plotstf(it,2,i))/mesh%rho(src%srcelem(i))   ! get direction and strength of source
                        srcArray(:,2,i,it)=srcInt(:,i)*(src%srcangle_force(2,i)*plotstf(it,2,i))/mesh%rho(src%srcelem(i))
                        srcArray(:,3,i,it)=srcInt(:,i)*(src%srcangle_force(3,i)*plotstf(it,2,i))/mesh%rho(src%srcelem(i))
                        iv=mesh%ibool(:,src%srcelem(i))
                        srcArray(:,1,i,it)=matmul(invmass,srcArray(:,1,i,it))/mesh%jacobian(iv)                        ! get source term for calculation
                        srcArray(:,2,i,it)=matmul(invmass,srcArray(:,2,i,it))/mesh%jacobian(iv)
                        srcArray(:,3,i,it)=matmul(invmass,srcArray(:,3,i,it))/mesh%jacobian(iv)
                    else if (src%srctype(i)==1) then
                        srcArrayM(:,1,i,it)=srcInt(:,i)*(ms(1,1,i)*plotstf(it,2,i))                              ! moment tensor calculation similar to point source
                        srcArrayM(:,2,i,it)=srcInt(:,i)*(ms(2,2,i)*plotstf(it,2,i))
                        srcArrayM(:,3,i,it)=srcInt(:,i)*(ms(3,3,i)*plotstf(it,2,i))
                        srcArrayM(:,4,i,it)=srcInt(:,i)*(ms(1,2,i)*plotstf(it,2,i))
                        srcArrayM(:,5,i,it)=srcInt(:,i)*(ms(2,3,i)*plotstf(it,2,i))
                        srcArrayM(:,6,i,it)=srcInt(:,i)*(ms(1,3,i)*plotstf(it,2,i))
                        iv=mesh%ibool(:,src%srcelem(i))
                        srcArrayM(:,1,i,it)=matmul(invmass,srcArrayM(:,1,i,it))/mesh%jacobian(iv)
                        srcArrayM(:,2,i,it)=matmul(invmass,srcArrayM(:,2,i,it))/mesh%jacobian(iv)
                        srcArrayM(:,3,i,it)=matmul(invmass,srcArrayM(:,3,i,it))/mesh%jacobian(iv)
                        srcArrayM(:,4,i,it)=matmul(invmass,srcArrayM(:,4,i,it))/mesh%jacobian(iv)
                        srcArrayM(:,5,i,it)=matmul(invmass,srcArrayM(:,5,i,it))/mesh%jacobian(iv)
                        srcArrayM(:,6,i,it)=matmul(invmass,srcArrayM(:,6,i,it))/mesh%jacobian(iv)
                    end if
                end do
            end if
        enddo
        do it=1,par%nt   ! start timeloop
            if (mesh%has_rec) then      ! interpolate receivers
                do r=1,rec%nrec
                    do j=1,Np
                        iglob=mesh%ibool(j,rec%recelem(r))
                        plotux(r,it) = plotux(r,it) + ux(iglob)*recint(j,r)             ! save in every timestep displacement, velocity and acceleration to plot seismograms
                        plotuy(r,it) = plotuy(r,it) + uy(iglob)*recint(j,r)             ! recint gives interpolation of receiver by every point in an element
                        plotuz(r,it) = plotuz(r,it) + uz(iglob)*recint(j,r)
                        plotvx(r,it)  = plotvx(r,it) + q(iglob,7)*recint(j,r)
                        plotvy(r,it)  = plotvy(r,it) + q(iglob,8)*recint(j,r)
                        plotvz(r,it)  = plotvz(r,it) + q(iglob,9)*recint(j,r)
                        plotax(r,it) = plotax(r,it) + ax(iglob)*recint(j,r)
                        plotay(r,it) = plotay(r,it) + ay(iglob)*recint(j,r)
                        plotaz(r,it) = plotaz(r,it) + az(iglob)*recint(j,r)
                    end do
                end do
            end if

            do irk=1,nrk ! start runge-kutta loop
                qm=q
                if (irk==1) qn=q
                if (par%attenuation) then
                    thetam=theta
                    if (irk==1) thetan=theta
                end if
                if (par%set_pml) then
                    fprimem=fprime
                    gprimem=gprime
                    hprimem=hprime
                    if (irk==1) then
                        fprimen=fprime
                        gprimen=gprime
                        hprimen=hprime
                    endif
                endif

                call SendRecVar(mesh,qm,qi,q_send,q_rec)      ! MPI send and receive. qsend and qrec just to use already allocated fields

                do ie=1,mesh%nelem                       ! start element loop
                    iv=mesh%ibool(:,ie)                ! get start point numbers of interpolation points in element _ie_
                    call elasticFluxes(qm(iv,:),elasticfluxvar(ie,:),ftemp,gtemp,htemp)   ! get elastic fluxes

                    if (pmlcheck(ie)) then          ! for PML elements modify flux terms [1, eq: 2.68]
                        do i=1,9
                            ftemp(:,i)=(fprimem(iv,i)+ftemp(:,i))/kx(iv)
                            gtemp(:,i)=(gprimem(iv,i)+gtemp(:,i))/ky(iv)
                            htemp(:,i)=(hprimem(iv,i)+htemp(:,i))/kz(iv)
                        enddo
                    endif

                    do i=1,9                        ! comp weak deriverate
                        dFdr = matmul(mesh%Dr,ftemp(:,i))    ! Dr=M*D_zeta*M^-T from [1, eq. 2.129] -> dFdr = A*Dr*u (same equation)
                        dFds = matmul(mesh%Ds,ftemp(:,i))
                        dFdt = matmul(mesh%Dt,ftemp(:,i))
                        dGdr = matmul(mesh%Dr,gtemp(:,i))
                        dGds = matmul(mesh%Ds,gtemp(:,i))
                        dGdt = matmul(mesh%Dt,gtemp(:,i))
                        dHdr = matmul(mesh%Dr,htemp(:,i))
                        dHds = matmul(mesh%Ds,htemp(:,i))
                        dHdt = matmul(mesh%Dt,htemp(:,i))
                        rq(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds + mesh%tx(iv)*dFdt) + &        ! add d(zeta)/d(x) etc. [1, eq. 2.219] and get total acceleration term without numeric fluxes
                                   (mesh%ry(iv)*dGdr + mesh%sy(iv)*dGds + mesh%ty(iv)*dGdt) + &
                                   (mesh%rz(iv)*dHdr + mesh%sz(iv)*dHds + mesh%tz(iv)*dHdt)
                    end do

                    call computeExactRiemannSF(ie,myrank,flux,qm,qi,mesh%neighbor(:,ie),VT(ie,:,:,:), VTfree(ie,:,:,:)&  ! compute Riemann fluxes on the surface
                        ,mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_connection,mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie), mesh%mpi_icon &
                        ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie))

                    if (par%attenuation) then
                        call anelasticFluxes(qm(iv,:),mesh%wl,a_ftemp,a_gtemp,a_htemp)    ! get fluxes for anelastic case (A*q, B*q)
                        do i=1,6*nMB
                            dFdr = matmul(mesh%Dr,a_ftemp(:,i))                              ! since anelastic matrices are just added to elastic A [1, eq. 2.87 for 2D]
                            dFds = matmul(mesh%Ds,a_ftemp(:,i))                              ! they also have to be incorporated in matrices of eq. 2.219
                            dFdt = matmul(mesh%Dt,a_ftemp(:,i))
                            dGdr = matmul(mesh%Dr,a_gtemp(:,i))
                            dGds = matmul(mesh%Ds,a_gtemp(:,i))
                            dGdt = matmul(mesh%Dt,a_gtemp(:,i))
                            dHdr = matmul(mesh%Dr,a_htemp(:,i))
                            dHds = matmul(mesh%Ds,a_htemp(:,i))
                            dHdt = matmul(mesh%Dt,a_htemp(:,i))
                            a_rq(iv,i) = (mesh%rx(iv)*dFdr + mesh%sx(iv)*dFds + mesh%tx(iv)*dFdt) + &  ! get anelastic part of total acceleration
                                         (mesh%ry(iv)*dGdr + mesh%sy(iv)*dGds + mesh%ty(iv)*dGdt) + &
                                         (mesh%rz(iv)*dHdr + mesh%sz(iv)*dHds + mesh%tz(iv)*dHdt)
                        end do

                        call computeExactRiemannSFAnelastic(ie,myrank,aflux,qm,qi,mesh%neighbor(:,ie),aVT(ie,:,:,:), aVTfree(ie,:,:,:)&  ! compute anelastic part of Riemann fluxes
                            ,mesh%face(:,ie),mesh%mpi_interface(:,:,ie),mesh%mpi_connection,mesh%mpi_ibool(:,ie), mesh%mpi_ibt(:,:,ie), mesh%mpi_icon &
                            ,mesh%ibt(:,:,ie),mesh%ibn(:,:,ie))
                        c=1
                        do i=1,nMB
                            do j=1,6
                                a_rQ(iv,c) = a_rQ(iv,c) + matmul(mesh%lift,mesh%fscale(:,ie)*mesh%wl(i)*aflux(:,j)/2.0) - mesh%wl(i) * thetam(iv,j,i) ! added equations for anelasticity [1, eq. 2.84 6*nMb last equations (2D case) ]
                                c=c+1
                            end do
                        end do
                        do j=1,nMB  ! first equations of [1, eq. 2.84 2D case]
                            rq(iv,1) = rq(iv,1) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+3) * thetam(iv,1,j) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+1) * (thetam(iv,2,j) + thetam(iv,3,j))  ! first 3 equations of eq. 2.84
                            rq(iv,2) = rq(iv,2) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+3) * thetam(iv,2,j) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+1) * (thetam(iv,3,j) + thetam(iv,1,j))
                            rq(iv,3) = rq(iv,3) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+3) * thetam(iv,3,j) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+1) * (thetam(iv,1,j) + thetam(iv,2,j))
                            rq(iv,4) = rq(iv,4) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+2) * thetam(iv,4,j)
                            rq(iv,5) = rq(iv,5) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+2) * thetam(iv,5,j)
                            rq(iv,6) = rq(iv,6) + anelasticvar(3*nMB*(ie-1)+3*(j-1)+2) * thetam(iv,6,j)
                        end do
                    end if

                    do i=1,9        ! lift surface integral
                        rq(iv,i) = rq(iv,i) + matmul(mesh%lift,mesh%fscale(:,ie)*flux(:,i)/2.0)
                    end do

                    do i=1,src%nsrc
                        if (src%srcelem(i)==ie) then      ! add source term
                            if (src%srctype(i) == 0) then
                                rq(iv,7:9)=rq(iv,7:9)+srcArray(:,1:3,i,it)
                            elseif (src%srctype(i) == 1) then
                                rq(iv,1:6)=rq(iv,1:6)+srcArrayM(:,1:6,i,it)
                            endif
                        endif
                    enddo
                enddo


                call RungeKuttaSteps(par,wrk,irk,mesh%nelem,mesh%ibool,q,qn,rq,deltat,theta,thetan,a_rq,fprime,gprime,hprime,fprimen,gprimen,hprimen,&
                    alphax,alphay,alphaz,ddx,ddy,ddz,kx,ky,kz,ftemp,gtemp,htemp,pmlcheck,onehalf,threefour,onefour,twothree,onethree)

                if (irk==nrk) then
                    ux(:) = ux(:) + deltat*q(:,7)        ! Euler integration to obtain displacement
                    uy(:) = uy(:) + deltat*q(:,8)
                    uz(:) = uz(:) + deltat*q(:,9)
                    ax(:) = rq(:,7)                       ! Save acceleration
                    ay(:) = rq(:,8)
                    az(:) = rq(:,9)
                    uplot(:) = sqrt(ux(:)**2+uy(:)**2+uz(:)**2)     ! norm of displacement
                    if (par%set_pml) then
                        do ie=1,mesh%nelem
                            iv=mesh%ibool(:,ie)
                            temp1=1./(4*mesh%mu(ie)**2+6*mesh%mu(ie)*mesh%lambda(ie))
                            temp2=2*mesh%mu(ie)
                            e(iv,1)=(2*(mesh%lambda(ie)+mesh%mu(ie))*q(iv,1)-mesh%lambda(ie)*(q(iv,2)+q(iv,3)))/temp1      ! Strain for energy calculation
                            e(iv,2)=(2*(mesh%lambda(ie)+mesh%mu(ie))*q(iv,2)-mesh%lambda(ie)*(q(iv,1)+q(iv,3)))/temp1
                            e(iv,3)=(2*(mesh%lambda(ie)+mesh%mu(ie))*q(iv,3)-mesh%lambda(ie)*(q(iv,1)+q(iv,1)))/temp1
                            e(iv,4)=q(iv,4)/temp2
                            e(iv,5)=q(iv,5)/temp2
                            e(iv,6)=q(iv,6)/temp2
                        enddo
                    endif
                endif
                if (par%nproc > 1) call sync_mpi()
            end do                      ! runge-kutta loop


            if(par%set_pml) then                                ! STA-LTA trigger for unstable PML (not necessary for MPML???) AL
                do ie=1,mesh%nelem               !Energy loop
                    iv=mesh%ibool(:,ie)
                        energy_loc(it) = energy_loc(it) + 0.5 * mesh%jacobian(iv(1)) * &
                                ( mesh%rho(ie) * dot_product(q(iv,7),matmul(mesh%mass,q(iv,7))) + mesh%rho(ie) * dot_product(q(iv,8),matmul(mesh%mass,q(iv,8))) + mesh%rho(ie) * dot_product(q(iv,9),matmul(mesh%mass,q(iv,9))) &
                                + dot_product(q(iv,1),matmul(mesh%mass,e(iv,1))) + dot_product(q(iv,2),matmul(mesh%mass,e(iv,2))) + dot_product(q(iv,3),matmul(mesh%mass,e(iv,3))) &
                                + 2*(dot_product(q(iv,4),matmul(mesh%mass,e(iv,4))) + dot_product(q(iv,5),matmul(mesh%mass,e(iv,5))) + dot_product(q(iv,6),matmul(mesh%mass,e(iv,6)))) )
                enddo
                call sum_real_all_DP(energy_loc(it),all_energy(it))    ! energy from all ranks
                if (pmlcrit) then
                    if (it>timecrit) then                   ! actual time bigger than source acting time
                        avg_energy1=0.0
                        do i=1,par%avg_window1                  ! long term average
                            if (it-i>0) then
                                avg_energy1=avg_energy1+all_energy(it-i)/avg_window1
                            endif
                        enddo
                        avg_energy2=0.0
                        do i=1,par%avg_window2                  ! short term average
                            if (it-i>0) then
                                avg_energy2=avg_energy2+all_energy(it-i)/avg_window2
                            endif
                        enddo
                        sta_lta=avg_energy2/avg_energy1     ! get trigger
                        if (abs(sta_lta-1)>par%sta_lta_trigger) then    ! trigger larger than treshold?
                            write(*,*) ''
                            write(*,*) '!! Energy grows too much, pml switched off !!'
                            write(*,*) ''
                            pmlcrit=.false.                 ! don't check again for trigger
                            pmlcheck(:)=.false.             ! switch off PML
                            pmlswitch=it                    ! save timestep at which switched off
                        endif
                    endif
                endif
            endif

            if (mod(it,n_output)==0 .or. it==par%nt) then
                if (myrank == 0 ) then                      ! print simulation information
                    call  system_clock(t2, rate, count_max) ! get time
                    dt1=real(t2 - t1, kind=8)/real(rate, kind=8)        ! get time difference from last check
                    dtges=dtges + dt1                       ! save total simulation run time
                end if
                if (par%nproc > 1) then
                    call maxval_real(maxval(uplot),maxu,CUSTOM_REAL)        ! get total norm
                else
                    maxu=maxval(uplot)
                endif
                if (par%log .and. (myrank == 0)) then           ! print information
                    write(*,*) "--------------------------------------------------------------------------------------"
                    if (par%set_pml) then
                        write(*,*) " timestep :", it, "norm u :", maxu, "Energy: ", all_energy(it)
                    else
                        write(*,*) " timestep :", it, "norm u :", maxu
                    endif
                    write(*,*) " time in simulation: ", stf(it), " seconds"
                    write(*,*) " time for ",n_output," timesteps: ", dt1, " seconds"
                    write(*,*) " simulation will ready in approx ", (par%nt-it)*(dt1/n_output)/60/60, "hours =", (par%nt-it)*(dt1/n_output)/60, "minutes"
                    write(*,*) " simulation runs for ", dtges, " seconds"
                    write(*,*) "--------------------------------------------------------------------------------------"
                    if (par%set_pml) then
                        if (it==par%nt.and.pmlswitch>0) then
                            write(*,*) ''
                            write(*,*) 'PML was switched off at time step ', pmlswitch
                        endif
                    endif
                end if
            end if

            if ((mod(it,par%frame)==0).and.(par%movie)) then
                write(filename,"('/moviedata',i6.6,'_it',i7.7,'.bin')") myrank+1,it     ! write files for vtk movies
                filename=trim(outpath)//trim(filename)
                open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
                write(27) uplot
                close(27)
            endif

            if (mod(it,n_output)==0) then
                if (myrank == 0) then
                    call  system_clock(t1, rate, count_max)                             !reset timer
                end if
            end if
        end do !timeloop

        if (mesh%has_rec) then                                                   ! print seismograms
            do r=1,rec%nrec
                write(filename,"('/seismo.x.',i7.7,'.sdu')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotux(r,:),filename)
                write(filename,"('/seismo.y.',i7.7,'.sdu')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotuy(r,:),filename)
                write(filename,"('/seismo.z.',i7.7,'.sdu')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotuz(r,:),filename)
                write(filename,"('/seismo.x.',i7.7,'.sdv')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotvx(r,:),filename)
                write(filename,"('/seismo.y.',i7.7,'.sdv')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotvy(r,:),filename)
                write(filename,"('/seismo.z.',i7.7,'.sdv')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotvz(r,:),filename)
                write(filename,"('/seismo.x.',i7.7,'.sda')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotax(r,:),filename)
                write(filename,"('/seismo.y.',i7.7,'.sda')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotay(r,:),filename)
                write(filename,"('/seismo.z.',i7.7,'.sda')") rec%recnr(r)
                filename=trim(outpath)//trim(filename)
                call plotPoints2D(stf,plotaz(r,:),filename)

            ! write in SEM name convention
            !!$          write(filename,"('/X',i2.2,'.DB.FXX.semd')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotux(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXY.semd')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotuy(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXZ.semd')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotuz(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXX.semv')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotvx(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXY.semv')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotvy(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXZ.semv')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotvz(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXX.sema')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotax(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXY.sema')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotay(r,:),filename)
            !!$          write(filename,"('/X',i2.2,'.DB.FXZ.sema')") recnr(r)
            !!$          filename=trim(outpath)//trim(filename)
            !!$          call plotPoints2D(stf(:),plotaz(r,:),filename)
            end do
        endif

        if (myrank==0) then
            write(filename,"('/energy')")
            filename=trim(outpath)//trim(filename)
            call plotPoints2D(stf,real(all_energy,CUSTOM_REAL),filename)
        endif

        call deallocateForwardVar(elasticfluxvar,APA,T,invT,VT,VTfree,ux,uy,uz,ax,ay,az,uplot,rQ,Q,Qn,Qm,&
                                    pmlcheck,all_energy,energy_loc,e,&
                                    stf,srcelemV)
        if (par%nproc > 1) call deallocateMPIVar(q_send,q_rec,qi)
        if (par%attenuation) call deallocateAttenuationVar(a_rQ,theta,thetan,thetam,anelasticvar,aVT,aVTfree,aT,aAPA)
        if (par%set_pml) call deallocatePMLVar(fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem,pmlloc, &
                                  ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz)
        if (mesh%has_src) then
            call deallocateSrcVar(srcInt,t0,srcArray,srcArrayM,plotstf,Ms)
            call deallocSrcVar(src)
        endif
        if (mesh%has_rec) then
            call deallocateRecVar(recInt,recTemp,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz)
            call deallocRecVar(rec)
        endif
        call deallocMeshVar(mesh)
    end subroutine timeloop3D






     subroutine writeEmptyLineToScreen(myrank,nproc)
        integer :: myrank
        integer :: nproc

        if (nproc > 1) call sync_mpi()
        if (myrank==0) write(*,*)
        if (nproc > 1) call sync_mpi()

     end subroutine writeEmptyLineToScreen


 end module timeloopMod
