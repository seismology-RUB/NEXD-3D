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
module allocateMod
    use constantsMod
    use typeMod
    implicit none
!
contains

     subroutine allocateForwardVar(mesh,par,elasticfluxvar,APA,T,invT,VT,VTfree,ux,uy,uz,ax,ay,az,uplot,rQ,Q,Qn,Qm,&
                                    pmlcheck,all_energy,energy_loc,e,&
                                    stf,srcelemV)

        type(meshVar) :: mesh
        type(parameterVar) :: par
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: T,invT,VT,VTfree
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: APA
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: elasticfluxvar,rQ,Q,Qn,Qm
        real(kind=SIZE_DOUBLE), dimension(:,:), allocatable :: e
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uy,uz,ax,ay,az,uplot,stf
        real(kind=SIZE_DOUBLE), dimension(:), allocatable :: all_energy,energy_loc
        integer, dimension(:), allocatable :: srcelemV
        logical, dimension(:), allocatable :: pmlcheck

        allocate(elasticfluxvar(mesh%nelem,4))
        allocate(APA(mesh%nelem,9,9))
        allocate(T(mesh%nelem,4,9,9), invT(mesh%nelem,4,9,9))
        allocate(VT(mesh%nelem,4,9,9), VTfree(mesh%nelem,4,9,9))
        allocate(ux(mesh%nglob),uy(mesh%nglob),uz(mesh%nglob))
        allocate(ax(mesh%nglob),ay(mesh%nglob),az(mesh%nglob))
        allocate(uplot(mesh%nglob))
        allocate(rQ(mesh%nglob,9),Q(mesh%nglob,9),Qn(mesh%nglob,9),Qm(mesh%nglob,9))
        allocate(pmlcheck(mesh%nelem))
        allocate(all_energy(par%nt), energy_loc(par%nt), e(mesh%nglob,6))
        allocate(stf(par%nt))
        allocate(srcelemV(mesh%nelem))

     end subroutine allocateForwardVar

     subroutine deallocateForwardVar(elasticfluxvar,APA,T,invT,VT,VTfree,ux,uy,uz,ax,ay,az,uplot,rQ,Q,Qn,Qm,&
                                    pmlcheck,all_energy,energy_loc,e,&
                                    stf,srcelemV)

        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: T,invT,VT,VTfree
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: APA
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: elasticfluxvar,rQ,Q,Qn,Qm
        real(kind=SIZE_DOUBLE), dimension(:,:), allocatable :: e
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: ux,uy,uz,ax,ay,az,uplot,stf
        real(kind=SIZE_DOUBLE), dimension(:), allocatable :: all_energy,energy_loc
        integer, dimension(:), allocatable :: srcelemV
        logical, dimension(:), allocatable :: pmlcheck

        deallocate(elasticfluxvar)
        deallocate(APA)
        deallocate(T, invT)
        deallocate(VT, VTfree)
        deallocate(ux,uy,uz)
        deallocate(ax,ay,az)
        deallocate(uplot)
        deallocate(rQ,Q,Qn,Qm)
        deallocate(pmlcheck)
        deallocate(all_energy, energy_loc, e)
        deallocate(stf)
        deallocate(srcelemV)

     end subroutine deallocateForwardVar

     subroutine allocateMPIVar(mesh,q_send,q_rec,qi)

        type(meshVar) :: mesh
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send,q_rec

        allocate(q_send(NpF*mesh%mpi_ne*9,mesh%mpi_nn))
        allocate(q_rec(NpF*mesh%mpi_ne*9,mesh%mpi_nn))
        allocate(qi(NpF,9,mesh%mpi_ne,mesh%mpi_nn))

     end subroutine allocateMPIVar

     subroutine deallocateMPIVar(q_send,q_rec,qi)

        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: qi
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: q_send,q_rec

        deallocate(q_send,q_rec,qi)

     end subroutine deallocateMPIVar

     subroutine allocateAttenuationVar(mesh,a_rQ,theta,thetan,thetam,anelasticvar,aVT,aVTfree,aT,aAPA)

        type(meshVar) :: mesh
        real(kind=CUSTOM_REAL),dimension(:,:,:,:), allocatable :: aVT,aVTfree,aT
        real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: theta,thetan,thetam,aAPA
        real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: a_rQ
        real(kind=CUSTOM_REAL),dimension(:), allocatable :: anelasticvar

        allocate(a_rQ(mesh%nglob,6*nMB),theta(mesh%nglob,6,nMB),thetan(mesh%nglob,6,nMB),thetam(mesh%nglob,6,nMB))
        allocate(anelasticvar(mesh%nelem*3*nMB))
        allocate(aVT(mesh%nelem,4,6,9), aVTfree(mesh%nelem,4,6,9),aT(mesh%nelem,4,6,6),aAPA(mesh%nelem,6,9))

     end subroutine allocateAttenuationVar

     subroutine deallocateAttenuationVar(a_rQ,theta,thetan,thetam,anelasticvar,aVT,aVTfree,aT,aAPA)

        real(kind=CUSTOM_REAL),dimension(:,:,:,:), allocatable :: aVT,aVTfree,aT
        real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: theta,thetan,thetam,aAPA
        real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: a_rQ
        real(kind=CUSTOM_REAL),dimension(:), allocatable :: anelasticvar

        deallocate(a_rQ,theta,thetan,thetam)
        deallocate(anelasticvar)
        deallocate(aVT, aVTfree,aT,aAPA)

     end subroutine deallocateAttenuationVar

     subroutine allocatePMLVar(mesh,fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem,pmlloc, &
                                  ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz)

        type(MeshVar) :: mesh
        real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem
        real(kind=CUSTOM_REAL),dimension(:), allocatable :: ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz
        integer, dimension(:,:), allocatable :: pmlloc

        allocate(fprime(mesh%nglob,9), gprime(mesh%nglob,9), hprime(mesh%nglob,9))
        allocate(fprimen(mesh%nglob,9), gprimen(mesh%nglob,9), hprimen(mesh%nglob,9))
        allocate(fprimem(mesh%nglob,9), gprimem(mesh%nglob,9), hprimem(mesh%nglob,9))
        allocate(pmlloc(mesh%nelem,3))
        allocate(ddx(mesh%nglob), ddy(mesh%nglob), ddz(mesh%nglob))
        allocate(alphax(mesh%nglob), alphay(mesh%nglob), alphaz(mesh%nglob))
        allocate(kx(mesh%nglob), ky(mesh%nglob), kz(mesh%nglob))

     end subroutine allocatePMLVar

     subroutine deallocatePMLVar(fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem,pmlloc, &
                                  ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz)

        real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: fprime,gprime,hprime,fprimen,gprimen,hprimen,fprimem,gprimem,hprimem
        real(kind=CUSTOM_REAL),dimension(:), allocatable :: ddx,ddy,ddz,alphax,alphay,alphaz,kx,ky,kz
        integer, dimension(:,:), allocatable :: pmlloc

        deallocate(fprime, gprime, hprime)
        deallocate(fprimen, gprimen, hprimen)
        deallocate(fprimem, gprimem, hprimem)
        deallocate(pmlloc)
        deallocate(ddx, ddy, ddz)
        deallocate(alphax, alphay, alphaz)
        deallocate(kx,ky,kz)

     end subroutine deallocatePMLVar

     subroutine allocateSrcVar(src,par,srcInt,t0,srcArray,srcArrayM,plotstf,Ms)

        type(SrcVar) :: src
        type(ParameterVar) :: par
        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: srcArray,srcArrayM
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf,Ms
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: t0

        allocate(srcInt(Np,src%nsrc))
        allocate(t0(src%nsrc))
        allocate(srcArray(Np,3,src%nsrc,par%nt))
        allocate(srcArrayM(Np,6,src%nsrc,par%nt))
        allocate(plotstf(par%nt,2,src%nsrc))
        allocate(Ms(3,3,src%nsrc))

     end subroutine allocateSrcVar

     subroutine deallocateSrcVar(srcInt,t0,srcArray,srcArrayM,plotstf,Ms)

        real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: srcArray,srcArrayM
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: plotstf,Ms
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: srcInt
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: t0

        deallocate(srcInt)
        deallocate(t0)
        deallocate(srcArray)
        deallocate(srcArrayM)
        deallocate(plotstf)
        deallocate(Ms)

     end subroutine deallocateSrcVar

     subroutine allocateRecVar(rec,par,recInt,recTemp,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz)

        type(RecVar) :: rec
        type(ParameterVar) :: par
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: recTemp

        allocate(recInt(Np,rec%nrec),recTemp(Np))
        allocate(plotux(rec%nrec,par%nt),plotuy(rec%nrec,par%nt),plotuz(rec%nrec,par%nt))
        allocate(plotvx(rec%nrec,par%nt),plotvy(rec%nrec,par%nt),plotvz(rec%nrec,par%nt))
        allocate(plotax(rec%nrec,par%nt),plotay(rec%nrec,par%nt),plotaz(rec%nrec,par%nt))

     end subroutine allocateRecVar

     subroutine deallocateRecVar(recInt,recTemp,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz)

        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: recInt,plotux,plotuy,plotuz,plotvx,plotvy,plotvz,plotax,plotay,plotaz
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: recTemp

        deallocate(recInt,recTemp)
        deallocate(plotux,plotuy,plotuz)
        deallocate(plotvx,plotvy,plotvz)
        deallocate(plotax,plotay,plotaz)

     end subroutine deallocateRecVar

     subroutine deallocMeshvar(this)
        type(meshVar) :: this
        if (allocated(this%matval)) deallocate(this%matval)
        if (allocated(this%vx)) deallocate(this%vx)
        if (allocated(this%vy)) deallocate(this%vy)
        if (allocated(this%vz)) deallocate(this%vz)
        if (allocated(this%rx)) deallocate(this%rx)
        if (allocated(this%ry)) deallocate(this%ry)
        if (allocated(this%rz)) deallocate(this%rz)
        if (allocated(this%sx)) deallocate(this%sx)
        if (allocated(this%sy)) deallocate(this%sy)
        if (allocated(this%sz)) deallocate(this%sz)
        if (allocated(this%tx)) deallocate(this%tx)
        if (allocated(this%ty)) deallocate(this%ty)
        if (allocated(this%tz)) deallocate(this%tz)
        if (allocated(this%nx)) deallocate(this%nx)
        if (allocated(this%ny)) deallocate(this%ny)
        if (allocated(this%nz)) deallocate(this%nz)
        if (allocated(this%sJ)) deallocate(this%sJ)
        if (allocated(this%fscale)) deallocate(this%fscale)
        if (allocated(this%nnormal)) deallocate(this%nnormal)
        if (allocated(this%stangential)) deallocate(this%stangential)
        if (allocated(this%ttangential)) deallocate(this%ttangential)
        if (allocated(this%jacobian)) deallocate(this%jacobian)
        if (allocated(this%mat)) deallocate(this%mat)
        if (allocated(this%ibool)) deallocate(this%ibool)
        if (allocated(this%coord)) deallocate(this%coord)
        if (allocated(this%elem)) deallocate(this%elem)
        if (allocated(this%neighbor)) deallocate(this%neighbor)
        if (allocated(this%face)) deallocate(this%face)
        if (allocated(this%vp)) deallocate(this%vp)
        if (allocated(this%vs)) deallocate(this%vs)
        if (allocated(this%rho)) deallocate(this%rho)
        if (allocated(this%mu)) deallocate(this%mu)
        if (allocated(this%mu)) deallocate(this%vpu)
        if (allocated(this%mu)) deallocate(this%vsu)
        if (allocated(this%mu)) deallocate(this%qp)
        if (allocated(this%mu)) deallocate(this%qs)
        if (allocated(this%mu)) deallocate(this%muu)
        if (allocated(this%mu)) deallocate(this%lambdau)
        if (allocated(this%lambda)) deallocate(this%lambda)
        if (allocated(this%loc2glob_nodes)) deallocate(this%loc2glob_nodes)
        if (allocated(this%mpi_interface)) deallocate(this%mpi_interface)
        if (allocated(this%mpi_neighbor)) deallocate(this%mpi_neighbor)
        if (allocated(this%mpi_connection)) deallocate(this%mpi_connection)
        if (allocated(this%mpi_ninterface)) deallocate(this%mpi_ninterface)
        if (allocated(this%mpi_ibool)) deallocate(this%mpi_ibool)
        if (allocated(this%mpi_ibt)) deallocate(this%mpi_ibt)
        if (allocated(this%mpi_icon)) deallocate(this%mpi_icon)
        if (allocated(this%pml)) deallocate(this%pml)
        if (allocated(this%vol)) deallocate(this%vol)
    end subroutine deallocMeshvar

    subroutine deallocSrcVar(this)
      type(srcVar) :: this
      if (associated(this%srcxyz)) deallocate(this%srcxyz)
      if (associated(this%srcrst)) deallocate(this%srcrst)
      if (associated(this%srcelem)) deallocate(this%srcelem)
      if (associated(this%srci)) deallocate(this%srci)
      if (associated(this%srctype)) deallocate(this%srctype)
      if (associated(this%srcstf)) deallocate(this%srcstf)
      if (associated(this%srcf0)) deallocate(this%srcf0)
      if (associated(this%srcfactor)) deallocate(this%srcfactor)
      if (associated(this%srcangle_force)) deallocate(this%srcangle_force)
      if (associated(this%srcm)) deallocate(this%srcm)
    end subroutine deallocSrcVar

    subroutine deallocRecVar(this)
      type(recVar) :: this
      if (associated(this%recxyz)) deallocate(this%recxyz)
      if (associated(this%recrst)) deallocate(this%recrst)
      if (associated(this%recelem)) deallocate(this%recelem)
      if (associated(this%reci)) deallocate(this%reci)
      if (associated(this%recnr)) deallocate(this%recnr)
    end subroutine deallocRecVar

end module allocateMod
