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
module sourceReceiverMod
	use constantsMod
	use parameterMod
	use nodesMod
	use tetTrafoMod
	use vandermondeMod
	use matrixMod
	use geometricFactorsMod
	use allocateMod
	use typeMod
	use stfMod
	use mpiMod
  implicit none

contains

  subroutine initSource(par,this)
    ! INPUT: general parameter container par
    ! DOES: reads in source parameters
    ! OUTPUT: source parameter container this
    type (srcVar) :: this
    type (parameterVar) :: par
    character(len=80) :: filename
    integer :: ier, isrc, pos
    real(kind=CUSTOM_REAL):: abs_vector

    filename=trim('data/source')
    open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
    if (par%log) then
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*, "(a26, a14, a12, a28)")&
           "|                          ","Begin reading ", filename, "...                        |"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
    endif

    call readIntPar(this%nsrc, "nsrc", filename, 0)
    if (this%nsrc==0) then
        write(*,*) 'CAUTION!! There is no source in your simulation!'
    else if (this%nsrc<0) then
        write(*,*) 'Number of source is negative. ABORT!'
        stop
    endif
    allocate(this%srcxyz(3,this%nsrc))
    allocate(this%srcrst(3,this%nsrc))
    allocate(this%srcelem(this%nsrc))
    allocate(this%srci(this%nsrc))
    allocate(this%srctype(this%nsrc))
    allocate(this%srcstf(this%nsrc))
    allocate(this%srcf0(this%nsrc))
    allocate(this%srcfactor(this%nsrc))
    allocate(this%srcangle_force(3,this%nsrc))
    allocate(this%srcM(6,this%nsrc))

    ! Cycle to read the parameters form the source file. The order of appearence of the parameters is not important
    do isrc=1,this%nsrc
        pos = setFilePosition("source", filename, isrc) + 11     ! define position in source file to search for parameters
        call readFloatPar(this%srcxyz(1, isrc), "xsource", filename, pos)    ! read in source parameters
        call readFloatPar(this%srcxyz(2, isrc), "ysource", filename, pos)
        call readFloatPar(this%srcxyz(3, isrc), "zsource", filename, pos)
        call readIntPar(this%srctype(isrc), "sourcetype", filename, pos)
        call readIntPar(this%srcstf(isrc), "stf", filename, pos)
        call readFloatPar(this%srcf0(isrc), "f0", filename, pos)
        if (this%srcf0(isrc)<=0.) then
            write(*,*) 'Source frequency of source ', isrc, ' is zero or negative. ABORT!'
            stop
        endif
        call readFloatPar(this%srcfactor(isrc), "factor", filename, pos)
        call readFloatPar(this%srcangle_force(1,isrc), "xangle_force", filename, pos)
        call readFloatPar(this%srcangle_force(2,isrc), "yangle_force", filename, pos)
        call readFloatPar(this%srcangle_force(3,isrc), "zangle_force", filename, pos)
        abs_vector=sqrt(this%srcangle_force(1,isrc)**2+this%srcangle_force(2,isrc)**2+this%srcangle_force(3,isrc)**2)
        this%srcangle_force(:,isrc)=this%srcangle_force(:,isrc)/abs_vector
        call readFloatPar(this%srcM(1, isrc), "Mxx", filename, pos)
        call readFloatPar(this%srcM(2, isrc), "Myy", filename, pos)
        call readFloatPar(this%srcM(3, isrc), "Mzz", filename, pos)
        call readFloatPar(this%srcM(4, isrc), "Mxy", filename, pos)
        call readFloatPar(this%srcM(5, isrc), "Mxz", filename, pos)
        call readFloatPar(this%srcM(6, isrc), "Myz", filename, pos)
    end do
    close(19)

    !print found values
    if (par%log) then
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*,"(a27, a13, a13, a27)") &
         "|                          ","Done reading ", filename, "                          |"     
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*,"(a40, i10, a30)")   "|                    Number of sources: ", this%nsrc, "                             |"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        do isrc = 1, this%nsrc
             write (*,"(a40, f10.2, a30)")&
               "|                X-value of the source: ", this%srcxyz(1, isrc), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|                Y-value of the source: ", this%srcxyz(2, isrc), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|                Z-value of the source: ", this%srcxyz(3, isrc), "                             |"
            write (*,"(a40, i10, a30)")&
               "|                       Type of source: ", this%srctype(isrc), "                             |" 
            write (*,"(a40, i10, a30)")&
               "|           source-time-function (stf): ", this%srcstf(isrc), "                             |"
            write (*,"(a40, f10.1, a30)")&
               "|              Center frequency of stf: ", this%srcf0(isrc), "                            |"
            write (*,"(a40, f10.1, a30)")&
               "|                    Factor of the stf: ", this%srcfactor(isrc), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|          X-component of force action: ", this%srcangle_force(1,isrc), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|          Y-component of force action: ", this%srcangle_force(2,isrc), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|          Z-component of force action: ", this%srcangle_force(3,isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mxx: ", this%srcM(1, isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Myy: ", this%srcM(2, isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mzz: ", this%srcM(3, isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mxy: ", this%srcM(4, isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Mxz: ", this%srcM(5, isrc), "                             |"
            write (*,"(a40, f10.3, a30)")&
               "|                     Momenttensor Myz: ", this%srcM(6, isrc), "                             |"
            write (*, "(a80)") "--------------------------------------------------------------------------------"
        enddo
    endif

  end subroutine initSource
!
  subroutine initReceiver(par,this)
    ! INPUT: general parameter container par
    ! DOES: reads in receiver parameters
    ! OUTPUT: receiver parameter container this
    type(parameterVar) :: par
    type(recVar) :: this
    integer :: irec
    character(len=256) :: filename

    filename=trim('data/stations')
    open(unit=19,file=trim(filename))
    if (par%log) then
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*, "(a26, a14, a12, a28)")&
         "|                          ","Begin reading ", filename, "...                        |"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
    endif

    read(19,*) this%nrec                ! get number of receivers

    allocate(this%recxyz(3,this%nrec))
    allocate(this%recrst(3,this%nrec))
    allocate(this%recelem(this%nrec))
    allocate(this%reci(this%nrec))
    allocate(this%recnr(this%nrec))
    do irec=1,this%nrec                 ! for all receivers: read their numbers and positions
        read(19,*) this%recnr(irec),this%recxyz(1,irec),this%recxyz(2,irec),this%recxyz(3,irec)
    end do
    close(19)
    this%recrst=0
    this%recelem=0
    this%reci=0
    if (par%log) then
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*,"(a27, a13, a13, a27)")&
         "|                          ","Done reading ", filename, "                          |"     
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        write (*, "(a40, i10, a30)")   "|                  Number of receivers: ", this%nrec, "                             |"
        write (*, "(a80)") "--------------------------------------------------------------------------------"
        do irec=1,this%nrec
            write (*,"(a40, i10, a30)") &
               "|                             receiver: ", this%recnr(irec), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|              X-value of the receiver: ", this%recxyz(1, irec), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|              Y-value of the receiver: ", this%recxyz(2, irec), "                             |"
            write (*,"(a40, f10.2, a30)")&
               "|              Z-value of the receiver: ", this%recxyz(3, irec), "                             |"
            write (*, "(a80)") "--------------------------------------------------------------------------------"
       enddo
    endif
  end subroutine initReceiver
!
!
  subroutine findSource(par,this,vx,vy,vz,nglob,nelem,ibool,coord,elem,dr,ds,dt)
    type(parameterVar) :: par
    type(srcVar) :: this
    integer, dimension(:,:), intent(in) :: elem
    real(kind=CUSTOM_REAL), dimension(:,:) :: coord
    integer, dimension(:,:), intent(in) :: ibool
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vx,vy,vz
    integer :: nglob, nelem
    real(kind=CUSTOM_REAL), dimension(Np) :: x,y,z,r,s,t
    logical, dimension(nelem,this%nsrc) :: checkelem
    integer :: isrc, i, ie, j, iglob, besave
    integer, dimension(NP) :: iv
    real(kind=CUSTOM_REAL) :: epsilon
    real(kind=CUSTOM_REAL) :: src_r, src_s, src_t
    real(kind=CUSTOM_REAL), dimension(1) :: src_rr, src_ss, src_tt
    real(kind=CUSTOM_REAL), dimension(3) :: xyztemp
    integer :: num_iter=4
    real(kind=CUSTOM_REAL), dimension(3) :: xyz
    real(kind=CUSTOM_REAL), dimension(Np) :: rx,ry,rz,sx,sy,sz,tx,ty,tz, jacobian
    real(kind=CUSTOM_REAL) :: dx,dy,dz
    real(kind=CUSTOM_REAL), dimension(Np) :: srctemp, srcInt 
    real(kind=CUSTOM_REAL), dimension(Np,NP) :: vdm, VdmTinv
    real(kind=CUSTOM_REAL), dimension(:,:) :: dr,ds,dt
    real(kind=CUSTOM_REAL) :: drdx,drdy,drdz,dsdx,dsdy,dsdz,dtdx,dtdy,dtdz
    real(kind=CUSTOM_REAL)  :: ddr,dds,ddt


    checkelem=.true. 
    ! get local element
    call nodes3D(balpha(NGLL),x,y,z)    ! get interpolation points in tetrahedron
    call xyzToRst(x,y,z,r,s,t)          ! coordinate transformation xyz -> rst
    if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
    call vdm3d(vdm,r,s,t)               ! get van-der-monde matrix, transpose, invert
    vdmTinv=vdm
    vdmTinv=transpose(vdmTinv)
    call invert(vdmTinv)

    ! do initial guess for source points
    isrc=1
    besave=1
    do while (isrc <= this%nsrc)
        epsilon=1e7
        do ie=1,nelem
            do i=1,Np
                iglob=ibool(i,ie)
                if ((sqrt((vx(iglob)-this%srcxyz(1,isrc))**2+(vy(iglob)-this%srcxyz(2,isrc))**2+(vz(iglob)-this%srcxyz(3,isrc))**2) < epsilon).and.checkelem(ie,isrc)) then  !check how close every point coordinate is to the searched source point
                    epsilon= sqrt((vx(iglob)-this%srcxyz(1,isrc))**2+(vy(iglob)-this%srcxyz(2,isrc))**2+(vz(iglob)-this%srcxyz(3,isrc))**2)     ! find smallest distance
                    this%srci=iglob ! save closest point as initial guess (final result is found after loop through all elements is finished)
                    src_r=r(i)      ! save coordinates in reference element
                    src_s=s(i)
                    src_t=t(i)
                    this%srcelem(isrc)=ie   ! save source element
                end if
            end do
        end do
        xyztemp(1) = this%srcxyz(1,isrc)    ! initial set of source coordinate are the coordinates from source file
        xyztemp(2) = this%srcxyz(2,isrc)
        xyztemp(3) = this%srcxyz(3,isrc)
        ie=this%srcelem(isrc)               ! save element number of source and all its point numbers
        iv = ibool(:,ie)
        call geometricFactors3d(rx,sx,tx,ry,sy,ty,rz,sz,tz,jacobian,vx(iv),vy(iv),vz(iv),dr,ds,dt)  ! get derivatives of r to x, etc.

        do i=1, num_iter
            xyz(1) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(1,elem(1,ie)) + (1.0+src_r) * coord(1,elem(2,ie)) +&  ! x,y,z get coordinates in element
               (1+src_s) * coord(1,elem(3,ie)) + (1+src_t) * coord(1,elem(4,ie)))
            xyz(2) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(2,elem(1,ie)) + (1.0+src_r) * coord(2,elem(2,ie)) +&
               (1+src_s) * coord(2,elem(3,ie)) + (1+src_t) * coord(2,elem(4,ie)))
            xyz(3) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(3,elem(1,ie)) + (1.0+src_r) * coord(3,elem(2,ie)) +&
               (1+src_s) * coord(3,elem(3,ie)) + (1+src_t) * coord(3,elem(4,ie)))

            src_rr(1)=src_r
            src_ss(1)=src_s
            src_tt(1)=src_t

            call vdm3D(srcTemp(:),src_rr,src_ss,src_tt) ! get vdm for current source guess
            srcInt(:)=matmul(vdmTinv,srcTemp(:))

            drdx=0
            drdy=0
            drdz=0
            dsdx=0
            dsdy=0
            dsdz=0
            dtdx=0
            dtdy=0
            dtdz=0
            do j=1,NP                   ! calculate correction
                drdx=drdx+srcInt(j)*rx(j)
                drdy=drdy+srcInt(j)*ry(j)
                drdz=drdz+srcInt(j)*rz(j)
                dsdx=dsdx+srcInt(j)*sx(j)
                dsdy=dsdy+srcInt(j)*sy(j)
                dsdz=dsdz+srcInt(j)*sz(j)
                dtdx=dtdx+srcInt(j)*tx(j)
                dtdy=dtdy+srcInt(j)*ty(j)
                dtdz=dtdz+srcInt(j)*tz(j)
            end do

            dx = (xyztemp(1)-xyz(1))
            dy = (xyztemp(2)-xyz(2))
            dz = (xyztemp(3)-xyz(3))
            ddr= drdx*dx+drdy*dy+drdz*dz
            dds= dsdx*dx+dsdy*dy+dsdz*dz
            ddt= dtdx*dx+dtdy*dy+dtdz*dz
            src_r = src_r + ddr         ! improve src position
            src_s = src_s + dds
            src_t = src_t + ddt

            if (i==num_iter) then
                if (src_r/=src_r .or. src_s/=src_s .or. src_t/=src_t) then  ! check if a variable is NaN
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else if (src_r < -1.001) then
                    checkelem(ie,isrc) = .false.        ! check if right element was found. If not flag element, restart
                    besave=besave+1
                else if (src_r > 1.001) then
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else if (src_s < -1.001) then
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else if (src_s > 1.001) then
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else if (src_t < -1.001) then
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else if (src_t > -0.999-src_r-src_s) then
                    checkelem(ie,isrc) = .false.
                    besave=besave+1
                else
                    if (src_r > 1.) src_r = 1.0
                    if (src_r < -1.) src_r = -1.0
                    if (src_s > 1.) src_s = 1.0
                    if (src_s < -1.) src_s = -1.0
                    if (src_t < -1.) src_t = -1.0
                    if (src_t > -1.-src_r-src_s) src_t = -1.0-src_r-src_s
                    this%srcrst(1,isrc) = src_r     ! save local coordinates
                    this%srcrst(2,isrc) = src_s
                    this%srcrst(3,isrc) = src_t
                    xyz(1) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(1,elem(1,ie)) + (1.0+src_r) * coord(1,elem(2,ie)) +&  ! x,y,z get coordinates in element
                        (1+src_s) * coord(1,elem(3,ie)) + (1+src_t) * coord(1,elem(4,ie)))
                    xyz(2) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(2,elem(1,ie)) + (1.0+src_r) * coord(2,elem(2,ie)) +&
                        (1+src_s) * coord(2,elem(3,ie)) + (1+src_t) * coord(2,elem(4,ie)))
                    xyz(3) = 0.5 * ( -(1.0+src_r+src_s+src_t) * coord(3,elem(1,ie)) + (1.0+src_r) * coord(3,elem(2,ie)) +&
                        (1+src_s) * coord(3,elem(3,ie)) + (1+src_t) * coord(3,elem(4,ie)))
                    if (par%log) write(*,*) "Found source   ",isrc, "at global coordinates", xyz(1),xyz(2),xyz(3)
                    if (par%log) write(*,*) "Final guess for source ",isrc, "in global element", this%srcelem(isrc), " and local coordinates", this%srcrst(1,isrc),this%srcrst(2,isrc),this%srcrst(3,isrc)
                    isrc=isrc+1                     ! go to next source
                    besave=1                        ! reset
                end if
                if (besave > 20) then !20 guesses, than abbort
                    if(par%log) write(*,*) "ERROR, could not find source", isrc,"at global coordinates", xyz(1),xyz(2),xyz(3)
                    stop
                end if
            end if
        end do
    end do !isrc
  end subroutine findSource
!
  subroutine findReceiver(par,this,vx,vy,vz,nglob,nelem,ibool,coord,elem,dr,ds,dt)
    type(parameterVar) :: par
    type(recVar), intent(inout) :: this
    integer, dimension(:,:), intent(in) :: elem
    real(kind=CUSTOM_REAL), dimension(:,:) :: coord
    integer, dimension(:,:), intent(in) :: ibool
    real(kind=CUSTOM_REAL), dimension(:), intent(in) :: vx,vy,vz
    real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: dr,ds,dt
    integer :: nglob, nelem
    real(kind=CUSTOM_REAL), dimension(Np) :: x,y,z,r,s,t

    logical, dimension(nelem,this%nrec) :: checkelem
    integer :: irec, i, ie, j, iglob, besave
    integer, dimension(NP) :: iv
    real(kind=CUSTOM_REAL) :: epsilon
    real(kind=CUSTOM_REAL) :: rec_r, rec_s, rec_t
    real(kind=CUSTOM_REAL), dimension(1) :: rec_rr, rec_ss, rec_tt
    real(kind=CUSTOM_REAL), dimension(3) :: xyztemp
    integer :: num_iter=4
    real(kind=CUSTOM_REAL), dimension(3) :: xyz
    real(kind=CUSTOM_REAL), dimension(Np) :: rx,ry,rz,sx,sy,sz,tx,ty,tz, jacobian
    real(kind=CUSTOM_REAL) :: dx,dy,dz
    real(kind=CUSTOM_REAL), dimension(Np) :: rectemp, recInt 
    real(kind=CUSTOM_REAL), dimension(Np,NP) :: vdm, VdmTinv
    real(kind=CUSTOM_REAL) :: drdx,drdy,drdz,dsdx,dsdy,dsdz,dtdx,dtdy,dtdz
    real(kind=CUSTOM_REAL)  :: ddr,dds,ddt

    write(*,*) "in find rec"

    checkelem=.true. 
    ! get local element
    call nodes3D(balpha(NGLL),x,y,z)
    call xyzToRst(x,y,z,r,s,t)
    if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
    call vdm3d(vdm,r,s,t)
    vdmTinv=vdm
    vdmTinv=transpose(vdmTinv)
    call invert(vdmTinv)

    ! do initial guess for receiver points  (This subroutine works analog to the source subroutine!!)
    irec=1
    besave=1
    do while (irec <= this%nrec)
        epsilon=1e7
        do ie=1,nelem
            do i=1,Np
                iglob=ibool(i,ie)
                if ((sqrt((vx(iglob)-this%recxyz(1,irec))**2+(vy(iglob)-this%recxyz(2,irec))**2+(vz(iglob)-this%recxyz(3,irec))**2) < epsilon).and.checkelem(ie,irec)) then
                    epsilon= sqrt((vx(iglob)-this%recxyz(1,irec))**2+(vy(iglob)-this%recxyz(2,irec))**2+(vz(iglob)-this%recxyz(3,irec))**2)
                    this%reci=iglob !initial guess
                    rec_r=r(i)
                    rec_s=s(i)
                    rec_t=t(i)
                    this%recelem(irec)=ie
                end if
            end do
        end do
        xyztemp(1) = this%recxyz(1,irec)
        xyztemp(2) = this%recxyz(2,irec)
        xyztemp(3) = this%recxyz(3,irec)
        ie=this%recelem(irec)
        iv = ibool(:,ie)
        call geometricFactors3d(rx,sx,tx,ry,sy,ty,rz,sz,tz,jacobian,vx(iv),vy(iv),vz(iv),dr,ds,dt)
        do i=1, num_iter
            xyz(1) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(1,elem(1,ie)) + (1.0+rec_r) * coord(1,elem(2,ie)) +&
               (1+rec_s) * coord(1,elem(3,ie)) + (1+rec_t) * coord(1,elem(4,ie)))
            xyz(2) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(2,elem(1,ie)) + (1.0+rec_r) * coord(2,elem(2,ie)) +&
               (1+rec_s) * coord(2,elem(3,ie)) + (1+rec_t) * coord(2,elem(4,ie)))
            xyz(3) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(3,elem(1,ie)) + (1.0+rec_r) * coord(3,elem(2,ie)) +&
               (1+rec_s) * coord(3,elem(3,ie)) + (1+rec_t) * coord(3,elem(4,ie)))

            rec_rr(1)=rec_r
            rec_ss(1)=rec_s
            rec_tt(1)=rec_t
            call vdm3D(recTemp(:),rec_rr,rec_ss,rec_tt)
            recInt(:)=matmul(vdmTinv,recTemp(:))

            drdx=0
            drdy=0
            drdz=0
            dsdx=0
            dsdy=0
            dsdz=0
            dtdx=0
            dtdy=0
            dtdz=0
            do j=1,NP
                drdx=drdx+recInt(j)*rx(j)
                drdy=drdy+recInt(j)*ry(j)
                drdz=drdz+recInt(j)*rz(j)
                dsdx=dsdx+recInt(j)*sx(j)
                dsdy=dsdy+recInt(j)*sy(j)
                dsdz=dsdz+recInt(j)*sz(j)
                dtdx=dtdx+recInt(j)*tx(j)
                dtdy=dtdy+recInt(j)*ty(j)
                dtdz=dtdz+recInt(j)*tz(j)
            end do

            dx = (xyztemp(1)-xyz(1))
            dy = (xyztemp(2)-xyz(2))
            dz = (xyztemp(3)-xyz(3))
            ddr= drdx*dx+drdy*dy+drdz*dz
            dds= dsdx*dx+dsdy*dy+dsdz*dz
            ddt= dtdx*dx+dtdy*dy+dtdz*dz
            rec_r = rec_r + ddr
            rec_s = rec_s + dds
            rec_t = rec_t + ddt

            if (i==num_iter) then
                if (rec_r/=rec_r .or. rec_s/=rec_s .or. rec_t/=rec_t) then  ! check if a variable is NaN
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_r < -1.001) then  ! check is found point is in element
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_r > 1.001) then
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_s < -1.001) then
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_s > 1.001) then
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_t < -1.001) then
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else if (rec_t > -.999-rec_r-rec_s) then
                    irec=irec
                    checkelem(ie,irec) = .false.
                    besave=besave+1
                else
                    if (rec_r > 1.) rec_r = 1.0
                    if (rec_r < -1.) rec_r = -1.0
                    if (rec_s > 1.) rec_s = 1.0
                    if (rec_s < -1.) rec_s = -1.0
                    if (rec_t < -1.) rec_t = -1.0
                    if (rec_t > -1.-rec_r-rec_s) rec_t = -1.0-rec_r-rec_s
                    this%recrst(1,irec) = rec_r
                    this%recrst(2,irec) = rec_s
                    this%recrst(3,irec) = rec_t
                    xyz(1) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(1,elem(1,ie)) + (1.0+rec_r) * coord(1,elem(2,ie)) +&
                        (1+rec_s) * coord(1,elem(3,ie)) + (1+rec_t) * coord(1,elem(4,ie)))
                    xyz(2) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(2,elem(1,ie)) + (1.0+rec_r) * coord(2,elem(2,ie)) +&
                        (1+rec_s) * coord(2,elem(3,ie)) + (1+rec_t) * coord(2,elem(4,ie)))
                    xyz(3) = 0.5 * ( -(1.0+rec_r+rec_s+rec_t) * coord(3,elem(1,ie)) + (1.0+rec_r) * coord(3,elem(2,ie)) +&
                        (1+rec_s) * coord(3,elem(3,ie)) + (1+rec_t) * coord(3,elem(4,ie)))
                    if (par%log) write(*,*) "Found receiver   ",irec, "at global coordinates", xyz(1),xyz(2),xyz(3)
                    if (par%log) write(*,*) "Final guess for receiver ",irec, "in global element", this%recelem(irec), " and local coordinates", this%recrst(1,irec),this%recrst(2,irec),this%recrst(3,irec)
                    irec=irec+1
                    besave=1
                end if
                if (besave > 20) then
                    if(par%log) write(*,*) "ERROR, could not find receiver", irec,"at global coordinates", xyz(1),xyz(2),xyz(3)
                    stop 'An error occured'
                end if
            end if
        end do
    end do !irec
  end subroutine findReceiver

! change src element
!
subroutine getlocalSrcElement(this,i,element)
  type(srcVar) :: this
  integer :: element,i
  this%srcelem(i) = element
end subroutine getlocalSrcElement
!"--------------------------------------------------------------------------------------"
! change rec element

subroutine getlocalRecElement(this,i,element)
  type(recVar) :: this
  integer :: element,i
  this%recelem(i) = element
end subroutine getlocalRecElement
!"--------------------------------------------------------------------------------------"
! prepare src recalculate
!
subroutine prepareRecalcSrc(this,nsrc,srcxyz,srctype,srcstf,srcf0,srcfactor,srcangle_force,srcM)
  type(srcVar) :: this
  integer :: nsrc
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcxyz, srcM
  real(kind=CUSTOM_REAL), dimension(:), pointer :: srcf0 ,srcfactor
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcangle_force
  integer, dimension(:), pointer :: srctype, srcstf
  real(kind=CUSTOM_REAL), dimension(3,nsrc) :: tempxyz
  real(kind=CUSTOM_REAL), dimension(6,nsrc) :: tempm
  real(kind=CUSTOM_REAL), dimension(nsrc) :: tempf0,tempfactor
  real(kind=CUSTOM_REAL), dimension(3,nsrc) :: tempangle_force
  integer, dimension(nsrc) :: temptype, tempstf 
  tempxyz = srcxyz
  temptype = srctype
  tempstf = srcstf
  tempf0 = srcf0
  tempfactor = srcfactor
  tempangle_force = srcangle_force
  tempm = srcm

  ! create new pointer to avoid haning pointers
  allocate(this%srcxyz(3,nsrc))
  allocate(this%srcrst(3,nsrc))
  allocate(this%srcelem(nsrc))
  allocate(this%srci(nsrc)) 
  allocate(this%srctype(nsrc))
  allocate(this%srcstf(nsrc))
  allocate(this%srcf0(nsrc))
  allocate(this%srcfactor(nsrc))
  allocate(this%srcangle_force(3,nsrc))
  allocate(this%srcM(6,nsrc))
  this%nsrc = nsrc
  this%srcxyz = tempxyz
  this%srctype = temptype
  this%srcstf = tempstf
  this%srcf0 = tempf0
  this%srcfactor = tempfactor
  this%srcangle_force = tempangle_force
  this%srcM = tempm
end subroutine prepareRecalcSrc
!"--------------------------------------------------------------------------------------"
! prepare rec recalculate
!
subroutine prepareRecalcRec(this,nrec,recxyz,recnr)
  type(recVar) :: this
  integer :: nrec
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recxyz
  integer, dimension(:), pointer :: recnr
  real(kind=CUSTOM_REAL), dimension(3,nrec), target :: tempxyz
  integer, dimension(nrec), target :: tempnr

  tempxyz = recxyz
  tempnr = recnr
! produce new pointer to avoid haning pointers
  allocate(this%recxyz(3,nrec))
  allocate(this%recrst(3,nrec))
  allocate(this%recelem(nrec))
  allocate(this%reci(nrec))
  allocate(this%recnr(nrec))
  this%nrec = nrec
  this%recxyz = tempxyz
  this%recnr = tempnr

end subroutine prepareRecalcRec
!"--------------------------------------------------------------------------------------"
! write srcVar
!
subroutine writeSrcVar(this,filename)
  type(srcVar) :: this
  character(len=*) :: filename
  write(*,*) "write scrVar ", filename
  open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
  write(27) this%nsrc
  write(27) this%srcxyz
  write(27) this%srcrst
  write(27) this%srcelem
  write(27) this%srci
  write(27) this%srctype
  write(27) this%srcstf
  write(27) this%srcf0
  write(27) this%srcfactor
  write(27) this%srcangle_force
  write(27) this%srcM
  close(27)
end subroutine writeSrcVar
!"--------------------------------------------------------------------------------------"
! write recVar
!
subroutine writeRecVar(this,filename)
  type(recVar) :: this
  character(len=*) :: filename
  write(*,*) "write scrVar ", filename
  open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
  write(27) this%nrec
  write(27) this%recxyz
  write(27) this%recrst
  write(27) this%recelem
  write(27) this%reci
  write(27) this%recnr
  close(27)
end subroutine writeRecVar

!"--------------------------------------------------------------------------------------"
! read srcVar
!
subroutine readSrcVar(this,filename)
  type(srcVar) :: this
  character(len=*) :: filename
  write(*,*) "read source variables from ", filename
  open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
  read(27) this%nsrc
  allocate(this%srcxyz(3,this%nsrc))
  allocate(this%srcrst(3,this%nsrc))
  allocate(this%srcelem(this%nsrc))
  allocate(this%srci(this%nsrc))
  allocate(this%srctype(this%nsrc))
  allocate(this%srcstf(this%nsrc))
  allocate(this%srcf0(this%nsrc))
  allocate(this%srcfactor(this%nsrc))
  allocate(this%srcangle_force(3,this%nsrc))
  allocate(this%srcM(6,this%nsrc))
  read(27) this%srcxyz
  read(27) this%srcrst
  read(27) this%srcelem
  read(27) this%srci
  read(27) this%srctype
  read(27) this%srcstf
  read(27) this%srcf0
  read(27) this%srcfactor
  read(27) this%srcangle_force
  read(27) this%srcM
  close(27)
end subroutine readSrcVar

!"--------------------------------------------------------------------------------------"
! read recVar
!
subroutine readRecVar(this,filename)
  type(recVar) :: this
  character(len=*) :: filename
  write(*,*) "read receiver variables from ", filename
  open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
  read(27) this%nrec
  allocate(this%recxyz(3,this%nrec))
  allocate(this%recrst(3,this%nrec))
  allocate(this%recelem(this%nrec))
  allocate(this%reci(this%nrec))
  allocate(this%recnr(this%nrec))
  read(27) this%recxyz
  read(27) this%recrst
  read(27) this%recelem
  read(27) this%reci
  read(27) this%recnr
  close(27)
end subroutine readRecVar


     subroutine readSrc(myrank,src)

        integer:: myrank
        logical:: file_exists
        type(SrcVar) :: src
        character(len=256) ::filename

            write(filename,"('/srcVar',i6.6)") myrank+1
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readSrcVar(src,filename)
            else
                write(*,*) "error in srcVar, files not existing"
                call stop_mpi()
            end if

     end subroutine readSrc

     subroutine readRec(myrank,rec)

        integer:: myrank
        logical:: file_exists
        type(RecVar) :: rec
        character(len=256) ::filename

            write(filename,"('/recVar',i6.6)") myrank+1
            filename=trim(outpath)//trim(filename)
            inquire(file=trim(filename), exist=file_exists)
            if (file_exists) then
                call readRecVar(rec,filename)
            else
                write(*,*) "error in srcVar, files not existing"
                call stop_mpi()
            end if

     end subroutine readRec

     subroutine prepareSources(mesh,src,par,srcInt,vdmTinv,t0,srcelemv,Ms,plotstf,deltat)

        type(MeshVar) :: mesh
        type(SrcVar) :: src
        type(ParameterVar) :: par
        real(kind=CUSTOM_REAL), dimension(:,:,:) :: Ms, plotstf
        real(kind=CUSTOM_REAL), dimension(:,:) :: srcInt
        real(kind=CUSTOM_REAL), dimension(:) :: t0
        real(kind=CUSTOM_REAL), dimension(1) :: r_v, s_v, t_v
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: vdmTinv
        real(kind=CUSTOM_REAL), dimension(Np) :: srcTemp
        real(kind=CUSTOM_REAL), dimension(3,3) :: M
        real(kind=CUSTOM_REAL) :: deltat
        integer, dimension(:) :: srcelemv
        integer :: i, it, j, k
        character(len=256) ::filename

        if (mesh%has_src) then
            do i=1,src%nsrc
                r_v(1) = src%srcrst(1,i)                ! save r,s,t coordinates of source
                s_v(1) = src%srcrst(2,i)
                t_v(1) = src%srcrst(3,i)
                call vdm3D(srcTemp,r_v,s_v,t_v)    ! to transform source
                srcInt(:,i)=matmul(vdmTinv,srcTemp)
                t0(i)=1.2/src%srcf0(i)                  ! Offset time for source
                srcelemv(src%srcelem(i)) = i            ! save source number
                if (src%srctype(i) == 1) then          ! if source is a momenttensor (works similar to single force)
                    M(1,1)=src%srcM(1,i)                    !Mxx
                    M(1,2)=src%srcM(4,i)                    !Mxy
                    M(1,3)=src%srcM(5,i)                    !Mxz
                    M(2,1)=src%srcM(4,i)                    !Myx
                    M(2,2)=src%srcM(2,i)                    !Myy
                    M(2,3)=src%srcM(6,i)                    !Myz
                    M(3,1)=src%srcM(5,i)                    !Mzx
                    M(3,2)=src%srcM(6,i)                    !Mzy
                    M(3,3)=src%srcM(3,i)                    !Mzz
                    do j=1,3
                        do k=1,3                        ! ??? AL
                            ms(j,k,i)=0.5*(M(j,k)+M(k,j))
                        end do
                    end do
                end if
                write(filename,"('/stf',i6.6)") i
                filename=trim(outpath)//trim(filename)
                write(*,*) 'write source time function to: ', filename
                open(unit=28,file=trim(filename),status='unknown')
                if (src%srcstf(i) == 3) then
                    do it=1,par%nt
                       plotstf(it,1,i)=(float(it)-1.)*deltat
                    enddo
                    call stfReadFile(plotstf(:,1,i),plotstf(:,2,i))
                else
                    do it=1,par%nt
                        plotstf(it,1,i)=(float(it)-1.)*deltat   ! time axis of stf
                        if (src%srcstf(i) == 1) then
                            plotstf(it,2,i)=-stfGauss(plotstf(it,1,i),src%srcf0(i),t0(i),src%srcfactor(i))     ! get Gaussian as stf
                        else if (src%srcstf(i) == 2) then
                            plotstf(it,2,i)=-stfRicker(plotstf(it,1,i),src%srcf0(i),t0(i),src%srcfactor(i))    ! get Ricker wavelet as stf
                        end if
                    end do
                endif
                do it=1,par%nt
                    write(28,*) plotstf(it,1,i),plotstf(it,2,i) ! print stf to file
                enddo
                if (src%srctype(i) == 1) then
                    call DerivateSTF(plotstf(:,1,i),plotstf(:,2,i),par%nt)
                endif
                close(28)
            enddo
        endif

     end subroutine prepareSources

     subroutine DerivateSTF(time, stf, nt)
        real(kind=CUSTOM_REAL), dimension(:) :: time, stf
        real(kind=CUSTOM_REAL), dimension(:), allocatable ::  stf_temp
        integer :: nt, i
        real(kind=CUSTOM_REAL) :: a,b

        allocate(stf_temp(nt))

        a=( (stf(1)-stf(2))/(time(1)-time(2)) - (stf(2)-stf(3))/(time(2)-time(3)) ) / (time(1)-time(3))
        b=(stf(1)-stf(2))/(time(1)-time(2)) - a * (time(1)+time(2))
        stf_temp(1) = 2*a*time(1)+b
        a=( (stf(nt-2)-stf(nt-1))/(time(nt-2)-time(nt-1)) - (stf(nt-1)-stf(nt))/(time(nt-1)-time(nt)) ) / (time(nt-2)-time(nt))
        b=(stf(nt-2)-stf(nt-1))/(time(nt-2)-time(nt-1)) - a * (time(nt-2)+time(nt-1))
        stf_temp(nt) = 2*a*time(nt)+b

        do i=2,nt-1
            stf_temp(i)=(stf(i+1)-stf(i-1))/(time(i+1)-time(i-1))
        enddo

        stf=stf_temp

        deallocate(stf_temp)

     end subroutine DerivateSTF
end module sourceReceiverMod
