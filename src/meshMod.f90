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
module meshMod
    ! module to partition the mesh provided by cubit
    use constantsMod
    use parameterMod
    use tetraNeighborMod
    use vtkMod
    use nodesMod
    use tetTrafoMod
    use vandermondeMod
    use dmatricesMod
    use liftMod
    use geometricFactorsMod
    use normalsMod
    use sourceReceiverMod
    use logoMod
    use anelasticMod
    use waveMod
    use typeMod
    use allocateMod
    use mpiMod

    implicit none

contains

    subroutine createMesh(par)

        ! INPUT: container of parameters read from parfile (uses also constants from constans.h)
        ! DOES:  Divides mesh in different partitions for every rank and provides files containing databases for every rank
        ! RETURNS: Nothing

        ! types
        type(meshVar) :: this                   ! global container for mesh variables
        type(parameterVar) :: par               ! container for parameters from parfile
        type(meshVar), dimension(:), allocatable :: db ! container for mesh variables of every rank for MPI run
        type(srcVar) :: src, tempsrc            ! container for source variables
        type(recVar) :: rec, temprec            ! container for receiver variables
        ! characters
        character(len=256) :: filename          ! filename for output to files
        character(len=16) :: screeninfo
        logical :: checked
        ! tempvars, counter and dummies
        integer :: i,j,k,c,d,e,is,ie,in,l,m,ee,iproc     ! integers for loops and other things
        integer :: temp1, temp2
        real(kind=CUSTOM_REAL) :: rtemp1, rtemp2,rtemp3, rtemp4,min1,min2
        integer :: dummy
        integer, dimension(:,:), allocatable :: neighbortemp
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: tempM1,tempM2
        logical, dimension(4) :: check = .false.
        integer :: nb_neigh
        integer , dimension(4) :: nb_temp
        integer :: num_glob, num_part
        integer, dimension(:), allocatable :: tempv
        !MPI and global/local numbering
        integer :: mpi_ti, mpi_ni                                           ! to save point numbers for permutation
        integer, dimension(:),allocatable :: part                           ! contains to which rank an element belongs
        integer, dimension(:), allocatable :: glob2loc_elmnts               ! translates global to local rank numbering of elements
        integer, dimension(:,:), allocatable :: loc2glob_elemnts            ! translates local numbering back to global numbering
        integer, dimension(:), allocatable :: num_loc
        integer, dimension(:), allocatable :: nodes_elmnts, nnodes_elmnts  ! global numbering of nodes, number of elements per node
        integer, dimension(:), allocatable :: glob2loc_nodes_nparts        ! number of a node in overall sum _size_glob2loc_nodes_
        integer, dimension(:), allocatable :: glob2loc_nodes_parts         ! rank where node of overall numbering is located
        integer, dimension(:), allocatable :: glob2loc_nodes,glob2loc_nodes2               ! local numbering of nodes in each rank depending on global overall numbering
        integer :: size_glob2loc_nodes                                     ! counter for overall sum of nodes if nodes belonging to multiple ranks counted multiple
        integer, dimension(:), allocatable :: part_nodes, num_parts        ! shows if a node belongs to a rank, number of nodes in a rank
        integer, dimension(4) :: loc_nodes                                 ! local number of a node in its rank
        integer, dimension(:,:), allocatable :: icom                       ! matrix with entry 1 if two ranks connected, 0 otherwise
        ! parVars
        real(kind=CUSTOM_REAL) :: cfl, v_long_edge                           ! cfl number
        ! surfaces
        integer :: nfree, nabs                                              ! number of free/absorbing surface elements
        integer, dimension(:,:), allocatable :: elem_free, elem_absorb      ! contains element number and node number of free/absorbing surface
        integer*4, dimension(Npf) :: perm                                   ! to permute point numbering on a face
        real(kind=CUSTOM_REAL), dimension(:),allocatable :: epsil                              ! maximum deviation to obtain permutation
        integer, dimension(4) :: side                                       ! to predefine faces of elements
        !metis
        integer, dimension(:), allocatable :: elmnts        ! vector containing nodes of all elements
        integer, dimension(:), allocatable  :: xadj         ! contains start indices for every element in adjncy
        integer, dimension(:), allocatable  :: adjncy       ! contains all neighbor elements sequentially
        integer, dimension(0:4) :: options                  ! necessary to call metis function
        integer :: edgecut
        ! nodes
        real(kind=CUSTOM_REAL), dimension(Np) :: x,y,z                  ! node coordinates in tetrahedron
        real(kind=CUSTOM_REAL), dimension(Np) :: r,s,t                  ! node coordinates in reference tetrahedron
        real(kind=CUSTOM_REAL), dimension(Np) :: rx,sx,tx,ry,sy,ty,rz,sz,tz,jacobian    ! derivatives of r with respect to x, etc.; Jacobian
        real(kind=CUSTOM_REAL) :: minGLL                                ! smallest gll distance
        ! bools
        integer, dimension(Np) :: iv,ivn                                ! save point numbers of an element
        !source
        integer, dimension(:), pointer :: locsrctype => null()              ! source type in rank
        integer, dimension(:), pointer :: locsrcstf => null()               ! stf in rank
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locsrcf0 => null() ! center frequency in rank
        real(kind=CUSTOM_REAL), dimension(:), pointer :: locsrcfactor => null() ! amplitude of stf in rank
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locsrcangle_force => null()  ! angle of force action in rank
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locsrcM => null()    ! moment tensor in rank
        integer :: locnsrc                                                      ! number of sources in rank
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locsrcxyz => null()  ! local source coordinates
        integer, dimension(:,:), allocatable :: srcnr                           ! global source numbering
        integer, dimension(:), allocatable :: locsrcnr                          ! local source numbering
        !receiver
        integer :: locnrec                                                      ! local number of receiver
        real(kind=CUSTOM_REAL), dimension(:,:), pointer :: locrecxyz => null()
        integer, dimension(:,:), allocatable :: recnum
        integer, dimension(:), pointer :: locrecnr => null()                    ! local receiver numbering
        integer, dimension(:), pointer :: locrecnum => null()                   ! local receiver in global numbering
        ! anelastic attenuation
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: wl
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: yl,ym
        ! Volume
        !    real(kind=CUSTOM_REAL) :: teta,tetb,tet,tetap,tetbp,tetcp,fa,fb,fc,tetdel
 
        side(4)=1
        side(3)=2
        side(1)=3
        side(2)=4

        ! get information from Parfile
        if (par%log) call writeLogo()
        if (par%log) write(*,*)
        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "Start mesher"
        if (par%log) write(*,*)

        allocate(db(par%nproc))

        ! get values from mesh files to allocate main arrays
        filename=trim('cubit/coord')
        open(unit=19,file=trim(filename))
        read(19,*) this%ncoord
        close(19)

        filename=trim('cubit/mesh')
        open(unit=19,file=trim(filename))
        read(19,*) this%nelem
        close(19)

        !free surface elements?
        filename=trim('cubit/free')
        open(unit=19,file=trim(filename))
        read(19,*) nfree
        close(19)

        !absorb surface elements?
        filename=trim('cubit/absorb')
        open(unit=19,file=trim(filename))
        read(19,*) nabs
        close(19)

        ! number of global interpolation points
        this%nglob = this%nelem*Np

        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "spatial order of simulation       :",ORDER
        if (par%log) write(*,*) "number of point per element       :",Np
        if (par%log) write(*,*) "number of global coordinates      :",this%ncoord
        if (par%log) write(*,*) "number of global elements         :",this%nelem
        if (par%log) write(*,*) "number of free surface elements   :",nfree
        if (par%log) write(*,*) "number of absorb surface elements :",nabs
        if (par%log) write(*,*) "number of global points           :",this%nglob

        ! allocate main arrays for the mesh
        allocate(this%coord(3,this%ncoord))         ! coordinates of all mesh points (x,y and z)
        allocate(this%elem(4,this%nelem))           ! contains the 4 mesh points an element is build of
        allocate(this%neighbor(4,this%nelem))       ! contains neighbor elements at every element face
        allocate(neighbortemp(4,this%nelem))        ! dummy for neighbors
        allocate(this%mat(this%nelem))              ! contains material number of every element
        ! physical parameters of every element (velocities, density, attenuation, Lamé parameters, etc.)
        allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%qp(this%nelem),this%qs(this%nelem),this%mu(this%nelem),this%lambda(this%nelem))
        allocate(this%vpu(this%nelem),this%vsu(this%nelem),this%muu(this%nelem),this%lambdau(this%nelem))
        allocate(this%ylambda(nMB,this%nelem),this%ymu(nMB,this%nelem), this%wl(nMB), this%wk(2*nMB-1))
        allocate(this%pml(this%nelem))              ! contains PML flag for every element
        allocate(this%face(4,this%nelem))
        allocate(elem_free(4,nfree))
        allocate(elem_absorb(4,nabs))
        allocate(this%vol(this%nelem))

        ! calculate min and max of setup coordinates for use of PML
        this%pmlxmin=1.E10              ! define start values for size search
        this%pmlxmax=0.0
        this%pmlymin=1.E10
        this%pmlymax=0.0
        this%pmlzmin=1.E10
        this%pmlzmax=0.0

        !get minmum and maximum dimensions of setup
        filename=trim('cubit/coord')
        open(unit=19,file=trim(filename))
        read(19,*) dummy
        do i=1,this%ncoord
            read(19,*) dummy,this%coord(1,i),this%coord(2,i),this%coord(3,i)
            if (this%coord(1,i) < this%pmlxmin) this%pmlxmin=this%coord(1,i)
            if (this%coord(1,i) > this%pmlxmax) this%pmlxmax=this%coord(1,i)
            if (this%coord(2,i) < this%pmlymin) this%pmlymin=this%coord(2,i)
            if (this%coord(2,i) > this%pmlymax) this%pmlymax=this%coord(2,i)
            if (this%coord(3,i) < this%pmlzmin) this%pmlzmin=this%coord(3,i)
            if (this%coord(3,i) > this%pmlzmax) this%pmlzmax=this%coord(3,i)
        end do
        close(19)

        !print size of setup
        if (par%log) write(*,*) 'x axis ranges from ', this%pmlxmin, ' to ', this%pmlxmax
        if (par%log) write(*,*) 'y axis ranges from ', this%pmlymin, ' to ', this%pmlymax
        if (par%log) write(*,*) 'z axis ranges from ', this%pmlzmin, ' to ', this%pmlzmax
    
        !read mesh file
        filename=trim('cubit/mesh')
        open(unit=19,file=trim(filename))
        read(19,*) dummy
        do i=1,this%nelem
            read(19,*) dummy, this%elem(1,i),this%elem(2,i),this%elem(3,i),this%elem(4,i),this%vp(i),this%vs(i),this%rho(i), this%qp(i), this%qs(i)
            if (this%vp(i)<=0. .or. this%vs(i)<=0.) then
                write(*,*) 'ERROR! vp or vs is zero or negative for at least one element. Stop!', i, this%vp(i), this%vs(i)
                stop
            endif
            if (this%vp(i)<=this%vs(i)) then
                write(*,*) 'ERROR! vp is not larger than vs for at least one element. Stop!', i, this%vp(i), this%vs(i)
                stop
            endif
            if (this%vp(i)<=sqrt(2.)*this%vs(i).and..not.checked) then
                write(*,*) "Poisson number is zero or negative for at least one element. Do you want to continue? (type '0': yes, anything else: no)"
                read(*,*) screeninfo
                if (trim(screeninfo)/='0') then
                    write(*,*) 'Abort!'
                    stop
                else
                    write(*,*) 'Continue!'
                    checked=.true.
                endif
            endif
        end do
        close(19)

        !get neighbors
        call tet_mesh_neighbor_tets(4,this%nelem,this%elem,neighbortemp)
        !need to reorder
        this%neighbor(1,:)=neighbortemp(4,:)
        this%neighbor(2,:)=neighbortemp(3,:)
        this%neighbor(3,:)=neighbortemp(1,:)
        this%neighbor(4,:)=neighbortemp(2,:)
 
        ! material values
        filename=trim('cubit/matprop')
        open(unit=19,file=trim(filename))
        read(19,*) this%nmat
        allocate(this%matval(this%nmat,11))
        do i=1,this%nmat
            ! vp,vs,rho,qp,qs
            read(19,*) dummy,dummy,this%matval(i,2),this%matval(i,3),this%matval(i,1),this%matval(i,6),this%matval(i,7)
            !compute mu and lambda
            this%matval(i,4) = this%matval(i,3)**2 * this%matval(i,1) !mu
            this%matval(i,5) = this%matval(i,2)**2 * this%matval(i,1) - 2 * this%matval(i,4) !lambda
        end do
        close(19)
        do i=1,this%nelem
            this%mu(i)=this%rho(i) * this%vs(i)**2
            this%lambda(i)=this%rho(i) * this%vp(i)**2 - 2 * this%mu(i)
        end do
        ! anelastic precalculation
        allocate(wl(nMB))
        allocate(yl(nMB,this%nmat))
        allocate(ym(nMB,this%nmat))
        if (par%attenuation) then
            if (par%log) write(*,*) "start anelastic precalculation"
            m = 2*nMB-1                          ! number of used frequencies
            call anelasticPreCalc(par%f0_att,par%fr_att,this%wl,this%wk)
            !do i=1,this%nmat                     ! for all materials
            do ie=1,this%nelem
                call calcAnelasticCoefficients(this%qp(ie),this%qs(ie),this%vp(ie),this%vs(ie),this%rho(ie),par%f0_att,par%fr_att,this%wl,this%wk,this%ylambda(:,ie),this%ymu(:,ie),this%lambdau(ie),this%muu(ie))
                this%vpu(ie) = sqrt((this%lambdau(ie)+2*this%muu(ie))/this%rho(ie)) ! get vp_inf and vs_inf
                this%vsu(ie) = sqrt(this%muu(ie)/this%rho(ie))
            end do
            if (par%log) write(*,*) "end anelastic pre calculation"
        end if

    
        ! materials array
        filename=trim('cubit/mat')
        open(unit=19,file=trim(filename))
        do i=1,this%nelem
            read(19,*) c,this%mat(c), this%pml(c)
        end do
        close(19)

        ! free surface array
        filename=trim('cubit/free')
        open(unit=19,file=trim(filename))
        read(19,*) dummy
        do i=1,nfree
            read(19,*) elem_free(1,i),elem_free(2,i),elem_free(3,i),elem_free(4,i)
        end do
        close(19)

        ! absorb surface array
        filename=trim('cubit/absorb')
        open(unit=19,file=trim(filename))
        read(19,*) dummy
        do i=1,nabs
            read(19,*) elem_absorb(1,i),elem_absorb(2,i),elem_absorb(3,i),elem_absorb(4,i)
        end do
        close(19)

        ! write meshfile
        if (par%vtk) then
            filename="/mesh.vtk"
            filename=trim(outpath)//trim(filename)
            write(*,*) "write ", filename
            call writeVtkTetraMeshIntdata(filename, this%elem, this%coord, this%mat)

            filename="/vp.vtk"
            filename=trim(outpath)//trim(filename)
            write(*,*) "write ", filename
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vp)

            filename="/vs.vtk"
            filename=trim(outpath)//trim(filename)
            write(*,*) "write ", filename
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vs)

            filename="/rho.vtk"
            filename=trim(outpath)//trim(filename)
            write(*,*) "write ", filename
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%rho)
        end if

        ! get Volume of elements

        !    do i=1,this%nelem
        !        teta=(this%coord(1,this%elem(1,i))-this%coord(1,this%elem(2,i)))**2+(this%coord(2,this%elem(1,i))-this%coord(2,this%elem(2,i)))**2+(this%coord(3,this%elem(1,i))-this%coord(3,this%elem(2,i)))**2
        !        tetb=(this%coord(1,this%elem(2,i))-this%coord(1,this%elem(3,i)))**2+(this%coord(2,this%elem(2,i))-this%coord(2,this%elem(3,i)))**2+(this%coord(3,this%elem(2,i))-this%coord(3,this%elem(3,i)))**2
        !        tetc=(this%coord(1,this%elem(1,i))-this%coord(1,this%elem(3,i)))**2+(this%coord(2,this%elem(1,i))-this%coord(2,this%elem(3,i)))**2+(this%coord(3,this%elem(1,i))-this%coord(3,this%elem(3,i)))**2
        !        tetap=(this%coord(1,this%elem(3,i))-this%coord(1,this%elem(4,i)))**2+(this%coord(2,this%elem(3,i))-this%coord(2,this%elem(4,i)))**2+(this%coord(3,this%elem(3,i))-this%coord(3,this%elem(4,i)))**2
        !        tetbp=(this%coord(1,this%elem(1,i))-this%coord(1,this%elem(4,i)))**2+(this%coord(2,this%elem(1,i))-this%coord(2,this%elem(4,i)))**2+(this%coord(3,this%elem(1,i))-this%coord(3,this%elem(4,i)))**2
        !        tetcp=(this%coord(1,this%elem(2,i))-this%coord(1,this%elem(4,i)))**2+(this%coord(2,this%elem(2,i))-this%coord(2,this%elem(4,i)))**2+(this%coord(3,this%elem(2,i))-this%coord(3,this%elem(4,i)))**2
        !        fa=tetb+tetbp+tetc+tetcp-teta-tetap
        !        fb=teta+tetap+tetc+tetcp-tetb-tetbp
        !        fc=teta+tetap+tetb+tetbp-tetc-tetcp
        !        tetdel=teta*tetb*tetc+teta*tetbp*tetcp+tetap*tetb*tetcp+tetap*tetbp*tetc
        !        this%vol(i)=sqrt(teta*tetap*fa+tetb*tetbp*fb+tetc*tetcp*fc-tetdel)/12.
        !    enddo

        !"--------------------------------------------------------------------------------------"
        !"--------------------------------------------------------------------------------------"
        ! create partions with metis
        !"--------------------------------------------------------------------------------------"
        !"--------------------------------------------------------------------------------------"

        ! get epsilon (small value) to calculate points at the same position in adjacent elements
        allocate(epsil(this%nelem))
        v_long_edge=1.E10
        do ie=1,this%nelem
            rtemp1=0.
            do i=1,4
                do j=i+1,4
                    rtemp1=rtemp1+sqrt((this%coord(1,this%elem(i,ie))-this%coord(1,this%elem(j,ie)))**2+(this%coord(2,this%elem(i,ie))-this%coord(2,this%elem(j,ie)))**2+(this%coord(3,this%elem(i,ie))-this%coord(3,this%elem(j,ie)))**2)
                enddo
            enddo
            if (this%vs(ie)/rtemp1*6. < v_long_edge) v_long_edge=this%vs(ie)/rtemp1*6.    ! find smallest ration of v_s to edge length
            epsil(ie)=rtemp1/6./30.    ! 1 percent of the average edge length
        enddo



        if (par%nproc > 0) then
            ! build element vector
            allocate(elmnts(this%nelem*4))
            k=1
            do i=1,this%nelem
                do j=1,4
                    elmnts(k)=this%elem(j,i)
                    k=k+1
                end do
            end do

            ! create compact neighbor array
            allocate(xadj(this%nelem+1))
            allocate(adjncy(max_neighbor*this%nelem))
            adjncy(:) = 0
            xadj(:)=0
            xadj(1)=1
            k=0
            do i=1,this%nelem
                nb_neigh=0
                nb_temp(:)=0
                do j=1,4
                    if (this%neighbor(j,i) > 0) then
                        nb_neigh=nb_neigh+1
                        nb_temp(nb_neigh)=this%neighbor(j,i)
                    end if
                end do
                l=1
                do j=k+1,k+nb_neigh
                    adjncy(j) = nb_temp(l)
                    l=l+1
                end do
                k=k+nb_neigh
                xadj(i+1)=k+1
            end do

            ! create metis partition
            allocate(part(this%nelem))
            options(:)=0

            call METIS_PartGraphRecursive(this%nelem, xadj, adjncy, 0, 0, 0, 1, par%nproc,options, edgecut, part)

            ! create glob2loc_elem
            ! inspired by the specfem partition (SPECFEM3D 2.0 in subroutine part_decompose_mesh_SCOTCH.f90
            ! written by Komatitsch at al. published under GPL taken from www.geodynamics.org)
            ! be careful, metis gives parts begining from 0 or from 1 depending on the version, so I compile a local version, here 4.0.3
            allocate(glob2loc_elmnts(this%nelem))
            glob2loc_elmnts(:) = 0

            allocate(num_loc(par%nproc))
            num_loc(:) = 1

            ! build a vector with the local numbering of the elements
            do num_glob=1,this%nelem
                num_part=part(num_glob)
                glob2loc_elmnts(num_glob) = num_loc(num_part)
                num_loc(num_part) = num_loc(num_part) + 1
            end do

            ! create local node numbering
            ! 2 vectors, nnodes with the positions of the elements in the vector nodes_elements, similar to the metis numbering
            allocate(nnodes_elmnts(this%ncoord))
            allocate(nodes_elmnts(this%ncoord*nsize))
            nnodes_elmnts(:)=0
            nodes_elmnts(:)=0
            do i=1, 4*this%nelem
                ! nodes is like a matrix with nodes as rows and nsize elements as colums
                nodes_elmnts((elmnts(i)-1)*nsize + nnodes_elmnts(elmnts(i))+1) = 1+(i-1)/4
                nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) +1
            end do

            ! create the local node numbering
            allocate(glob2loc_nodes_nparts(this%ncoord+1))
            allocate(part_nodes(par%nproc), num_parts(par%nproc))

            size_glob2loc_nodes = 1
            part_nodes(:) = 0

            do in=1,this%ncoord
                glob2loc_nodes_nparts(in) = size_glob2loc_nodes
                do ie = 1, nnodes_elmnts(in)
                    ! shows to which partitions a node belongs
                    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
                end do

                do num_part=1,par%nproc
                    if (part_nodes(num_part) == 1) then
                        ! get number of nodes in the global array, there might be some nodes which count double at the interfaces
                        size_glob2loc_nodes = size_glob2loc_nodes +1
                        part_nodes(num_part) = 0
                    end if
                end do
                glob2loc_nodes_nparts(this%ncoord+1) = size_glob2loc_nodes
            end do

            allocate(glob2loc_nodes_parts(glob2loc_nodes_nparts(this%ncoord+1)-1))
            allocate(glob2loc_nodes(glob2loc_nodes_nparts(this%ncoord+1)-1))
            allocate(glob2loc_nodes2(this%ncoord))

            glob2loc_nodes(:) = 1
            part_nodes(:) = 0
            num_parts(:)=1
            size_glob2loc_nodes = 1

            do in=1,this%ncoord
                do ie = 1, nnodes_elmnts(in)
                    part_nodes(part(nodes_elmnts((ie)+nsize*(in-1))))=1
                end do

                do num_part=1,par%nproc
                    if (part_nodes(num_part) == 1) then
                        ! build arrays with local nodes ordering
                        glob2loc_nodes_parts(size_glob2loc_nodes) = num_part
                        glob2loc_nodes(size_glob2loc_nodes) = num_parts(num_part)
                        size_glob2loc_nodes = size_glob2loc_nodes +1
                        num_parts(num_part) = num_parts(num_part) +1
                        part_nodes(num_part) = 0
                    end if
                end do
            end do
        else
            STOP "STOP. NexD 3D works only as MPI version. The calculations are too expensive for running them on just one processor."
        end if


        !"--------------------------------------------------------------------------------------"
        ! build boundarys for external mesh
        !"--------------------------------------------------------------------------------------"
    
        do i=1, this%nelem
            do j=1,4
                if (this%neighbor(j,i) <=0) then
                    this%neighbor(j,i) = 0
                end if
            end do
        end do

        ! find connections of faces
        this%face(:,:) = 0
        do i=1,this%nelem
            do j=1,4
                c=this%neighbor(j,i)
                if (c>0) then
                    do k=1,4
                        d=this%neighbor(k,c)
                        if (d==i) then
                            this%face(j,i)=k
                        end if
                    end do
                end if
            end do
        end do

        ! assign free faces, free=-1
        do i=1, nfree
            c = elem_free(1,i)
            check=.false.
            do j=1,4
                do k=1,3
                    if (this%elem(j,c)==elem_free(k+1,i)) then
                        check(j) = .true.
                    end if
                end do
            end do
            do j=1,4
                if (.not.check(j)) then
                    this%face(side(j),c)=-1
                end if
            end do
        end do

        ! assign absorbing faces, absorb=-2
        checked=.false.
        do i=1, nabs
            c = elem_absorb(1,i)
            check=.false.
            do j=1,4
                do k=1,3
                    if (this%elem(j,c)==elem_absorb(k+1,i)) then
                        check(j) = .true.
                    end if
                end do
            end do
            do j=1,4
                if (.not.check(j)) then
                    ! check for double assigned faces
                    if (this%face(side(j),c) == -1) then
                        write(*,*) 'Face ',j,' of element ',c,' is assigned as free and absorbing boundary at the same time!'
                        checked=.true.
                    endif
                    this%face(side(j),c)=-2
                end if
            end do
            if (sum(this%face(:,c))==-8) then
                write(*,*) 'Found single element with only absorbing boundaries!'
                stop
            end if
        end do

        ! check if all domain boundaries assigned

        do i=1, this%nelem
            do j=1,4
                if (this%neighbor(j,i) == 0) then
                    if (this%face(j,i) > -1) then
                        write(*,*) 'Face ',j,' of element ',i,' is a domain boundary but not assigned as free or absorbing boundary'
                        checked=.true.
                    endif
                end if
            end do
        end do
        if (checked) then
            write(*,*) 'ABORT!'
            stop
        endif
        !"--------------------------------------------------------------------------------------"
        ! start to create local elements
        !"--------------------------------------------------------------------------------------"

        call nodes3D(balpha(order),x,y,z)
        call xyzToRst(x,y,z,r,s,t)

        minGLL = abs(r(1)-r(2))
        if (par%log) write(*,*) "min gll point distance:           : ", minGLL

        call vdm3d(this%vdm,r,s,t)
        call invVdm3D(this%vdm,this%invvdm,0)
        this%mass = matmul(transpose(this%invvdm),this%invvdm)
        call dmatrices3d(this%Dr,this%Ds,this%Dt,r,s,t,this%vdm)
        call gradVdm3d(this%vdmr,this%Vdms,this%Vdmt,r,s,t)

        tempM1=matmul(this%vdm,transpose(this%vdm))
        tempM2=matmul(this%vdm,transpose(this%vdmr))
        call invert(tempM1)
        this%Drw=matmul(tempM2,tempM1)

        tempM2=matmul(this%vdm,transpose(this%vdms))
        this%Dsw=matmul(tempM2,tempM1)

        tempM2=matmul(this%vdm,transpose(this%vdmt))
        this%Dtw=matmul(tempM2,tempM1)
    
    
        ! find boundary nodes in the local element
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(t(i)+1.0) < EPS) then
                this%fmask(d,1)=c
                d=d+1
            end if
        end do
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(s(i)+1.0) < EPS) then
                this%fmask(d,2)=c
                d=d+1
            end if
        end do
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(1.0+r(i)+s(i)+t(i)) < EPS) then
                this%fmask(d,3)=c
                d=d+1
            end if
        end do
        c=0
        d=1
        do i=1,Np
            c=c+1
            if (abs(1.0+r(i)) < EPS) then
                this%fmask(d,4)=c
                d=d+1
            end if
        end do

        call lift3d(this%lift,this%fmask,this%vdm,r,s,t)


        !"--------------------------------------------------------------------------------------"
        ! build Databases
        !"--------------------------------------------------------------------------------------"

        allocate(loc2glob_elemnts(par%nproc,maxval(num_loc(:)-1)))
        do iproc=1,par%nproc

            ! build local to global numbering
            do ie=1,this%nelem
                if (part(ie)== iproc) then
                    loc2glob_elemnts(iproc,glob2loc_elmnts(ie))=ie
                endif
            enddo

            !prepare database for every partition
            db(iproc)%nglob = (num_loc(iproc)-1)*Np
            db(iproc)%nelem = num_loc(iproc)-1
            db(iproc)%nmat = this%nmat
            db(iproc)%ncoord = num_parts(iproc)-1

            db(iproc)%invVdm = this%invVdm
            db(iproc)%Vdm = this%Vdm
            db(iproc)%Dr = this%Dr
            db(iproc)%Ds = this%Ds
            db(iproc)%Dt = this%Dt
            db(iproc)%Drw = this%Drw
            db(iproc)%Dsw = this%Dsw
            db(iproc)%Dtw = this%Dtw
            db(iproc)%Vdmr = this%Vdmr
            db(iproc)%Vdms = this%Vdms
            db(iproc)%Vdmt = this%Vdmt
            db(iproc)%mass = this%mass
            db(iproc)%lift = this%lift
            db(iproc)%fmask = this%fmask
            db(iproc)%pml_delta=this%pml_delta

            allocate(db(iproc)%matval(db(iproc)%nmat,11))
            allocate(db(iproc)%mat(db(iproc)%nelem))
            allocate(db(iproc)%coord(3,db(iproc)%ncoord))
            allocate(db(iproc)%elem(4,db(iproc)%nelem))
            allocate(db(iproc)%neighbor(4,db(iproc)%nelem))
            allocate(db(iproc)%vp(db(iproc)%nelem),db(iproc)%vs(db(iproc)%nelem),db(iproc)%rho(db(iproc)%nelem),&
                db(iproc)%mu(db(iproc)%nelem),db(iproc)%lambda(db(iproc)%nelem))
            allocate(db(iproc)%vpu(db(iproc)%nelem),db(iproc)%vsu(db(iproc)%nelem),db(iproc)%qp(db(iproc)%nelem),db(iproc)%qs(db(iproc)%nelem),&
                db(iproc)%muu(db(iproc)%nelem),db(iproc)%lambdau(db(iproc)%nelem))
            allocate(db(iproc)%ylambda(nMB,db(iproc)%nelem),db(iproc)%ymu(nMB,db(iproc)%nelem),db(iproc)%wl(nMB))
            allocate(db(iproc)%pml(db(iproc)%nelem))
            allocate(db(iproc)%vx(db(iproc)%nglob),db(iproc)%vy(db(iproc)%nglob),db(iproc)%vz(db(iproc)%nglob))
            allocate(db(iproc)%rx(db(iproc)%nglob), db(iproc)%ry(db(iproc)%nglob),db(iproc)%rz(db(iproc)%nglob), &
                db(iproc)%sx(db(iproc)%nglob),db(iproc)%sy(db(iproc)%nglob), db(iproc)%sz(db(iproc)%nglob), &
                db(iproc)%tx(db(iproc)%nglob),db(iproc)%ty(db(iproc)%nglob), db(iproc)%tz(db(iproc)%nglob))
            allocate(db(iproc)%jacobian(db(iproc)%nglob))
            allocate(db(iproc)%ibool(Np,db(iproc)%nelem))
            allocate(db(iproc)%nx(4*NpF,db(iproc)%nelem),db(iproc)%ny(4*NpF,db(iproc)%nelem), db(iproc)%nz(4*NpF,db(iproc)%nelem), &
                db(iproc)%sj(4*NpF,db(iproc)%nelem),db(iproc)%fscale(4*NpF,db(iproc)%nelem))
            allocate(db(iproc)%nnormal(3,4,db(iproc)%nelem))
            allocate(db(iproc)%stangential(3,4,db(iproc)%nelem))
            allocate(db(iproc)%ttangential(3,4,db(iproc)%nelem))
            allocate(db(iproc)%mpi_interface(4,4,db(iproc)%nelem))
            allocate(db(iproc)%face(4,db(iproc)%nelem))
            allocate(db(iproc)%ibt(NpF,4,db(iproc)%nelem),db(iproc)%ibn(NpF,4,db(iproc)%nelem))
            allocate(db(iproc)%mpi_ibt(NpF,4,db(iproc)%nelem))
            allocate(db(iproc)%loc2glob_nodes(db(iproc)%ncoord))
            !       allocate(db(iproc)%vol(db(iproc)%nelem)

            ! build local arrays and coordinates
            do i=1, this%ncoord
                do j= glob2loc_nodes_nparts(i), glob2loc_nodes_nparts(i+1)-1
                    if (glob2loc_nodes_parts(j) == iproc) then
                        k=glob2loc_nodes(j)
                        db(iproc)%coord(1,k)=this%coord(1,i)
                        db(iproc)%coord(2,k)=this%coord(2,i)
                        db(iproc)%coord(3,k)=this%coord(3,i)
                        db(iproc)%loc2glob_nodes(k)=i
                    end if
                end do
            end do

            ! build local elements
            do i=1, this%nelem
                if (part(i) == iproc) then
                    do j=1,4
                        l=elmnts((i-1)*4+j)
                        do k=glob2loc_nodes_nparts(l), glob2loc_nodes_nparts(l+1)-1
                            if (glob2loc_nodes_parts(k) == iproc) then
                                loc_nodes(j) = glob2loc_nodes(k)
                            end if
                        end do
                    end do
                    k = glob2loc_elmnts(i)
                    db(iproc)%elem(1,k) = loc_nodes(1)
                    db(iproc)%elem(2,k) = loc_nodes(2)
                    db(iproc)%elem(3,k) = loc_nodes(3)
                    db(iproc)%elem(4,k) = loc_nodes(4)
                end if
            end do

            !"--------------------------------------------------------------------------------------"
            ! find mpi neighbors and the neighbor element
            !"--------------------------------------------------------------------------------------"

            db(iproc)%neighbor(:,:) = 0
            db(iproc)%mpi_interface(:,:,:) = 0
            db(iproc)%pinterfaces=0
            do ie=1,this%nelem
                do j=1,4
                    k=this%neighbor(j,ie)
                    if (k>0) then
                        if (part(ie) == iproc) then
                            if (part(k) == iproc) then
                                l = glob2loc_elmnts(ie)
                                m = glob2loc_elmnts(k)
                                db(iproc)%neighbor(j,l) = m
                            else
                                l = glob2loc_elmnts(ie)
                                m = glob2loc_elmnts(k)
                                db(iproc)%neighbor(j,l) = -1 ! means mpi neighbor
                                db(iproc)%mpi_interface(1,j,l) = part(k) ! neighbor partition
                                db(iproc)%mpi_interface(2,j,l) = m ! element in partition
                                db(iproc)%pinterfaces=db(iproc)%pinterfaces+1
                            end if
                        end if
                    end if
                end do
            end do

               !"--------------------------------------------------------------------------------------"
               ! transform elements to reference element
               !"--------------------------------------------------------------------------------------"

            db(iproc)%face(:,:) = 0
            c=0
            do i=1,db(iproc)%nelem
                do j=1,NP
                    c=c+1
                    ! get coordinates of interpolation points
                    db(iproc)%vx(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * db(iproc)%coord(1,db(iproc)%elem(1,i)) + (1.0+r(j)) * db(iproc)%coord(1,db(iproc)%elem(2,i)) +&
                        (1+s(j)) * db(iproc)%coord(1,db(iproc)%elem(3,i)) + (1+t(j)) * db(iproc)%coord(1,db(iproc)%elem(4,i)))
                    db(iproc)%vy(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * db(iproc)%coord(2,db(iproc)%elem(1,i)) + (1.0+r(j)) * db(iproc)%coord(2,db(iproc)%elem(2,i)) +&
                        (1+s(j)) * db(iproc)%coord(2,db(iproc)%elem(3,i)) + (1+t(j)) * db(iproc)%coord(2,db(iproc)%elem(4,i)))
                    db(iproc)%vz(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * db(iproc)%coord(3,db(iproc)%elem(1,i)) + (1.0+r(j)) * db(iproc)%coord(3,db(iproc)%elem(2,i)) +&
                        (1+s(j)) * db(iproc)%coord(3,db(iproc)%elem(3,i)) + (1+t(j)) * db(iproc)%coord(3,db(iproc)%elem(4,i)))
                    ! set up ibool (transforms local numbering of interpolation points in each element to a global numbering)
                    db(iproc)%ibool(j,i)=c
                end do
            end do


            ! get geometric factors
            do i=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,i)
                call geometricFactors3d(rx,sx,tx,ry,sy,ty,rz,sz,tz,jacobian,db(iproc)%vx(iv),db(iproc)%vy(iv),db(iproc)%vz(iv),db(iproc)%Dr,db(iproc)%Ds,db(iproc)%Dt)
                db(iproc)%rx(iv)=rx
                db(iproc)%sx(iv)=sx
                db(iproc)%tx(iv)=tx
                db(iproc)%ry(iv)=ry
                db(iproc)%sy(iv)=sy
                db(iproc)%ty(iv)=ty
                db(iproc)%rz(iv)=rz
                db(iproc)%sz(iv)=sz
                db(iproc)%tz(iv)=tz
                db(iproc)%jacobian(iv)=jacobian
            end do

            ! get normals
            do i=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,i)
                call normals3d(db(iproc)%nx(:,i),db(iproc)%ny(:,i),db(iproc)%nz(:,i),db(iproc)%sj(:,i),db(iproc)%Dr,db(iproc)%Ds,db(iproc)%Dt,db(iproc)%vx(iv),db(iproc)%vy(iv),db(iproc)%vz(iv),db(iproc)%fmask)
            end do

            ! calculate n: normal, s: tangential, t: tangential vector of one face
            do ie=1,db(iproc)%nelem
                do is=1,4
                    if (is == 1) then
                        c=1
                        d=2
                        e=0*NpF+1
                    else if (is == 2 ) then
                        c=1
                        d=2
                        e=1*NpF+1
                    else if (is == 3 ) then
                        c=2
                        d=3
                        e=2*NpF+1
                    else if (is == 4 ) then
                        c=1
                        d=3
                        e=3*NpF+1
                    end if
                    db(iproc)%nnormal(1,is,ie) =db(iproc)%nx(e,ie)
                    db(iproc)%nnormal(2,is,ie) =db(iproc)%ny(e,ie)
                    db(iproc)%nnormal(3,is,ie) =db(iproc)%nz(e,ie)
                    rtemp1 = db(iproc)%coord(1,db(iproc)%elem(d,ie)) - db(iproc)%coord(1,db(iproc)%elem(c,ie))
                    rtemp2 = db(iproc)%coord(2,db(iproc)%elem(d,ie)) - db(iproc)%coord(2,db(iproc)%elem(c,ie))
                    rtemp3 = db(iproc)%coord(3,db(iproc)%elem(d,ie)) - db(iproc)%coord(3,db(iproc)%elem(c,ie))
                    rtemp4= sqrt(rtemp1**2+rtemp2**2+rtemp3**2)
                    db(iproc)%stangential(1,is,ie)=rtemp1/rtemp4
                    db(iproc)%stangential(2,is,ie)=rtemp2/rtemp4
                    db(iproc)%stangential(3,is,ie)=rtemp3/rtemp4

                    rtemp1 = db(iproc)%nnormal(2,is,ie) * db(iproc)%stangential(3,is,ie) - db(iproc)%nnormal(3,is,ie)*db(iproc)%stangential(2,is,ie)
                    rtemp2 = db(iproc)%nnormal(3,is,ie) * db(iproc)%stangential(1,is,ie) - db(iproc)%nnormal(1,is,ie)*db(iproc)%stangential(3,is,ie)
                    rtemp3 = db(iproc)%nnormal(1,is,ie) * db(iproc)%stangential(2,is,ie) - db(iproc)%nnormal(2,is,ie)*db(iproc)%stangential(1,is,ie)
                    rtemp4= sqrt(rtemp1**2+rtemp2**2+rtemp3**2)
                    db(iproc)%ttangential(1,is,ie)=rtemp1/rtemp4
                    db(iproc)%ttangential(2,is,ie)=rtemp2/rtemp4
                    db(iproc)%ttangential(3,is,ie)=rtemp3/rtemp4
                end do
            end do

            ! make scaling vector
            do i=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,i)
                c=1
                do k=1,4
                    do j=1,NpF
                        d=iv(db(iproc)%fmask(j,k))
                        db(iproc)%fmaskv(c)=db(iproc)%fmask(j,k)
                        db(iproc)%fscale(c,i)=db(iproc)%sj(c,i)/db(iproc)%jacobian(d)
                        c=c+1
                    end do
                end do
            end do

            ! add physical parameters to the databases
            do ie=1, this%nelem
                if (part(ie)== iproc) then
                    l = glob2loc_elmnts(ie)
                    db(iproc)%face(:,l) = this%face(:,ie)
                    db(iproc)%mat(l)=this%mat(ie)
                    db(iproc)%rho(l)=this%rho(ie)
                    db(iproc)%vp(l)=this%vp(ie)
                    db(iproc)%vs(l)=this%vs(ie)
                    db(iproc)%qp(l)=this%qp(ie)
                    db(iproc)%qs(l)=this%qs(ie)
                    db(iproc)%vpu(l)=this%vpu(ie)
                    db(iproc)%vsu(l)=this%vsu(ie)
                    db(iproc)%lambdau(l)=this%lambdau(ie)
                    db(iproc)%muu(l)=this%muu(ie)
                    db(iproc)%ylambda(:,l)=this%ylambda(:,ie)
                    db(iproc)%ymu(:,l)=this%ymu(:,ie)
                    db(iproc)%wl(:)=this%wl(:)
                    db(iproc)%pml(l)=this%pml(ie)
                    !                db(iproc)%vol(l)=this%vol(ie)
                    if (par%attenuation) then
                        if (this%vs(l) .gt. 0) then
                            db(iproc)%mu(l)=this%muu(ie)
                        else
                            db(iproc)%mu(l)=0.0
                        end if
                        db(iproc)%lambda(l)=this%lambdau(ie)
                    else
                        if (this%vs(l) .gt. 0) then
                            db(iproc)%mu(l)=this%mu(ie)
                        else
                            db(iproc)%mu(l)=0.0
                        end if
                        db(iproc)%lambda(l)=this%lambda(ie)
                    endif
                end if
            end do
            if (.not.par%attenuation) then
                db(iproc)%vpu=db(iproc)%vp
                db(iproc)%vsu=db(iproc)%vs
                db(iproc)%lambdau=db(iproc)%lambda
                db(iproc)%muu=db(iproc)%mu
            endif

            !"--------------------------------------------------------------------------------------"
            ! recreate boundray nodes association with local numbering
            !"--------------------------------------------------------------------------------------"

            db(iproc)%ibn(:,:,:) =0
            do ie=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,ie)
                do is=1,4
                    in=db(iproc)%neighbor(is,ie)
                    if (in>0) then
                        ivn=db(iproc)%ibool(:,in)
                        j=db(iproc)%face(is,ie)
                        db(iproc)%ibt(:,is,ie)=iv(db(iproc)%fmask(:,is))
                        db(iproc)%ibn(:,is,ie)=ivn(db(iproc)%fmask(:,j))
                    else if( in<=0) then
                        ! for mpi boundaries we must also reorder (later)
                        db(iproc)%ibt(:,is,ie)=iv(db(iproc)%fmask(:,is))
                    end if
                end do
            end do

            ! we have to reorder some boundary arrays to match the points on the face
            do ie=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,ie)
                do is=1,4
                    in=db(iproc)%neighbor(is,ie)
                    if (in>0) then
                        do j=1,4
                            do i=1,Npf
                                if ( (abs(db(iproc)%vx(db(iproc)%ibt(i,is,ie)) - db(iproc)%vx(db(iproc)%ibn(i,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))&
                                    .and. (abs(db(iproc)%vy(db(iproc)%ibt(i,is,ie)) - db(iproc)%vy(db(iproc)%ibn(i,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie)))) &
                                    .and. (abs(db(iproc)%vz(db(iproc)%ibt(i,is,ie)) - db(iproc)%vz(db(iproc)%ibn(i,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))) then
                                    perm(i)=i
                                else
                                    c=1
                                    do k=1,Npf
                                        if ( (abs(db(iproc)%vx(db(iproc)%ibt(i,is,ie)) - db(iproc)%vx(db(iproc)%ibn(k,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))&
                                            .and. (abs(db(iproc)%vy(db(iproc)%ibt(i,is,ie)) - db(iproc)%vy(db(iproc)%ibn(k,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie)))) &
                                            .and. (abs(db(iproc)%vz(db(iproc)%ibt(i,is,ie)) - db(iproc)%vz(db(iproc)%ibn(k,is,ie)) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))) then
                                            perm(i)=c
                                            c=c+1
                                        else
                                            c=c+1
                                        end if
                                    end do
                                end if
                            end do
                            ! do permutation if not correct
                            db(iproc)%ibn(:,is,ie)=db(iproc)%ibn(perm(:),is,ie)
                        end do
                    end if
                end do
            end do
        end do ! iproc

        !"--------------------------------------------------------------------------------------"
        ! create mpi interfaces
        !"--------------------------------------------------------------------------------------"
        allocate(icom(par%nproc,par%nproc))
        icom(:,:) = 0
        do iproc=1,par%nproc
            do ie=1,db(iproc)%nelem
                do i=1,4
                    if (db(iproc)%mpi_interface(1,i,ie) > 0) then
                        icom(iproc,db(iproc)%mpi_interface(1,i,ie))=1
                    end if
                end do
            end do
            db(iproc)%mpi_nn=sum(icom(iproc,:))
        end do

        ! set up mpi_neighbor array
        do iproc=1,par%nproc
            allocate(db(iproc)%mpi_neighbor(db(iproc)%mpi_nn))
            c=1
            do i=1,par%nproc
                if (icom(iproc,i) == 1) then
                    db(iproc)%mpi_neighbor(c)=i
                    c=c+1
                end if
            end do
        end do

        ! howmany mpi interfaces do we have?
        do iproc=1,par%nproc
            allocate(tempv(db(iproc)%mpi_nn))
            allocate(db(iproc)%mpi_ninterface(db(iproc)%mpi_nn))
            tempv(:)=0
            do i=1,db(iproc)%mpi_nn
                in=db(iproc)%mpi_neighbor(i)
                do ie=1,db(iproc)%nelem
                    do is=1,4
                        if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
                            tempv(i)=tempv(i)+1
                        end if
                    end do
                end do
!                write(*,*) "proc",iproc," found", tempv(i), "elements on interface between", iproc ,"and", in
                db(iproc)%mpi_ninterface(i)=tempv(i)
            end do
            deallocate(tempv)
        end do

        if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
        if (par%log) write(*,*) "found total", sum(db(:)%mpi_nn), "mpi interfaces"

       
        temp2=0
        do iproc=1,par%nproc
            temp1=db(iproc)%nelem
            if (temp2<temp1) temp2=temp1
        end do

        ! set up mpi_ibool with max mpi neighbors
        do iproc=1,par%nproc
            db(iproc)%mpi_nemax=temp2
            allocate(db(iproc)%mpi_ibool(4,db(iproc)%mpi_nemax))
            db(iproc)%mpi_ibool(:,:)=0
        end do

        temp2=0
        do iproc=1,par%nproc
            temp1=maxval(db(iproc)%mpi_ninterface(:))
            if (temp2<temp1) temp2=temp1
        end do

        !"--------------------------------------------------------------------------------------"
        ! create main mpi reconnection arrays.
        !"--------------------------------------------------------------------------------------"
        do iproc=1,par%nproc
            db(iproc)%mpi_nnmax=maxval(db(:)%mpi_nn)
            allocate(db(iproc)%mpi_icon(db(iproc)%mpi_nnmax))
            db(iproc)%mpi_icon(:) = 0
            db(iproc)%mpi_ne=temp2
            allocate(db(iproc)%mpi_connection(db(iproc)%mpi_nn,db(iproc)%mpi_ne,2))
            db(iproc)%mpi_connection(:,:,:)=0
            do i=1,db(iproc)%mpi_nn
                in=db(iproc)%mpi_neighbor(i)
                c=1
                do ie=1,db(iproc)%nelem
                    do is=1,4
                        if ( (db(iproc)%neighbor(is,ie) == -1) .and. (db(iproc)%mpi_interface(1,is,ie) == in)) then
                            db(iproc)%mpi_connection(i,c,1)=ie ! mpi interface to global element
                            db(iproc)%mpi_connection(i,c,2)=is ! which side of the element is the interface
                            db(iproc)%mpi_interface(3,is,ie)=c ! global element to local mpi interface
                            db(iproc)%mpi_interface(4,is,ie)=i ! which local interface?
                            c=c+1
                        end if
                    end do
                end do
            end do
        end do

        do iproc=1,par%nproc
            do i=1,db(iproc)%mpi_nn
                l=db(iproc)%mpi_neighbor(i)
                do ie=1,db(iproc)%mpi_ne
                    if ( db(iproc)%mpi_connection(i,ie,1) >0) then
                        is=db(iproc)%mpi_connection(i,ie,2)
                        ee=db(iproc)%mpi_connection(i,ie,1)
                        k=db(iproc)%mpi_interface(2,is,ee)
                        c = db(l)%mpi_interface(3,db(iproc)%face(is,ee),k)
                        db(iproc)%mpi_icon(i)=db(l)%mpi_interface(4,db(iproc)%face(is,ee),k)
                        ! correct mpi boundarys
                        db(iproc)%mpi_ibool(is,ee)=c
                        perm=0
                        do j=1,NpF
                            mpi_ti = db(iproc)%ibt(j,is,ee)
                            mpi_ni = db(l)%ibt(j,db(iproc)%face(is,ee),k)
                            if ( (abs(db(iproc)%vx(mpi_ti) - db(l)%vx(mpi_ni) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie)))) &
                                .and. (abs(db(iproc)%vy(mpi_ti) - db(l)%vy(mpi_ni) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))) ) &
                                .and. (abs(db(iproc)%vz(mpi_ti) - db(l)%vz(mpi_ni))  < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))) then
                                perm(j)=j
                            else
                                d=1
                                do e=1,Npf
                                    mpi_ti = db(iproc)%ibt(e,is,ee)
                                    if ( (abs(db(iproc)%vx(mpi_ti) - db(l)%vx(mpi_ni) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie)))) &
                                        .and. (abs(db(iproc)%vy(mpi_ti) - db(l)%vy(mpi_ni) ) < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))) ) &
                                        .and. (abs(db(iproc)%vz(mpi_ti) - db(l)%vz(mpi_ni))  < epsil(loc2glob_elemnts(iproc,glob2loc_elmnts(ie))))) then
                                        perm(j)=d
                                        d=d+1
                                    else
                                        d=d+1
                                    end if
                                end do
                            end if
                        end do
                        ! do permutation if not correct
                        if (any(perm==0)) then
                            write(*,*) 'Error in permutation of GLL points at MPI boundaries'
                            stop
                        endif
                        db(l)%mpi_ibt(:,db(iproc)%face(is,ee),k)=perm(:)
                    end if
                end do
            end do
        end do

        ! find min dtfactor
        rtemp1=1e7
        do iproc=1,par%nproc
            do ie=1,db(iproc)%nelem
                iv=db(iproc)%ibool(:,ie)
                min1=minval(db(iproc)%jacobian(iv))
                min2=minval(db(iproc)%sj(:,ie))
                if ((min1/min2)<rtemp1) then
                    rtemp1=min1/min2
                end if
            end do
        end do
        rtemp3=1e7
        do iproc=1,par%nproc
            do ie=1,db(iproc)%nelem
                min1=2.0/3.0 * minGLL * (rtemp1/db(iproc)%vpu(ie))
                if ((min1)<rtemp3) then
                    rtemp3=min1
                end if
            end do
        end do
        db(:)%dtfactor = rtemp3
        db(:)%minGLL_jac = minGLL*rtemp1
        if (par%log) write(*,*) "overall min characteristic length  :",rtemp1
        if (par%log) write(*,*) "overall min dt factor (clf=1)      :",rtemp3
        if (par%log) write(*,*) "used time step with clf factor     :",cfl*rtemp3
        if (par%log) write(*,*) "Maximum Frequency for the used mesh:",NGLL/5.*v_long_edge

        ! Save outer boundaries of the mesh for PML
        db(:)%pmlxmin=this%pmlxmin
        db(:)%pmlxmax=this%pmlxmax
        db(:)%pmlymin=this%pmlymin
        db(:)%pmlymax=this%pmlymax
        db(:)%pmlzmin=this%pmlzmin
        db(:)%pmlzmax=this%pmlzmax

        !"--------------------------------------------------------------------------------------"
        ! find source and receivers
        !"--------------------------------------------------------------------------------------"


        call initSource(par,src)
        call initReceiver(par,rec)

        ! we need global coordinages to find the right source and receiver points
        allocate(this%vx(this%nglob),this%vy(this%nglob),this%vz(this%nglob))
        allocate(this%ibool(Np,this%nelem))
        this%ibool(:,:) = 0
        c=0
        do i=1,this%nelem
            do j=1,NP
                c=c+1
                this%vx(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * this%coord(1,this%elem(1,i)) + (1.0+r(j)) * this%coord(1,this%elem(2,i)) +&
                    (1+s(j)) * this%coord(1,this%elem(3,i)) + (1+t(j)) * this%coord(1,this%elem(4,i)))
                this%vy(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * this%coord(2,this%elem(1,i)) + (1.0+r(j)) * this%coord(2,this%elem(2,i)) +&
                    (1+s(j)) * this%coord(2,this%elem(3,i)) + (1+t(j)) * this%coord(2,this%elem(4,i)))
                this%vz(c) = 0.5 * ( -(1.0+r(j)+s(j)+t(j)) * this%coord(3,this%elem(1,i)) + (1.0+r(j)) * this%coord(3,this%elem(2,i)) +&
                    (1+s(j)) * this%coord(3,this%elem(3,i)) + (1+t(j)) * this%coord(3,this%elem(4,i)))
                ! set up ibool
                this%ibool(j,i)=c
            end do
        end do

        !find sources in global mesh
        write(*,*)
        write(*,*) "start finding sources in global mesh"
        write(*,*)
        call findSource(par,src,this%vx,this%vy,this%vz,this%nglob,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, this%Dt)
        ! get global src elements
        allocate(tempv(par%nproc))
        allocate(srcnr(par%nproc,size(src%srcelem)))
        srcnr(:,:)=0
        tempv(:)=0
        j=1
        do i=1,size(src%srcelem)
            db(part(src%srcelem(i)))%has_src=.true.
            tempv(part(src%srcelem(i)))=tempv(part(src%srcelem(i)))+1
            srcnr(part(src%srcelem(i)),tempv(part(src%srcelem(i))))=j
            j=j+1
        end do

        do i=1,size(src%srcelem)
            call getlocalSrcElement(src,i,glob2loc_elmnts(src%srcelem(i)))
        end do
        ! create local source var for proc that carries the source
        j=1
        write(*,*)
        write(*,*) "start searching sources in partitions"
        do iproc=1,par%nproc
            if (db(iproc)%has_src) then
                ! recreate src arrays
                allocate(locsrcxyz(3,tempv(iproc)))
                allocate(locsrcnr(tempv(iproc)))
                allocate(locsrctype(tempv(iproc)))
                allocate(locsrcstf(tempv(iproc)))
                allocate(locsrcf0(tempv(iproc)))
                allocate(locsrcfactor(tempv(iproc)))
                allocate(locsrcangle_force(3,tempv(iproc)))
                allocate(locsrcm(6,tempv(iproc)))
                locnsrc=tempv(iproc)
                do i=1,tempv(iproc)
                    locsrcnr(i) = srcnr(iproc,i)
                    locsrcxyz(:,i) = src%srcxyz(:,locsrcnr(i))
                    locsrctype(i) = src%srctype(locsrcnr(i))
                    locsrcstf(i) = src%srcstf(locsrcnr(i))
                    locsrcf0(i) = src%srcf0(locsrcnr(i))
                    locsrcfactor(i) = src%srcfactor(locsrcnr(i))
                    locsrcangle_force(:,i) = src%srcangle_force(:,locsrcnr(i))
                    locsrcM(:,i) = src%srcM(:,locsrcnr(i))
                    j=j+1
                end do
                call prepareRecalcSrc(tempsrc,locnsrc,locsrcxyz,locsrctype,locsrcstf,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM)
                call findSource(par,tempsrc,db(iproc)%vx,db(iproc)%vy,db(iproc)%vz,db(iproc)%nglob,db(iproc)%nelem,db(iproc)%ibool, &
                    db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds,db(iproc)%dt)
                write(filename,"('/srcVar',i6.6)") iproc
                filename=trim(outpath)//trim(filename)
                call writeSrcVar(tempsrc,filename)
                deallocate(locsrcnr)
                deallocate(locsrcxyz)
                deallocate(locsrctype,locsrcstf,locsrcf0,locsrcfactor,locsrcangle_force,locsrcM)
                call deallocSrcVar(tempsrc)
            end if
        end do
        write(*,*) "end source search"
        write(*,*)

        deallocate(tempv)
        deallocate(srcnr)

        ! find receiver, works as for sources!
        write(*,*) "start finding receiver in global mesh"
        write(*,*)
        call findReceiver(par,rec,this%vx,this%vy,this%vz,this%nglob,this%nelem,this%ibool,this%coord,this%elem, this%Dr, this%Ds, this%Dt)
        write(*,*) "end receiver search"
        ! get global rec elements
        allocate(tempv(par%nproc))
        allocate(recnum(par%nproc,size(rec%recelem)))
        recnum(:,:)=0
        tempv(:)=0
        j=1
        do i=1,size(rec%recelem)
            db(part(rec%recelem(i)))%has_rec=.true.
            tempv(part(rec%recelem(i)))=tempv(part(rec%recelem(i)))+1
            recnum(part(rec%recelem(i)),tempv(part(rec%recelem(i))))=j
            j=j+1
        end do
    
        !change element number from global to local numerbing
        do i=1,size(rec%recelem)
            call getlocalRecElement(rec,i,glob2loc_elmnts(rec%recelem(i)))
        end do
        ! create local source var for proc that carries the source
        j=1
        write(*,*)
        write(*,*) "start searching receivers in partitions"
        do iproc=1,par%nproc
            if (db(iproc)%has_rec) then
                ! recreate rec arrays
                allocate(locrecxyz(3,tempv(iproc)))
                allocate(locrecnr(tempv(iproc)))
                allocate(locrecnum(tempv(iproc)))
                locnrec=tempv(iproc)
                do i=1,tempv(iproc)
                    locrecnum(i) = recnum(iproc,i)
                    locrecxyz(:,i) = rec%recxyz(:,locrecnum(i))
                    locrecnr(i) = rec%recnr(locrecnum(i))
                    j=j+1
                end do
                call prepareRecalcRec(temprec,locnrec,locrecxyz,locrecnr)
                call findReceiver(par,temprec,db(iproc)%vx,db(iproc)%vy,db(iproc)%vz,db(iproc)%nglob,db(iproc)%nelem,db(iproc)%ibool,&
                    db(iproc)%coord,db(iproc)%elem,db(iproc)%dr,db(iproc)%ds,db(iproc)%dt)
                write(filename,"('/recVar',i6.6)") iproc
                filename=trim(outpath)//trim(filename)
                call writeRecVar(temprec,filename)
                deallocate(locrecnum)
                deallocate(locrecxyz)
                deallocate(locrecnr)
                call deallocRecVar(temprec)
            end if
        end do
        deallocate(tempv)

        ! write mesh partition
        ! write meshfiles to vtk
        if (par%vtk) then
            do iproc=1,par%nproc
                write(filename,"('/mesh',i6.6,'.vtk')") iproc
                filename=trim(outpath)//trim(filename)
                call writeVtkTetraMeshIntdata(filename, db(iproc)%elem, db(iproc)%coord, db(iproc)%mat)
            end do

            ! write nodefile to vtk
            do iproc=1,par%nproc
                write(filename,"('/node',i6.6,'.vtk')") iproc
                filename=trim(outpath)//trim(filename)
                call writeVtkNodes(filename, db(iproc)%vx,db(iproc)%vy,db(iproc)%vz)
            end do

            ! write src to vtk
            write(filename,"('/src.vtk')")
            filename=trim(outpath)//trim(filename)
            call writeVtkNodes(filename, src%srcxyz(1,:),src%srcxyz(2,:),src%srcxyz(3,:))

            ! write receiver to vtk
            write(filename,"('/rec.vtk')")
            filename=trim(outpath)//trim(filename)
            call writeVtkNodes(filename, rec%recxyz(1,:),rec%recxyz(2,:),rec%recxyz(3,:))

            write(filename,"('/vp_model.vtk')")
            filename=trim(outpath)//trim(filename)
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vp)

            write(filename,"('/vs_model.vtk')")
            filename=trim(outpath)//trim(filename)
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%vs)

            write(filename,"('/rho_model.vtk')")
            filename=trim(outpath)//trim(filename)
            call writeVtkTetraMeshRealdata(filename, this%elem, this%coord, this%rho)

        end if

        ! write databases
        do iproc=1,par%nproc
            write(filename,"('/meshVar',i6.6)") iproc
            filename=trim(outpath)//trim(filename)
            call writeMeshVar(db(iproc),filename)
        end do

        !"--------------------------------------------------------------------------------------"
        ! deallocate arrays
        !"--------------------------------------------------------------------------------------"

        call  deallocMeshvar(this)
        do iproc=1,par%nproc
            call deallocMeshvar(db(iproc))
        end do
        call deallocSrcVar(src)
        call deallocRecVar(rec)

        deallocate(elem_free)
        deallocate(elem_absorb)
        deallocate(elmnts)
        deallocate(xadj)
        deallocate(adjncy)
        deallocate(part)
        deallocate(glob2loc_elmnts)
        deallocate(num_loc)
        deallocate(nnodes_elmnts)
        deallocate(nodes_elmnts)
        deallocate(glob2loc_nodes_nparts)
        deallocate(part_nodes)
        deallocate(num_parts)
        deallocate(glob2loc_nodes_parts)
        deallocate(glob2loc_nodes)
        deallocate(glob2loc_nodes2)
        deallocate(icom)
        deallocate(epsil)

    end subroutine createMesh




    !"--------------------------------------------------------------------------------------"
    ! write meshVar to file
    !
    subroutine writeMeshVar(this,filename)
        implicit none
        type(meshVar) :: this
        character(len=80) filename

        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        write(27) this%nglob
        write(27) this%nelem
        write(27) this%nmat
        write(27) this%ncoord
        write(27) this%pinterfaces
        write(27) this%mpi_nn
        write(27) this%mpi_nnmax
        write(27) this%mpi_ne
        write(27) this%mpi_nemax
        write(27) this%nproc
        write(27) this%has_src
        write(27) this%has_rec
        write(27) this%dtfactor
        write(27) this%mingll_jac
        write(27) this%invVdm
        write(27) this%Vdm
        write(27) this%Dr
        write(27) this%Ds
        write(27) this%Dt
        write(27) this%Drw
        write(27) this%Dsw
        write(27) this%Dtw
        write(27) this%Vdmr
        write(27) this%Vdms
        write(27) this%Vdmt
        write(27) this%mass
        write(27) this%lift
        write(27) this%fmask
        write(27) this%fmaskv

        write(27) this%vx
        write(27) this%vy
        write(27) this%vz
        write(27) this%matval
        write(27) this%rx
        write(27) this%ry
        write(27) this%rz
        write(27) this%sx
        write(27) this%sy
        write(27) this%sz
        write(27) this%tx
        write(27) this%ty
        write(27) this%tz
        write(27) this%jacobian
        write(27) this%mat
        write(27) this%ibool
        write(27) this%nx
        write(27) this%ny
        write(27) this%nz
        write(27) this%sj
        write(27) this%fscale
        write(27) this%nnormal
        write(27) this%stangential
        write(27) this%ttangential
        write(27) this%coord
        write(27) this%elem
        write(27) this%neighbor
        write(27) this%face
        write(27) this%vp
        write(27) this%vs
        write(27) this%rho
        write(27) this%mu
        write(27) this%lambda
      
        write(27) this%qp
        write(27) this%qs
        write(27) this%vpu
        write(27) this%vsu
        write(27) this%muu
        write(27) this%lambdau
        write(27) this%ylambda
        write(27) this%ymu
        write(27) this%wl

        write(27) this%pml
        write(27) this%pmlxmin
        write(27) this%pmlxmax
        write(27) this%pmlymin
        write(27) this%pmlymax
        write(27) this%pmlzmin
        write(27) this%pmlzmax
    
        write(27) this%ibt
        write(27) this%ibn
        write(27) this%loc2glob_nodes
        write(27) this%mpi_interface
        write(27) this%mpi_neighbor
        write(27) this%mpi_connection
        write(27) this%mpi_ninterface
        write(27) this%mpi_ibool
        write(27) this%mpi_ibt
        write(27) this%mpi_icon

        write(27) this%pml_delta

        !    write(27) this%vol
        close(27)
    end subroutine writeMeshVar
    !"--------------------------------------------------------------------------------------"
    ! read meshVar from file
    !
    subroutine readMeshVar(this,filename)
        type(meshVar) :: this
        character(len=80) :: filename

        open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
        read(27) this%nglob
        read(27) this%nelem
        read(27) this%nmat
        read(27) this%ncoord
        read(27) this%pinterfaces
        read(27) this%mpi_nn
        read(27) this%mpi_nnmax
        read(27) this%mpi_ne
        read(27) this%mpi_nemax
        read(27) this%nproc
        read(27) this%has_src
        read(27) this%has_rec
        read(27) this%dtfactor
        read(27) this%mingll_jac
        read(27) this%invVdm
        read(27) this%Vdm
        read(27) this%Dr
        read(27) this%Ds
        read(27) this%Dt
        read(27) this%Drw
        read(27) this%Dsw
        read(27) this%Dtw
        read(27) this%Vdmr
        read(27) this%Vdms
        read(27) this%Vdmt
        read(27) this%mass
        read(27) this%lift
        read(27) this%fmask
        read(27) this%fmaskv

        allocate(this%vx(this%nglob),this%vy(this%nglob),this%vz(this%nglob))
        allocate(this%matval(this%nmat,11))
        allocate(this%rx(this%nglob), this%ry(this%nglob), this%rz(this%nglob))
        allocate(this%sx(this%nglob), this%sy(this%nglob), this%sz(this%nglob))
        allocate(this%tx(this%nglob), this%ty(this%nglob), this%tz(this%nglob))
        allocate(this%jacobian(this%nglob))
        allocate(this%mat(this%nelem))
        allocate(this%ibool(Np,this%nelem))
        allocate(this%nx(4*NpF,this%nelem),this%ny(4*NpF,this%nelem), this%nz(4*NpF,this%nelem), this%sj(4*NpF,this%nelem),&
            this%fscale(4*NpF,this%nelem))
        allocate(this%nnormal(3,4,this%nelem))
        allocate(this%stangential(3,4,this%nelem))
        allocate(this%ttangential(3,4,this%nelem))
        allocate(this%coord(3,this%ncoord))
        allocate(this%elem(4,this%nelem))
        allocate(this%neighbor(4,this%nelem))
        allocate(this%face(4,this%nelem))
        allocate(this%vp(this%nelem),this%vs(this%nelem),this%rho(this%nelem),this%mu(this%nelem),this%lambda(this%nelem))
        allocate(this%vpu(this%nelem),this%vsu(this%nelem),this%muu(this%nelem),this%lambdau(this%nelem))
        allocate(this%ylambda(nMB,this%nelem),this%ymu(nMB,this%nelem), this%wl(nMB))
        allocate(this%pml(this%nelem))
        allocate(this%ibt(NpF,4,this%nelem),this%ibn(NpF,4,this%nelem))
        allocate(this%loc2glob_nodes(this%ncoord))
        allocate(this%mpi_interface(4,4,this%nelem))
        allocate(this%mpi_neighbor(this%mpi_nn))
        allocate(this%mpi_connection(this%mpi_nn,this%mpi_ne,2))
        allocate(this%mpi_ninterface(this%mpi_nn))
        allocate(this%mpi_ibool(4,this%mpi_nemax))
        allocate(this%mpi_ibt(NpF,4,this%nelem))
        allocate(this%mpi_icon(this%mpi_nnmax))
        !    allocate(this%vol(this%nelem))
        read(27) this%vx
        read(27) this%vy
        read(27) this%vz
        read(27) this%matval
        read(27) this%rx
        read(27) this%ry
        read(27) this%rz
        read(27) this%sx
        read(27) this%sy
        read(27) this%sz
        read(27) this%tx
        read(27) this%ty
        read(27) this%tz
        read(27) this%jacobian
        read(27) this%mat
        read(27) this%ibool
        read(27) this%nx
        read(27) this%ny
        read(27) this%nz
        read(27) this%sj
        read(27) this%fscale
        read(27) this%nnormal
        read(27) this%stangential
        read(27) this%ttangential
        read(27) this%coord
        read(27) this%elem
        read(27) this%neighbor
        read(27) this%face
        read(27) this%vp
        read(27) this%vs
        read(27) this%rho
        read(27) this%mu
        read(27) this%lambda

        read(27) this%qp
        read(27) this%qs
        read(27) this%vpu
        read(27) this%vsu
        read(27) this%muu
        read(27) this%lambdau
        read(27) this%ylambda
        read(27) this%ymu
        read(27) this%wl

        read(27) this%pml
        read(27) this%pmlxmin
        read(27) this%pmlxmax
        read(27) this%pmlymin
        read(27) this%pmlymax
        read(27) this%pmlzmin
        read(27) this%pmlzmax
    
        read(27) this%ibt
        read(27) this%ibn
        read(27) this%loc2glob_nodes
        read(27) this%mpi_interface
        read(27) this%mpi_neighbor
        read(27) this%mpi_connection
        read(27) this%mpi_ninterface
        read(27) this%mpi_ibool
        read(27) this%mpi_ibt
        read(27) this%mpi_icon

        read(27) this%pml_delta
        !    read(27) this%vol
        close(27)
    end subroutine readMeshVar

    subroutine readmesh(myrank,mesh)
        integer:: myrank
        logical:: file_exists
        type(meshVar) :: mesh
        character(len=256) ::filename

        write(filename,"('/meshVar',i6.6)") myrank+1
        filename=trim(outpath)//trim(filename)
        inquire(file=trim(filename), exist=file_exists)
        if (file_exists) then
            call readMeshVar(mesh,filename)
        else
            write(*,*) "error in databases, files not existing"
            call stop_mpi()
        end if

     end subroutine readmesh
end module meshMod
