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
module typeMod
    use constantsMod

    implicit none


    type :: meshVar        ! container for mesh variables
        sequence
        ! Mesh variables
        integer, dimension(:,:), allocatable :: elem        ! contains mesh points to build an element
        integer, dimension(:,:), allocatable :: neighbor    ! contains neighbor elements of each element for every face
        integer :: nglob, nelem, nmat, ncoord           ! number of: global points, global elements, global materials, global coordinates
        integer(kind=CUSTOM_REAL), dimension(:), allocatable :: mat ! material numbers of simulation setup
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: coord    ! coordinates of element nodes
        integer, dimension(:,:), allocatable :: face        ! face number of neighbor element, where both elements are connected
        integer, dimension(:,:), allocatable :: ibool       ! increasing number for all points in a rank (assigns a for each element-point combination a number)
        integer, dimension(:,:,:), allocatable :: ibn,ibt   ! point numbers in NEIGHBOR element at face to neighbor; point numbers in element at face to a neighbor
        integer, dimension(:), allocatable :: loc2glob_nodes    ! translates local nodes back to global numbering
        real(kind=CUSTOM_REAL) :: minGLL_jac            ! minimal GLL distance times jacobian and normals for calculation of timestep
        ! PML variables
        integer, dimension(:), allocatable :: pml       ! 1: element is located in a PML, 0 otherwise
        real(kind=CUSTOM_REAL) :: pmlxmin, pmlxmax, pmlymin, pmlymax, pmlzmin, pmlzmax  ! coordinates of setup boundaries, important for PML calculation
        ! Physical variables
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: matval       ! array for physical parameters, for definition see read in
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: vp ,vs ,rho ,mu ,lambda ,qp ,qs    ! physical parameters
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: vpu,vsu,muu,lambdau    ! physical parameters for infinite frequency in attenuation case
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: ylambda,ymu      ! attenuation parameters
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: wl,wk          ! relaxation frequencies
        real(kind=CUSTOM_REAL) :: dtfactor          ! length of a timestep
        integer :: nt                       ! Number of timesteps
        ! MPI variables
        integer :: pinterfaces                  ! number of elements at mpi interface
        integer :: mpi_nn                   ! contains to how many other ranks a rank is connected
        integer :: mpi_nnmax                    ! maximum number of mpi interfaces in all ranks
        integer :: mpi_ne                   ! maximum number of interface faces in all ranks
        integer :: mpi_nemax                    ! maximum number of elements in a rank
        integer :: nproc                    ! number of processors
        logical :: has_src=.false.              ! boolean of a rank has a source
        logical :: has_rec=.false.              ! boolean of a rank has a receiver
        integer, dimension(:,:,:), allocatable :: mpi_interface        ! 1: partition number of neighbor element, 2: element number in other partition, 3: global to local mpi numbering of interface elements, 4: face number to 3
        integer, dimension(:), allocatable :: mpi_neighbor             ! contains all neighbor ranks of a rank
        integer, dimension(:,:,:), allocatable :: mpi_connection       ! 1: local to global numbering of interface elements, 2: at which face of these elements is the interface
        integer, dimension(:), allocatable :: mpi_ninterface           ! number of interfaces of a rank
        integer, dimension(:,:), allocatable :: mpi_ibool              ! local numbering of points in elements for a rank
        integer, dimension(:,:,:), allocatable :: mpi_ibt              ! points of face to a neighbor
        integer, dimension(:), allocatable :: mpi_icon                 ! neighbor interface to rank
        !Method variables
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: invVdm,Vdm,Dr,Ds,Dt ! inverse van-der-monde (vdm) matrix, vdm matrix, differentation matrices with respect to r,s and t
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: Drw,Dsw,Dtw,Vdmr,Vdms,Vdmt ! weak matrices with respect to r,s and t, gradient of vdm matrices with respect to r,s,t
        real(kind=CUSTOM_REAL), dimension(Np,Np) :: mass                ! mass matrix
        real(kind=CUSTOM_REAL), dimension(Np,4*NpF) :: lift             ! lift matrix
        integer, dimension(Npf,4) :: fmask                              ! mask for surface nodes
        integer, dimension(4*Npf) :: fmaskv
        real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: nx ,ny ,nz ,sj , fscale     ! normal components, length of normal, scaling
        real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nnormal ,stangential  ,ttangential           ! normal components in r,s,t
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: vx ,vy  ,vz     ! x,y,z coordinates of all points in rank in simulation space
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: rx ,ry ,rz  ! derivatives from r to x, etc.
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: sx ,sy ,sz
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: tx ,ty ,tz
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: jacobian
        real(kind=CUSTOM_REAL), dimension(:), allocatable :: vol
        real(kind=CUSTOM_REAL) :: pml_delta
    end type meshVar

  type :: parameterVar
!     private
     sequence
     logical :: log         !log to screen?
     logical :: external    !read external model?
     integer :: nproc
     logical :: vtk         !save vtk files?
     logical :: movie       !save movie?
     logical :: attenuation
     real(kind=CUSTOM_REAL) :: f0_att
     real(kind=CUSTOM_REAL) :: fr_att
     integer :: frame       !set frames for movie
     integer :: timeint     !which timeintegration?
     integer :: sft         !which sft sould we take
     integer :: nt          !how many timesteps
!
     real(kind=CUSTOM_REAL) :: cfl            !which cfl value for dt
     real(kind=CUSTOM_REAL) :: f0             !center frequency of sft
     real(kind=CUSTOM_REAL) :: factor         !factor of the stf
     real(kind=CUSTOM_REAL), dimension(3) :: angle_force    !angle of the force action in the media
!
     integer :: nsrc        !number of sources
     integer :: nrec        !number of receivers

     ! PML
     logical:: set_pml
     real(kind=CUSTOM_REAL) :: pml_delta, pml_rc, pml_kmax, pml_afac ! Width of PML, Relection coefficient, kmax, factor of amax
     integer :: avg_window1, avg_window2 ! Width of windows to get average
     real(kind=CUSTOM_REAL) sta_lta_trigger ! Threshold for energy growth
!
     integer :: sourcetype !force or moment tensor
     real(kind=CUSTOM_REAL) :: xsource,ysource,zsource !location of the source
     real(kind=CUSTOM_REAL) :: dt=0. !timestep
     real(kind=CUSTOM_REAL) :: t0 !timeshift of the source
     real(kind=CUSTOM_REAL) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz !Moment tensor
  end type parameterVar

  type :: srcVar
   integer :: nsrc
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcxyz => null()
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcrst => null()
   integer, dimension(:), pointer :: srcelem => null()
   integer, dimension(:), pointer :: srci => null()
!
   integer, dimension(:), pointer :: srctype => null()
   integer, dimension(:), pointer :: srcstf => null()
   real(kind=CUSTOM_REAL), dimension(:), pointer :: srcf0 => null()
   real(kind=CUSTOM_REAL), dimension(:), pointer :: srcfactor => null()
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcangle_force => null()
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: srcM => null()
end type srcVar
!
type :: recVar
   integer :: nrec
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recxyz => null()
   real(kind=CUSTOM_REAL), dimension(:,:), pointer :: recrst => null()
   integer, dimension(:), pointer :: recelem => null()
   integer, dimension(:), pointer :: reci => null()
   integer, dimension(:), pointer :: recnr => null()
end type recVar
end module typeMod
