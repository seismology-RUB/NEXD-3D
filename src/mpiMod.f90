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
module mpiMod
!  use mpi !
	use mpiincludeMod
	use constantsMod
	use typeMod

 implicit none
! module to access the mpi routines
contains

  subroutine init_mpi()
    implicit none
    
    integer :: ierr

    call MPI_Init(ierr)
    call check_ierr(ierr,"Error in Initializing MPI")
  end subroutine init_mpi

!!!
! call MPI_COMM_SIZE to get the world
!
  subroutine comm_size(size)
    implicit none
    
    integer :: size
    integer :: ierr
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ierr)
    call check_ierr(ierr,"Error in MPI_COMM_SIZE")
  end subroutine comm_size

!!!
! call MPI_COMM_RANK to geht the Rank of the node
!
  subroutine comm_rank(rank)
    implicit none
    integer :: ierr
    integer :: rank
    
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call check_ierr(ierr,"Error in MPI_COMM_RANK")
  end subroutine comm_rank
  
!!!
! call MPI_BARRIER to sync the processes
!
  subroutine sync_mpi()
    implicit none
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call check_ierr(ierr,"Error in MPI_BARRIER")
  end subroutine sync_mpi

!!!
! call MPI_Isend to send the data (vektor)
!
  subroutine isendV_real(buf,count,dest,tag,req,type_size)
    implicit none
    
    integer :: count, dest,tag,req,type_size
    real(kind=CUSTOM_REAL), dimension(count) :: buf
    integer :: ierr

    if (type_size==SIZE_REAL) then
        call MPI_ISSEND(buf(1),count,MPI_REAL,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_ISSEND")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_ISSEND(buf(1),count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_ISSEND")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_ISSEND")
    endif

  end subroutine isendV_real
!!!
! call MPI_Irecv to receive the data (vektor)
!
  subroutine irecV_real(buf,count,dest,tag,req,type_size)
    implicit none
    
    integer :: count, dest,tag,req,type_size
    real(kind=CUSTOM_REAL), dimension(count) :: buf
    integer :: ierr

    if (type_size==SIZE_REAL) then
        call MPI_IRECV(buf(1),count,MPI_REAL,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_IREC")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_IRECV(buf(1),count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD,req,ierr)
        call check_ierr(ierr,"Error in MPI_IREC")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_IREC")
    endif

  end subroutine irecV_real
!!!
! call mpi_wait to wait for the request
!
  subroutine wait_req(req)
    implicit none
    integer :: req
    integer, dimension(MPI_STATUS_SIZE) :: req_status
    integer :: ierr
    
    call mpi_wait(req,req_status,ierr)
    call check_ierr(ierr,"Error in MPI_WAIT")
  end subroutine wait_req
!!!
! call MPI_reduce to receive the maximum of a dataset
!
  subroutine maxval_real(send, rec,type_size)
    implicit none
    
    real(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_REDUCE(send,rec,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_REDUCE(send,rec,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    endif

  end subroutine maxval_real

! call MPI_reduce to receive the maximum of a dataset from all processes
!
  subroutine maxval_real_all(send, rec,type_size)
    implicit none

    real(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real_all")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_ALLREDUCE(send,rec,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real_all")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real_all")
    endif

  end subroutine maxval_real_all

! call MPI_reduce to receive the minimum of a dataset from all processes
!
  subroutine minval_real_all(send, rec,type_size)
    implicit none

    real(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_ALLREDUCE, minval_real_all")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_ALLREDUCE(send,rec,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_ALLREDUCE, minval_real_all")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_ALLREDUCE, minval_real_all")
    endif

  end subroutine minval_real_all
!!!
! call MPI_reduce to receive the sum of a dataset
!
  subroutine sum_real(send, rec,type_size)
    implicit none
    
    real(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_REDUCE(send,rec,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_REDUCE(send,rec,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_REDUCE, maxval_real")
    endif

  end subroutine sum_real

! call MPI_reduce to receive the sum of a dataset
!
  subroutine sum_real_all(send, rec,type_size)
    implicit none

    real(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_real_all")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_ALLREDUCE(send,rec,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_real_all")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_real_all")
    endif

  end subroutine sum_real_all

  subroutine sum_real_all_DP(send, rec)
    implicit none

    real(kind=SIZE_DOUBLE) :: send, rec
    integer :: ierr

    call MPI_ALLREDUCE(send,rec,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    call check_ierr(ierr,"Error in MPI_REDUCE, sum_real_all")


  end subroutine sum_real_all_DP
!!!
! call MPI_reduce to receive the sum of a dataset
!
  subroutine sum_int(send, rec,type_size)
    implicit none
    
    integer(kind=CUSTOM_REAL) :: send, rec
    integer :: ierr,type_size

    if (type_size==SIZE_REAL) then
        call MPI_REDUCE(send,rec,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_int")
    else if (type_size==SIZE_DOUBLE) then
        call MPI_REDUCE(send,rec,1,MPI_INTEGER8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_int")
    else
        write(*,*) 'Size of CUSTOM_REAL is not supported, see constants.h file!'
        call check_ierr(ierr,"Error in MPI_REDUCE, sum_int")
    endif

  end subroutine sum_int


!!!
! call MPI_finalize to stop MPI
!
  subroutine finalize_mpi()
    implicit none
    
    integer :: ierr

    call MPI_Finalize(ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*) "Error in finalizing MPI"
    end if
  end subroutine finalize_mpi

!!!
! call MPI_about to stop all comunications
!
  subroutine stop_mpi()
    implicit none
    integer :: ierr
    integer :: errcode

    call MPI_Abort(MPI_COMM_WORLD, errcode,ierr)
    stop 'error, program stop in stop_mpi'
  end subroutine stop_mpi
!!!
!
!
  subroutine check_ierr(ierr,msg)
    implicit none
    integer :: ierr
    character(len=*) :: msg

    if (ierr /= MPI_SUCCESS) then
       write(*,*) ierr,trim(msg)
       call stop_mpi()
    end if
  end subroutine check_ierr


  subroutine SendRecVar(mesh,qm,qi,q_send,q_rec)

      type(MeshVar) :: mesh
      real(kind=CUSTOM_REAL), dimension(:,:,:,:) :: qi
      real(kind=CUSTOM_REAL), dimension(:,:) :: qm,q_send, q_rec

      integer :: i, ie, k, j, c, dest, req, reqrec

      q_send(:,:) = 0
      do i=1,mesh%mpi_nn               ! build send buffer for MPI interfaces
          do ie=1,mesh%mpi_ne          ! loop over interface elements
              do k=1,9
                  do j=1,NpF
                      if ( mesh%mpi_connection(i,ie,1) >0) then
                          q_send((ie-1)*9*NpF + (k-1)*NpF + j,i)= qm(mesh%ibt(j,mesh%mpi_connection(i,ie,2),mesh%mpi_connection(i,ie,1)),k)
                      end if
                  end do
              end do
          end do
      end do
      do i=1,mesh%mpi_nn       !send and receive (MPI)
          dest=mesh%mpi_neighbor(i)-1
          call isendV_real(q_send(:,i),(mesh%mpi_ne*9*NpF),dest,0,req,CUSTOM_REAL)
          call irecV_real(q_rec(:,i),(mesh%mpi_ne*9*NpF),dest,0,reqrec,CUSTOM_REAL)
          call wait_req(req)
          call wait_req(reqrec)
      end do
      do i=1,mesh%mpi_nn       !prepare received mpi buffers from other ranks
          c=1
          do ie=1,mesh%mpi_ne
              do k=1,9
                  do j=1,NpF
                      if ( mesh%mpi_connection(i,ie,1) >0) then
                          qi(j,k,ie,i)=q_rec(c,i)
                      end if
                      c=c+1
                  end do
              end do
          end do
      end do              ! mpi end

  end subroutine SendRecVar

end module mpiMod
