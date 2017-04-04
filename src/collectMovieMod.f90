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
module collectMovieMod
	use parameterMod
	use meshMod
	use logoMod
	use vtkMod

  implicit none

contains

  subroutine collectMovie(par)
    type(parameterVar) :: par
    type(meshVar),dimension(:), allocatable :: db
    integer :: iproc
    integer :: it,nframe,c,d,e,i,j
    integer, dimension(:), allocatable :: ndata,ndata2,ndata3
    ! file
    logical :: file_exists
    character(len=80) ::filename
    ! data
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: uplot
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: xyzplot
    real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tempdata
    integer :: glob_nelem, glob_nglob, glob_ncoord
    real(kind=CUSTOM_REAL) , dimension(:,:), allocatable :: glob_coord, glob_coord2
    integer , dimension(:,:), allocatable :: glob_elem


    if (par%log) call writeLogo()
    if (par%log) write(*,*) "--------------------------------------------------------------------------------------"
    if (par%log) write(*,*) "generating movie"
    if (par%log) write(*,*) "--------------------------------------------------------------------------------------"

    ! load database
    
    nframe=par%nt/par%frame

    allocate(db(par%nproc))
    allocate(ndata(par%nproc))
    allocate(ndata2(par%nproc))
    allocate(ndata3(par%nproc))

    do iproc=1,par%nproc
       write(filename,"('/meshVar',i6.6)") iproc
       filename=trim(outpath)//trim(filename)
       inquire(file=trim(filename), exist=file_exists)
       if (file_exists) then
          call readMeshVar(db(iproc),filename)
       else
          write(*,*) "error in databases, files not existing"
          stop
       end if
    end do

    glob_nelem=0
    glob_nglob=0
    glob_ncoord=0
    do iproc=1,par%nproc
       glob_nelem=glob_nelem+db(iproc)%nelem
       glob_nglob=glob_nglob+db(iproc)%nglob
       glob_ncoord=glob_ncoord+db(iproc)%ncoord
       ndata(iproc) = db(iproc)%nglob
       ndata2(iproc) = db(iproc)%ncoord
       ndata3(iproc) = db(iproc)%nelem
    end do

    allocate(uplot(glob_nglob,nframe))
    allocate(tempdata(glob_nelem,nframe))
    allocate(xyzplot(3,glob_nglob))

    ! read data
    d=1
    e=ndata(1)
    do iproc=1,par%nproc
       c=1
       do it=par%frame,par%nt,par%frame
          write(filename,"('/moviedata',i6.6,'_it',i7.7,'.bin')") iproc,it
          filename=trim(outpath)//trim(filename)
          open(unit=27,file=trim(filename),form = "UNFORMATTED",status='unknown')
          read(27) uplot( d:e ,c )
          close(27)
          c=c+1
       end do
       if (iproc < par%nproc) then
          d=d+ndata(iproc)
          e=e+ndata(iproc+1)
       end if
    end do

    ! read coords
    d=1
    e=ndata(1)
    do iproc=1,par%nproc
       xyzplot(1,d:e) = db(iproc)%vx(:)
       xyzplot(2,d:e) = db(iproc)%vy(:)
       xyzplot(3,d:e) = db(iproc)%vz(:)
       if (iproc < par%nproc) then
          d=d+ndata(iproc)
          e=e+ndata(iproc+1)
       end if
    end do


    c=0
    d=0
    tempdata(:,:) = 0
    do iproc=1,par%nproc
       do it=1,nframe
          do i=1,db(iproc)%nelem
             do j=1,Np
                tempdata(i+c,it)=tempdata(i+c,it)+uplot(db(iproc)%ibool(j,i)+d,it)
             end do
          end do
       end do
       if (iproc < par%nproc) then
          c=c+ndata3(iproc)
          d=d+ndata(iproc)
       end if
    end do
    tempdata=tempdata/Np

    !build global mesh
    allocate(glob_elem(4,glob_nelem))
    allocate(glob_coord(3,glob_ncoord))

    c=1
    do iproc=1,par%nproc
       do i=1,db(iproc)%nelem
          glob_elem(1,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(1,i))
          glob_elem(2,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(2,i))
          glob_elem(3,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(3,i))
          glob_elem(4,c)=db(iproc)%loc2glob_nodes(db(iproc)%elem(4,i))
          c=c+1
       end do
    end do
    
    glob_coord=-1
    do iproc=1,par%nproc
       do i=1,db(iproc)%ncoord
          glob_coord(1,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(1,i)
          glob_coord(2,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(2,i)
          glob_coord(3,db(iproc)%loc2glob_nodes(i)) = db(iproc)%coord(3,i)
       end do
    end do
    
    c=0
    do i=1,glob_ncoord
       if (glob_coord(1,i)>-1) then
          c=c+1
       end if
    end do
    allocate(glob_coord2(3,c))    
    c=0
    do i=1,glob_ncoord
       if (glob_coord(1,i)>-1) then
          c=c+1
          glob_coord2(:,c)=glob_coord(:,i)
       end if
    end do

    if(par%log) write(*,*) "There are ",glob_nelem, "elements and ",c,"coordinates and ",glob_nglob," points and ",nframe," frames"
    do it=1,nframe
       write(filename,"('/movie_elem_it',i7.7,'.vtk')") it
       filename=trim(outpath)//trim(filename)
       call writeVtkTetraMeshRealdata(filename, glob_elem, glob_coord2, tempdata(:,it))
    end do

    ! It is possible to write a 3D VTK with the wave field at all GLL points. But be careful. This needs a lot of time and storage!
    !do it=1,nframe
    !   write(filename,"('/movie__points_it',i7.7,'.vtk')") it
    !   filename=trim(outpath)//trim(filename)
    !   call writeVtkNodesRealData(filename, xyzplot(1,:),xyzplot(2,:),xyzplot(3,:), uplot(:,it))
    !end do

    deallocate(uplot)
    deallocate(db)
    deallocate(ndata)
    deallocate(ndata2)
    deallocate(ndata3)
    deallocate(tempdata)
    deallocate(xyzplot)
    deallocate(glob_elem)
    deallocate(glob_coord)

  end subroutine collectMovie
end module collectMovieMod
