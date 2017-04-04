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
module parameterMod
! module to read a parameter file to set up the simulation

  use constantsMod
  use fileParameterMod
  use typeMod
  implicit none

contains

  subroutine readParfile(this,myrank)
    implicit none
    type (parameterVar) :: this

    character(len=256) :: filename
    integer :: ier
    integer :: myrank
    logical :: log = .true.

    filename=trim('data/parfile')
    open(unit=19, file=trim(filename), status='old', action='read', iostat=ier)
    if (myrank==0) then
      write (*, "(a80)") "--------------------------------------------------------------------------------"
      write (*, "(a26, a14, a12, a28)") "|                          ","Begin reading ", filename, "...                        |"
      write (*, "(a80)") "--------------------------------------------------------------------------------"
    endif



    call readLogicalPar(this%log, "log", filename, 0)
    call readLogicalPar(this%external, "external", filename, 0)
    call readIntPar(this%nproc, "nproc", filename, 0)
    call readLogicalPar(this%vtk, "vtk", filename, 0)
    call readLogicalPar(this%movie, "movie", filename, 0)
    call readIntPar(this%frame, "frame", filename, 0)
    call readIntPar(this%timeint, "timeint", filename, 0)
    call readIntPar(this%nt, "nt", filename, 0)
    call readFloatPar(this%cfl, "cfl", filename, 0)
    call readLogicalPar(this%attenuation, "attenuation", filename, 0)
    call readFloatPar(this%f0_att,"f0_att", filename, 0)
    call readFloatPar(this%fr_att,"fr_att", filename, 0)
    call readLogicalPar(this%set_pml, "set_pml", filename, 0)
    call readFloatPar(this%pml_delta,"pml_delta", filename, 0)
    call readFloatPar(this%pml_rc,"pml_rc", filename, 0)
    call readFloatPar(this%pml_kmax,"pml_kmax", filename, 0)
    call readFloatPar(this%pml_afac,"pml_afac", filename, 0)
    call readIntPar(this%avg_window1,"avg_window1", filename, 0)
    call readIntPar(this%avg_window2,"avg_window2", filename, 0)
    call readFloatPar(this%sta_lta_trigger,"sta_lta_trigger", filename, 0)
    close(19)

    log = this%log

    if (log.and.myrank==0) then

         write (*,"(a40, l10, a30)")   "|                           Create log: ", this%log, "                             |"
         write (*,"(a40, l10, a30)")   "|                   Use external model: ", this%external, "                             |"
         write (*,"(a40, i10, a30)")   "|                 Number of processors: ", this%nproc, "                             |"
         write (*,"(a40, l10, a30)")   "|                       Save vtk files: ", this%vtk, "                             |"
         write (*,"(a40, l10, a30)")   "|                           Save movie: ", this%movie, "                             |"
         write (*,"(a40, i10, a30)")   "|   Number of timesteps for movieframe: ", this%frame, "                             |"
         write (*,"(a40, i10, a30)")   "|                      Timeintegartion: ", this%timeint, "                             |"
         write (*,"(a40, i10, a30)")   "|                  Number of timesteps: ", this%nt, "                             |" 
         write (*,"(a40, f10.1, a30)") "|                            cfl value: ", this%cfl, "                             |"
         write (*,"(a40, l10, a30)")   &
              "|                          Attenuation: ", this%attenuation, "                             |"
         write (*,"(a40, f10.1, a30)") "|                                 f0_att: ", this%f0_att, "                             |"
         write (*,"(a40, f10.1, a30)") "|                                 fr_att: ", this%fr_att, "                             |"
         write (*,"(a40, l10, a30)")   "|                                 PML: ", this%set_pml, "                             |"
         write (*,"(a40, f10.1, a30)") "|                           pml_delta: ", this%pml_delta, "                             |"
         write (*,"(a40, f10.1, a30)") "|                              pml_rc: ", this%pml_rc, "                             |"
         write (*,"(a40, f10.1, a30)") "|                            pml_kmax: ", this%pml_kmax, "                             |"
         write (*,"(a40, f10.1, a30)") "|                            pml_afac: ", this%pml_afac, "                             |"
         write (*,"(a40, i10, a30)")   "|                          LTA Window: ", this%avg_window1, "                             |"
         write (*,"(a40, i10, a30)")   "|                          STA Window: ", this%avg_window2, "                             |"
         write (*,"(a40, f10.1, a30)") "|                   STA_LTA Threshold: ", this%sta_lta_trigger, "                             |"

    end if
    if (myrank==0) then
      write (*, "(a80)") "--------------------------------------------------------------------------------"
      write (*,"(a27, a13, a13, a27)") "|                          ","Done reading ", trim(filename), "                          |"
      write (*, "(a80)") "--------------------------------------------------------------------------------"
    endif
  end subroutine readParfile
   
end module parameterMod
