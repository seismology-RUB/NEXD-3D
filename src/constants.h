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

! -*- f90 -*-
  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

  integer, parameter :: CUSTOM_REAL = SIZE_REAL
  integer, parameter :: ORDER = 4
  integer, parameter :: NGLL = ORDER + 1 !5
  integer, parameter :: Np = (ORDER + 1) * (ORDER + 2) *(ORDER + 3) / 6
  integer, parameter :: NpF = (ORDER + 1) * (ORDER + 2) / 2
  real(kind=CUSTOM_REAL), parameter :: PI = 3.141592653589793 
  real(kind=CUSTOM_REAL), parameter :: EPS = 1.0e-6
  real(kind=CUSTOM_REAL), dimension(15), parameter :: balpha=(/0.0000, 0.0000, 0.0000, 0.1002,&
       & 1.1332, 1.5608, 1.3413, 1.2577, 1.1603, 1.10153, 0.6080, 0.4523, 0.8856, 0.8717, 0.9655/)

! attenuation
  integer, parameter :: nMB = 3 ! Number of SLS

  integer, parameter :: max_neighbor = 8  ! maximum number of neighbors per element
  integer, parameter :: nsize = 40 ! elements share one node
  character(len=256), parameter :: outpath = "out"

