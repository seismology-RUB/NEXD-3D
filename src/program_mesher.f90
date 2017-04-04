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
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::          
!:::;::::::::::::::::::::::::::::::::::::::::::::::::::::::::          
!:::::                                                  :::::
!:::::   Program to solve the elastic wave equation     :::::
!:::::   with the Nodal Discontinous Galerkin Method    :::::
!:::::   written by                                     :::::
!:::::                                                  :::::
!:::::   Lasse Lambrecht                                :;:::
!:::::   Ruhr-Universitaet Bochum                       :::::
!:::::                                                  :::::
!:::::::;::::::::::::::::::::::::::::::::::::::::::::::::::::          
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::          
!::::::::::::;::;::;::;::;::;::;::;::;::;::;::;::;::;::;::;::          
!::::;::::::::wwwwwwwwwas::::wwwc::::=vwwa::awaauwwwa/;::::::          
!:::::::;::;::QQQPT?T?$QQm;;:QQQf::::=jQQZ::##::::.-)Qs::::::          
!:::::::::::::QQQf::.::QWQL;:QQQf::::=3QQ#::ZW::::::<W"::::::          
!:::;:::::::::QQQL=isawWQQ(::QQQf::::=3QQZ::#QowwwmQQ>:::::::          
!:::::::::::::QQQf;9WWQ@?^:::QQQL::::<dQQZ::Z#:::-::+$g;:::::          
!:::::::;::;::QQQf::)BQWmc:::4QQQas<ayQQ@+::#W:::::::mQ:;::::          
!:::::::::::::QQQf;:::)QQWw::;?VQQQQQW@?~;::UmwawawwU?:::::::          
!::::;::::::::::::::::::::.:::::::::::::;:::::.-..-::::::::::          
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::          

program mesher
!program to mesh the given cubit files to build a MPI version of the code
	use typeMod
	use parameterMod
	use meshMod
 implicit none

  type (parameterVar) :: par            ! contains all parameters from parameter file (parfile)

  call readParfile(par,0)               ! Read information from parfile (0 is here just a dummy)                     ! get log from container par
  if (par%log) then
     write(*,*) "--------------------------------------------------------------------------------------"
     Write(*,*) "start mesher for mpi version of dg3d"
     write(*,*) "--------------------------------------------------------------------------------------"
  end if
  call createMesh(par)                  ! create the mesh and write it to database files

end program mesher
