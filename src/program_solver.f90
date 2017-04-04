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

program solver
	use parameterMod
	use timeloopMod
	use collectMovieMod
	use typeMod

    type(parameterVar) :: par
    integer :: world_size

    integer :: myrank

    ! program routine to solve the elastic wave equation with the nodal DG method

    call init_mpi()
    call comm_size(world_size)
    call comm_rank(myrank)

    if (myrank==0) then
        call writeLogo()
    end if
    call sync_mpi()
    call readParfile(par, myrank)
    if (par%log.and.myrank==0) then
        write(*,*) "-------------------------------------------------------------------------------"
        write(*,*) "start solver for mpi version of dg3d"
        write(*,*) "-------------------------------------------------------------------------------"
    end if
    call sync_mpi()
    call timeloop3d(par, myrank)
    call finalize_mpi()
end program solver
