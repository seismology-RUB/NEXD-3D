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
module logoMod
contains
  subroutine writeLogo()
    write(*,*)
    write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    write(*,*) ':::;::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    write(*,*) ':::::                                                  :::::'
    write(*,*) ':::::   Program to solve the (an)elastic wave equation :::::'
    write(*,*) ':::::   with the Nodal Discontinous Galerkin Method    :::::'
    write(*,*) ':::::   written by                                     :::::'
    write(*,*) ':::::                                                  :::::'
    write(*,*) ':::::   Lasse Lambrecht                                :;:::'
    write(*,*) ':::::   Andre Lamert                                   :;:::'
    write(*,*) ':::::   Ruhr-Universitaet Bochum 2011 - 2017           :::::'
    write(*,*) ':::::                                                  :::::'
    write(*,*) ':::::::;::::::::::::::::::::::::::::::::::::::::::::::::::::'
    write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    write(*,*) '::::::::::::;::;::;::;::;::;::;::;::;::;::;::;::;::;::;::;::'
    write(*,*) '::::;::::::::wwwwwwwwwas::::wwwc::::=vwwa::awaauwwwa/;::::::'
    write(*,*) ':::::::;::;::QQQPT?T?$QQm;;:QQQf::::=jQQZ::##::::.-)Qs::::::'
    write(*,*) ':::::::::::::QQQf::.::QWQL;:QQQf::::=3QQ#::ZW::::::<W"::::::'
    write(*,*) ':::;:::::::::QQQL=isawWQQ(::QQQf::::=3QQZ::#QowwwmQQ>:::::::'
    write(*,*) ':::::::::::::QQQf;9WWQ@?^:::QQQL::::<dQQZ::Z#:::-::+$g;:::::'
    write(*,*) ':::::::;::;::QQQf::)BQWmc:::4QQQas<ayQQ@+::#W:::::::mQ:;::::'
    write(*,*) ':::::::::::::QQQf;:::)QQWw::;?VQQQQQW@?~;::UmwawawwU?:::::::'
    write(*,*) '::::;::::::::::::::::::::.:::::::::::::;:::::.-..-::::::::::'
    write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
    write(*,*)
  end subroutine writeLogo
end module logoMod
