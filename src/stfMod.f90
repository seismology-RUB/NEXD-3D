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
module stfMod
  use constantsMod
  implicit none

contains
  function stfRicker(t,f0,t0,factor)
    ! INPUT: t: time, f0: central frequency of stf, t0: time shift, factor: amplitude
    ! DOES: calculates Ricker wavelet as source time function
    ! RETURNS: source time function at time t
    real(kind=CUSTOM_REAL) :: stfRicker
    real(kind=CUSTOM_REAL) :: f0,t0,factor,t
    real(kind=CUSTOM_REAL) :: aval
    aval = pi*pi*f0*f0
    stfRicker = - factor * (1.-2.*aval*(t-t0)**2.) * exp(-aval*(t-t0)**2.)
  end function stfRicker

  function stfGauss(t,f0,t0,factor)
    ! INPUT: t: time, f0: central frequency of stf, t0: time shift, factor: amplitude
    ! DOES: calculates Gaussian wavelet as source time function
    ! RETURNS: source time function at time t
    real(kind=CUSTOM_REAL) :: stfGauss
    real(kind=CUSTOM_REAL) :: f0,t0,factor,t
    real(kind=CUSTOM_REAL) :: aval
    aval = pi*pi*f0*f0
    stfGauss = -factor * exp(-aval*(t-t0)**2.)! * sqrt(pi) * f0
  end function stfGauss

  subroutine stfReadFile(t,stfFile)
    ! INPUT: t: time
    ! DOES: calculates Gaussian wavelet as source time function
    ! RETURNS: source time function at time t
    real(kind=CUSTOM_REAL) :: dummy
    real(kind=CUSTOM_REAL), dimension(:) :: t, stfFile
    real(kind=SIZE_DOUBLE), dimension(:), allocatable :: stfFile_tmp
    real(kind=SIZE_DOUBLE), dimension(:), allocatable :: time_old, stf_old
    integer :: j, ier, nt_old

    allocate(stfFile_tmp(size(stfFile)))
    open(unit=27,file='data/stf',status='old')
    j=0
    ier=0
    do while (.not. is_iostat_end(ier))
        read(27,*,iostat=ier) dummy, dummy
        if (is_iostat_end(ier)) exit
        j = j + 1
    enddo
    close(27)
    nt_old=j
    allocate(time_old(nt_old), stf_old(nt_old))
    open(unit=27,file='data/stf',status='old')
    do j=1,nt_old
        read(27,*) time_old(j), stf_old(j)
        time_old(j)=time_old(j)*1000
    enddo
    call cubicInterpolation(time_old, stf_old, dble(t), stfFile_tmp)
    stfFile=real(stfFile_tmp)
    close(27)
    deallocate(time_old,stf_old,stfFile_tmp)

  end subroutine stfReadFile

  subroutine cubicInterpolation(t_old, y_old, t_new, y_new)
    real(kind=8), dimension(:) :: t_old, t_new
    real(kind=8), dimension(:) :: y_old, y_new
    real(kind=8) ::a,b,c,d,y12,t23,t34,t24,y23,t12,t13,y34,t1t1,t1t2,t2t3,t3t3,t2t2,t3t4,t4t4, t1pt2
    integer :: i, ind, length_old, length_new

    length_old=size(t_old)
    length_new=size(t_new)
    do i=1,length_new
        if (i <= length_old) then
            if (abs(t_new(i)-t_old(i)) < 1.E-6) then
                y_new(i)=y_old(i)
            else
                call findindex(real(t_old,kind=4),real(t_new(i),kind=4),ind)
                if (ind==1) then
                    ind=2
                else if (length_old-ind==1) then
                    ind=length_old-2
                else if (length_old-ind==0) then
                    ind=length_old-2
                endif
                y12=y_old(ind-1)-y_old(ind)
                t23=t_old(ind)-t_old(ind+1)
                t34=t_old(ind+1)-t_old(ind+2)
                t24=t_old(ind)-t_old(ind+2)
                y23=y_old(ind)-y_old(ind+1)
                t12=t_old(ind-1)-t_old(ind)
                t13=t_old(ind-1)-t_old(ind+1)
                y34=y_old(ind+1)-y_old(ind+2)
                t1t1=t_old(ind-1)*t_old(ind-1)
                t1t2=t_old(ind-1)*t_old(ind)
                t2t3=t_old(ind)*t_old(ind+1)
                t3t3=t_old(ind+1)*t_old(ind+1)
                t2t2=t_old(ind)*t_old(ind)
                t3t4=t_old(ind+1)*t_old(ind+2)
                t4t4=t_old(ind+2)*t_old(ind+2)
                t1pt2=t_old(ind-1)+t_old(ind)
                a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
                  ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
                b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
                c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
                d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
                y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
            endif
        else
            call findindex(real(t_old,kind=4),real(t_new(i),kind=4),ind)
            if (ind==1) then
                ind=2
            else if (length_old-ind==1) then
                ind=length_old-2
            else if (length_old-ind==0) then
                ind=length_old-2
            endif
            y12=y_old(ind-1)-y_old(ind)
            t23=t_old(ind)-t_old(ind+1)
            t34=t_old(ind+1)-t_old(ind+2)
            t24=t_old(ind)-t_old(ind+2)
            y23=y_old(ind)-y_old(ind+1)
            t12=t_old(ind-1)-t_old(ind)
            t13=t_old(ind-1)-t_old(ind+1)
            y34=y_old(ind+1)-y_old(ind+2)
            t1t1=t_old(ind-1)*t_old(ind-1)
            t1t2=t_old(ind-1)*t_old(ind)
            t2t3=t_old(ind)*t_old(ind+1)
            t3t3=t_old(ind+1)*t_old(ind+1)
            t2t2=t_old(ind)*t_old(ind)
            t3t4=t_old(ind+1)*t_old(ind+2)
            t4t4=t_old(ind+2)*t_old(ind+2)
            t1pt2=t_old(ind-1)+t_old(ind)
            a=(y12*t23*t34*t24-y23*t12*t34*(t24+t13)+y34*t23*t13*t12)/(t12*t23*t34)/&
              ((t1t1+t1t2-t2t3-t3t3)*t24-(t2t2+t2t3-t3t4-t4t4)*t13)
            b=(y23*t34-y34*t23)/(t23*t34*t24)-a*(t2t2+t2t3-t3t4-t4t4)/t24
            c=y12/t12-a*(t1t1+t1t2+t2t2)-b*t1pt2
            d=y_old(ind)-a*t_old(ind)**3-b*t_old(ind)**2-c*t_old(ind)
            y_new(i)=a*t_new(i)**3+b*t_new(i)**2+c*t_new(i)+d
        endif
    enddo
  end subroutine cubicInterpolation

  subroutine findindex(t_old,t_new,ind)
    real(kind=4), dimension(:) :: t_old
    real(kind=4) :: t_new
    integer :: ind, left, right, middle, size_t, i

    ind=-1
    left=1
    size_t=size(t_old)
    right=size_t
    if (t_new >= t_old(right)) then
        left=right
        ind=right
        if (t_new > 2*t_old(right)-t_old(right-1)) then
            write(*,*) 't_new value outside of t_old range', t_new, t_old(right)
            stop
        endif
    endif
    do while (right-left > 5)
        middle=(right+left)/2

        if (t_old(middle) > t_new) then
            right=middle
        else if (t_old(middle) < t_new) then
            left=middle
        else
            ind=middle
            exit
        endif
    enddo
    do i=left,right-1
        if (t_old(i)<=t_new .and. t_old(i+1)>t_new) ind=i
    enddo
    if (ind==-1) then
        write(*,*) 'index not found', t_new, t_old(left), t_old(right)
        stop
    endif
    if (ind < size_t) then
        if (.not. t_new>=t_old(ind).and.t_new<t_old(ind+1)) then
            write(*,*) 'Error in index', ind, t_new, t_old(ind), t_old(ind+1)
            stop
        endif
    endif
  end subroutine findindex

end module stfMod
