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

module fileParameterMod
	use constantsMod
   implicit none

 contains

    function getParameterName(line) result(par)
        ! INPUT: _line_: line of text in form "parameter = value #comment"
        ! DOES: finds the name of a parameter at the beginning of a line
        ! RETURNS: _par_: parameter name
        character (len=*), intent(in) :: line
        character (len=len(line)) :: tmp
        character (len=50) :: par
        
        integer, parameter :: tab = 9, space = 32, eq = 61
        integer :: i, ASCII_char, n
        
        n = 0
        par = " "
        tmp = line(1:index(line, "="))
        do i=1, len(tmp)
            ASCII_char = iachar(tmp(i:i))
            if (ASCII_char /= tab .and. ASCII_char /= space .and. ASCII_char /= eq) then
                n = n + 1
                par(n:n) = tmp(i:i)
            end if
        end do
    end function getParameterName
    
    function getParameterValue(name, filename, pos) result(value)
        ! INPUT: _name_:parameter to search for, _filename_: filename to search in, _pos_: start position to start search
        ! DOES: finds value of paramater _name_ in file _filename_
        ! RETURNS: _value_: value of parameter _name_
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: pos
        
        character (len=80) :: tmp, par, value 
        character (len=255) :: line
        
        integer :: j, k, ier = 0
        integer, parameter :: seek_set = 0
        
        integer :: set_pos
        do while (.not. is_iostat_end(ier))     ! go through file and read line for line
            read(19, "(a50)", iostat=ier) line
            if (is_iostat_end(ier)) then        ! error message if end of file reached
                write (*,*) "Parameter not found in parameter file...", name
                stop
            endif
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle         ! cycle loop if line is a comment
            if (line(1:1) == "\n") cycle        ! cycle loop if line is emtpy
            par = getParameterName(line)
            par = trim(par)                     ! get paramter name of line
            j = scan(line, "=")                 ! get position of "="
            if (par /= name) then               ! test if searched parameter is found
            	cycle
            else
                do k = j, len_trim(line)
                    if (line(k:k) /= "=" .and. line(k:k) /= " ") then
                        tmp = trim(line(k:len_trim(line)))
                        read(tmp, *) value      ! read value if right one is found
                        value = trim(value)
                        exit
                    endif
                enddo
                exit
            endif
        enddo
        call fseek(19, pos, seek_set, ier)
        set_pos = ftell(19) !warum diese zeile hier hin muss verstehe ich nicht aber nur wenn sie da steht funktioniert die funktion. Es kann durchaus sein, dass dies ein Bug vom comlier bzgl. der funktion fseek ist
    end function getParameterValue
    
    subroutine readStringPar(par_to_read, name, filename, pos)
        !Function to read Parameters as text. (will probably not be used so much...)
        character (len=*) :: name, par_to_read, filename
        character (len=10) :: value!, tmp
        integer :: pos
        
        value = getParameterValue(name, filename, pos)  
        read(value, *) par_to_read 
    end subroutine readStringPar
    
    subroutine readIntPar(par_to_read, name, filename, pos)
        !function to read interger values from the parameter-file
        character (len=*), intent(in) :: name, filename
        integer :: par_to_read, pos
        character (len=10) :: value
        
        value = getParameterValue(name, filename, pos)
        read(value, *) par_to_read          
    end subroutine readIntPar
    
    subroutine readFloatPar(par_to_read, name, filename, pos)
        !Function to read floating point variables from the parameter-file
        character (len=*), intent(in) :: name, filename
        real(kind=CUSTOM_REAL) :: par_to_read
        character (len=10) :: value
        integer :: pos
        
        value = getParameterValue(name, filename, pos)
        read(value, *) par_to_read                  
    end subroutine readFloatPar
    
    subroutine readLogicalPar(par_to_read, name, filename, pos)
        !Function to read logical variables from the parameter-file
        character (len=*), intent(in) :: name, filename
        logical :: par_to_read
        character (len=10) :: value
        integer :: pos
        
        value = getParameterValue(name, filename, pos)  
        read(value, *) par_to_read              
    end subroutine readLogicalPar
    
    function setFilePosition(name, filename, nr) result(seek_pos)
        ! INPUT: _name_: variable name to search position for, _filename_: name of file to search in, _nr_: number of source
        ! DOES: searches start position to search for parameters of each source
        ! RETURNS: position of _name_ in _filename_
        character (len=*), intent(in) :: name, filename
        integer, intent(in) :: nr
        
        character (len=80) :: tmp, par, value 
        character (len=255) :: line
        integer :: seek_pos, ier, tmp_int, j, k
    
        integer, parameter :: seek_set = 0
            
        do          
            read(19, "(a50)", iostat=ier) line  ! go through file and read line for line
            if (is_iostat_end(ier)) then        ! error message if end of file reached
                write (*,*) "Parameter ",name," not found in file..."
                exit
            endif
            line = trim(adjustl(line))
            if (line(1:1) == "#") cycle         ! cycle loop if line is a comment
            if (line(1:1) == "\n") cycle        ! cycle loop if line is empty
            par = getParameterName(line)        ! get name of parameter in line
            par = trim(par)
            j = scan(line, "=")                 ! find position of "="
            if (par /= name) cycle              ! test if searched parameter is found
            if (par == name) then
                do k = j, len_trim(line)
                    if (line(k:k) /= "=" .and. line(k:k) /= " ") then
                        tmp = trim(line(k:len_trim(line)))
                        read(tmp, *) value      ! read rest of line
                        value = trim(value)
                        read(value,*) tmp_int   ! get value of parameter
                        exit
                    endif
                enddo
            endif
            if (tmp_int == nr) then             ! if searched number of source is found
                seek_pos = ftell(19)
                exit                            !Exit the loop if a matching source has been found
            endif
        enddo
    end function setFilePosition

end module fileParameterMod
