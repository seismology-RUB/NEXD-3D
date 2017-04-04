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
MODULE sortMod
!!! this code is from the rosetta code wiki
 
IMPLICIT NONE
 
CONTAINS
 
RECURSIVE SUBROUTINE Qsort(a)
 
  INTEGER, INTENT(IN OUT) :: a(:)
  INTEGER :: split
 
  IF(size(a) > 1) THEN
     CALL Partition(a, split)
     CALL Qsort(a(:split-1))
     CALL Qsort(a(split:))
  END IF
 
END SUBROUTINE Qsort
 
SUBROUTINE Partition(a, marker)
 
  INTEGER, INTENT(IN OUT) :: a(:)
  INTEGER, INTENT(OUT) :: marker
  INTEGER :: left, right, pivot, temp
 
  pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
  left = 0                         ! behavior with sorted or reverse sorted data
  right = size(a) + 1
 
  DO WHILE (left < right)
     right = right - 1
     DO WHILE (a(right) > pivot)
        right = right-1
     END DO
     left = left + 1
     DO WHILE (a(left) < pivot)
        left = left + 1
     END DO
     IF (left < right) THEN 
        temp = a(left)
        a(left) = a(right)
        a(right) = temp
     END IF
  END DO
 
  IF (left == right) THEN
     marker = left + 1
  ELSE
     marker = left
  END IF
 
END SUBROUTINE Partition
 
END MODULE sortMod
