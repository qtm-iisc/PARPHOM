! Package: PARPHOM
! Authors: Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
! License: GPL-3.0
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!> \file distribute_arrays.F90
!> \brief Routines for distributing arrays across processes in PARPHOM.
!>
!> This file contains routines for distributing and gathering arrays needed for parallel phonon calculations.
!>
!> - Handles data distribution for force constants, dynamical matrices, and eigenvectors
!> - Ensures correct mapping between global and local indices
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \warning
!>   Ensure that all arrays are allocated and initialized before distribution.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> 
subroutine distribution_length(tot_len, iprc, nprc, start, fin)
  
  ! -----------------------------------------------
  ! Distribute an array among different processors
  ! 
  ! tot_len (input) : total length of the array
  ! iprc    (input) : process id 
  ! nprc    (input) : total number of processes
  ! start   (input) : starting index of the distributed array in iprc
  ! fin     (input) : last index of the distributed array in iprc
  ! 
  ! ------------------------------------------------

  implicit none
  integer, intent(in) :: tot_len, iprc, nprc
  integer, intent(out) :: start, fin
  integer :: rem_each, size_each
  if (tot_len.ge.nprc) then
      size_each = tot_len/nprc
      rem_each = mod(tot_len,nprc)
      if (iprc.gt.(nprc-rem_each-1)) then
          size_each = size_each+1
          start = iprc*size_each-(nprc-rem_each)+1
          fin = (iprc+1)*size_each-(nprc-rem_each)
      else
          start = iprc*size_each+1
          fin = (iprc+1)*size_each
      end if
  else
      rem_each = nprc - tot_len - 1
      if ((iprc-rem_each).lt.0) then
          start = 0
          fin = start
      else
          start = iprc-rem_each
          fin = start
      end if
  end if
  return
end subroutine
