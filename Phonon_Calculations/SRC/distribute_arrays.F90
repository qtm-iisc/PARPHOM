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
