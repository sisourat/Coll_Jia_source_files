program print_mo_overlap
  implicit none
  integer :: i,j,k,l
  double precision :: h1e(mo_num,mo_num)


  h1e = mo_one_e_integrals
  open(unit=10,file='h1emo.txt')
  write(10,*)mo_num
  do i = 1, mo_num
   do j = 1, mo_num
     write(10,*)h1e(j,i)
   enddo
  enddo
  close(10)

end


