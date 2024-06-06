program cippres_two_e_int
  use general
  use map_module
 ! create a routine for two-e int mo_two_e_integral(l,k,j,i)

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_prop_collision solves the TDSE using the matrix elements from cippres_collision
  END_DOC

  logical :: exists
  integer :: i, j, ib, k, l, it
  !get_two_e_integral(i,j,k,l,mo_integrals_map)
  double precision :: twoeint
  double precision :: mo_two_e_integral,get_two_e_integral
  !double precision :: mo_two_e_integral(mo_num,mo_num,mo_num,mo_num)
  !double precision :: mo_one_e_integrals(mo_num,mo_num), mo_overlap(mo_num,mo_num) 
  double precision :: tmp, temp
  

  temp = big_array_coulomb_integrals(mo_num,mo_num,mo_num)

  !PROVIDE ezfio_filename !HF_bitmask mo_coef

  write(*,*) "mo_two_e_integral : prop 0"
  open(unit=10,file='twoe_int_test.txt')
  !open(unit=11,file='twoe_int_test1.txt')
       tmp=0.0d0
       do i=1, mo_num
         do j=1, mo_num
           do k=1, mo_num
             do l=1, mo_num
               !twoeint=get_two_e_integral(i,j,k,l,mo_integrals_map)
               twoeint=get_two_e_integral(l,k,j,i,mo_integrals_map)
               write(10,'(4(i5,1X),2(f30.16,3X))') l,k,j,i,twoeint,tmp
               !write(11,'(4(i5,1X),2(f30.16,3X))') l,k,j,i,mo_two_e_integral(l,k,j,i),tmp
             enddo
           enddo
         enddo
       enddo

       !do i=1, mo_num
         !do j=1, mo_num
           !do k=1, mo_num
             !do l=1, mo_num
               !twoeint=get_two_e_integral(i,j,k,l,mo_integrals_map)
               !twoeint=get_two_e_integral(l,k,j,i,mo_integrals_map)
               !write(10,'(4(i5,1X),2(f30.16,3X))') i,j,k,l,twoeint,tmp
               !write(11,'(4(i5,1X),2(f30.16,3X))') i,j,k,l,mo_two_e_integral(i,j,k,l),tmp
             !enddo
           !enddo
         !enddo
       !enddo

     close(10)
     close(11)

!  write(*,*) "mo_one_e_integral : prop 0"
!  open(unit=20,file='onee_int_test.txt')
     !do it = 1, n_time !! for ZGrid
     !  if(it==1 .and. ib==1) then
     !    open(unit=20,file='twoe_int_test.txt',status='replace')
     !    write(20,*) bgrid(ib)
     !  else
     !    open(unit=20,file='twoe_int_test.txt',status='old',access='append')
     !    if(it==1 .and. ib/=1) write(20,*) bgrid(ib)
     !  endif
     
     !  write(20,*) zgrid(it)
!       tmp=0.0d0
!       do i=1, mo_num
!         do j=1, mo_num
!           write(20,'(2(i5,1X),4(f30.16,3X))') j,i,mo_one_e_integrals(j,i),tmp,&
!                                               &mo_overlap(j,i),tmp
!         enddo
!       enddo
!     !end do
!     close(20)


end program cippres_two_e_int


