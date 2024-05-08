subroutine integralread(h1emott, h1emopp, coll_w1e_mo, coll_ov1e_mo, coll_r12_mo, impab, iz_total, mo_num_t, mo_num_p)
  implicit none

  double precision :: zg, bg
  double precision :: w1ea,w1eb,ov1ea,ov1eb
  double precision :: r12a,r12b
  integer :: i, j, k, l, ib, iz
  integer :: mm,nn,oo,pp
 
  double precision, intent(in) :: impab 
  integer, intent(in) :: iz_total, mo_num_t, mo_num_p
  double complex, dimension(mo_num_t+mo_num_p,mo_num_t+mo_num_p,iz_total), intent(in out) :: coll_w1e_mo, coll_ov1e_mo
  double complex, dimension(mo_num_t,mo_num_t), intent(out) :: h1emott
  double complex, dimension(mo_num_p,mo_num_p), intent(out) :: h1emopp
  double complex, dimension(mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p,iz_total), intent(in out) :: coll_r12_mo
 

  character(len=100) :: filett,filetp,filept,filepp,fileread
  character(len=100) :: newfilett,newfilepp
  character(len=100) :: filetwo1,filetwo2,filetwo3,filetwo4
  character(len=100) :: impactb_value,impab_filename, impab_name
  integer :: ilen

  write(*,*) "ccjia:","read coll_w1e_mo, coll_ov1e_mo"

 impab_name="impactb_"
 write(impactb_value,"(f10.2)") impab
 impab_filename= Trim(AdjustL(impab_name))//Trim(AdjustL(impactb_value))
 ilen=index(impab_filename,' ')


 newfilett= 'ints/'//'h1emo_tt.txt'
 open(unit=11,file=newfilett,status='old')
  read(11,*)nn
  do i = 1, mo_num_t
    do j = 1, mo_num_t
     read(11,*)h1emott(j,i)
    enddo
  enddo 
 close(11)
 newfilepp= 'ints/'//'h1emo_pp.txt'
 open(unit=11,file=newfilepp,status='old')
  read(11,*)nn
  do i = 1, mo_num_p
    do j = 1, mo_num_p
     read(11,*)h1emopp(j,i)
    enddo 
  enddo 
 close(11)

 h1emott(:,:) = 0d0 !nico
 h1emopp(:,:) = 0d0 !nico

 filett= 'ints/'//impab_filename(1:ilen-1)//'/'//'onee_int_tt.txt'
 filetp= 'ints/'//impab_filename(1:ilen-1)//'/'//'onee_int_tp.txt'
 filept= 'ints/'//impab_filename(1:ilen-1)//'/'//'onee_int_pt.txt'
 filepp= 'ints/'//impab_filename(1:ilen-1)//'/'//'onee_int_pp.txt'

 open(unit=11,file=filett,status='old')
 open(unit=12,file=filetp,status='old')
 open(unit=13,file=filept,status='old')
 open(unit=14,file=filepp,status='old')
!write(*,*) "now read one e int", impab, bg
!!!
!  do ib = ib_count, n_bimp
    read(11,*) bg
    read(12,*) bg
    read(13,*) bg
    read(14,*) bg
    
    do iz = 1, iz_total
      read(11,*) zg
      read(12,*) zg
      read(13,*) zg
      read(14,*) zg

      !write(1113,*) "b_coll= ",bg," ib=",ib,","," it/iz=",iz,",","zgrid=",zg,"
      !output for coll_w1e_mo(j,i,iz,ib)"
      do i = 1, mo_num_t
        do j = 1, mo_num_t
          read(11,*) k,l,w1ea,w1eb,ov1ea,ov1eb
          coll_w1e_mo(j,i,iz)= dcmplx(w1ea,w1eb) + h1emott(j,i)
          coll_ov1e_mo(j,i,iz)= dcmplx(ov1ea,ov1eb)!!mcoup(izgrid,j,i),movl(izgrid,j,i)
          !write(1113,"(A,2I6,4F20.16)") "j=mo_num,
          !i=mo_num,coll_w1e_and_ov1e__mo(j,i,iz,ib)=",j,i,coll_w1e_mo(j,i,iz,ib),coll_ov1e_mo(j,i,iz,ib)
        enddo
      enddo

      !  do i = 1, npmo
      !    do j = 1, ntmo !! code-2e-Junwen
      !        write(10,'(2(i5,1X),4(f20.16,1X))')j,i,mcoup(izgrid,j,i+ntmo),movl(izgrid,j,i+ntmo)
      !    enddo
      !  enddo
      do i = 1, mo_num_p
        do j = 1, mo_num_t
          read(12,*) k,l,w1ea,w1eb,ov1ea,ov1eb! (j,i)=(TP)-element
          coll_w1e_mo(j,i+mo_num_t,iz)= dcmplx(w1ea,w1eb)
          coll_ov1e_mo(j,i+mo_num_t,iz)= dcmplx(ov1ea,ov1eb)
        enddo
      enddo

      !  do i = 1, ntmo
      !    do j = 1, npmo !! code-2e-Junwen
      !        write(11,'(2(i5,1X),4(f20.16,1X))')j,i,mcoup(izgrid,j+ntmo,i),movl(izgrid,j+ntmo,i)
      !    enddo
      !  enddo
      do i = 1, mo_num_t
        do j = 1, mo_num_p
          read(13,*) k,l,w1ea,w1eb,ov1ea,ov1eb ! (j,i)=(PT)-element
          coll_w1e_mo(j+mo_num_t,i,iz)= dcmplx(w1ea,w1eb)
          coll_ov1e_mo(j+mo_num_t,i,iz)= dcmplx(ov1ea,ov1eb)
        enddo
      enddo

      do i = 1, mo_num_p
        do j = 1, mo_num_p
          read(14,*) k,l,w1ea,w1eb,ov1ea,ov1eb
          coll_w1e_mo(j+mo_num_t,i+mo_num_t,iz)= dcmplx(w1ea,w1eb) + h1emopp(j,i)
          coll_ov1e_mo(j+mo_num_t,i+mo_num_t,iz)= dcmplx(ov1ea,ov1eb)
        enddo
      enddo

    end do
!  end do

  print*,'read one e int ok'
  close(11)
  close(12)
  close(13)
  close(14)

! test r12mo rep.
! 
!write(*,*) "now read two e int", impab, bg
open(unit=11,file='ints/twoe_int_tttt.txt',status='old')
!    read(11,*)bg
    do iz = 1, iz_total
!      read(11,*)zg
      if(iz==1) then
      do i=1, mo_num_t
        do j=1, mo_num_t
          do k=1, mo_num_t
            do l=1, mo_num_t
              read(11,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k,j,i,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
      else
      do i=1, mo_num_t
        do j=1, mo_num_t
          do k=1, mo_num_t
            do l=1, mo_num_t
              !read(11,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k,j,i,iz)= coll_r12_mo(l,k,j,i,1)!dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
      end if
    end do ! for iz=zgrid(iz)
close(11)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_ptpt.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_tppt.txt'
filetwo3= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_pttp.txt'
filetwo4= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_tptp.txt'
open(unit=12,file=filetwo1,status='old')
open(unit=13,file=filetwo2,status='old')
open(unit=14,file=filetwo3,status='old')
open(unit=15,file=filetwo4,status='old')
!open(unit=12,file='ints/twoe_int_ptpt.txt',status='old')
!open(unit=13,file='ints/twoe_int_tppt.txt',status='old')
!open(unit=14,file='ints/twoe_int_pttp.txt',status='old')
!open(unit=15,file='ints/twoe_int_tptp.txt',status='old')
    read(12,*)bg
    read(13,*)bg
    read(14,*)bg
    read(15,*)bg
    do iz = 1, iz_total
      read(12,*)zg
      read(13,*)zg
      read(14,*)zg
      read(15,*)zg
      do i=1, mo_num_p
        do j=1, mo_num_t
          do k=1, mo_num_p
            do l=1, mo_num_t
              read(12,*) mm,nn,oo,pp, r12a,r12b !dcmplx(0.0d0,0.0d0)
              coll_r12_mo(k+mo_num_t,l,i+mo_num_t,j,iz)= dcmplx(r12a,r12b)
              read(13,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k+mo_num_t,i+mo_num_t,j,iz)= dcmplx(r12a,r12b)
              read(14,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(k+mo_num_t,l,j,i+mo_num_t,iz)= dcmplx(r12a,r12b)
              read(15,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k+mo_num_t,j,i+mo_num_t,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
    end do ! for iz=zgrid(iz)
close(12)
close(13)
close(14)
close(15)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_pttt.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_tptt.txt'
open(unit=16,file=filetwo1,status='old')
open(unit=17,file=filetwo2,status='old')
!open(unit=16,file='ints/twoe_int_pttt.txt',status='old')
!open(unit=17,file='ints/twoe_int_tptt.txt',status='old')
    read(16,*)bg
    read(17,*)bg
    do iz = 1, iz_total
      read(16,*)zg
      read(17,*)zg
      do i=1, mo_num_t
        do j=1, mo_num_t
          do k=1, mo_num_p
            do l=1, mo_num_t
              read(16,*) mm,nn,oo,pp, r12a,r12b!dcmplx(0.0d0,0.0d0)
              coll_r12_mo(k+mo_num_t,l,i,j,iz)= dcmplx(r12a,r12b)
              read(17,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k+mo_num_t,i,j,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
    end do ! for iz=zgrid(iz)
close(16)
close(17)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_ttpt.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_tttp.txt'
open(unit=18,file=filetwo1,status='old')
open(unit=19,file=filetwo2,status='old')
!open(unit=18,file='ints/twoe_int_ttpt.txt',status='old')
!open(unit=19,file='ints/twoe_int_tttp.txt',status='old')
    read(18,*)bg
    read(19,*)bg
    do iz = 1, iz_total
      read(18,*)zg
      read(19,*)zg
      do i=1, mo_num_t
        do j=1, mo_num_p
          do k=1, mo_num_t
            do l=1, mo_num_t
              read(18,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k,j+mo_num_t,i,iz)= dcmplx(r12a,r12b)
              read(19,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k,i,j+mo_num_t,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
    end do ! for iz=zgrid(iz)
close(18)
close(19)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_pptt.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_ttpp.txt'
open(unit=20,file=filetwo1,status='old')
open(unit=21,file=filetwo2,status='old')
!open(unit=20,file='ints/twoe_int_pptt.txt',status='old')
!open(unit=21,file='ints/twoe_int_ttpp.txt',status='old')
    read(20,*)bg
    read(21,*)bg
    do iz = 1, iz_total
      read(20,*)zg
      read(21,*)zg
      do i=1, mo_num_p
        do j=1, mo_num_p
          do k=1, mo_num_t
            do l=1, mo_num_t
              read(20,*) mm,nn,oo,pp, r12a,r12b !dcmplx(0.0d0,0.0d0)
              coll_r12_mo(j+mo_num_t,i+mo_num_t,l,k,iz)= dcmplx(r12a,r12b)
              read(21,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k,j+mo_num_t,i+mo_num_t,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo !!??? 
    end do ! for iz=zgrid(iz)
close(20)
close(21)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_ptpp.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_tppp.txt'
open(unit=22,file=filetwo1,status='old')
open(unit=23,file=filetwo2,status='old')
!open(unit=22,file='ints/twoe_int_ptpp.txt',status='old')
!open(unit=23,file='ints/twoe_int_tppp.txt',status='old')
    read(22,*)bg
    read(23,*)bg
    do iz = 1, iz_total
      read(22,*)zg
      read(23,*)zg
      do i=1, mo_num_p
        do j=1, mo_num_p
          do k=1, mo_num_p
            do l=1, mo_num_t
              read(22,*) mm,nn,oo,pp, r12a,r12b!dcmplx(0.0d0,0.0d0)
              coll_r12_mo(k+mo_num_t,l,i+mo_num_t,j+mo_num_t,iz)= dcmplx(r12a,r12b)
              read(23,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l,k+mo_num_t,j+mo_num_t,i+mo_num_t,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo!! ???
    end do ! for iz=zgrid(iz)
close(22)
close(23)
!!!!####
!!!!####
filetwo1= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_pppt.txt'
filetwo2= 'ints/'//impab_filename(1:ilen-1)//'/'//'twoe_int_pptp.txt'
open(unit=24,file=filetwo1,status='old')
open(unit=25,file=filetwo2,status='old')
!open(unit=24,file='ints/twoe_int_pppt.txt',status='old')
!open(unit=25,file='ints/twoe_int_pptp.txt',status='old')
    read(24,*)bg
    read(25,*)bg
    do iz = 1, iz_total
      read(24,*)zg
      read(25,*)zg
      do i=1, mo_num_p
        do j=1, mo_num_t
          do k=1, mo_num_p
            do l=1, mo_num_p
              read(24,*) mm,nn,oo,pp, r12a,r12b !dcmplx(0.0d0,0.0d0)
              coll_r12_mo(l+mo_num_t,k+mo_num_t,i+mo_num_t,j,iz)=dcmplx(r12a,r12b)
              read(25,*) mm,nn,oo,pp, r12a,r12b
              coll_r12_mo(l+mo_num_t,k+mo_num_t,j,i+mo_num_t,iz)=dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
    end do ! for iz=zgrid(iz)
close(24)
close(25)
!!!!####
!!!!####
open(unit=26,file='ints/twoe_int_pppp.txt',status='old')
!    read(26,*)bg
    do iz = 1, iz_total
!      read(26,*)zg
      if(iz==1) then
      do i=1, mo_num_p
        do j=1, mo_num_p
          do k=1, mo_num_p
            do l=1, mo_num_p
              read(26,*) mm,nn,oo,pp, r12a,r12b! dcmplx(0.0d0,0.0d0)
              coll_r12_mo(l+mo_num_t,k+mo_num_t,j+mo_num_t,i+mo_num_t,iz)= dcmplx(r12a,r12b)
!              !coll_r12_mo(k+mo_num_t,l+mo_num_t,j+mo_num_t,i+mo_num_t,iz)= dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
      else
      do i=1, mo_num_p
        do j=1, mo_num_p
          do k=1, mo_num_p
            do l=1, mo_num_p
              !read(26,*) mm,nn,oo,pp, r12a,r12b! dcmplx(0.0d0,0.0d0)
              coll_r12_mo(l+mo_num_t,k+mo_num_t,j+mo_num_t,i+mo_num_t,iz)= &
             &coll_r12_mo(l+mo_num_t,k+mo_num_t,j+mo_num_t,i+mo_num_t,1) !dcmplx(r12a,r12b)
!              !coll_r12_mo(k+mo_num_t,l+mo_num_t,j+mo_num_t,i+mo_num_t,iz)=
!              dcmplx(r12a,r12b)
            enddo
          enddo
        enddo
      enddo
      end if
    end do ! for iz=zgrid(iz)
close(26)

!  coll_r12_mo(:,:,:,:,:) = 0d0 !nico
  print*,'read two e int ok'

end subroutine integralread
