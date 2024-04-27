use bitmasks ! you need to include the bitmasks_module.f90 features
use general
use propdyn, only : esta

 BEGIN_PROVIDER [integer, n_sta_coll_max]
  implicit none
!!!  doc:  maximum number of CSFs for a given CI run (given by the python script generate_csfs.py via the parser.txt file)
!!!  here probably wrong. because we use the total (T + P) det, then should be number of all det. 
     if(n_csf_max<20000) then
        n_sta_coll_max = n_csf_max
     else
        n_sta_coll_max = 20000
     endif
 END_PROVIDER 

 BEGIN_PROVIDER [integer, ib_coll]
  implicit none
   ib_coll = 1
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, b_coll]
  implicit none
!!!! here we give the impact parameter b=0.0
!!!! then we will use: the 1st b equals the 1st b 
!!!! in your input file.
!!!! remember to make the orignal .xml consistent with
!!!! Junwens integral input .xml file.
   b_coll = 0d0
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, v_coll]
  implicit none
     call ezfio_get_cippres_v_coll(v_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, i_state_coll]
  implicit none
     call ezfio_get_cippres_i_state_coll(i_state_coll)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_bound]
  implicit none
     call ezfio_get_cippres_stamin_bound(stamin_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_bound]
  implicit none
     call ezfio_get_cippres_stamax_bound(stamax_bound)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_si]
  implicit none
     call ezfio_get_cippres_stamin_si(stamin_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_si]
  implicit none
     call ezfio_get_cippres_stamax_si(stamax_si)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamin_di]
  implicit none
     call ezfio_get_cippres_stamin_di(stamin_di)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, stamax_di]
  implicit none
     call ezfio_get_cippres_stamax_di(stamax_di)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_time]
  implicit none
     call ezfio_get_cippres_n_time(n_time)
     write(*,*) "ccjia:ezfio_get_cippres_n_time(n_time)",n_time
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_bimp]
  implicit none
     call ezfio_get_cippres_n_bimp(n_bimp)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, n_pcenter]
  implicit none
     call ezfio_get_cippres_n_pcenter(n_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, charge_pcenter, (n_pcenter)]
  implicit none
     call ezfio_get_cippres_charge_pcenter(charge_pcenter)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, bgrid, (n_bimp)]
  implicit none
     call ezfio_get_cippres_bgrid(bgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, zgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_zgrid(zgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [double precision, tgrid, (n_time)]
  implicit none
     call ezfio_get_cippres_tgrid(tgrid)
     !call ezfio_get_cippres_tgrid(tgrid)
 END_PROVIDER 

 BEGIN_PROVIDER [integer, ntdet]
&BEGIN_PROVIDER [integer, ntsta]
&BEGIN_PROVIDER [integer, npdet]
&BEGIN_PROVIDER [integer, npsta]
&BEGIN_PROVIDER [integer, mo_num_t]
&BEGIN_PROVIDER [integer, mo_num_p]
&BEGIN_PROVIDER [integer, Ndet_total]
&BEGIN_PROVIDER [integer, elec_alpha_num_t]
&BEGIN_PROVIDER [integer, elec_beta_num_t]
&BEGIN_PROVIDER [integer, elec_alpha_num_p]
&BEGIN_PROVIDER [integer, elec_beta_num_p]
&BEGIN_PROVIDER [integer, Ndet_Bound]
&BEGIN_PROVIDER [integer, Ndet_SE]
&BEGIN_PROVIDER [integer, Ndet_DE]
&BEGIN_PROVIDER [integer, Ndet_SC]
&BEGIN_PROVIDER [integer, Ndet_DC]

!!! here only for the number of np nt .....
  integer :: i, j
  ntdet = 0
  npdet = 0
  write(*,*) "ccjia: ntdet,npdet",ntdet,npdet
  open(unit=10,file='ints/tcistates_det.txt')
   read(10,*) mo_num_t,ntdet,ntsta
   read(10,*) elec_alpha_num_t, elec_beta_num_t
   write(*,*) "ccjia: mo_num_t,ntdet,ntsta",mo_num_t,ntdet,ntsta
   write(*,*) "elec_alpha_num_t, elec_beta_num_t",elec_alpha_num_t, elec_beta_num_t
  close(10)
  open(unit=10,file='ints/pcistates_det.txt')
   read(10,*) mo_num_p,npdet,npsta
   read(10,*) elec_alpha_num_p, elec_beta_num_p
   write(*,*) "ccjia: mo_num_p,npdet,npsta",mo_num_p,npdet,npsta
   write(*,*) "elec_alpha_num_p, elec_beta_num_p",elec_alpha_num_p, elec_beta_num_p
  close(10)
  !write(*,*) "ccjia:ntdet,npdet,00"
  
  open(unit=10,file='ints/det_total.dat')
  read(10,*) Ndet_total,i,j
  read(10,*) Ndet_Bound, Ndet_SE, Ndet_DE, Ndet_SC, Ndet_DC
  !read(10,*) Ndet_Bound, Ndet_SE, Ndet_SC
  close(10)
  write(*,*) "ccjia:Ndet_total=",Ndet_total
  !write(*,*) "ccjia:Ndet_Bound, SE, SC=",Ndet_Bound, Ndet_SE, Ndet_SC
  write(*,*) "ccjia:Ndet_Bound, SE, DE, SC, DC=",Ndet_Bound, Ndet_SE, Ndet_DE, Ndet_SC, Ndet_DC
  !Ndet_total= 1+ (mo_num_t + mo_num_p)*2 - elec_alpha_num_t-elec_beta_num_t

 END_PROVIDER


 BEGIN_PROVIDER [integer, detalpha, (elec_alpha_num_t+elec_alpha_num_p+1,Ndet_total)]
&BEGIN_PROVIDER [integer, detbeta, (elec_beta_num_t+elec_beta_num_p+1,Ndet_total)]
&BEGIN_PROVIDER [integer, tdeta, (elec_alpha_num_t+1,ntdet)]
&BEGIN_PROVIDER [integer, tdetb, (elec_beta_num_t+1,ntdet)]
&BEGIN_PROVIDER [double precision, ctdet, (ntdet)]
&BEGIN_PROVIDER [double precision, tci_e, (ntsta)]
&BEGIN_PROVIDER [double precision, tci_sta, (ntdet,ntsta)]
&BEGIN_PROVIDER [integer, pdeta, (elec_alpha_num_p+1,npdet)]
&BEGIN_PROVIDER [integer, pdetb, (elec_beta_num_p+1,npdet)]
&BEGIN_PROVIDER [double precision, cpdet, (npdet)]
&BEGIN_PROVIDER [double precision, pci_e, (npsta)]
&BEGIN_PROVIDER [double precision, pci_sta, (npdet,npsta)]
 integer :: i, j, k, l
 integer :: alp, bet


 detalpha(:,:) = 0
 detbeta(:,:) = 0
 tdeta(:,:) = 0
 tdetb(:,:) = 0
 ctdet(:) = 0.0d0
 pdeta(:,:) = 0
 pdetb(:,:) = 0
 cpdet(:) = 0.0d0


 open(unit=10,file='ints/tcistates_det.txt')
    read(10,*)i,j,k
    read(10,*)i,j
  do i = 1, ntdet
    write(*,*) "ntdet=",ntdet
    !write(*,*) "ccjia:1",ctdet(i), (tdeta(j,i),j=1,elec_alpha_num_t), (tdetb(j,i),j=1,elec_beta_num_t)
    read(10,*)ctdet(i), (tdeta(j,i),j=1,elec_alpha_num_t), (tdetb(j,i),j=1,elec_beta_num_t)
    !write(*,*) "ccjia:2",ctdet(i), (tdeta(j,i),j=1,elec_alpha_num_t), (tdetb(j,i),j=1,elec_beta_num_t)
  enddo 
  do i = 1, ntsta
    read(10,*)tci_e(i)
    read(10,*)(tci_sta(j,i),j=1,ntdet)
  enddo 
 close(10)

 open(unit=10,file='ints/pcistates_det.txt')
    read(10,*)i,j,k
    read(10,*)i,j
  do i = 1, npdet
    read(10,*) cpdet(i), (pdeta(j,i),j=1,elec_alpha_num_p), (pdetb(j,i),j=1,elec_beta_num_p)
  enddo
  read(10,*) 
  do i = 1, npsta
    read(10,*)pci_e(i)
    read(10,*)(pci_sta(j,i),j=1,npdet)
  enddo 
 close(10)

!!!! to define the total determinants based on the mo_t+mo_p
  open(unit=10,file='ints/det_total.dat')
  read(10,*) i,alp,bet
  write(*,*) "Ndet_total,alp_ele,bet_ele",i,alp,bet
  read(10,*) 
  do i = 1, Ndet_total
    read(10,*) (detalpha(j,i),j=1,alp), (detbeta(j,i),j=1,bet)
    write(*,*) (detalpha(j,i),j=1,alp), (detbeta(j,i),j=1,bet)
  enddo
  close(10)

 END_PROVIDER 

! BEGIN_PROVIDER [double complex, coll_couplings, (Ndet_total,Ndet_total,n_time)]
! BEGIN_PROVIDER [double complex, allocatable :: coll_couplings, (Ndet_total,Ndet_total,n_time)]
 BEGIN_PROVIDER [double complex, coll_couplings, (Ndet_total,Ndet_total,n_time)]
!&BEGIN_PROVIDER [double complex, coll_sta, (Ndet_total)]
!&BEGIN_PROVIDER [double complex, coll_sta, (Ndet_total,n_time)]
 !!  to use the Ndet_total as the all states. remove the n_csf_max. 
 use general
 use SlaterDeterminant
 implicit none
 integer :: i, j, k, l, imo
 integer :: it
 
 double complex, allocatable :: coll_csf_mat_M(:,:), coll_csf_mat_S(:,:)!, csf_mat(:)
 double complex, allocatable :: coll_csf_mat(:,:)
 double complex, allocatable :: energy_mat(:,:)

 double complex, allocatable :: coll_w1e_mo(:,:,:), coll_ov1e_mo(:,:,:)
 double complex, allocatable :: coll_r12_mo(:,:,:,:,:)

 double complex, dimension(mo_num_t,mo_num_t) :: h1emott
 double complex, dimension(mo_num_p,mo_num_p) :: h1emopp
 double complex, dimension(mo_num_t+mo_num_p,mo_num_t+mo_num_p) :: w1e !ccjia
 double complex, dimension(mo_num_t+mo_num_p,mo_num_t+mo_num_p) :: ovmo !ccjia
 double complex, dimension(mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p) :: r12mo !ccjia
 double complex :: ov, h1e, r12

 integer :: ne, nea, neb, n_mo
 integer :: nsta, ncsf

 double precision :: t1, t2
 logical :: exists, file_e

 integer :: r12test

 integer :: LWORK, INFO
 integer, dimension(:), allocatable :: IPIV
 double complex, dimension(:), allocatable ::  WORK
 double complex, allocatable :: matS1(:,:)
 !double complex, allocatable :: MMM1(:,:), MMM2(:,:)
 !double precision, allocatable :: Rwork(:)

 PROVIDE ezfio_filename !HF_bitmask mo_coef

!! TO BE DONE:  read all integrals and store them in single matrices 1 -> (ntdet+npdet)
!!              read eigevec from cistates_det.txt which should be written manually to incorporate both target and proj. MOs
!!              initiate psi with target CI coeff.
!!              test lowdin rules code (ok for tt-1e ints.)

  print*,'Computing coll_couplings', b_coll 
  call cpu_time(t1)
  !coll_sta(:)= dcmplx(0.0D0,0.0D0)
  coll_couplings(:,:,:) = dcmplx(0.0D0,0.0D0)
  
  !allocate(coll_couplings(Ndet_total,Ndet_total,n_time))


  allocate(coll_csf_mat_M(1:Ndet_total,1:Ndet_total))
  allocate(coll_csf_mat_S(1:Ndet_total,1:Ndet_total))
  allocate(coll_csf_mat(1:Ndet_total,1:Ndet_total))
  allocate(energy_mat(1:Ndet_total,1:Ndet_total))
  allocate(coll_w1e_mo(mo_num_t+mo_num_p,mo_num_t+mo_num_p,n_time))
  allocate(coll_ov1e_mo(mo_num_t+mo_num_p,mo_num_t+mo_num_p,n_time))
  allocate(coll_r12_mo(mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p,mo_num_t+mo_num_p,n_time))
  allocate(IPIV(1:Ndet_total))
  allocate(matS1(1:Ndet_total,1:Ndet_total))
  LWORK = 3*Ndet_total
  allocate(WORK(LWORK))
  WORK=dcmplx(0.0D0,0.0D0)

  coll_w1e_mo = dcmplx(0.0D0,0.0D0)
  coll_ov1e_mo = dcmplx(0.0D0,0.0D0)
  coll_r12_mo = dcmplx(0.0D0,0.0D0)

  nsta = Ndet_total
  if(nsta>20000) then
    write(*,*) "n total state or det. =",nsta, "nsta > n_sta_coll_max=20000, I stop"
    print*, "use less states......"
    stop
  endif
  write(*,*) "nsta, ntsta, npsta(wrong)",nsta, ntsta, npsta
  ncsf = Ndet_total

  if(dabs(b_coll-0.0d0) < 0.0001d0 ) goto 1001 !!(to avoid the bcoll==0.0 problem?)
  call integralread(h1emott, h1emopp, coll_w1e_mo, coll_ov1e_mo, coll_r12_mo, b_coll, n_time, mo_num_t, mo_num_p)

  nea = elec_alpha_num_t+elec_alpha_num_p !!! consider the T and P=0 ????
  neb = elec_beta_num_t+elec_beta_num_p !!! consider the T and P=0 ????
  ne = nea + neb
  n_mo = mo_num_t+mo_num_p !+mo_num  !! consider the T and P
  h1e = dcmplx(0.0D0,0.0D0)!0.0d0
  ov  = dcmplx(0.0D0,0.0D0)!0.0d0
  r12 = dcmplx(0.0D0,0.0D0)!0.0d0
  w1e(:,:) = coll_w1e_mo(:,:,1) 
!  do i = 1, n_mo
!    w1e(i,i) = w1e(i,i) - 2d0/zgrid(1)
!  enddo
  ovmo(:,:) = coll_ov1e_mo(:,:,1)
  r12mo(:,:,:,:) = coll_r12_mo(:,:,:,:,1)
  energy_mat(:,:) = 0d0
  do i=1,Ndet_total  ! 
    call lowdin(ne,nea,neb,n_mo,ovmo,w1e,r12mo,detalpha(:,i),detalpha(:,i),detbeta(:,i),detbeta(:,i),ov,h1e,r12)
    esta(i) = h1e + r12
    energy_mat(i,i) = esta(i)
    esta(i) = 0d0 !nico
  enddo

   do i = 1, mo_num_t
!     do j = 1, mo_num_t
     j = i
!       coll_w1e_mo(j,i,:) = coll_w1e_mo(j,i,:) - h1emott(j,i)
!     enddo
   enddo

   do i = 1, mo_num_p
!    do j = 1, mo_num_p
     j = i
!     coll_w1e_mo(j+mo_num_t,i+mo_num_t,:) = coll_w1e_mo(j+mo_num_t,i+mo_num_t,:) - h1emopp(j,i)
!    enddo
   enddo

! here, 'it' means the zGrid.
  do it = 1, n_time

    w1e(:,:) = coll_w1e_mo(:,:,it)
    ovmo(:,:) = coll_ov1e_mo(:,:,it)
    r12mo(:,:,:,:) = coll_r12_mo(:,:,:,:,it)
!    r12mo(:,:,:,:) = 0d0
!    r12mo(1:mo_num_t,:1:mo_num_t,1:mo_num_t,1:mo_num_t) = coll_r12_mo(1:mo_num_t,1:mo_num_t,1:mo_num_t,1:mo_num_t,it)
!    r12mo(1:mo_num_t,mo_num_t+1:mo_num_t+mo_num_p,1:mo_num_t,mo_num_t+1:mo_num_t+mo_num_p) = &
!    coll_r12_mo(1:mo_num_t,mo_num_t+1:mo_num_t+mo_num_p,1:mo_num_t,mo_num_t+1:mo_num_t+mo_num_p,it)

    coll_csf_mat_M(:,:) = dcmplx(0.0D0,0.0D0)
    coll_csf_mat_S(:,:) = dcmplx(0.0D0,0.0D0)
    coll_csf_mat(:,:) = dcmplx(0.0D0,0.0D0)
    !csf_mat(:) = 0d0

    ! touch ntdet

    nea = elec_alpha_num_t+elec_alpha_num_p !!! consider the T and P=0 ????
    neb = elec_beta_num_t+elec_beta_num_p !!! consider the T and P=0 ????
    ne = nea + neb
    n_mo = mo_num_t+mo_num_p !+mo_num  !! consider the T and P
    h1e = dcmplx(0.0D0,0.0D0)!0.0d0
    ov  = dcmplx(0.0D0,0.0D0)!0.0d0
    r12 = dcmplx(0.0D0,0.0D0)!0.0d0
    write(52,*)"###################################################################"
    write(52,*)"###################################################################"
    write(52,*)zgrid(it)
    write(52,*)"###################################################################"
    write(52,*)"###################################################################"
    do i=1,Ndet_total  ! 
      do j=1,Ndet_total   !

        write(52,*)detalpha(:,i),detbeta(:,i),detalpha(:,j),detbeta(:,j)
        call lowdin(ne,nea,neb,n_mo,ovmo,w1e,r12mo,detalpha(:,i),detalpha(:,j),detbeta(:,i),detbeta(:,j),ov,h1e,r12)
        coll_csf_mat_M(j,i) = h1e + r12 !real(h1e)!*DBLE(i+j)*0.01D0
        coll_csf_mat_S(j,i) = ov !real(ov)!*DBLE(i+j)*0.02D0
        write(1111,*) i,j,coll_csf_mat_M(j,i),coll_csf_mat_S(j,i)
        write(1112,*) i,j,h1e,r12,ov
        !coll_csf_mat_M(i,j) = dconjg(coll_csf_mat_M(j,i))
        !coll_csf_mat_S(i,j) = dconjg(coll_csf_mat_S(j,i))
      enddo
    enddo

    !!! compute M - ES
    !!matS1(:,:) = coll_csf_mat_M(:,:) - matmul(energy_mat(:,:),coll_csf_mat_S(:,:))
!nico    do i=1,Ndet_total  ! 
!nico      coll_csf_mat_M(i,i) = coll_csf_mat_M(i,i) - energy_mat(i,i)*coll_csf_mat_S(i,i)
!nico    enddo
    IPIV=0
   
    !!! compute the S^(-1) and also the S^(-1) * H
    matS1(:,:) = coll_csf_mat_S(:,:)

    CALL ZGETRF( Ndet_total, Ndet_total, matS1, Ndet_total, IPIV, INFO )

    CALL ZGETRI( Ndet_total, matS1, Ndet_total, IPIV, WORK, LWORK, INFO )
    
    coll_csf_mat(:,:)= matmul(matS1(:,:),coll_csf_mat_M(:,:))
    coll_couplings(1:Ndet_total,1:Ndet_total,it) = coll_csf_mat(1:Ndet_total,1:Ndet_total) !nico

  end do  !! (do it = 1, n_time) for ZGrid

  deallocate(matS1,WORK,IPIV,energy_mat)
  deallocate(coll_w1e_mo, coll_ov1e_mo, coll_r12_mo)
  !deallocate
  call cpu_time(t2)
  print*,t2-t1
  print*,' '
  deallocate(coll_csf_mat,coll_csf_mat_S,coll_csf_mat_M)

1001 continue

END_PROVIDER

