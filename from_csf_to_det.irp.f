
BEGIN_PROVIDER [integer, n_sta_cistate_analysis]
implicit none
   call ezfio_get_cippres_n_sta_cistate_analysis(n_sta_cistate_analysis)
END_PROVIDER

subroutine myprint_det(string,Nint,output) 
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant using the '+-' notation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(2)

  call bitstring_to_str( output(1), string(1,1), Nint )
  call bitstring_to_str( output(2), string(1,2), Nint )
!  print *,  trim(output(1))
!  print *,  trim(output(2))

end
 

subroutine from_csf_to_det(irun)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in) :: irun
 integer :: i,j,k,idet_tmp,istate
 double precision :: hij, accu

 character*(2048)   :: output(2)

 print*,'*************************************'
 print*,'n_det_max_csf, n_csf_max',n_det_max_csf, n_csf_max
 print*,'*************************************'

 idet_tmp = 0
 do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
    idet_tmp += 1
   enddo 
  enddo

 N_det = idet_tmp 
! N_states = n_csf_max
 N_states = n_sta_cistate_analysis
 psi_det_size =  idet_tmp
 touch psi_det_size
 touch N_det N_states
 print*,'N_det,psi_det_size',N_det,psi_det_size
 idet_tmp = 0
 do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
    idet_tmp += 1
    psi_det(:,:,idet_tmp) = csf_basis(:,:,k,i,irun) 
    do istate = 1, n_sta_cistate_analysis
     psi_coef(idet_tmp,istate) = coef_det_csf_basis(k,i,irun) * eigvectors_cippres(i,istate,irun)
    enddo
   enddo 
  enddo
 print*,'idet_tmp          = ',idet_tmp
 print*,'N_det,N_states    = ',N_det,N_states
 touch psi_det psi_coef 
 print*,'N_det_alpha_unique',N_det_alpha_unique
 print*,'N_det_beta_unique ',N_det_beta_unique
 print*,'*N_det,N_states    = ',N_det,N_states

!NICO double precision, allocatable :: coef_alpha_beta(:,:,:)
!NICO integer(bit_kind), allocatable :: psi_alpha_uniq_tmp(:,:), psi_beta_uniq_tmp(:,:)
!NICO allocate(psi_alpha_uniq_tmp(N_int,N_det_alpha_unique),psi_beta_uniq_tmp(N_int,N_det_beta_unique))

!NICO psi_alpha_uniq_tmp(:,1:N_det_alpha_unique) = psi_det_alpha_unique(:,1:N_det_alpha_unique)
!NICO psi_beta_uniq_tmp(:,1:N_det_beta_unique) = psi_det_beta_unique(:,1:N_det_beta_unique)

! allocate ( coef_alpha_beta(N_det_alpha_unique,N_det_beta_unique,N_states)  )
 
!NICO integer :: get_index_in_psi_det_alpha_unique,get_index_in_psi_det_beta_unique
!NICO integer :: n_alpha_tmp,n_beta_tmp
!NICO n_alpha_tmp = N_det_alpha_unique
!NICO n_beta_tmp = N_det_beta_unique
!NICO coef_alpha_beta = 0.d0
!NICO do k=1,N_det
!NICO   i = get_index_in_psi_det_alpha_unique(psi_det(1,1,k),N_int)
!NICO   ASSERT (i>0)
!NICO   ASSERT (i<=N_det_alpha_unique)
!NICO 
!NICO   j = get_index_in_psi_det_beta_unique (psi_det(1,2,k),N_int)
!NICO   ASSERT (j>0)
!NICO   ASSERT (j<=N_det_beta_unique)
!NICO   do istate = 1, N_states
!NICO    coef_alpha_beta(i,j,istate) += psi_coef(k,istate)
!NICO   enddo
!NICO enddo
!NICO print*,'N_det_alpha_unique',N_det_alpha_unique
!NICO print*,'N_det_beta_unique ',N_det_beta_unique

!NICON_det = N_det_alpha_unique * N_det_beta_unique 
!NICOtouch N_det 
!NICOidet_tmp = 0
!NICOdo j = 1, n_alpha_tmp 
!NICO do i = 1, n_beta_tmp 
!NICO  idet_tmp += 1
!NICO  psi_det(:,1,idet_tmp) = psi_alpha_uniq_tmp(:,j)
!NICO  psi_det(:,2,idet_tmp) = psi_beta_uniq_tmp(:,i)
!NICO   print*,'idet = ', idet_tmp
!NICO   call print_det(psi_det(1,1,idet_tmp),N_int)
!NICO   print*,''
!NICO  do istate = 1, N_states
!NICO   psi_coef(idet_tmp,istate) = coef_alpha_beta(j,i,istate)
!NICO  enddo
!NICO enddo
!NICOenddo
!NICO
!NICOtouch psi_det psi_coef 
!NICO
 open(unit=10,file='cistates_det.txt')
 write(10,*)mo_num,N_det,N_states
 do i = 1, N_det
  write(10,*)i
   call myprint_det(psi_det(1,1,i),N_int,output)
   write(10,'((a))')  trim(output(1))
   write(10,'((a))')  trim(output(2))
 enddo
 
  do istate = 1, N_states
  accu = 0.d0
  do i = 1, N_det
   do j = 1, N_det
    call i_H_j_one_e(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
    accu += hij * psi_coef(j,istate) * psi_coef(i,istate) 
   enddo
  enddo
 
    write(10,'(3(f30.16))') eigvalues_cippres(istate,irun), accu, nuclear_repulsion
    !write(10,*)eigvalues_cippres(istate,irun), accu, dip_couplings_cippres(1,istate), nuclear_repulsion
    do i = 1, N_det
     write(10,"(f30.16)") psi_coef(i,istate)
    enddo  
  enddo
  close(10)
  stop
!NICO do istate = 1, N_states
!NICO  double precision :: accu,hij
!NICO  accu = 0.d0
!NICO  do i = 1, N_det
!NICO   do j = 1, N_det
!NICO!    print*,'i,j',i,j
!NICO!    call print_det(psi_det(1,1,i),N_int)
!NICO!    call print_det(psi_det(1,1,j),N_int)
!NICO    call i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
!NICO    accu += hij * psi_coef(j,istate) * psi_coef(i,istate) 
!NICO   enddo
!NICO  enddo
!NICO  print*,'accu              = ',accu+nuclear_repulsion
!NICO  print*,'eigvalues_cippres = ',eigvalues_cippres(istate,irun)
!NICO enddo
 
end

