program cippres_prop_collision
  use general
  use propdyn
 ! create a routine for one-e general matrix with determinants

  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_prop_collision solves the TDSE using the matrix elements from cippres_collision
  END_DOC

  integer :: nsta_bound
  double precision :: t1, t2, tdyn, tmp
  logical :: exists, file_e
  integer :: zeroVp, oneVp, twoVp
  integer :: i, j, ib, k, l, it, ni, nf
  integer :: ii,jj,kk
  integer :: energytest, energytype 

  PROVIDE ezfio_filename !HF_bitmask mo_coef

  if (mpi_master) then
   call ezfio_has_cippres_n_time(exists)
    if (exists) then
      call ezfio_has_cippres_n_bimp(exists)
      if (exists) then
       call ezfio_has_cippres_tgrid(exists)
       if (exists) then
         call ezfio_has_cippres_zgrid(exists)
         if (exists) then
           call ezfio_has_cippres_bgrid(exists)
            if (exists) then
             call ezfio_has_cippres_v_coll(exists)
             if (exists) then
               call ezfio_has_cippres_i_state_coll(exists)
               if (exists) then
                 call ezfio_has_cippres_stamin_bound(exists)
                 if (exists) then
                   call ezfio_has_cippres_stamax_bound(exists)
                 endif
               endif
             endif
           endif
         endif
       endif
     endif
   endif
 endif


  if (exists) then
   call ezfio_get_cippres_n_time(n_time)
   call ezfio_get_cippres_n_bimp(n_bimp)
   call ezfio_get_cippres_zgrid(zgrid)
   call ezfio_get_cippres_tgrid(tgrid)
   call ezfio_get_cippres_bgrid(bgrid)
   call ezfio_get_cippres_stamin_bound(stamin_bound)
   call ezfio_get_cippres_stamax_bound(stamax_bound)
   call ezfio_get_cippres_stamin_si(stamin_si)
   call ezfio_get_cippres_stamax_si(stamax_si)
   call ezfio_get_cippres_stamin_di(stamin_di)
   call ezfio_get_cippres_stamax_di(stamax_di)
   call ezfio_get_cippres_i_state_coll(i_state_coll)
   call ezfio_get_cippres_v_coll(v_coll)
   print*, "" 
   print*, "Collision info"
   print*, "" 
   print*,"n_bimp,n_time,i_state_coll,v_coll =", n_bimp, n_time, i_state_coll, v_coll
   print*,"stamin, stamax, nsta_bound =", stamin_bound, stamax_bound, stamax_bound-stamin_bound+1
   print*,"impact. param. =", bgrid
!   print*,"tgrid = ", tgrid
!   print*,"zgrid = ", zgrid
   print*, "" 


!   call ezfio_get_cippres_ici1(ici1)
   call ezfio_get_cippres_n_csf_cippres(n_csf_cippres)

   ntime= n_time
   !nsave_time= n_time !! nico
   !nsave_time=1  !ccjia test, change to remove the loop in propdyn-ABM

   nsta_bound = stamax_bound-stamin_bound+1
   nsi  = stamax_si-stamin_si+1
   ndi  = stamax_di-stamin_di+1
   ni = stamin_bound
   nf = stamax_di
   if(stamax_di == 0) then  
      nsi = 0
      ndi = 0 
      ni = stamin_bound
      nf = stamax_bound
   endif
!!!
  nsta = Ndet_total
  !Det_Tar = ntdet
  vproj = v_coll

  energytest=1
  energytype=1
  vprojtest=1
  nsave_time=n_time
  output_cmd=0

  allocate(mcoup(ntime,nsta,nsta),timegrid(ntime),esta(nsta))
  !allocate(mcoup(ntime,nsta,nsta),timegrid(ntime),esta(nsta),mat(ntime,nsta,nsta))
  allocate(psi(nsta))
  allocate(rmat2intrp(ntime,nsta,nsta),cmat2intrp(ntime,nsta,nsta))
  allocate(matintrp(nsta,nsta))

   !allocate(g_si(1:ntime,1:nsta,1:nsta),g_ddi(1:ntime,1:nsta,1:nsta),g_sdi(1:ntime,1:nsta))
  !allocate(g_si_intrp(1:ntime,1:nsta,1:nsta),g_ddi_intrp(1:ntime,1:nsta,1:nsta),g_sdi_intrp(1:ntime,1:nsta))

   timegrid(1:ntime) = tgrid(1:n_time)
   esta=0.0d0

!!!
   zeroVp = Ndet_Bound + Ndet_SE + Ndet_DE
   oneVp  = Ndet_Bound + Ndet_SE + Ndet_DE + Ndet_SC
   twoVp  = Ndet_Bound + Ndet_SE + Ndet_DE + Ndet_SC + Ndet_DC
!!!
   !write(21,*) "n_tsta, n_psta, n_tpsta",n_tsta, n_psta, n_tpsta
   write(*,*) "zeroVp, oneVp, twoVp", zeroVp, oneVp, twoVp

   open(unit=20,file='Prop_collision.out')
   write(20,*) n_bimp, Ndet_total
   close(20)
   print*,'start dyn'
   do ib = 1, n_bimp

    b_coll = bgrid(ib)
    ib_coll = ib
    touch b_coll

    mcoup(:,:,:) = dcmplx(0.0D0,0.0D0)
    do it = 1, ntime !! for ZGrid
      mcoup(it,1:nsta,1:nsta) = coll_couplings(1:nsta,1:nsta,it)
    end do
    FREE coll_couplings
    
!    if(output_cmd==1) then
      do it = 1, ntime !! for ZGrid
        write(997,'(9(f20.12,1X))')zgrid(it),mcoup(it,1,1),mcoup(it,1,2),mcoup(it,1,3)
        write(998,'(9(f20.12,1X))')zgrid(it),mcoup(it,2,1),mcoup(it,2,2),mcoup(it,2,3)
        write(999,'(9(f20.12,1X))')zgrid(it),mcoup(it,3,1),mcoup(it,3,2),mcoup(it,3,3)
      end do
!        do i=1, nsta
!          do j=1, nsta
!            !mcoup(it,j,i) = coll_couplings(j,i,it)
!            write(999,*) "it,nsta,nsta,j,i,mcoup(it,j,i)",it,nsta,nsta,j,i,mcoup(it,j,i)
!          end do
!        end do
!    end if
!!!!!!!!!!!!!!!!!
    n_tsta = Ndet_Bound + Ndet_SE + Ndet_DE
    n_tpsta = Ndet_SC
    n_psta = Ndet_DC

    do i=1, nsta
!      esta(i) = mcoup(1,i,i)
      write(*,*)"esta:",i,esta(i)
    end do
!!!!!!!!!!!!!!!!

    psi(:) = 0.0d0
    psi(1) = 1.0d0 ! only true for CIS, otherwise use lines below (Nico 22.04.2024)
    !open(unit=111,file='det_coeff.dat',status='old')
!    do j=1, n_tsta
!      psi(j)=dcmplx(tci_sta(j,1), 0.0D0)
!      write(*,*) "psi(j)=dcmplx(tci_sta(j,1), 0.0D0)",psi(j),tci_sta(j,1)
!    end do

! propagation 
    call cpu_time(t1)
    call dyn
    call cpu_time(t2)
    tdyn = t2-t1
    write(*,*)'DYN takes',tdyn

    write(*,*) bgrid(ib),(cdabs(psi(i))**2.0d0,i=1,nsta), sum(cdabs(psi(:))**2.0d0),v_coll
    open(unit=20,file='Prop_collision.out',status='old',access='append')
    write(20,'(20000(f15.8,1X))')bgrid(ib),(cdabs(psi(i))**2.0d0,i=1,nsta), sum(cdabs(psi(:))**2.0d0),v_coll
    close(20)

   enddo 

   !close(20)

   deallocate(mcoup,timegrid,esta)
   deallocate(matintrp,psi,rmat2intrp,cmat2intrp)

  else

    print*, "Z and/or b grids are not setup correctly"
    stop

  endif

end program cippres_prop_collision
