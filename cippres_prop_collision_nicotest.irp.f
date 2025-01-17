program cippres_prop_collision_nicotest
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
   print*,"Target HF Energy =", t_ehf
   do i = 1, tmo_num
     print*,i,tmo_eig(i)
   enddo
   print*,"Projectile HF Energy =", p_ehf
   do i = 1, pmo_num
     print*,i,pmo_eig(i)
   enddo
!   print*,"tgrid = ", tgrid
!   print*,"zgrid = ", zgrid
   print*, "" 

   ntime= n_time

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
  nstatest = 2
  touch nstatest
  nsta = nstatest
    n_tsta = 1
    n_tpsta = 0
    n_psta = 1
  vproj = v_coll

  allocate(mcoup(ntime,nsta,nsta),timegrid(ntime),esta(nsta))
  allocate(psi(nsta))
  allocate(rmat2intrp(ntime,nsta,nsta),cmat2intrp(ntime,nsta,nsta))
  allocate(matintrp(nsta,nsta))

  timegrid(1:ntime) = tgrid(1:n_time)
  esta=0.0d0
  !do i = 1, tmo_num
  !  esta(i) = tmo_eig(i)
  !enddo
  !do i = 1, pmo_num
  !  esta(i+tmo_num) = pmo_eig(i)
  !enddo

  !!tmo_eig(2) = -0.078925734200000! -0.237691366041379
  tmo_eig(2) = -0.237691366041379
  pmo_eig(1) = -0.499859669283984

  esta(1) = 0.0d0
  esta(2) = 0.0d0 !!!- tmo_eig(2) + pmo_eig(1)


!!!

  print*,'start dyn'
  do ib = 1, n_bimp

    b_coll = bgrid(ib)
    ib_coll = ib
    TOUCH b_coll

    mcoup(:,:,:) = dcmplx(0.0D0,0.0D0)
    do it = 1, ntime !! for ZGrid
      mcoup(it,1:nsta,1:nsta) = coll_couplings_nicotest(1:nsta,1:nsta,it)
      !print*,mcoup(it,1,2),'nico'
    end do
    FREE coll_couplings_nicotest
    
!!!!!!!!!!!!!!!!!
   if(ib==1) then
     open(unit=20,file='Prop_collision.out')
      write(20,*) n_bimp, nsta
      do i = 1, nsta
        write(*,*)"esta:",i,real(esta(i))
        write(20,*)real(esta(i))
      enddo
     close(20)
   endif
!!!!!!!!!!!!!!!!

    psi(:) = 0.0d0
    psi(1) = 1.0d0 ! only true for CIS, otherwise use lines below (Nico 22.04.2024)

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

end program cippres_prop_collision_nicotest
