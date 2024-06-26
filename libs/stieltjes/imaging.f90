subroutine imaging(npt,e_point,g_point,nmax,maxord,eorder,gorder)
implicit none

 integer :: nmax

 integer, parameter :: QR_K = selected_real_kind (32) ! use quadruple precision
 double precision, parameter :: overmax=1d0 ! set arbitrary to 1 to determine the maximum order of pol.

 integer, intent(in) :: npt
 real (kind=QR_K), dimension(npt), intent(in) :: e_point, g_point

 real (kind=QR_K), dimension(0:nmax,npt) :: qpol
 real (kind=QR_K), dimension(nmax) :: acoef
 real (kind=QR_K), dimension(0:nmax) :: bcoef

 real (kind=QR_K), dimension(nmax) :: diag, offdiag
 real (kind=QR_K), dimension(nmax,nmax) :: abvec
 real (kind=QR_K) :: asum, bprod, qnorm, qoverlap
 integer :: iord, ierr, maxord, min, max

 real (kind=QR_K), dimension(nmax) :: enew, gnew
 real (kind=QR_K), dimension(nmax,nmax) :: eorder, gorder

 integer :: i, j

! initiate the recursive computation of the a,b coefficients and the orthogonal 
! polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,npt
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,npt
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,npt
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

! calculate the higher-order coefficients and polynomials recursively
! up to the (NMAX-1)th order (total of NMAX polynomials)

       asum=acoef(1)
       do i=3,nmax

          asum=asum+acoef(i-1)

          do j=1,npt
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)-bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,npt
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)/(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,npt
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

       end do

! calculate the nmax-th order polynomial just for the orthogonality check 
       do j=1,npt
          qpol(nmax,j)=(1.q0/e_point(j)-acoef(nmax))*qpol(nmax-1,j)-bcoef(nmax-1)*qpol(nmax-2,j)
       end do

!       write(*,*)"a",acoef(1:nmax)
!       write(*,*)"b",bcoef(0:nmax)
!       write(*,*)"qpol",qpol(:,:)
!       stop

! check the orthogonality of the polynomials to define the maximal approximation order 
! if the orthogonality is preserved for all orders, MAXORD is set to NMAX
       maxord=nmax
       qnorm=bcoef(0)
       do i=1,nmax
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,npt
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
          if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          
          if (qnorm/qabs(qoverlap).le.overmax) then
! MAXORD=I-1 is appropriate since the polynomial failing 
! the orthogonality check should not be used
             maxord=i-1
             print*, ' MAXORD=',maxord
             go to 10
          end if
       end do
 10    continue
!!       maxord=nmax

! look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:'
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
          print*, ' MAXORD=',maxord
       end if

! perform the gamma calculation using the successive approximations 
! n=5,...,nmax

   do iord=5,maxord

     write(*,*)"Performs Stieltjes at order",iord 

! fill the coefficients matrix
       do i=1,iord
          diag(i)=acoef(i)
!          write(*,*)"diag",i,diag(i)
       end do
       do i=2,iord
          offdiag(i)=-qsqrt(bcoef(i-1))
!          write(*,*)"offdiag",i,offdiag(i)
       end do
!stop

! diagonalize the coefficients matrix
! initialize the arrays
       do i=1,nmax
          do j=1,nmax
             abvec(i,j)=0.q0
          end do
          abvec(i,i)=1.q0
       end do
       call tql2(iord,iord,diag(1:iord),offdiag(1:iord),abvec(1:iord,1:iord),ierr)
       if (ierr.ne.0) then
          print*, '***WARNING*** Stieltjes:'
          print*, ' the eigenvalue no. ',ierr,' failed to converge'
       end if

! fill the Stieltjes energy and gamma arrays
! note that the eigenvalues are inverse energies and are given in ascending order 

!       do i = 1, iord
!        print*,diag(i)
!        print*,i,abvec(i,1:iord)
!        print*,""
!       enddo
!       stop
       do i=1,iord
!          print*,"nico",i,diag(iord+1-i),abvec(1,iord+1-i)**2
          enew(i)=1.q0/diag(iord+1-i)
          gnew(i)=bcoef(0)*abvec(1,iord+1-i)**2
!          print*,'new',i,enew(i),gnew(i)
       end do
!          print*,''
!       stop

       call eigsrtnico(enew,gnew,iord,nmax)

! calculate the gamma's by simple numerical differentiation at the middle 
! point of each [ENEW(I),ENEW(I+1)] interval
       do i=1,iord-1
          eorder(i,iord)=0.5d0*(enew(i)+enew(i+1))
          gorder(i,iord)=0.5d0*(gnew(i+1)+gnew(i))/(enew(i+1)-enew(i))
       end do

   enddo ! loop over iord  


end subroutine
