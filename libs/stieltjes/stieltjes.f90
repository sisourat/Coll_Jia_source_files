program stieltjes
use interpolation
implicit none

integer :: npt
integer, parameter :: QR_K = 16 !selected_real_kind (32)
real (kind=QR_K), allocatable, dimension(:) :: e, g

integer :: nmin=5, nmax = 25 ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
real (kind=QR_K), allocatable, dimension(:) :: sk, gord
real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1
real (kind=QR_K), allocatable, dimension(:) :: eallord,gallord
real (kind=QR_K) :: shift1 = 0.1d0
integer :: imax1, imax2
integer :: inmax, ishift
integer :: i, j, k, ichan, ijob
integer :: exit_cycle
character(len=60) :: fname

! Tsveta
real (kind=QR_K) :: g_
real (kind=QR_K) :: temp, pi, gav, stadev, gtot

pi = dacos(-1d0)

! reads and sorts energy and matrix elements

read(*,*)ijob
read(*,*)npt
allocate(e(npt),g(npt))

do i = 1, npt
  read(*,*)e(i),g(i)
enddo
!g(:) = abs(g(:))**2

if(sum(g)==0d0) then
 write(*,*)"All matrix elements are zero, I stop"
 stop
endif

call sort2(npt,e,g)
open(unit=10,file='input.sorted.txt')
do i = 1, npt
  write(10,'(2(f20.16,1X))')e(i),g(i)
enddo
close(10)
shift1 = 0d0
if(ijob/=2) then
 ishift = 2
 shift1 = ishift * 0.1d0
 shift1 = shift1 + abs(e(1))
endif
!shift1=35d0
!print*,'SHIFT=',shift1
e(:) = e(:) + shift1
allocate(e1(nmax,nmax), g1(nmax,nmax))       
call imaging(npt,e,g,nmax,nmax,e1,g1)
e(:) = e(:) - shift1

open(unit=10,file='stieltjes.info.txt')
 write(10,*)nmin,nmax,shift1
close(10)

k =1
do i = nmin+2, nmax
  k = k + 1
  write(fname, '(A16,I2,A4)')"stieltjes.order.",i,".txt"
  write(*, '(A16,I2,A4)')"stieltjes.order.",i,".txt"
  open(238,file=fname)
  do  j = nmin, i-1
     write(238,'(2(f20.10,1X))')e1(j,i)-shift1,g1(j,i)
  enddo
  if(ijob==1) then
    call interp(e1(nmin:i-1,i)-shift1,g1(nmin:i-1,i),i-nmin,0q0,g_)
    g_ = 2.0d0*pi*g_
    write(*,'(I3,A3,F28.15,A)')i," ",g_*27211,' in meV'
  endif
  close(238)
enddo

deallocate(e1,g1)
deallocate(e,g)

end program stieltjes
