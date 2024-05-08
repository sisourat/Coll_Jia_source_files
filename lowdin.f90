!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function compute_det(n,mat)
implicit none

integer :: n, lda, info, nperm
double complex, dimension(n,n) :: mat
double complex :: det
double complex :: compute_det
integer, dimension(n) :: ipiv
integer :: i

 if (n==0) then
  compute_det = 1d0 !nico
  !compute_det = 0d0 ! ccjia test
  return
 endif 
lda = n
call ZGETRF(n,n,mat,lda,ipiv,info)

! number of perm
nperm=0
det=1d0
do i = 1, n
 if(ipiv(i)/=i) nperm=nperm+1
 det=det*mat(i,i)
enddo
compute_det = (-1d0)**nperm*(det)
return
end function compute_det

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine lowdin(ne,nea,neb,nmo,ovmo,h1emo,r12mo,cdet1,cdet2,deta1,deta2,detb1,detb2,ov,h1e,r12)
subroutine lowdin(ne,nea,neb,nmo,ovmo,h1emo,r12mo,deta1,deta2,detb1,detb2,ov,h1e,r12)
implicit none

!! compute the matrix elements between two Slater determinants (det1 and
!det2) formed with non-orthogonal MOs
!! for efficiency, the matrix elements in the MO basis (ovmo, h1emo and
!r12mo) are global variables
!! ne is the number of electrons
!! ov, h1e, r12 are the matrix elements between the two determinants
!det1 and det2

!! FOR EFFICIENCY, r12mo SHOULD BE CHANGED AS GLOBAL VARIABLES (AVOID
!COPYING LARGE ARRAYS)

integer, intent(in) :: ne, nea, neb, nmo
double complex, dimension(nmo,nmo), intent(in) :: ovmo, h1emo
double complex, dimension(nmo,nmo,nmo,nmo), intent(in) :: r12mo

!double precision, intent(in) :: cdet1, cdet2
integer, dimension(nea), intent(in) :: deta1, deta2
integer, dimension(neb), intent(in) :: detb1, detb2

double complex, intent(out) :: ov, h1e, r12

double complex, dimension(ne) :: vectemp
double complex, dimension(ne,ne) :: ovmat, ovstore
double complex, dimension(ne-1,ne-1) :: comat
double complex, dimension(ne-2,ne-2) :: comat2

double complex, external :: compute_det

integer :: i, j , k, l, ia, ib, ja, jb, ka, kb, la, lb

!! compute and store overlap matrix between the determinants
ovmat(:,:) = dcmplx(0.0D0,0.0D0)
do i = 1, nea
 do j = 1, nea
  ja =  deta2(j)
  ia =  deta1(i)
  !write(*,*) "ccjia: ja,ia",ja,ia
  ovmat(j,i) = ovmo(ja,ia)
 enddo
enddo

do i = 1, neb
 do j = 1, neb
  jb =  detb2(j)
  ib =  detb1(i)
  ovmat(j+nea,i+nea) = ovmo(jb,ib)
 enddo
enddo

ovstore(:,:) = ovmat(:,:)

!! compute the overlap between the two determinants
ov = compute_det(ne,ovmat)

!! compute the 1e hamiltonian between the two determinants
h1e = dcmplx(0.0D0,0.0D0)
do i = 1, nea
 do j = 1, nea
  ovmat(:,:) = ovstore(:,:)
  ovmat(j,:) = ovmat(ne,:)
  ovmat(:,i) = ovmat(:,ne)
  comat(:,:) = ovmat(1:ne-1,1:ne-1)
 
  ja =  deta2(j)
  ia =  deta1(i)
  h1e = h1e + (-1d0)**(i+j)*h1emo(ja,ia)*compute_det(ne-1,comat)
!  write(*,*)ja,ia,h1emo(ja,ia)
!  write(*,*)compute_det(ne-1,comat)
!  write(*,*)"h1e a",h1e
 enddo
enddo

do i = 1, neb
 do j = 1, neb
  ovmat(:,:) = ovstore(:,:)
  ovmat(j+nea,:) = ovmat(ne,:)
  ovmat(:,i+nea) = ovmat(:,ne)
  comat(:,:) = ovmat(1:ne-1,1:ne-1)

  jb =  detb2(j)
  ib =  detb1(i)
  h1e = h1e + (-1d0)**(nea+i+nea+j)*h1emo(jb,ib)*compute_det(ne-1,comat)
!  write(*,*)"h1e b",h1e
 enddo
enddo

!! compute the 2e repulsion between the two determinants
r12 = dcmplx(0.0D0,0.0D0)
!! alpha el for e1 and e2
do i = 1, nea  !e1
 do j = 1, nea  !e1
  ja =  deta2(j)
  ia =  deta1(i)

   do k = i+1, nea  !e2
     do l = j+1, nea  !e2
       la =  deta2(l)
       ka =  deta1(k)

       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j) = ovmat(:,ne)
       ovmat(:,i) = ovmat(:,ne-1)
       ovmat(l,:) = ovmat(ne,:)
       ovmat(k,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

!nico indices to be checked
       r12 = r12 + r12mo(la,ka,ja,ia)*(-1d0)**(i+j+k+l)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

!! alpha el for e1 and beta el for e2
do i = 1, nea  !e1
 do j = 1, nea  !e1
  ja =  deta2(j)
  ia =  deta1(i)

   do k = 1, neb  !e2
     do l = 1, neb  !e2
       lb =  detb2(l)
       kb =  detb1(k)
       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j) = ovmat(:,ne)
       ovmat(:,i) = ovmat(:,ne-1)
       ovmat(l+nea,:) = ovmat(ne,:)
       ovmat(k+nea,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

!nico indices checked
!nico       r12 = r12 + r12mo(lb,kb,ja,ia)*(-1d0)**(i+j+k+l+nea+nea)*compute_det(ne-2,comat2)
       r12 = r12 + r12mo(ia,kb,ja,lb)*(-1d0)**(i+j+k+l+nea+nea)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

!! beta el for e1 and alpha el for e2
!do i = 1, neb  !e1
! do j = 1, neb  !e1
!  jb =  detb2(j)
!  ib =  detb1(i)
!
!   do k = 1, nea  !e2
!     do l = 1, nea  !e2
!       la =  deta2(l)
!       ka =  deta1(k)
!       ovmat(:,:) = ovstore(:,:)
!       ovmat(:,j+nea) = ovmat(:,ne)
!       ovmat(:,i+nea) = ovmat(:,ne-1)
!       ovmat(l,:) = ovmat(ne,:)
!       ovmat(k,:) = ovmat(ne-1,:)
!       comat2(:,:) = ovmat(1:ne-2,1:ne-2)
!
!       r12 = r12 + r12mo(la,ka,jb,ib)*(-1d0)**(i+j+k+l+nea+nea)*compute_det(ne-2,comat2)
!
!     enddo
!   enddo
!
! enddo
!enddo

!! beta el for e1 and beta el for e2
do i = 1, neb  !e1
 do j = 1, neb  !e1
  jb =  detb2(j)
  ib =  detb1(i)

   do k = i+1, neb  !e2
     do l = j+1, neb  !e2
       lb =  detb2(l)
       kb =  detb1(k)
       ovmat(:,:) = ovstore(:,:)
       ovmat(:,j+nea) = ovmat(:,ne)
       ovmat(:,i+nea) = ovmat(:,ne-1)
       ovmat(l+nea,:) = ovmat(ne,:)
       ovmat(k+nea,:) = ovmat(ne-1,:)
       comat2(:,:) = ovmat(1:ne-2,1:ne-2)

!nico indices to be checked
      r12 = r12 + r12mo(lb,kb,jb,ib)*(-1d0)**(i+j+k+l+2*nea+2*nea)*compute_det(ne-2,comat2)

     enddo
   enddo

 enddo
enddo

end subroutine lowdin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
