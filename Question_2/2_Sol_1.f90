program dynamo
	implicit none
! Dynamic allocaiton
INTEGER, ALLOCATABLE, DIMENSION(:,:):: A, trans
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: inv
integer :: siz,i,ok,colsiz, trace

open(unit = 10, file = 'input_matrix_2.txt')
read(10,*) siz, colsiz
if (siz == colsiz) then
	ALLOCATE(A(siz,siz),STAT=ok)
	ALLOCATE(trans(siz,siz),STAT=ok)

	print *, 'Input matrix read from file is'
	do i = 1,siz
		read(10,*) A(i,1:siz)
		print *, A(i,1:siz)
	end do

	call transposecalc(A,trans)
	call tracecalc(A,trace,siz)
	! print *, trans
else
	print *, 'Input matrix is not a square matrix. Please check the input file'
end if
end program

subroutine transposecalc(B, C)
		implicit none
	integer, intent(in) :: B(3,3)
	integer, intent(out) :: C(3,3)
	integer :: i, j

	C = 0
	do i = 1,3
		do j = 1,3
			C(i,j) = B(j,i)
		end do
	end do

	print *, 'Transpose calculated using a subroutine is'
	do i = 1 , 3
		print *, C(i,1), C(i,2), C(i,3)
	end do

end subroutine

subroutine tracecalc(B, tr, siz)
		implicit none
	integer, intent(in) :: siz
	integer, intent(in) :: B(siz,siz)
	integer, intent(out) :: tr
	integer :: i

	tr = 0
	do i = 1,siz
		tr = tr + B(i,i)
	end do
	print *, "The trace calculated using subroutine is ", tr

end subroutine

subroutine detcalc(B, det, siz)
		implicit none
	integer, intent(in) :: siz
	integer, intent(in) :: B(siz,siz)
	integer, intent(out) :: det
	integer :: i, j, k, sumdet = 0
	integer :: epsilon(3,3,3)

	epsilon = 0
	epsilon(1,2,3) = 1
	epsilon(2,3,1) = 1
	epsilon(3,1,2) = 1
	epsilon(1,3,2) = -1
	epsilon(2,1,3) = -1
	epsilon(3,2,1) = -1
	do i = 1,3
		do j = 1,3
			do k = 1, 3
				sumdet = sumdet + epsilon(i,j,k)*B(i,1)*B(j,2)*B(k,3)
			end do
		end do
	end do

	print *, 'Determinant of the matrix is', sumdet

end subroutine

subroutine inverse(a,c,n)
implicit none 
integer n, a(n,n)
double precision c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

