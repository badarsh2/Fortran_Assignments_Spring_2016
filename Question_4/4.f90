! Program to solve Linear equations in n variables
! PLEASE GO THROUGH THE INPUT_MATRIX_2.txt FILE SO THAT U CAN UNDERSTAND THE CODE BETTER
program dynamo
	implicit none
	! Dynamic allocaiton
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: A, inv, B, C
	integer :: siz,i,ok,colsiz, trace

	open(unit = 10, file = 'input_matrix_4.txt')
	! Reading size of the system
	read(10,*) siz
	! Ax = B
	ALLOCATE(A(siz,siz),STAT=ok)
	ALLOCATE(B(siz,1),STAT=ok)
	ALLOCATE(C(siz,1),STAT=ok)
	ALLOCATE(inv(siz,siz),STAT=ok)

	! Reading matrix A
	print *, 'Matrix A in Ax = B is'
	do i = 1,siz
		read(10,*) A(i,1:siz)
		print *, A(i,1:siz)
	end do

	! Reading matrix B
	print *, 'Matrix B in Ax = B is'
	do i = 1,siz
		read(10,*) B(i,1)
		print *, B(i,1)
	end do

	! print *, trans
	call inverse(A, inv, siz)
	call multiply(inv, B,  C, siz)
end program

subroutine inverse(a, inv, n)
	integer, intent(in) :: n
	DOUBLE PRECISION, intent(in) :: a(n,n)
	DOUBLE PRECISION, intent(out) :: inv(n,n)
	DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: aug
	DOUBLE PRECISION :: tempvar
	integer :: ok, i, j, k
	ALLOCATE(aug(2*n,2*n),STAT=ok)
	do i = 1,2*n
		aug(i,1:2*n) = 0.0
	end do
	do i = 1,n
		aug(i,1:n) = a(i,1:n)
	end do
	do i = 1, n
		do j = 1, 2*n
			if (j==i+n) then
				aug(i,j) = 1
			end if
		end do
	end do
	do i = n,2,-1
		if (aug(i-1,1) < aug(i,1)) then
			do j = 1, 2*n
				tempvar = aug(i,j)
				aug(i,j) = aug(i-1, j)
				aug(i-1,j) = tempvar
			end do
		end if
	end do

	do i = 1,n
		! print *, aug(i,1:2*n)
	end do

	tempvar = 0.0

	do i = 1,n
		do j = 1, 2*n
			if (j == i) then
				! print *, "Lol"
			else
				tempvar = aug(j,i)/aug(i,i)
				! print *, "Tempvar is", tempvar
				do k = 1, 2*n
					aug(j,k) = aug(j,k) - (aug(i,k) * tempvar)
					! print *, i, j, k, aug(j,k)
				end do
				tempvar = 0.0
			end if
		end do
	end do

	do i = 1,n
		tempvar = aug(i,i)
		do j = 1, 2*n
			aug(i,j) = aug(i,j) / tempvar
		end do
	end do

	inv(1:n,1:n) = aug(1:n, n+1:2*n)
end subroutine

subroutine multiply(a, b, c, n)
	integer, intent(in) :: n
	DOUBLE PRECISION, intent(in) :: a(n,n)
	DOUBLE PRECISION, intent(in) :: b(n,1)
	DOUBLE PRECISION, intent(out) :: c(n,1)
	integer ::  i, k
	do i=1, n
		do k=1, n
			c(i,1) = c(i,1) + (a(i,k)*b(k,1))
		end do
	end do
	print *, 'Approx value of B is:'
	do i = 1, n
		print *, c(i,1)
	end do
end subroutine
					





