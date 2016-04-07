! Program to multiply two matrices
! PLEASE GO THROUGH THE INPUT_MATRIX_2.txt FILE SO THAT U CAN UNDERSTAND THE CODE BETTER

program matmult
	implicit none
! Dynamic allocaiton
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: A, B, C
integer :: A_rowsize, A_colsize, B_rowsize, B_colsize, i, j, k, ok

! Reading matrix A row and column size
open(unit = 10, file = 'input_matrix_2.txt')
read(10,*) A_rowsize, A_colsize
ALLOCATE(A(A_rowsize, A_colsize),STAT=ok)

! Reading elements of Matrix A
print *, 'Matrix A:'
do i = 1, A_rowsize
	read(10,*) A(i,1:A_colsize)
	print *, A(i,1:A_colsize)
end do

! Reading matrix B row and column size
read(10,*) B_rowsize, B_colsize
ALLOCATE(B(B_rowsize, B_colsize),STAT=ok)

! Reading elements of Matrix B
print *, 'Matrix B:'
do i = 1, B_rowsize
	read(10,*) B(i,1:B_colsize)
	print *, B(i,1:B_colsize)
end do

! Matrices can be multiplied only if A's colsize and B's rowsize are equal
if (A_colsize == B_rowsize) then
	ALLOCATE(C(A_rowsize, B_colsize),STAT=ok)
	do i=1, A_rowsize
		do j=1, B_colsize
			C(i,j) = 0
			do k=1, A_colsize
				C(i,j) = C(i,j) + (A(i,k)*B(k,j))
			end do
		end do
	end do
	print *, 'Matrix multiplication result:'
	do i = 1, A_rowsize
		print *, C(i,1:B_colsize)
	end do
else
	print *, 'Matrices cannot be multiplied'
end if
! print *, trans
end program
