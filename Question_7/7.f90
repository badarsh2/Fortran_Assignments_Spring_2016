! Program to multiply two matrices
! PLEASE GO THROUGH THE INPUT_MATRIX_2.txt FILE SO THAT U CAN UNDERSTAND THE CODE BETTER

program matmult
	implicit none
! Dynamic allocaiton
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:):: A, B, C
integer :: A_rowsize, A_colsize, B_rowsize, B_colsize, i, j, k, l, ok, delta(3,3), bigdelta(3,3,3,3) = 0
double precision :: C11, C12, C44, stiffness(3,3,3,3) = 0.0, sigma(3,3) = 0.0
double precision :: aux

delta = RESHAPE((/1, 0, 0, 0, 1, 0, 0, 0, 1/),(/3, 3/))
bigdelta(1,1,1,1) = 1
bigdelta(2,2,2,2) = 1
bigdelta(3,3,3,3) = 1

! Reading matrix A row and column size
open(unit = 10, file = 'input_matrix_2.txt')
read(10,*) C11, C12, C44

print *, "C11 is"
print *, C11
print *, "C12 is"
print *, C12
print *, "C44 is"
print *, C44

! Reading matrix B row and column size
read(10,*) B_rowsize, B_colsize
ALLOCATE(B(B_rowsize, B_colsize),STAT=ok)

! Reading elements of Matrix B
do i = 1, B_rowsize
	read(10,*) B(i,1:B_colsize)
end do

do i = 1,3
	do j = 1,3
		do k = 1,3
			do l = 1,3
				stiffness(i,j,k,l) = stiffness(i,j,k,l) + C12*delta(i,j)*delta(k,l)
				stiffness(i,j,k,l) = stiffness(i,j,k,l) + (C44*delta(i,k)*delta(j,l)) 
                                stiffness(i,j,k,l) = stiffness(i,j,k,l) + (C44*delta(i,l)*delta(j,k)) 
                                stiffness(i,j,k,l) = stiffness(i,j,k,l) + (C11 - C12 - 2*C44)*bigdelta(i,j,k,l)
			end do
		end do
	end do
end do

print *, 'Strain Tensor:'
do i = 1, B_rowsize
	print *, B(i,1:B_colsize)
end do

do i = 1,3
	do j = 1,3
		aux = 0.0
		do k = 1,3
			do l = 1,3
				! print *,stiffness(i,j,k,l)
				aux = aux + stiffness(i,j,k,l)*B(k,l) 
				! print *, i,j,k,l,aux
			end do
		end do
		sigma(i,j) =aux
	end do
end do

print *, 'The Stress Tensor (all values in GPa) is:'
	do i = 1, 3
		print *, sigma(i,1:3)
	end do
! print *, trans
end program
