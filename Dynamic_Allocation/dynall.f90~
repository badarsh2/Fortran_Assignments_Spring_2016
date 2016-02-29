program dynamo
	implicit none
! Dynamic allocaiton
INTEGER, ALLOCATABLE, DIMENSION(:,:):: A, trans
integer :: siz,i,ok

open(unit = 10, file = 'input_matrix_2.txt')
read(10,*) siz
ALLOCATE(A(siz,siz),STAT=ok)
ALLOCATE(trans(siz,siz),STAT=ok)

print *, 'Input matrix read from file is'
do i = 1,siz
	read(10,*) A(i,1:siz)
	print *, A(i,1:siz)
end do

call transposecalc(A,trans)
! print *, trans 

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
