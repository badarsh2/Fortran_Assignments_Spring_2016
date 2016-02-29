program trace
		implicit none
	! To calculate the trace of a 3x3 matrix by reading from input_matrix.txt and using DO loops
	INTEGER :: A(3,3), trace1=0, i, j

	! Unit just is an ID reference
	open(unit=10, file='input_matrix.txt')
	do i = 1,3
		read(10,*) A(i,1),A(i,2),A(i,3)
	end do

	! DO loop for iterating and calculaing trace
	do i=1,3
		trace1 = trace1 + A(i,i)
	end do

	print *, "The trace calculated using DO LOOP is ", trace1
end program
