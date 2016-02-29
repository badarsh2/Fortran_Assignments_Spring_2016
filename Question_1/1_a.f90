program trace
		implicit none
	! To calculate the trace of a 3x3 matrix by hardcoding
	INTEGER :: A(3,3), B(3,3), C(3,3), trace1=0, sumofsquares = 0, i, j, determinant

	! Method 1
	A = RESHAPE((/3, 4, 8, 1, 5, 6, 9, 4, 7/),(/3, 3/))

	! Method 2
	A(1,1) = 3
	A(1,2) = 4
	A(1,3) = 8
	A(2,1) = 1
	A(2,2) = 5
	A(2,3) = 6
	A(3,1) = 9
	A(3,2) = 4
	A(3,3) = 7

	trace1 = A(1,1) + A(2,2) + A(3,3)
	print *, "The trace is ", trace1
end program

