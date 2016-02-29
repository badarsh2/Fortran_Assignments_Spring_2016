program trace
		implicit none
	! To calculate the trace of a 3x3 matrix by reading from input_matrix.txt and using a subroutine
	INTEGER :: A(3,3), trace1=0, i

	! Unit just is an ID reference
	open(unit=10, file='input_matrix.txt')
	do i = 1,3
		read(10,*) A(i,1),A(i,2),A(i,3)
	end do

	! Calling the subroutine for calculaing trace
	call tracecalc(A, trace1)
end program

! Subroutine definition
! (Subroutines are similar to functions in C language)
subroutine tracecalc(B, tr)
		implicit none
	integer, intent(in) :: B(3,3)
	integer, intent(out) :: tr
	integer :: i

	tr = 0
	do i = 1,3
		tr = tr + B(i,i)
	end do
	print *, "The trace calculated using subroutine is ", tr

end subroutine

