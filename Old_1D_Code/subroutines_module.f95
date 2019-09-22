module subroutines_module

use constants_module

contains

	subroutine CRSmatvec(n, m, NZE, val, col_ind, row_ptr, b, R)
		!This subroutine performs the multiplication A*b, where 'A' is an n by m matrix
		!stored in Compressed Row Storage (CRS) format and 'b' is a full m by 1 vector.
		!The input arguments are:
		!	n   --> number of rows of A
		!	m   --> number of columns of A
		!	NZE --> number of non-zero entries in A
		!	val, col_ind & row_ptr --> arrays that define the matrix A in CRS format
		!	b   --> right vector for the multiplication A*b
		!	R   --> Variable to store the result of the multiplication (size n by 1)

		implicit none
		integer, intent (in) :: n, m, NZE
		real (kind=dp), dimension(1:NZE), intent (in) :: val
		integer, dimension(1:NZE), intent (in) :: col_ind
		integer, dimension(1:n+1), intent (in) :: row_ptr
		real (kind=dp), dimension(1:m), intent (in) :: b
		real (kind=dp), dimension(1:n), intent (out) :: R
		integer :: i, j, k

		R = 0.0D0
		k = 1

		do i = 1, n
			do j = 1, (row_ptr(i+1)-row_ptr(i))
				R(i) = R(i) + val(k)*b(col_ind(k))
				k = k + 1
			end do
		end do
	end subroutine CRSmatvec


	subroutine Modified_TDMA(N, a, b, c, d, R)
		!a
		implicit none
		integer, intent (in) :: N
		real (kind=dp), dimension(1:N), intent(in) :: a, b, c, d
		real (kind=dp), dimension(1:N), intent(out) :: R
		real (kind=dp), dimension(1:N) :: bp, cp, dprime, u, v, q, y
		real (kind=dp) :: temp
		integer :: i

		bp = b; cp = c; dprime = d
		u = 0.0D0; v = 0.0D0

		u(N) = cp(N)
		u(1) = -bp(1)

		v(N) = -a(1)/bp(1)
		v(1) = 1.0D0

		bp(N) = bp(N) + (a(1)*cp(N))/bp(1)
		bp(1) = 2.0D0*bp(1)

		cp(1) = cp(1)/bp(1)
		dprime(1) = dprime(1)/bp(1)
		u(1) = u(1)/bp(1)

		do i = 2, N, 1
			temp = bp(i) - a(i) * cp(i-1)
			if (i<N) then
				cp(i) = cp(i)/temp
			end if
		dprime(i) = (dprime(i) - a(i)*dprime(i-1))/temp
		u(i) = (u(i) - a(i)*u(i-1))/temp
		end do

		y(N) = dprime(N)
		q(N) = u(N)
		do i = N-1, 1, -1
			y(i) = dprime(i) - cp(i)*y(i+1)
			q(i) = u(i) - cp(i)*q(i+1)
		end do

		R = y - (dot_product(v,y)/(1+dot_product(v,q)))*q

	end subroutine Modified_TDMA


	subroutine TDMA(N, a, b, c, d, R)
		!a
		implicit none
		integer, intent (in) :: N
		real (kind=dp), dimension(1:N), intent(in) :: a, b, c, d
		real (kind=dp), dimension(1:N), intent(out) :: R
		real (kind=dp), dimension(1:N) :: cp, dprime
		real (kind=dp) :: temp
		integer :: i

		dprime(1) = d(1)/b(1)
		cp(1) = c(1)/b(1)

		do i = 2, N
			temp = b(i) - a(i) * cp(i-1)
			if (i<N) then
				cp(i) = c(i)/temp
			end if
			dprime(i) = (d(i) - a(i)*dprime(i-1))/temp
		end do

		R(N) = dprime(N)

		do i = n-1, 1, -1
			R(i) = dprime(i) - cp(i)*R(i+1)
		end do

	end subroutine TDMA


end module subroutines_module
