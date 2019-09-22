program ConvDifReac
! This program solves the Convection-Diffusion-Reaction equation
! with constant velocity flux and constant Diffusion coefficient
use constants_module
use subroutines_module
implicit none
integer :: i, j, n_cells, n_t, NZE
real (kind=dp), dimension(:), allocatable :: a1, b1, c1, vec_der, C
real (kind=dp), dimension(:), allocatable :: val, x
integer, dimension(:), allocatable :: col_ind, row_ptr
real (kind=dp) :: delta_t, delta_x, tf
real (kind=dp) :: D, u, k, Cref, L, theta
real (kind=dp) :: S, CFL

! OPENING FILES
open (unit=300, file='inputs.dat', status='old')
open (unit=305, file='log.dat')
open (unit=310, file='results.dat')

! READING INPUTS
read (unit=300, fmt=*) D
read (unit=300, fmt=*) u
read (unit=300, fmt=*) k
read (unit=300, fmt=*) Cref
read (unit=300, fmt=*) L
read (unit=300, fmt=*) n_cells
read (unit=300, fmt=*) delta_t
read (unit=300, fmt=*) tf
read (unit=300, fmt=*) theta


! CALCULATING SOME THINGS...
 delta_x = L/dble(n_cells)

 S = (D*delta_t)/(delta_x**2)
 CFL = (u*delta_t)/delta_x

 NZE = 3*n_cells - 2

 n_t = floor(tf/delta_t) + 1

! ALLOCATING MEMORY
allocate(a1(1:n_cells))
allocate(b1(1:n_cells))
allocate(c1(1:n_cells))
allocate(vec_der(1:n_cells))
allocate(C(1:n_cells))
allocate(val(1:NZE))
allocate(col_ind(1:NZE))
allocate(row_ptr(1:n_cells+1))
allocate(x(1:n_cells))

! CREATING THE x VECTOR
x(1) = delta_x/2
do i = 2, n_cells
	x(i) = x(i-1) + delta_x
end do

! WRITING LOG
write (unit=305, fmt=100) 'delta_y =', delta_x
write (unit=305, fmt=100) 'delta_t =', delta_t
write (unit=305, fmt=100) '      S =', S
write (unit=305, fmt=100) '    CFL =', CFL
100 format(A10, 1x, f16.11)

! MAKING THE MATRICES
 a1 = -((1.0D0-theta)*S + (1.0D0-theta)*0.5D0*CFL)
 a1(1) = 0.0D0
 b1 = 1.0D0 + (1.0D0-theta)*delta_t*k + 2.0D0*(1.0D0-theta)*S
 b1(1) = b1(1) - ((1.0D0-theta)*S + (1.0D0-theta)*0.5D0*CFL)*(1.0D0 - delta_x*u/D) !by Robin BC
 b1(n_cells) = b1(n_cells) - (1.0D0-theta)*S + (1.0D0-theta)*0.5D0*CFL !by Neumann BC
 c1 = -(1.0D0-theta)*S + (1.0D0-theta)*0.5D0*CFL
 c1(n_cells) = 0.0D0

col_ind(1) = 1
col_ind(2) = 2
do i = 1, n_cells-2
	col_ind(3*i) = i
	col_ind(3*i+1) = i+1
	col_ind(3*i+2) = i+2
end do
col_ind(NZE-1) = n_cells-1
col_ind(NZE) = n_cells

row_ptr(1) = 0
row_ptr(2) = 2
do i = 3, n_cells
	row_ptr(i) = row_ptr(i-1) + 3
end do
row_ptr(n_cells+1) = row_ptr(n_cells) + 2

val(1) = 1.0D0 - theta*delta_t*k - 2.0D0*theta*S + (theta*S + theta*0.5D0*CFL)*(1-delta_x*u/D) !by Robin BC
val(2) = theta*S - theta*0.5D0*CFL
do i = 1, n_cells-2
	val(3*i) = theta*S + theta*0.5D0*CFL
	val(3*i+1) = 1.0D0 - theta*delta_t*k - 2.0D0*theta*S
	val(3*i+2) = theta*S - theta*0.5D0*CFL
end do
val(NZE-1) = theta*S + theta*0.5D0*CFL
val(NZE) = 1.0D0 - theta*delta_t*k - 2.0D0*theta*S + theta*S - theta*0.5D0*CFL !by Neumann BC

! SOLVING THE SYSTEM

 C = 0.008D0 !initial conditions for C

write (unit=310, fmt='(A13)') '# Delta No. 0'	
write (unit=310, fmt='(A7)') '# t = 0'
do j = 1, n_cells
	write (unit=310, fmt=105) x(j), C(j)
end do
write (unit=310, fmt=*)
write (unit=310, fmt=*)

do i = 1, n_t

	call CRSmatvec(n_cells, n_cells, NZE, val, col_ind, row_ptr, C, vec_der)

	! here we add the only component of the BC's vector (due to the Robin BC at the first boundary)
	vec_der(1) = vec_der(1) + (theta*S + theta*0.5D0*CFL)*(Cref*delta_x*u/D) + &
	((1.0D0-theta)*S + (1.0D0-theta)*0.5D0*CFL)*(Cref*delta_x*u/D)

	call TDMA(n_cells, a1, b1, c1, vec_der, C)

	write (unit=310, fmt='(A12, i5)') '# Delta No. ', i	
	write (unit=310, fmt='(A5, f15.12)') '# t =', i*delta_t
	do j = 1, n_cells
		write (unit=310, fmt=105) x(j), C(j)
		105 format(f18.14, 3x, f18.14)
	end do
	write (unit=310, fmt=*)
	write (unit=310, fmt=*)

end do


deallocate(a1)
deallocate(b1)
deallocate(c1)
deallocate(vec_der)
deallocate(C)
deallocate(val)
deallocate(col_ind)
deallocate(row_ptr)
deallocate(x)

end program ConvDifReac
