module constants_module
implicit none

integer, parameter :: sp = selected_real_kind( 6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qd = selected_real_kind(30, 291)
real (kind=dp), parameter :: pi = 3.141592653589793238D0

end module constants_module
