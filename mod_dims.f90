module mod_dims
!
!   Purpose:
!       To define a derived type containing information about the dimension (upper bound
!       and lower bound) of an array
!
    use iso_fortran_env, only: int32, real64, int32, real64
    implicit none
    
    type :: dims
        character(len=5) :: name
        integer(int32) :: lb, ub
    end type dims
contains
    
end module mod_dims