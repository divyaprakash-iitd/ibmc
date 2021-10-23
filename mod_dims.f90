module mod_dims
    use iso_fortran_env, only: int32, real64, int32, real64
    implicit none
    
    type :: dims
        character(len=5) :: name
        integer(int32) :: lb, ub
    end type dims
contains
    
end module mod_dims