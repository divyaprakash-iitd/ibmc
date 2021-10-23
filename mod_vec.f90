module mod_vec
    use iso_fortran_env, only: int32, real64, real64
    implicit none
   
    private
    public :: vec

    type :: vec
        real(real64) :: x
        real(real64) :: y 
    end type 

    interface vec
        module procedure :: vec_constructor
    end interface vec

contains

    pure type(vec) function vec_constructor() result(self)
        self%x = 0.0d0
        self%y = 0.0d0
    end function
end module mod_vec
