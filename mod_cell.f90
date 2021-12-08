module mod_cell
    use iso_fortran_env, only: int32, real64, real64
    use mod_vec
    use mod_particle
    implicit none
   
    private
    public :: cell

    type particleptr
        type(particle), pointer :: pptr
        integer(int32) :: ciliaId(2)
    end type particleptr

    type :: cell
        type(vec)                           :: loc                  ! Location of the Bottom-Left point
        type(vec)                           :: L                    ! Dimension of the square cell
        type(particleptr)                   :: Nlist(50)    ! Pointers of all the particles in a cell
        integer(int32)                      :: NN ! No. of particles in the cell 
    end type 

    interface cell
        module procedure :: cell_constructor
    end interface cell

contains

    pure type(cell) function cell_constructor() result(self)

        integer(int32) :: i 

        self%loc%x = 0.0d0
        self%loc%y = 0.0d0
        self%L%x = 0.0d0
        self%L%y = 0.0d0
        self%NN = 0
        do i = 1,size(self%Nlist)
            nullify(self%Nlist(i)%pptr)
            self%Nlist(i)%ciliaId = 0
        end do
    end function
end module mod_cell