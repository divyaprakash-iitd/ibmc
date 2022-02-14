module mod_cell
!
!   Purpose:
!       The define the derived type, particleptr and cell.
!       particleptr: It is used to store the pointers of the particle data type which lie in a cell
!                    along with an id for the cilia which differentiates between the types of the cilia
!                    arrays (Particle or Cilia) and the numbers the cilia within each array.
!       cell:        It is used to create square cells of the given dimensions and collect the particles which
!                    lie witihin the cell in it's particleptr array.
!                      
    use iso_fortran_env, only: int32, real64, real64
    use mod_vec
    use mod_particle
    implicit none
   
    private
    public :: cell
    
    ! Define a type to store an array of pointers in Fortran
    type particleptr
        type(particle), pointer :: pptr         ! A pointer to the 'particle derived data type'
        integer(int32)          :: ciliaId(2)   ! Cilia id to differentiate 
                                                ! (1) The type of cilia (particle or cilia)
                                                ! (2) The cilium within the array of cilia or particle
    end type particleptr

    type :: cell
        type(vec)         :: loc           ! Location of the Bottom-Left point
        type(vec)         :: L             ! Dimension of the square cell
        type(particleptr) :: Nlist(100)    ! Pointers of all the particles in a cell
        integer(int32)    :: NN            ! No. of particles in the cell
        integer(int32)    :: cellId        ! The tag to identify a cell 
    end type 

    interface cell
        module procedure :: cell_constructor
    end interface cell

contains

    pure type(cell) function cell_constructor() result(self)
    !
    !   Purpose: 
    !       To initialize the attributes of the cell and particleptr data type 
    !
        integer(int32) :: i 

        self%loc%x = 0.0d0
        self%loc%y = 0.0d0
        self%L%x = 0.0d0
        self%L%y = 0.0d0
        self%NN = 0
        self%cellId = 0
        do i = 1,size(self%Nlist)
            nullify(self%Nlist(i)%pptr)
            self%Nlist(i)%ciliaId = 0
        end do
    end function
end module mod_cell