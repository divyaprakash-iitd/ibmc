module mod_field
    use iso_fortran_env, only: int32, real32, int64, real64
    use mod_dims
    use mod_diff, only: diffx_real => diffx, diffy_real => diffy, &
                        diff2x_real => diff2x, diff2y_real => diff2y
    implicit none

    ! private
    ! public :: Field, diffx, diffy, diff2x, diff2y, write_field

    type :: Field

        character(len=:), allocatable   :: name
        type(dims)                      :: x, y
        real(real64), allocatable       :: data(:,:)

        contains
            procedure, public, pass(self)   :: diffx
            procedure, public, pass(self)   :: apply_boundary_conditions
            procedure, public, pass(self)   :: print_field
            procedure, public, pass(self)   :: write_field
    end type Field

    interface Field
        module procedure :: field_constructor
    end interface Field

contains

    pure type(Field) function field_constructor(name,x,y) result(self)
            character(*), intent(in) :: name
            type(dims), intent(in) :: x, y

            self%name = name
            self%x = x
            self%y = y

            allocate(self%data(x%lb:x%ub,y%lb:y%ub))
            self%data = 0.0d0
    end function field_constructor

    subroutine apply_boundary_conditions(self, xmin, xmax, ymin, ymax)
        class(Field), intent(in out) :: self
        real(real32) :: xmin, xmax, ymin, ymax

        self%data(:,self%y%lb) = ymin
        self%data(:,self%y%ub) = ymax
        self%data(self%x%lb,:) = xmin
        self%data(self%x%ub,:) = xmax
    end subroutine apply_boundary_conditions

    subroutine print_field(self)
        class(Field), intent(in):: self

        print *, self%data
    end subroutine print_field

    function diffx(self)
        class(Field), intent(in) :: self
        real(real64) :: diffx(self%x%lb:self%x%ub,self%y%lb:self%y%ub)

        diffx = diffx_real(self%data)
    end function diffx
        
    function diffy(self)
        class(Field), intent(in) :: self
        real(real64) :: diffy(self%x%lb:self%x%ub,self%y%lb:self%y%ub)
    
        diffy = diffy_real(self%data)
    end function diffy

    function diff2x(self)
        class(Field), intent(in) :: self
        real(real64) :: diff2x(self%x%lb:self%x%ub,self%y%lb:self%y%ub)
    
        diff2x = diff2x_real(self%data)
    end function diff2x

    function diff2y(self)
        class(Field), intent(in) :: self
        real(real64) :: diff2y(self%x%lb:self%x%ub,self%y%lb:self%y%ub)
    
        diff2y = diff2y_real(self%data)
    end function diff2y

    subroutine set(self,data)
        class(Field), intent(in) :: self
        real(real64) :: data(:,:)
    end subroutine

    subroutine write_field(self, time)
        ! Writes the field in a text file
        class(Field), intent(in) :: self
        integer(int32), intent(in) :: time
        
        integer(int32) :: i,j
        character(len=100) :: filename, timestr

        write(timestr, '(i4.4)') time
        filename = self%name // '_' // trim(timestr) // '.txt'

        open(unit=75, file=filename, ACTION="write", STATUS="replace")
        do i = self%x%lb, self%x%ub
            write(75, '(*(F14.7))')( real(self%data(i,j)) ,j = self%y%lb, self%y%ub)
        end do

    end subroutine write_field

end module mod_field