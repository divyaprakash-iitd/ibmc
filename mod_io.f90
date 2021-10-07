module mod_io
    use iso_fortran_env, only: int32, real32, real64
    use mod_vec
    implicit none
    
    private
    public :: write_field
contains

    subroutine write_field(F,fieldname,timestep)
        real(real64), intent(in) :: F(:,:)
        character(len=1), intent(in) :: fieldname
        integer(int32), intent(in) :: timestep

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        character(len=8) :: itnumber
        integer(int32) :: i,j

        write(itnumber,"(I8.8)") timestep
        filename = fieldname // '_' // itnumber // '.txt'
        
        open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
        do j = 1,size(F,2)
            write(fileunit, '(*(F14.7))')( real(F(i,j)) , i = 1,size(F,1))
        end do
        close(fileunit)
    end subroutine write_field

    subroutine write_mesh(phi_mesh,c)
        class(vec), intent(in) :: phi_mesh(:,:)
        character(len=1), intent(in) :: c

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: i, j

        filename = c // '_' // 'x' // 'mesh' //'.txt'
        open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
        do j = 1,size(phi_mesh,2)
            write(fileunit, '(*(F14.7))')( real(phi_mesh(i,j)%x) , i = 1,size(phi_mesh,1))
        end do
        close(fileunit)

    end subroutine

end module mod_io