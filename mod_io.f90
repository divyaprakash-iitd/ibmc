module mod_io
    use iso_fortran_env, only: int32, real32, real64
    use mod_mesh
    implicit none
    
    private
    public :: write_field, write_mesh
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

    subroutine write_mesh(M,c)
        class(mesh), intent(in) :: M
        character(len=1), intent(in) :: c

        integer(int32) :: fileunit = 8
        character(len=:), allocatable :: filename
        integer(int32) :: i, j

        if (c.eq.'u') then 
            filename = c // '_x_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yu%lb,M%yu%ub 
                write(fileunit, '(*(F14.7))')(M%u_mesh(i,j)%x , i = M%xu%lb,M%xu%ub)
            end do
            close(fileunit)

            filename = c // '_y_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yu%lb,M%yu%ub
                write(fileunit, '(*(F14.7))')(M%u_mesh(i,j)%y , i = M%xu%lb,M%xu%ub)
            end do
            close(fileunit)
        end if

        if (c.eq.'v') then 
            filename = c // '_x_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yv%lb,M%yv%ub 
                write(fileunit, '(*(F14.7))')(M%v_mesh(i,j)%x , i = M%xv%lb,M%xv%ub)
            end do
            close(fileunit)

            filename = c // '_y_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yv%lb,M%yv%ub
                write(fileunit, '(*(F14.7))')(M%v_mesh(i,j)%y , i = M%xv%lb,M%xv%ub)
            end do
            close(fileunit)
        end if

        if (c.eq.'p') then
            filename = c // '_x_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yp%lb,M%yp%ub 
                write(fileunit, '(*(F14.7))')(M%p_mesh(i,j)%x , i = M%xp%lb,M%xp%ub)
            end do
            close(fileunit)

            filename = c // '_y_' // 'mesh' //'.txt'
            open(unit=fileunit, file=filename, ACTION="write", STATUS="replace")
            do j = M%yp%lb,M%yp%ub
                write(fileunit, '(*(F14.7))')(M%p_mesh(i,j)%y , i = M%xp%lb,M%xp%ub)
            end do
            close(fileunit)
        end if

    end subroutine write_mesh

end module mod_io