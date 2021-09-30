module mod_helper
    use iso_fortran_env, only: int32, int64, real32, real64
    implicit none
    save

    private
    public :: generate_crs_matrix

    contains

    subroutine generate_crs_matrix()
        use iso_c_binding, only: c_int, c_double, c_loc
        use mod_global, only: A,Nx,Ny,col_ind,row_ptr,val
        
        integer(int32) :: i,j,flag,L,R,T,B,row,col
        real(real64) :: zero
        real(real64) :: AA(1:Nx,1:Ny,1:5)

        ! Neighbours of a grid point
        integer(int32), parameter :: BOTTOM = 1
        integer(int32), parameter :: LEFT   = 2
        integer(int32), parameter :: CENTER = 3
        integer(int32), parameter :: RIGHT  = 4
        integer(int32), parameter :: TOP    = 5
        integer(int32), parameter :: NN     = 5

        ! Initialize empty arrays
        val     = [real ::]
        col_ind = [integer ::]
        row_ptr = [integer ::]

        AA = A
        zero = 1e-14

        do j = 1,Ny
            do i = 1,Nx
                flag = 0
                L = i-1
                R = i+1
                T = j+1
                B = j-1 

                row = id(i,j)
                ! Bottom
                if (abs(AA(i,j,BOTTOM)).gt.zero) then
                    col = id(i,B)

                    val     = [val, AA(i,j,BOTTOM)]
                    col_ind = [col_ind, col]
                    row_ptr = [row_ptr, ubound(col_ind)]
                    flag = 1 
                end if

                ! Left
                if (abs(AA(i,j,LEFT)).gt.zero) then
                    col = id(L,j)

                    val     = [val, AA(i,j,LEFT)]
                    col_ind = [col_ind, col]
                    if (flag == 0) then
                        row_ptr = [row_ptr, ubound(col_ind)]
                        flag = 1
                    end if
                end if

                !Center
                if (abs(AA(i,j,CENTER)).gt.zero) then
                    col = id(i,j)
                    
                    val     = [val, AA(i,j,CENTER)]
                    col_ind = [col_ind, col]
                    if (flag == 0) then
                        row_ptr = [row_ptr, ubound(col_ind)]
                        flag = 1
                    end if
                end if

                ! Right
                if (abs(AA(i,j,RIGHT)).gt.zero) then
                    col = id(R,j)
                    
                    val     = [val, AA(i,j,RIGHT)]
                    col_ind = [col_ind, col]
                end if

                ! Top
                if (abs(AA(i,j,TOP)).gt.zero) then
                    col = id(i,T)
                    
                    val     = [val, AA(i,j,TOP)]
                    col_ind = [col_ind, col]
                end if

            end do
        end do
        row_ptr = row_ptr - 1
        col_ind = col_ind - 1        

        row_ptr = [row_ptr, size(val)+1]

       !open(unit=85, file='L_CRS.txt', ACTION="write", STATUS="replace")
       !do i=1,size(lap,1)
       !    write(85, '(*(F14.7))')( real(lap(i,j)) ,j=1,size(lap,2))
       !end do
       
       open(unit=3, file='val.txt', action="write", status="replace")
       do i=1,size(val)
           write(3, '(*(f14.7))') real(val(i))
       end do
       close(unit=3)

       open(unit=3, file='col_ind.txt', action="write", status="replace")
       do i=1,size(col_ind)
           write(3, '(*(f14.7))') real(col_ind(i))
       end do
       close(unit=3)

       open(unit=3, file='row_ptr.txt', action="write", status="replace")
       do i=1,size(row_ptr)
           write(3, '(*(f14.7))') real(row_ptr(i))
       end do
       close(unit=3)
        
        contains 

        function id(ix,jy)  
            integer(int32),intent(in) :: ix, jy
            integer(int32) :: id

            ! Elements are mapped column by column along a single row
            id = (jy-1)*Nx + ix
        end function id
    end subroutine generate_crs_matrix

end module mod_helper
