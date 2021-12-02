module mod_amgx
    use iso_fortran_env, only: int32, int32, real64, real64
    use iso_c_binding, only: c_int, c_double, c_loc
    use nvtx
    implicit none
    save

    private
    public :: calculate_pressure_amgx

    ! CRS data (Common across all the subroutines in this module)
    real(c_double), allocatable, target :: val(:)
    integer(c_int), allocatable, target :: col_ind(:)
    integer(c_int), allocatable, target :: row_ptr(:)
    real(c_double), allocatable, target :: rhsv(:)
    real(c_double), allocatable, target :: sol(:)
    integer(c_int),dimension(4), target :: crs_data
    integer(c_int) :: NNZ
    integer(int32) :: init_am, solve_am ! Return Variables
    
    ! Domain data
    integer(int32) :: Nx, Ny

    contains

    subroutine calculate_pressure_amgx(A,x,b,init_status)

        real(real64), intent(in)        :: A(:,:,:)    ! Coefficient matrix
        real(real64), intent(in out)    :: x(:,:)        ! Pressure matrix
        real(real64), intent(in)        :: b(:,:)        ! Divergence matrix 
        logical, intent(in out) :: init_status
        
        if (init_status .eqv..False.) then

            call nvtxStartRange('Generate CRS matrix')
            call generate_crs_matrix(A)
            call nvtxEndRange
            NNZ = size(col_ind)

            ! Get the CRS data to pass to AmgX
            crs_data = [Nx*Ny, NNZ, 1, 1]
            allocate(sol(Nx*Ny))
            sol = 0.0d0

            ! Reshape divergence matrix
            rhsv = reshape(b,[Nx*Ny])

            call nvtxStartRange('AmgX: Initialize and solve')
            call init_amgx()
            call nvtxEndRange
            init_status = .True.
        else
            rhsv = reshape(b,[Nx*Ny])
            sol = reshape(x,[Nx*Ny])
            call nvtxStartRange('AmgX: Solve')
            call solve_amgx()
            call nvtxEndRange
        end if

        ! Reshape the solution to a matrix
        x = reshape(sol,[Nx,Ny])
    end subroutine calculate_pressure_amgx

    subroutine generate_crs_matrix(A)
       
        real(real64), intent(in) :: A(:,:,:)

        integer(int32) :: i,j,flag,L,R,T,B,row,col
        real(real64) :: zero

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

        zero = 1e-14

        Nx = size(A,1)
        Ny = size(A,2)

        do j = 1,Ny
            do i = 1,Nx
                flag = 0
                L = i-1
                R = i+1
                T = j+1
                B = j-1 

                row = id(i,j)
                ! Bottom
                if (abs(A(i,j,BOTTOM)).gt.zero) then
                    col = id(i,B)

                    val     = [val, A(i,j,BOTTOM)]
                    col_ind = [col_ind, col]
                    row_ptr = [row_ptr, ubound(col_ind)]
                    flag = 1 
                end if

                ! Left
                if (abs(A(i,j,LEFT)).gt.zero) then
                    col = id(L,j)

                    val     = [val, A(i,j,LEFT)]
                    col_ind = [col_ind, col]
                    if (flag == 0) then
                        row_ptr = [row_ptr, ubound(col_ind)]
                        flag = 1
                    end if
                end if

                !Center
                if (abs(A(i,j,CENTER)).gt.zero) then
                    col = id(i,j)
                    
                    val     = [val, A(i,j,CENTER)]
                    col_ind = [col_ind, col]
                    if (flag == 0) then
                        row_ptr = [row_ptr, ubound(col_ind)]
                        flag = 1
                    end if
                end if

                ! Right
                if (abs(A(i,j,RIGHT)).gt.zero) then
                    col = id(R,j)
                    
                    val     = [val, A(i,j,RIGHT)]
                    col_ind = [col_ind, col]
                end if

                ! Top
                if (abs(A(i,j,TOP)).gt.zero) then
                    col = id(i,T)
                    
                    val     = [val, A(i,j,TOP)]
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
       
    !    open(unit=3, file='val.txt', action="write", status="replace")
    !    do i=1,size(val)
    !        write(3, '(*(f14.7))') real(val(i))
    !    end do
    !    close(unit=3)

    !    open(unit=3, file='col_ind.txt', action="write", status="replace")
    !    do i=1,size(col_ind)
    !        write(3, '(*(f14.7))') real(col_ind(i))
    !    end do
    !    close(unit=3)

    !    open(unit=3, file='row_ptr.txt', action="write", status="replace")
    !    do i=1,size(row_ptr)
    !        write(3, '(*(f14.7))') real(row_ptr(i))
    !    end do
    !    close(unit=3)
        
        contains 

        function id(ix,jy)  
            integer(int32),intent(in) :: ix, jy
            integer(int32) :: id

            ! Elements are mapped column by column along a single row
            id = (jy-1)*Nx + ix
        end function id
    end subroutine generate_crs_matrix

    subroutine init_amgx()
        use ftn_c

        init_am = initialize_amgx(c_loc(crs_data), c_loc(val), c_loc(col_ind), c_loc(row_ptr),c_loc(rhsv),c_loc(sol))
        write(*,*) 'AmgX Initialization done!'
 
    end subroutine init_amgx
   
    subroutine solve_amgx()
        use ftn_c

        solve_am = solveamg(c_loc(rhsv),c_loc(sol))
        
    end subroutine solve_amgx

end module mod_amgx
