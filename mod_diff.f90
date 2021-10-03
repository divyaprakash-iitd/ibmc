module mod_diff
    ! Finite Difference Functions
    use iso_fortran_env, only: int32, real32, int64, real64
    implicit none

    private
    public :: diffx, diffy, diff2x, diff2y 
    
contains

    pure function diffx(M) 
        ! Centered Finite Difference in x
        real(real64), intent(in) :: M(:,:)
        real(real64), allocatable :: diffx(:,:)
        integer(int32) :: i, j

        allocate(diffx(size(M,1),size(M,2)))

        diffx = 0.0d0
        do concurrent (j = 1:size(M,2), i = 2:size(M,1)-1)
            diffx(i,j) = 0.5d0*(M(i+1,j) - M(i-1,j))
        end do
    end function diffx 
    
    pure function diffy(M) 
        ! Centered Finite Difference in y
        real(real64), intent(in) :: M(:,:)
        real(real64), allocatable :: diffy(:,:)
        integer(int32) :: i, j

        allocate(diffy(size(M,1),size(M,2)))

        diffy = 0.0d0
        do concurrent (j = 2:size(M,2)-1, i = 1:size(M,1))
            diffy(i,j) = 0.5d0*(M(i,j+1) - M(i,j-1))
        end do
    end function diffy
    
    pure function diff2x(M) 
        ! Centered Finite Difference in x
        real(real64), intent(in) :: M(:,:)
        real(real64), allocatable :: diff2x(:,:)
        integer(int32) :: i, j

        allocate(diff2x(size(M,1),size(M,2)))

        diff2x = 0.0d0
        do concurrent (j = 1:size(M,2), i = 2:size(M,1)-1)
            diff2x(i,j) = (M(i-1,j) - 2*M(i,j) + M(i+1,j))
        end do
    end function diff2x 
    
    pure function diff2y(M) 
        ! Centered Finite Difference in y
        real(real64), intent(in) :: M(:,:)
        real(real64), allocatable :: diff2y(:,:)
        integer(int32) :: i, j

        allocate(diff2y(size(M,1),size(M,2)))

        diff2y = 0.0d0
        do concurrent (j = 2:size(M,2)-1, i = 1:size(M,1))
            diff2y(i,j) = (M(i,j-1) - 2*M(i,j) + M(i,j+1))
        end do
    end function diff2y 
end module mod_diff