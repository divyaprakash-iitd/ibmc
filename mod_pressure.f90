module mod_pressure
    use iso_fortran_env, only: int32, int64, real32, real64 
    implicit none

    private
    public :: generate_laplacian_sparse, calculate_pressure_sparse

contains

subroutine generate_laplacian_sparse(idx,idy,A,imin,imax,jmin,jmax)


    integer(int32), intent(in)      :: imin,imax,jmin,jmax
    real(real32),   intent(in)      :: idx, idy
    real(real64),   intent(in out)  :: A(imin:imax,jmin:jmax,5)
    
    
    real(real64)                    :: idx2, idy2 ! Squared inverse of dx and dy
    ! Indices
    integer(int32) :: i,j,iminL,jminL,imaxL,jmaxL
    integer(int32), parameter :: BOTTOM = 1
    integer(int32), parameter :: LEFT   = 2
    integer(int32), parameter :: CENTER = 3
    integer(int32), parameter :: RIGHT  = 4
    integer(int32), parameter :: TOP    = 5
    integer(int32), parameter :: NN     = 5

    idx2 = idx**2
    idy2 = idy**2

    iminL = imin
    imaxL = imax
    jminL = jmin
    jmaxL = jmax

    ! Initialize
    A = 0.0d0
        
    ! Center Nodes
    do concurrent (j = jminL+1:jmaxL-1)
        do concurrent (i = iminL+1:imaxL-1)
            
            ! d^2P/dx^2
            A(i,j,CENTER)   = -2*idx2      ! Center
            A(i,j,RIGHT)    = idx2         ! Right
            A(i,j,LEFT)     = idx2         ! Left
            
            ! d^2P/dy^2
            A(i,j,CENTER)   = A(i,j,CENTER) + (-2*idy2)       ! Center
            A(i,j,TOP)      = idy2                            ! Top
            A(i,j,BOTTOM)   = idy2                            ! Bottom
        end do
    end do
    
    ! Left Boundary Nodes (Excluding corners)
    do concurrent (j = jminL+1:jmaxL-1)
        do concurrent (i = iminL:iminL)
            
            ! d^2P/dx^2
            A(i,j,LEFT)    = 0                                ! Left
            A(i,j,RIGHT)   = idx2                             ! Right
            A(i,j,CENTER)  = -idx2                            ! Center
            
            ! d^2P/dy^2 
            A(i,j,CENTER)  = A(i,j,CENTER) + (-2*idy2)        ! Center
            A(i,j,BOTTOM)  = idy2                             ! Bottom
            A(i,j,TOP)     = idy2                             ! Top                    
        end do
    end do
    
    ! Right Boundary Nodes (Excluding corners)
    do concurrent (j = jminL+1:jmaxL-1)
        do i = imaxL,imaxL
            
            ! d^2P/dx^2
            A(i,j,RIGHT)    = 0                               ! Right
            A(i,j,LEFT)     = idx2                            ! Left
            A(i,j,CENTER)   = -idx2                           ! Center
            
            ! d^2P/dy^2
            A(i,j,CENTER)   = A(i,j,CENTER) + (-2*idy2)       ! Center
            A(i,j,TOP)      = idy2                            ! Top
            A(i,j,BOTTOM)   = idy2                            ! Bottom
        end do
    end do
    
    ! Top Boundary Nodes (Excluding corners)
    do concurrent (j = jmaxL:jmaxL)
        do concurrent (i = iminL+1:imaxL-1)
            
            ! d^2P/dy^2
            A(i,j,TOP)      = 0                               ! Top
            A(i,j,BOTTOM)   = idy2                            ! Bottom
            A(i,j,CENTER)   = -idy2                           ! Center
            
            ! d^2P/dx^2
            A(i,j,CENTER)   = A(i,j,CENTER) + (-2*idx2)       ! Center
            A(i,j,RIGHT)    = idx2                            ! Right
            A(i,j,LEFT)     = idx2                            ! Left
        end do
    end do
        
    ! Bottom Boundary Nodes (Excluding corners)
    do concurrent (j = jminL:jminL)
        do concurrent (i = iminL+1:imaxL-1)
            
            ! d^2P/dy^2
            A(i,j,BOTTOM)   = 0                                ! Bottom
            A(i,j,TOP)      = idy2                             ! Top
            A(i,j,CENTER)   = -idy2                            ! Center
            
            ! d^2P/dx^2
            A(i,j,CENTER)   = A(i,j,CENTER) + (-2*idx2)        ! Center
            A(i,j,RIGHT)    = idx2                             ! Right
            A(i,j,LEFT)     = idx2                             ! Left                    
        end do
    end do
    
    ! Corner Node (Top-Left Boundary)
    i = iminL
    j = jmaxL
    ! d^2P/dy^2
    A(i,j,TOP)      = 0                                        ! Top
    A(i,j,BOTTOM)   = idy2                                     ! Bottom
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,LEFT)     = 0                                        ! Right
    A(i,j,RIGHT)    = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Top-Right Boundary)
    i = imaxL
    j = jmaxL
    ! d^2P/dy^2
    A(i,j,TOP)      = 0                                        ! Top
    A(i,j,BOTTOM)   = idy2                                     ! Bottom
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,RIGHT)    = 0                                        ! Right
    A(i,j,LEFT)     = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Bottom-Left Boundary)
    i = iminL
    j = jminL
    ! d^2P/dy^2
    A(i,j,BOTTOM)   = 0                                        ! Bottom
    A(i,j,TOP)      = idy2                                     ! Top
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,LEFT)     = 0                                        ! Right
    A(i,j,RIGHT)    = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Bottom-Right Boundary)
    i = imaxL
    j = jminL
    ! d^2P/dy^2
    A(i,j,BOTTOM)   = 0                                        ! Bottom
    A(i,j,TOP)      = idy2                                     ! Top
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,RIGHT)    = 0                                        ! Right
    A(i,j,LEFT)     = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center

    ! Set the pressure at the corner cell
    A = -A
    A(imin,jmin,:) = 0.0d0
    A(imin,jmin,CENTER) = 1.0d0

end subroutine generate_laplacian_sparse


subroutine calculate_pressure_sparse(A,PIN,R)

    real(real64),       intent(in)      :: R(:,:), A(:,:,:)
    real(real64),       intent(in out)  :: PIN(:,:)

    integer(int32)      :: imin, imax, jmin, jmax
    real(real64), allocatable :: P(:,:)
    integer(int32) :: Nx, Ny

    ! Indices 
    integer(int64) :: i,j,ITER

    ! Parameters for SOR
    real(real32),   parameter :: OMEGA      = 1.99
    real(real32),   parameter :: TOLERANCE  = 1E-10
    integer(int64), parameter :: MAXITER    = 1E5
    integer(int32), parameter :: BOTTOM     = 1
    integer(int32), parameter :: LEFT       = 2
    integer(int32), parameter :: CENTER     = 3
    integer(int32), parameter :: RIGHT      = 4
    integer(int32), parameter :: TOP        = 5
    integer(int32), parameter :: NN         = 5
    
    real(real64) :: ERROR
    real(real64) :: PTEMP

    real(real64), allocatable :: RSDL(:,:)     ! Residual Matrix
    real(real64), allocatable :: RSDLV(:)    ! Residual Vector
    
    imin = lbound(A,1)
    imax = ubound(A,1)
    jmin = lbound(A,2)
    jmax = ubound(A,2)
   
    ! Domain Size
    Nx = size(PIN,1)
    Ny = size(PIN,2)

    ! Allocate variables
    allocate(P(imin-1:imax+1,jmin-1:jmax+1))
    allocate(RSDL(Nx,Ny))

    ! Pressure array with zero padding on the boundaries
    P = 0.0d0
    PIN = 0.0d0
    P(imin:imax,jmin:jmax) = PIN

    allocate(RSDLV(Nx*Ny))        
    RSDL = 0.0d0
    ERROR = 1.0d0

    ITER = 0
    do while ((ITER.lt.MAXITER).and.(ERROR.gt.TOLERANCE))
        do j = jmin,jmax
            do i = imin,imax
                PTEMP = R(i,j) - &
                        (A(i,j,BOTTOM)   *   P(i,j-1)    +   &
                        A(i,j,LEFT)     *   P(i-1,j)    +   &
                        A(i,j,RIGHT)    *   P(i+1,j)    +   &
                        A(i,j,TOP)      *   P(i,j+1))

                P(i,j) = (1-OMEGA) * P(i,j) + OMEGA/A(i,j,CENTER) * PTEMP
            end do
        end do

        ! Calculate Residual
        do concurrent (j = jmin:jmax)
            do concurrent (i = imin:imax)
                RSDL(i,j) = A(i,j,BOTTOM)   *   P(i,j-1)    +   &
                            A(i,j,LEFT)     *   P(i-1,j)    +   &
                            A(i,j,CENTER)   *   P(i,j)      +   &
                            A(i,j,RIGHT)    *   P(i+1,j)    +   &
                            A(i,j,TOP)      *   P(i,j+1)    -   &
                            R(i,j)
            end do
        end do

        RSDLV = reshape(RSDL, [Nx*Ny])
        ERROR = norm2(abs(RSDLV))/Nx/Ny

        ITER = ITER + 1
    end do
    PIN = P(imin:imax,jmin:jmax)
    print *, 'Error = ', ERROR
    print *, 'Iteration = ', ITER
end subroutine calculate_pressure_sparse

end module mod_pressure
