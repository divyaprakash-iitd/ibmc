module mod_pressure
    use iso_fortran_env, only: int32, int32, real64, real64 
    use mod_mesh
    implicit none

    private
    public :: generate_laplacian_sparse, calculate_pressure_sparse

contains

subroutine calculate_pressure_sparse(A,PIN,R)

    real(real64),       intent(in)      :: R(:,:), A(:,:,:)
    real(real64),       intent(in out)  :: PIN(:,:)

    integer(int32)      :: imin, imax, jmin, jmax
    real(real64), allocatable :: P(:,:)
    integer(int32) :: Nx, Ny

    ! Indices 
    integer(int32) :: i,j,ITER

    ! Parameters for SOR
    real(real64),   parameter :: OMEGA      = 1.99
    real(real64),   parameter :: TOLERANCE  = 1E-8
    integer(int32), parameter :: MAXITER    = 1E5
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
    ! PIN = 0.0d0
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

subroutine generate_laplacian_sparse(A,dx,dy)

    real(real64),   intent(in)      :: dx, dy
    real(real64),   intent(in out)  :: A(:,:,:)
    
    
    real(real64)                    :: idx2, idy2 ! Squared inverse of dx and dy
    ! Indices
    integer(int32) :: i,j
    integer(int32), parameter :: BOTTOM = 1
    integer(int32), parameter :: LEFT   = 2
    integer(int32), parameter :: CENTER = 3
    integer(int32), parameter :: RIGHT  = 4
    integer(int32), parameter :: TOP    = 5
    integer(int32), parameter :: NN     = 5
    integer(int32) :: imin, imax, jmin, jmax

    ! Define limits
    imin = lbound(A,1)
    imax = ubound(A,1)
    jmin = lbound(A,2)
    jmax = ubound(A,2)

    idx2 = (1.0d0/dx)**2
    idy2 = (1.0d0/dy)**2

    ! Initialize
    A = 0.0d0
        
    ! Center Nodes
    do concurrent (j = jmin+1:jmax-1, i = imin+1:imax-1)
        ! d^2P/dx^2
        A(i,j,CENTER)   = -2*idx2      ! Center
        A(i,j,RIGHT)    = idx2         ! Right
        A(i,j,LEFT)     = idx2         ! Left
        
        ! d^2P/dy^2
        A(i,j,CENTER)   = A(i,j,CENTER) + (-2*idy2)       ! Center
        A(i,j,TOP)      = idy2                            ! Top
        A(i,j,BOTTOM)   = idy2                            ! Bottom
    end do
    
    ! Left Boundary Nodes (Excluding corners)
    do concurrent (j = jmin+1:jmax-1)
        do concurrent (i = imin:imin)
            
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
    do concurrent (j = jmin+1:jmax-1)
        do i = imax,imax
            
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
    do concurrent (j = jmax:jmax)
        do concurrent (i = imin+1:imax-1)
            
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
    do concurrent (j = jmin:jmin)
        do concurrent (i = imin+1:imax-1)
            
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
    i = imin
    j = jmax
    ! d^2P/dy^2
    A(i,j,TOP)      = 0                                        ! Top
    A(i,j,BOTTOM)   = idy2                                     ! Bottom
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,LEFT)     = 0                                        ! Right
    A(i,j,RIGHT)    = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Top-Right Boundary)
    i = imax
    j = jmax
    ! d^2P/dy^2
    A(i,j,TOP)      = 0                                        ! Top
    A(i,j,BOTTOM)   = idy2                                     ! Bottom
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,RIGHT)    = 0                                        ! Right
    A(i,j,LEFT)     = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Bottom-Left Boundary)
    i = imin
    j = jmin
    ! d^2P/dy^2
    A(i,j,BOTTOM)   = 0                                        ! Bottom
    A(i,j,TOP)      = idy2                                     ! Top
    A(i,j,CENTER)   = -idy2                                    ! Center
    ! d^2P/dx^2
    A(i,j,LEFT)     = 0                                        ! Right
    A(i,j,RIGHT)    = idx2                                     ! Left
    A(i,j,CENTER)   = A(i,j,CENTER) + (-idx2)                  ! Center
    
    ! Corner Node (Bottom-Right Boundary)
    i = imax
    j = jmin
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
    A(imax-1,jmin:jmax,RIGHT) = 0.0d0;
    A(imin+1,jmin:jmax,LEFT) = 0.0d0;
    ! A(imin,jmin,:) = 0.0d0
    ! A(imin,jmin,CENTER) = 1.0d0

end subroutine generate_laplacian_sparse

end module mod_pressure
