program ibmc
    use iso_fortran_env, only: int32, real32, real64
    use mod_pressure, only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx, only: calculate_pressure_amgx
    use mod_mesh
    implicit none

    ! Execution time
    real(real64) :: start, finish

    ! Computational Domain
    real(real32)    :: Lx = 1.0
    real(real32)    :: Ly = 1.0
    integer(int32)  :: Nx = 20
    integer(int32)  :: Ny = 20

    ! Mesh Paramaters
    real(real32) :: dx, dxi
    real(real32) :: dy, dyi

    ! Simulation Paramaters
    real(real32) :: tsim    = 1
    real(real32) :: dt      = 0.001
    real(real32) :: t

    ! Physical Constants
    real(real32) :: nu = 1.0/100.0
    real(real32) :: rho = 1.0d0

    ! Index ranges
    integer(int32) :: imin, imax, jmin, jmax
    integer(int32) :: i,j

    ! Boundary values
    real(real32) :: utop, vtop, ubottom, vbottom, &
                    uleft, vleft, uright, vright

    ! Matrices to store fields
    real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), P(:,:), A(:,:,:)

    ! Temporary/Miscellaneous variable
    real(real64) :: ucenter, vcenter
    integer(int32) :: NN
    logical :: init_status

    ! Mesh
    type(mesh) :: M

    call cpu_time(start)

    ! Construct Mesh
    M = mesh('M',Lx,Ly,Nx,Ny)

    ! Define boundary conditions
    utop    = 1.0
    vtop    = 0.0
    ubottom = 0.0
    vbottom = 0.0
    uleft   = 0.0
    vleft   = 0.0
    uright  = 0.0
    vright  = 0.0

    ! Allocate matrices ofr u,v,us,vs,rhs
    allocate(u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    NN = 5
    allocate(A(1:Nx,1:Ny,NN))

    ! Initialize
    u   = 0.0d0
    v   = 0.0d0
    us  = 0.0d0
    vs  = 0.0d0
    R   = 0.0d0
    A   = 0.0d0

    ! Mesh values
    dx  = Lx/Nx
    dy  = Ly/Ny
    dxi = 1.0/dx
    dyi = 1.0/dy

    ! Generate Laplacian matrix
    call generate_laplacian_sparse(A,dxi,dyi)
    ! Start time loop
    t = 0.0d0
    init_status = .False.
    do while (t.lt.tsim)
        t = t+dt 

        ! Apply boundary conditions
        u(:,M%yu%lb) = 2*ubottom - u(:,M%yu%lb+1);
        u(:,M%yu%ub) = 2*utop    - u(:,M%yu%ub-1);
        v(M%xv%lb,:) = 2*vleft   - v(M%xv%lb+1,:);
        v(M%yv%ub,:) = 2*vright  - v(M%yv%ub-1,:);
    

        u(M%xu%lb,:)    = uleft;
        u(M%xu%ub,:)    = uright;
        v(:,M%yv%ub)    = vtop;
        v(:,M%yv%lb)    = vbottom;

        ! Perform predictor step
        ! us 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
            us(i,j) = u(i,j) + dt* &
                        ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                        + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                        - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                        - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi)
        end do

        ! vs 
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
            vs(i,j) = v(i,j) + dt* &
                        ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                        + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                        - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                        - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)
        end do

        ! Form the right hand side of the pressure poisson equation
        ! (Pressure cell)
        do concurrent (j = M%yp%lb:M%yp%ub, i = M%xp%lb:M%xp%ub)
            R(i,j) = -rho/dt* &
                    ( (us(i+1,j) - us(i,j))*dxi &
                    + (vs(i,j+1) - vs(i,j))*dyi)
        end do

        ! Solve for presssure
        call calculate_pressure_sparse(A,P,R)

        !call calculate_pressure_amgx(A,P,R,init_status)

        ! Perform the corrector steps
        ! u 
        ! (u-velocity cell)
        do concurrent (j = M%yu%lb+1:M%yu%ub-1, i = M%xu%lb+1:M%xu%ub-1)
            u(i,j) = us(i,j) - dt/rho * (p(i,j) - p(i-1,j)) * dxi
        end do
        
        ! v
        ! (v-velocity cell)
        do concurrent (j = M%yv%lb+1:M%yv%ub-1, i = M%xv%lb+1:M%xv%ub-1)
            v(i,j) = vs(i,j) - dt/rho * (p(i,j) - p(i,j-1)) * dyi
        end do

        print *, 'time = ', t
    end do

    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start 

    open(unit=55, file='u.txt', ACTION="write", STATUS="replace")
    do i=M%xu%lb,M%xu%ub
        write(55, '(*(F14.7))')( real(u(i,j)) ,j=M%yu%lb,M%yu%ub)
    end do

    open(unit=65, file='v.txt', ACTION="write", STATUS="replace")
    do i=M%xv%lb,M%xv%ub
        write(65, '(*(F14.7))')( real(v(i,j)) ,j=M%yv%lb,M%yv%ub)
    end do

    ! open(unit=75, file='p.txt', ACTION="write", STATUS="replace")
    ! do i=imin-1,imax+1
    !     write(75, '(*(F14.7))')( real(p(i,j)) ,j=jmin-1,jmax+1)
    ! end do

end program ibmc
