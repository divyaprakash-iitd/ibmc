program ibmc
    use iso_fortran_env,    only: int32, real32, real64
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    use mod_boundary
    use mod_time_stepping
    use mod_io
    use mod_ib
    use mod_ibm
    implicit none

    ! Code execution time
    real(real64)    :: start, finish
    ! Computational Domain
    real(real32)    :: Lx       = 1.0
    real(real32)    :: Ly       = 1.0
    ! Mesh Paramaters
    integer(int32)  :: Nx       = 30
    integer(int32)  :: Ny       = 30
    ! Simulation time Paramaters
    real(real32)    :: tsim     = 10
    real(real32)    :: dt       = 0.001
    real(real32)    :: t
    ! Physical Constants
    real(real32)    :: nu       = 1.0/100.0
    real(real32)    :: rho      = 1.0d0
    real(real64)    :: ks       = 1.0d0
    real(real64)    :: Rl       = 0.2d0
    real(real64)    :: kb       = 1.0d0
    real(real64)    :: theta    = 1.0472D0
    ! Boundary values
    real(real32)    :: utop, vtop, ubottom, vbottom, &
                       uleft, vleft, uright, vright
    ! Matrices to store fields
    real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), &
                                 P(:,:), A(:,:,:), Fx(:,:), Fy(:,:)
    ! Temporary/Miscellaneous variable
    integer(int32)  :: it, NN
    logical         :: init_status
    ! Mesh
    type(mesh)      :: M
    ! Immersed boundary
    type(ib)        :: ptcle

    !---------------------- Begin Calculations ------------------------------------!
    call cpu_time(start)

    ! Construct Mesh
    M = mesh('M',Lx,Ly,Nx,Ny)
    ! Allocate matrices ofr u,v,us,vs,rhs
    allocate(u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    NN = 5
    allocate(A(1:Nx,1:Ny,NN))
    ! Initialize
    u   = 0.0d0
    v   = 0.0d0
    us  = 0.0d0
    vs  = 0.0d0
    R   = 0.0d0
    A   = 0.0d0
    Fx  = 0.0d0
    Fy  = 0.0d0
    ! Define boundary conditions for velocity
    utop    = 1.0
    vtop    = 0.0
    ubottom = 0.0
    vbottom = 0.0
    uleft   = 0.0
    vleft   = 0.0
    uright  = 0.0
    vright  = 0.0

    ! Create the IB structure
    ptcle = ib('spring_array',3)
    call initialize_ib(ptcle)

    ! Generate Laplacian matrix
    call generate_laplacian_sparse(A,M%dx,M%dy)

    ! Start time loop
    t = 0.0d0
    it = 0
    init_status = .False. ! AmgX initialization status
    do while (t.lt.tsim)
        t = t+dt 
        it = it + 1

        ! Apply velocity boundary conditions
        call apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

        ! Calculate forces in the immersed boundary structure
        call calculate_spring_force(ptcle,Ks,Rl)
        call calculate_torsional_spring_force(ptcle,kb,theta)

        ! Spread force from the immersed boundary
        call spread_force(M,ptcle,Fx,Fy)

        ! Calculate intermediate/predicted velocity
        ! call predictor(M,u,v,us,vs,nu,dt) 
        call euler(M,u,v,us,vs,nu,dt,Fx,Fy)
        ! call euler(M,u,v,us,vs,nu,dt)
        ! call RK2(M,u,v,us,vs,nu,dt,Fx,Fy)
        ! call RK4(M,u,v,us,vs,nu,dt,Fx,Fy)
        
        ! Form the right hand side of the pressure poisson equation
        call calculate_rhs(M,us,vs,R,rho,dt)

        ! Solve for presssure
        ! call calculate_pressure_sparse(A,P,R)
        call calculate_pressure_amgx(A,P,R,init_status)

        ! Perform the corrector steps to obtain the velocity
        call corrector(M,u,v,us,vs,p,rho,dt)

        ! Interpolate the Eulerian grid velocity to the Lagrangian structure
        call interpolate_velocity(M,ptcle,u,v)

        ! Update the Immersed Boundary
        call update_ib(ptcle,dt) 

        print *, 'time = ', t

        ! Write files every 10th timestep
        if (mod(it,10).eq.0) then 
            call write_field(u,'u',it) 
            call write_location(ptcle,it)
        end if
    end do

    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

    call write_mesh(M,'u')
    call write_mesh(M,'v')
    call write_mesh(M,'p')

    print *, ptcle%boundary(1)%x
    print *, ptcle%boundary(1)%y
end program ibmc
