program ibmc
    use iso_fortran_env,    only: int32, int64, real32, real64
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    use mod_boundary
    use mod_time_stepping
    use mod_io
    use mod_ib
    use mod_ibm
    use mod_vec
    use mod_cilia
    use mod_cilia_array
    implicit none

    ! Parameters
    integer(int64), parameter :: PI = 3.141592653589793
    ! Code execution time
    real(real64)    :: start, finish
    ! Computational Domain
    real(real32)    :: Lx       = 1.0
    real(real32)    :: Ly       = 1.0
    ! Mesh Paramaters
    integer(int32)  :: Nx       = 30
    integer(int32)  :: Ny       = 30
    ! Simulation time Paramaters
    real(real32)    :: tsim     = 20
    real(real32)    :: dt       = 0.0001
    real(real32)    :: t
    ! Physical Constants
    real(real32)    :: nu       = 1.0/100
    real(real32)    :: rho      = 1.0d0
    real(real64)    :: ks       = 0.1d0
    real(real64)    :: kb       = 1.5d0
    real(real64)    :: theta    = 3.1416
    real(real64)    :: Rl
    real(real32)    :: FP(2)
    real(real64)    :: SP(3)
    ! Boundary values
    real(real32)    :: utop, vtop, ubottom, vbottom, &
                       uleft, vleft, uright, vright
    real(real32)    :: BC(8)
    ! Matrices to store fields
    real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), &
                                 P(:,:), A(:,:,:), Fx(:,:), Fy(:,:)
    ! Temporary/Miscellaneous variable
    integer(int32)  :: it, NN, il, ip
    logical         :: init_status
    real(real64)    :: tp = 2.0d0 ! Time period
    real(real64)    :: Ftip = 0.01 ! Tip force
    ! Mesh
    type(mesh)      :: M
    ! Immersed boundary
    type(ib)        :: ptcle
    character(1)    :: btype        ! Boundary type (Open or Closed)
    type(vec)       :: origin       ! Origin of the immersed boundary
    real(real64)    :: ibL          ! Length of the immersed boundary
    integer(int32)  :: np           ! Number of particles in the immersed boundary
    integer(int32)  :: nl           ! Number of layers in the cilia
    real(real64)    :: Wbl          ! Width of cilia
    real(real64)    :: dc           ! Distance between cilium in an array
    integer(int32)  :: nc       ! Total number of cilia
    ! Cilia
    type(cilia)     :: Cil
    ! Cilia array
    type(cilia_array) :: CA
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

    BC = [utop,vtop,ubottom,vbottom,uleft,vleft,uright,vright]
    FP = [nu,rho]
    SP = [ks,Rl,Ftip]
    
    ! Generate Laplacian matrix
    call generate_laplacian_sparse(A,M%dx,M%dy)

    init_status = .False. ! AmgX initialization status
    call write_mesh(M,'u')
    call write_mesh(M,'v')
    call write_mesh(M,'p')

    ! Create cilia
    nl = 2
    np = 5
    ibl = 0.3
    wbl = 0.05
    Rl = ibL/(np-1) 
    origin = vec(0.4,0.05)

    nc = 3 ! Number of cilia
    dc = 16*M%dx ! Distance between the cilia structures

    ! CA = cilia_array(nc,nl,np)

    ! call create_cilia_array(CA,ibl,wbl,dc,origin)

    call time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,A,P,R,tsim,dt)

    call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start

end program ibmc
