program ibmc
    use iso_fortran_env,    only: int32, real64
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    use mod_time_stepping
    use mod_boundary
    use mod_io
    use mod_ib
    use mod_ibm
    use mod_vec
    use mod_cilia
    use mod_cilia_array
    use mod_closed_cilia
    implicit none

    ! Namelists for input
    namelist /time/ dt, it_save, tsim
    namelist /grid/ Nx, Ny, Lx, Ly
    namelist /flow/ nu, rho, uleft
    namelist /ciliaprop/ nc, np, ko, kd
    namelist /particleprop/ nparticles, radius, npparticles, kop, kod
  
    ! Parameters
    integer(int32), parameter :: PI = 3.141592653589793
    ! Code execution time
    real(real64)    :: start, finish
    ! Computational Domain
    real(real64)    :: Lx       = 5.0d0
    real(real64)    :: Ly       = 1.0d0
    ! Mesh Paramaters
    integer(int32)  :: Nx       = 500
    integer(int32)  :: Ny       = 100
    ! Simulation time Paramaters
    real(real64)    :: tsim     = 1.000d0
    real(real64)    :: dt       = 0.001d0
    real(real64)    :: t
    integer(int32)  :: it_save  = 100
    ! Physical Constants
    real(real64)    :: nu       = 1.0d0/50.0d0
    real(real64)    :: rho      = 1.0d0
    real(real64)    :: ko       = 1.0d0
    real(real64)    :: kd       = 0.5d0
    real(real64)    :: kb       = 10.0d0
    real(real64)    :: theta    = 3.1416d0
    real(real64)    :: Rl
    real(real64)    :: FP(2)
    real(real64)    :: SP(4)
    ! Boundary values
    real(real64)    :: utop, vtop, ubottom, vbottom, &
                       uleft, vleft, uright, vright
    real(real64)    :: BC(8)
    ! Matrices to store fields
    real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), &
                                 P(:,:), A(:,:,:), Fx(:,:), Fy(:,:)
    ! Temporary/Miscellaneous variable
    integer(int32)  :: it, NN, il, ip, err
    logical         :: init_status
    real(real64)    :: tp = 2.0d0 ! Time period
    real(real64)    :: Ftip = -1.0d0 ! Tip force
    ! Mesh
    type(mesh)      :: M
    ! Immersed boundary
    type(ib)        :: ptcle
    character(1)    :: btype        ! Boundary type (Open or Closed)
    type(vec)       :: origin       ! Origin of the immersed boundary
    real(real64)    :: ibL          ! Length of the immersed boundary
    integer(int32)  :: np           ! Number of particles in the immersed boundary
    integer(int32)  :: nl           ! Number of layers in the cilia
    real(real64)    :: dp           ! Spacing between two particles in a cilia
    real(real64)    :: Ll           ! Length of a layer of cilia
    real(real64)    :: Wbl          ! Width of cilia
    real(real64)    :: dc           ! Distance between cilium in an array
    integer(int32)  :: nc       ! Total number of cilia
    ! Cilia
    type(cilia)     :: Cil
    ! Cilia array
    type(cilia_array) :: CA
    type(cilia_array) :: CAmid

    ! Closed cilia parameters
    type(cilia_array) :: CAP
    real(real64) :: radius
    type(vec) :: originP
    real(real64) :: kop = 1.0d0
    real(real64) :: kod = 0.50d0
    integer(int32) :: nparticles
    integer(int32) :: npparticles
    !---------------------- Begin Calculations ------------------------------------!
    call cpu_time(start)

    ! Read input data from file
    open(1004,file="input_params.txt",form='formatted')
    READ(unit=1004,nml=grid,iostat=err)
    close(1004)

    ! Construct and write Mesh data
    M = mesh('M',Lx,Ly,Nx,Ny)
    call write_mesh(M,'u')
    call write_mesh(M,'v')
    call write_mesh(M,'p')

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
    utop    = 0.0d0
    vtop    = 0.0d0
    ubottom = 0.0d0
    vbottom = 0.0d0
    uleft   = 0.1d0
    vleft   = 0.0d0
    uright  = 0.0d0
    vright  = 0.0d0

    
    ! Generate Laplacian matrix
    call generate_laplacian_sparse(A,M%dx,M%dy)

    ! AmgX initialization status
    init_status = .False.

    ! Create cilia
    nl      = 2                     ! No. of Layers/Cilia (More or less than 2 is not supported for now!)
    Rl      = 4*M%dx                ! Resting Length of Spring
    dp      = Rl                    ! Spacing between two particles
    np      = 6                     ! No. of Particles/Layer
    wbl     = Rl                    ! Width/Distance between two Layers
    dc      = 3*Rl                    ! Distance between two Cilia
    nc      = 15                    ! Number of cilia
    origin  = vec(Lx/4,0.1d0)      ! Location of the first Cilium (Bottom-Left Particle)
    radius = 0.02*Lx
    originP = vec(Lx/9,2*Ly/3)
    nparticles = 1
    npparticles = 8
    ! originP = vec(Lx/3,2.25*Ly/3)

    ! Read input data from file
    open(1004,file="input_params.txt",form='formatted')
    READ(unit=1004,nml=time,iostat=err)
    READ(unit=1004,nml=flow,iostat=err)
    READ(unit=1004,nml=ciliaprop,iostat=err)
    READ(unit=1004,nml=particleprop,iostat=err)
    close(1004)

    write(*,'(2(A,I8),2(A,1p1e15.6))') "Nx = ",Nx, " Ny = ",Ny, &
    " Lx = ",Lx, " Ly = ",Ly
    
    write(*,'(2(A,I8),2(A,1p1e15.6))') "nc = ",nc, " np = ",np

    ! Create cilia and particle arrays
    CAP = cilia_array(npparticles,nl,npparticles)
    call create_closed_loop_array(CAP,0.5d0*radius,radius,originP)
    CA = cilia_array(nc,nl,np)
    call create_cilia_array(CA,wbl,dc,dp,origin)

    ! call write_field(u,'u',1) 
    ! call write_field(v,'v',1) 
    ! call write_location_cilia(CA,1,'c')
    ! call write_location_cilia(CAP,1,'p')
    ! call write_location_cilia_force(CA,1,'c')
    ! call write_location_cilia_force(CAP,1,'p')
    ! call write_location_cilia_velocity(CA,1,'c')
    ! call write_location_cilia_velocity(CAP,1,'p')

    call apply_parabolic_initialization(M,u,uleft)

    ! Arrays to transfer data
    BC = [utop,vtop,ubottom,vbottom,uleft,vleft,uright,vright]
    FP = [nu,rho]
    SP = [ko,kd,Rl,Ftip]
    call time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,CAP,A,P,R,tsim,dt,it_save)

    ! call write_location_cilia(CAP,10)
    ! call write_field(u,'u',10) 
    ! call write_field(v,'v',10) 
    call cpu_time(finish)
    print '("Time = ",f15.10," seconds.")',finish-start

end program ibmc
