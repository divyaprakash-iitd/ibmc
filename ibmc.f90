program ibmc
    use iso_fortran_env, only: int32, real32, real64
    use mod_pressure, only: generate_laplacian_sparse, calculate_pressure_sparse
    implicit none

    ! Computational Domain
    real(real32)    :: Lx = 1.0
    real(real32)    :: Ly = 1.0
    integer(int32)  :: Nx = 30
    integer(int32)  :: Ny = 30

    ! Mesh Paramaters
    real(real32) :: dx, dxi
    real(real32) :: dy, dyi

    ! Simulation Paramaters
    real(real32) :: tsim    = 1.0
    real(real32) :: dt      = 0.0005
    real(real32) :: t

    ! Physical Constants
    real(real32) :: nu = 1.0/200.0
    real(real32) :: rho = 1.0d0

    ! Index ranges
    integer(int32) :: imin, imax, jmin, jmax
    integer(int32), dimension(:), allocatable :: iu,ju,iv,jv,ip,jp
    integer(int32), dimension(:), allocatable :: iuc,juc,ivc,jvc,ipc,jpc
    integer(int32) :: i,j

    ! Boundary values
    real(real32) :: utop, vtop, ubottom, vbottom, &
                    uleft, vleft, uright, vright

    ! Matrices to store fields
    real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:), P(:,:), Pv(:)

    ! Limits of arrays
    integer(int32) :: iuS, iuE, juS, juE, ivS, ivE, jvS, jvE, ipS, ipE, jpS, jpE  
    integer(int32) :: iucS, iucE, jucS, jucE, ivcS, ivcE, jvcS, jvcE  

    ! Temporary/Miscellaneous variable
    real(real64) :: ucenter, vcenter
    integer(int32) :: n, NN

    ! Declaration for subroutines
    real(real64), allocatable :: A(:,:,:)

    ! Define ranges
    imin = 2
    imax = imin + Nx - 1
    jmin = 2
    jmax = jmin + Ny - 1

    ! u: u-complete
    iuS = imin
    iuE = imax+1
    juS = jmin-1
    juE = jmax+1

    ! v: v-complete
    ivS = imin-1
    ivE = imax+1
    jvS = jmin
    jvE = jmax+1

    ! uc: u-calculated     
    iucS = imin+1
    iucE = imax
    jucS = jmin
    jucE = jmax

    ! vc: v-calculated
    ivcS = imin
    ivcE = imax
    jvcS = jmin+1
    jvcE = jmax

    ! p: p-complete
    ipS = imin
    ipE = imax
    jpS = jmin
    jpE = jmax

    ! ! Complete index range
    ! iu = [  (i, i = imin      ,   imax+1) ]
    ! ju = [  (i, i = jmin-1    ,   jmax+1) ]
    ! iv = [  (i, i = imin-1    ,   imax+1) ]
    ! jv = [  (i, i = jmin      ,   jmax+1) ]
    ! ip = [  (i, i = imin      ,   imax)   ]
    ! jp = [  (i, i = jmin      ,   jmax)   ]

    ! ! Indices range for unknowns/calculated cells
    ! iuc = [ (i, i = jmin+1  ,   imax)   ]
    ! juc = [ (i, i = jmin    ,   jmax)   ]
    ! ivc = [ (i, i = imin    ,   imax)   ]
    ! jvc = [ (i, i = jmin+1  ,   jmax)   ]
    ! ipc = [ (i, i = imin    ,   imax)   ]
    ! jpc = [ (i, i = jmin    ,   jmax)   ]



    ! Allocate matrices ofr u,v,us,vs,rhs
    allocate(u(iuS:iuE,juS:juE))    
    allocate(v(ivS:ivE,jvS:jvE))
    allocate(us(iuS:iuE,juS:juE))    
    allocate(vs(ivS:ivE,jvS:jvE))
    allocate(R(Nx*Ny))
    allocate(P(ipS:ipE,jpS:jpE))
    allocate(Pv(Nx*Ny))

    ! Allocation for sub-routines
    allocate(A(imin:imax,jmin:jmax,NN))

    ! Initialize
    u   = 0.0d0
    v   = 0.0d0
    us  = 0.0d0
    vs  = 0.0d0
    R   = 0.0d0

    ! Mesh values
    dx  = Lx/Nx
    dy  = Ly/Ny
    dxi = 1/dx
    dyi = 1/dy

    ! Generate Laplacian matrix
    call generate_laplacian_sparse(dxi,dyi,A,imin,imax,jmin,jmax)

    ! Start time loop
    t = 0.0d0
    do while (t.lt.tsim)
        t = t+dt 

        ! Apply boundary conditions
        u(:,juS) = 2*ubottom - u(:,jucS);
        u(:,juE) = 2*utop    - u(:,jucE);
        v(ivS,:) = 2*vleft   - v(ivcS,:);
        v(ivE,:) = 2*vright  - v(ivcE,:);
    

        u(iuS,:)    = uleft;
        u(iuE,:)    = uright;
        v(:,jvE)    = vtop;
        v(:,jvS)    = vbottom;


        ! Perform predictor step
        ! us 
        ! (u-velocity cell)
        do j = jucS,jucE
            do i = iucS,iucE
                vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1));
                us(i,j) = u(i,j) + dt* &
                            ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
                            + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
                            - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
                            - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi)
            end do
        end do

        ! vs 
        ! (v-velocity cell)
        do j = jvcS,jvcE
            do i = ivcS,ivcE
                ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j));
                vs(i,j) = v(i,j) + dt* &
                            ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
                            + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
                            - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
                            - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)
            end do
        end do

        ! Form the right hand side of the pressure poisson equation
        ! (Pressure cell)
        n = 0;
        do j = jpS,jpE
            do i = ipS,ipE
                n = n + 1
                R(n) = -rho/dt* &
                        ( (us(i+1,j) - us(i,j))*dxi &
                        + (vs(i,j+1) - vs(i,j))*dyi)
            end do
        end do

        ! Solve for presssure   
        call calculate_pressure_sparse(imin,imax,jmin,jmax,R,A,P)

        ! ! Convert p to mesh representation
        ! ! (Pressure cell)
        ! n = 0
        ! P = 0.0d0
        ! do j = jpS,jpE
        !     do i = ipS,ipE
        !         n = n + 1;
        !         P(i,j) = Pv(n);
        !     end do
        ! end do

        ! Perform the corrector steps
        ! u 
        ! (u-velocity cell)
        do j = jucS,jucE
            do i = iucS,iucE
                u(i,j) = us(i,j) - dt/rho * (p(i,j) - p(i-1,j)) * dxi;
            end do
        end do
        
        ! v
        ! (v-velocity cell)
        do j = jvcS,jvcE
            do i = ivcS,ivcE
                v(i,j) = vs(i,j) - dt/rho * (p(i,j) - p(i,j-1)) * dyi;
            end do
        end do

    end do

end program ibmc
