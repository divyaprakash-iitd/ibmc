program ibmc
    use iso_fortran_env, only: int32, real32, real64
    use mod_pressure, only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx, only: calculate_pressure_amgx
    use mod_dims, only: dims
    use mod_field, only: Field
    implicit none
    
    ! Execution time
    real(real64) :: start, finish

    ! Define fields
    type(dims) :: ui, uj, vi, vj
    type(Field) :: u, v

    ! Computational Domain
    real(real32)    :: Lx = 1.0
    real(real32)    :: Ly = 1.0
    integer(int32)  :: Nx = 3
    integer(int32)  :: Ny = 3

    ! Mesh Paramaters
    real(real32) :: dx, dxi
    real(real32) :: dy, dyi

    ! Simulation Paramaters
    real(real32) :: tsim    = 20
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
    ! real(real64), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), P(:,:), A(:,:,:)

    ! Limits of arrays
    integer(int32) :: iuS, iuE, juS, juE, ivS, ivE, jvS, jvE, ipS, ipE, jpS, jpE  
    integer(int32) :: iucS, iucE, jucS, jucE, ivcS, ivcE, jvcS, jvcE  

    ! Temporary/Miscellaneous variable
    ! real(real64) :: ucenter, vcenter
    ! integer(int32) :: NN
    ! logical :: init_status

    call cpu_time(start)

    ! Define boundary conditions
    utop    = 1.0
    vtop    = 0.0
    ubottom = 0.0
    vbottom = 0.0
    uleft   = 0.0
    vleft   = 0.0
    uright  = 0.0
    vright  = 0.0

    ! Define ranges
    imin = 2
    imax = imin + Nx - 1
    jmin = 2
    jmax = jmin + Ny - 1

    ! u: u-complete
    ! iuS = imin
    ! iuE = imax+1
    ! juS = jmin-1
    ! juE = jmax+1

    ! v: v-complete
    ! ivS = imin-1
    ! ivE = imax+1
    ! jvS = jmin
    ! jvE = jmax+1

    ! uc: u-calculated     
    ! ! iucS = imin+1
    ! ! iucE = imax
    ! ! jucS = jmin
    ! ! jucE = jmax

    ! ! ! vc: v-calculated
    ! ! ivcS = imin
    ! ! ivcE = imax
    ! ! jvcS = jmin+1
    ! ! jvcE = jmax

    ! ! p: p-complete
    ! ipS = imin
    ! ipE = imax
    ! jpS = jmin
    ! jpE = jmax

    ! Allocate matrices ofr u,v,us,vs,rhs
    ! allocate(u(iuS:iuE,juS:juE))    
    ! allocate(v(ivS:ivE,jvS:jvE))
    ! allocate(us(iuS:iuE,juS:juE))    
    ! allocate(vs(ivS:ivE,jvS:jvE))
    ! allocate(R(ipS:ipE,jpS:jpE))
    ! allocate(P(ipS:ipE,jpS:jpE))
    ! NN = 5
    ! allocate(A(imin:imax,jmin:jmax,NN))

    ! Initialize
    ! u   = 0.0d0
    ! v   = 0.0d0
    ! us  = 0.0d0
    ! vs  = 0.0d0
    ! R   = 0.0d0
    ! A   = 0.0d0

    ! Mesh values
    dx  = Lx/Nx
    dy  = Ly/Ny
    dxi = 1.0/dx
    dyi = 1.0/dy

    ui = dims('x',imin,imax+1)
    uj = dims('y',jmin-1,jmax+1)
    vi = dims('x',imin-1,imax+1)
    vj = dims('y',jmin,jmax+1)

    u = Field('u',ui,uj)
    v = Field('v',vi,vj)

    ! Generate Laplacian matrix
    ! call generate_laplacian_sparse(A,dxi,dyi)
    ! Start time loop
    t = 0.0d0
    ! init_status = .False.
    do while (t.lt.tsim)
        t = t+dt 

        ! Apply boundary conditions

        u%data(:,u%y%lb) = 2*ubottom - u%data(:,u%y%lb+1)
        u%data(:,u%y%ub) = 2*utop - u%data(:,u%y%ub-1)
        u%data(u%x%lb,:) = uleft
        u%data(u%x%ub,:) = uright

        v%data(v%x%lb,:) = 2*vleft - v%data(v%x%lb+1,:)
        v%data(v%x%ub,:) = 2*vright - v%data(v%x%ub-1,:)
        v%data(:,v%y%lb) = vbottom
        v%data(:,v%y%ub) = vtop

        ! call u%write_field(10)


        ! Predictor Step
        u%data = u%data + dt* (nu*(u%diff2x()/dxi**2 + u%diff2y()/dyi**2) - (u%diffx()/dx + u%diffy()/dy))

    !     ! Perform predictor step
    !     ! us 
    !     ! (u-velocity cell)
    !     do concurrent (j = jucS:jucE, i = iucS:iucE)
    !         vcenter = 0.25*(v(i-1,j) + v(i-1,j+1) + v(i,j) + v(i,j+1))
    !         us(i,j) = u(i,j) + dt* &
    !                     ( nu*(u(i-1,j) - 2*u(i,j) + u(i+1,j))*dxi**2 &
    !                     + nu*(u(i,j-1) -2*u(i,j) + u(i,j+1))*dyi**2 &
    !                     - u(i,j)*(u(i+1,j) - u(i-1,j))*0.5*dxi &
    !                     - vcenter*(u(i,j+1)-u(i,j-1))*0.5*dyi)
    !     end do

    !     ! vs 
    !     ! (v-velocity cell)
    !     do concurrent (j = jvcS:jvcE, i = ivcS:ivcE)
    !         ucenter = 0.25*(u(i,j-1)+u(i,j)+u(i+1,j-1)+u(i+1,j))
    !         vs(i,j) = v(i,j) + dt* &
    !                     ( nu*(v(i-1,j) - 2*v(i,j) + v(i+1,j))*dxi**2 &
    !                     + nu*(v(i,j-1) - 2*v(i,j) + v(i,j+1))*dyi**2 &
    !                     - ucenter*(v(i+1,j) - v(i-1,j))*0.5*dxi &
    !                     - v(i,j)*(v(i,j+1)-v(i,j-1))*0.5*dyi)
    !     end do

    !     ! Form the right hand side of the pressure poisson equation
    !     ! (Pressure cell)
    !     do concurrent (j = jpS:jpE, i = ipS:ipE)
    !         R(i,j) = -rho/dt* &
    !                 ( (us(i+1,j) - us(i,j))*dxi &
    !                 + (vs(i,j+1) - vs(i,j))*dyi)
    !     end do

    !     ! Solve for presssure
    !     call calculate_pressure_sparse(A,P,R)

    !     !call calculate_pressure_amgx(A,P,R,init_status)

    !     ! Perform the corrector steps
    !     ! u 
    !     ! (u-velocity cell)
    !     do concurrent (j = jucS:jucE, i = iucS:iucE)
    !         u(i,j) = us(i,j) - dt/rho * (p(i,j) - p(i-1,j)) * dxi
    !     end do
        
    !     ! v
    !     ! (v-velocity cell)
    !     do concurrent (j = jvcS:jvcE, i = ivcS:ivcE)
    !         v(i,j) = vs(i,j) - dt/rho * (p(i,j) - p(i,j-1)) * dyi
    !     end do

    !     print *, 'time = ', t
    ! end do

    ! call cpu_time(finish)
    ! print '("Time = ",f6.3," seconds.")',finish-start 

    ! open(unit=55, file='u.txt', ACTION="write", STATUS="replace")
    ! do i=imin,imax+1
    !     write(55, '(*(F14.7))')( real(u(i,j)) ,j=jmin-1,jmax+1)
    ! end do

    ! open(unit=65, file='v.txt', ACTION="write", STATUS="replace")
    ! do i=imin-1,imax+1
    !     write(65, '(*(F14.7))')( real(v(i,j)) ,j=jmin,jmax+1)
    ! end do

    ! ! open(unit=75, file='p.txt', ACTION="write", STATUS="replace")
    ! ! do i=imin-1,imax+1
    ! !     write(75, '(*(F14.7))')( real(p(i,j)) ,j=jmin-1,jmax+1)
    end do

end program ibmc
