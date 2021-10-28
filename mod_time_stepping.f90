module mod_time_stepping
    use iso_fortran_env, only: int32, real64
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse
    use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    use mod_boundary
    use mod_io
    use mod_ib
    use mod_ibm
    use mod_vec
    use mod_cilia
    use mod_cilia_array

    implicit none
    
contains

    pure subroutine euler(M,u,v,us,vs,nu,dt,Fx,Fy)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu, dt
    
        
        call cdu(M,u,v,us,nu,Fx)
        call cdv(M,u,v,vs,nu,Fy) 

        us = u + dt*us
        vs = v + dt*vs

    end subroutine euler
    
    subroutine RK2(M,u,v,us,vs,nu,dt,Fx,Fy)
        class(mesh), intent(in) :: M
        real(real64), intent(inout) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu, dt

        real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: k1u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k1v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        
        real(real64) :: umid(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vmid(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        
        ! ! Step 1
        ! call cdu(M,u,v,k1u,nu,Fx)
        ! du1 = k1u*dt/2
        ! call cdv(M,u,v,k1v,nu,Fy)
        ! dv1 = k1v*dt/2

        ! u = u + du1
        ! v = v + dv1

        ! ! print *, 'Sum = ', sum(u)

        ! ! Step 2
        ! call cdu(M,u,v,k2u,nu,Fx)
        ! call cdv(M,u,v,k2v,nu,Fy)

        ! du2 =k2u*dt
        ! dv2 =k2v*dt

        ! us = u + du2
        ! vs = v + dv2

        ! call cdu(M,u,v,us,nu,Fx)
        ! call cdv(M,u,v,vs,nu,Fy)
        
        umid = u + dt/2*cdu_f(M,u,v,nu,Fx)
        vmid = v + dt/2*cdv_f(M,u,v,nu,Fy)
        
        us = u + dt*cdu_f(M,umid,vmid,nu,Fx)
        vs = v + dt*cdv_f(M,umid,vmid,nu,Fx)

    end subroutine RK2

    subroutine RK4(M,u,v,us,vs,nu,dt,Fx,Fy)
        class(mesh), intent(in) :: M
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in out) :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(in) :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in) :: nu, dt

        real(real64) :: du1(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv1(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du2(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv2(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: du3(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), dv3(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: k1u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k1v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: k2u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k2v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: k3u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k3v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64) :: k4u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), k4v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Step 1
        call cdu(M,u,v,k1u,nu,Fx)
        du1 = k1u*dt/2
        call cdv(M,u,v,k1v,nu,Fy)
        dv1 = k1v*dt/2
        
        ! Step 2
        call cdu(M,(u+du1),(v+dv1),k2u,nu,Fx)
        du2 = k2u*dt/2
        call cdv(M,(u+du1),(v+dv1),k2v,nu,Fy)
        dv2 = k2v*dt/2

        ! Step 3
        call cdu(M,(u+du2),(v+dv2),k3u,nu,Fx)
        du3 = dt*k3u
        call cdv(M,(u+du2),(v+dv2),k3v,nu,Fy)
        dv3 = dt*k3v

        ! Step 4
        call cdu(M,(u+du3),(v+dv3),k4u,nu,Fx)
        call cdv(M,(u+du3),(v+dv3),k4v,nu,Fy)

        us = u + 1.0d0/6*dt*(k1u + 2*k2u + 2*k3u + k4u)
        vs = v + 1.0d0/6*dt*(k1v + 2*k2v + 2*k3v + k4v)

        ! us = us/2
        ! vs = vs/2

        ! call cdu(M,u,v,us,nu)
        ! call cdv(M,u,v,vs,nu) 

        ! us = u + dt*us
        ! vs = v + dt*vs
    end subroutine RK4

!    subroutine time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,A,P,R,tsim,dt)
!        real(real64), intent(in)          :: FP(:)
!        real(real64), intent(in)          :: BC(:)
!        class(mesh), intent(inout)        :: M
!        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
!        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
!        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
!        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
!        real(real64), intent(in)          :: SP(:)
!        class(cilia_array), intent(inout) :: CA
!        real(real64), intent(in)          :: A(:,:,:)
!        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
!        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
!        real(real64), intent(in)          :: tsim
!        real(real64), intent(in)          :: dt


!        ! Iteration 
!        integer(int32) :: it

!        ! Time
!        real(real64) :: t

!        ! Fluid properties
!        real(real64) :: nu, rho        

!        ! Boundary conditions
!        real(real64)    :: utop, vtop, ubottom, vbottom, &
!                               uleft, vleft, uright, vright
!        ! Spring properties
!        real(real64)    :: ks, Rl, Ftip
       
!        ! Assign boundary values
!        utop    = BC(1)
!        vtop    = BC(2)
!        ubottom = BC(3)
!        vbottom = BC(4)
!        uleft   = BC(5)
!        vleft   = BC(6)
!        uright  = BC(7)
!        vright  = BC(8)

!        ! Get spring parameters
!        ks      = SP(1)
!        Rl      = SP(2)
!        Ftip    = SP(3)

!        ! Get fluid properties
!        nu  = FP(1)
!        rho = FP(2)

!        ! Start time loop
!        t = 0.0d0
!        it = 0

!        do while (t.lt.tsim)
!            t = t + dt
!            it = it + 1

!            ! Apply velocity boundary conditions
!            call apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

!            ! Calculate forces in the cilia
!            ! call calculate_cilia_array_force(CA,ks,Rl)

!            ! Apply tip force
!            ! call apply_tip_force_cilia_array(CA,Ftip,t)

!            ! Spread force from the cilia to the fluid
!            Fx = 0.0d0
!            Fy = 0.0d0
!            ! call spread_force_cilia_array(M,CA,Fx,Fy)

!            ! Calculate intermediate velocity (Implement RK4 here)
!            call euler(M,u,v,us,vs,nu,dt,Fx,Fy)

!            ! Form the RHS of the pressure poisson equation
!            call calculate_rhs(M,us,vs,R,rho,dt)

!            ! Solve for pressure
!            call calculate_pressure_sparse(A,P,R)

!            ! Perform the corrector step to obtain the velocity
!            call corrector(M,u,v,us,vs,p,rho,dt)

!            ! Interpolate the fluid velocity to the cilia
!            ! call initialize_velocity_cilia_array(CA)
!            ! call interpolate_velocity_cilia_array(M,CA,u,v)

!            ! Update the cilia node locations
!            ! call update_cilia_array(CA,dt)

!            print *, 'time = ', t
           
!            ! Write files every 10th timestep
!            if (mod(it,50).eq.0) then 
!                call write_field(u,'u',it) 
!                call write_field(v,'v',it) 
!                ! call write_location(CA,it)
!            end if

!        end do
!    end subroutine time_loop

    subroutine time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,A,P,R,tsim,dt)
        real(real64), intent(in)          :: FP(:)
        real(real64), intent(in)          :: BC(:)
        class(mesh), intent(inout)        :: M
        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in)          :: SP(:)
        class(cilia_array), intent(inout) :: CA
        ! class(cilia_array), intent(inout) :: CAmid
        real(real64), intent(in)          :: A(:,:,:)
        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in)          :: tsim
        real(real64), intent(in)          :: dt

        ! Intermediate cilia
        type(cilia_array) :: CAmid

        ! Intermediate velocities
        real(real64)       :: umid(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vmid(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Iteration 
        integer(int32) :: it

        ! Time
        real(real64) :: t

        ! Fluid properties
        real(real64) :: nu, rho        

        ! Boundary conditions
        real(real64)    :: utop, vtop, ubottom, vbottom, &
                               uleft, vleft, uright, vright
        ! Spring properties
        real(real64)    :: ks, Rl, Ftip
        
        ! AmgX
        logical         :: init_status

        ! Create tempeorary cilia array
        CAmid = cilia_array(CA%nc,CA%array(1)%nl,CA%array(1)%np)


        ! Assign boundary values
        utop    = BC(1)
        vtop    = BC(2)
        ubottom = BC(3)
        vbottom = BC(4)
        uleft   = BC(5)
        vleft   = BC(6)
        uright  = BC(7)
        vright  = BC(8)

        ! Get spring parameters
        ks      = SP(1)
        Rl      = SP(2)
        Ftip    = SP(3)

        ! Get fluid properties
        nu  = FP(1)
        rho = FP(2)

        ! Start time loop
        t = 0.0d0
        it = 0
    
        init_status = .False. ! AmgX initialization status

        do while (t.lt.tsim)
            t = t + dt
            it = it + 1
            

            ! Calculate forces in the immersed boundary structure
            call calculate_cilia_array_force(CA,ks,Rl)

            ! Apply tip force for the first 1 second
            ! if (t.lt.0.5) then
            call apply_tip_force_cilia_array(CA,Ftip,t)
            ! end if
            
            call copy_cilia(CA,CAmid)
 
            ! RK2: Step 1
            ! Apply velocity boundary conditions
            call apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)

            ! us = u + 0.5d0*dt*cdu_f(M,u,v,nu,Fx)
            ! vs = v + 0.5d0*dt*cdv_f(M,u,v,nu,Fy)
            call cdu(M,u,v,us,nu,Fx)
            us = u + us*dt/2
            call cdv(M,u,v,vs,nu,Fy)
            vs = v + vs*dt/2

            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,0.5d0*dt)

            ! Solve for pressure
            ! call calculate_pressure_sparse(A,P,R)
            call calculate_pressure_amgx(A,P,R,init_status)

            ! Perform the corrector step to obtain the velocity
            call corrector(M,umid,vmid,us,vs,P,rho,0.5d0*dt)

            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CAmid)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CAmid,umid,vmid)

            ! Update the Immersed Boundary
            call update_cilia_array(CAmid,dt/2)

            ! RK2: Step 2
            
            ! Calculate forces in the immersed boundary structure
            call calculate_cilia_array_force(CAmid,ks,Rl)

            ! Apply tip force for the first 1 second
            ! if (t.lt.0.5) then
            call apply_tip_force_cilia_array(CAmid,Ftip,t)
            ! end if

            call apply_boundary(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)
            
            ! us = u + dt*cdu_f(M,umid,vmid,nu,Fx)
            ! vs = v + dt*cdv_f(M,umid,vmid,nu,Fx)

            call cdu(M,u,v,us,nu,Fx)
            us = u + us*dt
            call cdv(M,u,v,vs,nu,Fy)
            vs = v + vs*dt
            
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,dt)

            ! Solve for pressure
            ! call calculate_pressure_sparse(A,P,R)
            call calculate_pressure_amgx(A,P,R,init_status)

            ! Perform the corrector step to obtain the velocity
            call corrector(M,u,v,us,vs,p,rho,dt)

            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CA)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CA,u,v)

            ! Update the Immersed Boundary
            call update_cilia_array(CA,dt)

            print *, 'time = ', t
            
            ! Write files every Nth timestep
            if (mod(it,25).eq.0) then 
                call write_field(u,'u',it) 
                call write_field(v,'v',it) 
                call write_location_cilia(CA,it)
            end if

        end do
    end subroutine time_loop

    subroutine copy_cilia(CA1,CA2)
        class(cilia_array), intent(in) :: CA1
        class(cilia_array), intent(inout) :: CA2

        integer(int32) :: nc, nl, np
        integer(int32) :: ic, il, ip

        nc = CA1%nc
        nl = CA1%array(1)%nl
        np = CA1%array(1)%layers(1)%np

        do ic = 1,nc
            do il = 1,nl
                do ip = 1,np
                    CA2%array(ic)%layers(il)%boundary(ip)%x = CA1%array(ic)%layers(il)%boundary(ip)%x
                    CA2%array(ic)%layers(il)%boundary(ip)%y = CA1%array(ic)%layers(il)%boundary(ip)%y
                    CA2%array(ic)%layers(il)%boundary(ip)%Ux = CA1%array(ic)%layers(il)%boundary(ip)%Ux
                    CA2%array(ic)%layers(il)%boundary(ip)%Uy = CA1%array(ic)%layers(il)%boundary(ip)%Uy
                    CA2%array(ic)%layers(il)%boundary(ip)%Fx = CA1%array(ic)%layers(il)%boundary(ip)%Fx
                    CA2%array(ic)%layers(il)%boundary(ip)%Fy = CA1%array(ic)%layers(il)%boundary(ip)%Fy
                    ! CA2%array(ic)%layers(il)%boundary(ip)%tag = CA1%array(ic)%layers(il)%boundary(ip)%tag
                end do
            end do
        end do

    end subroutine copy_cilia

end module mod_time_stepping

    !     do while (t.lt.tsim)
    !         t = t + dt
    !         it = it + 1

    !         ! RK2: Step 1
    !         ! Apply velocity boundary conditions
    !         call apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

    !         ! Calculate forces in the cilia
    !         ! call calculate_cilia_array_force(CA,ks,Rl)

    !         ! Apply tip force
    !         ! call apply_tip_force_cilia_array(CA,Ftip,t)

    !         ! Spread force from the cilia to the fluid
    !         Fx = 0.0d0
    !         Fy = 0.0d0
    !         ! call spread_force_cilia_array(M,CA,Fx,Fy)

    !         ! Calculate intermediate velocity (Implement RK4 here)
    !         ! call euler(M,u,v,us,vs,nu,dt,Fx,Fy)
    !         ! call RK2(M,u,v,us,vs,nu,dt,Fx,Fy)
    !         ! call RK4(M,u,v,us,vs,nu,dt,Fx,Fy)

    !         us = u + dt/2*cdu_f(M,u,v,nu,Fx)
    !         vs = v + dt/2*cdv_f(M,u,v,nu,Fy)

    !         ! Form the RHS of the pressure poisson equation
    !         call calculate_rhs(M,us,vs,R,rho,dt/2)

    !         ! Solve for pressure
    !         call calculate_pressure_sparse(A,P,R)
    !         ! call calculate_pressure_amgx(A,P,R,init_status)

    !         ! Perform the corrector step to obtain the velocity
    !         call corrector(M,umid,vmid,us,vs,p,rho,dt/2)

    !         ! print *, sum(r)
    !         ! Interpolate the fluid velocity to the cilia
    !         ! call initialize_velocity_cilia_array(CA)
    !         ! call interpolate_velocity_cilia_array(M,CA,u,v)

    !         ! Update the cilia node locations
    !         ! call update_cilia_array(CA,dt)

    !         ! RK2: Step 2

    !         call apply_boundary(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            
    !         us = u + dt*cdu_f(M,umid,vmid,nu,Fx)
    !         vs = v + dt*cdv_f(M,umid,vmid,nu,Fx)

    !         ! Form the RHS of the pressure poisson equation
    !         call calculate_rhs(M,us,vs,R,rho,dt)

    !         ! Solve for pressure
    !         call calculate_pressure_sparse(A,P,R)
    !         ! call calculate_pressure_amgx(A,P,R,init_status)

    !         ! Perform the corrector step to obtain the velocity
    !         call corrector(M,u,v,us,vs,p,rho,dt)

    !         ! print *, sum(r)
    !         ! Interpolate the fluid velocity to the cilia
    !         ! call initialize_velocity_cilia_array(CA)
    !         ! call interpolate_velocity_cilia_array(M,CA,u,v)

    !         ! Update the cilia node locations
    !         ! call update_cilia_array(CA,dt)

    !         print *, 'time = ', t
            
    !         ! Write files every 10th timestep
    !         if (mod(it,50).eq.0) then 
    !             call write_field(u,'u',it) 
    !             call write_field(v,'v',it) 
    !             ! call write_location(CA,it)
    !         end if

    !     end do
    ! end subroutine time_loop
