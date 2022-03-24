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
    use mod_closed_cilia
    use nvtx
    use mod_cell
    use mod_inter_particle_force
    implicit none
    
contains

    subroutine time_loop_particle(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,CAP,A,P,R,tsim,dt,it_save)
        real(real64), intent(in)          :: FP(:)
        real(real64), intent(in)          :: BC(:)
        class(mesh), intent(inout)        :: M
        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in)          :: SP(:)
        class(cilia_array), intent(inout) :: CA
        class(cilia_array), intent(inout) :: CAP
        ! class(cilia_array), intent(inout) :: CAmid
        real(real64), intent(in)          :: A(:,:,:)
        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in)          :: tsim
        real(real64), intent(in)          :: dt
        integer(int32), intent(in)        :: it_save
        ! Intermediate cilia
        type(cilia_array) :: CAmid, CAPmid

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
                               uleft, vleft, uright, vright, ulid
        ! Spring properties
        real(real64)    :: ko,kd,Rl,tp,kop,kod
        
        ! AmgX
        logical         :: init_status

        ! Resting lengths
        real(real64) :: RLH, RLD, RLV
    
        ! Create an array of cells
        type(cell), allocatable, target :: cell_array(:,:)
        integer(int32) :: Nxcell, Nycell

        RLH = 0.0; RLD = 0.0; RLV = 0.0


        ! Create temporary cilia array
        CAmid = cilia_array(CA%nc,CA%array(1)%nl,CA%array(1)%np)
        CAPmid = cilia_array(CAP%nc,CAP%array(1)%nl,CAP%array(1)%np)
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
        ko      = SP(1)
        kd      = SP(2)
        kop     = SP(3)
        kod     = SP(4)
        Rl      = SP(5)
        tp    = SP(6)

        ! Get fluid properties
        nu  = FP(1)
        rho = FP(2)

        ! Start time loop
        t = 0.0d0
        it = 0
    
        init_status = .False. ! AmgX initialization status

        ! Neghbour's list computations
        Nxcell = 10
        Nycell = 10
        allocate(cell_array(Nxcell,Nycell))
        ! Create cells to assign particles to it
        call create_cells(M,cell_array)

        ! Lid velocity
        ulid = utop

        do while (t.lt.tsim)
            ! ulid is the amplitude
            utop = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
            !utop = sin(2*3.1415*t/tp) * ulid
            call nvtxStartRange("Time Loop")
            t = t + dt
            it = it + 1
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("RK2:Step-1")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            ! call assign_particles_to_cells(cell_array,CA)
            call assign_particles_to_cells(cell_array,CAP)
            call nvtxEndRange
           
            call nvtxStartRange("Calculate cilia forces")
            ! call calculate_cilia_array_force(CA,ko,kd,Rl)
            call calculate_closed_loop_array_force(CAP,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange
           
            call nvtxStartRange("Copy cilia")
            ! call copy_cilia(CA,CAmid)
            call copy_cilia(CAP,CAPmid)
            call nvtxEndRange
 
            ! RK2: Step 1
            ! Apply velocity boundary conditions
            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange
        
            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            ! call spread_force_cilia_array(M,CAmid,Fx,Fy)
            call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,u,v,us,nu,Fx)
            us = u + us*dt/2
            call cdv(M,u,v,vs,nu,Fy)
            vs = v + vs*dt/2
            call nvtxEndRange

            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange


            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,umid,vmid,us,vs,P,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            ! call initialize_velocity_cilia_array(CAmid)
            call initialize_velocity_cilia_array(CAPmid)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            ! call interpolate_velocity_cilia_array(M,CAmid,umid,vmid)
            call interpolate_velocity_cilia_array(M,CAPmid,umid,vmid)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            ! call update_cilia_array(CAmid,dt/2)
            call update_cilia_array(CAPmid,dt/2)
            call nvtxEndrange


            call nvtxEndRange
            ! RK2: Step 2
            
            call nvtxStartRange("RK2:Step-2")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            ! call assign_particles_to_cells(cell_array,CAmid)
            call assign_particles_to_cells(cell_array,CAPmid)
            call nvtxEndRange
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("Calculate cilia forces")
            ! call calculate_cilia_array_force(CAmid,ko,kd,Rl)
            call calculate_closed_loop_array_force(CAPmid,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange

            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            ! call spread_force_cilia_array(M,CAmid,Fx,Fy)
            call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,umid,vmid,us,nu,Fx)
            us = u + us*dt
            call cdv(M,umid,vmid,vs,nu,Fy)
            vs = v + vs*dt
            call nvtxEndRange
            
            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,u,v,us,vs,p,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            ! call initialize_velocity_cilia_array(CA)
            call initialize_velocity_cilia_array(CAP)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            ! call interpolate_velocity_cilia_array(M,CA,u,v)
            call interpolate_velocity_cilia_array(M,CAP,u,v)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            ! call update_cilia_array(CA,dt)
            call update_cilia_array(CAP,dt)
            call nvtxEndRange
        
            call nvtxEndRange
            write(*,'(A,F10.5)') 'time = ', t
            
        call nvtxEndRange
            ! Write files every Nth timestep
            if (mod(it,it_save).eq.0) then 
                call write_field(u,'u',it) 
                call write_field(v,'v',it) 
                ! call write_location_cilia(CA,it,'c')
                call write_location_cilia(CAP,it,'p')
                ! call write_location_cilia_force(CA,it,'c')
                call write_location_cilia_force(CAP,it,'p')
                ! call write_location_cilia_velocity(CA,it,'c')
                call write_location_cilia_velocity(CAP,it,'p')
                call write_cell_data(cell_array)
                call write_particle_data(cell_array)
            end if

        end do
    end subroutine time_loop_particle

    subroutine time_loop_cilia(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,CAP,A,P,R,tsim,dt,it_save)
        real(real64), intent(in)          :: FP(:)
        real(real64), intent(in)          :: BC(:)
        class(mesh), intent(inout)        :: M
        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in)          :: SP(:)
        class(cilia_array), intent(inout) :: CA
        class(cilia_array), intent(inout) :: CAP
        ! class(cilia_array), intent(inout) :: CAmid
        real(real64), intent(in)          :: A(:,:,:)
        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in)          :: tsim
        real(real64), intent(in)          :: dt
        integer(int32), intent(in)        :: it_save
        ! Intermediate cilia
        type(cilia_array) :: CAmid, CAPmid

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
                               uleft, vleft, uright, vright, ulid
        ! Spring properties
        real(real64)    :: ko,kd,Rl,tp,kop,kod
        
        ! AmgX
        logical         :: init_status

        ! Resting lengths
        real(real64) :: RLH, RLD, RLV
    
        ! Create an array of cells
        type(cell), allocatable, target :: cell_array(:,:)
        integer(int32) :: Nxcell, Nycell

        RLH = 0.0; RLD = 0.0; RLV = 0.0


        ! Create temporary cilia array
        CAmid = cilia_array(CA%nc,CA%array(1)%nl,CA%array(1)%np)
        CAPmid = cilia_array(CAP%nc,CAP%array(1)%nl,CAP%array(1)%np)
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
        ko      = SP(1)
        kd      = SP(2)
        kop     = SP(3)
        kod     = SP(4)
        Rl      = SP(5)
        tp    = SP(6)

        ! Get fluid properties
        nu  = FP(1)
        rho = FP(2)

        ! Start time loop
        t = 0.0d0
        it = 0
    
        init_status = .False. ! AmgX initialization status

        ! Neghbour's list computations
        Nxcell = 10
        Nycell = 10
        allocate(cell_array(Nxcell,Nycell))
        ! Create cells to assign particles to it
        call create_cells(M,cell_array)

        ! Lid velocity
        ulid = utop

        do while (t.lt.tsim)
            ! ulid is the amplitude
            utop = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
            !utop = sin(2*3.1415*t/tp) * ulid
            call nvtxStartRange("Time Loop")
            t = t + dt
            it = it + 1
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("RK2:Step-1")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            call assign_particles_to_cells(cell_array,CA)
            ! call assign_particles_to_cells(cell_array,CAP)
            call nvtxEndRange
           
            call nvtxStartRange("Calculate cilia forces")
            call calculate_cilia_array_force(CA,ko,kd,Rl)
            ! call calculate_closed_loop_array_force(CAP,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange
           
            call nvtxStartRange("Copy cilia")
            call copy_cilia(CA,CAmid)
            ! call copy_cilia(CAP,CAPmid)
            call nvtxEndRange
 
            ! RK2: Step 1
            ! Apply velocity boundary conditions
            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange
        
            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)
            ! call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,u,v,us,nu,Fx)
            us = u + us*dt/2
            call cdv(M,u,v,vs,nu,Fy)
            vs = v + vs*dt/2
            call nvtxEndRange

            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange


            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,umid,vmid,us,vs,P,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CAmid)
            ! call initialize_velocity_cilia_array(CAPmid)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CAmid,umid,vmid)
            ! call interpolate_velocity_cilia_array(M,CAPmid,umid,vmid)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            call update_cilia_array(CAmid,dt/2)
            ! call update_cilia_array(CAPmid,dt/2)
            call nvtxEndrange


            call nvtxEndRange
            ! RK2: Step 2
            
            call nvtxStartRange("RK2:Step-2")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            call assign_particles_to_cells(cell_array,CAmid)
            ! call assign_particles_to_cells(cell_array,CAPmid)
            call nvtxEndRange
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("Calculate cilia forces")
            call calculate_cilia_array_force(CAmid,ko,kd,Rl)
            ! call calculate_closed_loop_array_force(CAPmid,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange

            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)
            ! call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,umid,vmid,us,nu,Fx)
            us = u + us*dt
            call cdv(M,umid,vmid,vs,nu,Fy)
            vs = v + vs*dt
            call nvtxEndRange
            
            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,u,v,us,vs,p,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CA)
            ! call initialize_velocity_cilia_array(CAP)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CA,u,v)
            ! call interpolate_velocity_cilia_array(M,CAP,u,v)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            call update_cilia_array(CA,dt)
            ! call update_cilia_array(CAP,dt)
            call nvtxEndRange
        
            call nvtxEndRange
            write(*,'(A,F10.5)') 'time = ', t
            
        call nvtxEndRange
            ! Write files every Nth timestep
            if (mod(it,it_save).eq.0) then 
                call write_field(u,'u',it) 
                call write_field(v,'v',it) 
                call write_location_cilia(CA,it,'c')
                ! call write_location_cilia(CAP,it,'p')
                call write_location_cilia_force(CA,it,'c')
                ! call write_location_cilia_force(CAP,it,'p')
                call write_location_cilia_velocity(CA,it,'c')
                ! call write_location_cilia_velocity(CAP,it,'p')
                call write_cell_data(cell_array)
                call write_particle_data(cell_array)
            end if

        end do
    end subroutine time_loop_cilia

    subroutine time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,SP,CA,CAP,A,P,R,tsim,dt,it_save)
        real(real64), intent(in)          :: FP(:)
        real(real64), intent(in)          :: BC(:)
        class(mesh), intent(inout)        :: M
        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in)          :: SP(:)
        class(cilia_array), intent(inout) :: CA
        class(cilia_array), intent(inout) :: CAP
        ! class(cilia_array), intent(inout) :: CAmid
        real(real64), intent(in)          :: A(:,:,:)
        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in)          :: tsim
        real(real64), intent(in)          :: dt
        integer(int32), intent(in)        :: it_save
        ! Intermediate cilia
        type(cilia_array) :: CAmid, CAPmid

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
                               uleft, vleft, uright, vright, ulid
        ! Spring properties
        real(real64)    :: ko,kd,Rl,tp,kop,kod
        
        ! AmgX
        logical         :: init_status

        ! Resting lengths
        real(real64) :: RLH, RLD, RLV

        ! Deformation force
        real(real64) :: fdef

        ! Create an array of cells
        type(cell), allocatable, target :: cell_array(:,:)
        integer(int32) :: Nxcell, Nycell

        RLH = 0.0; RLD = 0.0; RLV = 0.0


        ! Create temporary cilia array
        CAmid = cilia_array(CA%nc,CA%array(1)%nl,CA%array(1)%np)
        CAPmid = cilia_array(CAP%nc,CAP%array(1)%nl,CAP%array(1)%np)
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
        ko      = SP(1)
        kd      = SP(2)
        kop     = SP(3)
        kod     = SP(4)
        Rl      = SP(5)
        tp    = SP(6)

        ! Get fluid properties
        nu  = FP(1)
        rho = FP(2)

        ! Start time loop
        t = 0.0d0
        it = 0
    
        init_status = .False. ! AmgX initialization status

        ! Neghbour's list computations
        Nxcell = 10
        Nycell = 10
        allocate(cell_array(Nxcell,Nycell))
        ! Create cells to assign particles to it
        call create_cells(M,cell_array)

        ! Lid velocity
        ulid = utop

        do while (t.lt.tsim)
            ! ulid is the amplitude
            utop = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
            
            fdef = 2*3.1415/tp*0.0000001*cos(2*3.1415*t/tp)
            !utop = sin(2*3.1415*t/tp) * ulid
            call nvtxStartRange("Time Loop")
            t = t + dt
            it = it + 1
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("RK2:Step-1")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            call assign_particles_to_cells(cell_array,CA)
            call assign_particles_to_cells(cell_array,CAP)
            call nvtxEndRange
           
            call nvtxStartRange("Calculate cilia forces")
            call calculate_cilia_array_force(CA,ko,kd,Rl)
            call calculate_closed_loop_array_force(CAP,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange

            ! Apply particle deformation forces
            call dynamic_deform(CAP,fdef)
           
            call nvtxStartRange("Copy cilia")
            call copy_cilia(CA,CAmid)
            call copy_cilia(CAP,CAPmid)
            call nvtxEndRange
 
            ! RK2: Step 1
            ! Apply velocity boundary conditions
            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange
        
            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)
            call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,u,v,us,nu,Fx)
            us = u + us*dt/2
            call cdv(M,u,v,vs,nu,Fy)
            vs = v + vs*dt/2
            call nvtxEndRange

            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange


            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,umid,vmid,us,vs,P,rho,0.5d0*dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CAmid)
            call initialize_velocity_cilia_array(CAPmid)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CAmid,umid,vmid)
            call interpolate_velocity_cilia_array(M,CAPmid,umid,vmid)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            call update_cilia_array(CAmid,dt/2)
            call update_cilia_array(CAPmid,dt/2)
            call nvtxEndrange


            call nvtxEndRange
            ! RK2: Step 2
            
            call nvtxStartRange("RK2:Step-2")
            
            call nvtxStartRange("Assign particles to cells")
            ! Reinitialize the no. of particles before assignment
            cell_array%NN = 0
            ! Assign particles to cells
            call assign_particles_to_cells(cell_array,CAmid)
            call assign_particles_to_cells(cell_array,CAPmid)
            call nvtxEndRange
            
            ! Calculate forces in the immersed boundary structure
            call nvtxStartRange("Calculate cilia forces")
            call calculate_cilia_array_force(CAmid,ko,kd,Rl)
            call calculate_closed_loop_array_force(CAPmid,kop,kod,RLV,RLH,RLD)
            call nvtxEndRange
            
            call nvtxStartRange("Calculate inter-particle forces")
            call calculate_forces(cell_array,M)
            call nvtxEndRange

            ! Apply particle deformation forces
            call dynamic_deform(CAPmid,fdef)
            
            call nvtxStartRange("Apply BC")
            call apply_boundary_channel(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,u,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Spread Forces")
            ! Spread force from the immersed boundary
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0
            call spread_force_cilia_array(M,CAmid,Fx,Fy)
            call spread_force_cilia_array(M,CAPmid,Fx,Fy)
            call nvtxEndRange

            call nvtxStartRange("Predictor")
            call cdu(M,umid,vmid,us,nu,Fx)
            us = u + us*dt
            call cdv(M,umid,vmid,vs,nu,Fy)
            vs = v + vs*dt
            call nvtxEndRange
            
            call nvtxStartRange("Apply BC")
            ! Apply velocity boundary conditions to us and vs
            call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            ! call apply_parabolic_inlet(M,us,uleft)
            ! call apply_pulsating_inlet(M,u,uleft,t)
            call nvtxEndRange

            call nvtxStartRange("Calculate RHS")
            ! Form the RHS of the pressure poisson equation
            call calculate_rhs(M,us,vs,R,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Calculate pressure")
            ! Solve for pressure
            ! call calculate_pressure_amgx(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:),init_status)
            call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
            ! call calculate_pressure_sparse(A(1:M%Nx-1,:,:),P(1:M%Nx-1,:),R(1:M%Nx-1,:))
            call nvtxEndRange

            call nvtxStartRange("Corrector")
            ! Perform the corrector step to obtain the velocity
            call corrector(M,u,v,us,vs,p,rho,dt)
            call nvtxEndRange

            call nvtxStartRange("Interpolate velocity")
            ! Initialize the velocity at every time-step
            call initialize_velocity_cilia_array(CA)
            call initialize_velocity_cilia_array(CAP)
            ! Interpolate the Eulerian grid velocity to the Lagrangian structure
            call interpolate_velocity_cilia_array(M,CA,u,v)
            call interpolate_velocity_cilia_array(M,CAP,u,v)
            call nvtxEndRange

            call nvtxStartRange("Update cilia locations")
            ! Update the Immersed Boundary
            call update_cilia_array(CA,dt)
            call update_cilia_array(CAP,dt)
            call nvtxEndRange
        
            call nvtxEndRange
            write(*,'(A,F10.5)') 'time = ', t
            
        call nvtxEndRange
            ! Write files every Nth timestep
            if (mod(it,it_save).eq.0) then 
                call write_field(u,'u',it) 
                call write_field(v,'v',it) 
                call write_location_cilia(CA,it,'c')
                call write_location_cilia(CAP,it,'p')
                call write_location_cilia_force(CA,it,'c')
                call write_location_cilia_force(CAP,it,'p')
                call write_location_cilia_velocity(CA,it,'c')
                call write_location_cilia_velocity(CAP,it,'p')
                call write_cell_data(cell_array)
                call write_particle_data(cell_array)
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
                CA2%array(ic)%layers(il)%t = CA1%array(ic)%layers(il)%t
                do ip = 1,np
                    CA2%array(ic)%layers(il)%boundary(ip)%x = CA1%array(ic)%layers(il)%boundary(ip)%x
                    CA2%array(ic)%layers(il)%boundary(ip)%y = CA1%array(ic)%layers(il)%boundary(ip)%y
                    CA2%array(ic)%layers(il)%boundary(ip)%xo = CA1%array(ic)%layers(il)%boundary(ip)%xo
                    CA2%array(ic)%layers(il)%boundary(ip)%yo = CA1%array(ic)%layers(il)%boundary(ip)%yo
                    CA2%array(ic)%layers(il)%boundary(ip)%Ux = CA1%array(ic)%layers(il)%boundary(ip)%Ux
                    CA2%array(ic)%layers(il)%boundary(ip)%Uy = CA1%array(ic)%layers(il)%boundary(ip)%Uy
                    CA2%array(ic)%layers(il)%boundary(ip)%Fx = CA1%array(ic)%layers(il)%boundary(ip)%Fx
                    CA2%array(ic)%layers(il)%boundary(ip)%Fy = CA1%array(ic)%layers(il)%boundary(ip)%Fy
                end do
            end do
        end do

    end subroutine copy_cilia

end module mod_time_stepping
