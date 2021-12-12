module mod_inter_particle_force
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use mod_cilia_array
    use mod_cell
    use mod_vec
    use mod_particle
    implicit none
   
    private
    public :: create_cells, assign_particles_to_cells, write_cell_data, write_particle_data, calculate_forces 
contains
    subroutine create_cells(M,cell_array)
        class(mesh), intent(in) :: M
        class(cell), target, intent(inout) :: cell_array(:,:)

        integer(int32) :: Nxcell, Nycell
        type(vec) :: L
        type(vec) :: origin
        integer(int32) :: i, j

        Nxcell = size(cell_array,1)
        Nycell = size(cell_array,2)

        L%x = M%Lx/Nxcell
        L%y = M%Ly/Nycell

        ! Divide the domain into cells and assign the locations and dimensions to each cell
        
        origin%x = 0.0d0
        do i = 1,Nxcell
            origin%y = 0.0d0
            do j = 1,Nycell
                cell_array(i,j)%cellId = i + (j-1)*Nxcell
                cell_array(i,j)%loc%x = origin%x
                cell_array(i,j)%loc%y = origin%y
                cell_array(i,j)%L%x = L%x
                cell_array(i,j)%L%y = L%y
                origin%y= origin%y+L%y
            end do
            origin%x = origin%x+L%x
        end do         

    end subroutine create_cells

    subroutine assign_particles_to_cells(cell_array,cilia_all)
        class(cilia_array), target, intent(in) :: cilia_all
        class(cell), target, intent(inout) :: cell_array(:,:)

        integer(int32) :: Nxcell, Nycell

        integer(int32) :: i, j, k, q, r, t
        real(real64) :: px, py
       
        ! Check if the assumption that dim=1 is along is correct or wrong
        Nxcell = size(cell_array,1)
        Nycell = size(cell_array,2)

        ! Check if cilia or particle
        if (cilia_all%array(1)%layers(1)%t.eq.'o') then
            t = 1 ! Open cilia or Cilia
        elseif (cilia_all%array(1)%layers(1)%t.eq.'c') then
            t = 2 ! Closed cilia or particle
        end if

        ! Loop over the cilia arrays and sort them into the appropriate cell
        do i = 1,cilia_all%nc ! loop over the no. of cilia
            do j = 1,cilia_all%array(i)%nl ! loop over the no. layers per cilia
                do k = 1,cilia_all%array(i)%layers(j)%np ! loop over the no. of nodes/particles per layer
                    px = cilia_all%array(i)%layers(j)%boundary(k)%x
                    py = cilia_all%array(i)%layers(j)%boundary(k)%y
                    ! print *, px, py
                    ! find a box for this particle
                    outer: do q = 1,Nxcell
                        do r = 1,Nycell
                            if (     (px.ge.cell_array(q,r)%loc%x) &
                                .and.(py.ge.cell_array(q,r)%loc%y) &
                                .and.(px.lt.(cell_array(q,r)%loc%x + cell_array(q,r)%L%x)) &
                                .and.(py.lt.(cell_array(q,r)%loc%y + cell_array(q,r)%L%y))) then
                                    ! Increment the No. of neighbours for that cell
                                    ! print *, Nxcell, Nycell
                                    cell_array(q,r)%NN = cell_array(q,r)%NN+1 
                                    ! Assign the particle's pointer to the Neighbour list of that cell
                                    cell_array(q,r)%Nlist(cell_array(q,r)%NN)%pptr => cilia_all%array(i)%layers(j)%boundary(k)
                                    ! Attach a cilia Id to pointer data type
                                    cell_array(q,r)%Nlist(cell_array(q,r)%NN)%ciliaId(1) = t 
                                    cell_array(q,r)%Nlist(cell_array(q,r)%NN)%ciliaId(2) = i 
                                    exit outer
                            end if
                        end do
                    end do outer
                end do
            end do
        end do 
        
    end subroutine assign_particles_to_cells

    subroutine write_cell_data(cell_array)
        class(cell), target, intent(inout) :: cell_array(:,:)
       
        integer(int32) :: i, j, fileunit

        fileunit = 699
        open(unit=fileunit, file='cell.txt', ACTION="write", STATUS="replace")
        do j = 1,size(cell_array,2)
            write(fileunit, '(F10.5,F10.5,F10.5,F10.5)')( cell_array(i,j)%loc%x, cell_array(i,j)%loc%y,&
                                          cell_array(i,j)%L%x, cell_array(i,j)%L%y , i = 1,size(cell_array,1))
        end do
        close(fileunit)
    end subroutine write_cell_data

    subroutine write_particle_data(cell_array)
        class(cell), intent(inout) :: cell_array(:,:)

        integer(int32) :: i, j, k, fileunit

        fileunit = 699
        open(unit=fileunit, file='cell_particles.txt', ACTION="write", STATUS="replace")
        do j = 1,size(cell_array,2)
            do i = 1,size(cell_array,1)
                if (cell_array(i,j)%NN.gt.0) then
                    write(fileunit, '(*(F14.7))')( cell_array(i,j)%Nlist(k)%pptr%x , k = 1,cell_array(i,j)%NN)
                    write(fileunit, '(*(F14.7))')( cell_array(i,j)%Nlist(k)%pptr%y , k = 1,cell_array(i,j)%NN)
                end if
            end do
        end do
        close(fileunit)
    end subroutine write_particle_data

    subroutine calculate_forces(cell_array)
        class(cell), intent(inout) :: cell_array(:,:)

        integer(int32) :: i, j, NNx, NNy

        NNx = size(cell_array,1)
        NNy = size(cell_array,2)

        do j = 1,NNy
            do i = 1,NNx
                ! ! Calculate same cell forces
                call calculate_force_neighbouring_cells(cell_array(i,j),cell_array(i,j))

                ! Right
                if ((i+1).le.NNx) then
                    call calculate_force_neighbouring_cells(cell_array(i,j),cell_array(i+1,j))
                end if
                
                ! Top
                if ((j+1).le.NNy) then
                    call calculate_force_neighbouring_cells(cell_array(i,j),cell_array(i,j+1))
                end if

                ! Top Right
                if (((i+1).le.NNx).and.((j+1).le.NNy)) then
                    call calculate_force_neighbouring_cells(cell_array(i,j),cell_array(i+1,j+1))
                end if
                
                ! Top Left
                if (((i-1).ge.1).and.((j+1).le.NNy)) then
                    call calculate_force_neighbouring_cells(cell_array(i,j),cell_array(i-1,j+1))
                end if

            end do
        end do
    end subroutine calculate_forces

    subroutine calculate_force_neighbouring_cells(masterC, slaveC)
        class(cell), intent(inout) :: masterC
        class(cell), intent(inout) :: slaveC

        integer(int32) :: i, j
        real(real64) :: xm, ym, xsl, ysl, d, Fmx, Fmy, Fslx, Fsly, dcutoff, k

        ! Spring constants
        k = 0.1

        ! Cut-off distnace
        dcutoff = 0.1


        do i = 1,masterC%NN
                ! Master particle location
                xm = masterC%Nlist(i)%pptr%x
                ym = masterC%Nlist(i)%pptr%y
            do j = 1,slaveC%NN
                ! Check for cilia type
                if (masterC%Nlist(i)%ciliaId(1).ne.slaveC%Nlist(j)%ciliaId(1) &
                ! Check for cilia Id    
                .or.(masterC%Nlist(i)%ciliaId(2).ne.slaveC%Nlist(j)%ciliaId(2))) then 
                    ! Calculate the distance between the two particles

                    ! print *, masterC%NN
                    ! Slave particle location
                    xsl = slaveC%Nlist(j)%pptr%x
                    ysl = slaveC%Nlist(j)%pptr%y

                    ! Calculate the distance between the master and cell particle
                    d = norm2([(xsl-xm),(ysl-ym)])

                    ! Check for cut-off distance
                    if (abs(d).le.dcutoff) then
                        ! Calculate forces (Master particle)
                        if (((xsl-xm).lt.1e-10).and.((ysl-ym).lt.1e-10)) then
                            Fmx = 0.0d0
                            Fmy = 0.0d0
                        else
                            Fmx = k*(1.0d0-dcutoff/d)*(xsl-xm)
                            Fmy = k*(1.0d0-dcutoff/d)*(ysl-ym)
                        end if
                        ! Calculate forces (Slave particle)
                        Fslx = -Fmx
                        Fsly = -Fmy
                        ! print *, xsl-xm
                        ! Assign the forces to the respective pointers
                        ! Master cell pointers
                        masterC%Nlist(i)%pptr%Fx = masterC%Nlist(i)%pptr%Fx + Fmx
                        masterC%Nlist(i)%pptr%Fy = masterC%Nlist(i)%pptr%Fy + Fmy
                        ! Slave cell pointers
                        slaveC%Nlist(j)%pptr%Fx = slaveC%Nlist(j)%pptr%Fx + Fslx
                        slaveC%Nlist(j)%pptr%Fy = slaveC%Nlist(j)%pptr%Fy + Fsly
                    end if
                end if                
            end do
        end do
    end subroutine calculate_force_neighbouring_cells

    subroutine calculate_inter_particle_force(cell_array)
        class(cell), intent(inout) :: cell_array(:,:)

        integer(int32) :: i, j, p, q, m, n
        real(real64) :: Fx, Fy

        ! For every cell, we look at the (i+1,j) and (i,j+1) cells as its neighbours
        do j = 1,size(cell_array,2)
            do i = 1,size(cell_array,1)
                ! Loop over every particle
                do p = 1,cell_array(i,j)%NN
                    do q = 1,cell_array(i,j)%NN
                        ! Check if the particles isn't being compared against itself
                        if (p.ne.q) then
                            ! Check for Cilia type
                            if ((cell_array(i,j)%Nlist(p)%ciliaId(1).ne.cell_array(i,j)%Nlist(q)%ciliaId(1)) & 
                            ! Check for Cilia Id
                            .or.(cell_array(i,j)%Nlist(p)%ciliaId(2).ne.cell_array(i,j)%Nlist(q)%ciliaId(2))) then
                                ! Calculate pair force and assign to the particle
                                Fx = 0.0d0
                                Fy = 0.0d0
                                cell_array(i,j)%Nlist(p)%pptr%Fx = cell_array(i,j)%Nlist(p)%pptr%Fx + Fx 
                                cell_array(i,j)%Nlist(p)%pptr%Fy = cell_array(i,j)%Nlist(p)%pptr%Fy + Fy
                                Fx = -Fx 
                                Fy = -Fy
                                cell_array(i,j)%Nlist(q)%pptr%Fx = cell_array(i,j)%Nlist(q)%pptr%Fx + Fx 
                                cell_array(i,j)%Nlist(q)%pptr%Fy = cell_array(i,j)%Nlist(q)%pptr%Fy + Fy
                            end if
                            ! Also check in neighbouring cells
                            if ((i+1).le.size(cell_array,1)) then
                                ! Check for Cilia type
                                if ((cell_array(i+1,j)%Nlist(p)%ciliaId(1).ne.cell_array(i+1,j)%Nlist(q)%ciliaId(1)) & 
                                ! Check for Cilia Id
                                .or.(cell_array(i+1,j)%Nlist(p)%ciliaId(2).ne.cell_array(i+1,j)%Nlist(q)%ciliaId(2))) then
                                    ! Calculate pair force and assign to the particle
                                    Fx = 0.0d0
                                    Fy = 0.0d0
                                    cell_array(i+1,j)%Nlist(p)%pptr%Fx = cell_array(i+1,j)%Nlist(p)%pptr%Fx + Fx 
                                    cell_array(i+1,j)%Nlist(p)%pptr%Fy = cell_array(i+1,j)%Nlist(p)%pptr%Fy + Fy
                                    Fx = -Fx 
                                    Fy = -Fy
                                    cell_array(i+1,j)%Nlist(q)%pptr%Fx = cell_array(i+1,j)%Nlist(q)%pptr%Fx + Fx 
                                    cell_array(i+1,j)%Nlist(q)%pptr%Fy = cell_array(i+1,j)%Nlist(q)%pptr%Fy + Fy
                                end if
                            end if

                            if ((j+1).le.size(cell_array,2)) then
                                ! Check for Cilia type
                                if ((cell_array(i,j+1)%Nlist(p)%ciliaId(1).ne.cell_array(i,j+1)%Nlist(q)%ciliaId(1)) & 
                                ! Check for Cilia Id
                                .or.(cell_array(i,j+1)%Nlist(p)%ciliaId(2).ne.cell_array(i,j+1)%Nlist(q)%ciliaId(2))) then
                                    ! Calculate pair force and assign to the particle
                                    Fx = 0.0d0
                                    Fy = 0.0d0
                                    cell_array(i,j+1)%Nlist(p)%pptr%Fx = cell_array(i,j+1)%Nlist(p)%pptr%Fx + Fx 
                                    cell_array(i,j+1)%Nlist(p)%pptr%Fy = cell_array(i,j+1)%Nlist(p)%pptr%Fy + Fy
                                    Fx = -Fx 
                                    Fy = -Fy
                                    cell_array(i,j+1)%Nlist(q)%pptr%Fx = cell_array(i,j+1)%Nlist(q)%pptr%Fx + Fx 
                                    cell_array(i,j+1)%Nlist(q)%pptr%Fy = cell_array(i,j+1)%Nlist(q)%pptr%Fy + Fy
                                end if
                            end if

                        end if
                    end do
                end do 
            end do
        end do

    end subroutine calculate_inter_particle_force

end module mod_inter_particle_force