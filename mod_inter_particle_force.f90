module mod_inter_particle_force
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use mod_cilia_array
    use mod_cell
    use mod_vec
    use mod_particle
    implicit none
   
    private
    public :: create_cells, assign_particles_to_cells, write_cell_data 
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
                cell_array(i,j)%loc = origin
                cell_array(i,j)%L = L
                origin%y= origin%y+L%y
            end do
            origin%x = origin%x+L%x
        end do         

    end subroutine create_cells

    subroutine assign_particles_to_cells(cell_array,cilia_all,particles_all)
        class(cilia_array), target, intent(in) :: cilia_all, particles_all
        class(cell), target, intent(inout) :: cell_array(:,:)

        integer(int32) :: Nxcell, Nycell

        integer(int32) :: i, j, k, q, r
        real(real64) :: px, py
       
        ! Check if the assumption that dim=1 is along is correct or wrong
        Nxcell = size(cell_array,1)
        Nycell = size(cell_array,2)

        ! Loop over the cilia arrays and sort them into the appropriate cell
        do i = 1,cilia_all%nc
            do j = 1,cilia_all%array(i)%nl
                do k = 1,cilia_all%array(i)%layers(j)%np
                    px = cilia_all%array(i)%layers(j)%boundary(k)%x
                    py = cilia_all%array(i)%layers(j)%boundary(k)%y
                    print *, px, py
                    ! find a box for this particle
                    outer: do q = 1,Nxcell
                        do r = 1,Nycell
                            if (     (px.ge.cell_array(q,r)%loc%x) &
                                .and.(py.ge.cell_array(q,r)%loc%y) &
                                .and.(px.lt.(cell_array(q,r)%loc%x + cell_array(q,r)%L%x)) &
                                .and.(py.lt.(cell_array(q,r)%loc%y + cell_array(q,r)%L%y))) then
                                    ! Increment the No. of neighbours for that cell
                                    print *, Nxcell, Nycell
                                    cell_array(q,r)%NN = cell_array(q,r)%NN+1 
                                    ! Assign the particle's pointer to the Neighbour list of that cell
                                    cell_array(q,r)%Nlist(cell_array(q,r)%NN)%pptr => cilia_all%array(i)%layers(j)%boundary(k)
                                    ! Attach a cilia Id to pointer data type
                                    cell_array(q,r)%Nlist(cell_array(q,r)%NN)%ciliaId = i 
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

        integer(int32) :: i, j, fileunit

        fileunit = 699
        open(unit=fileunit, file='cell_particles.txt', ACTION="write", STATUS="replace")
        do j = 1,size(cell_array,2)
            do i = 1,size(cell_array,1)
            write(fileunit, '(F10.5,F10.5,F10.5,F10.5)') cell_array(i,j)%loc%x, cell_array(i,j)%loc%y,&
                                                         cell_array(i,j)%L%x, cell_array(i,j)%L%y
            end do
        end do
        close(fileunit)
    end subroutine write_particle_data

end module mod_inter_particle_force