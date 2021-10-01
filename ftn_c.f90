module ftn_c
    interface
        integer (c_int) function initialize_amgx(crs_data, datam, col_ind, row_ptr, rhs, sol) bind(c, name='initialize_amgx')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: crs_data
            type (c_ptr), value :: datam
            type (c_ptr), value :: col_ind
            type (c_ptr), value :: row_ptr
            type (c_ptr), value :: rhs
            type (c_ptr), value :: sol
        end function initialize_amgx
        
        integer (c_int) function solveamg(rhs, sol) bind(c, name='solveamg')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: rhs
            type (c_ptr), value :: sol
        end function solveamg
        
        integer (c_int) function destroy_amgx() bind(c, name='destroy_amgx')
            use iso_c_binding
            implicit none
        end function destroy_amgx
    end interface    
end module ftn_c
