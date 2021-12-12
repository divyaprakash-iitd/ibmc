CC=gcc
FORT= gfortran
#FORT= /opt/nvidia/hpc_sdk/Linux_x86_64/2021/compilers/bin/nvfortran
#OPTS = -O3 -g -fbounds-check
#OPTS = -O3 -g -fcheck=bounds
OPTS = -O3 -fcheck=bounds
#OPTS = -O3 -g -stdpar -Minfo=accel
#OPTS = -fast -stdpar=gpu -Minfo=accel
#OPTS = -fast

OBJECTS= ibmc.o mod_pressure.o mod_amgx.o ftn_c.o mod_dims.o mod_mesh.o mod_time.o mod_boundary.o mod_time_stepping.o mod_io.o mod_particle.o mod_ib.o mod_vec.o mod_ibm.o mod_cilia.o mod_cilia_array.o mod_closed_cilia.o nvtx.o mod_cell.o mod_inter_particle_force.o

LIB_DIR= -L/home/divyaprakash/wrappers/amgx_code/axb_amgx

INC_DIR=

LIB_L= -laxbamgx -lnvToolsExt

ibmc: $(OBJECTS) 
	$(FORT) $(OPTS) -o $@ $(OBJECTS) $(LIB_DIR) $(LIB_L)

nvtx.o: nvtx.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_pressure.o: mod_pressure.f90 mod_mesh.o Makefile
	$(FORT) -c $(OPTS) $<

mod_amgx.o: mod_amgx.f90 ftn_c.o nvtx.o Makefile
	$(FORT) -c $(OPTS) $<

ftn_c.o: ftn_c.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_dims.o: mod_dims.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_mesh.o: mod_mesh.f90 mod_dims.o mod_vec.o Makefile
	$(FORT) -c $(OPTS) $<

mod_time.o: mod_time.f90 mod_mesh.o Makefile
	$(FORT) -c $(OPTS) $<

mod_boundary.o: mod_boundary.f90 mod_mesh.o Makefile
	$(FORT) -c $(OPTS) $<

mod_time_stepping.o: mod_time_stepping.f90 mod_pressure.o mod_amgx.o mod_mesh.o mod_dims.o mod_time.o mod_boundary.o mod_io.o mod_particle.o mod_ib.o mod_vec.o mod_ibm.o mod_cilia.o mod_cilia_array.o mod_closed_cilia.o nvtx.o mod_cell.o mod_inter_particle_force.o Makefile
	$(FORT) -c $(OPTS) $<

mod_io.o: mod_io.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_particle.o: mod_particle.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_ib.o: mod_ib.f90 mod_particle.o Makefile
	$(FORT) -c $(OPTS) $<

mod_vec.o: mod_vec.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_ibm.o: mod_ibm.f90 mod_mesh.o mod_ib.o mod_cilia.o mod_cilia_array.o Makefile
	$(FORT) -c $(OPTS) $<

mod_cilia.o: mod_cilia.f90 mod_ib.o Makefile
	$(FORT) -c $(OPTS) $<

mod_cilia_array.o: mod_cilia_array.f90 mod_cilia.o Makefile
	$(FORT) -c $(OPTS) $<

mod_closed_cilia.o: mod_closed_cilia.f90 mod_mesh.o mod_ib.o mod_vec.o mod_cilia.o mod_cilia_array.o mod_ibm.o Makefile
	$(FORT) -c $(OPTS) $<

mod_cell.o: mod_cell.f90 mod_vec.o Makefile
	$(FORT) -c $(OPTS) $<

mod_inter_particle_force.o: mod_inter_particle_force.f90 mod_mesh.o mod_cilia_array.o mod_cell.o mod_vec.o Makefile
	$(FORT) -c $(OPTS) $<

ibmc.o: ibmc.f90 mod_pressure.o mod_amgx.o mod_mesh.o mod_dims.o mod_time.o mod_boundary.o mod_time_stepping.o mod_io.o mod_particle.o mod_ib.o mod_vec.o mod_ibm.o mod_cilia.o mod_cilia_array.o mod_closed_cilia.o mod_cell.o mod_inter_particle_force.o Makefile
	$(FORT) -c $(OPTS) $<

clean:
	rm -f ibmc *.mod *.bak *~ $(OBJECTS) 
