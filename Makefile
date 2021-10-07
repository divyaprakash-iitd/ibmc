CC=gcc
#FORT= gfortran
FORT= /opt/nvidia/hpc_sdk/Linux_x86_64/2021/compilers/bin/nvfortran
#OPTS = -O3 -g -fbounds-check
#OPTS = -O3 -g -fcheck=bounds
OPTS = -O3 -g -stdpar -Minfo=accel

OBJECTS= ibmc.o mod_pressure.o mod_amgx.o ftn_c.o mod_dims.o mod_mesh.o mod_time.o mod_boundary.o mod_time_stepping.o mod_io.o mod_particle.o mod_ib.o

LIB_DIR= -L/home/divyaprakash/wrappers/amgx_code/axb_amgx


INC_DIR=

LIB_L= -laxbamgx


ibmc: $(OBJECTS) 
	$(FORT) $(OPTS) -o $@ $(OBJECTS) $(LIB_DIR) $(LIB_L)

mod_pressure.o: mod_pressure.f90 Makefile
	$(FORT) -c $(OPTS) $<

mod_amgx.o: mod_amgx.f90 ftn_c.o Makefile
	$(FORT) -c $(OPTS) $<

ftn_c.o: ftn_c.f90 Makefile
	$(FORT) -c $<

mod_dims.o: mod_dims.f90 Makefile
	$(FORT) -c $<

mod_mesh.o: mod_mesh.f90 mod_dims.o Makefile
	$(FORT) -c $<

mod_time.o: mod_time.f90 mod_mesh.o Makefile
	$(FORT) -c $<

mod_boundary.o: mod_boundary.f90 mod_mesh.o Makefile
	$(FORT) -c $<

mod_time_stepping.o: mod_time_stepping.f90 Makefile
	$(FORT) -c $<

mod_io.o: mod_io.f90 Makefile
	$(FORT) -c $<

mod_particle.o: mod_particle.f90 Makefile
	$(FORT) -c $<

mod_ib.o: mod_ib.f90 mod_particle.o Makefile
	$(FORT) -c $<

ibmc.o: ibmc.f90 mod_pressure.o mod_amgx.o mod_mesh.o mod_dims.o mod_time.o mod_boundary.o mod_time_stepping.o mod_io.o mod_particle.o mod_ib.o Makefile
	$(FORT) -c $(OPTS) $<

clean:
	rm -f ibmc *.mod *.bak *~ $(OBJECTS) 
