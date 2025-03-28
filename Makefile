CC=gcc
FORT= gfortran
#FORT= /opt/nvidia/hpc_sdk/Linux_x86_64/2021/compilers/bin/nvfortran
#OPTS = -O3 -g -fbounds-check
OPTS = -O3 -g -fcheck=bounds
#OPTS = -O3 -g -stdpar -Minfo=accel

OBJECTS= ibmc.o mod_pressure.o mod_amgx.o ftn_c.o

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

ibmc.o: ibmc.f90 mod_pressure.o mod_amgx.o Makefile
	$(FORT) -c $(OPTS) $<

clean:
	rm -f ibmc *.mod *.bak *~ $(OBJECTS) 
