CC=gcc
#FORT= gfortran
FORT= /opt/nvidia/hpc_sdk/Linux_x86_64/2021/compilers/bin/nvfortran
#OPTS = -O3 -g -fbounds-check
#OPTS = -O3 -g -fcheck=bounds
OPTS = -O3 -g

OBJECTS=ibmc.o mod_pressure.o

LIB_DIR=

INC_DIR=

LIB_L=

ibmc: $(OBJECTS) 
	$(FORT) $(OPTS) -o $@ $(OBJECTS) $(LIB_DIR) $(LIB_L)

mod_pressure.o: mod_pressure.f90 Makefile
	$(FORT) -c $(OPTS) $<

ibmc.o: ibmc.f90 mod_pressure.o Makefile
	$(FORT) -c $(OPTS) $<

clean:
	rm -f ibmc *.mod *.bak *~ $(OBJECTS) 
