# Makefile for building Green function
# 
# Modified by Shuying SUN for EXL 2019/05/24 
#
MACHINE := $(shell uname)
VPATH =../src
RUN_DIR = ../run

ifeq ($(MACHINE), Darwin)
  FC = gfortran
  FFLAGS = -g
  LDFLAGS =
  LIB =
endif
ifeq ($(MACHINE),Linux)
  FC = pgf90
  FFLAGS = -g
  LDFLAGS = -Mlfs
  LIB =
endif

SOURCES = modules.f90 subroutines.f90 main.f90 getset.f90 calcgfunc.f90 \
					writegfunc.f90

OBJECTS=$(SOURCES:.f90=.o)

gfun:$(OBJECTS) 
	$(FC) $(LDFLAGS) -o $@ $^ $(LIB)
	cp $@ $(RUN_DIR)
	@echo "GFUN BUILT"

%.o:%.f90
	$(FC) $(FFLAGS) -c $< $(LIB)

$(filter-out modules.o,$(OBJECTS)):modules.o

main.o: main.f90 subroutines.o getset.o calcgfunc.o writegfunc.o
	$(FC) $(FFLAGS) -c $< $(LIB)

#modules.o: modules.f90
#	$(FC) $(FFLAGS) -c $< $(LIB)
#
#subroutines.o: subroutines.f90 module.o
#	$(FC) $(FFLAGS) -c $< $(LIB)
#
#getset.o: getset.f90 module.o
#	$(FC) $(FFLAGS) -c $< $(LIB)
#
#calcgfunc.o: calcgfunc.f90 module.o
#	$(FC) $(FFLAGS) -c $< $(LIB)
#
#writegfunc.o: writegfunc.f90 module.o
#	$(FC) $(FFLAGS) -c $< $(LIB)

clean:
	rm *.o *.mod gfun
	@echo "CLEAN DONE"

